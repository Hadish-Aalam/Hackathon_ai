"""Power-invariant orthogonal two-phase dq model (no three-phase 3/2 factor)."""
import math
import numpy as np

class PIController:
    def __init__(self, kp, ki, kaw):
        self.kp, self.ki, self.kaw, self.integral = kp, ki, kaw, 0.0

    def output(self, error):
        return self.kp * error + self.integral

    def update(self, error, requested, applied, dt):
        self.integral += dt * (self.ki * error + self.kaw * (applied - requested))

class Motor:
    def __init__(self, parameters, ambient):
        self.p, self.ambient = parameters, ambient

    def resistance(self, temperature):
        p = self.p
        return p.resistance_ohm * (1 + p.copper_alpha_per_c * (temperature-p.reference_temp_c)) + p.cable_ohm_per_m*p.cable_length_m

    def torque(self, id_, iq):
        p = self.p
        return p.pole_pairs*(p.flux_wb + (p.ld_h-p.lq_h)*id_)*iq

    def derivative(self, state, voltage, load, brake):
        id_, iq, omega, theta, temp = state
        vd, vq = voltage
        p = self.p
        we, r = p.pole_pairs*omega, self.resistance(temp)
        torque = self.torque(id_, iq)
        friction = (p.coulomb_nm + brake)*math.tanh(omega/0.1) + p.viscous_nm_s*omega
        copper = r*(id_*id_+iq*iq)
        return np.array([(vd-r*id_+we*p.lq_h*iq)/p.ld_h,
                         (vq-r*iq-we*(p.ld_h*id_+p.flux_wb))/p.lq_h,
                         (torque-load-friction)/p.inertia_kg_m2, omega,
                         (copper-(temp-self.ambient)/p.thermal_resistance_c_per_w)/p.thermal_capacity_j_per_c])

    def advance(self, state, voltage, load, brake, dt):
        # Fixed control period, zero-order-held inverter voltage, RK4 plant integration.
        f = lambda x: self.derivative(x, voltage, load, brake)
        a = f(state); b = f(state+dt*a/2); c = f(state+dt*b/2); d = f(state+dt*c)
        return state + dt*(a+2*b+2*c+d)/6

class FieldWeakening:
    def __init__(self, motor, control):
        self.motor, self.c = motor, control

    def reference(self, omega, iq, resistance):
        c, p = self.c, self.motor.p
        if not (c.field_weakening_controller and c.field_weakening_motor):
            return 0.0
        # Exact steady-state voltage ellipse: (R id-we Lq iq)^2 +
        # (R iq+we(Ld id+psi))^2 <= Vavailable^2.
        we = p.pole_pairs*omega
        a = resistance**2 + (we*p.ld_h)**2
        b = 2*(-resistance*we*p.lq_h*iq + we*p.ld_h*(resistance*iq+we*p.flux_wb))
        cc = (we*p.lq_h*iq)**2+(resistance*iq+we*p.flux_wb)**2-c.voltage_limit**2
        if cc <= 0:
            return 0.0
        disc = b*b-4*a*cc
        candidate = (-b+math.sqrt(disc))/(2*a) if disc >= 0 else -b/(2*a)
        return float(np.clip(candidate, -c.max_negative_id_a, 0))

class DriveController:
    def __init__(self, motor, control):
        self.motor, self.c = motor, control
        self.speed = PIController(control.speed_kp, control.speed_ki, control.antiwindup_gain_per_s)
        self.d = PIController(control.current_kp, control.current_ki, control.antiwindup_gain_per_s)
        self.q = PIController(control.current_kp, control.current_ki, control.antiwindup_gain_per_s)
        self.fw = FieldWeakening(motor, control)

    def step(self, state, reference, dt):
        id_, iq, omega, _, temp = state
        c, p = self.c, self.motor.p
        limit = c.max_current_a * float(np.clip((c.shutdown_c-temp)/(c.shutdown_c-c.derating_start_c), 0, 1))
        error = reference-omega
        raw_iq = self.speed.output(error)
        iq_ref = float(np.clip(raw_iq, -limit, limit))
        r = self.motor.resistance(temp)
        # Reference speed enables anticipatory weakening; magnitude chooses direction from actual motion.
        fw_speed = reference if abs(reference) > abs(omega) else omega
        id_ref = max(-limit, self.fw.reference(fw_speed, iq_ref, r))
        # Conservative L1 reference bound specified in the document.
        iq_limit = max(0, limit-abs(id_ref))
        iq_ref = float(np.clip(iq_ref, -iq_limit, iq_limit))
        self.speed.update(error, raw_iq, iq_ref, dt)
        ed, eq = id_ref-id_, iq_ref-iq
        we = p.pole_pairs*omega
        raw_d = self.d.output(ed) - we*p.lq_h*iq
        raw_q = self.q.output(eq) + we*(p.ld_h*id_+p.flux_wb)
        scale = min(1, c.voltage_limit/max(math.hypot(raw_d,raw_q), 1e-12))
        vd, vq = raw_d*scale, raw_q*scale
        self.d.update(ed, raw_d, vd, dt)
        self.q.update(eq, raw_q, vq, dt)
        return np.array([vd,vq]), id_ref, iq_ref, limit

class Simulation:
    columns = ['time_s','speed_ref_rad_s','speed_rad_s','position_rad','id_a','iq_a',
               'id_ref_a','iq_ref_a','vd_v','vq_v','torque_nm','load_nm','temperature_c',
               'electrical_power_w','copper_loss_w','mechanical_power_w','current_limit_a']

    def __init__(self, config):
        self.config = config
        self.motor = Motor(config.motor, config.simulation.ambient_c)
        self.controller = DriveController(self.motor, config.control)

    def run(self):
        s, c = self.config.simulation, self.config.control
        state = np.array([0.,0.,s.initial_speed_rad_s,s.initial_position_rad,s.initial_temp_c])
        rows = []
        count = math.ceil(s.duration_s/s.step_s)
        for k in range(count+1):
            t = min(k*s.step_s, s.duration_s)
            ref = s.speed_profile[0][1]
            for when, value in s.speed_profile:
                if t >= when: ref = value
            if s.mode == 'position':
                ref = float(np.clip(c.position_kp*(s.target_position_rad-state[3]), -s.max_speed_rad_s, s.max_speed_rad_s))
            dt = min(s.step_s, max(0, s.duration_s-t))
            voltage, dr, qr, limit = self.controller.step(state, ref, dt)
            load = s.load_torque_nm if t >= s.load_step_s else 0.
            brake = s.brake_torque_nm if t >= s.brake_start_s else 0.
            id_, iq, omega, theta, temp = state
            torque = self.motor.torque(id_,iq)
            if k % s.sample_every == 0 or k == count:
                rows.append([t,ref,omega,theta,id_,iq,dr,qr,*voltage,torque,load,temp,
                             voltage[0]*id_+voltage[1]*iq,self.motor.resistance(temp)*(id_*id_+iq*iq),torque*omega,limit])
            if dt: state = self.motor.advance(state, voltage, load, brake, dt)
            if not np.all(np.isfinite(state)):
                raise RuntimeError('Simulation diverged; check parameters and integration step')
        return {name: np.array(rows)[:,i] for i,name in enumerate(self.columns)}


class TwoPhaseTransform:
    """Orthogonal Park transform; two-phase alpha/beta require no Clarke stage."""
    @staticmethod
    def to_dq(alpha, beta, electrical_angle):
        c, s = math.cos(electrical_angle), math.sin(electrical_angle)
        return np.array([c*alpha+s*beta, -s*alpha+c*beta])

    @staticmethod
    def to_phase(d, q, electrical_angle):
        c, s = math.cos(electrical_angle), math.sin(electrical_angle)
        return np.array([c*d-s*q, s*d+c*q])
