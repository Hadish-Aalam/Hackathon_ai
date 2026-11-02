import unittest
from copy import deepcopy
from pathlib import Path
import numpy as np
from motor_sim.config import Configuration
from motor_sim.model import Motor, FieldWeakening, Simulation, TwoPhaseTransform

ROOT = Path(__file__).resolve().parents[1]

class ModelTests(unittest.TestCase):
    def setUp(self):
        self.cfg = Configuration.load(ROOT/'config.json')

    def test_park_roundtrip_and_power(self):
        for angle in np.linspace(-7,7,19):
            u = np.array([3.,-4.]); i = np.array([2.,5.])
            dq = TwoPhaseTransform.to_dq(*u,angle)
            np.testing.assert_allclose(TwoPhaseTransform.to_phase(*dq,angle),u,atol=1e-12)
            self.assertAlmostEqual(float(u@i),float(dq@TwoPhaseTransform.to_dq(*i,angle)))

    def test_instantaneous_energy_balance(self):
        motor = Motor(self.cfg.motor,25)
        p = motor.p
        state = np.array([-3.,2.,180.,4.,40.]); voltage = np.array([-8.,24.])
        load, brake = .08,.1
        dx = motor.derivative(state,voltage,load,brake)
        id_,iq,w,_,temp = state
        stored_rate = p.ld_h*id_*dx[0]+p.lq_h*iq*dx[1]+p.inertia_kg_m2*w*dx[2]
        losses = motor.resistance(temp)*(id_*id_+iq*iq)
        friction = (p.coulomb_nm+brake)*np.tanh(w/.1)+p.viscous_nm_s*w
        self.assertAlmostEqual(float(voltage@state[:2]),stored_rate+losses+(load+friction)*w,places=9)

    def test_field_weakening_voltage_solution_and_flags(self):
        motor = Motor(self.cfg.motor,25); c = self.cfg.control; p = motor.p
        fw = FieldWeakening(motor,c); r = motor.resistance(25)
        self.assertEqual(fw.reference(20,1,r),0)
        id_ = fw.reference(260,1,r); we = 260*p.pole_pairs
        self.assertLess(id_,0); self.assertGreaterEqual(id_,-c.max_negative_id_a)
        self.assertAlmostEqual(np.hypot(r*id_-we*p.lq_h, r+we*(p.ld_h*id_+p.flux_wb)),c.voltage_limit)
        c.field_weakening_motor = False
        self.assertEqual(fw.reference(260,1,r),0)

    def test_closed_loop_limits_and_timestep_convergence(self):
        cfg = self.cfg; cfg.simulation.duration_s = .6
        a = Simulation(cfg).run()
        fine = deepcopy(cfg); fine.simulation.step_s /= 2; fine.simulation.sample_every *= 2
        b = Simulation(fine).run()
        np.testing.assert_allclose(a['speed_rad_s'],b['speed_rad_s'],atol=.15,rtol=.003)
        self.assertLessEqual(np.hypot(a['vd_v'],a['vq_v']).max(),cfg.control.voltage_limit+1e-10)
        self.assertTrue(np.all(abs(a['id_ref_a'])+abs(a['iq_ref_a']) <= a['current_limit_a']+1e-10))
        self.assertAlmostEqual(a['speed_rad_s'][-1],80,delta=1)

    def test_brake_and_thermal_derating(self):
        from motor_sim.model import DriveController
        motor = Motor(self.cfg.motor,25)
        state = np.array([0.,1.,20.,0.,25.])
        free = motor.derivative(state,[0.,0.],0,0)
        braked = motor.derivative(state,[0.,0.],0,.3)
        self.assertLess(braked[2], free[2])
        state[4] = self.cfg.control.shutdown_c
        controller = DriveController(motor,self.cfg.control)
        _,dr,qr,limit = controller.step(state,100,.0001)
        self.assertEqual((dr,qr,limit),(0,0,0))

    def test_position_loop(self):
        cfg = self.cfg
        cfg.simulation.mode = 'position'
        cfg.simulation.target_position_rad = 1
        cfg.simulation.duration_s = 1
        cfg.simulation.load_step_s = 2
        data = Simulation(cfg).run()
        self.assertAlmostEqual(data['position_rad'][-1],1,delta=.03)

    def test_invalid_configuration(self):
        self.cfg.control.voltage_reserve_fraction = 1
        with self.assertRaises(ValueError): self.cfg.validate()

if __name__ == '__main__': unittest.main()
