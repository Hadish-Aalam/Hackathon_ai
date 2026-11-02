"""Typed configuration; all numerical values use SI units unless named otherwise."""
from dataclasses import dataclass
import json
import math
from pathlib import Path

@dataclass
class MotorParameters:
    resistance_ohm: float
    ld_h: float
    lq_h: float
    flux_wb: float
    pole_pairs: int
    inertia_kg_m2: float
    viscous_nm_s: float
    coulomb_nm: float
    copper_alpha_per_c: float
    reference_temp_c: float
    cable_ohm_per_m: float
    cable_length_m: float
    thermal_capacity_j_per_c: float
    thermal_resistance_c_per_w: float

@dataclass
class ControlParameters:
    dc_link_v: float
    voltage_reserve_fraction: float
    max_current_a: float
    max_negative_id_a: float
    current_kp: float
    current_ki: float
    speed_kp: float
    speed_ki: float
    position_kp: float
    antiwindup_gain_per_s: float
    derating_start_c: float
    shutdown_c: float
    field_weakening_controller: bool
    field_weakening_motor: bool

    @property
    def voltage_limit(self):
        return self.dc_link_v / 2**0.5 * (1 - self.voltage_reserve_fraction)

@dataclass
class SimulationParameters:
    duration_s: float
    step_s: float
    sample_every: int
    ambient_c: float
    initial_temp_c: float
    initial_speed_rad_s: float
    initial_position_rad: float
    mode: str
    target_position_rad: float
    max_speed_rad_s: float
    speed_profile: list
    load_step_s: float
    load_torque_nm: float
    brake_start_s: float
    brake_torque_nm: float

@dataclass
class Configuration:
    motor: MotorParameters
    control: ControlParameters
    simulation: SimulationParameters

    @classmethod
    def load(cls, path):
        raw = json.loads(Path(path).read_text())
        cfg = cls(MotorParameters(**raw['motor']), ControlParameters(**raw['control']),
                  SimulationParameters(**raw['simulation']))
        cfg.validate()
        return cfg

    def validate(self):
        m, c, s = self.motor, self.control, self.simulation
        def finite(value):
            if isinstance(value, (int, float)) and not math.isfinite(value):
                raise ValueError('All numerical parameters must be finite')
            if isinstance(value, list):
                for item in value: finite(item)
        for section in (m, c, s):
            for value in vars(section).values(): finite(value)
        if int(m.pole_pairs) != m.pole_pairs:
            raise ValueError('pole_pairs must be an integer')
        if any(value < 0 for value in (m.viscous_nm_s, m.coulomb_nm, m.cable_ohm_per_m,
                                      m.cable_length_m, c.current_kp, c.current_ki,
                                      c.speed_kp, c.speed_ki, c.position_kp, c.antiwindup_gain_per_s)):
            raise ValueError('Loss parameters and controller gains must be nonnegative')
        if 1 + m.copper_alpha_per_c*(min(s.initial_temp_c,s.ambient_c)-m.reference_temp_c) <= 0:
            raise ValueError('Temperature produces nonphysical resistance')
        for value in (m.resistance_ohm, m.ld_h, m.lq_h, m.flux_wb, m.pole_pairs,
                      m.inertia_kg_m2, m.thermal_capacity_j_per_c, m.thermal_resistance_c_per_w,
                      c.dc_link_v, c.max_current_a, s.step_s, s.duration_s):
            if value <= 0:
                raise ValueError('Motor, supply and integration scales must be positive')
        if not 0 <= c.voltage_reserve_fraction < 1:
            raise ValueError('Voltage reserve must lie in [0, 1)')
        if not 0 <= c.max_negative_id_a <= c.max_current_a:
            raise ValueError('Negative d-current limit must lie within current rating')
        if c.shutdown_c <= c.derating_start_c:
            raise ValueError('Shutdown temperature must exceed derating start')
        if s.mode not in ('speed', 'position') or s.sample_every < 1 or int(s.sample_every) != s.sample_every:
            raise ValueError('Use speed/position mode and a positive integer sample_every')
        times = [p[0] for p in s.speed_profile]
        if not times or times[0] != 0 or any(b <= a for a,b in zip(times, times[1:])):
            raise ValueError('Speed profile must start at zero with strictly increasing times')
        if s.step_s > min(m.ld_h, m.lq_h) / m.resistance_ohm / 10:
            raise ValueError('Reduce integration step to resolve electrical dynamics')
