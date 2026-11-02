# Electromechanical motor simulation in Python

An object-oriented, averaged two-phase permanent-magnet motor simulation reconstructed from **Simulink Model Description with FW.pdf**. It includes cascaded speed/current PI control, optional position control, field weakening, electrical and mechanical dynamics, cable resistance, temperature-dependent resistance, thermal derating, load, and braking.

**This is an executable engineering reconstruction, not a numerically verified conversion of the original Simulink model.** The supplied PDF has no `.slx`/`.mdl` file, complete block equations, parameter values, initial conditions, controller sampling details, or reference traces. Every numerical value in `config.json` is an illustrative assumption, not a manufacturer specification. Replace them with measured/datasheet values before evaluating a real drive. This model must not be used as hardware protection firmware.

## Run

Python 3.10+; validated here with Python 3.12.14, NumPy 2.3.5 and Matplotlib 3.10.8.

```bash
cd /workspace/Hackathon_ai
python -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
.venv/bin/python -m motor_sim --config config.json --output outputs
.venv/bin/python -m unittest discover -s tests -v
```

If the dependencies already exist in your active Python environment, use `python` directly. No MATLAB, GPU, external service, API key, or proprietary library is needed. NumPy handles state/vector calculations; Matplotlib generates headless PNG plots. RK4 is implemented explicitly so the controller's fixed sampling period and zero-order-held voltage are unambiguous.

## Configure parameters

Edit **config.json**; typed dataclasses and validation are in **motor_sim/config.py**. All angles are radians, time seconds, mechanical speed rad/s, inductance H, resistance ohms, current A, torque Nm, flux Wb, inertia kg m², and temperature °C. RPM = rad/s × 60/(2π). Pole pairs convert mechanical speed to electrical speed. Currents/voltages use an orthonormal two-phase dq convention, not a three-phase RMS convention; do not paste catalog RMS or line-to-line constants without conversion.

- `motor`: R, Ld, Lq, magnet flux, pole pairs, inertia, friction, copper temperature coefficient, cable resistance per metre and length, lumped thermal capacity/resistance. `cable_ohm_per_m` represents the entire configured phase wiring path.
- `control`: DC-link supply, reserve fraction, total current reference limit, negative d-current limit, PI gains, anti-windup gain, thermal thresholds, and independent controller/motor field-weakening flags. The model derives the torque constant from pole pairs × magnet flux for energy consistency; it does not independently specify incompatible torque/back-EMF constants.
- `simulation`: duration, integration/control step, output decimation, ambient/initial temperature, initial speed/position, load step, brake step, and command profile. Profile entries `[time_s, speed_rad_s]` are held until the next entry; they are not interpolated. For position control set `mode` to `position` and configure `target_position_rad` and `max_speed_rad_s`. The position loop is proportional around the PI speed and current loops; the PDF does not define an implementable position PI law.

The default scenario accelerates to 80 rad/s at 0.2 s, commands 260 rad/s at 0.8 s, adds 0.08 Nm load at 1.5 s, then commands 100 rad/s at 2.4 s. The brake starts outside the default simulation at 4 s; move `brake_start_s` inside the duration to exercise it. Positive configured load is a signed opposing torque for forward operation; reverse load conventions must be configured deliberately.

## Classes and equations

`Motor` owns electrical, thermal, and mechanical dynamics; `PIController` owns integrator state; `FieldWeakening` calculates negative d-current; `DriveController` coordinates limiting, derating and cascaded control; `Simulation` integrates states and records signals. `TwoPhaseTransform` implements Park/inverse Park for phase reconstruction. The state is `[id, iq, omega, theta, temperature]`.

For electrical speed `we = pole_pairs * omega`:

```
did/dt = (vd - R*id + we*Lq*iq) / Ld
diq/dt = (vq - R*iq - we*(Ld*id + psi_pm)) / Lq
torque = pole_pairs * (psi_pm + (Ld-Lq)*id) * iq
domega/dt = (torque - load - viscous*omega
             - (coulomb + brake)*tanh(omega/0.1)) / inertia
dtheta/dt = omega
R = R0 * (1 + alpha*(temperature-reference_temperature)) + cable_resistance
P_copper = R*(id² + iq²)
P_electrical = vd*id + vq*iq
P_mechanical = torque*omega
dtemperature/dt = (P_copper - (temperature-ambient)/thermal_resistance) / thermal_capacity
```

The lumped thermal node includes cable Joule loss for simplicity. Core, switching, inverter, and detailed spatial thermal losses are omitted. Electrical input also supplies stored magnetic energy; it is not equal to copper loss.

Field weakening solves the steady-state voltage ellipse analytically for the negative d-current nearest zero:

```
(R*id - we*Lq*iq)² + (R*iq + we*(Ld*id + psi_pm))² <= V_available²
V_available = dc_link/sqrt(2) * (1-reserve_fraction)
```

If no real feasible root exists, it selects the voltage-minimizing d-current and applies the configured limit. The root uses the speed request (or actual speed if its magnitude is higher) and preliminary q-current. q-current is subsequently reduced to satisfy `abs(id_ref)+abs(iq_ref) <= temperature_derated_current_limit`. This conservative sequence can leave additional voltage headroom; it is not an optimal torque-per-volt solver. Both field-weakening flags must be enabled. Current PI controllers include back-EMF/cross-coupling feedforward, a shared vector voltage clamp, and back-calculation anti-windup. The speed PI also tracks its saturated current command.

The 25% reserve is interpreted as a fraction of `dc_link/sqrt(2)`. The PDF mixes peak/RMS conventions and does not specify reserve ordering; this assumption is configurable. Negative id reduces **effective q-axis flux linkage** `psi_pm + Ld*id`; permanent magnet flux itself is constant. These choices resolve the PDF's inconsistent flux arrow and its mislabeling of electrical input as copper loss.

## Outputs and validation

- `outputs/comparison.png`: commanded/actual speed, d/q currents and references, torque, voltage and temperature, with field weakening configured versus disabled.
- `outputs/power_position.png`: electrical input, copper loss, mechanical power, and shaft position.
- `outputs/configured.csv`, `without_fw.csv`: sampled signals with units in column names.
- `outputs/parameters_used.json`: exact input configuration.
- `outputs/summary.json`: numerical metrics. The high-speed mean uses 2.1–2.4 s for the default experiment; customize that analysis window if changing the profile.

Seven automated tests cover Park power invariance and round-trip, instantaneous electromechanical energy balance, field-weakening voltage feasibility/enable flags, closed-loop tracking/reference limits/voltage saturation/time-step refinement, invalid configuration rejection, braking/thermal derating, and position tracking. Full-run step halving changes speed by at most 0.065 rad/s. No Simulink parity claim is possible without the original model and traces.

Default result: mean high-speed operation approximately **260 rad/s with field weakening versus 151.7 rad/s without**, with applied voltage limited to **25.456 V**. Negative id reaches about **−8 A** and increases heating.

**Reference current limiting is not a hard physical current clamp.** Actual current can overshoot under transients and voltage saturation: the default runs peak at approximately 10.02 A (configured) and 10.77 A (disabled) in the L1 norm despite a 10 A reference limit. The plant state is not artificially projected onto a current boundary; actual overcurrent protection, demagnetization maps and inverter shutdown dynamics require additional hardware-specific information. Temperature derating and controller saturation do not establish hardware safety.

No PWM switching, inverter dead-time, encoder quantization, magnetic saturation, detent/cogging torque, or vendor-specific stepper behavior is modeled. The default Ld=Lq is a nonsalient two-phase PM approximation, not an identified EMMT-ST motor. A switching model, vendor parameter identification, or exact Simulink trace matching is separate work requiring source data.
