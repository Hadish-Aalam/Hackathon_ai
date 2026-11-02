import argparse
from copy import deepcopy
import json
from pathlib import Path
import numpy as np
import os
os.environ.setdefault('MPLCONFIGDIR', '/tmp/electromechanical-matplotlib')
os.environ.setdefault('XDG_CACHE_HOME', '/tmp/electromechanical-cache')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from .config import Configuration
from .model import Simulation

def main():
    parser = argparse.ArgumentParser(description='Two-phase FOC electromechanical simulation')
    parser.add_argument('--config', default='config.json')
    parser.add_argument('--output', default='outputs')
    args = parser.parse_args()
    cfg = Configuration.load(args.config)
    out = Path(args.output); out.mkdir(parents=True, exist_ok=True)
    results = {}
    for name, enabled in [('configured', True), ('without_fw', False)]:
        case = deepcopy(cfg)
        if not enabled: case.control.field_weakening_controller = False
        data = Simulation(case).run()
        results[name] = data
        np.savetxt(out/f'{name}.csv', np.column_stack(list(data.values())),
                   delimiter=',', header=','.join(data), comments='')
    colors = {'configured':'tab:blue', 'without_fw':'tab:orange'}
    fig, axs = plt.subplots(3, 2, figsize=(13,10), sharex=True, constrained_layout=True)
    specifications = [('speed_rad_s','Speed (rad/s)'),('id_a','d-axis current (A)'),
                      ('iq_a','q-axis current (A)'),('torque_nm','Motor torque (Nm)'),
                      ('voltage','Voltage magnitude (V)'),('temperature_c','Temperature (°C)')]
    for ax, (key, label) in zip(axs.flat, specifications):
        for name,d in results.items():
            y = np.hypot(d['vd_v'],d['vq_v']) if key == 'voltage' else d[key]
            ax.plot(d['time_s'],y,label=name,color=colors[name])
        if key == 'speed_rad_s': ax.plot(data['time_s'],data['speed_ref_rad_s'],'k--',label='reference')
        if key in ('id_a','iq_a'):
            ax.plot(results['configured']['time_s'], results['configured'][key.replace('_a','_ref_a')], '--',color='tab:green',label='configured reference')
        if key == 'voltage': ax.axhline(cfg.control.voltage_limit,color='k',ls='--',label='available voltage')
        ax.set_ylabel(label); ax.grid(alpha=.3); ax.legend(fontsize=8)
    for ax in axs[-1]: ax.set_xlabel('Time (s)')
    fig.suptitle('Averaged two-phase FOC — illustrative parameters')
    fig.savefig(out/'comparison.png', dpi=160); plt.close(fig)
    fig,axs = plt.subplots(2,1,figsize=(11,7),sharex=True,constrained_layout=True)
    d = results['configured']
    for key in ['electrical_power_w','copper_loss_w','mechanical_power_w']:
        axs[0].plot(d['time_s'],d[key],label=key)
    axs[0].set_ylabel('Power (W)'); axs[0].legend()
    axs[1].plot(d['time_s'],d['position_rad'],label='shaft position')
    axs[1].set_ylabel('Position (rad)'); axs[1].set_xlabel('Time (s)'); axs[1].legend()
    for ax in axs: ax.grid(alpha=.3)
    fig.savefig(out/'power_position.png',dpi=160); plt.close(fig)
    summary = {}
    for name,d in results.items():
        mask = (d['time_s'] >= 2.1) & (d['time_s'] < 2.4)
        summary[name] = {'peak_voltage_v': float(np.hypot(d['vd_v'],d['vq_v']).max()),
                         'minimum_id_a':float(d['id_a'].min()),
                         'peak_actual_l1_current_a':float((abs(d['id_a'])+abs(d['iq_a'])).max()),
                         'high_speed_mean_rad_s':float(d['speed_rad_s'][mask].mean()) if mask.any() else None,
                         'final_temperature_c':float(d['temperature_c'][-1])}
    (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    (out/'parameters_used.json').write_text(Path(args.config).read_text())
    print(json.dumps(summary,indent=2))
    print(f'Plots and CSV data saved to {out.resolve()}')

if __name__ == '__main__':
    main()
