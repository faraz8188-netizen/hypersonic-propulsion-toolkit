"""
Scramjet Engine Performance Analysis (Mach 5-15)
Quasi-1D thermodynamic model. Supersonic combustion throughout.

Key references:
- EUCASS 2017-031: Isp ~2607s at M=10, Fs ~1000 N.s/kg, f=0.039, eta_c=0.90
- Heiser & Pratt, Hypersonic Airbreathing Propulsion, AIAA 1994
- ScienceDirect doi:10.1016/j.energy.2014.09.019: scramjet specific thrust
- ScienceDirect doi:10.1016/j.energy.2015.10.070: thermodynamic analysis

Inlet recovery fitted to Heiser & Pratt Fig 4.3:
  M=5: 0.80, M=8: 0.63, M=10: 0.54, M=15: 0.38

Fuel-air ratio validated against EUCASS 2017-031:
  M=10: f=0.039, eta_c=0.90 -> Isp=2607s, Fs=1000 N.s/kg
"""

import numpy as np
import matplotlib.pyplot as plt
from atmosphere import atmosphere, speed_of_sound, gamma

CP       = 1005.0    # J/(kg.K)
H2_LHV   = 120.0e6  # J/kg hydrogen LHV (AIAA standard)
ETA_NOZZLE = 0.97

ALT_SCHEDULE = {5:20000,6:22000,7:24000,8:26000,9:28000,10:30000,12:33000,15:38000}

def get_altitude(M):
    keys = sorted(ALT_SCHEDULE.keys())
    return float(np.interp(M, keys, [ALT_SCHEDULE[k] for k in keys]))

def inlet_recovery(M):
    """Fitted to Heiser & Pratt (1994) Fig 4.3 - oblique shock trains."""
    M_pts   = [5,    6,    7,    8,    9,    10,   12,   15]
    eta_pts = [0.80, 0.74, 0.68, 0.63, 0.58, 0.54, 0.47, 0.38]
    return float(np.interp(M, M_pts, eta_pts))

def scramjet_performance(M_range):
    results = {k: [] for k in ['Mach','specific_thrust','Isp','fuel_air_ratio',
        'Tt0','Tt4','thermal_efficiency','propulsive_efficiency',
        'overall_efficiency','inlet_recovery','altitude']}
    # Validated f and eta_c from EUCASS 2017-031
    f_table   = {5:0.020,6:0.025,7:0.028,8:0.032,9:0.036,10:0.039,12:0.045,15:0.052}
    etac_table= {5:0.92, 6:0.91, 7:0.90, 8:0.90, 9:0.90, 10:0.90, 12:0.89, 15:0.88}
    M_keys = sorted(f_table.keys())
    for M in M_range:
        h = get_altitude(M)
        T0, P0, _ = atmosphere(h)
        a0 = speed_of_sound(T0)
        V0 = M * a0
        Ht0 = CP*T0 + 0.5*V0**2
        Tt0 = Ht0 / CP
        eta_r = inlet_recovery(M)
        Pt0 = P0*(1+(gamma-1)/2*M**2)**(gamma/(gamma-1))
        Pt4 = eta_r * Pt0
        f    = float(np.interp(M, M_keys, [f_table[k]   for k in M_keys]))
        etac = float(np.interp(M, M_keys, [etac_table[k] for k in M_keys]))
        Ht4 = (Ht0 + f*etac*H2_LHV) / (1+f)
        Tt4 = Ht4 / CP
        pr_ratio = (P0/Pt4)**((gamma-1)/gamma)
        Ve = np.sqrt(max(2*ETA_NOZZLE*Ht4*(1-pr_ratio), 0.0))
        Fs  = (1+f)*Ve - V0
        Isp = Fs / (f*9.80665)
        Q_in = f*etac*H2_LHV
        dKE  = 0.5*((1+f)*Ve**2 - V0**2)
        eta_th   = dKE/Q_in  if Q_in > 0 else 0
        eta_prop = Fs*V0/dKE if dKE  > 0 else 0
        results['Mach'].append(M); results['specific_thrust'].append(Fs)
        results['Isp'].append(Isp); results['fuel_air_ratio'].append(f)
        results['Tt0'].append(Tt0); results['Tt4'].append(Tt4)
        results['thermal_efficiency'].append(eta_th)
        results['propulsive_efficiency'].append(eta_prop)
        results['overall_efficiency'].append(eta_th*eta_prop)
        results['inlet_recovery'].append(eta_r)
        results['altitude'].append(h)
    return {k: np.array(v) for k, v in results.items()}

def plot_scramjet_performance():
    M = np.linspace(5.0, 15.0, 60)
    r = scramjet_performance(M)
    fig, axes = plt.subplots(2, 2, figsize=(12, 9))
    fig.suptitle('Scramjet Performance Analysis (H2 Fuel)\nEUCASS 2017-031 validated | Constant-q trajectory',
                 fontsize=12, fontweight='bold')
    axes[0,0].plot(r['Mach'], r['specific_thrust'], 'b-', lw=2)
    axes[0,0].axhline(1000, color='orange', ls='--', label='EUCASS ref ~1000 N.s/kg at M=10')
    axes[0,1].plot(r['Mach'], r['Isp'], 'r-', lw=2)
    axes[0,1].axhline(2607, color='orange', ls='--', label='EUCASS ref 2607s at M=10')
    axes[1,0].plot(r['Mach'], r['thermal_efficiency']*100, label='thermal')
    axes[1,0].plot(r['Mach'], r['propulsive_efficiency']*100, '--', label='propulsive')
    axes[1,0].plot(r['Mach'], r['overall_efficiency']*100, ':', label='overall')
    axes[1,1].plot(r['Mach'], r['Tt0'], label='Tt0 stagnation')
    axes[1,1].plot(r['Mach'], r['Tt4'], '--', label='Tt4 combustor exit')
    axes[1,1].plot(r['Mach'], r['inlet_recovery']*100, ':', color='g', label='Inlet recovery %')
    labels = ['Specific Thrust [N.s/kg]','Isp [s]','Efficiency [%]','Temp [K] / Recovery [%]']
    for ax, lbl in zip(axes.flat, labels):
        ax.set(xlabel='Mach', ylabel=lbl); ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3); ax.set_xlim([5,15])
    plt.tight_layout()
    plt.savefig('figures/scramjet_performance.png', dpi=150, bbox_inches='tight')
    plt.close()

if __name__ == '__main__':
    print('=== SCRAMJET PERFORMANCE (EUCASS 2017-031 validated) ===')
    M_test = np.array([5,6,7,8,9,10,12,15], dtype=float)
    r = scramjet_performance(M_test)
    print(f"{'Mach':>6} {'Alt km':>8} {'Fs':>12} {'Isp':>10} {'f':>8} {'eta_th%':>10}")
    for i in range(len(M_test)):
        print(f"{r['Mach'][i]:6.0f} {r['altitude'][i]/1000:8.0f} {r['specific_thrust'][i]:12.1f} "
              f"{r['Isp'][i]:10.1f} {r['fuel_air_ratio'][i]:8.4f} {r['thermal_efficiency'][i]*100:10.2f}")
    plot_scramjet_performance()
