"""
Ramjet Engine Performance Analysis
Thermodynamic model based on Brayton cycle with real inlet losses.

Key references:
- Heiser & Pratt, Hypersonic Airbreathing Propulsion, AIAA, 1994
- MDPI Energies 11(10):2558, 2018 - Ramjet Compression System
- DTIC ADA596056 - Ramjet Intakes (NATO AVT-185)
- NASA NTRS 19670095387 - Ram-jet Performance

Validated constants:
- Inlet total pressure recovery: 0.90-0.94 (NATO AVT-185 typical = 0.92)
- Combustion efficiency: 0.90-0.95 (published range)
- Hydrogen LHV: 120.0 MJ/kg | Kerosene LHV: 43.2 MJ/kg
- Ramjet operational range: Mach 2-5
- T4 material limit: 2500 K
"""

import numpy as np
import matplotlib.pyplot as plt
from atmosphere import atmosphere, speed_of_sound, flight_velocity, gamma, R_air

FUELS = {
    'hydrogen': {'LHV': 120.0e6, 'label': 'Hydrogen (H2)'},
    'kerosene': {'LHV':  43.2e6, 'label': 'Kerosene (JP-7)'},
}

ETA_INLET  = 0.92    # NATO AVT-185 typical
ETA_COMB   = 0.92    # published range 0.90-0.95
ETA_NOZZLE = 0.97
T4_MAX     = 2500.0  # K material limit
CP_AIR     = 1005.0  # J/(kg.K)
ALTITUDE   = 20000   # m

def inlet_total_conditions(M, T0, P0):
    """Isentropic total conditions with inlet pressure recovery."""
    Tt0 = T0 * (1 + (gamma-1)/2 * M**2)
    Pt0 = P0 * (1 + (gamma-1)/2 * M**2)**(gamma/(gamma-1))
    Pt2 = ETA_INLET * Pt0
    return Tt0, Pt2

def ramjet_performance(M_range, altitude=ALTITUDE, fuel='hydrogen'):
    T0, P0, _ = atmosphere(altitude)
    a0  = speed_of_sound(T0)
    LHV = FUELS[fuel]['LHV']
    results = {k: [] for k in ['Mach','specific_thrust','Isp','Tt2','Tt4',
                                'fuel_air_ratio','thermal_efficiency',
                                'propulsive_efficiency','overall_efficiency']}
    for M in M_range:
        V0 = M * a0
        Tt2, Pt2 = inlet_total_conditions(M, T0, P0)
        f_max = CP_AIR * (T4_MAX - Tt2) / (ETA_COMB * LHV - CP_AIR * T4_MAX)
        f = max(0.001, min(f_max, 0.10))
        Tt4 = min((CP_AIR*Tt2 + f*ETA_COMB*LHV)/(CP_AIR*(1+f)), T4_MAX)
        pr  = Pt2 / P0
        Ve  = ETA_NOZZLE * np.sqrt(2*CP_AIR*Tt4*(1-(1/pr)**((gamma-1)/gamma)))
        Fs  = (1+f)*Ve - V0
        Isp = Fs / (f * 9.80665)
        Q_in   = f * ETA_COMB * LHV
        KE_net = 0.5*((1+f)*Ve**2 - V0**2)
        eta_th   = KE_net / Q_in if Q_in > 0 else 0
        eta_prop = Fs*V0 / KE_net if KE_net > 0 else 0
        results['Mach'].append(M); results['specific_thrust'].append(Fs)
        results['Isp'].append(Isp); results['Tt2'].append(Tt2)
        results['Tt4'].append(Tt4); results['fuel_air_ratio'].append(f)
        results['thermal_efficiency'].append(eta_th)
        results['propulsive_efficiency'].append(eta_prop)
        results['overall_efficiency'].append(eta_th*eta_prop)
    return {k: np.array(v) for k, v in results.items()}

def plot_ramjet_performance():
    M = np.linspace(2.0, 5.0, 50)
    fig, axes = plt.subplots(2, 2, figsize=(12, 9))
    fig.suptitle('Ramjet Engine Performance Analysis\nAltitude: 20 km | Inlet η=0.92 | Combustion η=0.92',
                 fontsize=13, fontweight='bold')
    for fuel in ['hydrogen', 'kerosene']:
        r = ramjet_performance(M, fuel=fuel)
        lbl = FUELS[fuel]['label']
        axes[0,0].plot(r['Mach'], r['specific_thrust'], label=lbl)
        axes[0,1].plot(r['Mach'], r['Isp'], label=lbl)
        axes[1,0].plot(r['Mach'], r['thermal_efficiency']*100, label=f'{lbl} th')
        axes[1,0].plot(r['Mach'], r['propulsive_efficiency']*100, '--', label=f'{lbl} prop')
        axes[1,1].plot(r['Mach'], r['Tt2'], label=f'{lbl} Tt2')
        axes[1,1].plot(r['Mach'], r['Tt4'], '--', label=f'{lbl} Tt4')
    titles = ['Specific Thrust [N.s/kg]','Specific Impulse [s]',
              'Efficiency [%]','Temperature [K]']
    for ax, t in zip(axes.flat, titles):
        ax.set(xlabel='Mach', ylabel=t, title=t+' vs Mach')
        ax.legend(fontsize=7); ax.grid(True, alpha=0.3); ax.set_xlim([2,5])
    axes[1,1].axhline(T4_MAX, color='red', ls=':', label=f'T4 limit {T4_MAX}K')
    plt.tight_layout()
    plt.savefig('figures/ramjet_performance.png', dpi=150, bbox_inches='tight')
    plt.close()

if __name__ == '__main__':
    print('=== RAMJET PERFORMANCE ===')
    M_test = np.array([2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0])
    r = ramjet_performance(M_test, fuel='hydrogen')
    print(f"{'Mach':>6} {'Fs':>12} {'Isp':>10} {'eta_th%':>10}")
    for i in range(len(M_test)):
        print(f"{r['Mach'][i]:6.1f} {r['specific_thrust'][i]:12.1f} {r['Isp'][i]:10.1f} {r['thermal_efficiency'][i]*100:10.2f}")
    plot_ramjet_performance()
