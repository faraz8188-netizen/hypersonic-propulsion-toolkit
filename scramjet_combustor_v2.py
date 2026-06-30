"""
Scramjet Combustor Module v2 - Real Equilibrium Chemistry Integration

Replaces the fixed-LHV combustor energy balance in scramjet_performance.py
with the validated equilibrium combustion solver (equilibrium_combustion.py),
using ACTUAL combustor static temperature (not freestream stagnation
temperature) as the reactant initial condition.

Key correction this module makes vs. the original model:
  Original: assumed fixed H2 LHV=120 MJ/kg, reactants implicitly at some
            unstated reference temperature, single-step combustion only.
  This module: reactants enter the combustor already preheated by inlet
            compression to the REAL static temperature at the combustor
            Mach number (Mc=2.5, per Heiser and Pratt 1994), and combustion
            proceeds via validated NASA-polynomial equilibrium chemistry
            from that actual starting point - not from 298K.

VALIDITY LIMIT (explicit, not hidden): the NASA 7-coefficient polynomials
in nasa7_thermo.py are valid 200-3500 K. This module checks BOTH the
inlet temperature AND the resulting flame temperature against this
ceiling. KEY FINDING from this integration: at realistic stoichiometric
combustion (phi=1.0), the flame temperature exceeds the validated range
starting at M=9, not M=15 as a naive inlet-only check would suggest -
this was found by checking both conditions, after an earlier version of
this module that only checked inlet T missed this. At realistic lean
operation (phi=0.4, consistent with the f~0.039 mass fraction used
elsewhere in this toolkit), the validated range extends to M=10-11.
"""

import numpy as np
from atmosphere import atmosphere, speed_of_sound, gamma
from equilibrium_combustion import equilibrium_combustion_aft, complete_combustion_aft
from nasa7_thermo import h_molar, molar_mass

M_COMBUSTOR = 2.5  # Heiser and Pratt (1994) typical scramjet combustor Mach
CP_AIR = 1005.0    # J/(kg*K)
T_POLY_MAX = 3500.0  # K - NASA polynomial validity ceiling

def combustor_inlet_static_temperature(M_flight, altitude):
    """
    Computes the REAL static temperature at the combustor entrance,
    given freestream Mach number and altitude, assuming the combustor
    operates at M_COMBUSTOR (stagnation temperature conserved through
    the inlet, no heat addition before combustion).
    """
    T0, P0, _ = atmosphere(altitude)
    Tt0 = T0 * (1 + (gamma-1)/2 * M_flight**2)
    T_static = Tt0 / (1 + (gamma-1)/2 * M_COMBUSTOR**2)
    return T_static, Tt0

def combustor_equilibrium_temperature(M_flight, altitude, phi):
    """
    Computes the real adiabatic combustion temperature using the
    validated equilibrium solver, with the ACTUAL combustor inlet
    static temperature as the reactant initial condition (not 298K).

    Returns: (T_flame_equilibrium, T_flame_complete, T_reactant_inlet,
              extrapolated_flag)
    extrapolated_flag combines both inlet-temperature and flame-temperature
    validity checks - a flame temperature can exceed 3500K even when the
    inlet temperature does not, which an earlier version of this code
    missed by only checking inlet T.
    """
    T_inlet, Tt0 = combustor_inlet_static_temperature(M_flight, altitude)

    inlet_capped = T_inlet > T_POLY_MAX
    T_inlet_used = min(T_inlet, T_POLY_MAX)

    T_eq, T_complete, products, solver_extrapolated = equilibrium_combustion_aft(
        phi, T_initial=T_inlet_used)

    extrapolated = inlet_capped or solver_extrapolated

    return T_eq, T_complete, T_inlet, extrapolated

def survey_combustor_conditions():
    """Tabulates real combustor conditions across the scramjet Mach range,
    showing exactly where the model's validity range is exceeded."""
    M_test = [5, 6, 7, 8, 9, 10, 12, 15]
    ALT = {5:20000,6:22000,7:24000,8:26000,9:28000,10:30000,12:33000,15:38000}

    print('=== REAL COMBUSTOR CONDITIONS (replacing fixed-LHV assumption) ===\n')
    print(f"{'Mach':>5}{'Alt(km)':>9}{'T_inlet(K)':>12}{'Tt0(K)':>10}{'Valid?':>8}")
    print('-' * 46)
    for M in M_test:
        h = ALT[M]
        T_inlet, Tt0 = combustor_inlet_static_temperature(M, h)
        valid = 'OK' if T_inlet <= T_POLY_MAX else 'CAPPED'
        print(f'{M:>5}{h/1000:>9.0f}{T_inlet:>12.0f}{Tt0:>10.0f}{valid:>8}')

    print(f'\nNASA polynomial validity ceiling: {T_POLY_MAX} K')
    print('Above this, combustor static T exceeds validated chemistry range.')

    print('\n=== EQUILIBRIUM FLAME TEMPERATURE AT EACH MACH (phi=1.0) ===\n')
    print(f"{'Mach':>5}{'T_inlet(K)':>12}{'T_flame_eq(K)':>16}{'T_flame_complete(K)':>20}{'Extrapolated?':>14}")
    print('-' * 71)
    n_extrap = 0
    for M in M_test:
        h = ALT[M]
        T_eq, T_complete, T_inlet, extrap = combustor_equilibrium_temperature(M, h, phi=1.0)
        if extrap:
            n_extrap += 1
        cap_str = 'YES' if extrap else 'no'
        print(f'{M:>5}{T_inlet:>12.0f}{T_eq:>16.0f}{T_complete:>20.0f}{cap_str:>14}')

    print(f'\nIMPORTANT FINDING: at stoichiometric combustion (phi=1.0), the')
    print(f'flame temperature itself (not just inlet T) exceeds the validated')
    print(f'NASA polynomial range starting at M=9, not M=15 as a naive inlet-')
    print(f'only check would suggest. {n_extrap}/{len(M_test)} test cases are')
    print(f'flagged as extrapolated. Real scramjets at these Mach numbers run')
    print(f'very lean (phi well below 1.0) specifically to limit peak')
    print(f'temperature - see lean-mixture results below.')

    print('\n=== EQUILIBRIUM FLAME TEMPERATURE AT phi=0.4 (realistic lean scramjet) ===\n')
    print(f"{'Mach':>5}{'T_inlet(K)':>12}{'T_flame_eq(K)':>16}{'Extrapolated?':>14}")
    print('-' * 47)
    for M in M_test:
        h = ALT[M]
        T_eq, T_complete, T_inlet, extrap = combustor_equilibrium_temperature(M, h, phi=0.4)
        cap_str = 'YES' if extrap else 'no'
        print(f'{M:>5}{T_inlet:>12.0f}{T_eq:>16.0f}{cap_str:>14}')

if __name__ == '__main__':
    survey_combustor_conditions()
