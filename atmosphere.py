"""
US Standard Atmosphere 1976 Model
Source: NASA TM-X-74335, US Standard Atmosphere 1976
Valid altitude range: 0 - 86 km
"""

import numpy as np

# Constants (SI units)
R_air = 287.05   # J/(kg.K) - specific gas constant for air
gamma = 1.4      # ratio of specific heats for air (calorically perfect)
g0    = 9.80665  # m/s^2 - standard gravity

# Layer definitions: (base altitude m, base temp K, lapse rate K/m)
_LAYERS = [
    (0,     288.15, -0.0065),   # Troposphere
    (11000, 216.65,  0.0),      # Tropopause
    (20000, 216.65,  0.001),    # Stratosphere 1
    (32000, 228.65,  0.0028),   # Stratosphere 2
    (47000, 270.65,  0.0),      # Stratopause
    (51000, 270.65, -0.0028),   # Mesosphere 1
    (71000, 214.65, -0.002),    # Mesosphere 2
]
P0   = 101325.0  # Pa - sea level pressure
rho0 = 1.225     # kg/m^3 - sea level density

def atmosphere(h):
    """Compute T(K), P(Pa), rho(kg/m^3) at altitude h (m)."""
    h = np.atleast_1d(np.float64(h))
    T = np.zeros_like(h); P = np.zeros_like(h); rho = np.zeros_like(h)
    for i, alt in enumerate(h):
        alt = float(np.clip(alt, 0, 86000))
        T_base, P_base, h_base, L = 288.15, P0, 0.0, -0.0065
        for j, (h_bot, T_bot, lapse) in enumerate(_LAYERS):
            h_top = _LAYERS[j+1][0] if j+1 < len(_LAYERS) else 86000
            if j > 0:
                prev_h, prev_T, prev_L = _LAYERS[j-1]
                dh = h_bot - prev_h
                P_base = (P_base * np.exp(-g0*dh/(R_air*T_base)) if prev_L==0
                          else P_base*(T_base/(T_base+prev_L*dh))**(g0/(R_air*prev_L)))
                T_base = T_bot; h_base = h_bot; L = lapse
            if alt <= h_top: break
        dh = alt - h_base; T_i = T_base + L*dh
        P_i = (P_base*np.exp(-g0*dh/(R_air*T_base)) if L==0
               else P_base*(T_i/T_base)**(-g0/(R_air*L)))
        T[i], P[i], rho[i] = T_i, P_i, P_i/(R_air*T_i)
    return T.squeeze(), P.squeeze(), rho.squeeze()

def speed_of_sound(T):
    """Speed of sound (m/s) given temperature T (K)."""
    return np.sqrt(gamma * R_air * T)

def flight_velocity(M, T):
    """Flight velocity (m/s) from Mach number and temperature."""
    return M * speed_of_sound(T)
