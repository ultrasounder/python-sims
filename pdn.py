import numpy as np
import matplotlib.pyplot as plt

# --- 1. Define Physics Constants ---

freq = np.logspace(3, 7, 1000) # 1Khz steps of 10MHz
omega = 2 * np.pi * freq

# 2. define Components 
# VRM Model (Sinmplified as Inductor + Resistance)
L_vrm = 150e-9 #Typical for a switching regulator
R_vrm = 1e-3 # 1 mOhm

# Bad bulk capacitor 
C_bulk_bad = 47e-6
ESR_bulk = 5e-3 # 5 mOhm
ESL_bulk = 2e-9 # 2nH

#Good Bulk Capacitor (The "fix": 3x 10uF + 2x 1uF)
# For simplicity in this script, we lump them, but strictly you wil parallel them all
c_bulk_good = (3 * 10e-6) + (2 * 1e-6) # All of these should be discrete values in parallel
ESR_good = 2e-3 # Low ESR due to parallel caps
ESL_good = 1e-9 # Lower ESL

# 3. Calculate Impedance (z)
# z = R + jwL + 1/jwC

def get_impedance_cap(C, ESR, ESL, w):
    return ESR + 1j * (w * ESL - 1/(w * C))

def get_impedance_vrm(R, L, w):
    return R + 1j * w * L

#calculate individual impedances
Z_vrm_vals = get_impedance_vrm(R_vrm, L_vrm, omega)

Z_cap_bad_vals = get_impedance_cap(C_bulk_bad, ESL_bulk, ESL_bulk, omega)

Z_cap_good_vals = get_impedance_cap(c_bulk_good, ESR_good, ESL_good, omega)

# Calculate PDN Impedance (VRM || CAP)
# Parallel formula Z_total = (z1 * z2) / (z1 + z2)
Z_pdn_bad = (Z_vrm_vals * Z_cap_bad_vals) / (Z_vrm_vals + Z_cap_bad_vals)
Z_pdn_good = (Z_vrm_vals * Z_cap_good_vals) / (Z_vrm_vals + Z_cap_good_vals)

# 4. Target Impedance Logic 
# Based on a 0.9V PCIE Retimer power supply spec 2 A transient, 3% Tolerance
V_dd = 0.9
Allowed_Ripple = 0.03 # 3%
I_transient = 2.0
Z_target = (V_dd * Allowed_Ripple) / I_transient

# 5. Visualization
plt.figure(figsize=(10, 6))
plt.loglog(freq, np.abs(Z_pdn_bad), 'r', linewidth=2, label='Rev A (Fail): 1x 47uF')
plt.loglog(freq, np.abs(Z_pdn_good), 'g', linewidth=2, label='Rev B (Fix): 3x10uF + 2x1uF')
plt.axhline(y=Z_target, color='b', linestyle='--', label=f'Target Z ({Z_target*1000:.1f} mOhm)')
plt.title("PDN Impedance: The 'Rogue Wave' at 320 kHz")
plt.xlabel("Frequency (Hz)")
plt.ylabel("Impedance (Ohms)")
plt.grid(True, which="both", ls="-")
plt.legend()
plt.show()
