import numpy as np
import matplotlib.pyplot as plt

#10 uF MLCC parameters
C = 10e-6 # Farads
ESR = 5e-3 # Ohms
L = 0.7e-9 # Henries

freq = np.logspace(3, 10, 1000) # 1 KHz to 10 GHz
omega = 2 * np.pi * freq

# Impedance components
Z_R = ESR * np.ones_like(freq)
Z_L = 1j * omega * L
Z_C = -1j / (omega * C)

# Total impedance
Z_total = Z_R + Z_L + Z_C
Z_mag = np.abs(Z_total)

# self resonant frequency (SRF): where inductive and capacitive reactance cancel each other
f_SRF = 1 / (2 * np.pi * np.sqrt(L * C))
print(f"SRF for 10uF cap with 0.7nH ESL: {f_SRF/1e6: .2f} MHz")

# Plot
plt.figure(figsize=(12, 6))
plt.loglog(freq, Z_mag * 1e3, 'b-', linewidth=2.5, label='|Z_cap|')
plt.axhline(ESR * 1e3, color='r', linestyle='--', label=f'ESR = {ESR*1e3:.1f} mΩ')
plt.axvline(f_SRF, color='g', linestyle='--', label=f'SRF = {f_SRF/1e6:.1f} MHz')
plt.grid(True, which='both', alpha=0.3)
plt.xlabel('Frequency (Hz)')
plt.ylabel('|Z| (mΩ)')
plt.title('10µF Ceramic MLCC Impedance')
plt.legend()
plt.ylim([0.5, 500])
plt.show()