import numpy as np
import matplotlib.pyplot as plt

freq = np.logspace(3, 10, 2000)  # 1 kHz to 10 GHz
omega = 2 * np.pi * freq

# Bulk capacitor
C_bulk = 470e-6
ESR_bulk = 50e-3
L_bulk = 5e-9
Z_bulk = ESR_bulk + 1j * (omega * L_bulk - 1 / (omega * C_bulk))

# MLCC
C_mlcc = 10e-6
ESR_mlcc = 5e-3
L_mlcc = 0.7e-9
Z_mlcc = ESR_mlcc + 1j * (omega * L_mlcc - 1 / (omega * C_mlcc))

# Parallel combination
Z_para = 1 / (1/Z_bulk + 1/Z_mlcc)

# Plot
fig, axes = plt.subplots(2, 2, figsize=(14, 10))

# Plot 1: Individual impedances
ax = axes[0, 0]
ax.loglog(freq, np.abs(Z_bulk)*1e3, 'b-', linewidth=2, label='470µF Bulk')
ax.loglog(freq, np.abs(Z_mlcc)*1e3, 'r-', linewidth=2, label='10µF MLCC')
ax.loglog(freq, np.abs(Z_para)*1e3, 'k--', linewidth=3, label='Parallel')
ax.grid(True, which='both', alpha=0.3)
ax.set_ylabel('|Z| (mΩ)')
ax.set_title('PDN Impedance: Two Cap Banks')
ax.legend()
ax.set_ylim([1, 5000])

# Plot 2: Real part (resistance)
ax = axes[0, 1]
ax.loglog(freq, np.real(Z_bulk)*1e3, 'b-', linewidth=2, label='470µF Bulk')
ax.loglog(freq, np.real(Z_mlcc)*1e3, 'r-', linewidth=2, label='10µF MLCC')
ax.loglog(freq, np.real(Z_para)*1e3, 'k--', linewidth=3, label='Parallel')
ax.grid(True, which='both', alpha=0.3)
ax.set_ylabel('Re(Z) (mΩ)')
ax.set_title('Resistance Component')
ax.legend()

# Plot 3: Imaginary part (reactance)
ax = axes[1, 0]
ax.semilogx(freq, np.imag(Z_bulk)*1e3, 'b-', linewidth=2, label='470µF Bulk')
ax.semilogx(freq, np.imag(Z_mlcc)*1e3, 'r-', linewidth=2, label='10µF MLCC')
ax.semilogx(freq, np.imag(Z_para)*1e3, 'k--', linewidth=3, label='Parallel')
ax.axhline(0, color='k', linestyle='-', alpha=0.2)
ax.grid(True, which='both', alpha=0.3)
ax.set_ylabel('Im(Z) (mΩ)')
ax.set_title('Reactance (Inductive = +, Capacitive = −)')
ax.legend()

# Plot 4: Highlight anti-resonance
ax = axes[1, 1]
ax.loglog(freq, np.abs(Z_para)*1e3, 'k-', linewidth=2.5, label='Parallel Z')
ax.axhline(10, color='g', linestyle='--', label='Target Z_target = 10 mΩ')

# Find anti-resonance peak (between the two SRFs)
f_srf_bulk = 1 / (2*np.pi*np.sqrt(L_bulk*C_bulk))
f_srf_mlcc = 1 / (2*np.pi*np.sqrt(L_mlcc*C_mlcc))
idx_areason = np.argmax(np.abs(Z_para)) 
f_areason = freq[idx_areason]
z_areason = np.abs(Z_para)[idx_areason]

ax.plot(f_areason, z_areason*1e3, 'ro', markersize=10, label=f'Anti-resonance\n@ {f_areason/1e6:.1f} MHz')
ax.axvline(f_srf_bulk, color='b', linestyle=':', alpha=0.5, label=f'Bulk SRF = {f_srf_bulk/1e6:.2f} MHz')
ax.axvline(f_srf_mlcc, color='r', linestyle=':', alpha=0.5, label=f'MLCC SRF = {f_srf_mlcc/1e6:.2f} MHz')

ax.grid(True, which='both', alpha=0.3)
ax.set_ylabel('|Z| (mΩ)')
ax.set_xlabel('Frequency (Hz)')
ax.set_title('Anti-Resonance: THE FAILURE POINT')
ax.legend(fontsize=9)
ax.set_ylim([1, 200])

plt.tight_layout()
plt.savefig('anti_resonance_explanation.png', dpi=150, bbox_inches='tight')
plt.show()

print(f"\nBulk SRF: {f_srf_bulk/1e6:.2f} MHz")
print(f"MLCC SRF: {f_srf_mlcc/1e6:.2f} MHz")
print(f"Anti-resonance peak @ {f_areason/1e6:.2f} MHz: {z_areason*1e3:.2f} mΩ")
print(f"Peak is {z_areason/10:.1f}x higher than target impedance!")


# ## **Why Anti-Resonance Happens**

# Between the two SRFs (3.3 MHz and 19 MHz), something pathological occurs:

# 1. **Below bulk SRF (< 3.3 MHz):**
#    - Bulk cap dominates (capacitive)
#    - MLCC helps a bit
#    - Parallel impedance is low (good)

# 2. **Right at bulk SRF (≈ 3.3 MHz):**
#    - Bulk impedance = 50 mΩ (its minimum)
#    - MLCC impedance ≈ 100 mΩ (still somewhat inductive)
#    - They parallel to ≈ 33 mΩ

# 3. **Between bulk SRF and MLCC SRF (3.3 to 19 MHz) ← THE PROBLEM:**
#    - Bulk is now **INDUCTIVE** (looks like a small inductor)
#    - MLCC is still **CAPACITIVE** (still helping)
#    - Series inductance + parallel capacitance = **LC resonance tank**
#    - The parallel combination **PEAKS** at the anti-resonant frequency
#    - This peak can exceed **100 mΩ** or even **1 Ω** in bad designs!

# 4. **At MLCC SRF (≈ 19 MHz):**
#    - MLCC impedance = 5 mΩ (its minimum)
#    - Bulk is now very inductive ≈ 3 Ω
#    - They parallel to ≈ 5 mΩ (good again)

# 5. **Above MLCC SRF (> 100 MHz):**
#    - Both are inductive
#    - Parallel impedance ≈ min(inductive reactances)
#    - Begins to rise with frequency again

# ---

# ## **The Physics of Anti-Resonance**

# Think of it like a tuned LC circuit:
# ```
#     ┌─────────┐     ┌─────────┐
#     │ Bulk    │     │ MLCC    │
#     │ (L1+C1) │  || │ (L2+C2) │
#     └─────────┘     └─────────┘
    
# At some frequency between SRF1 and SRF2:
#   - L1 (inductive) of the bulk cap and 
#   - C2 (capacitive) of the MLCC 
#   form an LC tank circuit
  
# This tank "rings" at its own resonant frequency
# |Z_parallel| = max at the anti-resonance frequency
# ```

# The **Q factor** (sharpness of the peak) depends on ESR:
# - **High ESR:** Broad, gentle peak (bad cap, but easier to handle)
# - **Low ESR:** Sharp, needle-like peak (good cap, but deadly anti-resonance)

# Modern ultra-low-ESR MLCCs have ESR ≈ 1 mΩ, so anti-resonance peaks can reach **10−100 mΩ**, easily exceeding a 10 mΩ target impedance.

# ---

# ## **Why This Kills Signal Integrity**

# When SSO (simultaneous switching output) transient current at frequency f_anti hits this peak:
# ```
# V_noise = Z_PDN(f_anti) × I_transient

# If f_anti = 8 MHz and Z_peak = 100 mΩ, I_SSO = 10 A:
#   V_noise = 0.1 Ω × 10 A = 1.0 V!

# For a 1.0 V supply with 3% ripple tolerance (30 mV budget):
#   1.0 V ripple CONSUMES THE ENTIRE BUDGET and causes supply collapse!