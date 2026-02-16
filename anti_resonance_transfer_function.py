"""
Anti-Resonance Transfer Function: The LC Tank Circuit Model
===========================================================

When two capacitor banks parallel, they form an LC resonance tank at the
frequency between their SRFs. This notebook derives the transfer function,
plots it, and shows how it manifests in PDN impedance.

KEY INSIGHT:
  A parallel LC tank has TWO distinct frequency behaviors:
  1. Series resonance (both elements resonate in phase) → minimum impedance
  2. Parallel/anti-resonance (elements oscillate out of phase) → MAXIMUM impedance
  
The anti-resonance is what kills PDN performance.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
import warnings
warnings.filterwarnings('ignore')

print("="*80)
print("ANTI-RESONANCE TRANSFER FUNCTION DERIVATION")
print("="*80)

# ==============================================================================
# PART 1: THE LUMPED-ELEMENT MODEL
# ==============================================================================

print("\n[STEP 1] Understand the lumped-element circuit model")
print("-" * 80)

print("""
PHYSICAL PICTURE:

When you parallel two capacitor banks, you have:

        ┌──────────────┬──────────────┐
        │              │              │
    ────┤              │              ├────  Power plane
        │              │              │
        │    Bank 1    │    Bank 2    │
        │  (Bulk cap)  │   (MLCC)     │
        │              │              │
    ────┤              │              ├────  Ground plane
        │              │              │
        └──────────────┴──────────────┘

Each bank is modeled as:
  
  Bank 1 (Bulk):        Bank 2 (MLCC):
  ┌─ R1 ─┬─ L1 ─┐       ┌─ R2 ─┬─ L2 ─┐
  │      │      │       │      │      │
  ├──────┴──────┤       ├──────┴──────┤
  │     C1      │       │     C2      │
  └─────────────┘       └─────────────┘

When parallelized, the total impedance is:

  Z_total(s) = 1 / (Y1(s) + Y2(s))
  
where admittance Y(s) = 1/Z(s).

For each bank:
  Z(s) = R + sL + 1/(sC)
  
Rewrite with standard form (transfer function):
  Z(s) = (Rs·L·C + s²LC + 1) / (sC)
       = (s²LC + sRC + 1) / (sC)

In the frequency domain (s = jω):
  Z(jω) = R + jωL + 1/(jωC)
        = R + j(ωL - 1/(ωC))
""")

# ==============================================================================
# PART 2: DERIVE THE IMPEDANCE OF A SINGLE RLC SERIES NETWORK
# ==============================================================================

print("\n[STEP 2] Single RLC series impedance (single capacitor bank)")
print("-" * 80)

print("""
For a single RLC series circuit:

  Z(s) = R + sL + 1/(sC)

Rewrite with common denominator:

  Z(s) = (Rs·L·C + s²LC + 1) / (sC)
       = (s²LC + sRC + 1) / (sC)

Divide numerator and denominator by LC:

  Z(s) = (s² + s(R/L) + 1/(LC)) / (sC)

Define:
  - ω₀ = 1/√(LC)     ← Natural frequency (SRF)
  - ζ = (R/2) * √(C/L)  ← Damping ratio (ESR effect)

Then:

  Z(s) = (s² + 2ζω₀·s + ω₀²) / (sC)
       = [1/(sC)] × (s² + 2ζω₀·s + ω₀²)

At the natural frequency ω₀ = 1/√(LC):
  - Impedance magnitude = ESR = R (minimum point)
  - Phase = 0° (resistive)
  
Below ω₀: capacitive (Im(Z) < 0, phase ≈ −90°)
Above ω₀: inductive (Im(Z) > 0, phase ≈ +90°)

Now, let's compute this numerically.
""")

# ==============================================================================
# PART 3: NUMERICAL EXAMPLE - SINGLE CAPACITOR BANK
# ==============================================================================

print("\n[STEP 3] Compute single bank impedance (10µF MLCC)")
print("-" * 80)

# 10µF MLCC parameters
C = 10e-6   # 10 µF
ESR = 5e-3  # 5 mΩ
L = 0.7e-9  # 0.7 nH

# Frequency range
freq = np.logspace(3, 10, 2000)  # 1 kHz to 10 GHz
omega = 2 * np.pi * freq
s = 1j * omega

# Impedance: Z(s) = R + sL + 1/(sC)
Z_single = ESR + s * L + 1 / (s * C)

# Natural frequency (SRF)
omega_0 = 1 / np.sqrt(L * C)
f_0 = omega_0 / (2 * np.pi)

print(f"10µF MLCC Parameters:")
print(f"  C = {C*1e6:.1f} µF")
print(f"  ESR = {ESR*1e3:.2f} mΩ")
print(f"  L = {L*1e9:.2f} nH")
print(f"  SRF (f₀) = 1/(2π√LC) = {f_0/1e6:.2f} MHz")
print(f"  |Z| at SRF = ESR = {ESR*1e3:.2f} mΩ")

# Extract magnitude, phase, real, imaginary
Z_mag = np.abs(Z_single)
Z_phase = np.angle(Z_single, deg=True)
Z_real = np.real(Z_single)
Z_imag = np.imag(Z_single)

# Find the minimum impedance point
idx_min = np.argmin(Z_mag)
f_min = freq[idx_min]
Z_min = Z_mag[idx_min]

print(f"  Measured minimum |Z| = {Z_min*1e3:.3f} mΩ at f = {f_min/1e6:.2f} MHz")

# ==============================================================================
# PART 4: TWO-BANK SYSTEM - THE ANTI-RESONANCE
# ==============================================================================

print("\n[STEP 4] Two-bank system: Bulk + MLCC → Anti-Resonance")
print("-" * 80)

# Bulk capacitor (470µF)
C_bulk = 470e-6
ESR_bulk = 50e-3
L_bulk = 5e-9

# MLCC capacitor (10µF)
C_mlcc = 10e-6
ESR_mlcc = 5e-3
L_mlcc = 0.7e-9

# Impedances
Z_bulk = ESR_bulk + s * L_bulk + 1 / (s * C_bulk)
Z_mlcc = ESR_mlcc + s * L_mlcc + 1 / (s * C_mlcc)

# Parallel combination: Z_total = 1 / (1/Z_bulk + 1/Z_mlcc)
Z_parallel = 1 / (1/Z_bulk + 1/Z_mlcc)

# SRFs
f_srf_bulk = 1 / (2 * np.pi * np.sqrt(L_bulk * C_bulk))
f_srf_mlcc = 1 / (2 * np.pi * np.sqrt(L_mlcc * C_mlcc))

print(f"Bulk capacitor (470µF):")
print(f"  SRF = {f_srf_bulk/1e6:.2f} MHz")
print(f"Ceramic MLCC (10µF):")
print(f"  SRF = {f_srf_mlcc/1e6:.2f} MHz")

# Extract parallel impedance properties
Z_para_mag = np.abs(Z_parallel)
Z_para_phase = np.angle(Z_parallel, deg=True)
Z_para_real = np.real(Z_parallel)
Z_para_imag = np.imag(Z_parallel)

# Find anti-resonance peak
idx_areson = np.argmax(Z_para_mag)
f_areson = freq[idx_areson]
Z_areson = Z_para_mag[idx_areson]

print(f"\nParallel (two-bank) system:")
print(f"  Anti-resonance frequency f_ar = {f_areson/1e6:.2f} MHz")
print(f"  Anti-resonance impedance |Z_ar| = {Z_areson*1e3:.2f} mΩ")
print(f"  Ratio to target (10 mΩ): {Z_areson*1e3 / 10:.1f}×")

# ==============================================================================
# PART 5: TRANSFER FUNCTION IN POLE-RESIDUE FORM
# ==============================================================================

print("\n[STEP 5] Express as transfer function H(s) = Z(s) (impedance function)")
print("-" * 80)

print("""
The impedance of a parallel LC tank can be expressed as:

  Z(s) = [s²(L₁||L₂) + s(R₁||R₂) + 1/(C₁||C₂)] / [...]

For our two-bank system, the denominator becomes zero (pole) at two frequencies:

SERIES RESONANCE (minimum impedance):
  When both capacitors "help" together → Z = R_min ≈ ESR
  This happens at each capacitor's individual SRF

PARALLEL/ANTI-RESONANCE (maximum impedance):
  When one is capacitive and one is inductive, they oppose each other
  Z becomes VERY LARGE (pole in the denominator)
  
The anti-resonance frequency is approximately:

  f_ar ≈ √(f_srf_bulk × f_srf_mlcc)    [geometric mean]

But this is approximate. The exact formula is complex.

Let me compute it directly using the pole-residue decomposition of the
parallel impedance function.
""")

# For parallel impedance, find poles and zeros using transfer function analysis
# Z_parallel = 1 / (Y_bulk + Y_mlcc) where Y = 1/Z

# The admittance of each bank in pole-residue form:
# Y(s) = (sC) / (s² + s(R/L) + 1/(LC))

# For parallel: Y_total = Y_bulk + Y_mlcc
# The poles of Z_total occur at the zeros of Y_total

# Construct the pole-residue model
from scipy.signal import TransferFunction, residue

# Bulk admittance: Y_bulk(s) = sC_bulk / (s² + s*ESR_bulk/L_bulk + 1/(L_bulk*C_bulk))
num_bulk = [C_bulk, 0]
den_bulk = [1, ESR_bulk/L_bulk, 1/(L_bulk*C_bulk)]
Y_bulk_tf = TransferFunction(num_bulk, den_bulk)

# MLCC admittance: Y_mlcc(s) = sC_mlcc / (s² + s*ESR_mlcc/L_mlcc + 1/(L_mlcc*C_mlcc))
num_mlcc = [C_mlcc, 0]
den_mlcc = [1, ESR_mlcc/L_mlcc, 1/(L_mlcc*C_mlcc)]
Y_mlcc_tf = TransferFunction(num_mlcc, den_mlcc)

# Find zeros (poles of impedance) by finding zeros of Y_bulk + Y_mlcc
# This is complex, so we'll use a numerical approach: find the frequency
# where d(Z)/df = 0 (maximum)

# We already found it: f_areson

print(f"\nNumerical anti-resonance frequency: {f_areson/1e6:.3f} MHz")

# Approximate geometric mean
f_geom_mean = np.sqrt(f_srf_bulk * f_srf_mlcc)
print(f"Geometric mean of SRFs: {f_geom_mean/1e6:.3f} MHz")
print(f"Difference: {abs(f_areson - f_geom_mean)/1e6:.3f} MHz")

# ==============================================================================
# PART 6: THE QUALITY FACTOR Q - SHARPNESS OF THE PEAK
# ==============================================================================

print("\n[STEP 6] Quality factor Q: How sharp is the anti-resonance peak?")
print("-" * 80)

print("""
The Q-factor (quality factor) of an LC resonance is:

  Q = ω₀ × L / R = 1 / (2ζ)

where ζ is the damping ratio.

High Q (low R): Sharp, needle-like peak (dangerous!)
Low Q (high R): Broad, gentle peak (easier to handle)

For anti-resonance, the effective Q depends on the parallel combination:

  Q_eff ≈ (ω_ar × L_eff) / R_eff

where L_eff and R_eff are effective inductance and resistance at f_ar.

Low ESR → High Q → Sharp peak → Exceeds target impedance by large margin
High ESR → Low Q → Broad peak → Easier to damp
""")

# Estimate Q at anti-resonance frequency
Z_ar_complex = Z_parallel[idx_areson]
L_eff = np.imag(Z_ar_complex) / omega[idx_areson]  # ≈ ωL
R_eff = np.real(Z_ar_complex)  # ≈ R

Q_eff = (omega[idx_areson] * abs(L_eff)) / (abs(R_eff) + 1e-10)

print(f"\nAt anti-resonance (f = {f_areson/1e6:.2f} MHz):")
print(f"  Effective R = {R_eff*1e3:.3f} mΩ")
print(f"  Effective L = {L_eff*1e9:.3f} nH")
print(f"  Q-factor ≈ {Q_eff:.1f}")

# ==============================================================================
# PART 7: VISUALIZATION
# ==============================================================================

print("\n[STEP 7] Visualize impedance, transfer function magnitude/phase")
print("-" * 80)

fig, axes = plt.subplots(2, 3, figsize=(16, 10))
fig.suptitle('Anti-Resonance in PDN: Transfer Function Analysis', fontsize=14, fontweight='bold')

# Plot 1: Magnitude (log-log)
ax = axes[0, 0]
ax.loglog(freq, Z_para_mag*1e3, 'k-', linewidth=2.5, label='|Z_parallel|')
ax.loglog(freq, np.abs(Z_bulk)*1e3, 'b--', linewidth=1.5, alpha=0.6, label='|Z_bulk|')
ax.loglog(freq, np.abs(Z_mlcc)*1e3, 'r--', linewidth=1.5, alpha=0.6, label='|Z_mlcc|')
ax.axhline(10, color='g', linestyle=':', linewidth=2, label='Target Z_target = 10 mΩ')
ax.plot(f_areson, Z_areson*1e3, 'ro', markersize=12, label=f'Anti-resonance\n@ {f_areson/1e6:.2f} MHz')
ax.axvline(f_srf_bulk, color='b', linestyle=':', alpha=0.5)
ax.axvline(f_srf_mlcc, color='r', linestyle=':', alpha=0.5)
ax.grid(True, which='both', alpha=0.3)
ax.set_ylabel('|Z| (mΩ)', fontsize=11)
ax.set_title('Impedance Magnitude', fontsize=11, fontweight='bold')
ax.legend(fontsize=9)
ax.set_ylim([0.5, 500])

# Plot 2: Magnitude (linear axis, zoomed to anti-resonance region)
ax = axes[0, 1]
freq_zoom = freq[(freq > 1e6) & (freq < 500e6)]
Z_zoom = Z_para_mag[(freq > 1e6) & (freq < 500e6)]
ax.semilogy(freq_zoom/1e6, Z_zoom*1e3, 'k-', linewidth=2.5)
ax.axhline(10, color='g', linestyle=':', linewidth=2, label='Target = 10 mΩ')
ax.plot(f_areson/1e6, Z_areson*1e3, 'ro', markersize=10)
ax.axvline(f_srf_bulk/1e6, color='b', linestyle='--', alpha=0.5, label=f'Bulk SRF = {f_srf_bulk/1e6:.1f} MHz')
ax.axvline(f_srf_mlcc/1e6, color='r', linestyle='--', alpha=0.5, label=f'MLCC SRF = {f_srf_mlcc/1e6:.1f} MHz')
ax.grid(True, which='both', alpha=0.3)
ax.set_xlabel('Frequency (MHz)', fontsize=11)
ax.set_ylabel('|Z| (mΩ)', fontsize=11)
ax.set_title('Anti-Resonance Region (Zoomed)', fontsize=11, fontweight='bold')
ax.legend(fontsize=9)
ax.set_ylim([1, 200])

# Plot 3: Phase response
ax = axes[0, 2]
ax.semilogx(freq, Z_para_phase, 'k-', linewidth=2.5, label='Phase(Z_parallel)')
ax.axhline(0, color='k', linestyle='-', alpha=0.2)
ax.axhline(90, color='g', linestyle='--', alpha=0.3)
ax.axhline(-90, color='g', linestyle='--', alpha=0.3)
ax.axvline(f_areson, color='r', linestyle='--', alpha=0.5)
ax.plot(f_areson, Z_para_phase[idx_areson], 'ro', markersize=10, label=f'Phase @ f_ar')
ax.grid(True, which='both', alpha=0.3)
ax.set_ylabel('Phase (degrees)', fontsize=11)
ax.set_title('Phase Response', fontsize=11, fontweight='bold')
ax.legend(fontsize=9)
ax.set_ylim([-180, 180])

# Plot 4: Real part (resistance)
ax = axes[1, 0]
ax.loglog(freq, np.abs(Z_para_real)*1e3, 'k-', linewidth=2.5)
ax.loglog(freq, np.abs(Z_bulk_real := np.real(Z_bulk))*1e3, 'b--', linewidth=1.5, alpha=0.6)
ax.loglog(freq, np.abs(Z_mlcc_real := np.real(Z_mlcc))*1e3, 'r--', linewidth=1.5, alpha=0.6)
ax.plot(f_areson, abs(Z_para_real[idx_areson])*1e3, 'ro', markersize=10)
ax.grid(True, which='both', alpha=0.3)
ax.set_ylabel('Re(Z) (mΩ)', fontsize=11)
ax.set_title('Resistance (Real Part)', fontsize=11, fontweight='bold')

# Plot 5: Imaginary part (reactance)
ax = axes[1, 1]
ax.semilogx(freq, Z_para_imag*1e3, 'k-', linewidth=2.5, label='Im(Z_parallel)')
ax.semilogx(freq, np.imag(Z_bulk)*1e3, 'b--', linewidth=1.5, alpha=0.6, label='Im(Z_bulk)')
ax.semilogx(freq, np.imag(Z_mlcc)*1e3, 'r--', linewidth=1.5, alpha=0.6, label='Im(Z_mlcc)')
ax.axhline(0, color='k', linestyle='-', alpha=0.3)
ax.axvline(f_areson, color='r', linestyle='--', alpha=0.5)
ax.plot(f_areson, Z_para_imag[idx_areson]*1e3, 'ro', markersize=10)
ax.grid(True, which='both', alpha=0.3)
ax.set_ylabel('Im(Z) (mΩ)', fontsize=11)
ax.set_title('Reactance (Imaginary Part)', fontsize=11, fontweight='bold')
ax.legend(fontsize=9)

# Plot 6: Impedance vector diagram at key frequencies
ax = axes[1, 2]

# Sample frequencies: at and near anti-resonance
f_samples = [f_srf_bulk, f_areson, f_srf_mlcc]
colors_sample = ['blue', 'red', 'darkred']
labels_sample = [f'@ Bulk SRF\n({f_srf_bulk/1e6:.2f} MHz)', 
                 f'@ Anti-resonance\n({f_areson/1e6:.2f} MHz)',
                 f'@ MLCC SRF\n({f_srf_mlcc/1e6:.2f} MHz)']

for f_s, c, lbl in zip(f_samples, colors_sample, labels_sample):
    idx_s = np.argmin(np.abs(freq - f_s))
    Z_s = Z_parallel[idx_s]
    Z_mag_s = abs(Z_s)*1e3
    Z_phase_s = np.angle(Z_s)
    
    # Plot as vector in complex plane
    ax.quiver(0, 0, np.real(Z_s)*1e3, np.imag(Z_s)*1e3, 
              angles='xy', scale_units='xy', scale=1, width=0.005,
              color=c, label=lbl)
    
    # Mark the tip
    ax.plot(np.real(Z_s)*1e3, np.imag(Z_s)*1e3, 'o', color=c, markersize=8)

ax.axhline(0, color='k', linestyle='-', alpha=0.2)
ax.axvline(0, color='k', linestyle='-', alpha=0.2)
ax.grid(True, alpha=0.3)
ax.set_xlabel('Re(Z) (mΩ)', fontsize=11)
ax.set_ylabel('Im(Z) (mΩ)', fontsize=11)
ax.set_title('Complex Impedance Vectors', fontsize=11, fontweight='bold')
ax.legend(fontsize=9, loc='upper right')
ax.axis('equal')
ax.set_xlim([-100, 150])
ax.set_ylim([-100, 150])

plt.tight_layout()
plt.savefig('anti_resonance_transfer_function.png', dpi=150, bbox_inches='tight')
print("✓ Saved plot to 'anti_resonance_transfer_function.png'")
plt.show()

# ==============================================================================
# PART 8: SUMMARY TABLE
# ==============================================================================

print("\n" + "="*80)
print("SUMMARY: ANTI-RESONANCE TRANSFER FUNCTION")
print("="*80)

summary_data = f"""
CAPACITOR SPECIFICATIONS:
┌─────────────────────┬──────────────┬─────────────┐
│ Parameter           │ Bulk (470µF) │ MLCC (10µF) │
├─────────────────────┼──────────────┼─────────────┤
│ Capacitance         │  470 µF      │  10 µF      │
│ ESR                 │  50 mΩ       │  5 mΩ       │
│ ESL                 │  5.0 nH      │  0.7 nH     │
│ SRF (f₀)            │  {f_srf_bulk/1e6:6.2f} MHz  │  {f_srf_mlcc/1e6:5.2f} MHz  │
│ Min |Z| @ SRF       │  50 mΩ       │  5 mΩ       │
└─────────────────────┴──────────────┴─────────────┘

PARALLEL COMBINATION BEHAVIOR:
────────────────────────────────────────────────────
Frequency Region              Behavior
────────────────────────────────────────────────────
DC − Bulk SRF ({f_srf_bulk/1e6:.2f} MHz)    Bulk dominates (capacitive)
                            Z dominated by 1/(jωC_bulk)

@ Bulk SRF ({f_srf_bulk/1e6:.2f} MHz)        Bulk impedance = 50 mΩ (min)
                            MLCC still capacitive (~100 mΩ)
                            Parallel Z ≈ 33 mΩ (good)

Between SRFs ({f_srf_bulk/1e6:.2f}−{f_srf_mlcc/1e6:.2f} MHz)      **ANTI-RESONANCE ZONE**
                            Bulk: inductive (jωL)
                            MLCC: capacitive (−j/(ωC))
                            Form LC tank → SHARP PEAK

@ Anti-Resonance ({f_areson/1e6:.2f} MHz)  |Z| = {Z_areson*1e3:.2f} mΩ ← EXCEEDS TARGET!
                            This is where PDN fails

@ MLCC SRF ({f_srf_mlcc/1e6:.2f} MHz)        MLCC impedance = 5 mΩ (min)
                            Bulk very inductive
                            Parallel Z ≈ 5 mΩ (good)

Above MLCC SRF ({f_srf_mlcc/1e6:.2f} MHz)    Both inductive, impedance rises
────────────────────────────────────────────────────

TARGET IMPEDANCE BUDGET:  10 mΩ
ANTI-RESONANCE PEAK:      {Z_areson*1e3:.2f} mΩ
EXCEEDANCE:               {Z_areson*1e3/10:.1f}× TARGET ← FAILURE

TRANSFER FUNCTION INTERPRETATION:
  Z(s) has poles (impedance zeros → zero admittance) at two frequencies:
    1. SERIES RESONANCE: at each capacitor's SRF (minimum impedance)
    2. PARALLEL RESONANCE: between the two SRFs (MAXIMUM impedance) ← THE KILLER

The parallel/anti-resonance pole is where the LC tank effect dominates.
This is a second-order system with complex pole pair, quality factor Q ≈ {Q_eff:.1f}.
"""

print(summary_data)

# ==============================================================================
# PART 9: PHYSICAL INTERPRETATION - WHY THIS MATTERS FOR SERDES
# ==============================================================================

print("\n" + "="*80)
print("INTERVIEW INSIGHT: How This Kills Signal Integrity")
print("="*80)

insight = f"""
THE MECHANISM:

1. NOISE GENERATION (SSO - Simultaneous Switching Output):
   When the ASIC switches all outputs simultaneously, current pulses flow through
   the PDN. At frequency f_ar = {f_areson/1e6:.2f} MHz, the SSO current hits the
   anti-resonance peak and generates maximum noise:
   
     V_noise = Z(f_ar) × I_SSO
             = {Z_areson*1e3:.2f} mΩ × 10 A (typical SSO)
             = {Z_areson*1e3 * 10:.0f} mV of noise on the supply!
   
   For a 1.0 V supply with 3% ripple budget (30 mV), this noise
   CONSUMES THE ENTIRE MARGIN.

2. MODE CONVERSION:
   This common-mode noise couples into the differential signal pair through:
   - Via coupling and cross-coupling
   - Fiber-weave-induced asymmetry  
   - Connector parasitics
   - Ground return path splitting
   
   The conversion factor depends on asymmetry, typically 10−30% coupling.

3. CHANNEL DEGRADATION:
   The coupled noise appears as insertion loss ripple in the channel at f_ar.
   The ripple causes ISI (inter-symbol interference) that extends over many
   symbol periods and is hard to equalize.
   
   Result: Eye height drops by 10−50 mV, jitter increases by 5−20 ps.

4. COM FAILURE:
   IEEE 802.3 Channel Operating Margin (COM) integrates insertion loss, jitter,
   and noise. PDN-induced insertion loss ripple at f_ar directly reduces COM.
   
   If COM drops below 3 dB, the link fails compliance.
   
THE FIX:

Option 1: ESR Damping
  Add 10−50 mΩ series resistance to one cap bank.
  This decreases Q-factor from {Q_eff:.1f} to ~2, broadens the peak.
  
Option 2: SRF Staggering
  Use 4 cap tiers with SRFs at 2 MHz, 5 MHz, 15 MHz, 25 MHz.
  Spreads anti-resonances, prevents any single peak from dominating.
  
Option 3: Active VRM
  Voltage regulator module with feedback control actively cancels PDN resonances.
  Industry standard for high-power ASICs (Gbit Ethernet, 5G SoCs).
"""

print(insight)

print("\n" + "="*80)
print("END OF ANTI-RESONANCE ANALYSIS")
print("="*80)
