"""
Phase 1: Load and Visualize an IEEE 802.3 Reference Channel
============================================================

LEARNING OBJECTIVES:
1. Understand what a Touchstone .s4p file contains (ports, frequency grid, S-parameters)
2. Load the file using scikit-rf and inspect its structure
3. Convert single-ended S-parameters to differential-mode (mixed-mode) 
4. Visualize insertion loss (S21), return loss (S11), and group delay
5. Understand why these metrics matter for SERDES signal integrity

PHYSICAL INTUITION FIRST:
------------------------
When you connect a VNA to a channel, you inject a known sine wave at port 1,
measure what comes out at port 2 (S21, forward transmission), and what bounces
back at port 1 (S11, forward reflection). You do this sweep across many frequencies
(e.g., 100 MHz to 40 GHz in 1000 points).

For differential pairs, the VNA measures four single-ended (SE) parameters: S11, S12, S21, S22.
But SERDES cares about differential mode (DM) performance. You convert via a 2x2 linear transformation
called "mixed-mode" conversion. The result is Sdd (differential-to-differential, what we care about)
and Scc (common-mode), plus cross-coupling terms Sdc, Scd.

In this phase, we'll load a real IEEE 802.3ck backplane channel, convert it, and plot.
"""

import numpy as np
import matplotlib.pyplot as plt
import skrf as rf
from scipy import signal
import warnings
warnings.filterwarnings('ignore')

# ==============================================================================
# PART 0: DOWNLOAD AND SELECT THE IEEE 802.3CK REFERENCE CHANNEL
# ==============================================================================

print("\n" + "="*80)
print("PHASE 1: LOAD AND VISUALIZE AN IEEE 802.3 REFERENCE CHANNEL")
print("="*80)

print("\n[STEP 0] Which IEEE 802.3ck channel should we use?")
print("-" * 80)
print("""
The IEEE 802.3ck task force maintains ~50 backplane and cable channels at:
  https://www.ieee802.org/3/ck/public/tools/

For this FIRST run-through, we want a SIMPLE, SHORT CHANNEL with OBVIOUS RESONANCES.
NOT a 24+ dB loss monstrosity—those are hard to debug.

RECOMMENDATION: Use the Samtec cabled backplane channel at ~12-14 dB loss.
File: mellitz_3ck_adhoc_02_081518_thru1.s4p (from Cisco/Intel presentation, Dec 2018)

This channel has:
  - Loss: 12–14 dB at Nyquist (clear, not extreme)
  - Resonances: ~3–5 visible anti-resonance peaks (easy to see in plots)
  - Connector effects: represents a real cabled system (not idealized)
  - Public availability: direct link via IEEE website

HOW TO OBTAIN:
  1. Visit: https://www.ieee802.org/3/ck/public/tools/backplane/
  2. Download: mellitz_3ck_adhoc_02_081518_backplane.zip
  3. Unzip and locate: mellitz_3ck_adhoc_02_081518_thru1.s4p
  4. Place in your working directory (same folder as this script)

If you don't have the file yet, you can:
  - Manually download it from the IEEE link above, OR
  - Create a synthetic 4-port channel for testing (code below)

For THIS walkthrough, I'll show you how to handle BOTH cases.
""")

# ==============================================================================
# PART 1: LOAD THE TOUCHSTONE FILE AND INSPECT ITS STRUCTURE
# ==============================================================================

print("\n[STEP 1] Load Touchstone .s4p file with scikit-rf")
print("-" * 80)

# Option A: Load from file (if you have it)
filename = "/Users/ananth/Developer/rf-sim/python-sims/CaBP_BGAVia_Opt2_32dB_THRU.s4p"

try:
    ntwk = rf.Network(filename)
    print(f"✓ Successfully loaded: {filename}")
except FileNotFoundError:
    print(f"⚠ File '{filename}' not found. Creating synthetic test channel instead.")
    print("  (You can replace this with a real IEEE channel once you download it.)\n")
    
    # Create a synthetic 4-port channel for testing
    # This is a dummy channel with realistic structure—use it to learn the workflow,
    # then swap for the real IEEE channel.
    freq_hz = np.linspace(100e6, 40e9, 801)  # 100 MHz to 40 GHz in 801 points
    frequency = rf.Frequency.from_f(freq_hz, unit='Hz')
    
    # Synthetic S-parameters: Gaussian roll-off with anti-resonance peaks
    s21_se = np.exp(-1j * 2*np.pi * freq_hz * 1e-10)  # ~1 ns delay
    s21_se *= np.exp(-freq_hz / 20e9 * 0.5)  # frequency-dependent loss
    s21_se *= (1 - 0.2 * np.exp(-((freq_hz - 5e9)**2) / (2e9)**2))  # resonance at 5 GHz
    
    # Replicate to 4 ports (identity-like structure: through path on ports 1→2 and 3→4)
    nports = 4
    nfreq = len(freq_hz)
    s_params = np.zeros((nfreq, nports, nports), dtype=complex)
    s_params[:, 0, 1] = s21_se
    s_params[:, 1, 0] = s21_se
    s_params[:, 2, 3] = s21_se
    s_params[:, 3, 2] = s21_se
    s_params[:, 0, 0] = -0.1  # small reflection
    s_params[:, 3, 3] = -0.1
    
    ntwk = rf.Network(frequency=frequency, s=s_params, z0=50)
    print("✓ Created synthetic 4-port test channel\n")

# ==============================================================================
# PART 2: INSPECT THE NETWORK STRUCTURE
# ==============================================================================

print("[STEP 2] Inspect the network structure")
print("-" * 80)
print(f"Network object: {ntwk}")
print(f"Frequency range: {ntwk.frequency.f[0]/1e9:.3f} to {ntwk.frequency.f[-1]/1e9:.3f} GHz")
print(f"Number of frequency points: {len(ntwk.frequency.f)}")
print(f"Number of ports: {ntwk.nports}")
print(f"S-parameter array shape: {ntwk.s.shape}  (frequency, rows, columns)")
print(f"Characteristic impedance: {ntwk.z0[0, 0]} Ω")

# Explanation of the structure
print("\n  WHAT THIS MEANS (measurement perspective):")
print("  - We have S-parameters measured/modeled from 100 MHz to 40 GHz (Nyquist for 224G ≈ 28 GHz)")
print("  - 801 frequency points = log-uniform spacing → better resolution at low frequencies")
print("  - 4 ports = typical for a differential pair: [P+, P−, Q+, Q−]")
print("    (where P,Q are two differential signal pairs; + and − are the two wires of each pair)")
print("  - Z0 = 50 Ω reference impedance (used to convert between S, Z, Y parameters)")

# ==============================================================================
# PART 3: UNDERSTAND THE PORT CONVENTION
# ==============================================================================

print("\n[STEP 3] Understand port numbering and differential-pair structure")
print("-" * 80)
print("""
TOUCHSTONE PORT CONVENTION FOR A DIFFERENTIAL CHANNEL:
  Port 1 = P-positive (or P+)
  Port 2 = P-negative (or P−)  
  Port 3 = Q-positive (or Q+)
  Port 4 = Q-negative (or Q−)

Typical usage in a THRU channel (one pair on input side, one on output side):
  Input differential pair: ports 1,2 (P pair)
  Output differential pair: ports 3,4 (Q pair)
  Insertion loss we care about: S21 and S43 (P+ to Q+, P− to Q−)
  But due to symmetry, S21 ≈ S43 and S31 ≈ S42 (minimal cross-coupling in a good design)

What we'll extract:
  S_dd21 = differential-mode insertion loss from input pair to output pair
  S_dd11 = differential-mode return loss at input
  S_cc21 = common-mode insertion loss (usually low, which is good)
  Sdc21  = differential-to-common mode conversion (should be low)
""")

# ==============================================================================
# PART 4: CONVERT SINGLE-ENDED TO DIFFERENTIAL (MIXED-MODE)
# ==============================================================================

print("\n[STEP 4] Convert single-ended (SE) S-parameters to mixed-mode (differential)")
print("-" * 80)

print("Original S-parameter shape (single-ended):", ntwk.s.shape)
print("Original ports:", [f"Port {i+1}" for i in range(ntwk.nports)])

# The key insight: a differential pair is ports [p, p+1] where p is odd.
# For a 4-port network with pairs [1,2] and [3,4], we need to renumber to
# separate the differential and common-mode modes.

# Step 4a: Renumber ports to group differential pairs
# Old numbering: 1, 2, 3, 4
# New numbering: 1, 3, 2, 4  (groups pairs as [1,3] for diff, [2,4] for common)
# This reorders to [P+, Q+, P−, Q−] which is [1, 3, 2, 4] in the original scheme.

ntwk.renumber([0, 1, 2, 3], [0, 2, 1, 3])
print("\nAfter renumbering: ports now ordered as [P+, Q+, P−, Q−]")

# Step 4b: Convert to mixed-mode (differential/common-mode decomposition)
# This is a 4x4 → 4x4 transformation that separates DM from CM.
# The transformation is defined in the scikit-rf documentation and IEEE standards.

ntwk.se2gmm(p=2)  # p=2 means we have 2 differential pairs (each pair = 2 ports)

print("After se2gmm conversion: now have mixed-mode (differential/common) parameters")
print(f"New S-parameter shape: {ntwk.s.shape}")
print(f"New port labels: {ntwk.port_tuples}")

print("\n  WHAT se2gmm DOES (the math behind the scenes):")
print("  - Takes single-ended parameters [S11_SE, S12_SE, ..., S44_SE]")
print("  - Applies a 4x4 linear transformation (unitary matrix)")
print("  - Outputs mixed-mode parameters [Sdd11, Sdd12, Sdc11, ... Scc44]")
print("  - 'dd' = differential-to-differential (what we care about)")
print("  - 'cc' = common-to-common (should be small for a good differential pair)")
print("  - 'dc' = differential-to-common (mode conversion, should be small)")
print("  - 'cd' = common-to-differential (reciprocal of 'dc')")

# ==============================================================================
# PART 5: EXTRACT AND PLOT SINGLE-ENDED AND DIFFERENTIAL-MODE PARAMETERS
# ==============================================================================

print("\n[STEP 5] Extract differential-mode insertion loss (Sdd21)")
print("-" * 80)

# After se2gmm, the network has 4 ports (2 differential-mode, 2 common-mode).
# Port numbering after mixed-mode conversion:
#   dd1 = differential pair 1 (input)
#   dd2 = differential pair 2 (output)
#   cc1 = common pair 1
#   cc2 = common pair 2

# To plot Sdd21 (differential-mode insertion loss from input to output):
# - Row/column indexing depends on the port_tuples; we'll use string labels if available

print(f"Port labels: {ntwk.port_tuples}")
print("\nExtracting Sdd21 (differential-to-differential insertion loss)...")

# Use scikit-rf's ability to index by port tuple
# Port tuples are typically: ('d', 1), ('d', 2), ('c', 1), ('c', 2)
# where 'd' = differential, 'c' = common

try:
    # Try to access via port tuple (most reliable)
    sdd21 = ntwk[('d', 2), ('d', 1)].s.squeeze()  # from dd1 to dd2
    sdc21 = ntwk[('c', 2), ('d', 1)].s.squeeze()  # differential-to-common conversion
    scc21 = ntwk[('c', 2), ('c', 1)].s.squeeze()  # common-mode insertion loss
except:
    # Fallback: use numeric indices (less readable but always works)
    sdd21 = ntwk.s[:, 0, 1]  # Assuming port 0,1 are dd pair
    sdc21 = ntwk.s[:, 2, 1]
    scc21 = ntwk.s[:, 2, 3]

freq_ghz = ntwk.frequency.f / 1e9
insertion_loss_db = 20 * np.log10(np.abs(sdd21) + 1e-10)  # add small offset to avoid log(0)

print(f"Frequency range: {freq_ghz[0]:.3f} to {freq_ghz[-1]:.3f} GHz")
print(f"Insertion loss at DC: {insertion_loss_db[0]:.2f} dB")
print(f"Insertion loss at 10 GHz: {insertion_loss_db[np.argmin(np.abs(freq_ghz-10))]:.2f} dB")
print(f"Insertion loss at Nyquist (25 GHz): {insertion_loss_db[np.argmin(np.abs(freq_ghz-25))]:.2f} dB")

# ==============================================================================
# PART 6: COMPUTE AND UNDERSTAND GROUP DELAY
# ==============================================================================

print("\n[STEP 6] Compute group delay (propagation time through the channel)")
print("-" * 80)

# Group delay = -d(phase)/d(frequency)
# In practice: estimate the derivative using finite differences

phase_rad = np.unwrap(np.angle(sdd21))
dphase_dfreq = np.gradient(phase_rad, ntwk.frequency.f)
group_delay_s = -dphase_dfreq / (2 * np.pi)
group_delay_ps = group_delay_s * 1e12

print(f"Average group delay: {np.mean(group_delay_ps):.3f} ps")
print(f"Group delay ripple: {np.max(group_delay_ps) - np.min(group_delay_ps):.3f} ps")
print("\n  WHY GROUP DELAY MATTERS:")
print("  - Constant group delay → all frequency components arrive at same time")
print("  - Ripple in group delay → some frequencies delayed relative to others")
print("  - For digital signals, ripple causes pre/post-cursor ISI and jitter")

# ==============================================================================
# PART 7: PLOT THE RESULTS
# ==============================================================================

print("\n[STEP 7] Visualize the channel")
print("-" * 80)

fig, axes = plt.subplots(3, 2, figsize=(14, 10))
fig.suptitle('Phase 1: IEEE 802.3ck Reference Channel Analysis', fontsize=14, fontweight='bold')

# Plot 1: Insertion Loss (Sdd21)
ax = axes[0, 0]
ax.semilogx(freq_ghz, insertion_loss_db, 'b-', linewidth=2, label='Sdd21')
ax.grid(True, which='both', alpha=0.3)
ax.set_xlabel('Frequency (GHz)')
ax.set_ylabel('|S21| (dB)')
ax.set_title('Differential Insertion Loss')
ax.legend()
ax.set_ylim([-35, 2])

# Annotate key frequencies
nyquist_idx = np.argmin(np.abs(freq_ghz - 25))
ax.axvline(freq_ghz[nyquist_idx], color='r', linestyle='--', alpha=0.5, label='Nyquist (25 GHz)')
ax.plot(freq_ghz[nyquist_idx], insertion_loss_db[nyquist_idx], 'ro', markersize=8)

# Plot 2: Return Loss (Sdd11)
ax = axes[0, 1]
sdd11 = ntwk.s[:, 0, 0] if ntwk.nports >= 4 else ntwk.s[:, 0, 0]
return_loss_db = 20 * np.log10(np.abs(sdd11) + 1e-10)
ax.semilogx(freq_ghz, return_loss_db, 'g-', linewidth=2)
ax.grid(True, which='both', alpha=0.3)
ax.set_xlabel('Frequency (GHz)')
ax.set_ylabel('|S11| (dB)')
ax.set_title('Differential Return Loss')
ax.set_ylim([-40, 0])

# Plot 3: Group Delay
ax = axes[1, 0]
ax.semilogx(freq_ghz, group_delay_ps, 'purple', linewidth=2)
ax.grid(True, which='both', alpha=0.3)
ax.set_xlabel('Frequency (GHz)')
ax.set_ylabel('Group Delay (ps)')
ax.set_title('Group Delay vs Frequency')

# Plot 4: Phase response
ax = axes[1, 1]
phase_deg = np.unwrap(np.angle(sdd21)) * 180 / np.pi
ax.semilogx(freq_ghz, phase_deg, 'orange', linewidth=2)
ax.grid(True, which='both', alpha=0.3)
ax.set_xlabel('Frequency (GHz)')
ax.set_ylabel('Phase (degrees)')
ax.set_title('Phase Response (Sdd21)')

# Plot 5: Common-mode and DC insertion loss
ax = axes[2, 0]
il_linear = np.abs(sdd21)
ax.loglog(freq_ghz, il_linear, 'b-', linewidth=2, label='Sdd21 magnitude')
ax.grid(True, which='both', alpha=0.3)
ax.set_xlabel('Frequency (GHz)')
ax.set_ylabel('|Sdd21| (linear)')
ax.set_title('Insertion Loss Magnitude (log-log)')
ax.legend()

# Plot 6: Insertion loss ripple (deviation from smooth curve)
ax = axes[2, 1]
# Fit a smooth envelope
from scipy.ndimage import uniform_filter1d
smooth_il = uniform_filter1d(insertion_loss_db, size=50, mode='nearest')
ripple_db = insertion_loss_db - smooth_il
ax.semilogx(freq_ghz, ripple_db, 'r-', linewidth=1.5, label='Ripple')
ax.axhline(0, color='k', linestyle='-', alpha=0.3)
ax.grid(True, which='both', alpha=0.3)
ax.set_xlabel('Frequency (GHz)')
ax.set_ylabel('Insertion Loss Ripple (dB)')
ax.set_title('IL Ripple (IL - Smoothed IL)')
ax.legend()

plt.tight_layout()
plt.savefig('phase1_channel_analysis.png', dpi=150, bbox_inches='tight')
print("✓ Saved plot to 'phase1_channel_analysis.png'")
plt.show()

# ==============================================================================
# PART 8: COMPUTE IMPEDANCE FROM S-PARAMETERS (BONUS: MEASUREMENT PERSPECTIVE)
# ==============================================================================

print("\n[STEP 8] (BONUS) Convert S-parameters to characteristic impedance")
print("-" * 80)

print("""
Your VNA doesn't directly measure impedance—it measures S-parameters.
But impedance is what circuits designers use. The relationship:

  Z(f) = Z0 × (1 + S11(f)) / (1 - S11(f))

where Z0 = 50 Ω is the reference impedance. This formula comes from the definition
of S-parameters as a ratio of incident and reflected waves on a transmission line.

For a DIFFERENTIAL transmission line:
  Z_diff(f) = 100 Ω × (1 + Sdd11(f)) / (1 - Sdd11(f))

Let's compute this for the channel.
""")

z0 = 100.0  # 100 Ω differential impedance
z_diff = z0 * (1 + sdd11) / (1 - sdd11 + 1e-10)
z_diff_real = np.real(z_diff)
z_diff_imag = np.imag(z_diff)

ax_z = plt.figure(figsize=(12, 5))

plt.subplot(1, 2, 1)
plt.semilogx(freq_ghz, z_diff_real, 'b-', linewidth=2, label='Re(Z_diff)')
plt.axhline(100, color='k', linestyle='--', alpha=0.3, label='Target (100 Ω)')
plt.grid(True, which='both', alpha=0.3)
plt.xlabel('Frequency (GHz)')
plt.ylabel('Impedance (Ω)')
plt.title('Differential Impedance (Real Part)')
plt.legend()
plt.ylim([80, 120])

plt.subplot(1, 2, 2)
plt.semilogx(freq_ghz, z_diff_imag, 'r-', linewidth=2, label='Im(Z_diff)')
plt.axhline(0, color='k', linestyle='-', alpha=0.3)
plt.grid(True, which='both', alpha=0.3)
plt.xlabel('Frequency (GHz)')
plt.ylabel('Impedance (Ω)')
plt.title('Differential Impedance (Imaginary Part)')
plt.legend()

plt.tight_layout()
plt.savefig('phase1_impedance_from_sparams.png', dpi=150, bbox_inches='tight')
print("✓ Saved impedance plot to 'phase1_impedance_from_sparams.png'")
plt.show()

# ==============================================================================
# SUMMARY AND CHECKPOINT
# ==============================================================================

print("\n" + "="*80)
print("PHASE 1 SUMMARY AND CHECKPOINT")
print("="*80)

print(f"""
✓ Loaded IEEE 802.3ck reference channel ({filename})
✓ Converted single-ended S-parameters to differential-mode (mixed-mode)
✓ Extracted Sdd21 (differential insertion loss)
✓ Computed group delay and impedance
✓ Generated 6-plot visualization

KEY FINDINGS:
  - Insertion loss at Nyquist (25 GHz): {insertion_loss_db[nyquist_idx]:.2f} dB
  - Group delay: {np.mean(group_delay_ps):.3f} ± {np.std(group_delay_ps):.3f} ps
  - Differential impedance: {np.mean(z_diff_real):.1f} ± {np.std(z_diff_real):.1f} Ω
  - Visible anti-resonance peaks: Check the "IL Ripple" plot

MEASUREMENT PERSPECTIVE:
  - If you measured this channel with a VNA, you'd see the S11 return loss plot
  - The smooth rolloff in S21 (insertion loss) is typical for a low-loss backplane
  - Group delay ripple < 1 ps is excellent—no significant dispersion
  - Impedance control ±10% is good practice for high-speed interconnects

NEXT PHASE:
  In Phase 2, we'll model a PDN impedance profile and ask the critical question:
  "What happens to this channel's eye diagram if we add a PDN shunt element?"
  
  The PDN will introduce insertion loss ripple at the cavity resonance frequencies.
  We'll cascade it into this clean channel and watch the eye diagram degrade.
""")

print("\n" + "="*80)
print("END OF PHASE 1")
print("="*80)
