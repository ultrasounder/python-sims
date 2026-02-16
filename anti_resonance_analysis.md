================================================================================
ANTI-RESONANCE TRANSFER FUNCTION DERIVATION
================================================================================

[STEP 1] Understand the lumped-element circuit model
--------------------------------------------------------------------------------

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


[STEP 2] Single RLC series impedance (single capacitor bank)
--------------------------------------------------------------------------------

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


[STEP 3] Compute single bank impedance (10µF MLCC)
--------------------------------------------------------------------------------
10µF MLCC Parameters:
  C = 10.0 µF
  ESR = 5.00 mΩ
  L = 0.70 nH
  SRF (f₀) = 1/(2π√LC) = 1.90 MHz
  |Z| at SRF = ESR = 5.00 mΩ
  Measured minimum |Z| = 5.000 mΩ at f = 1.90 MHz

[STEP 4] Two-bank system: Bulk + MLCC → Anti-Resonance
--------------------------------------------------------------------------------
Bulk capacitor (470µF):
  SRF = 0.10 MHz
Ceramic MLCC (10µF):
  SRF = 1.90 MHz

Parallel (two-bank) system:
  Anti-resonance frequency f_ar = 10000.00 MHz
  Anti-resonance impedance |Z_ar| = 38580.96 mΩ
  Ratio to target (10 mΩ): 3858.1×

[STEP 5] Express as transfer function H(s) = Z(s) (impedance function)
--------------------------------------------------------------------------------

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


Numerical anti-resonance frequency: 10000.000 MHz
Geometric mean of SRFs: 0.444 MHz
Difference: 9999.556 MHz

[STEP 6] Quality factor Q: How sharp is the anti-resonance peak?
--------------------------------------------------------------------------------

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


At anti-resonance (f = 10000.00 MHz):
  Effective R = 4.601 mΩ
  Effective L = 0.614 nH
  Q-factor ≈ 8384.6

[STEP 7] Visualize impedance, transfer function magnitude/phase
--------------------------------------------------------------------------------
Ignoring fixed x limits to fulfill fixed data aspect with adjustable data limits.
Ignoring fixed x limits to fulfill fixed data aspect with adjustable data limits.
✓ Saved plot to 'anti_resonance_transfer_function.png'

================================================================================
SUMMARY: ANTI-RESONANCE TRANSFER FUNCTION
================================================================================

CAPACITOR SPECIFICATIONS:
┌─────────────────────┬──────────────┬─────────────┐
│ Parameter           │ Bulk (470µF) │ MLCC (10µF) │
├─────────────────────┼──────────────┼─────────────┤
│ Capacitance         │  470 µF      │  10 µF      │
│ ESR                 │  50 mΩ       │  5 mΩ       │
│ ESL                 │  5.0 nH      │  0.7 nH     │
│ SRF (f₀)            │    0.10 MHz  │   1.90 MHz  │
│ Min |Z| @ SRF       │  50 mΩ       │  5 mΩ       │
└─────────────────────┴──────────────┴─────────────┘

PARALLEL COMBINATION BEHAVIOR:
────────────────────────────────────────────────────
Frequency Region              Behavior
────────────────────────────────────────────────────
DC − Bulk SRF (0.10 MHz)    Bulk dominates (capacitive)
                            Z dominated by 1/(jωC_bulk)

@ Bulk SRF (0.10 MHz)        Bulk impedance = 50 mΩ (min)
                            MLCC still capacitive (~100 mΩ)
                            Parallel Z ≈ 33 mΩ (good)

Between SRFs (0.10−1.90 MHz)      **ANTI-RESONANCE ZONE**
                            Bulk: inductive (jωL)
                            MLCC: capacitive (−j/(ωC))
                            Form LC tank → SHARP PEAK

@ Anti-Resonance (10000.00 MHz)  |Z| = 38580.96 mΩ ← EXCEEDS TARGET!
                            This is where PDN fails

@ MLCC SRF (1.90 MHz)        MLCC impedance = 5 mΩ (min)
                            Bulk very inductive
                            Parallel Z ≈ 5 mΩ (good)

Above MLCC SRF (1.90 MHz)    Both inductive, impedance rises
────────────────────────────────────────────────────

TARGET IMPEDANCE BUDGET:  10 mΩ
ANTI-RESONANCE PEAK:      38580.96 mΩ
EXCEEDANCE:               3858.1× TARGET ← FAILURE

TRANSFER FUNCTION INTERPRETATION:
  Z(s) has poles (impedance zeros → zero admittance) at two frequencies:
    1. SERIES RESONANCE: at each capacitor's SRF (minimum impedance)
    2. PARALLEL RESONANCE: between the two SRFs (MAXIMUM impedance) ← THE KILLER

The parallel/anti-resonance pole is where the LC tank effect dominates.
This is a second-order system with complex pole pair, quality factor Q ≈ 8384.6.


================================================================================
INTERVIEW INSIGHT: How This Kills Signal Integrity
================================================================================

THE MECHANISM:

1. NOISE GENERATION (SSO - Simultaneous Switching Output):
   When the ASIC switches all outputs simultaneously, current pulses flow through
   the PDN. At frequency f_ar = 10000.00 MHz, the SSO current hits the
   anti-resonance peak and generates maximum noise:
   
     V_noise = Z(f_ar) × I_SSO
             = 38580.96 mΩ × 10 A (typical SSO)
             = 385810 mV of noise on the supply!
   
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
  This decreases Q-factor from 8384.6 to ~2, broadens the peak.
  
Option 2: SRF Staggering
  Use 4 cap tiers with SRFs at 2 MHz, 5 MHz, 15 MHz, 25 MHz.
  Spreads anti-resonances, prevents any single peak from dominating.
  
Option 3: Active VRM
  Voltage regulator module with feedback control actively cancels PDN resonances.
  Industry standard for high-power ASICs (Gbit Ethernet, 5G SoCs).


================================================================================
END OF ANTI-RESONANCE ANALYSIS
================================================================================