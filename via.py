import numpy as np

# --- 1. Physics Constants & Setup ---
freq_start = 10e6   # 10 MHz
freq_stop  = 20e9   # 20 GHz
num_points = 2000
freq = np.linspace(freq_start, freq_stop, num_points)
w = 2 * np.pi * freq

# Impedance Targets
Z0_line = 85.0      # Target Diff Impedance (Ohms)
Z0_ref  = 85.0      # Port Reference Impedance

# The "Bad" Via (Capacitive Discontinuity)
# A 25mil anti-pad on a thick PCB creates excess capacitance.
# Modeled as a shunt capacitor in the middle of the line.
C_via = 0.6e-12     # 0.6 pF (This causes the reflection)

# Transmission Line Properties (Megtron-6ish)
# Loss model: alpha = a_cond * sqrt(f) + a_diel * f
L_total_inch = 10.0
L_seg_meter  = (L_total_inch / 2.0) * 0.0254 # Split into two 5-inch segments

# Loss coefficients (Approx for stripline)
alpha_cond = 2.5e-7  # Skin effect
alpha_diel = 5.0e-11 # Dielectric loss

# --- 2. ABCD Matrix Calculation ---
# We calculate the matrix for: Line (5") -> Via (Cap) -> Line (5")

# Propagation Constant (Gamma)
alpha = alpha_cond * np.sqrt(freq) + alpha_diel * freq
beta  = w / (1.5e8) # Approx phase velocity (Dk ~ 4)
gamma = alpha + 1j * beta

# T-Line ABCD Matrix (Half Channel)
# A = cosh(gl), B = Z0*sinh(gl), C = sinh(gl)/Z0, D = cosh(gl)
cosh_gl = np.cosh(gamma * L_seg_meter)
sinh_gl = np.sinh(gamma * L_seg_meter)

A_t = cosh_gl
B_t = Z0_line * sinh_gl
C_t = sinh_gl / Z0_line
D_t = cosh_gl

# Via ABCD Matrix (Shunt Capacitor)
# A=1, B=0, C=jwc, D=1
A_v = np.ones_like(freq)
B_v = np.zeros_like(freq)
C_v = 1j * w * C_via
D_v = np.ones_like(freq)

# Cascade Matrices: [T_total] = [T_line] * [T_via] * [T_line]
# MatMul: [[A1, B1], [C1, D1]] * [[A2, B2], [C2, D2]]
# Step 1: Line * Via
A_tv = A_t*A_v + B_t*C_v
B_tv = A_t*B_v + B_t*D_v
C_tv = C_t*A_v + D_t*C_v
D_tv = C_t*B_v + D_t*D_v

# Step 2: (Line*Via) * Line
A_sys = A_tv*A_t + B_tv*C_t
B_sys = A_tv*B_t + B_tv*D_t
C_sys = C_tv*A_t + D_tv*C_t
D_sys = C_tv*B_t + D_tv*D_t

# --- 3. Convert ABCD to S-Parameters ---
# Formula for S21 and S11 from ABCD normalized to Z0_ref
denom = A_sys + B_sys/Z0_ref + C_sys*Z0_ref + D_sys

s21 = 2.0 / denom
s11 = (A_sys + B_sys/Z0_ref - C_sys*Z0_ref - D_sys) / denom
# Assume symmetry/reciprocity
s12 = s21
s22 = s11

# --- 4. Write Touchstone (.s2p) File ---
filename = "Bad_Via_Model.s2p"
with open(filename, 'w') as f:
    f.write("! Validation Capstone Generated Data\n")
    f.write("! Freq: 10MHz - 20GHz, 10 inch trace, 0.6pF Via Cap\n")
    f.write("# Hz S RI R 85.0\n") # RI = Real/Imag, Ref Imp = 85 Ohms
    
    for i in range(len(freq)):
        row = f"{freq[i]:.0f} {s11[i].real:.5f} {s11[i].imag:.5f} {s21[i].real:.5f} {s21[i].imag:.5f} {s12[i].real:.5f} {s12[i].imag:.5f} {s22[i].real:.5f} {s22[i].imag:.5f}\n"
        f.write(row)

print(f"Generated {filename}. Insertion Loss @ 8GHz: {20*np.log10(np.abs(s21[800])):.2f} dB")