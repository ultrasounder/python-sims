from matplotlib import pyplot as plt
import skrf as rf
from skrf.media import DefinedGammaZ0
import numpy as np

#1. create Frequency Vector
# PCIe Gen4 needs visibility up to Nyquist (8 GHz) and beyond. Lets go to 20 GHz
freq = rf.Frequency(0.01, 20, 1001, 'GHz')
omega = freq.w # Angular frequency (2*pi*f)

#2. Define the Transmission line media
# we simulate the Stripline in Megtron-6 (From PCIe Gen4 design Spec)
# Df = 0.002, Dk = 3.7
Dk = 3.7 # Dielectric constant of Megtron 6
Df = 0.002 # Dissipation factor (Low loss)
z0_target = 85 # PCI-E target impedance

#3. calculate the Phase constant (Beta)
#Physics: signal slows down by the sqrt(Dk)
c = 3e8 #Speed of light in  =vacuum
v_phase = c / np.sqrt(Dk)
beta = omega / v_phase

# 4. Calculate Attenuation factor (Alpha)
# A. Dielectric Loss (Scales lineaerly with frequency)
alpha_dielectric = (omega/(2 * v_phase)) * Df

#B. Conductor Loss (scales with sqrt of Frequency)
# Approximating 'k' to match typical stripline loss (~0.8 db/in at 8GHz)
alpha_conductor = 0.5 * np.sqrt(freq.f / 1e9)

#Total Alpha (Nepers/meter)
alpha_total = alpha_dielectric + alpha_conductor

#5. Create the complex Gamma 
# Gamma = Alpha + j*beta
gamma_complex = alpha_total + 1j * beta

#6. create media object
 
megtron6 = DefinedGammaZ0(freq, z0=z0_target,gamma=gamma_complex)
#Note. This is a simplified loss model. alpha controls the loss.

# 3. Create the components of the channel
# Trace 1: 6 inches
trace1 = megtron6.line(6, unit='in', name='Trace_6in')

#Trace 2: 4 inches
trace2 = megtron6.line(4, unit='in', name='Trace_4in')

# The Via (simulated as a small impedance mismatch for this demo)
# A via is often capacitive. Lets molde it as a short line with low z.
via = megtron6.line(50, unit='mil', z0=45, name='Via_Capacitive_Dip')

#4. Cascade the channels
#Chain: Trace 1 -> Via -> Trace2
channel = trace1 ** via ** trace2
# Create a copy of the channel to simulate the Lab Environment
channel_lab = channel.copy()

# # Critical Step: Connect 50-Ohm VNA cables to the 85-Ohm board
channel_lab.renormalize(50)

# --- 7. Plotting the Comparison ---
plt.figure(figsize=(10, 6))

# Plot 1: The Ideal Design (85 Ohm Reference) - What you simulated before
channel.plot_s_db(m=0, n=0, lw=2, color='orange', label='S11 (Design View - 85 Ohm Port)')

# Plot 2: The Lab Measurement (50 Ohm Reference) - THE RIPPLE
channel_lab.plot_s_db(m=0, n=0, lw=1, color='blue', alpha=0.7, linestyle='--', label='S11 (Lab View - 50 Ohm VNA)')

plt.title("Why Lab Measurements Have Ripples (50 vs 85 Ohms)")
plt.axvline(x=8e9, color='k', linestyle=':', label='Nyquist (8 GHz)')
plt.ylabel("Return Loss (dB)")
plt.xlabel("Frequency (Hz)")
plt.legend()
plt.grid(True)
plt.show()

# 5. Plot Insertion loss (S21)
plt.figure(figsize=(10,6))
channel.plot_s_db(m=1, n=0, lw=2) # S21 is port 2 to port 1
plt.title("PCIe Gen4 Channel Insertion Loss (simulated in Python)")
# Plot S11 (Return Loss) - The Reflection
# This is the line I described: it will be jagged and low
channel.plot_s_db(m=0, n=0, lw=2, label='S11 (Return Loss)')

plt.title("PCIe Gen4 Channel: S21 Vs S11 (python simulation)")
plt.axvline(x=8e9, color='k', linestyle=':', label='Nyquist (8 GHz)')
plt.ylabel("Magnitude (dB)")
plt.xlabel("Frequency (Hz)")


plt.legend()
plt.grid(True)
plt.show()

s21_at_8g = channel.s_db[np.argmin(np.abs(freq.f - 8e9)), 1, 0]
s11_at_8g = channel.s_db[np.argmin(np.abs(freq.f - 8e9)), 0, 0]

print(f"--- SIMULATION RESULTS AT 8 GHz ---")
print(f"S21 (Loss): {s21_at_8g:.2f} dB (Spec: > -28 dB) -> PASS")
print(f"S11 (Reflection): {s11_at_8g:.2f} dB (Goal: < -10 dB) -> PASS")

print(f"Loss at 8GHz: {channel.s_db[np.argmin(np.abs(freq.f - 8e9)), 1, 0]:.2f} dB")