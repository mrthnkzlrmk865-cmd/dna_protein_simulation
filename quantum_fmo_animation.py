"""
Computational Biophysics: Live Animation of Quantum Energy Transfer in FMO Complex
Author: Murathan Kizilirmak
"""

import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
from scipy.linalg import expm

# --- QUANTUM SYSTEM SETUP ---
# Hamiltonian Matrix (H) - Site energies and quantum coupling strengths
H = np.array([
    [200.0, -50.0,   0.0],  # Site 1 couples to Site 2
    [-50.0, 180.0, -80.0],  # Site 2 couples to Site 1 and Site 3
    [  0.0, -80.0, 120.0]   # Site 3 (Reaction Center Trap) couples to Site 2
])

# Initial State Vector (Psi_0): Exciton starts 100% at Site 1
Psi_0 = np.array([1.0, 0.0, 0.0], dtype=complex)
h_bar = 1.0

# Simulation Time Parameters
total_frames = 200
time_steps = np.linspace(0, 0.15, total_frames)

site_names = ["BChl Site 1 (Input)", "BChl Site 2 (Bridge)", "Reaction Center (Sink)"]
colors = ['#d62728', '#ff7f0e', '#1f77b4'] # Red, Orange, Blue

# --- PLOT & ANIMATION SETUP ---
fig, ax = plt.subplots(figsize=(9, 6))
ax.set_xlim(0, 0.15)
ax.set_ylim(-0.05, 1.05)
ax.set_title("Real-Time Quantum Coherence Dynamics (FMO Complex)", fontsize=12, fontweight='bold', pad=15)
ax.set_xlabel("Time Scale (Picoseconds)", fontsize=11)
ax.set_ylabel("Excitation Probability Population", fontsize=11)
ax.grid(True, linestyle=':', alpha=0.6)

# Initialize blank lines for the animation
lines = []
x_data = []
y_data = [[] for _ in range(3)]

for i in range(3):
    line, = ax.plot([], [], label=site_names[i], color=colors[i], linewidth=2.5 if i in [0,2] else 1.5)
    lines.append(line)

ax.legend(loc='upper right', frameon=True, shadow=True)

# --- ANIMATION CORE FUNCTIONS ---
def init():
    """Initializes the background of the animation."""
    for line in lines:
        line.set_data([], [])
    return lines

def animate(frame):
    """Calculates quantum states and updates the plot lines frame by frame."""
    t = time_steps[frame]
    
    # Solve Time-Dependent Schrödinger Equation for the current millisecond/picosecond
    evolution_operator = expm(-1j * H * t / h_bar)
    current_state = np.dot(evolution_operator, Psi_0)
    
    # Born Rule: Compute probability density P = |Psi|^2
    probabilities = np.abs(current_state) ** 2
    
    # Append current time and data point
    x_data.append(t)
    for i in range(3):
        y_data[i].append(probabilities[i])
        lines[i].set_data(x_data, y_data[i])
        
    return lines

# Create the live animation object
# interval=30 means a new frame will be rendered every 30 milliseconds
ani = animation.FuncAnimation(
    fig, animate, init_func=init, frames=total_frames, interval=30, blit=True, repeat=False
)

print("[Quantum Live Engine] Generating animated matrix graphics... Look at the active plot window!")
plt.tight_layout()
plt.show()
