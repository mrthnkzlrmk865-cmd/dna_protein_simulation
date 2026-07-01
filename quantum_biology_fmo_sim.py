"""
Computational Biophysics: Quantum Coherence & Energy Transfer in FMO Complex
Author: Murathan Kizilirmak
"""

import matplotlib.pyplot as plt
import numpy as np
from scipy.linalg import expm

def simulate_quantum_dynamics(hamiltonian, initial_state, time_steps, h_bar=1.0):
    """
    Simulates the time evolution of a quantum state using the Time-Dependent Schrödinger Equation.
    Psi(t) = exp(-i * H * t / h_bar) * Psi(0)
    """
    print("\n[Quantum Engine] Solving Schrödinger time-evolution matrix...")
    
    num_sites = len(initial_state)
    probabilities_over_time = []
    
    for t in time_steps:
        # Compute the unitary time-evolution operator: U(t) = exp(-i * H * t)
        # In quantum mechanics, i is represented by 1j in Python
        evolution_operator = expm(-1j * hamiltonian * t / h_bar)
        
        # Evolve the wave function state vector
        current_state = np.dot(evolution_operator, initial_state)
        
        # Calculate Probability Density: P(i) = |Psi(i)|^2
        probabilities = np.abs(current_state) ** 2
        probabilities_over_time.append(probabilities)
        
    return np.array(probabilities_over_time)

def plot_quantum_transfer(time_steps, probabilities, site_names):
    """Plots the population dynamics of quantum excitation transfer over time."""
    fig, ax = plt.subplots(figsize=(9, 6))
    
    for i, site in enumerate(site_names):
        # Line width settings to highlight the initial and final target sites
        linewidth = 2.5 if i in [0, 2] else 1.5
        ax.plot(time_steps, probabilities[:, i], label=site, linewidth=linewidth)
        
    ax.set_title("Quantum Coherence Exciton Energy Transfer (FMO Complex)", fontsize=12, fontweight='bold', pad=15)
    ax.set_xlabel("Time Scale (Picoseconds / Arbitrary Units)", fontsize=11)
    ax.set_ylabel("Excitation Probability Population", fontsize=11)
    ax.set_ylim(-0.05, 1.05)
    ax.grid(True, linestyle=':', alpha=0.6)
    ax.legend(loc='upper right', frameon=True, shadow=True)
    
    plt.tight_layout()
    plt.show()

def main():
    print("=" * 75)
    print("   BIOINPHYSICS SIMULATION: QUANTUM BIOLOGY LIGHT-HARVESTING ECOSYSTEM   ")
    print("=" * 75)
    
    # Simulating a simplified 3-site network within the FMO Complex:
    # Site 1: Bacteriochlorophyll where the light energy is absorbed (Initial state)
    # Site 2: Intermediate bridging chlorophyll molecule
    # Site 3: Target Reaction Center interface (Final destination)
    site_names = ["BChl Site 1 (Light Input)", "BChl Site 2 (Bridge)", "Reaction Center (Sink)"]
    
    # Define the Quantum Hamiltonian Matrix (H) in arbitrary energy units.
    # Diagonal elements: Site energies
    # Off-diagonal elements: Quantum coupling interaction coefficients between sites
    H = np.array([
        [200.0, -50.0,   0.0],  # Site 1 couples with Site 2
        [-50.0, 180.0, -80.0],  # Site 2 couples with both Site 1 and Site 3
        [  0.0, -80.0, 120.0]   # Site 3 couples with Site 2 (Lower energy to trap exciton)
    ])
    
    # Initial State Vector (Psi_0): Exciton starts 100% at Site 1
    # Complex state vector format: [1.0 + 0j, 0.0 + 0j, 0.0 + 0j]
    Psi_0 = np.array([1.0, 0.0, 0.0], dtype=complex)
    
    # Time vector for simulation tracking (0 to 0.15 picoseconds equivalent scale)
    t_space = np.linspace(0, 0.15, 300)
    
    # Execute biophysical quantum transport equation
    prob_matrix = simulate_quantum_dynamics(H, Psi_0, t_space)
    
    print("\nSimulation successful. Rendering quantum trajectory waveforms...")
    plot_quantum_transfer(t_space, prob_matrix, site_names)

if __name__ == "__main__":
    main()
