import matplotlib.pyplot as plt
import numpy as np
from qiskit import QuantumCircuit
from qiskit_aer import AerSimulator

def create_ansatz_circuit(theta):
    qc = QuantumCircuit(2)
    qc.x(0) 
    qc.ry(theta, 1)
    qc.cx(1, 0)
    return qc

def compute_expectation_value(theta, hamiltonian_weights):
    simulator = AerSimulator()
    
    qc_z = create_ansatz_circuit(theta)
    qc_z.measure_all()
    
    shots = 2048
    result = simulator.run(qc_z, shots=shots).result()
    counts = result.get_counts()
    
    p00 = counts.get('00', 0) / shots
    p01 = counts.get('01', 0) / shots
    p10 = counts.get('10', 0) / shots
    p11 = counts.get('11', 0) / shots
    
    exp_Z0 = (p00 + p10) - (p01 + p11)
    exp_Z1 = (p00 + p01) - (p10 + p11)
    
    qc_x = create_ansatz_circuit(theta)
    qc_x.h(0)
    qc_x.h(1)
    qc_x.measure_all()
    
    result_x = simulator.run(qc_x, shots=shots).result()
    counts_x = result_x.get_counts()
    
    px00 = counts_x.get('00', 0) / shots
    px01 = counts_x.get('01', 0) / shots
    px10 = counts_x.get('10', 0) / shots
    px11 = counts_x.get('11', 0) / shots
    
    exp_X0X1 = (px00 + px11) - (px01 + px10)
    
    c0, c1, c2, c3 = hamiltonian_weights
    total_energy = c0 + (c1 * exp_Z0) + (c2 * exp_Z1) + (c3 * exp_X0X1)
    
    return total_energy

def main():
    print("=" * 75)
    print("   QUANTUM COMPUTING: MOLECULAR HAMILTONIAN SIMULATION (QISKIT ENGINE)   ")
    print("=" * 75)
    
    distances = np.linspace(0.3, 2.5, 20)
    ground_state_energies = []
    
    print("[Quantum VQE] Mapping molecular potential energy surface (PES)...")
    
    for d in distances:
        c0 = -1.05 + (0.4 / d)
        c1 = 0.40 - (0.1 / d)
        c2 = 0.40 - (0.1 / d)
        c3 = 0.20 + (0.05 * d)
        weights = [c0, c1, c2, c3]
        
        angles = np.linspace(0, 2 * np.pi, 50)
        energies_for_distance = [compute_expectation_value(theta, weights) for theta in angles]
        
        min_energy = min(energies_for_distance)
        ground_state_energies.append(min_energy)
        
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.plot(distances, ground_state_energies, 'o-', color='#8c564b', linewidth=2, label="H2 Ground State (VQE)")
    
    optimal_idx = np.argmin(ground_state_energies)
    ax.scatter(distances[optimal_idx], ground_state_energies[optimal_idx], color='#d62728', s=100, zorder=5, 
               label=f"Stable Bond Length (~{distances[optimal_idx]:.2f} Å)")
    
    ax.set_title("H2 Molecule Potential Energy Surface via Quantum Simulation", fontsize=12, fontweight='bold', pad=15)
    ax.set_xlabel("Interatomic Distance (Ångström)", fontsize=11)
    ax.set_ylabel("Total Molecular Energy (Hartree)", fontsize=11)
    ax.grid(True, linestyle=':', alpha=0.6)
    ax.legend(frameon=True, shadow=True)
    
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    main()
