"""
Computational Biology: Digital Gel Electrophoresis Simulation
Author: Murathan Kizilirmak
"""

import matplotlib.pyplot as plt
import numpy as np

def simulate_gel_electrophoresis(dna_samples):
    """
    Simulates an agarose gel electrophoresis run.
    Larger DNA fragments stay near the top (wells), 
    while smaller fragments travel further down towards the positive electrode.
    """
    # Calculate sizes (lengths) of each DNA sample in base pairs (bp)
    sample_sizes = [len(dna) for dna in dna_samples]
    
    # Setup the plot area to mimic a real Agarose Gel matrix
    fig, ax = plt.subplots(figsize=(6, 8))
    ax.set_facecolor('#1a1a1a')  # Dark background representing the UV transilluminator field
    
    num_lanes = len(dna_samples)
    
    # Plot each DNA fragment as a fluorescent band on the gel matrix
    for lane_idx, size in enumerate(sample_sizes):
        # Biophysical model: Migration distance is inversely proportional to the log of the fragment size
        migration_distance = 100 - (np.log(size) * 15)
        
        # Draw the fluorescent DNA band using specific capstyle for visual realism
        ax.hlines(y=migration_distance, xmin=lane_idx + 0.6, xmax=lane_idx + 1.4, 
                  colors='#00FF66', linewidth=6, alpha=0.9, capstyle='round')
        
        # Display the fragment size metrics next to each distinct band
        ax.text(lane_idx + 1.45, migration_distance, f"{size} bp", 
                color='white', fontsize=9, va='center')

    # Render the sample loading wells at the negative terminal (top)
    for lane_idx in range(num_lanes):
        rect = plt.Rectangle((lane_idx + 0.7, 95), 0.6, 3, facecolor='#333333', edgecolor='#555555')
        ax.add_patch(rect)
        ax.text(lane_idx + 1.0, 92, f"Lane {lane_idx + 1}", color='#888888', 
                fontsize=10, ha='center', fontweight='bold')

    # Configure scientific chart attributes and multi-lane laboratory formatting
    ax.set_title("Agarose Gel Electrophoresis Simulation", color='white', fontsize=14, fontweight='bold', pad=20)
    ax.set_xlim(0, num_lanes + 1)
    ax.set_ylim(0, 100)
    
    # Y-axis indicates molecular migration velocity vector from Negative (-) to Positive (+)
    ax.set_ylabel("DNA Migration Direction (—) ---------> (+)", color='white', fontsize=12)
    ax.set_xticks([]) 
    ax.set_yticks([]) 
    
    # Formal boundaries for the simulated gel framework
    for spine in ax.spines.values():
        spine.set_color('#444444')
        spine.set_linewidth(2)

    plt.tight_layout()
    plt.show()

def main():
    print("=" * 50)
    print("   COMPUTATIONAL BIOLOGY: GEL ELECTROPHORESIS   ")
    print("=" * 50)
    
    # Generation of 4 distinct DNA matrices simulating diverse nucleotide molecular weights
    sample_1 = "ATGC" * 5    # 20 bp molecular weight
    sample_2 = "ATGC" * 50   # 200 bp molecular weight
    sample_3 = "ATGC" * 150  # 600 bp molecular weight
    sample_4 = "ATGC" * 300  # 1200 bp molecular weight
    
    samples = [sample_1, sample_2, sample_3, sample_4]
    
    print("Analyzing DNA sample sizes...")
    for i, sample in enumerate(samples):
        print(f"Lane {i+1} Molecular Weight: {len(sample)} base pairs (bp)")
        
    print("\nExecuting digital gel electrophoresis framework...")
    simulate_gel_electrophoresis(samples)

if __name__ == "__main__":
    main()
