"""
Computational Biology: DNA -> mRNA -> Protein Simulation with Visualization
Author: Murathan Kizilirmak
"""

import matplotlib.pyplot as plt

# Genetic Code Dictionary (RNA Codon Table)
RNA_CODON_TABLE = {
    'UUU': 'F', 'UUC': 'F', 'UUA': 'L', 'UUG': 'L',
    'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S',
    'UAU': 'Y', 'UAC': 'Y', 'UAA': 'STOP', 'UAG': 'STOP',
    'UGU': 'C', 'UGC': 'C', 'UGA': 'STOP', 'UGG': 'W',
    'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L',
    'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'CAU': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AUU': 'I', 'AUC': 'I', 'AUA': 'I', 'AUG': 'M',
    'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'AAU': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGU': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GUU': 'V', 'GUC': 'V', 'GUA': 'V', 'GUG': 'V',
    'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'GAU': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGU': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
}

def transcribe(dna_sequence):
    """Converts a DNA sequence into an mRNA sequence (T -> U)."""
    return dna_sequence.upper().replace('T', 'U')

def calculate_gc_content(dna_sequence):
    """Calculates the percentage of Guanine (G) and Cytosine (C) in the sequence."""
    dna = dna_sequence.upper()
    g_count = dna.count('G')
    c_count = dna.count('C')
    
    if len(dna) == 0:
        return 0.0
    
    return ((g_count + c_count) / len(dna)) * 100

def translate_mrna(mrna_sequence):
    """Translates mRNA sequence into an amino acid chain (protein) by codons."""
    mrna = mrna_sequence.upper()
    protein = []
    
    for i in range(0, len(mrna) - 2, 3):
        codon = mrna[i:i+3]
        amino_acid = RNA_CODON_TABLE.get(codon, '?')
        
        if amino_acid == 'STOP':
            protein.append("[STOP]")
            break
        else:
            protein.append(amino_acid)
            
    return "-".join(protein)

def plot_base_frequencies(dna_sequence):
    """Counts DNA bases and generates a bar chart visualization."""
    dna = dna_sequence.upper()
    
    # Count the occurrences of each base
    frequencies = {
        'Adenine (A)': dna.count('A'),
        'Thymine (T)': dna.count('T'),
        'Guanine (G)': dna.count('G'),
        'Cytosine (C)': dna.count('C')
    }
    
    # Graph Styling
    colors = ['#4CAF50', '#FF5722', '#2196F3', '#9C27B0'] # Custom colors for DNA bases
    plt.figure(figsize=(8, 5))
    
    # Create Bar Chart
    plt.bar(frequencies.keys(), frequencies.values(), color=colors, edgecolor='black')
    
    # Add Titles and Labels
    plt.title("DNA Base Frequency Analysis", fontsize=14, fontweight='bold')
    plt.xlabel("Nucleotide Bases", fontsize=12)
    plt.ylabel("Count / Frequency", fontsize=12)
    plt.grid(axis='y', linestyle='--', alpha=0.7) # Gridlines for better readability
    
    # Show the plot
    plt.tight_layout()
    plt.show()

def main():
    print("=" * 50)
    print("   COMPUTATIONAL BIOLOGY: SIMULATION TOOL   ")
    print("=" * 50)
    
    # A longer sample DNA sequence to make the graph look better
    sample_dna = "ATGGCCATTGTAATGGGCCGCTGAAAGGGTCCCAATTTTTAA"
    
    print(f"Input DNA Sequence : {sample_dna}")
    
    # 1. GC Content Analysis
    gc_percent = calculate_gc_content(sample_dna)
    print(f"GC Content Ratio   : {gc_percent:.2f}%")
    
    # 2. Transcription (DNA -> mRNA)
    mrna = transcribe(sample_dna)
    print(f"mRNA Sequence      : {mrna}")
    
    # 3. Translation (mRNA -> Protein)
    protein_sequence = translate_mrna(mrna)
    print(f"Synthesized Protein: {protein_sequence}")
    print("=" * 50)
    
    # 4. Data Visualization
    print("Generating visualization window...")
    plot_base_frequencies(sample_dna)

if __name__ == "__main__":
    main()
