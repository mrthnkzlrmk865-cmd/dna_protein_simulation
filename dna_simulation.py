"""
Computational Biology: DNA -> mRNA -> Protein Simulation
Author: Murathan Kizilirmak
"""

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
    'AUU': 'I', 'AUC': 'I', 'AUA': 'I', 'AUG': 'M',  # AUG is also the START codon.
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
    
    # Loop through the sequence in steps of 3 (Codons)
    for i in range(0, len(mrna) - 2, 3):
        codon = mrna[i:i+3]
        amino_acid = RNA_CODON_TABLE.get(codon, '?') # '?' for unknown codons
        
        if amino_acid == 'STOP':
            protein.append("[STOP]")
            break # STOP codon terminates translation
        else:
            protein.append(amino_acid)
            
    return "-".join(protein)

def main():
    print("=" * 50)
    print("   COMPUTATIONAL BIOLOGY: SIMULATION TOOL   ")
    print("=" * 50)
    
    # Sample DNA sequence (Feel free to change it)
    sample_dna = "ATGTACTCGGCAATCTACTTTGCAACCAATATTTAA"
    
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

if __name__ == "__main__":
    main()
