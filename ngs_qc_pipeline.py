"""
Computational Biology: NGS Data Analysis (Quality Control & Read Mapping Pipeline)
Author: murathan kizilirmak
"""

import matplotlib.pyplot as plt
import numpy as np

def phred_to_prob(ascii_char):
    """Converts a Phred-33 ASCII character to a numeric Quality Score (Q)."""
    q_score = ord(ascii_char) - 33
    # Statistical probability of incorrect base calling: P = 10^(-Q/10)
    error_prob = 10 ** (-q_score / 10)
    return q_score, error_prob

def process_fastq(fastq_data, q_threshold=28):
    """
    Parses FASTQ records, filters low-quality reads, 
    and tracks position-specific quality metrics.
    """
    print(f"\n[QC Filter] Starting Quality Control with threshold Q >= {q_threshold}...")
    
    high_quality_reads = []
    position_scores = []
    
    for read_id, sequence, quality_str in fastq_data:
        q_scores = [ord(char) - 33 for char in quality_str]
        avg_q = np.mean(q_scores)
        
        # Track quality per position for downstream analytics
        if not position_scores:
            position_scores = [[q] for q in q_scores]
        else:
            for i, q in enumerate(q_scores):
                if i < len(position_scores):
                    position_scores[i].append(q)
                    
        if avg_q >= q_threshold:
            high_quality_reads.append((read_id, sequence))
        else:
            print(f"-> Read {read_id} Rejected (Low Quality: Mean Q = {avg_q:.1f})")
            
    print(f"-> QC Complete. Passed: {len(high_quality_reads)}/{len(fastq_data)} reads.")
    return high_quality_reads, position_scores

def map_reads_to_reference(reads, reference_genome):
    """
    Simulates NGS Read Mapping (Alignment) against a reference genome template
    using a deterministic string-matching index matrix.
    """
    print("\n[Read Mapping] Aligning high-quality fragments to reference genome...")
    mapping_results = {}
    
    for read_id, sequence in reads:
        # Scan the reference genome for exact sub-string alignment matches
        position = reference_genome.find(sequence)
        if position != -1:
            mapping_results[read_id] = position
            print(f"-> Read {read_id} successfully mapped to Genome Position: {position}")
        else:
            print(f"-> Read {read_id} unmapped (No match found in reference template)")
            
    return mapping_results

def plot_quality_scores(position_scores):
    """Generates an industry-standard FastQC quality distribution chart across read lengths."""
    mean_qualities = [np.mean(positions) for positions in position_scores]
    positions = range(1, len(mean_qualities) + 1)
    
    fig, ax = plt.subplots(figsize=(8, 5))
    
    # Render background zones matching FastQC interface
    ax.axhspan(28, 40, color='#e2f0d9', alpha=0.6, label='Very Good Quality Zone')
    ax.axhspan(20, 28, color='#fff2cc', alpha=0.6, label='Warning Zone')
    ax.axhspan(0, 20, color='#fce4d6', alpha=0.6, label='Poor Quality Zone')
    
    # Plot metric analytics
    ax.plot(positions, mean_qualities, color='#1f77b4', marker='o', linewidth=2, label='Mean Q-Score')
    
    ax.set_title("NGS Sequencing Quality Distribution (FastQC Simulation)", fontsize=12, fontweight='bold')
    ax.set_xlabel("Position in Read (Nucleotide Base Index)", fontsize=10)
    ax.set_ylabel("Quality Score (Phred-33 Scale)", fontsize=10)
    ax.set_ylim(0, 42)
    ax.grid(True, linestyle='--', alpha=0.5)
    ax.legend(loc='lower left')
    
    plt.tight_layout()
    plt.show()

def main():
    print("=" * 70)
    print("   BIOINFORMATICS PIPELINE: NEXT-GENERATION SEQUENCING (NGS) ANALYSIS  ")
    print("=" * 70)
    
    # Synthetic target reference genome locus (e.g., a segment of a viral genome)
    reference_locus = "AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC"
    
    # Mock FASTQ dataset format: (Read_ID, Sequence, Phred-33 Quality String)
    # Note: 'I' = Q40 (excellent), '#' = Q2 (terrible raw error machine noise)
    mock_fastq = [
        ("Read_1", "TTTTCATTCTGA", "IIIIIIIIIIII"), # Perfect data
        ("Read_2", "CAACGGGCAATA", "IIII###IIIII"), # Fragmented noise inside read
        ("Read_3", "TGTGTGGATTAA", "IIIIIIIIIIII"), # Perfect data
        ("Read_4", "AGAGTGTCTGAT", "##II##II####")  # Failed sequencing cycle
    ]
    
    # Step 1: Quality Control Filtering & FastQC Analytics
    cleaned_reads, quality_metrics = process_fastq(mock_fastq, q_threshold=28)
    
    # Step 2: Read Mapping & Genomic Localization
    alignment_map = map_reads_to_reference(cleaned_reads, reference_locus)
    
    print(f"\nFinal Mapping Output Registry: {alignment_map}")
    
    # Step 3: FastQC Style Visualization
    print("\nGenerating quality metric reports...")
    plot_quality_scores(quality_metrics)

if __name__ == "__main__":
    main()
