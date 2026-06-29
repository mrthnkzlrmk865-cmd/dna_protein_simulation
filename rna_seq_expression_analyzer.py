"""
Computational Biology: RNA-Seq Gene Expression Profiles & RPKM Normalization
Author: Murathan Kizilirmak
"""

import matplotlib.pyplot as plt
import numpy as np

def calculate_rpkm(raw_counts, gen_lengths, total_reads):
    """
    Normalizes raw RNA-Seq read counts using the RPKM method.
    Formula: RPKM = (C * 10^9) / (L * N)
    Where:
    C = Number of reads mapped to a gene
    L = Gene length in base pairs (bp)
    N = Total mapped reads in the entire sequencing run
    """
    rpkm_values = {}
    for gene, count in raw_counts.items():
        length = gen_lengths[gene]
        # Mathematical RPKM matrix calculation
        rpkm = (count * 1e9) / (length * total_reads)
        rpkm_values[gene] = round(rpkm, 2)
    return rpkm_values

def analyze_differential_expression(normal_rpkm, tumor_rpkm, threshold=2.0):
    """
    Identifies up-regulated or down-regulated genes by computing the fold change.
    Fold Change = Tumor RPKM / Normal RPKM
    """
    print(f"\n[Differential Analysis] Scanning genes for Expression Fold Change >= {threshold}x...")
    analysis_report = {}
    
    for gene in normal_rpkm:
        norm_val = normal_rpkm[gene] if normal_rpkm[gene] > 0 else 0.1 # Avoid division by zero
        tumor_val = tumor_rpkm[gene]
        
        fold_change = tumor_val / norm_val
        analysis_report[gene] = fold_change
        
        if fold_change >= threshold:
            print(f"-> ALERT: {gene} is UP-REGULATED in tumor tissue! (Fold Change = {fold_change:.2f}x)")
        elif fold_change <= (1 / threshold):
            print(f"-> NOTE: {gene} is DOWN-REGULATED in tumor tissue! (Fold Change = {fold_change:.2f}x)")
            
    return analysis_report

def plot_expression_profiles(genes, normal_rpkm, tumor_rpkm):
    """Generates a comparative side-by-side bar chart of normalized gene expressions."""
    normal_values = [normal_rpkm[g] for g in genes]
    tumor_values = [tumor_rpkm[g] for g in genes]
    
    x = np.arange(len(genes))
    width = 0.35  # Width of the bars
    
    fig, ax = plt.subplots(figsize=(9, 6))
    
    # Render paired bar charts
    rects1 = ax.bar(x - width/2, normal_values, width, label='Healthy Control (Normal)', color='#2ca02c', alpha=0.85)
    rects2 = ax.bar(x + width/2, tumor_values, width, label='Patient Biopsy (Tumor)', color='#d62728', alpha=0.85)
    
    ax.set_title("RNA-Seq Comparative Transcriptome Analysis (Normalized RPKM Metrics)", fontsize=12, fontweight='bold', pad=15)
    ax.set_xlabel("Target Biomarker Genes", fontsize=11)
    ax.set_ylabel("Gene Expression Level (RPKM Score)", fontsize=11)
    ax.set_xticks(x)
    ax.set_xticklabels(genes, fontweight='bold')
    ax.legend()
    ax.grid(axis='y', linestyle='--', alpha=0.5)
    
    plt.tight_layout()
    plt.show()

def main():
    print("=" * 75)
    print("   BIOINFORMATICS PIPELINE: RNA-SEQ TRANSCRIPTOME EXPRESSION ANALYZER   ")
    print("=" * 75)
    
    # Target cancer-related panel genes and their genomic lengths in base pairs (bp)
    gene_lengths = {
        "TP53": 1200,   # Tumor suppressor
        "BRCA1": 5500,  # DNA repair gene (Very long)
        "MYC": 1400,    # Oncogene (Cancer promoter)
        "EGFR": 3800    # Growth factor receptor
    }
    genes = list(gene_lengths.keys())
    
    # Raw Read Counts received directly from the sequencing machine instruments
    # Notice how BRCA1 has high raw counts in both, partly because it is a massive gene
    raw_counts_normal = {"TP53": 450, "BRCA1": 1200, "MYC": 150, "EGFR": 300}
    raw_counts_tumor  = {"TP53": 110, "BRCA1": 1300, "MYC": 980, "EGFR": 1450}
    
    total_reads_normal = sum(raw_counts_normal.values())
    total_reads_tumor = sum(raw_counts_tumor.values())
    
    # Step 1: Execute RPKM Mathematical Normalization Matrix
    print("Executing transcriptomic depth and length normalization...")
    normal_rpkm = calculate_rpkm(raw_counts_normal, gene_lengths, total_reads_normal)
    tumor_rpkm = calculate_rpkm(raw_counts_tumor, gene_lengths, total_reads_tumor)
    
    print(f"-> Healthy Tissue Expression Registry (RPKM): {normal_rpkm}")
    print(f"-> Tumor Tissue Expression Registry (RPKM): {tumor_rpkm}")
    
    # Step 2: Run Differential Expression Analytics
    analyze_differential_expression(normal_rpkm, tumor_rpkm, threshold=2.5)
    
    # Step 3: Plot comparative figures
    print("\nRendering comparative transcriptional profiles...")
    plot_expression_profiles(genes, normal_rpkm, tumor_rpkm)

if __name__ == "__main__":
    main()
