"""
Computational Biology: Differential Gene Expression (DEG) & Volcano Plot Analytics
Author: Murathan Kizilirmak
"""

import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

def perform_deg_analysis(normal_matrix, tumor_matrix, gene_names):
    """
    Performs statistical differential expression analysis.
    Computes Log2 Fold Change and Student's t-test p-values for each gene.
    """
    print("\n[Statistical Engine] Running Student's t-test across transcriptomic matrices...")
    
    results = []
    for i, gene in enumerate(gene_names):
        normal_samples = normal_matrix[i]
        tumor_samples = tumor_matrix[i]
        
        # Calculate mean expression levels
        mean_normal = np.mean(normal_samples)
        mean_tumor = np.mean(tumor_samples)
        
        # Avoid division by zero or log issues by ensuring means are non-zero
        mean_normal = max(mean_normal, 0.1)
        mean_tumor = max(mean_tumor, 0.1)
        
        # Calculate Log2 Fold Change: log2(Tumor / Normal)
        log2_fc = np.log2(mean_tumor / mean_normal)
        
        # Perform Independent Student's t-test to check statistical significance
        t_stat, p_value = stats.ttest_ind(tumor_samples, normal_samples, equal_var=False)
        
        # Calculate -Log10(p-value) for Volcano Plot scaling
        p_value = max(p_value, 1e-10)
        neg_log10_p = -np.log10(p_value)
        
        results.append({
            "gene": gene,
            "log2FC": log2_fc,
            "p_value": p_value,
            "neg_log10_p": neg_log10_p
        })
        
    return results

def plot_volcano(deg_results):
    """Generates an industry-standard Volcano Plot marking significant biomolecules."""
    log2_fcs = [r["log2FC"] for r in deg_results]
    neg_log_ps = [r["neg_log10_p"] for r in deg_results]
    genes = [r["gene"] for r in deg_results]
    
    fig, ax = plt.subplots(figsize=(8, 7))
    
    # Statistical significance thresholds
    fc_threshold = 1.0  # Represents a 2-fold change line
    p_threshold_line = -np.log10(0.05)  # Significance cut-off line (p = 0.05)
    
    # Classify genes based on thresholds for conditional coloring
    for r in deg_results:
        x, y = r["log2FC"], r["neg_log10_p"]
        
        if y >= p_threshold_line and x >= fc_threshold:
            # Significant up-regulated genes (Oncogenes)
            ax.scatter(x, y, color='#d62728', s=50, edgecolors='black', alpha=0.8)
            ax.text(x + 0.05, y, r["gene"], fontsize=8, fontweight='bold', color='#d62728')
        elif y >= p_threshold_line and x <= -fc_threshold:
            # Significant down-regulated genes (Tumor Suppressors)
            ax.scatter(x, y, color='#1f77b4', s=50, edgecolors='black', alpha=0.8)
            ax.text(x - 0.25, y, r["gene"], fontsize=8, fontweight='bold', color='#1f77b4')
        else:
            # Non-significant genes (Background noise)
            ax.scatter(x, y, color='#7f7f7f', s=30, alpha=0.4)

    # Draw operational threshold divider lines
    ax.axhline(y=p_threshold_line, color='black', linestyle='--', linewidth=1.2, alpha=0.7)
    ax.axvline(x=fc_threshold, color='black', linestyle='--', linewidth=1.2, alpha=0.7)
    ax.axvline(x=-fc_threshold, color='black', linestyle='--', linewidth=1.2, alpha=0.7)
    
    # Labels and Titles
    ax.set_title("Differential Gene Expression Analysis (Volcano Plot Profile)", fontsize=12, fontweight='bold', pad=15)
    ax.set_xlabel("Log2 Fold Change (Tumor vs. Normal)", fontsize=11)
    ax.set_ylabel("-Log10 (p-value) Significance Index", fontsize=11)
    
    # Fixed the text bug here (Removed the invalid linestyle parameter)
    ax.text(ax.get_xlim()[0] + 0.2, p_threshold_line + 0.1, "p = 0.05 cutoff", color='black', fontsize=9)
    
    ax.grid(True, linestyle=':', alpha=0.5)
    plt.tight_layout()
    plt.show()

def main():
    print("=" * 75)
    print("   BIOINFORMATICS PIPELINE: STATISTICAL DEG ANALYZER & VOLCANO ENGINE   ")
    print("=" * 75)
    
    # Target Biomarker Genes names
    gene_panel = ["TP53", "BRCA1", "MYC", "EGFR", "PTEN", "MDM2", "VEGFA", "AKT1", "GAPDH", "ACTB"]
    
    # Simulating 5 healthy samples and 5 tumor samples for each gene
    np.random.seed(42) 
    
    normal_matrix = [
        np.random.normal(50, 5, 5),   # TP53
        np.random.normal(40, 4, 5),   # BRCA1
        np.random.normal(15, 3, 5),   # MYC
        np.random.normal(20, 3, 5),   # EGFR
        np.random.normal(60, 6, 5),   # PTEN
        np.random.normal(10, 2, 5),   # MDM2
        np.random.normal(25, 4, 5),   # VEGFA
        np.random.normal(30, 4, 5),   # AKT1
        np.random.normal(100, 5, 5),  # GAPDH
        np.random.normal(120, 6, 5)   # ACTB
    ]
    
    tumor_matrix = [
        np.random.normal(12, 2, 5),   # TP53
        np.random.normal(42, 5, 5),   # BRCA1
        np.random.normal(95, 8, 5),   # MYC
        np.random.normal(80, 7, 5),   # EGFR
        np.random.normal(15, 2, 5),   # PTEN
        np.random.normal(45, 5, 5),   # MDM2
        np.random.normal(70, 6, 5),   # VEGFA
        np.random.normal(35, 4, 5),   # AKT1
        np.random.normal(102, 5, 5),  # GAPDH
        np.random.normal(118, 6, 5)   # ACTB
    ]
    
    deg_report = perform_deg_analysis(normal_matrix, tumor_matrix, gene_panel)
    
    print("\nStatistical Analysis Complete. Rendering Volcano matrix plot...")
    plot_volcano(deg_report)

if __name__ == "__main__":
    main()
