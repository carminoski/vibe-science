"""
Differential Expression Analysis for UPR genes in GSE157153.
Statistical testing with multiple comparison correction.
"""

import pandas as pd
import numpy as np
from scipy import stats
from statsmodels.stats.multitest import multipletests
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

DATA_DIR = Path(__file__).parent

# Load data
print("Loading data...")
expr = pd.read_csv(DATA_DIR / "GSE157153_endometriosis_associated.txt", sep="\t", index_col=0)
metadata = pd.read_csv(DATA_DIR / "GSE157153_clinical_metadata.txt", sep="\t")
upr_expr = pd.read_csv(DATA_DIR / "upr_expression_matrix.csv", index_col=0)

# Map sample names to groups
sample_to_group = dict(zip(metadata['sample_title'], metadata['group']))

# Define comparisons
COMPARISONS = [
    ("ovarian carcinoma", "endometriosis"),           # Cancer vs Benign
    ("atypical endometriosis", "endometriosis"),      # Atypical vs Benign
    ("adjacent endometriosis", "endometriosis"),      # Adjacent vs Benign
    ("adjacent endometriosis", "atypical endometriosis"),  # Adjacent vs Atypical
    ("ovarian carcinoma", "adjacent endometriosis"),  # Cancer vs Adjacent
    ("ovarian carcinoma", "atypical endometriosis"),  # Cancer vs Atypical
]

def run_de_analysis(expr_df, group1_samples, group2_samples, group1_name, group2_name):
    """Run differential expression analysis between two groups."""
    results = []

    for gene in expr_df.index:
        g1_vals = expr_df.loc[gene, group1_samples].values
        g2_vals = expr_df.loc[gene, group2_samples].values

        # Calculate statistics
        mean1 = np.mean(g1_vals)
        mean2 = np.mean(g2_vals)
        log2fc = mean1 - mean2  # Already log2 data

        # t-test
        try:
            stat, pval = stats.ttest_ind(g1_vals, g2_vals)
        except:
            pval = 1.0

        results.append({
            'gene': gene,
            f'mean_{group1_name}': mean1,
            f'mean_{group2_name}': mean2,
            'log2FC': log2fc,
            'pvalue': pval
        })

    # Create dataframe
    results_df = pd.DataFrame(results)

    # Multiple testing correction
    _, padj, _, _ = multipletests(results_df['pvalue'].fillna(1), method='fdr_bh')
    results_df['padj'] = padj

    # Sort by p-value
    results_df = results_df.sort_values('pvalue')

    return results_df

# ============================================================
# RUN ALL COMPARISONS
# ============================================================
print("\n" + "=" * 70)
print("DIFFERENTIAL EXPRESSION ANALYSIS - UPR GENES")
print("=" * 70)

all_results = {}

for group1, group2 in COMPARISONS:
    print(f"\n--- {group1} vs {group2} ---")

    g1_samples = [s for s in upr_expr.columns if sample_to_group.get(s) == group1]
    g2_samples = [s for s in upr_expr.columns if sample_to_group.get(s) == group2]

    print(f"  Samples: {len(g1_samples)} vs {len(g2_samples)}")

    de_results = run_de_analysis(upr_expr, g1_samples, g2_samples, group1, group2)

    # Count significant genes
    sig_up = de_results[(de_results['padj'] < 0.05) & (de_results['log2FC'] > 0.5)].shape[0]
    sig_down = de_results[(de_results['padj'] < 0.05) & (de_results['log2FC'] < -0.5)].shape[0]
    print(f"  Significant (padj<0.05, |log2FC|>0.5): {sig_up} UP, {sig_down} DOWN")

    # Save
    comparison_name = f"{group1.replace(' ', '_')}_vs_{group2.replace(' ', '_')}"
    de_results.to_csv(DATA_DIR / f"DE_{comparison_name}.csv", index=False)
    all_results[comparison_name] = de_results

    # Show top genes
    print(f"\n  Top 5 UP-regulated (padj<0.05):")
    top_up = de_results[(de_results['padj'] < 0.05) & (de_results['log2FC'] > 0)].head(5)
    for _, row in top_up.iterrows():
        print(f"    {row['gene']}: log2FC={row['log2FC']:.2f}, padj={row['padj']:.2e}")

    print(f"\n  Top 5 DOWN-regulated (padj<0.05):")
    top_down = de_results[(de_results['padj'] < 0.05) & (de_results['log2FC'] < 0)].head(5)
    for _, row in top_down.iterrows():
        print(f"    {row['gene']}: log2FC={row['log2FC']:.2f}, padj={row['padj']:.2e}")

# ============================================================
# KEY COMPARISON: CANCER vs BENIGN ENDOMETRIOSIS
# ============================================================
print("\n" + "=" * 70)
print("FOCUS: OVARIAN CARCINOMA vs BENIGN ENDOMETRIOSIS")
print("=" * 70)

key_result = all_results["ovarian_carcinoma_vs_endometriosis"]

# Master UPR genes
master_genes = ["HSPA5", "ERN1", "EIF2AK3", "ATF6", "XBP1", "ATF4", "DDIT3", "CHAC1"]

print("\nMaster UPR Genes:")
print("-" * 60)
print(f"{'Gene':<12} {'log2FC':>10} {'pvalue':>12} {'padj':>12} {'Sig?':<6}")
print("-" * 60)

for gene in master_genes:
    row = key_result[key_result['gene'] == gene]
    if len(row) > 0:
        row = row.iloc[0]
        sig = "***" if row['padj'] < 0.001 else ("**" if row['padj'] < 0.01 else ("*" if row['padj'] < 0.05 else ""))
        print(f"{gene:<12} {row['log2FC']:>+10.2f} {row['pvalue']:>12.2e} {row['padj']:>12.2e} {sig:<6}")

# ============================================================
# PROGRESSION ANALYSIS
# ============================================================
print("\n" + "=" * 70)
print("PROGRESSION ANALYSIS: Trend across stages")
print("=" * 70)

# Define progression order
progression_order = ["endometriosis", "atypical endometriosis", "adjacent endometriosis", "ovarian carcinoma"]

# Calculate mean expression per stage for master genes
print("\nMean expression across progression stages:")
print("-" * 75)
header = f"{'Gene':<12}"
for stage in progression_order:
    header += f" {stage[:12]:>12}"
print(header)
print("-" * 75)

progression_data = []
for gene in master_genes:
    if gene in upr_expr.index:
        row_data = {'gene': gene}
        means = []
        for stage in progression_order:
            stage_samples = [s for s in upr_expr.columns if sample_to_group.get(s) == stage]
            mean_val = upr_expr.loc[gene, stage_samples].mean()
            row_data[stage] = mean_val
            means.append(mean_val)

        # Test for monotonic trend (Spearman correlation with stage order)
        rho, pval = stats.spearmanr([0, 1, 2, 3], means)
        row_data['trend_rho'] = rho
        row_data['trend_pval'] = pval
        progression_data.append(row_data)

        # Print
        row_str = f"{gene:<12}"
        for mean_val in means:
            row_str += f" {mean_val:>12.2f}"
        trend_dir = "↑" if rho > 0.5 else ("↓" if rho < -0.5 else "→")
        row_str += f"  {trend_dir} (ρ={rho:.2f}, p={pval:.3f})"
        print(row_str)

progression_df = pd.DataFrame(progression_data)
progression_df.to_csv(DATA_DIR / "UPR_progression_analysis.csv", index=False)

# ============================================================
# SUMMARY
# ============================================================
print("\n" + "=" * 70)
print("SUMMARY")
print("=" * 70)

# Count significant UPR genes in main comparison
main_de = all_results["ovarian_carcinoma_vs_endometriosis"]
sig_genes = main_de[main_de['padj'] < 0.05]
up_genes = sig_genes[sig_genes['log2FC'] > 0.5]['gene'].tolist()
down_genes = sig_genes[sig_genes['log2FC'] < -0.5]['gene'].tolist()

print(f"\nCancer vs Benign Endometriosis:")
print(f"  Total UPR genes tested: {len(main_de)}")
print(f"  Significant (padj<0.05): {len(sig_genes)}")
print(f"  Upregulated (log2FC>0.5): {len(up_genes)}")
print(f"  Downregulated (log2FC<-0.5): {len(down_genes)}")

if up_genes:
    print(f"\n  Top upregulated UPR genes:")
    for gene in up_genes[:10]:
        row = main_de[main_de['gene'] == gene].iloc[0]
        print(f"    {gene}: log2FC={row['log2FC']:.2f}, padj={row['padj']:.2e}")

print("\nAnalysis complete! Results saved.")
