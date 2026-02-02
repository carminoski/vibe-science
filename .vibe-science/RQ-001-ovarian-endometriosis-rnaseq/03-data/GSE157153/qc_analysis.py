"""
Quality Control Analysis for GSE157153.
Generates PCA, sample clustering, and UPR gene heatmap.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.stats import zscore
from pathlib import Path

# Settings
plt.style.use('seaborn-v0_8-whitegrid')
DATA_DIR = Path(__file__).parent
OUT_DIR = DATA_DIR / "figures"
OUT_DIR.mkdir(exist_ok=True)

# Load data
print("Loading data...")
expr = pd.read_csv(DATA_DIR / "GSE157153_endometriosis_associated.txt", sep="\t", index_col=0)
metadata = pd.read_csv(DATA_DIR / "GSE157153_clinical_metadata.txt", sep="\t")
upr_expr = pd.read_csv(DATA_DIR / "upr_expression_matrix.csv", index_col=0)

print(f"Expression matrix: {expr.shape}")
print(f"UPR genes: {upr_expr.shape}")
print(f"Samples: {len(metadata)}")

# Map sample names to groups
sample_to_group = dict(zip(metadata['sample_title'], metadata['group']))

# Define progression order (for visualization)
GROUP_ORDER = [
    "endometriosis",
    "atypical endometriosis",
    "adjacent endometriosis",
    "ovarian carcinoma"
]

# Color palette for groups
COLORS = {
    "endometriosis": "#2ecc71",            # Green - benign
    "atypical endometriosis": "#f1c40f",   # Yellow - atypical
    "adjacent endometriosis": "#e67e22",   # Orange - pre-malignant
    "ovarian carcinoma": "#e74c3c"         # Red - cancer
}

# ============================================================
# 1. PCA ANALYSIS
# ============================================================
print("\n" + "=" * 60)
print("1. PCA ANALYSIS")
print("=" * 60)

# Transpose for PCA (samples as rows)
X = expr.T.values
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)

# PCA
pca = PCA(n_components=10)
pca_result = pca.fit_transform(X_scaled)

# Create PCA dataframe
pca_df = pd.DataFrame({
    'PC1': pca_result[:, 0],
    'PC2': pca_result[:, 1],
    'PC3': pca_result[:, 2],
    'sample': expr.columns,
    'group': [sample_to_group.get(s, 'unknown') for s in expr.columns]
})

# PCA plot
fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# PC1 vs PC2
ax = axes[0]
for group in GROUP_ORDER:
    mask = pca_df['group'] == group
    ax.scatter(
        pca_df.loc[mask, 'PC1'],
        pca_df.loc[mask, 'PC2'],
        c=COLORS.get(group, 'gray'),
        label=group,
        alpha=0.7,
        s=80,
        edgecolors='black',
        linewidth=0.5
    )
ax.set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)')
ax.set_ylabel(f'PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)')
ax.set_title('PCA: All Genes')
ax.legend(fontsize=8)

# PC2 vs PC3
ax = axes[1]
for group in GROUP_ORDER:
    mask = pca_df['group'] == group
    ax.scatter(
        pca_df.loc[mask, 'PC2'],
        pca_df.loc[mask, 'PC3'],
        c=COLORS.get(group, 'gray'),
        label=group,
        alpha=0.7,
        s=80,
        edgecolors='black',
        linewidth=0.5
    )
ax.set_xlabel(f'PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)')
ax.set_ylabel(f'PC3 ({pca.explained_variance_ratio_[2]*100:.1f}%)')
ax.set_title('PCA: All Genes')
ax.legend(fontsize=8)

plt.tight_layout()
plt.savefig(OUT_DIR / "01_pca_all_genes.png", dpi=150, bbox_inches='tight')
plt.close()

print(f"Variance explained: PC1={pca.explained_variance_ratio_[0]*100:.1f}%, "
      f"PC2={pca.explained_variance_ratio_[1]*100:.1f}%, "
      f"PC3={pca.explained_variance_ratio_[2]*100:.1f}%")

# ============================================================
# 2. PCA ON UPR GENES ONLY
# ============================================================
print("\n" + "=" * 60)
print("2. PCA ON UPR GENES")
print("=" * 60)

X_upr = upr_expr.T.values
X_upr_scaled = scaler.fit_transform(X_upr)

pca_upr = PCA(n_components=10)
pca_upr_result = pca_upr.fit_transform(X_upr_scaled)

pca_upr_df = pd.DataFrame({
    'PC1': pca_upr_result[:, 0],
    'PC2': pca_upr_result[:, 1],
    'sample': upr_expr.columns,
    'group': [sample_to_group.get(s, 'unknown') for s in upr_expr.columns]
})

fig, ax = plt.subplots(figsize=(8, 6))
for group in GROUP_ORDER:
    mask = pca_upr_df['group'] == group
    ax.scatter(
        pca_upr_df.loc[mask, 'PC1'],
        pca_upr_df.loc[mask, 'PC2'],
        c=COLORS.get(group, 'gray'),
        label=group,
        alpha=0.7,
        s=80,
        edgecolors='black',
        linewidth=0.5
    )
ax.set_xlabel(f'PC1 ({pca_upr.explained_variance_ratio_[0]*100:.1f}%)')
ax.set_ylabel(f'PC2 ({pca_upr.explained_variance_ratio_[1]*100:.1f}%)')
ax.set_title('PCA: UPR Genes Only (N=121)')
ax.legend(fontsize=8)
plt.tight_layout()
plt.savefig(OUT_DIR / "02_pca_upr_genes.png", dpi=150, bbox_inches='tight')
plt.close()

print(f"UPR PCA variance: PC1={pca_upr.explained_variance_ratio_[0]*100:.1f}%, "
      f"PC2={pca_upr.explained_variance_ratio_[1]*100:.1f}%")

# ============================================================
# 3. HEATMAP OF KEY UPR GENES
# ============================================================
print("\n" + "=" * 60)
print("3. UPR GENE HEATMAP")
print("=" * 60)

# Key UPR genes for heatmap
key_upr = [
    # IRE1 branch
    "ERN1", "XBP1",
    # PERK branch
    "EIF2AK3", "ATF4", "DDIT3", "CHAC1", "PPP1R15A",
    # ATF6 branch
    "ATF6", "HSPA5", "CALR", "HSP90B1",
    # ERAD
    "EDEM1", "SEL1L", "SYVN1", "VCP",
    # Other markers
    "HERPUD1", "DNAJB9", "PDIA6", "HYOU1"
]

# Filter to available genes
key_upr_avail = [g for g in key_upr if g in upr_expr.index]
print(f"Key UPR genes available: {len(key_upr_avail)}/{len(key_upr)}")

# Prepare heatmap data
heatmap_data = upr_expr.loc[key_upr_avail].copy()

# Order samples by group
sample_order = []
for group in GROUP_ORDER:
    group_samples = [s for s in heatmap_data.columns if sample_to_group.get(s) == group]
    sample_order.extend(group_samples)

heatmap_data = heatmap_data[sample_order]

# Z-score normalize for visualization
heatmap_z = heatmap_data.apply(lambda x: (x - x.mean()) / x.std(), axis=1)

# Create heatmap
fig, ax = plt.subplots(figsize=(16, 10))

# Create color bar for groups
group_colors = [COLORS.get(sample_to_group.get(s, 'unknown'), 'gray') for s in sample_order]

# Plot heatmap
sns.heatmap(
    heatmap_z,
    cmap='RdBu_r',
    center=0,
    vmin=-2, vmax=2,
    xticklabels=False,
    yticklabels=True,
    ax=ax
)

# Add group color bar at top
for i, color in enumerate(group_colors):
    ax.add_patch(plt.Rectangle((i, len(key_upr_avail)), 1, 0.5, facecolor=color, edgecolor='none'))

ax.set_title('UPR Gene Expression (Z-score)\nOrdered by Progression: Endo → Atypical → Adjacent → Cancer')
ax.set_xlabel('Samples (ordered by group)')
ax.set_ylabel('UPR Genes')

plt.tight_layout()
plt.savefig(OUT_DIR / "03_upr_heatmap.png", dpi=150, bbox_inches='tight')
plt.close()

# ============================================================
# 4. BOXPLOTS OF KEY UPR GENES BY GROUP
# ============================================================
print("\n" + "=" * 60)
print("4. UPR GENE BOXPLOTS")
print("=" * 60)

# Master regulators
master_genes = ["HSPA5", "ERN1", "EIF2AK3", "ATF6"]
downstream_genes = ["XBP1", "ATF4", "DDIT3", "CHAC1"]

fig, axes = plt.subplots(2, 4, figsize=(16, 8))

for i, gene in enumerate(master_genes + downstream_genes):
    ax = axes.flatten()[i]

    data = []
    for group in GROUP_ORDER:
        group_samples = [s for s in upr_expr.columns if sample_to_group.get(s) == group]
        if gene in upr_expr.index:
            vals = upr_expr.loc[gene, group_samples].values
            for v in vals:
                data.append({'group': group, 'expression': v})

    plot_df = pd.DataFrame(data)

    sns.boxplot(data=plot_df, x='group', y='expression', palette=COLORS, ax=ax)
    ax.set_title(gene, fontweight='bold')
    ax.set_xlabel('')
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right', fontsize=8)
    ax.set_ylabel('Log2 Expression')

plt.suptitle('UPR Gene Expression Across Progression Stages', fontsize=14, fontweight='bold')
plt.tight_layout()
plt.savefig(OUT_DIR / "04_upr_boxplots.png", dpi=150, bbox_inches='tight')
plt.close()

# ============================================================
# 5. SUMMARY STATISTICS
# ============================================================
print("\n" + "=" * 60)
print("5. SUMMARY STATISTICS")
print("=" * 60)

# Calculate mean expression per group for key genes
summary = []
for gene in master_genes + downstream_genes:
    if gene in upr_expr.index:
        for group in GROUP_ORDER:
            group_samples = [s for s in upr_expr.columns if sample_to_group.get(s) == group]
            vals = upr_expr.loc[gene, group_samples]
            summary.append({
                'gene': gene,
                'group': group,
                'mean': vals.mean(),
                'std': vals.std(),
                'n': len(vals)
            })

summary_df = pd.DataFrame(summary)
summary_pivot = summary_df.pivot(index='gene', columns='group', values='mean')
summary_pivot = summary_pivot[GROUP_ORDER]  # Reorder columns
summary_pivot.to_csv(DATA_DIR / "upr_mean_expression_by_group.csv")

print("\nMean Expression by Group:")
print(summary_pivot.round(2).to_string())

# Calculate fold changes (cancer vs benign)
print("\n\nFold Changes (Cancer vs Benign Endometriosis):")
for gene in master_genes + downstream_genes:
    if gene in upr_expr.index:
        benign = [s for s in upr_expr.columns if sample_to_group.get(s) == "endometriosis"]
        cancer = [s for s in upr_expr.columns if sample_to_group.get(s) == "ovarian carcinoma"]

        mean_benign = upr_expr.loc[gene, benign].mean()
        mean_cancer = upr_expr.loc[gene, cancer].mean()

        # Log2 fold change (already log2 data, so subtract)
        log2fc = mean_cancer - mean_benign
        direction = "UP" if log2fc > 0 else "DOWN"

        print(f"  {gene}: log2FC = {log2fc:+.2f} ({direction})")

print(f"\nFigures saved to: {OUT_DIR}")
print("\nQC Analysis Complete!")
