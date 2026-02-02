"""
Check UPR gene coverage in GSE157153 expression data.
"""

import pandas as pd
from pathlib import Path

DATA_DIR = Path(__file__).parent

# UPR genes (Hallmark + Core)
HALLMARK_UPR = """ALDH18A1 ARFGAP1 ASNS ATF3 ATF4 ATF6 ATP6V0D1 BAG3 BANF1 CALR
CCL2 CEBPB CEBPG CHAC1 CKS1B CNOT2 CNOT4 CNOT6 CXXC1 DCP1A
DCP2 DCTN1 DDIT4 DDX10 DKC1 DNAJA4 DNAJB9 DNAJC3 EDC4 EDEM1
EEF2 EIF2AK3 EIF2S1 EIF4A1 EIF4A2 EIF4A3 EIF4E EIF4EBP1 EIF4G1
ERN1 ERO1A EXOC2 EXOSC1 EXOSC10 EXOSC2 EXOSC4 EXOSC5 EXOSC9
FKBP14 FUS GEMIN4 GOSR2 H2AX HERPUD1 HSP90B1 HSPA5 HSPA9 HYOU1
IARS1 IFIT1 IGFBP1 IMP3 KDELR3 KHSRP KIF5B LSM1 LSM4 MTHFD2
NFYA NFYB NHP2 NOLC1 NOP14 NOP56 NPM1 NABP1 PAIP1 PARN PDIA5
PDIA6 POP4 PREB PSAT1 RPS14 RRP9 SDAD1 SEC11A SEC31A SERP1
SHC1 MTREX SLC1A4 SLC30A5 SLC7A5 SPCS1 SPCS3 SRPRA SRPRB SSR1
STC2 TARS1 TATDN2 TSPYL2 SKIC3 TUBB2A VEGFA WFS1 WIPI1 XBP1
XPOT YIF1A YWHAZ ZBTB17""".split()

CORE_UPR = """HSPA5 ERN1 EIF2AK3 ATF6 XBP1 ATF4 DDIT3 CHAC1 PPP1R15A
CALR CANX PDIA3 PDIA6 HSP90B1 HYOU1 EDEM1 SEL1L SYVN1
DERL1 VCP BAX BBC3 TRIB3 HERPUD1 DNAJB9 DNAJC3""".split()

# Lab panel genes (from user)
LAB_PANEL = ["HSPA5", "XBP1", "STAT3", "IL6", "CXCL8", "LIFR", "RELA"]  # GRP78, XBP1, STAT3, IL-6, IL-8, LIFR, p65

# Aliases
ALIASES = {
    "GRP78": "HSPA5",
    "BiP": "HSPA5",
    "CHOP": "DDIT3",
    "GADD34": "PPP1R15A",
    "IL8": "CXCL8",
    "IL-8": "CXCL8",
    "IL-6": "IL6",
    "p65": "RELA",
    "GRP94": "HSP90B1",
    "H2AFX": "H2AX",  # alias
    "IARS": "IARS1",  # alias
    "TARS": "TARS1",  # alias
}

# Load expression data
print("Loading expression data...")
expr = pd.read_csv(DATA_DIR / "GSE157153_endometriosis_associated.txt", sep="\t", index_col=0)
print(f"Expression matrix: {expr.shape[0]} genes × {expr.shape[1]} samples")

# Get gene list from expression data
expr_genes = set(expr.index.str.upper())
print(f"\nTotal genes in dataset: {len(expr_genes)}")

# Check Hallmark UPR coverage
print("\n" + "=" * 60)
print("HALLMARK UPR GENE COVERAGE")
print("=" * 60)

hallmark_found = []
hallmark_missing = []

for gene in HALLMARK_UPR:
    gene_upper = gene.upper()
    if gene_upper in expr_genes or gene in expr.index:
        hallmark_found.append(gene)
    else:
        # Check aliases
        found_via_alias = False
        for alias, canonical in ALIASES.items():
            if canonical.upper() == gene_upper and alias.upper() in expr_genes:
                hallmark_found.append(gene)
                found_via_alias = True
                break
        if not found_via_alias:
            hallmark_missing.append(gene)

print(f"Found: {len(hallmark_found)}/{len(HALLMARK_UPR)} ({100*len(hallmark_found)/len(HALLMARK_UPR):.1f}%)")
print(f"Missing: {hallmark_missing}")

# Check Core UPR coverage
print("\n" + "=" * 60)
print("CORE UPR GENE COVERAGE")
print("=" * 60)

core_found = []
core_missing = []

for gene in CORE_UPR:
    gene_upper = gene.upper()
    if gene_upper in expr_genes or gene in expr.index:
        core_found.append(gene)
    else:
        core_missing.append(gene)

print(f"Found: {len(core_found)}/{len(CORE_UPR)} ({100*len(core_found)/len(CORE_UPR):.1f}%)")
print(f"Missing: {core_missing}")

# Check Lab Panel coverage
print("\n" + "=" * 60)
print("LAB PANEL GENE COVERAGE")
print("=" * 60)

lab_found = []
lab_missing = []

for gene in LAB_PANEL:
    gene_upper = gene.upper()
    if gene_upper in expr_genes or gene in expr.index:
        lab_found.append(gene)
    else:
        lab_missing.append(gene)

print(f"Found: {len(lab_found)}/{len(LAB_PANEL)} ({100*len(lab_found)/len(LAB_PANEL):.1f}%)")
if lab_missing:
    print(f"Missing: {lab_missing}")

# Sample check: look for key UPR genes in actual index
print("\n" + "=" * 60)
print("KEY UPR GENES - ACTUAL VALUES CHECK")
print("=" * 60)

key_genes = ["HSPA5", "XBP1", "ATF4", "ATF6", "ERN1", "EIF2AK3", "DDIT3", "CHAC1", "CALR"]
for gene in key_genes:
    if gene in expr.index:
        vals = expr.loc[gene]
        print(f"{gene}: min={vals.min():.2f}, max={vals.max():.2f}, mean={vals.mean():.2f}")
    else:
        print(f"{gene}: NOT FOUND")

# Save UPR gene expression subset
print("\n" + "=" * 60)
print("SAVING UPR GENE SUBSET")
print("=" * 60)

# Get all available UPR genes
all_upr = list(set(hallmark_found + core_found + LAB_PANEL))
upr_in_expr = [g for g in all_upr if g in expr.index]

upr_expr = expr.loc[upr_in_expr]
upr_expr.to_csv(DATA_DIR / "upr_expression_matrix.csv")
print(f"Saved {len(upr_in_expr)} UPR genes to upr_expression_matrix.csv")

# Summary
print("\n" + "=" * 60)
print("SUMMARY")
print("=" * 60)
print(f"Total genes in dataset: {len(expr_genes)}")
print(f"Hallmark UPR coverage: {len(hallmark_found)}/{len(HALLMARK_UPR)} ({100*len(hallmark_found)/len(HALLMARK_UPR):.1f}%)")
print(f"Core UPR coverage: {len(core_found)}/{len(CORE_UPR)} ({100*len(core_found)/len(CORE_UPR):.1f}%)")
print(f"Lab panel coverage: {len(lab_found)}/{len(LAB_PANEL)} ({100*len(lab_found)/len(LAB_PANEL):.1f}%)")
