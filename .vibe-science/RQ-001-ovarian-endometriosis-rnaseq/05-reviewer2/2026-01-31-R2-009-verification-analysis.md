# FINDING-009 Verification Analysis

**Date:** 2026-01-31
**Purpose:** Systematic verification of Reviewer 2 criticisms
**Status:** IN PROGRESS

---

## Executive Summary

After thorough investigation, **Reviewer 2's criticisms are largely valid**. Several claims in FINDING-009 are overinterpreted or based on inconsistent cross-dataset comparisons.

### Verdict by Issue

| # | Issue | R2 Correct? | Action Required |
|---|-------|-------------|-----------------|
| 1 | Platform error (Ion Torrent) | ✅ YES | Already documented in code |
| 2 | XBP1 splicing logical jump | ✅ **YES** | Major revision needed |
| 3 | N too small, IL6 outlier risk | ✅ YES | Document uncertainty |
| 4 | IL6 epithelial attribution | ⚠️ PARTIAL | IL6 also in stroma |
| 5 | Literature gap claim | ✅ YES | Refine wording |

---

## Issue 1: Platform Error - VERIFIED

GSE157153 uses **Ion Torrent Proton (GPL17303)**, not Illumina.

**Evidence:** `03-data/GSE157153/download_gse157153.py` line 4:
```
Platform: GPL17303 (Ion Torrent Proton RNA-seq)
```

**Documentation Status:** Already correct in code, but FINDINGS.md may reference "Illumina" incorrectly.

**Action:** Audit all documentation for platform references.

---

## Issue 2: XBP1 Splicing - CRITICAL FINDING

### The Problem

We claimed: "IRE1 → XBP1 splicing → IL6 transcription"

**But the data shows OPPOSITE patterns in bulk vs single-cell:**

| Dataset | Comparison | XBP1 log2FC | padj | Status |
|---------|------------|-------------|------|--------|
| Bulk GSE157153 | Mixed CCOC+EnOC vs Endo | **+1.74** | 0.0002 | **Significant** |
| scRNA CCOC | CCOC epithelium vs nOV | **-0.84** | 0.40 | **Not significant** |

### Root Cause

The bulk and single-cell analyses use **different comparisons**:

1. **Bulk:** CCOC (n=17) + EnOC (n=12) combined vs pure endometriosis (n=9)
2. **scRNA:** CCOC epithelium only (n=3) vs normal ovary

These cannot be directly compared.

### Branch-Specific Analysis (CCOC vs nOV, Epithelium)

| Branch | Gene | log2FC | padj | Status |
|--------|------|--------|------|--------|
| **PERK** | EIF2AK3 | **+4.08** | **0.002** | **UP** |
| **PERK** | DDIT3/CHOP | **+2.64** | **0.04** | **UP** |
| PERK | ATF4 | -0.28 | 0.76 | NS |
| PERK | ASNS | -7.57 | 0.046 | DOWN |
| IRE1 | ERN1 | +1.71 | 0.08 | Trending |
| **IRE1** | **XBP1** | **-0.84** | **0.40** | **NS** |
| IRE1 | DNAJB9 | -1.43 | 0.20 | NS |
| ATF6 | ATF6 | +0.37 | 0.70 | NS |
| ATF6 | HSPA5/GRP78 | +0.01 | 0.99 | NS |

### Conclusion

**The XBP1 claim is NOT supported by single-cell CCOC data.**

- PERK branch (EIF2AK3, DDIT3) IS activated
- IRE1/XBP1 branch is NOT activated in CCOC epithelium
- The mechanistic model must be revised

**Revised interpretation:**
- CCOC shows PERK branch activation, not IRE1/XBP1
- XBP1 +1.74 from bulk reflects mixed tumor effect, not CCOC-specific
- Cannot claim "XBP1 splicing drives IL6" without junction read evidence

---

## Issue 3: N Too Small - VERIFIED

### Sample Sizes

| Dataset | CCOC samples | Cells/sample |
|---------|--------------|--------------|
| scRNA | 3 | 330, 1008, 2554 |

### IL6 Statistics

From DGE_CCOC_vs_nOV.csv (Epithelium):
```
baseMean: 299.14
log2FC: +10.98
lfcSE: 2.07  ← HIGH UNCERTAINTY
stat: 5.31
padj: 8.66e-06
```

**lfcSE = 2.07** means 95% CI is approximately:
- Lower: 10.98 - 2×2.07 = 6.84
- Upper: 10.98 + 2×2.07 = 15.12

This is a **wide confidence interval** for a log2FC estimate.

### Outlier Risk - CONFIRMED

**Per-sample IL6 counts (Epithelium):**

| Sample | IL6 counts | % of CCOC total |
|--------|------------|-----------------|
| CCCOC_1 | 0 | 0% |
| CCCOC_2 | 42 | 3.7% |
| **CCCOC_3** | **1080** | **96.3%** |
| nOV (8 samples) | 0-1 | - |

**CONFIRMED: ONE SAMPLE DRIVES 96% OF THE IL6 SIGNAL.**

The log2FC of +10.98 is almost entirely attributable to CCCOC_3. If this sample were excluded:
- Remaining CCOC IL6: 0 + 42 = 42 counts (2 samples)
- log2FC would drop dramatically

This validates Reviewer 2's concern about instability with N=3.

---

## Issue 4: IL6 Epithelial Attribution - PARTIALLY VALID

### Cross-Compartment Analysis (CCOC vs nOV)

| Compartment | baseMean | log2FC | padj | Status |
|-------------|----------|--------|------|--------|
| **Epithelium** | 299 | **+10.98** | 8.7e-06 | **Significant** |
| **Stroma** | 185 | **+3.64** | 0.012 | **Significant** |
| Immune | 10 | +1.31 | 0.50 | NS |

### Interpretation

IL6 is elevated in BOTH epithelium AND stroma in CCOC. This suggests:

1. **Possibility A:** Tumor-intrinsic IL6 expression (real epithelial signal)
2. **Possibility B:** Paracrine signaling (stroma → epithelium or vice versa)
3. **Possibility C:** Technical artifact (ambient RNA, doublets)

### STAT3 Target Coherence (CCOC Epithelium)

If IL6→JAK/STAT3 is active, downstream targets should be up:

| Target | log2FC | padj | Status |
|--------|--------|------|--------|
| SOCS3 | +2.50 | 0.13 | Trending |
| STAT3 | +1.27 | 0.37 | NS |
| BCL2L1 | +1.04 | 0.58 | NS |
| MYC | -0.81 | 0.71 | NS |
| IL6R | +0.20 | 0.96 | NS |
| BIRC5 | -0.37 | 0.91 | NS |

**STAT3 targets show WEAK signal** - not strongly coherent with IL6 +10.98.

### Doublet Analysis

From doublet_report.csv:
```
CCCOC_1: 330 cells, 3 doublets (0.91%)
CCCOC_2: 1008 cells, 4 doublets (0.40%)
CCCOC_3: 2554 cells, 7 doublets (0.27%)
```

Doublet rates are low (0.27-0.91%), but this doesn't specifically address IL6+ cell doublet enrichment.

### Missing Analyses (Required)

1. [ ] % cells IL6+ per cluster (not mean expression)
2. [ ] UMI distribution for IL6 across cell types
3. [ ] Doublet scores specifically for IL6+ epithelial cells
4. [ ] Ambient RNA assessment (no SoupX/CellBender results found)

---

## Issue 5: Literature Gap Claim - VALID

### Original Claim
"UPR + EAOC = 2-3 papers only"

### R2 Response
Reviews on UPR in gynecological cancers exist (Tandfonline 2025).

### Refined Claim
"No published model compares EnOC vs CCOC UPR trajectories in the endometriosis-associated context with divergent downstream pathways."

This is more defensible and specific.

---

## Summary of Required Corrections

### FINDING-009 Must Be Revised

| Section | Current Claim | Required Revision |
|---------|---------------|-------------------|
| XBP1 mechanism | "IRE1→XBP1→IL6" | "PERK branch activation observed; IRE1/XBP1 role unconfirmed" |
| UPR activation | Both branches active | PERK dominant in CCOC; IRE1/XBP1 not significant |
| IL6 source | Epithelium-specific | Present in both epithelium and stroma |
| Confidence | HIGH | **MODERATE** - requires validation |

### Documentation Updates Required

1. **FINDINGS.md:**
   - Line 69: Change "IRE1/XBP1 branch: STRONGLY activated" to acknowledge scRNA contradiction
   - Add caveat about bulk vs scRNA comparison limitations

2. **FINDING-009 document:**
   - Remove "XBP1 splicing drives IL6" claim
   - Add cross-compartment IL6 data
   - Document STAT3 target incoherence
   - Lower confidence level

3. **STATE.md:**
   - Update status to reflect verification issues

### New Analyses Required

1. Per-sample IL6 expression (leave-one-out)
2. % IL6+ cells per cluster
3. External validation dataset search (TCGA-OV histotype split)

---

## Conclusion

Reviewer 2's critique is **scientifically rigorous and largely correct**. The divergent mechanism hypothesis (CCOC→IL6 vs EnOC→OxPhos) may still hold, but:

1. The mechanistic pathway (XBP1→IL6) is NOT supported by CCOC single-cell data
2. IL6 elevation is real but compartment attribution is unclear
3. STAT3 target activation is weak
4. Sample size limitations create high uncertainty

**Recommendation:** Downgrade FINDING-009 confidence from HIGH to MODERATE, revise mechanistic claims, and explicitly document limitations.

---

*Analysis completed: 2026-01-31*
*Reviewer 2 criticisms: 4/5 validated*
