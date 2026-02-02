# Reviewer #2 Critique - FINDING-009

**Received:** 2026-01-31
**Verdict:** MAJOR REVISION
**Severity:** HIGH - Multiple foundational issues

---

## Summary of Criticisms

| # | Issue | Severity | Status |
|---|-------|----------|--------|
| 1 | Platform error (GSE157153 = Ion Torrent, NOT Illumina) | CRITICAL | TO FIX |
| 2 | XBP1 splicing: logical jump (expression ≠ activity) | MAJOR | TO ADDRESS |
| 3 | N too small, FDR unstable, IL6 could be outlier | MAJOR | TO VERIFY |
| 4 | IL6 "epithelial": doublet/ambient RNA risk | MAJOR | TO VERIFY |
| 5 | Literature gap claim too broad | MINOR | TO REFINE |

---

## MAJOR ISSUE 1: Platform Error (CRITICAL)

### The Problem

> GSE157153 uses **GPL17303 = Ion Torrent Proton**, NOT Illumina.

This affects:
- Bias profiles
- Coverage patterns
- Quantification (especially low-count genes like IL6)
- Comparability with 10x scRNA-seq

### Action Required

- [ ] Correct ALL documentation (STATE.md, FINDINGS.md, CONSOLIDATED-SUMMARY.md)
- [ ] Add methodological note on Ion Torrent vs Illumina bias
- [ ] Run sensitivity analysis with stricter low-expression filters
- [ ] Verify IL6 remains significant with shrinkage estimators

### Files to Update

```
FINDINGS.md - Line mentioning "Illumina"
CONSOLIDATED-SUMMARY.md - Dataset description
01-discovery/2026-01-30-dataset-selection-criteria-analysis.md
03-data/GSE157153/*.py scripts
```

---

## MAJOR ISSUE 2: XBP1 Splicing Logical Jump

### The Problem

> "XBP1 active" requires **XBP1s (spliced)**, not just XBP1 expression.
> UPR is branch-specific: PERK/ATF4, ATF6, IRE1/XBP1 can move discordantly.
> Ciavattini 2018 showed ATF6/GRP78 ↑ but CHOP/XBP1 ↓ in endometrioid.

### Current Claim (WRONG)

"IRE1/XBP1 → IL6 → JAK/STAT3"

### Corrected Claim

"Pattern compatible with remodeled UPR branches; branch-specific activity requires validation"

### Action Required

- [ ] Calculate separate scores for PERK/ATF4, ATF6, IRE1/XBP1 branches
- [ ] Use Reactome/GO curated gene sets, not just Hallmark
- [ ] Check if bulk data has junction reads for XBP1s estimation
- [ ] Align interpretation with Ciavattini findings
- [ ] Remove causal language ("XBP1 splicing drives IL6")

### Branch-Specific Gene Sets Needed

| Branch | Key Genes |
|--------|-----------|
| IRE1/XBP1 | ERN1, XBP1, DNAJB9, EDEM1, SEC61A1 |
| PERK/ATF4 | EIF2AK3, ATF4, DDIT3, ASNS, TRIB3 |
| ATF6 | ATF6, HSPA5, CALR, PDIA4, XBP1 (unspliced target) |

---

## MAJOR ISSUE 3: N Too Small, Unstable Statistics

### The Problem

> CCOC n≈3. One outlier patient could invent "+10.98 log2FC" for IL6.
> FDR and TOP50 rankings are unstable with N this low.

### Action Required

- [ ] Report individual sample IL6 values (not just aggregate)
- [ ] Calculate leave-one-out sensitivity (remove each sample, recompute)
- [ ] Show forest plot or individual datapoints
- [ ] Find external validation dataset (bulk or scRNA)
- [ ] Acknowledge limitation explicitly in manuscript

### Validation Options

| Dataset | Type | Samples | Feasibility |
|---------|------|---------|-------------|
| GSE226870 | Bulk EAOC | ? | Check for CCOC/EnOC split |
| GSE189553 | Bulk CCOC | Yes | Already excluded as not scRNA |
| TCGA-OV | Bulk | Limited CCOC | Check histotype annotation |

---

## MAJOR ISSUE 4: IL6 Epithelial Attribution Risk

### The Problem

> IL6 classically comes from stroma/immune cells.
> "Epithelial IL6" could be:
> - Doublets (epithelial + immune)
> - Ambient RNA contamination
> - Overly permissive cluster annotation

### Action Required

- [ ] Show % IL6+ cells per cluster (not just mean expression)
- [ ] Show UMI distribution for IL6 across cell types
- [ ] Report doublet scores (SOLO/Scrublet) for IL6+ epithelial cells
- [ ] Apply ambient RNA correction (SoupX/CellBender) and re-check
- [ ] Verify with stringent epithelial markers (EPCAM+, KRT+, CDH1+)
- [ ] Check if STAT3 targets (SOCS3, IL6R) are coherent in same compartment

### Key Question

Is IL6 **tumor-intrinsic** or **TME-derived** in CCOC?

Literature (AACR) shows IL6-STAT3-HIF axis is active in OCCC, but doesn't confirm tumor-intrinsic source.

---

## MINOR ISSUE 5: Literature Gap Claim

### The Problem

> "UPR+EAOC = 2-3 papers" is too broad.
> UPR in gynecological/ovarian cancer has reviews (Tandfonline 2025).

### Corrected Claim

> "No published model compares **EnOC vs CCOC UPR trajectories** in the **endometriosis-associated context** with divergent downstream pathways."

---

## Reviewer's Minimum Requirements

| # | Requirement | Priority |
|---|-------------|----------|
| 1 | Platform correction + sensitivity analysis | IMMEDIATE |
| 2 | Branch-specific UPR scores (PERK, ATF6, IRE1) | HIGH |
| 3 | IL6 cellular attribution (doublets, ambient, cluster %) | HIGH |
| 4 | At least 1 external validation (bulk or scRNA) | HIGH |
| 5 | Alternative explanation audit (confounders) | MEDIUM |

### Confounders to Address

**Clinical:**
- Stage, grade
- Prior treatment
- Tumor purity

**Technical:**
- Library size normalization
- Percent mitochondrial
- Dissociation stress
- Batch effects

**Biological:**
- Endometrioma hypoxia
- Hemorrhage artifacts
- Macrophage infiltrate

---

## Response Plan

### Phase 1: Corrections (Immediate)

1. Fix platform error everywhere
2. Soften causal language
3. Refine literature gap claim

### Phase 2: Verification (Required)

1. Calculate branch-specific UPR scores
2. Individual sample IL6 analysis
3. Doublet/ambient RNA check for IL6+ epithelial cells

### Phase 3: Validation (If Possible)

1. Find external dataset with CCOC/EnOC split
2. Replicate UPR/IL6 pattern

### Phase 4: Resubmission

1. Update all documentation
2. Add methodological caveats
3. Resubmit to R2

---

## Key Lesson

> **"La direzione può diventare Q1-worthy, ma solo se smetti di chiudere il cerchio a parole e lo chiudi con i dati."**

The finding is interesting but over-sold relative to the evidence strength.

---

*Critique received: 2026-01-31*
*Response plan initiated: 2026-01-31*
