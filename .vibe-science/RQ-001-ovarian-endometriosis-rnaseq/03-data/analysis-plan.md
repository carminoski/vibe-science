# UPR Analysis Plan - GSE157153

## Dataset Overview

| Property | Value |
|----------|-------|
| GEO Accession | GSE157153 |
| Title | Gene expression profiling in endometriosis and ovarian cancer |
| Platform | GPL17303 (Ion Torrent Proton RNA-seq) |
| Samples | 66 |
| Genes covered | ~20,803 |
| SRA | SRP279373 |

## Sample Groups (6 conditions)

| Group | Description | Expected N |
|-------|-------------|------------|
| 1 | Endometriosis (benign) | ~11 |
| 2 | Atypical endometriosis | ~11 |
| 3 | Endometriosis adjacent to clear cell OC | ~11 |
| 4 | Endometriosis adjacent to endometrioid OC | ~11 |
| 5 | Clear cell ovarian carcinoma | ~11 |
| 6 | Endometrioid ovarian carcinoma | ~11 |

## Analysis Objectives

### Primary Objective
Characterize UPR pathway activation across the endometriosis-to-cancer progression spectrum.

### Hypotheses

**H1:** UPR genes are differentially expressed between benign endometriosis and EAOC.

**H2:** UPR activation increases progressively along the transformation spectrum (benign → atypical → adjacent → cancer).

**H3:** Specific UPR branches (IRE1/XBP1, PERK/ATF4, ATF6) show distinct patterns in different histotypes (clear cell vs endometrioid).

## Analysis Pipeline

### Phase 1: Data Acquisition

```bash
# Download from GEO
# Option A: Processed counts
wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE157nnn/GSE157153/suppl/GSE157153_RAW.tar

# Option B: Raw reads from SRA
prefetch SRP279373
fastq-dump --split-files SRR*
```

### Phase 2: Quality Control

1. **Sample QC:**
   - Check for outliers (PCA, hierarchical clustering)
   - Verify sample annotation matches data
   - Check for batch effects

2. **Expression QC:**
   - Filter low-expressed genes
   - Normalize (TMM or DESeq2 normalization)
   - Verify UPR genes are detected

### Phase 3: Differential Expression

**Comparisons:**

| Comparison | Biological Question |
|------------|---------------------|
| Cancer vs Benign endo | Overall transformation signature |
| Atypical vs Benign | Early transformation |
| Adjacent vs Atypical | Pre-malignant changes |
| Cancer vs Adjacent | Final transformation step |
| Clear cell vs Endometrioid | Histotype differences |

**Method:** DESeq2 with appropriate design matrix

**Thresholds:**
- |log2FC| > 1
- padj < 0.05

### Phase 4: UPR-Focused Analysis

1. **Targeted DE:**
   - Filter results for 113 Hallmark UPR genes
   - Create heatmap of UPR genes across groups

2. **Pathway Enrichment:**
   - GSEA with Hallmark UPR gene set
   - Test enrichment in each comparison

3. **Branch-Specific Analysis:**
   - IRE1 branch: ERN1, XBP1 (spliced form markers), target genes
   - PERK branch: EIF2AK3, ATF4, DDIT3, CHAC1
   - ATF6 branch: ATF6, HSPA5, CALR, PDIA genes

4. **Progression Analysis:**
   - Order groups by transformation stage
   - Test for linear trend in UPR activation
   - Identify genes with monotonic increase/decrease

### Phase 5: Validation Approaches

1. **Internal validation:**
   - Bootstrap confidence intervals
   - Leave-one-out cross-validation of findings

2. **External validation:**
   - Compare with GSE230956 (independent OCCC dataset)
   - Check consistency with Ciavattini 2018 results (ATF6, GRP78, CHOP, XBP1)

3. **Biological plausibility:**
   - Correlate UPR with oxidative stress markers (if available)
   - Check consistency with ARID1A mutation status (known in clear cell)

## Expected Outputs

1. **Figure 1:** Heatmap of UPR genes across 6 conditions
2. **Figure 2:** GSEA enrichment plots for UPR pathway
3. **Figure 3:** UPR branch-specific activation patterns
4. **Figure 4:** Progression analysis (UPR score vs transformation stage)
5. **Table 1:** Top differentially expressed UPR genes
6. **Table S1:** Full DE results for all UPR genes

## Statistical Considerations

- **Multiple testing:** BH correction within each comparison
- **Power:** With ~11 samples/group and expected effect sizes from cancer studies, power is adequate for pathway-level analysis
- **Confounders:** Age, batch (if applicable) should be included as covariates

## Timeline (Estimated)

| Phase | Task | Duration |
|-------|------|----------|
| 1 | Data download | 1 day |
| 2 | QC | 1 day |
| 3 | DE analysis | 2 days |
| 4 | UPR analysis | 2 days |
| 5 | Validation | 2 days |
| | **Total** | **~1 week** |

## Success Criteria

- [ ] UPR genes show significant enrichment in at least one comparison
- [ ] At least 10 UPR genes differentially expressed (padj < 0.05)
- [ ] Findings consistent between clear cell and endometrioid where applicable
- [ ] Results replicate Ciavattini 2018 direction for ATF6/GRP78/CHOP/XBP1
- [ ] External validation in GSE230956 shows concordant direction

## Kill Conditions

- UPR genes not differentially expressed in any comparison
- Results contradicted by external dataset
- Batch effects confound biological signal
