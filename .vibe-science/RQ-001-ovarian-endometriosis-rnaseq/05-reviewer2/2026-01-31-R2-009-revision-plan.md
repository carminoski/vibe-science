# FINDING-009 (REVISED): Divergent UPR Branch Engagement in CCOC vs EnOC
**Revision date:** 2026-01-31  
**Status:** MAJOR REVISION REQUIRED (post-verification)  
**Confidence:** MEDIUM (signal present, mechanism not yet proven)

---

## 0) What changed (non-negotiable corrections)

### A. Platform correction
- GEO **GPL17303** corresponds to **Ion Torrent Proton** (not Illumina).  
  *Implication:* quantify/bias and low-count behavior (e.g., IL6) must be treated carefully.

### B. Mechanistic chain removed
- **Removed:** “IRE1/XBP1 splicing → IL6 → JAK/STAT3” as a primary model.  
- **Reason:** in CCOC epithelium (scRNA), **XBP1 is not up** and not significant; bulk XBP1 signal is mixed-histotype.

### C. Outlier disclosure added
- IL6 +10.98 log2FC in CCOC epithelium is **outlier-driven** (≈96% from a single sample; N=3).  
  *This must be stated explicitly in the Executive Summary and Limitations.*

---

## 1) Revised claim (what is still defensible)

### Claim 1 — Branch-specific UPR in CCOC
CCOC shows a **PERK/eIF2α/ATF4–CHOP**-leaning UPR signature (e.g., EIF2AK3/PERK and DDIT3/CHOP up), whereas **IRE1/XBP1 evidence is weak in CCOC epithelium**.

### Claim 2 — Divergent inflammatory gene-set behavior
Inflammation/IL6-related **gene set** behavior diverges between histotypes (CCOC enriched; EnOC depleted), but **gene-level IL6 is not robust** given outlier sensitivity and small N.

### Claim 3 — Stronger global UPR in EnOC
Across gene sets, EnOC shows **stronger global UPR enrichment** than CCOC, and EnOC also exhibits **metabolic reprogramming signatures** (OxPhos/MYC).

---

## 2) Revised mechanistic hypothesis (clearly labeled as hypothesis)

### Model A — CCOC
**Stress context (oxidative/hypoxic/proteotoxic) → PERK activation → translational repression + ATF4/CHOP program**
- Downstream candidates (hypothesis, not proven here):
  - **NF-κB-linked inflammatory output** (PERK can couple to inflammatory programs)
  - **HIF axis coupling** (known in OCCC contexts)
  - Autophagy / amino-acid stress adaptation

### Model B — EnOC
**More global UPR engagement + metabolic reprogramming**
- Hypothesis: UPR in EnOC is more “broad/organizing” (proteostasis + metabolism), rather than a PERK-skewed stress response.

---

## 3) Why STAT3 evidence is currently insufficient
Even if IL6 gene sets suggest inflammation, a transcript-only readout does **not** guarantee STAT3 pathway activity:
- STAT3 activity is often **post-translational (pSTAT3)**.
- However, the absence/weakness of canonical STAT3 targets (SOCS3, MYC, BCL2L1) in the same compartment is a serious inconsistency.

**Current position:** “IL6→JAK/STAT3 active” is **not supported** by the present transcriptomic evidence.

---

## 4) What must be shown next (decisive validations)

### A. Bulk (Ion Torrent) robustness
1. Re-run DE with **robust methods** (Cook’s distance/outlier handling; leave-one-out) and report whether IL6 remains significant.
2. Report **effect sizes with CIs**, not just log2FC.
3. Confirm that results are not driven by library size/coverage artifacts.

### B. scRNA epithelium credibility
4. Verify **IL6 localization**:
   - % IL6+ cells per epithelial cluster
   - doublet detection
   - ambient RNA correction sensitivity
5. Use **patient-level pseudobulk** to avoid cell-count inflation.

### C. Branch-specific UPR scoring
6. Score PERK vs IRE1 vs ATF6 branches separately (curated gene sets; Reactome/GO).
7. If feasible, estimate **XBP1 splicing proxy** (junction evidence) in bulk; otherwise, do not claim “XBP1 splicing.”

### D. External replication
8. Replicate CCOC PERK-skew and EnOC global UPR in ≥1 independent dataset.

### E. Orthogonal biology
9. If you can access IHC/WB: p-eIF2α, ATF4, CHOP, GRP78; optionally pSTAT3.
10. If only transcriptomics: use multiple STAT3 target sets and show consistency or explicitly state “not supported.”

---

## 5) Suggested wording changes (for the manuscript / finding sheet)

- Replace “XBP1 splicing drives IL6” with:  
  **“CCOC epithelium shows a PERK/ATF4–CHOP-leaning ER-stress program; IL6-related inflammatory gene sets are enriched, but gene-level IL6 is sensitive to an outlier (N=3) and does not currently support a causal IL6→STAT3 model.”**

- Replace “STAT3 pathway active” with:  
  **“STAT3 activation is not supported at transcript level by canonical targets in the present data; phosphorylation-level validation would be required.”**

---

## 6) Minimal reference anchors (IDs only; add full refs in your bib)
- GEO platform GPL17303: Ion Torrent Proton.
- OCCC IL6–STAT3–HIF axis: Anglesio et al., Clin Cancer Res (PubMed: 21343371).
- ER-stress ↔ inflammation/NF-κB crosstalk: Schmitz et al. 2018 (PMC6027367); Tam et al. 2012 (PLOS ONE e0045078).
- PERK inflammatory induction examples: Sheshadri et al. 2021 (Cell Death Discov).
