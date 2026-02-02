---
id: RQ-001
created: 2026-01-30
status: active
priority: high
serendipity_origin: null
---

# Research Question

## Question

What are the unexploited research opportunities at the intersection of ovarian cancer, endometriosis, and transcriptomic analysis (single-cell and bulk RNA-seq)?

## Hypothesis

There exist publicly available RNA-seq datasets (GEO, ArrayExpress, paper repositories, GitHub) on ovarian cancer and/or endometriosis that have not been fully exploited with modern computational methods, creating opportunities for novel insights.

## Sub-Questions

1. **Data Landscape:** What scRNA-seq and bulk RNA-seq datasets exist for:
   - Ovarian cancer (subtypes: high-grade serous, clear cell, endometrioid, mucinous)
   - Endometriosis (eutopic vs ectopic endometrium)
   - Endometriosis-associated ovarian cancer (EAOC)

2. **Methodological Gaps:** What computational approaches have NOT been applied?
   - Trajectory inference
   - Cell-cell communication analysis
   - Deconvolution of bulk data
   - Integration of multiple datasets
   - Machine learning for biomarker discovery

3. **Biological Gaps:** What biological questions remain unanswered?
   - Transformation from endometriosis to cancer
   - Microenvironment differences
   - Cell type-specific expression patterns
   - Drug response signatures

## Success Criteria

- [ ] Catalog of at least 10 relevant GEO datasets with metadata
- [ ] Identification of 2-3 specific research gaps with available data
- [ ] At least one gap with no published papers (confirmed via Scopus)
- [ ] Data accessibility verified (can download, has sufficient samples)

## Data Requirements

**Must have:**
- Public RNA-seq datasets (GEO, ArrayExpress, TCGA)
- Sample annotations (disease state, tissue type)
- Raw or normalized counts (not just summary statistics)

**Nice to have:**
- Single-cell data
- Matched patient samples (endometriosis + ovarian cancer)
- Clinical metadata (survival, stage, treatment)

## Kill Conditions

- No public datasets found for this intersection
- All obvious analyses already published
- Data exists but is inaccessible (restricted access, missing raw data)
- Sample sizes too small for meaningful analysis (n < 10 per group)

## Search Strategy

### Phase 1: Literature Landscape
- Scopus: "endometriosis" AND "ovarian cancer" AND "RNA-seq"
- Scopus: "endometriosis-associated ovarian cancer" AND transcriptom*
- PubMed: endometriosis[MeSH] AND ovarian neoplasms[MeSH] AND RNA-seq

### Phase 2: Dataset Discovery
- GEO: "ovarian cancer" + "RNA-seq" + organism:human
- GEO: "endometriosis" + "RNA-seq" + organism:human
- TCGA: OV (ovarian serous cystadenocarcinoma)
- GitHub: search for analysis repositories

### Phase 3: Gap Identification
- Compare methods used in papers vs available methods
- Check if single-cell data exists and what's been done with it
- Look for datasets mentioned but not fully analyzed
