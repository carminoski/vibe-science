# UPR Gap Discovery

**Date:** 2026-01-30
**Cycle:** 1
**Type:** MAJOR FINDING

## Summary

Identified a significant research gap at the intersection of Unfolded Protein Response (UPR), endometriosis, and ovarian cancer. Despite substantial literature on UPR in ovarian cancer (187 papers) and moderate work on UPR in endometriosis (43 papers), the specific intersection of all three has only **2 papers**.

## Evidence

### Literature Landscape

| Search Query | Results |
|--------------|---------|
| UPR/ER stress + ovarian cancer | 187 papers |
| UPR/ER stress + endometriosis | 43 papers |
| **UPR + endometriosis + ovarian cancer** | **2 papers** |

### The Two Papers

#### 1. Ciavattini 2018 (Key Paper)
- **Title:** "Unfolded protein response, a link between endometrioid ovarian carcinoma and endometriosis: A pilot study"
- **Journal:** Oncology Letters
- **DOI:** 10.3892/ol.2018.9381
- **Citations:** 13
- **What they did:**
  - Retrospective pilot study
  - Compared UPR gene expression across tissue samples
  - Only 4 genes: ATF6, GRP78, CHOP, XBP1
- **Limitations:**
  - Small pilot study
  - Limited gene panel
  - **NO PUBLIC DATA**
  - No transcriptome-wide analysis

#### 2. Singla 2022
- **Title:** Review on nanomaterials for ER stress modulation in ovarian cancer treatment
- **Relevance:** Minimal - focuses on drug delivery, not mechanistic understanding
- Not useful for our purposes

### Available Datasets (Not Yet Analyzed for UPR)

| GEO ID | Title | Samples | Data Type |
|--------|-------|---------|-----------|
| GSE157153 | Gene expression in endometriosis progression to ovarian cancer | 66 | RNA-seq |
| GSE230956 | Clear-cell carcinoma with concurrent endometriosis | 16 | Bulk + small-RNA |
| GSE226575 | FOS in EAOC progression | 9 | RNA-seq |
| GSE291389 | IL-17C signaling in EAOC | 3 | Spatial transcriptomics |
| Nature Genetics 2023 | Single-cell atlas of endometriosis | 370,000 cells | scRNA-seq |

## Gap Analysis

### What Exists
- General UPR studies in ovarian cancer (mechanistic, drug response)
- General UPR studies in endometriosis (oxidative stress, granulosa cells)
- One pilot study (Ciavattini 2018) suggesting UPR as transformation link

### What's Missing
1. **Transcriptome-wide UPR analysis** in endometriosis → EAOC transformation
2. **Single-cell resolution** of UPR pathway in endometriosis/EAOC
3. **Public dataset** with UPR-focused analysis
4. **Systematic comparison** across disease stages
5. **UPR gene signature** for EAOC risk prediction

## Opportunity

**Proposed Analysis:**
Apply comprehensive UPR pathway analysis to existing public datasets:

1. **Use GSE157153** (66 samples across endometriosis → cancer spectrum)
   - Extract UPR pathway genes (full KEGG/Reactome pathways)
   - Differential expression analysis
   - Pathway enrichment
   - Stage-specific signatures

2. **Use Nature Genetics 2023 scRNA-seq data**
   - Cell type-specific UPR activation
   - Trajectory analysis with UPR overlay
   - Identify which cells show UPR dysregulation

3. **Integrate multiple datasets**
   - Meta-analysis of UPR genes across EAOC studies
   - Build UPR-based risk signature

## Confidence Assessment

**Confidence: HIGH**

Reasons:
- Gap confirmed by Scopus search (2 vs 187 papers)
- Ciavattini 2018 explicitly calls for follow-up studies
- Multiple public datasets available with RNA-seq data
- UPR is well-characterized pathway (gene lists available)
- Clear path to novel contribution

## Counter-Evidence Searched

- Searched for UPR analyses in existing EAOC papers: NOT FOUND
- Checked if Nature Genetics 2023 mentioned UPR: NOT MENTIONED in abstract
- Verified no recent 2024-2025 papers filled this gap: CONFIRMED STILL OPEN

## Next Steps

1. Retrieve full Ciavattini 2018 paper
2. Download GSE157153 dataset
3. Compile comprehensive UPR gene list
4. Verify data quality and sample annotations
5. **INVOKE REVIEWER 2** to challenge this finding
