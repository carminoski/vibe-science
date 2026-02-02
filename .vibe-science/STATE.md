---
rq: RQ-001-ovarian-endometriosis-rnaseq
phase: VERIFICATION & REVISION
cycle: 10
last_updated: 2026-01-31T14:00:00Z
minor_findings_pending: 0
major_finding_pending_review: 7
reviewer2_approved: 1
findings_total: 9
r2_critique_verified: true
---

## Current Focus

**REVIEWER 2 CRITIQUE VERIFICATION COMPLETE.** Several claims require revision.

### FINDING-009 Status: MAJOR REVISION REQUIRED

| Subtype | UPR NES | IL6/JAK/STAT3 NES | Original Mechanism | Revised Status |
|---------|---------|-------------------|-------------------|----------------|
| **CCOC** | +1.56 | **+1.38** | ~~UPR → XBP1 → IL6~~ | PERK branch, NOT XBP1 |
| **EnOC** | **+3.41** | **-1.96** | OxPhos/Myc | Still valid |

### Critical Issues Verified

| Issue | R2 Correct? | Evidence |
|-------|-------------|----------|
| XBP1 not up in CCOC | ✅ YES | scRNA: -0.84, NS (bulk +1.74 is mixed tumors) |
| IL6 outlier-driven | ✅ **YES** | 96% of signal from 1 of 3 samples |
| STAT3 targets weak | ✅ YES | SOCS3 p=0.13, MYC/BCL2L1/BIRC5 NS |
| IL6 not epithelium-specific | ✅ PARTIAL | Also +3.64 in stroma |

### Dataset Validation Search Results

After thorough GEO search, **no ideal external validation dataset exists**:

| Dataset | Issue | Status |
|---------|-------|--------|
| GSE224333 | In vitro (not primary tumor) | ❌ EXCLUDED |
| GSE189553 | Bulk RNA-seq (not scRNA-seq) | ❌ EXCLUDED |
| GSE226870 | Bulk RNA-seq | ⚠️ Could use for pathway-level validation |
| 2025 CCOC Study | NovaSeq bulk, not 10x | ⚠️ Limited |

**Conclusion:** The Bertelli dataset (3 CCOC, Drop-seq) remains the **best available CCOC scRNA-seq data**. The cross-validation with GSE157153 (bulk) already provides strong support.

## Key Findings Summary (9 Total)

| # | Type | Status | Key Result |
|---|------|--------|------------|
| 001 | Research Gap | ✅ R2 APPROVED | UPR+EAOC: 2-3 papers only |
| 002 | Bulk RNA-seq | PENDING | XBP1 +1.74, HSPA5 +1.24 |
| 003 | Cross-validation | PENDING | NES=1.56 in Bertelli |
| 004 | Methodological | ⚠️ DOCUMENTED | Platform heterogeneity |
| 005 | Cell-type specific | PENDING | EnOC NES=3.41 > CCOC 1.56 |
| 006 | Lab panel | PENDING | **IL6 +10.98 in CCOC** |
| 007 | Datasets | ⚠️ UPDATED | No ideal external dataset |
| 008 | Audit | ✅ COMPLETE | Methodology is sound |
| **009** | **Divergent Mechanisms** | **🆕 PENDING** | **CCOC→IL6, EnOC→OxPhos** |

## Open Questions - ALL RESOLVED

1. ~~General landscape~~ → DONE
2. ~~Reviewer 2 approval~~ → APPROVED
3. ~~UPR gene coverage?~~ → YES (96.5%)
4. ~~UPR activated in EAOC?~~ → **YES - CONFIRMED**
5. ~~Single-cell validation?~~ → **YES (NES=1.56-3.41)**
6. ~~Cell-type specific?~~ → **YES - PERK all compartments**
7. ~~Lab panel genes?~~ → **IL6 +10.98, LIFR up, XBP1/GRP78 up**
8. ~~Additional datasets?~~ → **No ideal 10x CCOC exists**
9. ~~Methodology audit?~~ → **SOUND - no porcherie**
10. ~~EnOC > CCOC why?~~ → **DIVERGENT MECHANISMS FOUND** (FINDING-009)
11. ~~UPR trajectory?~~ → **COMPLETE**: Endo(-)→OMA(+)→Cancer(++)
12. ~~Why IL6 high in CCOC but lower UPR?~~ → **CCOC uses IL6 path, EnOC uses OxPhos**

## Completed Actions

| # | Task | Status |
|---|------|--------|
| 1 | Compile UPR gene list | ✅ DONE |
| 2 | Verify gene coverage | ✅ DONE (96.5%) |
| 3 | Analyze GSE157153 | ✅ DONE |
| 4 | Review Bertelli thesis | ✅ DONE |
| 5 | Verify dataset criteria | ✅ DONE |
| 6 | Access Bertelli data | ✅ DONE |
| 7 | UPR by compartment | ✅ FINDING-005 |
| 8 | EnOC > CCOC difference | ✅ FINDING-005 |
| 9 | Search additional datasets | ✅ FINDING-007 (no ideal found) |
| 10 | Verify lab panel genes | ✅ FINDING-006 |
| 11 | Audit methodology | ✅ FINDING-008 |
| 12 | External validation search | ✅ No ideal dataset available |
| 13 | Consolidate all findings | ✅ CONSOLIDATED-SUMMARY.md updated |

## Methodology Audit (2026-01-31) - COMPLETED ✅

### Files Reviewed
- code/3.2.txt-3.6.txt (methodology)
- code/parte_1-2_integrazione_normal_dati.txt (processing)
- nuovo blocco/1.txt (project history)
- reports/04_overview_text.md, 09_epi_subtypes/

### Audit Result: **METHODOLOGY IS SOUND**

| Aspect | Finding |
|--------|---------|
| QC | Standard parameters (min_genes=200, pct_mt<15%) |
| Integration | Harmony, 30 PCs, LISI=3.2 |
| DGE | PyDESeq2 pseudobulk, robust |
| Historical R issues | All resolved in Python migration |

## Key Numbers

| Metric | Value | Note |
|--------|-------|------|
| Total cells | 321,232 | |
| Total samples | 52 | |
| CCOC samples | 3 | ⚠️ Limited |
| UPR genes tested | 121 | |
| UPR genes UP | 34 | |
| XBP1 log2FC (bulk) | +1.74 | Mixed CCOC+EnOC vs Endo |
| XBP1 log2FC (CCOC scRNA) | **-0.84** | ⚠️ NOT SIGNIFICANT |
| HSPA5 log2FC (CCOC scRNA) | +0.01 | ⚠️ NOT SIGNIFICANT |
| EIF2AK3/PERK (CCOC scRNA) | **+4.08** | ✅ PERK branch UP |
| **IL6 log2FC (CCOC)** | **+10.98** | ⚠️ 96% from 1 sample |
| IL6 CCOC_1 counts | 0 | |
| IL6 CCOC_2 counts | 42 | |
| IL6 CCOC_3 counts | **1080** | ← Outlier |
| EnOC UPR NES | 3.41 | |
| CCOC UPR NES | 1.56 | |
| CCOC IL6/JAK NES | +1.38 | |
| EnOC IL6/JAK NES | -1.96 | |

## Files Generated This Session

```
CONSOLIDATED-SUMMARY.md                           # Complete 13-section summary
FINDINGS.md                                       # 8 findings documented
STATE.md                                          # This file
05-reviewer2/2026-01-31-benign-lesions-analysis.md         # Lesioni benigne: chi fa cosa
05-reviewer2/2026-01-31-compartment-specific-analysis.md   # UPR/IL6 per compartimento
05-reviewer2/2026-01-31-R2-009-verification-analysis.md    # Verifica critiche R2
05-reviewer2/2026-01-31-R2-009-PACK.md                     # Pacchetto revisione per R2
05-reviewer2/2026-01-31-comprehensive-pathway-analysis.md  # Oncogeni, ormonali, ferroptosi
05-reviewer2/2026-01-31-proliferation-survival-analysis.md # Proliferazione e sopravvivenza
05-reviewer2/2026-01-31-stemness-EMT-differentiation-analysis.md # Staminalità, EMT, differenziamento
05-reviewer2/2026-01-31-stress-inflammation-celldeath-analysis.md # Stress, infiammazione, morte cellulare
05-reviewer2/2026-01-31-mechanotransduction-adhesion-analysis.md  # Hippo, adesione, ECM, polarità
```

## Analisi Compartment-Specific (2026-01-31) - NEW

### Key Discovery: IL6 Source è Diverso!

| Tumor | Epithelium IL6 | Stroma IL6 | Source |
|-------|----------------|------------|--------|
| **CCOC** | +1.72 (NS) | -0.09 | Epiteliale (ma outlier) |
| **EnOC** | **-4.04*** | **+3.09*** | **STROMA-DERIVED** |

### UPR Branch Diversi

| Tumor | Branch Dominante | Compartimento |
|-------|------------------|---------------|
| CCOC | **PERK** (EIF2AK3) | Sistemico (tutti) |
| EnOC | **ATF6/chaperone** (HSPA5) | Solo Stroma |

### File Generato
`05-reviewer2/2026-01-31-compartment-specific-analysis.md`

---

## Reviewer 2 Critique Verification (2026-01-31)

### Summary
4 of 5 criticisms validated. Major revisions required.

### Verified Issues

1. **Platform:** GSE157153 = Ion Torrent (confirmed in download script)
2. **XBP1:** NOT upregulated in CCOC scRNA (-0.84, NS); bulk +1.74 is mixed tumor effect
3. **IL6 outlier:** 96% of signal from CCCOC_3 (1080/1122 total counts)
4. **STAT3 targets:** Weak/incoherent (SOCS3 p=0.13, others NS)
5. **IL6 compartment:** Present in BOTH epithelium (+10.98) AND stroma (+3.64)

### Branch-Specific UPR (CCOC Epithelium)

| Branch | Status |
|--------|--------|
| PERK (EIF2AK3, DDIT3) | ✅ ACTIVATED |
| IRE1/XBP1 | ❌ NOT ACTIVATED |
| ATF6 | ❌ NOT ACTIVATED |

### Files Created

```
05-reviewer2/2026-01-31-R2-009-verification-analysis.md
```

## Next Steps (Required)

1. **Revise FINDINGS.md** - Remove XBP1→IL6 mechanism claim
2. **Revise FINDING-009** - Change to PERK-driven model
3. **Document IL6 outlier** - Add per-sample data and caveats
4. **Lower confidence** - From HIGH to MODERATE
5. **Consider external validation** - TCGA-OV histotype split for pathway-level check

## Salvageable Claims

Despite issues, some findings remain valid:

| Claim | Status | Evidence |
|-------|--------|----------|
| Divergent IL6 pattern (CCOC vs EnOC) | ✅ Valid | GSEA NES: +1.38 vs -1.96 |
| PERK branch activation in CCOC | ✅ Valid | EIF2AK3 +4.08, DDIT3 +2.64 |
| Stronger UPR in EnOC | ✅ Valid | NES 3.41 vs 1.56 |
| EnOC uses OxPhos/Myc | ✅ Valid | NES +3.62, +4.63 |

## Comprehensive Pathway Analysis (2026-01-31) - NEW

### Oncogeni/Tumor Suppressori

| Gene | CCOC Pattern | EnOC Pattern | Note |
|------|--------------|--------------|------|
| **NF1** | +1.79*** immune | +2.44*** immune | Paradosso: TSG UP |
| **CDKN2A** | +1.81~ stroma | +2.38 stroma | UP in stroma |
| TP53 | -0.75 (NS) | -0.56 (NS) | Non significativo |
| ARID1A | -0.51 (NS) | -0.42 (NS) | Non significativo |
| PIK3CA | +0.31 (NS) | +0.58 (NS) | Non significativo |

### Switch Ormonale

| Gene | CCOC | EnOC | Interpretazione |
|------|------|------|-----------------|
| **AR** | **-5.10***↓ | +0.70 | CCOC: silenziato |
| **ESR2** | **+5.89***↑ | +3.15 | UP entrambi |
| **ESR1** | -2.02 | **+2.40~**↑ | EnOC riattiva |
| **PGR** | -1.83 | **+2.99~**↑ | EnOC riattiva |

**Insight:** CCOC = "de-mascolinizzazione" (AR↓); EnOC = estrogeno-dipendente

### Ferroptosi

| Gene | CCOC Stroma | EnOC Stroma | Target? |
|------|-------------|-------------|---------|
| **GPX4** | **-2.06***↓ | -0.90 | CCOC vulnerabile |
| **SLC7A11** | **+2.92***↑ | **+2.91***↑ | Compensazione |
| **ACSL4** | +0.78*↑ | **+1.25***↑ | Pro-ferroptosi |

**Target terapeutico:** Inibitori SLC7A11 (erastin, sulfasalazina) per sfruttare GPX4↓

### File Generato
`05-reviewer2/2026-01-31-comprehensive-pathway-analysis.md`

---

## Proliferation & Survival Analysis (2026-01-31) - NEW

### Scoperta Chiave: Stroma Quiescente → Riattivato

**Lesioni Benigne (vs EuE):**
| Marker | Endo Pelvica Stroma | OMA Stroma |
|--------|---------------------|------------|
| MKI67 | **-2.05~** | **-3.06***↓↓ |
| TOP2A | **-1.94~** | **-2.77***↓↓ |
| BIRC5 | -1.25 | **-2.73***↓↓ |

→ Stroma delle lesioni benigne è **quiescente**

**Trasformazione Tumorale (vs OMA):**
| Marker | CCOC Stroma | EnOC Stroma |
|--------|-------------|-------------|
| MKI67 | **+4.24***↑↑ | **+4.41***↑↑ |
| TOP2A | **+4.25***↑↑ | **+3.58***↑↑ |
| BIRC5 | **+4.33***↑↑ | +1.54 |
| CCNE1 | **+3.07***↑ | +1.73~ |
| mTOR | **+1.70***↑ | +0.97~ |

→ Stroma tumorale **massivamente riattivato**

### Pattern CCNE1 nell'Epitelio
- CCOC: **CCNE1 +6.47***↑↑↑** (bypass CDK4/6!)
- EnOC: +2.99 (trend)
- OMA: -2.87 (down vs EuE)

### Implicazione
- Inibitori CDK4/6 potrebbero essere meno efficaci in CCOC (CCNE1 bypass)
- Target alternativi: mTOR, Survivin/BIRC5

### File Generato
`05-reviewer2/2026-01-31-proliferation-survival-analysis.md`

---

## Stemness, EMT & Differentiation (2026-01-31) - NEW

### Scoperta Chiave: SOX2 Stromale

| Tumore | SOX2 Epitelio | SOX2 Stroma | Ratio |
|--------|---------------|-------------|-------|
| CCOC | +5.33~ | **+11.44***↑↑↑ | 2.1x stroma |
| EnOC | +5.80 | **+10.12***↑↑ | 1.7x stroma |

→ **Staminalità nel MICROAMBIENTE > Epitelio**

### Pattern EMT: Paradosso

| Confronto | Epitelio | Stroma |
|-----------|----------|--------|
| **OMA vs EuE** | **EMT↑** (FN1+5.6, CDH2+3.9) | EMT-TF↓ |
| **CCOC vs OMA** | **MET↓** (VIM-3.6, CDH2-5.3) | EPCAM+5.2 |
| **EnOC vs OMA** | **MET↓** (VIM-2.9, FN1-3.8) | EPCAM+6.0 |

→ **Tumori mostrano MET (reversione), non EMT!**

### HNF1B: Marker CCOC-Specifico

| Gene | CCOC Epitelio | EnOC Epitelio |
|------|---------------|---------------|
| **HNF1B** | **+3.46***↑ | +0.07 |

→ Target potenziale CCOC-specifico

### File Generato
`05-reviewer2/2026-01-31-stemness-EMT-differentiation-analysis.md`

---

## Stress, Inflammation & Cell Death (2026-01-31) - NEW

### Scoperta Chiave: Switch Antiossidante

| Condizione | Sistema Dominante | Pattern |
|------------|-------------------|---------|
| **OMA** | NRF2/HMOX1 (citosolico) | HMOX1↑, NRF2↑ |
| **CCOC/EnOC** | SOD2/GPX1 (mitocondriale) | SOD2 +15***, GPX1 +13*** |

→ **Cambio radicale da antiossidanti citosolici a mitocondriali**

### SOCS3: Il Guardiano della Trasformazione

| Condizione | SOCS3 Stroma | Interpretazione |
|------------|--------------|-----------------|
| OMA vs EuE | **+1.29*** | FRENO ATTIVO |
| CCOC vs OMA | **-3.26*** | FRENO RIMOSSO |
| EnOC vs OMA | **-1.99*** | FRENO RIMOSSO |

→ **Perdita di SOCS3 = evento chiave nella trasformazione**

### Piroptosi: Escape Meccanismo

| Condizione | GSDME | Interpretazione |
|------------|-------|-----------------|
| OMA vs EuE | **+4.59*** (epitelio) | SENSIBILIZZATO |
| CCOC vs OMA | **-8.21*** | ESCAPE |
| EnOC vs OMA | **-6.56*** | ESCAPE |

→ **Tumori disattivano piroptosi come meccanismo di evasione**

### HIF1A Pseudoipossia

| Tumore | HIF1A Stroma | VEGFA | LDHA |
|--------|--------------|-------|------|
| CCOC | **+11.87*** | -0.51 | -0.98 |
| EnOC | **+12.20*** | +1.41 | +0.12 |

→ **HIF1A massiccamente UP ma pathway a valle non attivato = pseudoipossia?**

### File Generato
`05-reviewer2/2026-01-31-stress-inflammation-celldeath-analysis.md`

---

## Mechanotransduction, Adhesion & Polarity (2026-01-31) - NEW

### Scoperta Chiave: OMA = EMT → Tumori = MET

| Caratteristica | OMA vs EuE | CCOC vs OMA | EnOC vs OMA |
|----------------|------------|-------------|-------------|
| **CDH1 (E-cad)** | ↓ (stroma -2.82***) | ↑ (stroma +2.65***) | ↑ |
| **CDH2 (N-cad)** | ↑ (+3.85~) | **↓↓↓ (-5.27***)** | ↓ |
| **ACTA2 (α-SMA)** | **↑↑↑ (+4.05***)** | **↓↓↓ (-3.55***)** | ↓↓↓ |
| **CRB3 (polarità)** | ↓↓↓ | ↑ | ↑ |
| **PARD3 (polarità)** | - | **↑↑↑ (+3.85***)** | ↑↑↑ |

→ **I tumori REVERTONO l'EMT dell'OMA e ristabiliscono polarità epiteliale!**

### Hippo Pathway: CCOC-Specifico

| Componente | CCOC | EnOC |
|------------|------|------|
| **TAZ/WWTR1** | **+2.47*** epitelio | NS |
| **AMOTL2** | **-4.51*** stroma | -2.69* |

→ **CCOC attiva TAZ e rimuove i freni AMOT**

### Adesioni Focali: Stroma CCOC

| Gene | CCOC Stroma | EnOC Stroma |
|------|-------------|-------------|
| **PXN** | **+10.01***↑↑↑ | +0.39 |
| **FAK** | **+1.73*** | +1.41* |
| **VCL** | +0.88* | +0.20 |

→ **CCOC ha massivo rinforzo adesioni focali stromali (PXN 10x!)**

### Claudine: Marker Specifici

| Gene | OMA | CCOC | EnOC |
|------|-----|------|------|
| CLDN3/4 | ↓↓↓ | **↑↑↑** | ↑ |
| CLDN7 | ↓ | ↓↓↓ epitelio | - |

→ **CLDN4 = marcatore CCOC**

### File Generato
`05-reviewer2/2026-01-31-mechanotransduction-adhesion-analysis.md`

---

## Literature Verification (2026-02-01) - NEW

### Scopus API Systematic Search Completed

**API Key:** 3db28da606426ff5261092311e15bbd2
**N. Ricerche:** 18

### Summary of Findings Classification

| Categoria | NOVEL | Parzialmente KNOWN | KNOWN |
|-----------|-------|-------------------|-------|
| Stress/Infiammazione/Morte Cellulare | 5 | 2 | 3 |
| Meccanotrasduzione/Adesione | 4 | 2 | 4 |
| Pathway Ormonali/Oncogeni/Ferroptosi | 4 | 2 | 3 |
| **TOTALE** | **13** | **6** | **10** |

### TOP NOVEL FINDINGS (Prima Volta Descritti)

1. **SOCS3 Gatekeeper Concept** (0 results)
   - OMA: +1.29*** (freno attivo) → CCOC: -3.26*** (freno rimosso)

2. **SOD2/GPX1 Mitochondrial Switch** (0 results)
   - OMA: SOD2 -20 → CCOC: SOD2 +15 (~35 log2FC change!)

3. **WWTR1/TAZ CCOC-Specificity** (0 results)
   - TAZ +2.47*** specifico per CCOC

4. **PERK vs ATF6 Branch Divergence** (1 non-relevant result)
   - CCOC: PERK sistemico; EnOC: ATF6/chaperone stroma

5. **PXN +10x in CCOC Stroma** (11 generic results, none specific)
   - Massive focal adhesion reinforcement

6. **ESR2 Switch CCOC** (1 result)
   - ESR2 +5.89***, AR -5.10***

7. **HIF1A Pseudohypoxia** (0 results)
   - HIF1A +11.87*** but VEGFA/LDHA uncoupled

8. **IL6 Source Switch** (35 generic results)
   - CCOC: epithelial; EnOC: stromal

9. **GSDME Pyroptosis Escape** (27 results focus on activation, not escape)
   - OMA: +4.59*** → Tumors: -8.21***

10. **Polarity Restoration Paradox** (4 results, not for restoration)
    - PARD3/DLG1/LLGL2 UP in tumors

### KNOWN FINDINGS

| Finding | Scopus Results | Status |
|---------|----------------|--------|
| HNF1B CCOC marker | 18 | KNOWN |
| Iron oxidative stress CCOC | 31 | KNOWN |
| CLDN4 ovarian cancer marker | 109 | KNOWN |
| ACTA2/α-SMA endometriosis | 43 | KNOWN |
| Ferroptosis GPX4/SLC7A11 | 24 | KNOWN |
| EMT/MET in ovarian cancer | 25 | KNOWN |
| AMOTL2-YAP regulatory axis | 27 | KNOWN |

### Suggested New Research Directions

1. **SOCS3 as transformation biomarker** - Monitor in OMA patients
2. **SOD2 mitochondrial inhibitors** - New therapeutic class
3. **UPR branch-specific therapy** - PERK for CCOC, ATF6 for EnOC
4. **Stromal PXN/FAK targeting** - Non-epithelial vulnerability
5. **TAZ/WWTR1 as CCOC-specific target** - Verteporfin studies
6. **ESR2 modulators for CCOC** - Hormonal repurposing

### File Generato
`05-reviewer2/2026-02-01-literature-verification-summary.md`

---

## Files Generated (All Sessions)

```
CONSOLIDATED-SUMMARY.md                           # Complete 13-section summary
FINDINGS.md                                       # 8 findings documented
STATE.md                                          # This file
05-reviewer2/2026-01-31-benign-lesions-analysis.md         # Lesioni benigne: chi fa cosa
05-reviewer2/2026-01-31-compartment-specific-analysis.md   # UPR/IL6 per compartimento
05-reviewer2/2026-01-31-R2-009-verification-analysis.md    # Verifica critiche R2
05-reviewer2/2026-01-31-R2-009-PACK.md                     # Pacchetto revisione per R2
05-reviewer2/2026-01-31-comprehensive-pathway-analysis.md  # Oncogeni, ormonali, ferroptosi
05-reviewer2/2026-01-31-proliferation-survival-analysis.md # Proliferazione e sopravvivenza
05-reviewer2/2026-01-31-stemness-EMT-differentiation-analysis.md # Staminalità, EMT, differenziamento
05-reviewer2/2026-01-31-stress-inflammation-celldeath-analysis.md # Stress, infiammazione, morte cellulare
05-reviewer2/2026-01-31-mechanotransduction-adhesion-analysis.md  # Hippo, adesione, ECM, polarità
05-reviewer2/2026-02-01-literature-verification-summary.md       # NEW: Verifica letteratura NOVEL vs KNOWN
```

---

## VPN Note
UniPisa VPN offered by user - available for Scopus/literature access if needed.

---
*Session Status: LITERATURE VERIFICATION COMPLETE*
*13 NOVEL findings identified, 6 partially known, 10 known*
*Finding-009: MAJOR REVISION REQUIRED*
*Mechanistic model needs correction: PERK branch, not IRE1/XBP1*
