# AUDIT DEFINITIVO: Pipeline scRNA-seq Originale

**Data audit:** 2026-02-01
**Scopo:** Verificare cosa è stato realmente eseguito nella pipeline originale e identificare lacune metodologiche

---

## SOMMARIO ESECUTIVO

| Step | Eseguito? | Evidenza |
|------|-----------|----------|
| **Ingestione dati (52 campioni)** | ✅ SI | `FASE_1.py`, 52 `.h5ad` in `work/01_qc_h5ad/` |
| **Standardizzazione Ensembl ID** | ✅ SI | GENCODE v43 GTF mapping |
| **QC base (MT%, min_genes, min_cells)** | ✅ SI | `MT_PCT=20%, MIN_GENES=200, MIN_CELLS=3` |
| **Doublet removal (SOLO + Scrublet)** | ✅ SI | `doublet_consensus_report.csv` (52 samples) |
| **Harmony batch integration** | ✅ SI | Log: "Converged after 3 iterations" |
| **UMAP + Leiden clustering** | ✅ SI | `03_harmony_umap_leiden.h5ad` |
| **Cell type annotation** | ✅ SI | Marker-based (EPCAM, PTPRC, COL1A1, PECAM1) |
| **Pseudobulk DGE (PyDESeq2)** | ✅ SI | Formula: `~ state + study_id` |
| **Ambient RNA correction (SoupX/CellBender)** | ❌ **NO** | Nessun codice/output trovato |
| **LOOCV sensitivity analysis** | ❌ **NO** | Mai implementato |
| **Module scores (gene signatures)** | ❌ **NO** | Mai implementato |
| **Compositional analysis (scCODA)** | ❌ **NO** | Mai implementato |

---

## 1. FASI ESEGUITE CORRETTAMENTE

### 1.1 Ingestione Dati (FASE 1)
**File:** `FASE_1.py`

**Operazioni:**
- Caricamento multi-formato (10x MTX, H5, CSV/TSV)
- Conversione nomi geni → Ensembl ID (GENCODE v43)
- Collasso geni duplicati per somma
- QC preliminare per campione
- Output: 52 file `.h5ad` compressi

**Parametri QC:**
```python
MT_PCT    = 20.0   # Max % mitochondrial
MIN_GENES = 200    # Min geni per cellula
MIN_CELLS = 3      # Min cellule per gene
```

### 1.2 Merge e Pulizia (FASE 2)
**File:** `STORIA_0_01.md` (Cella 13)

**Operazioni:**
- Concatenazione 52 campioni (`ad.concat`, outer join)
- Collasso geni duplicati post-merge
- Mapping gruppi biologici (OMA, EuE, nOV, CCOC, EnOC, LGsOC)
- Filtro geni rari (< 20 cellule globali)
- Output: `02_merged_raw.h5ad`

### 1.3 Doublet Removal (FASE 2.5)
**Strategia:** Consenso SOLO + Scrublet

**Evidenza:** `work/doublet_consensus_report.csv`
```
sample_id,n,solo_rate,scru_rate,cons_rate,solo_threshold_used
oma_5_GSM6574533_sample29,17576,0.485,0.138,0.127,0.5
...
(52 campioni processati)
```

**Criterio di esclusione:**
- AND: Classificato come doublet da entrambi i metodi
- OR: Score > 0.90 in almeno un metodo

**Cellule rimosse:** ~10-15% per campione (variabile)

### 1.4 Integrazione Harmony (FASE 3)
**Log esecuzione:** (2025-09-08)
```
harmonypy - INFO - Computing initial centroids with sklearn.KMeans...
harmonypy - INFO - Iteration 1 of 20
harmonypy - INFO - Iteration 2 of 20
harmonypy - INFO - Iteration 3 of 20
harmonypy - INFO - Converged after 3 iterations
```

**Parametri:**
```python
sce.pp.harmony_integrate(
    adata,
    key="sample_id",
    basis="X_pca",
    adjusted_basis="X_pca_harmony",
    max_iter_harmony=20
)
sc.pp.neighbors(adata, use_rep="X_pca_harmony", n_neighbors=15)
```

**Output:** `03_harmony_umap_leiden.h5ad` (321,232 cellule × 33,596 geni)

### 1.5 Analisi Differenziale (PyDESeq2)
**Formula utilizzata:**
```python
formula = "~ state" + (f" + {cov}" if cov else "")
# cov = "study_id" se disponibile
```

**Interpretazione:**
- `state` = condizione biologica (OMA vs EuE, etc.)
- `study_id` = covariata per batch (dataset di origine)
- **NON include** sample_id (batch a livello campione)

**Soglie:**
- FDR < 0.05
- |log2FC| ≥ 1.0

---

## 2. STEP NON ESEGUITI (LACUNE CRITICHE)

### 2.1 Ambient RNA Correction ❌

**Ricerca effettuata:**
```bash
grep -ri "SoupX\|CellBender\|ambient\|decontX" ./
# Risultato: ZERO match (solo "ambiente" italiano)
```

**Impatto:**
- Geni a bassa espressione (IL6, CXCL8) potrebbero essere contaminazione
- Marker "source" vs "destination" non affidabili
- Particolarmente critico per cellule epiteliali che esprimono geni stromali/immuni

### 2.2 Leave-One-Donor-Out (LOOCV) ❌

**Mai implementato.** Nessun codice o output trovato.

**Impatto:**
- Robustezza dei risultati non verificata
- Donors outlier potrebbero guidare i risultati
- Particolarmente critico per CCOC (n=3 donatori)

### 2.3 Module Scores ❌

**Mai calcolati.** Analisi basata su singoli geni.

**Impatto:**
- Attivazione pathway inferita da singoli geni (es. EIF2AK3 mRNA ≠ PERK attivo)
- Nessuna aggregazione dei segnali biologici
- Overclaim su pathway attivati/inibiti

### 2.4 Compositional Analysis (scCODA) ❌

**Mai eseguita.** Nessuna normalizzazione compositional.

**Impatto:**
- Shift nelle proporzioni cellulari confondono i risultati DGE
- Non distinguibile: "gene up in tipo X" vs "più cellule di tipo X"

---

## 3. VALUTAZIONE IMPATTO

### 3.1 Criticità per le Conclusioni Pubblicate

| Finding | Affidabilità | Motivo |
|---------|-------------|--------|
| IL6 come marker OMA | ⚠️ BASSA | Low-count gene senza ambient RNA correction |
| SOCS3 "gatekeeper" | ⚠️ MODERATA | Concetto valido, quantificazione incerta |
| UPR activation (PERK) | ⚠️ BASSA | mRNA ≠ protein activity (kinase) |
| CCOC-specific signatures | ❌ MOLTO BASSA | n=3 donors, no LOOCV |
| Compartment-specific effects | ⚠️ MODERATA | No compositional adjustment |

### 3.2 Cosa è Salvabile

1. **Architettura generale** del microambiente (Epithelium/Stroma/Immune)
2. **Trend qualitativi** (infiammazione in OMA, stress in CCOC)
3. **Doublet removal** è stato eseguito correttamente
4. **Batch correction** (Harmony + study_id in formula) mitiga parzialmente

### 3.3 Cosa Richiede Reprocessing

1. **Ambient RNA correction** → SoupX o CellBender sui raw counts
2. **LOOCV** → Verificare robustezza, escludere outlier donors
3. **Module scores** → Non fare claim su pathway da singoli geni
4. **N donors per CCOC** → Riconoscere limitazione (n=3)

---

## 4. RACCOMANDAZIONI

### 4.1 Per la Tesi (Minima Correzione)

1. **Aggiungere disclaimer** nelle sezioni metodologiche:
   - "Ambient RNA correction was not applied; low-count gene results should be interpreted with caution"
   - "LOOCV sensitivity analysis was not performed"

2. **Downgrade claims:**
   - Da "IL6 is THE source marker" → "IL6 shows differential expression"
   - Da "PERK pathway is activated" → "PERK-related genes show altered expression"

3. **Aggiungere tabella N donors** per ogni confronto

### 4.2 Per Pubblicazione (Reprocessing Richiesto)

1. Run SoupX/CellBender sui raw data
2. Implement LOOCV for all major findings
3. Calculate module scores per pathway
4. Add compositional analysis (scCODA or similar)
5. Increase CCOC sample size if possible

---

## 5. FILE SORGENTE ESAMINATI

| File | Contenuto |
|------|-----------|
| `FASE_1.py` | Pipeline ingestione |
| `STORIA_0_01.md` - `STORIA_0_04.md` | ChatGPT conversation history |
| `STORIA_1.md`, `STORIA_2.md`, `STORIA_3.md` | Methodology documentation |
| `code/Interactive_Processing.ipynb` | Main analysis notebook |
| `work/doublet_consensus_report.csv` | Doublet removal report |
| `code/3.2.txt`, `code/3.5.txt` | Thesis methodology text |

---

## 6. CONCLUSIONE

**La pipeline originale è tecnicamente corretta per gli step eseguiti**, ma presenta **lacune metodologiche significative** che compromettono l'affidabilità di alcune conclusioni:

1. **Ambient RNA correction assente** → claim su geni low-count non affidabili
2. **LOOCV non eseguito** → robustezza non verificata
3. **Module scores non calcolati** → overclaim su pathway activation
4. **CCOC underpowered** → n=3 donors insufficiente

**Azione raccomandata:** Aggiornare il documento `2026-01-31-compartment-specific-analysis.md` con disclaimer appropriati e downgrade delle conclusioni non supportate.
