Ecco il file Markdown completo e strutturato per Canvas (o per la tua documentazione di tesi). Riassume l'intero percorso, le motivazioni tecniche, i problemi risolti e il codice definitivo.

***

# Migrazione Workflow scRNA-seq: da R (Seurat) a Python (Scanpy/scVI)

**Progetto:** Analisi Single-Cell RNA-seq su 52 campioni (Endometriosi, Tumori Ovarici, Tessuti Sani).
**Data:** 1 Febbraio 2026
**Stato:** Fase 2 completata (Merge & Cleaning). Fase 3 (Integrazione) pronta.

---

## 1. Contesto e Motivazione

Inizialmente, l'analisi è stata condotta in **R** utilizzando la suite **Seurat/Monocle3**. Tuttavia, il workflow si è interrotto durante la fase di integrazione e riduzione dimensionale (UMAP) a causa di limiti architetturali di R nella gestione della memoria RAM.
Nonostante l'utilizzo di una istanza AWS con **750GB di RAM**, la costruzione del grafo dei vicini (KNN) e l'integrazione di oltre **390.000 cellule** hanno saturato la memoria, causando crash sistematici.

### Soluzione Adottata
Abbiamo migrato l'intero workflow in **Python**, sfruttando librerie progettate per l'efficienza e la scalabilità:
1.  **Scanpy / AnnData:** Utilizzo di matrici sparse (CSR) e gestione efficiente della memoria (senza duplicazione di oggetti densi).
2.  **scvi-tools:** Integrazione dei dati tramite Deep Learning (Variational AutoEncoders) con training in *mini-batch*. Questo permette di scalare a milioni di cellule senza caricare tutto in RAM.
3.  **UMAP "Subset-and-Transform":** Strategia *fail-safe* che calcola l'embedding su un sottoinsieme di cellule e proietta le restanti, evitando il calcolo di una matrice di adiacenza NxN globale.

---

## 2. Fase 1: Pre-processing dei Singoli Campioni

Prima del merge, ogni campione (formati 10x MTX, H5 o CSV) è stato processato individualmente per uniformare i dati.

**Operazioni svolte:**
*   **Caricamento:** Lettura dei dati grezzi.
*   **Standardizzazione Geni:** Conversione di tutti i nomi geni in **Ensembl ID** (univoci) usando un GTF di riferimento (GENCODE v43).
*   **QC preliminare:** Calcolo metriche QC e salvataggio in formato `.h5ad`.
*   **Sparsità:** Conversione forzata delle matrici in formato `CSR` (Compressed Sparse Row) per risparmiare memoria.

*Risultato:* 52 file `.h5ad` salvati nella cartella `work/01_qc_h5ad`.

---

## 3. Fase 2: Merge, Pulizia Strutturale e Annotazione (Codice Definitivo)

Questa è la fase critica in cui i 52 campioni vengono uniti in un unico oggetto. Abbiamo dovuto risolvere diversi bug logici e strutturali rispetto al workflow R.

### Problemi affrontati e risolti (Bug Fixes)

1.  **Gestione Path Windows:**
    *   *Problema:* Errore `FileNotFoundError` nel caricamento del GTF.
    *   *Fix:* Utilizzo di `pathlib` con percorso assoluto corretto per l'ambiente Anaconda su Windows.

2.  **Duplicazione Geni post-merge:**
    *   *Problema:* `ad.concat` gestisce i geni non comuni creando duplicati con suffisso (es. `ENSG000001...-1`). Questo gonfiava la matrice e impediva analisi corrette.
    *   *Fix:* Implementato un blocco di codice che identifica i duplicati basandosi sull'ID radice e **collassa le colonne per somma** (summing counts).

3.  **Mapping dei Gruppi Biologici (Regex):**
    *   *Problema:* Il token "OMA" (Endometrioma) veniva trovato erroneamente dentro "Carcin**OMA**", classificando i tumori come endometriosi.
    *   *Fix:* Funzione `map_group` riscritta con priorità gerarchica (prima i carcinomi, poi endometriosi) e uso di Regex (`\bOMA\b`) per match esatti.

4.  **Tipi di Dato ed Efficienza:**
    *   *Problema:* Rischio di conversione in `float64` denso.
    *   *Fix:* Forzatura di `adata.layers["counts"]` a `int32` e mantenimento della matrice `X` come sparsa CSR.

### Codice Eseguito (Cella 13)

```python
# Cella 13: Fase 2 — Merge, Pulizia Strutturale e Annotazione (JupyterLab/Anaconda, Windows)

import os, glob, re, gc, gzip
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
from scipy import sparse

print("--- FASE 2: INIZIO MERGE E PULIZIA STRUTTURALE ---")

# 1) Gene map da GTF (Ensembl -> HGNC) con path assoluto richiesto
def load_gtf_gene_map(gtf_path: str) -> dict:
    mp = {}
    with gzip.open(gtf_path, "rt", encoding="utf-8") as f:
        for line in f:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9 or cols[2] != "gene":
                continue
            gid = re.search(r'gene_id "([^"]+)"', cols[8])
            gnm = re.search(r'gene_name "([^"]+)"', cols[8])
            if gid and gnm:
                ens = gid.group(1).split(".")[0]
                if ens not in mp:
                    mp[ens] = gnm.group(1)
    return mp

# FIX: Path assoluto per Windows
GTF_PATH = Path(r"C:\Users\Utente\Desktop\Tesi_Python_scRNA\data\gtf\gencode.v43.annotation.gtf.gz")
assert GTF_PATH.exists(), f"File non trovato: {GTF_PATH}"
gene_map = load_gtf_gene_map(GTF_PATH.as_posix())
print(f"Gene map caricata: {len(gene_map):,} geni.")

# 2) Carica e concatena tutti i campioni (.h5ad) usando dict + label='sample_id'
INPUT_DIR = "../work/01_qc_h5ad"
files = sorted(glob.glob(os.path.join(INPUT_DIR, "*.h5ad")))
print(f"Trovati {len(files)} campioni. Caricamento...")

adatas = {}
for f in files:
    sid = os.path.basename(f)[:-5]  # rimuove ".h5ad"
    a = sc.read_h5ad(f)
    if not a.var_names.is_unique:
        a.var_names_make_unique(join='-')
    adatas[sid] = a  # 'sid' sarà usato per costruire obs['sample_id'] da ad.concat

print("Concatenazione...")
# FIX: label="sample_id" crea automaticamente la colonna metadata corretta
adata = ad.concat(adatas, join="outer", label="sample_id",
                  fill_value=0, index_unique="-")
del adatas; gc.collect()
print(f"Merge completato. Forma iniziale: {adata.shape}")

# 3) Collasso geni duplicati (FIX CRITICO: somma colonne con suffisso -1, -2...)
if not adata.var_names.is_unique:
    print("Geni duplicati post-merge trovati. Collasso per somma...")
    base_ids = adata.var_names.to_series().str.replace(r"-\d+$", "", regex=True)
    groups = base_ids.groupby(base_ids).indices  # dict: base_id -> lista indici

    new_X = sparse.lil_matrix((adata.n_obs, len(groups)), dtype=adata.X.dtype)
    new_var_names = []
    for j, (gid, idxs) in enumerate(groups.items()):
        new_var_names.append(gid)
        s = adata.X[:, list(idxs)].sum(axis=1)
        if sparse.issparse(s):
            new_X[:, j] = s
        else:
            new_X[:, j] = sparse.csr_matrix(s).T

    adata = ad.AnnData(new_X.tocsr(), obs=adata.obs.copy())
    adata.var_names = pd.Index(new_var_names)
    print(f"Geni duplicati collassati. Nuova forma: {adata.shape}")

# 4) Ricostruzione .var (Ensembl, gene_symbol, mitocondriali)
adata.var["ensembl_id"]  = adata.var_names
adata.var["gene_symbol"] = adata.var_names.map(gene_map).astype("object")
sym = adata.var["gene_symbol"].fillna("").astype(str).str.upper()
adata.var["mt"] = sym.str.startswith("MT-")
adata.var_names_make_unique()

# 5) Mapping gruppi biologici (FIX: Regex e priorità per evitare errori carcinoma/oma)
def map_group(sample_id: str) -> str:
    sid = sample_id.upper()
    # Priorità 1: Carcinomi
    if "LGSOC" in sid or "LG-SOC" in sid: return "LGsOC"
    if "CCOC"  in sid or "CLEARCELL" in sid: return "CCOC"
    if "ENOC_S" in sid or "ENOC-S" in sid or "OVCA" in sid: return "EnOC_S"
    if re.search(r"\bENOC\b", sid): return "EnOC"
    # Priorità 2: Endometrioma (evita match parziale con carcinOMA)
    if "ENDOMETRIOMA" in sid or re.search(r"\bOMA\b", sid): return "OMA"
    # Priorità 3: Tessuti sani o eutopici
    if "EUTOPIC" in sid or "EUE" in sid: return "EuE"
    if "NORMAL" in sid or "NOV" in sid: return "nOV"
    if "ENDOMETRIOSIS" in sid or "ENDOMETRIOSI" in sid or "ENDO_" in sid: return "Endo"
    return "Unknown"

adata.obs["sample_id"] = adata.obs["sample_id"].astype("category")
adata.obs["princ_etic_1"] = adata.obs["sample_id"].map(map_group).astype("category")
print("Gruppi princ_etic_1:\n", adata.obs["princ_etic_1"].value_counts())

# 6) Filtro geni rari
print(f"Geni prima del filtro: {adata.n_vars}")
sc.pp.filter_genes(adata, min_cells=20)
print(f"Geni dopo filtro (>=20 cellule): {adata.n_vars}")

# 7) Standardizzazione Tipi & Layer
if not sparse.issparse(adata.X):
    adata.X = sparse.csr_matrix(adata.X)
adata.layers["counts"] = adata.X.astype(np.int32)

# 8) Sanity check + log
mapped_ratio = adata.var["gene_symbol"].notna().mean()
print(f"Gene symbols mappati: {mapped_ratio:.2%}")
adata.obs.groupby("sample_id").size().to_csv("../work/per_sample_counts.csv")
adata.obs.groupby("princ_etic_1").size().to_csv("../work/group_counts.csv")

# 9) Salvataggi finali
# Zarr è raccomandato per accessi futuri in modalità backed
adata.write_h5ad("../work/02_merged_raw.h5ad")
adata.write_zarr("../work/02_merged_raw.zarr", chunk_size=(10_000, 2_000))
print("\n--- FASE 2 COMPLETATA (Hardened) ---")
```

---

## 4. Prossimi Passaggi: Integrazione e Analisi (Fase 3)

Una volta creato il file `02_merged_raw.h5ad`, procederemo con l'analisi vera e propria.

**Punti chiave della Fase 3:**
1.  **HVG Batch-Aware:** Selezione dei 4000 geni più variabili considerando i campioni come batch (evita di selezionare geni "batch-specifici").
2.  **scVI Integration:** Integrazione su layer `counts` usando GPU (se disp.) o CPU. Questo sostituisce `Seurat::IntegrateData`.
3.  **UMAP Scalabile:** Invece di calcolare la UMAP su tutte le 390k cellule (crash RAM), calcoleremo la mappa su un subset di 80.000 cellule e proietteremo le restanti sopra di esse (*subset-and-transform*).
4.  **Clustering & DGE:** Leiden clustering e analisi differenziale robusta (con downsampling per velocità).

Il codice per la Fase 3 è già stato predisposto e sarà eseguito nella Cella 14.