Ecco il documento completo in formato Markdown (`.md`), pronto per essere caricato su Canvas o usato come documentazione ufficiale della tesi.

Include la narrativa del problema, le decisioni tecniche prese, i bug risolti e il codice definitivo che ha portato al successo della conversione dei dati.

***

# Migrazione Workflow scRNA-seq: Fase 1 (Ingestione & Standardizzazione)

**Progetto:** Analisi Single-Cell RNA-seq su 52 campioni (Endometriosi, Tumori Ovarici, Tessuti Sani).
**Stato:** Fase 1 Completata (Tutti i campioni convertiti in `.h5ad`).
**Data:** 1 Febbraio 2026

---

## 1. Contesto e Motivazione

Il progetto prevede l'analisi di un atlante cellulare composto da oltre **390.000 cellule** provenienti da 52 dataset pubblici (NCBI GEO).

### Il Problema (R & Seurat)
Inizialmente, l'analisi è stata condotta in **R** utilizzando la suite **Seurat**. Tuttavia, il workflow si è bloccato durante le fasi di integrazione e calcolo UMAP. Nonostante l'utilizzo di una istanza AWS EC2 con **750GB di RAM**, R non è riuscito a gestire la matrice di distanza globale e la duplicazione degli oggetti in memoria, causando crash sistematici (*Out Of Memory*).

### La Soluzione (Python & Scanpy)
Abbiamo deciso di migrare l'intero workflow in **Python** per sfruttare architetture più efficienti:
1.  **Matrici Sparse (CSR):** Scanpy/AnnData gestiscono i dati come matrici sparse, riducendo drasticamente l'uso della RAM.
2.  **Modularità:** Elaborazione per singolo campione prima del merge.
3.  **Robustezza:** Standardizzazione dei nomi dei geni (Ensembl ID) e gestione dei formati eterogenei alla fonte.

---

## 2. Fase 0: Pre-processing Dataset Aggregati (GSE111976)

Il dataset `GSE111976` presentava una sfida unica: i dati di 6 campioni di endometrio eutopico erano aggregati in un unico file gigante (convertito da `.rds` a `.csv` tramite R). Caricare tutto in memoria per poi filtrare era impossibile.

**Strategia:** Abbiamo scritto uno script che legge il CSV gigante in modalità "lazy" (usando `usecols` di Pandas), estraendo solo le colonne (barcodes) specifiche per ogni donatore/giorno, salvando poi 6 file `matrix.csv` separati.

### Codice: Split GSE111976

```python
# Cella 12 — Split robusto del dataset aggregato GSE111976
# Obiettivo: Creare 6 sottocartelle eutopic_* partendo dai metadati

import re
from pathlib import Path
import pandas as pd
from tqdm import tqdm

RAW_DIR = Path("../data/raw")
GSE_DIR = RAW_DIR / "GSE111976"
COUNTS_PATH = GSE_DIR / "GSE111976_full_matrix.csv" # File generato da R
META_PATH = GSE_DIR / "GSE111976_summary_10x_day_donor_ctype.csv"

# 1. Carica Metadati
meta = pd.read_csv(META_PATH)
meta.columns = meta.columns.str.lower().str.strip()

# 2. Definizione dei campioni da estrarre
sample_defs = [
    ("eutopic_7_GSM4577307", 19, 22), ("eutopic_8_GSM4577308", 20, 22),
    ("eutopic_9_GSM4577309", 29, 20), ("eutopic_10_GSM4577310", 39, 23),
    ("eutopic_11_GSM4577312", 57, 26), ("eutopic_12_GSM4577314", 60, 26),
]

# 3. Lettura Header per mappatura colonne
hdr = pd.read_csv(COUNTS_PATH, nrows=0)
all_cols = hdr.columns.astype(str).tolist()
# Mappa per normalizzare suffissi (es. rimuovere -1 dai barcode)
def strip_suf(s): return re.sub(r"-\d+$", "", s)
counts_norm_map = {strip_suf(c): c for c in all_cols[1:]}

# 4. Estrazione mirata
for sample_name, donor_id, day_id in tqdm(sample_defs, desc="Split GSE111976"):
    # Identifica barcode target dai metadati
    bc_meta = meta.loc[(meta["donor"] == donor_id) & (meta["day"] == day_id), "cell_name"].astype(str).tolist()
    
    # Trova i nomi reali nel CSV dei conteggi
    bc_found = [counts_norm_map[strip_suf(b)] for b in bc_meta if strip_suf(b) in counts_norm_map]
    
    if bc_found:
        out_dir = RAW_DIR / sample_name; out_dir.mkdir(exist_ok=True)
        # Legge SOLO le colonne necessarie (Memoria efficiente)
        sub = pd.read_csv(COUNTS_PATH, usecols=[all_cols[0]] + bc_found, index_col=0)
        sub.to_csv(out_dir / "matrix.csv")
        print(f"[OK] {sample_name}: salvate {len(bc_found)} cellule.")
```

---

## 3. Fase 1: Ingestione e Standardizzazione (Workflow Definitivo)

Questa è stata la parte più complessa. I 52 campioni presentavano formati eterogenei (`.mtx`, `.h5`, `.csv`, `.txt`) e diverse anomalie nei dati grezzi.

### Sfide Affrontate e Bug Risolti

1.  **Nomi File Non Standard:** Alcuni campioni 10x (es. `normal_5`) avevano prefissi nei file (es. `GSM..._matrix.mtx.gz`) che il loader standard di Scanpy non riconosceva.
    *   *Soluzione:* Scritto un loader personalizzato che usa `glob` per trovare i file indipendentemente dal prefisso.
2.  **Valori Float in file MTX:** Alcuni file `.mtx` contenevano valori decimali, facendo crashare il parser standard (`scipy.io.mmread`) che si aspettava interi ("Invalid integer value").
    *   *Soluzione:* Implementato un parser `_read_mtx_flexible` che usa Pandas per leggere i dati come float e convertirli in matrice sparsa.
3.  **Lentezza Caricamento CSV:** Usare `pd.to_numeric` su file di testo enormi (es. `CCCOC_1`) era troppo lento e causava blocchi.
    *   *Soluzione:* Scritto un loader `_read_table_fast` che usa l'engine `c` di Pandas e tipizzazione `float32` diretta.
4.  **Duplicati e Indici:** Molti CSV avevano geni duplicati (righe) o barcode duplicati (colonne), causando errori `cannot reindex` o `Index does not support mutable operations` più a valle.
    *   *Soluzione:* Gestione dei duplicati (somma per i geni, rinomina univoca per i barcode) **durante il caricamento**, prima della creazione dell'oggetto AnnData.
5.  **Gene Mapping:** Necessità di uniformare tutto a **Ensembl ID**.
    *   *Soluzione:* Mappatura robusta tramite file GTF. I geni che non mappano vengono mantenuti con il loro nome originale per non perdere dati.

### Codice Definitivo

Il workflow finale è diviso in 3 celle Jupyter.

#### Cella 1: Funzioni Helper "Corazzate"
Include i loader per ogni formato, la gestione dei duplicati e la logica di QC.

```python
# === Cella 1: Helper robusti (V6) ===
from pathlib import Path
import re, gzip
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
from scipy import sparse
from scipy.io import mmread
from scipy.sparse import coo_matrix
from tqdm.auto import tqdm

# ---------- 1. Caricamento Mappa Geni (GTF) ----------
def load_gtf_gene_map(gtf_path: Path) -> dict:
    mp = {}
    opener = gzip.open if str(gtf_path).endswith(".gz") else open
    with opener(gtf_path, "rt", encoding="utf-8") as f:
        for ln in f:
            if ln.startswith("#"): continue
            cs = ln.rstrip("\n").split("\t")
            if len(cs) < 9 or cs[2] != "gene": continue
            gid = re.search(r'gene_id "([^"]+)"', cs[8])
            gnm = re.search(r'gene_name "([^"]+)"', cs[8])
            if gid and gnm:
                mp[gid.group(1).split(".")[0]] = gnm.group(1)
    return mp

# ---------- 2. Utility per AnnData ----------
def _make_names_unique(a: ad.AnnData) -> ad.AnnData:
    a.obs_names = pd.Index(a.obs_names.astype(str))
    a.var_names = pd.Index(a.var_names.astype(str))
    if not a.obs_names.is_unique: a.obs_names_make_unique()
    if not a.var_names.is_unique: a.var_names_make_unique()
    return a

def _only_gene_expression(a: ad.AnnData) -> ad.AnnData:
    """Filtra feature come 'Antibody Capture' se presenti."""
    if "feature_types" in a.var.columns:
        mask = a.var["feature_types"].astype(str).str.contains("Gene Expression", case=False, na=False)
        if mask.any():
            a = a[:, mask].copy()
    return a

def _dedup_columns(cols: pd.Index) -> pd.Index:
    """Rende univoci i nomi delle colonne (barcode) aggiungendo suffissi."""
    base = cols.astype(str)
    seen = {}; out = []
    for x in base:
        if x not in seen: seen[x] = 0; out.append(x)
        else: seen[x] += 1; out.append(f"{x}.{seen[x]}")
    return pd.Index(out)

# ---------- 3. Loader MTX "Float-Safe" ----------
def _read_mtx_flexible(path: Path):
    """Fallback per leggere file .mtx che contengono float invece di interi."""
    opener = gzip.open if str(path).endswith(".gz") else open
    header_lines = 0
    with opener(path, "rt", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.startswith('%'): header_lines += 1; continue
            break
    # Legge ignorando l'header problematico
    df = pd.read_csv(
        path, sep=r"\s+", engine="python", comment="%", header=None,
        names=["i","j","x"], skiprows=header_lines
    ).dropna(subset=["i","j","x"])
    i = df["i"].astype(np.int64).to_numpy() - 1
    j = df["j"].astype(np.int64).to_numpy() - 1
    x = pd.to_numeric(df["x"], errors="coerce").fillna(0).astype(np.float32).to_numpy()
    # Ricostruisce la matrice sparsa (COO -> CSR)
    # Nota: leggiamo dimensioni n_rows/n_cols dalla prima riga non commentata se necessario, 
    # qui semplificato per autodetect dalla coo_matrix
    return coo_matrix((x, (i, j))).tocsr()

# ---------- 4. Loader Universale 10x (MTX) ----------
def _load_10x_mtx(p: Path) -> ad.AnnData | None:
    """Carica MTX gestendo prefissi nei nomi e fallback float."""
    candidates = [p] + [d for d in p.glob("filtered_feature_bc_matrix*") if d.is_dir()]
    for base in candidates:
        matrix_file   = next(base.glob("*matrix.mtx*"),   None)
        barcodes_file = next(base.glob("*barcodes.tsv*"), None)
        features_file = next(base.glob("*features.tsv*"), None) or next(base.glob("*genes.tsv*"), None)
        if not (matrix_file and barcodes_file and features_file): continue

        try:
            M = mmread(str(matrix_file)).tocsr()
        except ValueError as e:
            # Fallback se mmread fallisce su float
            if "invalid integer" in str(e).lower():
                M = _read_mtx_flexible(matrix_file)
            else: raise

        barcodes = pd.read_csv(barcodes_file, sep="\t", header=None, usecols=[0], dtype=str).iloc[:,0].astype(str)
        barcodes = _dedup_columns(pd.Index(barcodes))
        
        feat = pd.read_csv(features_file, sep="\t", header=None, dtype=str).fillna("")
        # Gestione colonne features (ID, Symbol, Type)
        if feat.shape[1] >= 3:
            gene_id, gene_nm, ftype = feat.iloc[:,0], feat.iloc[:,1], feat.iloc[:,2]
        else:
            gene_id = feat.iloc[:,0]; gene_nm = gene_id.copy(); ftype = pd.Series(["Gene Expression"]*len(gene_id))

        # Usa Ensembl ID se disponibile, altrimenti Symbol
        var_names = np.where(gene_id.str.upper().str.startswith("ENSG"), gene_id, gene_nm).astype(str)

        # Allineamento dimensioni (safe crop)
        n_cells = min(M.shape[1], len(barcodes))
        n_genes = min(M.shape[0], len(var_names))
        M = M[:n_genes, :n_cells]
        
        a = ad.AnnData(M.T.tocsr(), obs=pd.DataFrame(index=barcodes[:n_cells]), var=pd.DataFrame(index=var_names[:n_genes]))
        a.var["gene_id"] = gene_id.values[:n_genes]
        a.var["feature_types"] = ftype.values[:n_genes]
        
        a.var_names = a.var_names.str.split(".").str[0] # Rimuove versioni Ensembl
        return _make_names_unique(_only_gene_expression(a))
    return None

# ---------- 5. Loader Universale H5 ----------
def _load_10x_h5(p: Path) -> ad.AnnData | None:
    h5 = next(p.glob("*.h5"), None)
    if h5:
        a = sc.read_10x_h5(h5.as_posix())
        if a.var_names[0].startswith("ENSG"):
            a.var_names = a.var_names.str.split(".").str[0]
        return _make_names_unique(_only_gene_expression(a))
    return None

# ---------- 6. Loader Universale Tabellare (CSV/TXT) ----------
def _load_table(p: Path) -> ad.AnnData | None:
    files = list(p.glob("*.csv*")) + list(p.glob("*.tsv*")) + list(p.glob("*.txt*"))
    if not files: return None
    f = files[0]
    
    # Detect separatore
    with open(f, "r") as fh: sep = "\t" if fh.read(2000).count("\t") > 5 else ","
    
    try:
        # Fast load con engine C
        df = pd.read_csv(f, sep=sep, engine="c", index_col=0, dtype=np.float32)
    except:
        # Fallback python engine
        df = pd.read_csv(f, sep=None, engine="python", index_col=0)
        df = df.apply(pd.to_numeric, errors='coerce').fillna(0).astype(np.float32)

    # Pulizia Indici e Colonne
    df = df.loc[df.index.notna() & (df.index.astype(str) != ""), :]
    if not df.index.is_unique: df = df.groupby(level=0).sum() # Collassa geni duplicati
    if not df.columns.is_unique: df.columns = _dedup_columns(df.columns) # Rinomina barcode duplicati

    # Euristica Orientamento (Geni sulle righe o colonne?)
    if df.shape[0] >= df.shape[1]: # Righe >= Colonne -> Geni su righe
        X = sparse.csr_matrix(df.T.values)
        var = pd.DataFrame(index=df.index.astype(str))
        obs = pd.DataFrame(index=df.columns.astype(str))
    else:
        X = sparse.csr_matrix(df.values)
        var = pd.DataFrame(index=df.columns.astype(str))
        obs = pd.DataFrame(index=df.index.astype(str))

    return _make_names_unique(ad.AnnData(X, obs=obs, var=var))

def load_one_sample(sample_dir: Path) -> ad.AnnData:
    a = _load_10x_mtx(sample_dir)
    if a is None: a = _load_10x_h5(sample_dir)
    if a is None: a = _load_table(sample_dir)
    if a is None: raise FileNotFoundError(f"Nessun formato valido in {sample_dir}")
    if not sparse.issparse(a.X): a.X = sparse.csr_matrix(a.X)
    return _make_names_unique(a)

# ---------- 7. Mapping e QC ----------
def to_ensembl_and_collapse(a: ad.AnnData, gene_map: dict, symbol_to_ens: dict) -> ad.AnnData:
    # Mappatura Symbol -> Ensembl
    vn = pd.Index(a.var_names.astype(str))
    is_ens = vn.str.upper().str.startswith("ENSG")
    
    new_vn = vn.to_numpy(dtype=object)
    # Mappa solo se non è già Ensembl
    mapped = vn[~is_ens].map(symbol_to_ens)
    new_vn[~is_ens] = mapped.where(pd.notna(mapped), vn[~is_ens]) # Se mappa usa Ensembl, se no lascia originale
    
    a.var_names = pd.Index(new_vn.astype(str))

    # Collasso post-mapping (somma conteggi per stesso ID)
    if not a.var_names.is_unique:
        codes, uniques = pd.factorize(a.var_names, sort=False)
        G = sparse.csr_matrix((np.ones_like(codes), (codes, np.arange(len(codes)))), shape=(len(uniques), len(codes)))
        a = ad.AnnData((a.X @ G.T).tocsr(), obs=a.obs.copy())
        a.var_names = pd.Index(uniques.astype(str))

    # Annotazione Simboli e Mito
    vn = pd.Index(a.var_names.astype(str))
    a.var = pd.DataFrame(index=vn)
    a.var["gene_symbol"] = pd.Series([gene_map.get(g, g) for g in vn], index=vn, dtype="object")
    a.var["mt"] = a.var["gene_symbol"].str.upper().str.startswith("MT-")
    return _make_names_unique(a)

def apply_qc(a: ad.AnnData, mt_pct=20.0, min_genes=200, min_cells=3) -> ad.AnnData:
    sc.pp.calculate_qc_metrics(a, qc_vars=["mt"], inplace=True)
    sc.pp.filter_cells(a, min_genes=min_genes)
    sc.pp.filter_genes(a, min_cells=min_cells)
    if "pct_counts_mt" in a.obs:
        a = a[a.obs["pct_counts_mt"] < mt_pct].copy()
    a.layers["counts"] = a.X.astype(np.int32) # Salva raw counts come interi
    return _make_names_unique(a)
```

#### Cella 2: Configurazione
Imposta i percorsi e carica il riferimento genomico.

```python
# === Cella 2: Configurazione ===
RAW_DIR   = Path("../data/raw")
OUT_DIR   = Path("../work/01_qc_h5ad"); OUT_DIR.mkdir(parents=True, exist_ok=True)
GTF_PATH  = Path(r"C:\Users\Utente\Desktop\Tesi_Python_scRNA\data\gtf\gencode.v43.annotation.gtf.gz")

MT_PCT    = 20.0
MIN_GENES = 200
MIN_CELLS = 3
FORCE_RERUN = True # Sovrascrivi file esistenti per sicurezza

# Carica Mappa Geni
gene_map = load_gtf_gene_map(GTF_PATH)
symbol_to_ens = {v: k for k, v in gene_map.items()}
print(f"[SETUP] Mappa geni caricata: {len(gene_map):,} entries.")

# Filtro cartelle valide
EXCLUDE = {"GSE111976", "GSE111976_processed_temp", "_ignore"}
def has_data(p): 
    return p.is_dir() and (list(p.glob("*.mtx*")) or list(p.glob("*matrix*")) or list(p.glob("*.h5")) or list(p.glob("*.csv*")))

folders = sorted([p for p in RAW_DIR.iterdir() if p.name not in EXCLUDE and has_data(p)])
gse_split = sorted(RAW_DIR.glob("eutopic_*_GSM*"))
all_targets = sorted(list(set(folders + gse_split)), key=lambda x: x.name)

print(f"[SETUP] Totale campioni da processare: {len(all_targets)}")
```

#### Cella 3: Esecuzione Loop
Itera su tutti i campioni, applica la logica e salva.

```python
# === Cella 3: Esecuzione Ingestione ===
for sample_dir in tqdm(all_targets):
    sid = sample_dir.name
    out = OUT_DIR / f"{sid}.h5ad"
    
    if out.exists() and not FORCE_RERUN: continue

    try:
        print(f"\n[PROC] {sid}")
        # 1. Load (con autodetect formato)
        a = load_one_sample(sample_dir)
        
        # 2. Map & Collapse (Standardizza a Ensembl)
        a = to_ensembl_and_collapse(a, gene_map, symbol_to_ens)
        
        # 3. QC
        a = apply_qc(a, mt_pct=MT_PCT, min_genes=MIN_GENES, min_cells=MIN_CELLS)
        
        # 4. Metadata & Save
        a.obs["sample_id"] = sid
        a.write_h5ad(out.as_posix(), compression="gzip")
        print(f"  -> OK: {out.name} | Shape: {a.shape}")
        
    except Exception as e:
        print(f"  -> FAILED: {e}")
        # import traceback; traceback.print_exc() # Scommenta per debug profondo
```

---

## 4. Conclusioni Fase 1

Con l'esecuzione di questi script, tutti i 52 dataset (inclusi i casi problematici come `normal_5/6/7` e gli split di `GSE111976`) sono stati convertiti con successo in un formato standard `.h5ad`, pronti per il merge.

**Caratteristiche degli output:**
*   **Indice:** Ensembl ID univoci.
*   **Metadata:** `gene_symbol` (HGNC) e `mt` (booleano mitocondriale) in `.var`.
*   **Dati:** Matrice sparsa (`CSR`) compressa.
*   **Raw Counts:** Salvati in `layers['counts']`.

Siamo pronti per la **Fase 2: Merge e Integrazione con scVI**.