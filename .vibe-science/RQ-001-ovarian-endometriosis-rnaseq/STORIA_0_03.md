Hai perfettamente ragione. Ho sbagliato a sintetizzare troppo nell'ultimo passaggio, pensando di fare un favore "riassuntivo", ma ho ottenuto l'effetto contrario perdendo i dettagli e la roadmap che avevamo faticosamente costruito. **Scusa.**

Rimedio subito.

Qui sotto trovi i **due file completi e non censurati**.
1.  **File Narrativa:** Contiene lo stato attuale ESATTO e, soprattutto, la **SCALETTA DEFINITIVA** dei prossimi step (4, 5, 6, ecc.) che abbiamo concordato e che dobbiamo ancora eseguire.
2.  **File Codice Master:** Contiene **TUTTE** le celle definitive, nell'ordine esatto di esecuzione, con tutte le patch (loader MTX flessibile, collasso duplicati sparso, consenso SOLO+Scrublet, Harmony).

Niente riassunti. Ecco la documentazione rigorosa.

---

### **FILE 1: NARRATIVA & ROADMAP DI PROGETTO**
*(Salva come `01_Narrativa_e_Roadmap.md`)*

# Analisi Bioinformatica scRNA-seq: Stato dell'Arte e Roadmap

**Data:** 1 Febbraio 2026
**Progetto:** Atlante integrato Endometriosi/Ovaio (52 dataset, >330k cellule).
**Obiettivo Tesi:** Identificare correlazioni molecolari tra endometriosi e carcinomi ovarici.

---

## PARTE 1: LAVORO SVOLTO (Methodology)

### 1. Ingestione Dati (Data Wrangling)
Abbiamo creato un workflow di ingestione custom per gestire l'eterogeneità dei 52 campioni GEO.
*   **Problema:** Formati misti (`.mtx` con prefissi non standard, `.h5`, `.csv`), file corrotti (float in file integer), dataset aggregati.
*   **Soluzione:**
    *   Script di **Split** per il dataset aggregato `GSE111976`.
    *   **Loader Universale (V6)**: Riconoscimento automatico del formato, parser `mtx` tollerante ai float, parser `csv` ottimizzato (engine C) per evitare timeout.
    *   **Standardizzazione:** Mappatura immediata a **Ensembl ID** (tramite GTF Gencode v43) per evitare ambiguità sui simboli genici.

### 2. Merge e Sanificazione
Dopo l'ingestione, i file `.h5ad` singoli sono stati concatenati.
*   **Problema:** Duplicazione dei geni (stesso Ensembl ID presente più volte) e perdita di metadati.
*   **Soluzione:** Algoritmo di **collasso sparso**: le colonne con lo stesso nome gene sono state sommate (non mediate) senza mai densificare la matrice in RAM.
*   **Esito:** Matrice pulita `(331.703 cellule x 33.596 geni)`.

### 3. Rimozione Doppietti (Strategia "Consenso")
Per garantire la massima pulizia per la tesi, abbiamo usato due metodi ortogonali:
1.  **SOLO (scVI-tools):** Approccio Deep Learning generativo.
2.  **Scrublet:** Approccio simulativo (kNN) eseguito *per-campione*.
*   **Regola di Filtro:** Abbiamo rimosso le cellule considerate doppietti da **entrambi** i metodi (intersezione calibrata) O quelle con punteggi di probabilità estremi (>0.90/0.95).
*   **File:** `02d_merged_nodoublets_consensus.h5ad` (321.232 cellule rimaste).

### 4. Normalizzazione e Riduzione Dimensionale
*   **Norm:** Scaling a 10.000 conteggi/cella + Log1p.
*   **HVG:** Selezione 2000 geni variabili (metodo Seurat v3) calcolata *per batch* (sample_id) per robustezza.
*   **PCA:** Calcolata sui soli HVG, senza centratura (`zero_center=False`) per efficienza computazionale.

### 5. Integrazione Batch (Harmony)
Per rimuovere l'effetto dei 52 diversi campioni senza distorcere la biologia, abbiamo applicato **Harmony** sulle prime 30 PC.
*   **Output:** UMAP e Clustering (Leiden res 0.5 e 1.0) calcolati sullo spazio armonizzato.
*   **File Corrente:** `../work/03_harmony_umap_leiden.h5ad`.

---

## PARTE 2: ROADMAP DEFINITIVA (Prossimi Step)

Questa è la scaletta concordata per completare l'analisi scientifica.

1.  **STEP 4: Annotazione Macro-Tipi (Marker-Based)**
    *   **Obiettivo:** Assegnare le classi principali (Epithelial, Stromal, Immune, Endothelial) usando pannelli di geni noti.
    *   **Metodo:** Calcolo `score_genes` per ogni firma. Assegnazione al cluster Leiden basata sullo score medio predominante. Verifica con dotplot.

2.  **STEP 5: Sub-clustering Epiteliale (Focus Tesi)**
    *   **Obiettivo:** Distinguere le popolazioni epiteliali sane (Ciliate, Secretorie) da quelle tumorali o intermedie.
    *   **Metodo:** Subset delle sole cellule epiteliali -> Nuova selezione HVG -> Nuova PCA/Harmony -> Clustering fine.
    *   **Marker:** FOXJ1/TPPP3 (Ciliate), PAX8/MUC16 (Secretorie), KRT7/MSLN.

3.  **STEP 6: Analisi Differenziale (Pseudobulk DGE)**
    *   **Obiettivo:** Trovare geni differenzialmente espressi tra Endometriosi (OMA/EuE) e Tumori.
    *   **Metodo:** Aggregazione dei conteggi ("pseudobulk") a livello di *Campione*.
    *   **Confronti:**
        *   OMA vs EuE (Eutopico)
        *   OMA vs CCOC/EnOC (Tumori)
    *   **Tool:** Test statistico robusto (T-test su LogCPM o PyDESeq2).

4.  **STEP 7: Pathway Enrichment (GSEA)**
    *   **Obiettivo:** Dare un significato biologico alle liste di geni (es. infiammazione, ciclo cellulare).
    *   **Tool:** `gseapy` su database GO_Biological_Process, KEGG, Reactome.

5.  **STEP 8: Traiettorie (Trajectory Inference)**
    *   **Obiettivo:** Vedere se esiste una transizione continua tra stati epiteliali (es. Secretorie -> Ciliate o Normale -> Tumorale).
    *   **Tool:** PAGA (Partition-based Graph Abstraction) + DPT (Diffusion Pseudotime).

6.  **STEP 9: Interazioni Cellula-Cellula (Opzionale/Premium)**
    *   **Obiettivo:** Capire come lo stroma o il sistema immune parlano con l'epitelio.
    *   **Tool:** LIANA o Squidpy.

---
---

### **FILE 2: MASTER CODE (Notebook Completo)**
*(Salva come `02_Master_Pipeline_Code.py` o incolla nelle celle di Jupyter)*

Questo file contiene **tutto** il codice validato fino ad ora, incluse le funzioni helper V6 corrette.

#### **CELLA 1: LIBRERIE E HELPER V6 (Loader, QC, Fix)**
*Contiene il loader MTX con fallback float-safe e la gestione duplicati.*

```python
# === CELLA 1: HELPER V6 (DEFINITIVI) ===
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

# --- 1. GTF Loader ---
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

# --- 2. Utils AnnData ---
def _make_names_unique(a: ad.AnnData) -> ad.AnnData:
    a.obs_names = pd.Index(a.obs_names.astype(str))
    a.var_names = pd.Index(a.var_names.astype(str))
    if not a.obs_names.is_unique: a.obs_names_make_unique()
    if not a.var_names.is_unique: a.var_names_make_unique()
    return a

def _only_gene_expression(a: ad.AnnData) -> ad.AnnData:
    if "feature_types" in a.var.columns:
        mask = a.var["feature_types"].astype(str).str.contains("Gene Expression", case=False, na=False)
        if mask.any(): a = a[:, mask].copy()
    return a

def _dedup_columns(cols: pd.Index) -> pd.Index:
    base = cols.astype(str); seen = {}; out = []
    for x in base:
        if x not in seen: seen[x] = 0; out.append(x)
        else: seen[x] += 1; out.append(f"{x}.{seen[x]}")
    return pd.Index(out)

# --- 3. Reader MTX Float-Safe ---
def _read_mtx_flexible(path: Path):
    opener = gzip.open if str(path).endswith(".gz") else open
    header_lines = 0
    with opener(path, "rt", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.startswith('%'): header_lines += 1; continue
            parts = line.strip().split()
            if len(parts) < 3: raise ValueError("Header MTX malformato")
            n_rows, n_cols, _ = map(int, parts[:3])
            header_lines += 1; break
    df = pd.read_csv(
        path, sep=r"\s+", engine="python", comment="%", header=None,
        names=["i","j","x"], skiprows=header_lines
    ).dropna(subset=["i","j","x"])
    i = df["i"].astype(np.int64).to_numpy() - 1
    j = df["j"].astype(np.int64).to_numpy() - 1
    x = pd.to_numeric(df["x"], errors="coerce").fillna(0).astype(np.float32).to_numpy()
    return coo_matrix((x, (i, j)), shape=(n_rows, n_cols)).tocsr()

# --- 4. Loader 10x MTX (con prefissi e fallback) ---
def _load_10x_mtx(p: Path) -> ad.AnnData | None:
    candidates = [p] + [d for d in p.glob("filtered_feature_bc_matrix*") if d.is_dir()]
    for base in candidates:
        matrix_file = next(base.glob("*matrix.mtx*"), None)
        barcodes_file = next(base.glob("*barcodes.tsv*"), None)
        features_file = next(base.glob("*features.tsv*"), None) or next(base.glob("*genes.tsv*"), None)
        if not (matrix_file and barcodes_file and features_file): continue
        
        try: M = mmread(str(matrix_file)).tocsr()
        except ValueError as e:
            if "invalid integer" in str(e).lower(): M = _read_mtx_flexible(matrix_file)
            else: raise

        barcodes = pd.read_csv(barcodes_file, sep="\t", header=None, usecols=[0], dtype=str).iloc[:,0].astype(str)
        barcodes = _dedup_columns(pd.Index(barcodes))
        feat = pd.read_csv(features_file, sep="\t", header=None, dtype=str).fillna("")
        
        if feat.shape[1] >= 3:
            gene_id, gene_nm, ftype = feat.iloc[:,0], feat.iloc[:,1], feat.iloc[:,2]
        else:
            gene_id = feat.iloc[:,0]; gene_nm = gene_id.copy(); ftype = pd.Series(["Gene Expression"]*len(gene_id))

        var_names = np.where(gene_id.str.upper().str.startswith("ENSG"), gene_id, gene_nm).astype(str)
        n_cells = min(M.shape[1], len(barcodes)); n_genes = min(M.shape[0], len(var_names))
        M = M[:n_genes, :n_cells].tocsr()
        
        a = ad.AnnData(M.T.tocsr(), obs=pd.DataFrame(index=barcodes[:n_cells]), var=pd.DataFrame(index=var_names[:n_genes]))
        a.var["gene_id_raw"] = gene_id.values[:n_genes]
        a.var["feature_types"] = ftype.values[:n_genes]
        a.var_names = a.var_names.str.split(".").str[0]
        return _make_names_unique(_only_gene_expression(a))
    return None

# --- 5. Loader H5 ---
def _load_10x_h5(p: Path) -> ad.AnnData | None:
    h5 = next(p.glob("*.h5"), None)
    if h5:
        a = sc.read_10x_h5(h5.as_posix())
        if a.var_names[0].startswith("ENSG"): a.var_names = a.var_names.str.split(".").str[0]
        return _make_names_unique(_only_gene_expression(a))
    return None

# --- 6. Loader Table (CSV/TXT) ---
def _load_table(p: Path) -> ad.AnnData | None:
    files = list(p.glob("*.csv*")) + list(p.glob("*.tsv*")) + list(p.glob("*.txt*"))
    if not files: return None
    f = files[0]
    with open(f, "r") as fh: sep = "\t" if fh.read(2000).count("\t") > 5 else ","
    try: df = pd.read_csv(f, sep=sep, engine="c", index_col=0, dtype=np.float32)
    except: df = pd.read_csv(f, sep=None, engine="python", index_col=0).apply(pd.to_numeric, errors='coerce').fillna(0)
    
    df = df.loc[df.index.notna(), :]; df = df.loc[:, df.columns.notna()]
    if not df.index.is_unique: df = df.groupby(level=0).sum()
    if not df.columns.is_unique: df.columns = _dedup_columns(df.columns)
    
    if df.shape[0] >= df.shape[1]: X = sparse.csr_matrix(df.T.values); var = df.index; obs = df.columns
    else: X = sparse.csr_matrix(df.values); var = df.columns; obs = df.index
    return _make_names_unique(ad.AnnData(X, obs=pd.DataFrame(index=obs.astype(str)), var=pd.DataFrame(index=var.astype(str))))

def load_one_sample(sample_dir: Path) -> ad.AnnData:
    a = _load_10x_mtx(sample_dir); 
    if a is None: a = _load_10x_h5(sample_dir)
    if a is None: a = _load_table(sample_dir)
    if a is None: raise FileNotFoundError(f"Nessun formato in {sample_dir}")
    if not sparse.issparse(a.X): a.X = sparse.csr_matrix(a.X)
    return _make_names_unique(a)

# --- 7. Mapping & QC ---
def to_ensembl_and_collapse(a: ad.AnnData, gene_map: dict, symbol_to_ens: dict) -> ad.AnnData:
    vn = pd.Index(a.var_names.astype(str)); is_ens = vn.str.upper().str.startswith("ENSG")
    new_vn = vn.to_numpy(dtype=object); mapped = vn[~is_ens].map(symbol_to_ens)
    new_vn[~is_ens] = mapped.where(pd.notna(mapped), vn[~is_ens])
    a.var_names = pd.Index(new_vn.astype(str))
    
    if not a.var_names.is_unique:
        codes, uniques = pd.factorize(a.var_names, sort=False)
        G = sparse.csr_matrix((np.ones_like(codes), (codes, np.arange(len(codes)))), shape=(len(uniques), len(codes)))
        a = ad.AnnData((a.X @ G.T).tocsr(), obs=a.obs.copy())
        a.var_names = pd.Index(uniques.astype(str))
        
    vn = pd.Index(a.var_names.astype(str))
    a.var = pd.DataFrame(index=vn)
    a.var["gene_symbol"] = pd.Series([gene_map.get(g, g) for g in vn], index=vn, dtype="object")
    a.var["mt"] = a.var["gene_symbol"].str.upper().str.startswith("MT-")
    return _make_names_unique(a)

def apply_qc(a: ad.AnnData, mt_pct=20.0, min_genes=200, min_cells=3) -> ad.AnnData:
    sc.pp.calculate_qc_metrics(a, qc_vars=["mt"], inplace=True)
    sc.pp.filter_cells(a, min_genes=min_genes); sc.pp.filter_genes(a, min_cells=min_cells)
    if "pct_counts_mt" in a.obs: a = a[a.obs["pct_counts_mt"] < mt_pct].copy()
    a.layers["counts"] = a.X.astype(np.int32)
    return _make_names_unique(a)
```

#### **CELLA 2: CONFIGURAZIONE & INGESTIONE (Loop)**

```python
# === CELLA 2: CONFIG E INGESTIONE ===
RAW_DIR = Path("../data/raw"); OUT_DIR = Path("../work/01_qc_h5ad"); OUT_DIR.mkdir(parents=True, exist_ok=True)
GTF_PATH = Path(r"C:\Users\Utente\Desktop\Tesi_Python_scRNA\data\gtf\gencode.v43.annotation.gtf.gz")

gene_map = load_gtf_gene_map(GTF_PATH)
symbol_to_ens = {v: k for k, v in gene_map.items()}

EXCLUDE = {"GSE111976", "GSE111976_processed_temp", "_ignore"}
folders = sorted([p for p in RAW_DIR.iterdir() if p.name not in EXCLUDE and has_known_data(p)])
gse_split = sorted(RAW_DIR.glob("eutopic_*_GSM*"))
all_targets = sorted(list(set(folders + gse_split)), key=lambda x: x.name)

print(f"Campioni da processare: {len(all_targets)}")

for sample_dir in tqdm(all_targets):
    sid = sample_dir.name; out = OUT_DIR / f"{sid}.h5ad"
    if out.exists(): continue
    try:
        a = load_one_sample(sample_dir)
        a = to_ensembl_and_collapse(a, gene_map, symbol_to_ens)
        a = apply_qc(a)
        a.obs["sample_id"] = sid
        a.write_h5ad(out.as_posix(), compression="gzip")
        print(f"OK: {sid} {a.shape}")
    except Exception as e: print(f"FAIL {sid}: {e}")
```

#### **CELLA 3: MERGE E SANIFICAZIONE FINALE**

```python
# === CELLA 3: MERGE E COLLASSO DUPLICATI ===
files = sorted(OUT_DIR.glob("*.h5ad"))
adatas = {f.stem: sc.read_h5ad(f) for f in files}
adata = ad.concat(adatas, join="outer", label="sample_id", fill_value=0, index_unique="-")
print("Merged shape:", adata.shape)

# Collasso Finale Duplicati (Sparse)
adata.var_names = adata.var_names.astype(str).str.split(".").str[0]
if not adata.var_names.is_unique:
    codes, uniques = pd.factorize(adata.var_names, sort=False)
    G = sparse.csr_matrix((np.ones_like(codes), (codes, np.arange(len(codes)))), shape=(len(uniques), len(codes)))
    adata = ad.AnnData((adata.X @ G.T).tocsr(), obs=adata.obs.copy())
    adata.var_names = pd.Index(uniques.astype(str))

# Ricostruzione .var e counts
vn = pd.Index(adata.var_names.astype(str))
adata.var = pd.DataFrame(index=vn)
adata.var["gene_symbol"] = pd.Series([gene_map.get(g, g) for g in vn], index=vn, dtype="object")
adata.var["mt"] = adata.var["gene_symbol"].str.upper().str.startswith("MT-")
adata.layers["counts"] = adata.X.astype(np.int32)

adata.write_h5ad("../work/02_merged_raw.h5ad", compression="gzip")
print("Saved 02_merged_raw.h5ad")
```

#### **CELLA 4: DOUBLET CONSENSUS (SOLO + SCRUBLET)**

```python
# === CELLA 4: DOUBLET CONSENSUS ===
# Nota: richiede scvi-tools, scrublet
IN_MERGED = "../work/02_merged_raw.h5ad"
OUT_NODBL = "../work/02d_merged_nodoublets_consensus.h5ad"

adata = sc.read_h5ad(IN_MERGED)

# 1. SOLO
import scvi
from scvi.external import SOLO
scvi.model.SCVI.setup_anndata(adata, layer="counts", batch_key="sample_id")
vae = scvi.model.SCVI(adata, n_latent=20)
vae.train(max_epochs=40, batch_size=2048, accelerator="cpu") # Aumenta epoch per finale
solo = SOLO.from_scvi_model(vae)
solo.train(max_epochs=20, batch_size=2048, accelerator="cpu")
preds = solo.predict(soft=True)
adata.obs["doublet_score_solo"] = preds[:, 1] if preds.ndim==2 else preds

# 2. Scrublet (Per Sample)
adata.obs["predicted_doublet_scr"] = False
for s in adata.obs["sample_id"].unique():
    idx = adata.obs["sample_id"] == s
    if idx.sum() < 200: continue
    sub = adata[idx].copy(); sub.X = sub.layers["counts"]
    sc.pp.scrublet(sub, verbose=False)
    adata.obs.loc[idx, "predicted_doublet_scr"] = sub.obs["predicted_doublet"]

# 3. Consenso e Filtro
# Soglia SOLO dinamica (quantile 95% per sample)
solo_thr = adata.obs.groupby("sample_id")["doublet_score_solo"].transform(lambda x: np.quantile(x, 0.95))
solo_flag = (adata.obs["doublet_score_solo"] > np.maximum(solo_thr, 0.5))
cons = (solo_flag & adata.obs["predicted_doublet_scr"]) | (adata.obs["doublet_score_solo"] > 0.95)

adata = adata[~cons].copy()
adata.write_h5ad(OUT_NODBL, compression="gzip")
print(f"Saved Nodoublets: {adata.shape}")
```

#### **CELLA 5: ANALISI BASE (Norm, HVG, PCA, Harmony, UMAP)**

```python
# === CELLA 5: ANALISI BASE & INTEGRAZIONE ===
import scanpy.external as sce

IN = "../work/02d_merged_nodoublets_consensus.h5ad"
OUT = "../work/03_harmony_umap_leiden.h5ad"

adata = sc.read_h5ad(IN)

# Norm & Log
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata

# HVG & PCA
sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=2000, batch_key="sample_id")
sc.tl.pca(adata, n_comps=30, use_highly_variable=True, zero_center=False)

# Harmony Integration
sce.pp.harmony_integrate(adata, key="sample_id", basis="X_pca", adjusted_basis="X_pca_harmony")

# Embedding & Clustering
sc.pp.neighbors(adata, use_rep="X_pca_harmony", n_neighbors=15)
sc.tl.umap(adata, min_dist=0.5)
sc.tl.leiden(adata, resolution=0.5, key_added="leiden05")
sc.tl.leiden(adata, resolution=1.0, key_added="leiden10")

adata.write_h5ad(OUT, compression="gzip")
print("Analisi Base Completata e Salvata.")
```