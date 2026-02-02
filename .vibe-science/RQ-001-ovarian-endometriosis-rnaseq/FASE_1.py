# ==============================================================================
# PIPELINE DI INGESTIONE scRNA-seq (FASE 1)
# ==============================================================================
# Descrizione: 
# Script robusto per convertire dati raw eterogenei (10x MTX, H5, CSV/TXT) 
# in oggetti AnnData (.h5ad) standardizzati.
# Include gestione per:
# - Nomi file non standard (prefissi)
# - MTX con valori float (fallback reader)
# - Mapping Geni Ensembl/HGNC
# - QC di base
# ==============================================================================

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

# ------------------------------------------------------------------------------
# 1. FUNZIONI HELPER DI CARICAMENTO E UTILITIES
# ------------------------------------------------------------------------------

def load_gtf_gene_map(gtf_path: Path) -> dict:
    """Crea un dizionario Ensembl ID -> Gene Symbol dal file GTF."""
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

def _make_names_unique(a: ad.AnnData) -> ad.AnnData:
    """Garantisce che i nomi di osservazioni e variabili siano stringhe univoche."""
    a.obs_names = pd.Index(a.obs_names.astype(str))
    a.var_names = pd.Index(a.var_names.astype(str))
    if not a.obs_names.is_unique: a.obs_names_make_unique()
    if not a.var_names.is_unique: a.var_names_make_unique()
    return a

def _only_gene_expression(a: ad.AnnData) -> ad.AnnData:
    """Filtra solo le feature di tipo 'Gene Expression' se presenti."""
    if "feature_types" in a.var.columns:
        mask = a.var["feature_types"].astype(str).str.contains("Gene Expression", case=False, na=False)
        if mask.any():
            a = a[:, mask].copy()
    return a

def _dedup_columns(cols: pd.Index) -> pd.Index:
    """De-duplica i nomi delle colonne aggiungendo suffissi numerici."""
    base = cols.astype(str)
    seen = {}
    out = []
    for x in base:
        if x not in seen:
            seen[x] = 0; out.append(x)
        else:
            seen[x] += 1; out.append(f"{x}.{seen[x]}")
    return pd.Index(out)

def has_known_data(p: Path) -> bool:
    """Verifica se una cartella contiene dati riconoscibili."""
    if not p.is_dir(): return False
    # Check 10x patterns
    if list(p.glob("filtered_feature_bc_matrix*/matrix.mtx*")): return True
    if (p/"matrix.mtx.gz").exists() or (p/"matrix.mtx").exists(): return True
    # Check H5
    if list(p.glob("*.h5")): return True
    # Check Tables
    if list(p.glob("*.csv*")) or list(p.glob("*.tsv*")) or list(p.glob("*.txt*")): return True
    return False

# ------------------------------------------------------------------------------
# 2. LOADERS SPECIALIZZATI (MTX, H5, TABLE)
# ------------------------------------------------------------------------------

def _read_mtx_flexible(path: Path):
    """
    Legge MTX ignorando l'header integer se i dati sono float.
    Risolve il bug 'Invalid integer value'.
    """
    opener = gzip.open if str(path).endswith(".gz") else open
    header_lines = 0
    with opener(path, "rt", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.startswith('%'):
                header_lines += 1; continue
            parts = line.strip().split()
            if len(parts) < 3: raise ValueError("Header MTX malformato")
            n_rows, n_cols, _ = map(int, parts[:3])
            header_lines += 1
            break
    
    # Leggi con pandas engine python (lento ma flessibile per i float)
    df = pd.read_csv(
        path, sep=r"\s+", engine="python", comment="%", header=None,
        names=["i","j","x"], skiprows=header_lines
    ).dropna(subset=["i","j","x"])
    
    i = df["i"].astype(np.int64).to_numpy() - 1
    j = df["j"].astype(np.int64).to_numpy() - 1
    x = pd.to_numeric(df["x"], errors="coerce").fillna(0).astype(np.float32).to_numpy()
    return coo_matrix((x, (i, j)), shape=(n_rows, n_cols)).tocsr()

def _load_10x_mtx(p: Path) -> ad.AnnData | None:
    """
    Loader 10x MTX capace di gestire prefissi non standard e float.
    """
    # Cerca nella sottocartella standard O nella cartella corrente
    candidates = [p] + [d for d in p.glob("filtered_feature_bc_matrix*") if d.is_dir()]
    
    for base in candidates:
        matrix_file   = next(base.glob("*matrix.mtx*"),   None)
        barcodes_file = next(base.glob("*barcodes.tsv*"), None)
        features_file = next(base.glob("*features.tsv*"), None) or next(base.glob("*genes.tsv*"), None)
        
        if not (matrix_file and barcodes_file and features_file):
            continue

        # Prova lettura standard, fallback se fallisce
        try:
            M = mmread(str(matrix_file)).tocsr()
        except ValueError as e:
            if "invalid integer" in str(e).lower() or "invalid literal" in str(e).lower():
                M = _read_mtx_flexible(matrix_file)
            else:
                raise

        barcodes = pd.read_csv(barcodes_file, sep="\t", header=None, usecols=[0], dtype=str, compression="infer").iloc[:,0].astype(str)
        barcodes = _dedup_columns(pd.Index(barcodes))
        feat = pd.read_csv(features_file, sep="\t", header=None, dtype=str, compression="infer").fillna("")
        
        # Gestione colonne features.tsv (a volte mancano colonne)
        if feat.shape[1] >= 3:
            gene_id, gene_nm, ftype = feat.iloc[:,0].astype(str), feat.iloc[:,1].astype(str), feat.iloc[:,2].astype(str)
        else:
            gene_id = feat.iloc[:,0].astype(str); gene_nm = gene_id.copy()
            ftype = pd.Series(["Gene Expression"] * len(gene_id))

        var_names = np.where(gene_id.str.upper().str.startswith("ENSG"), gene_id, gene_nm).astype(str)

        # Allineamento dimensioni (safe crop se i file non combaciano)
        n_cells = min(M.shape[1], len(barcodes))
        n_genes = min(M.shape[0], len(var_names))
        if (M.shape[1] != len(barcodes)) or (M.shape[0] != len(var_names)):
            M = M[:n_genes, :n_cells].tocsr()
            barcodes = pd.Index(barcodes[:n_cells])
            var_names = var_names[:n_genes]
            ftype = ftype.iloc[:n_genes]

        a = ad.AnnData(M.T.tocsr(), obs=pd.DataFrame(index=pd.Index(barcodes, name="barcode")), var=pd.DataFrame(index=pd.Index(var_names, name="gene")))
        a.var["gene_id_raw"] = gene_id.values[:a.n_vars]
        a.var["feature_types"] = ftype.values[:a.n_vars]
        a.var_names = a.var_names.str.split(".").str[0] # Rimuovi versioni ENSG
        a = _only_gene_expression(a)
        
        sc.pp.filter_cells(a, min_counts=1)
        sc.pp.filter_genes(a, min_counts=1)
        return _make_names_unique(a)
    return None

def _load_10x_h5(p: Path) -> ad.AnnData | None:
    h5 = next(p.glob("*.h5"), None)
    if h5:
        a = sc.read_10x_h5(h5.as_posix())
        a.var_names = a.var_names.astype(str).str.split(".").str[0]
        a = _only_gene_expression(a)
        sc.pp.filter_cells(a, min_counts=1)
        sc.pp.filter_genes(a, min_counts=1)
        return _make_names_unique(a)
    return None

def _load_table(p: Path) -> ad.AnnData | None:
    files = list(p.glob("*.csv*")) + list(p.glob("*.tsv*")) + list(p.glob("*.txt*"))
    if not files: return None
    
    # Detect separator
    f = files[0]
    with open(f, "r", encoding="utf-8", errors="ignore") as fh: s = fh.read(20000)
    sep = "\t" if s.count("\t") > s.count(",") else ","

    try:
        df = pd.read_csv(f, sep=sep, engine="c", index_col=0, dtype=np.float32, low_memory=False)
    except Exception:
        df = pd.read_csv(f, sep=None, engine="python", index_col=0, on_bad_lines="skip")
        for c in df.columns: df[c] = pd.to_numeric(df[c], errors="coerce")
        df = df.fillna(0).astype(np.float32)

    # Pulizia
    df = df.loc[df.index.astype(str) != "", :]
    df = df.loc[:, df.columns.astype(str) != ""]
    if not df.index.is_unique: df = df.groupby(level=0).sum()
    if not df.columns.is_unique: df.columns = _dedup_columns(df.columns)

    # Orientamento
    if df.shape[0] >= df.shape[1]: # Righe=Geni
        X = sparse.csr_matrix(df.T.values, dtype=np.float32)
        var, obs = pd.Index(df.index.astype(str)), pd.Index(df.columns.astype(str))
    else: # Righe=Cellule
        X = sparse.csr_matrix(df.values, dtype=np.float32)
        var, obs = pd.Index(df.columns.astype(str)), pd.Index(df.index.astype(str))

    a = ad.AnnData(X, obs=obs, var=var)
    sc.pp.filter_cells(a, min_counts=1)
    sc.pp.filter_genes(a, min_counts=1)
    return _make_names_unique(a)

def load_one_sample(sample_dir: Path) -> ad.AnnData:
    a = _load_10x_mtx(sample_dir)
    if a is None: a = _load_10x_h5(sample_dir)
    if a is None: a = _load_table(sample_dir)
    if a is None: raise FileNotFoundError(f"Nessun formato valido in {sample_dir}")
    if not sparse.issparse(a.X): a.X = sparse.csr_matrix(a.X)
    return _make_names_unique(a)

# ------------------------------------------------------------------------------
# 3. MAPPING & QC
# ------------------------------------------------------------------------------

def to_ensembl_and_collapse(a: ad.AnnData, gene_map: dict, symbol_to_ens: dict) -> ad.AnnData:
    vn = pd.Index(a.var_names.astype(str))
    is_ens = vn.str.upper().str.startswith("ENSG")
    
    # Mapping sicuro senza in-place
    mapped_vals = pd.Series(vn, index=vn).map(symbol_to_ens)
    new_names = np.where((~is_ens) & pd.notna(mapped_vals), mapped_vals, vn)
    a.var_names = pd.Index(new_names)

    # Collasso duplicati (somma sparsa)
    if not a.var_names.is_unique:
        codes, uniques = pd.factorize(a.var_names, sort=False)
        G = sparse.csr_matrix((np.ones_like(codes), (codes, np.arange(len(codes)))), shape=(len(uniques), len(codes)))
        a = ad.AnnData((a.X @ G.T).tocsr(), obs=a.obs.copy())
        a.var_names = pd.Index(uniques.astype(str))

    # Ricostruzione .var
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
    if not sparse.issparse(a.X): a.X = sparse.csr_matrix(a.X)
    a.layers["counts"] = a.X.astype(np.int32)
    return _make_names_unique(a)

# ------------------------------------------------------------------------------
# 4. ESECUZIONE (CONFIGURAZIONE E LOOP)
# ------------------------------------------------------------------------------

# Configurazione Percorsi
RAW_DIR   = Path("../data/raw")
OUT_DIR   = Path("../work/01_qc_h5ad"); OUT_DIR.mkdir(parents=True, exist_ok=True)
GTF_PATH  = Path(r"C:\Users\Utente\Desktop\Tesi_Python_scRNA\data\gtf\gencode.v43.annotation.gtf.gz")

# Parametri QC
MT_PCT    = 20.0
MIN_GENES = 200
MIN_CELLS = 3
FORCE_RERUN = False # Metti True per sovrascrivere

# Caricamento GTF
if not GTF_PATH.exists(): raise FileNotFoundError(f"GTF non trovato: {GTF_PATH}")
gene_map = load_gtf_gene_map(GTF_PATH)
symbol_to_ens = {v: k for k, v in gene_map.items()}
print(f"[Setup] GTF caricato: {len(gene_map):,} geni.")

# Logica di scoperta cartelle (esclude temp e spazzatura)
def discover_folders() -> list[Path]:
    candidates = []
    EXCLUDE_NAMES = {"GSE111976", "GSE111976_processed_temp", "C1", "_ignore"}
    
    # 1. Scansiona cartelle standard
    for p in RAW_DIR.iterdir():
        if not p.is_dir(): continue
        if p.name in EXCLUDE_NAMES or p.name.startswith("_"): continue
        if has_known_data(p): candidates.append(p)
    
    # 2. Aggiungi split di GSE111976 (se presenti)
    for p in RAW_DIR.glob("eutopic_*_GSM*"):
        if p.is_dir() and has_known_data(p): candidates.append(p)

    # 3. Deduplica
    seen = set(); unique = []
    for p in sorted(candidates, key=lambda x: x.name):
        if p.name not in seen:
            seen.add(p.name); unique.append(p)
    return unique

# Esecuzione Loop
folders = discover_folders()
print(f"[INFO] Campioni da processare: {len(folders)}")

for sample_dir in tqdm(folders):
    sid = sample_dir.name
    out = OUT_DIR / f"{sid}.h5ad"
    
    if out.exists() and not FORCE_RERUN:
        continue

    try:
        print(f"\n[PROC] {sid}")
        # 1. Caricamento (Smart Loader)
        a = load_one_sample(sample_dir)
        
        # 2. Standardizzazione (Ensembl + Collasso)
        a = to_ensembl_and_collapse(a, gene_map, symbol_to_ens)
        
        # 3. QC e Metadati
        a = apply_qc(a, mt_pct=MT_PCT, min_genes=MIN_GENES, min_cells=MIN_CELLS)
        a.obs["sample_id"] = sid
        
        # 4. Salvataggio
        a.write_h5ad(out.as_posix(), compression="gzip")
        print(f"  -> OK: {out.name}  shape={a.shape}")
        
    except Exception as e:
        import traceback
        print(f"  -> FAILED: {e}")
        # traceback.print_exc(limit=1) # Scommenta per debug dettagliato

print("\n=== FASE 1 COMPLETATA ===")