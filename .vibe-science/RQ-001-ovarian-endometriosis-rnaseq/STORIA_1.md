Relazione Tecnica: Migrazione Workflow scRNA-seq (R Seurat → Python Scanpy/scVI)1. Contesto e MotivazioniIl Progetto: Analisi di un atlante single-cell RNA-seq comprendente 52 campioni (Endometrio Eutopico, Endometriosi/OMA, Ovaio Sano, e vari istotipi di Carcinoma Ovarico: CCOC, EnOC, LGsOC) per un totale di oltre 390.000 cellule.Il Problema (R/Seurat):Il workflow originale in R (Seurat, Monocle3) si è interrotto durante le fasi di integrazione e riduzione dimensionale (UMAP). Nonostante l'utilizzo di una istanza AWS con 750 GB di RAM, R saturava la memoria.Causa: R carica interamente le matrici dense in memoria. La costruzione del grafo dei vicini (Nearest Neighbor) e l'integrazione con IntegrateData su >300k cellule richiedono risorse esponenziali.La Soluzione (Python/Scanpy + scVI):Abbiamo migrato l'analisi in Python per sfruttare:Sparsità: Gestione nativa di matrici sparse (CSR) tramite AnnData, riducendo l'uso della RAM del 90%.Scalabilità (scVI): Utilizzo di scvi-tools (Variational AutoEncoder) per l'integrazione batch-aware, che utilizza mini-batching e non richiede il caricamento totale dei dati in memoria.UMAP Ottimizzata: Strategia "Subset + Transform" (calcolo UMAP su un sottoinsieme e proiezione del resto) per evitare il calcolo di grafi giganti.2. Fase 1: Ingestion e Pre-processing (Completata)Stato: I 52 campioni grezzi (formati misti: .mtx, .h5, .csv) sono stati processati individualmente.Operazioni Svolte:Caricamento dati grezzi.Standardizzazione nomi geni (Mappatura Ensembl ID univoci).QC preliminare per campione (filtro cellule morte/vuote).Salvataggio di 52 file .h5ad nella cartella work/01_qc_h5ad.3. Fase 2: Merge, Annotazione e Pulizia Strutturale (Codice Definitivo)Questa è la fase critica appena conclusa. L'obiettivo era unire 52 file in un unico oggetto AnnData senza esplodere la RAM e correggendo errori strutturali.Problemi Riscontrati e Bug FixatiBug "GTF Path": Su Windows/Anaconda, il percorso relativo del file GTF (necessario per rimappare i geni) causava FileNotFoundError.Fix: Utilizzo di pathlib e percorsi assoluti verificati.Bug "Collisione Geni Duplicati": Durante il merge, alcuni geni apparivano duplicati con suffissi (es. GeneA, GeneA-1).Fix: Implementata una logica post-merge che somma le colonne dei geni duplicati e ricostruisce la matrice sparsa.Bug "Mapping Gruppi Biologici": La logica precedente classificava erroneamente i carcinomi (es. "carcinOMA") come "OMA" (Endometrioma) a causa di una regex troppo permissiva. Inoltre c'era un typo "CCCOC".Fix: Riscritta la funzione map_group con regex rigorose e priorità ai carcinomi.Codice Definitivo (Cella 13)Questo script unisce i file, pulisce i duplicati, applica le annotazioni corrette e salva in formato .h5ad e .zarr (ottimizzato per letture parziali).# Cella 13: Fase 2 — Merge, Pulizia Strutturale e Annotazione (Versione Finale Windows/Anaconda)

import os, glob, re, gc, gzip
from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
from scipy import sparse

print("--- FASE 2: INIZIO MERGE E PULIZIA STRUTTURALE ---")

# ---------------------------------------------------------
# 1. Caricamento Mappa Geni (GTF) - Fix Path Assoluto
# ---------------------------------------------------------
def load_gtf_gene_map(gtf_path: str) -> dict:
    mp = {}
    # Usa gzip per leggere il file .gz compresso
    with gzip.open(gtf_path, "rt", encoding="utf-8") as f:
        for line in f:
            if line.startswith("#"): continue
            cols = line.rstrip("\n").split("\t")
            # Parsing del GTF per estrarre gene_id (Ensembl) e gene_name (Symbol)
            if len(cols) < 9 or cols[2] != "gene": continue
            gid = re.search(r'gene_id "([^"]+)"', cols[8])
            gnm = re.search(r'gene_name "([^"]+)"', cols[8])
            if gid and gnm:
                ens = gid.group(1).split(".")[0] # Rimuove versione (es. ENSG000.1 -> ENSG000)
                if ens not in mp:
                    mp[ens] = gnm.group(1)
    return mp

# Impostazione Path Assoluto per Windows
GTF_PATH = Path(r"C:\Users\Utente\Desktop\Tesi_Python_scRNA\data\gtf\gencode.v43.annotation.gtf.gz")
assert GTF_PATH.exists(), f"File non trovato: {GTF_PATH}"
gene_map = load_gtf_gene_map(GTF_PATH.as_posix())
print(f"Gene map caricata: {len(gene_map):,} geni.")

# ---------------------------------------------------------
# 2. Caricamento e Concatenazione Campioni
# ---------------------------------------------------------
INPUT_DIR = "../work/01_qc_h5ad"
files = sorted(glob.glob(os.path.join(INPUT_DIR, "*.h5ad")))
print(f"Trovati {len(files)} campioni. Caricamento...")

adatas = {}
for f in files:
    sid = os.path.basename(f)[:-5]  # rimuove ".h5ad" per ottenere il Sample ID
    a = sc.read_h5ad(f)
    # Gestione preliminare duplicati pre-merge
    if not a.var_names.is_unique:
        a.var_names_make_unique(join='-')
    adatas[sid] = a  # Il dizionario preserva l'ID per il merge

print("Concatenazione in corso...")
# Merge outer: mantiene tutti i geni presenti in almeno un campione
adata = ad.concat(adatas, join="outer", label="sample_id",
                  fill_value=0, index_unique="-")
del adatas; gc.collect() # Pulizia RAM immediata
print(f"Merge completato. Forma iniziale: {adata.shape}")

# ---------------------------------------------------------
# 3. Fix Strutturale: Collasso Geni Duplicati
# ---------------------------------------------------------
# Se il merge ha creato duplicati (es. GENE-1), li sommiamo
if not adata.var_names.is_unique:
    print("Geni duplicati post-merge trovati. Collasso per somma...")
    # Rimuove i suffissi numerici per trovare la radice del gene
    base_ids = adata.var_names.to_series().str.replace(r"-\d+$", "", regex=True)
    groups = base_ids.groupby(base_ids).indices  # Mappa: Gene -> [colonna1, colonna2...]

    new_X = sparse.lil_matrix((adata.n_obs, len(groups)), dtype=adata.X.dtype)
    new_var_names = []
    
    # Itera sui gruppi di geni duplicati e somma le colonne
    for j, (gid, idxs) in enumerate(groups.items()):
        new_var_names.append(gid)
        s = adata.X[:, list(idxs)].sum(axis=1)
        if sparse.issparse(s):
            new_X[:, j] = s
        else:
            new_X[:, j] = sparse.csr_matrix(s).T

    # Ricrea l'oggetto AnnData pulito
    adata = ad.AnnData(new_X.tocsr(), obs=adata.obs.copy())
    adata.var_names = pd.Index(new_var_names)
    print(f"Geni duplicati collassati. Nuova forma: {adata.shape}")

# ---------------------------------------------------------
# 4. Ricostruzione Annotazioni Geniche (.var)
# ---------------------------------------------------------
adata.var["ensembl_id"]  = adata.var_names
adata.var["gene_symbol"] = adata.var_names.map(gene_map).astype("object")
# Annotazione geni mitocondriali
sym = adata.var["gene_symbol"].fillna("").astype(str).str.upper()
adata.var["mt"] = sym.str.startswith("MT-")
adata.var_names_make_unique()

# ---------------------------------------------------------
# 5. Annotazione Gruppi Biologici (Mapping Robusto)
# ---------------------------------------------------------
def map_group(sample_id: str) -> str:
    sid = sample_id.upper()
    # Priorità ai carcinomi per evitare falsi match con 'OMA'
    if "LGSOC" in sid or "LG-SOC" in sid: return "LGsOC"
    if "CCOC"  in sid or "CLEARCELL" in sid: return "CCOC"
    if "ENOC_S" in sid or "ENOC-S" in sid or "OVCA" in sid: return "EnOC_S"
    if re.search(r"\bENOC\b", sid): return "EnOC"
    
    # Endometrioma / OMA (Regex boundary \b per evitare 'carcinOMA')
    if "ENDOMETRIOMA" in sid or re.search(r"\bOMA\b", sid): return "OMA"
    
    # Tessuti sani / Eutopici / Endometriosi generica
    if "EUTOPIC" in sid or "EUE" in sid: return "EuE"
    if "NORMAL" in sid or "NOV" in sid: return "nOV"
    if "ENDOMETRIOSIS" in sid or "ENDOMETRIOSI" in sid or "ENDO_" in sid: return "Endo"
    return "Unknown"

# Casting a category per risparmiare memoria
adata.obs["sample_id"] = adata.obs["sample_id"].astype("category")
adata.obs["princ_etic_1"] = adata.obs["sample_id"].map(map_group).astype("category")
print("Distribuzione Gruppi princ_etic_1:\n", adata.obs["princ_etic_1"].value_counts())

# ---------------------------------------------------------
# 6. Pulizia Finale e Salvataggio
# ---------------------------------------------------------
# Filtro geni rari (meno di 20 cellule) per alleggerire scVI
print(f"Geni prima del filtro: {adata.n_vars}")
sc.pp.filter_genes(adata, min_cells=20)
print(f"Geni dopo filtro: {adata.n_vars}")

# Assicura formato sparso CSR e Interi per i conteggi raw
if not sparse.issparse(adata.X):
    adata.X = sparse.csr_matrix(adata.X)
adata.layers["counts"] = adata.X.astype(np.int32)

# Sanity Check finale
mapped_ratio = adata.var["gene_symbol"].notna().mean()
print(f"Gene symbols mappati correttamente: {mapped_ratio:.2%}")

# Export CSV per controllo manuale
adata.obs.groupby("sample_id").size().to_csv("../work/per_sample_counts.csv")
adata.obs.groupby("princ_etic_1").size().to_csv("../work/group_counts.csv")

# Salvataggio
adata.write_h5ad("../work/02_merged_raw.h5ad")
# Zarr è consigliato per letture parziali future
adata.write_zarr("../work/02_merged_raw.zarr", chunk_size=(10_000, 2_000))
print("\n--- FASE 2 COMPLETATA CON SUCCESSO ---")
4. Fase 2.5: Rimozione Doppietti (Next Step)Obiettivo: Utilizzare scvi-tools (modello SOLO) per identificare le droplet contenenti due cellule invece di una. Eseguiamo questo prima dell'integrazione finale per pulire il dataset.Codice (Cella 13.5)# Cella 13.5: Rilevamento Doppietti (SOLO)
import scvi

print("\n--- FASE 2.5: RILEVAMENTO DOPPIETTI ---")
# Caricamento del dataset unito
adata = sc.read_h5ad("../work/02_merged_raw.h5ad")

# Setup scVI: usa i conteggi raw e il sample_id come batch
scvi.model.SCVI.setup_anndata(adata, layer="counts", batch_key="sample_id")

# 1. Training modello base (VAE)
print("Training modello base SCVI...")
base_model = scvi.model.SCVI(adata, n_latent=20)
base_model.train(max_epochs=70, early_stopping=True)

# 2. Training modello SOLO (classificatore doppietti)
print("Training modello SOLO...")
solo_model = scvi.external.SOLO.from_scvi_model(base_model)
solo_model.train(max_epochs=50, early_stopping=True)

# 3. Predizione e Filtro
pred = solo_model.predict(as_dataframe=True)
adata.obs["doublet_score"] = pred["doublet_score"].values
adata.obs["doublet_prediction"] = pred["prediction"].values

print("Risultati SOLO:\n", adata.obs["doublet_prediction"].value_counts())

# Manteniamo solo i singlets
adata_singlets = adata[adata.obs["doublet_prediction"] == "singlet"].copy()
print(f"Rimossi {adata.n_obs - adata_singlets.n_obs} doppietti.")

adata_singlets.write_h5ad("../work/02b_merged_nodoublets.h5ad")
print("File pulito salvato: 02b_merged_nodoublets.h5ad")
5. Fase 3: Integrazione, UMAP e Clustering (Next Step)Obiettivo: Creare l'atlante finale.Strategia Anti-Crash:HVG Batch-Aware: Selezionare solo i 4000 geni più informativi considerando la varianza tra campioni.Integrazione scVI: Correzione batch effect senza densificazione.UMAP Subset+Transform: Invece di calcolare il grafo su 300k cellule (che satura la RAM), calcoliamo l'UMAP su 80k cellule e "proiettiamo" matematicamente le restanti.Codice (Cella 14)# Cella 14: Fase 3 — Integrazione, Clustering e Analisi

import scvi
import numpy as np
import umap

print("\n--- FASE 3: INTEGRAZIONE E CLUSTERING ---")
scvi.settings.seed = 42
adata = sc.read_h5ad("../work/02b_merged_nodoublets.h5ad")

# 1. Selezione HVG (riduce drasticamente la memoria necessaria)
print("Selezione HVG...")
sc.pp.highly_variable_genes(
    adata, n_top_genes=4000, flavor="seurat_v3",
    batch_key="sample_id", subset=True, layer="counts"
)

# 2. Integrazione (scVI)
print("Training scVI per integrazione batch...")
scvi.model.SCVI.setup_anndata(adata, layer="counts", batch_key="sample_id")
model = scvi.model.SCVI(adata, n_latent=30, n_layers=2, gene_likelihood="nb")
model.train(max_epochs=100, early_stopping=True, early_stopping_patience=10)
# Salviamo l'embedding latente (spazio corretto dai batch)
adata.obsm["X_scVI"] = model.get_latent_representation()

# 3. UMAP (Tecnica Fail-Safe: Subset + Transform)
print("Calcolo UMAP (Subset 80k)...")
X = adata.obsm["X_scVI"]
rng = np.random.default_rng(42)
idx = rng.permutation(adata.n_obs)
n_sub = min(80000, adata.n_obs) # Usa subset per evitare crash RAM
sub, rest = idx[:n_sub], idx[n_sub:]

# Fit UMAP solo sul subset
um = umap.UMAP(n_neighbors=15, min_dist=0.4, metric="cosine", random_state=42)
emb = np.zeros((adata.n_obs, 2), np.float32)
emb[sub] = um.fit_transform(X[sub])
# Proiezione del resto
if rest.size > 0:
    emb[rest] = um.transform(X[rest])
adata.obsm["X_umap"] = emb

# 4. Clustering Leiden
print("Clustering Leiden...")
sc.pp.neighbors(adata, use_rep="X_scVI", n_neighbors=15)
sc.tl.leiden(adata, resolution=0.8, key_added="leiden")

# 5. Marcatori (Differential Expression) con Downsampling
print("Calcolo Marcatori...")
# Usiamo max 500 cellule per cluster per velocizzare il test statistico
keep = []
for g in adata.obs["leiden"].cat.categories:
    idx = adata.obs_names[adata.obs["leiden"]==g]
    keep.extend(rng.choice(idx, size=min(500, len(idx)), replace=False))

adata_ds = adata[keep].copy()
# Normalizzazione al volo solo per il calcolo statistico
sc.pp.normalize_total(adata_ds, target_sum=1e4)
sc.pp.log1p(adata_ds)
sc.tl.rank_genes_groups(adata_ds, "leiden", method="wilcoxon")

# Salvataggio Risultati
sc.get.rank_genes_groups_df(adata_ds, None).to_csv("../work/markers_leiden.csv", index=False)
adata.write_h5ad("../work/03_integrated.h5ad")
model.save("../work/scvi_model", overwrite=True)

print("--- ANALISI COMPLETATA: Atlante integrato salvato in 03_integrated.h5ad ---")


