# SCOPUS QUERY PLAN — UPR in EAOC

Gap claim da blindare: **"Nessuno studio transcriptomico sistematico ha caratterizzato l'attivazione dell'UPR nella progressione endometriosi → cancro ovarico"**

---

## PACCHETTO A — Core Claim (UPR/ER stress + EAOC)

```
A1  TITLE-ABS-KEY(("unfolded protein response" OR "UPR") AND (endometriosis OR endometrioma) AND ("ovarian cancer" OR "ovarian carcinoma" OR "ovarian neoplasm*"))

A2  TITLE-ABS-KEY(("ER stress" OR "endoplasmic reticulum stress") AND (endometriosis OR endometrioma) AND ("ovarian cancer" OR "ovarian carcinoma"))

A3  TITLE-ABS-KEY((GRP78 OR BiP OR HSPA5) AND (endometriosis) AND ("ovarian cancer" OR "ovarian carcinoma" OR EAOC))

A4  TITLE-ABS-KEY((XBP1 OR "X-box binding protein") AND (endometriosis) AND ("ovarian cancer" OR "ovarian carcinoma"))

A5  TITLE-ABS-KEY((CHOP OR DDIT3 OR "C/EBP homologous protein") AND (endometriosis) AND ("ovarian cancer" OR "ovarian carcinoma"))

A6  TITLE-ABS-KEY((ATF4 OR ATF6 OR "activating transcription factor") AND (endometriosis) AND ("ovarian cancer" OR "ovarian carcinoma"))

A7  TITLE-ABS-KEY((IRE1 OR ERN1 OR PERK OR EIF2AK3) AND (endometriosis) AND ("ovarian cancer" OR "ovarian carcinoma"))

A8  TITLE-ABS-KEY(("ER stress" OR "UPR" OR "unfolded protein response") AND ("clear cell" W/2 (carcinoma OR cancer OR ovarian)) AND (endometriosis OR endometrioma))

A9  TITLE-ABS-KEY(("ER stress" OR "UPR") AND ("endometrioid" W/2 (carcinoma OR cancer OR ovarian)) AND (endometriosis))

A10 TITLE-ABS-KEY(("endometriosis-associated ovarian cancer" OR "EAOC" OR "endometriosis associated ovarian") AND ("ER stress" OR "UPR" OR "unfolded protein" OR GRP78 OR XBP1 OR CHOP))

A11 TITLE-ABS-KEY((endometriosis) AND ("malignant transformation" OR "carcinogenesis" OR "oncogenic transformation") AND ("ER stress" OR "UPR" OR "proteostasis" OR "protein homeostasis"))

A12 TITLE-ABS-KEY(("ovarian clear cell" OR OCCC) AND ("ER stress" OR "endoplasmic reticulum stress" OR GRP78 OR BiP))
```

**Perché questo pacchetto:** Ricerca diretta del gap claim con tutte le varianti terminologiche UPR + istotipi EAOC.

---

## PACCHETTO B — Competitor Scan (studi transcriptomici esistenti)

```
B1  TITLE-ABS-KEY((endometriosis) AND ("ovarian cancer" OR "ovarian carcinoma") AND ("RNA-seq" OR "RNA sequencing" OR "transcriptom*"))

B2  TITLE-ABS-KEY((endometriosis) AND ("ovarian cancer" OR "ovarian carcinoma") AND ("single-cell" OR "scRNA-seq" OR "single cell RNA"))

B3  TITLE-ABS-KEY((endometriosis) AND ("ovarian cancer") AND ("gene expression profiling" OR "expression profiling" OR microarray))

B4  TITLE-ABS-KEY(("EAOC" OR "endometriosis-associated ovarian cancer") AND (transcriptom* OR "RNA-seq" OR "gene expression"))

B5  TITLE-ABS-KEY(("clear cell ovarian" OR OCCC) AND (endometriosis) AND (transcriptom* OR "RNA-seq" OR "gene expression"))

B6  TITLE-ABS-KEY((endometriosis) AND ("malignant transformation") AND (transcriptom* OR "expression profil*" OR "pathway analysis"))

B7  TITLE-ABS-KEY((endometriosis) AND ("ovarian cancer") AND ("differential expression" OR "DEG" OR "differentially expressed gene*"))

B8  TITLE-ABS-KEY((endometriosis) AND ("ovarian cancer") AND ("pathway enrichment" OR "GSEA" OR "gene set enrichment" OR "pathway analysis"))
```

**Perché questo pacchetto:** Identificare TUTTI gli studi transcriptomici su EAOC per verificare se qualcuno ha già analizzato UPR.

---

## PACCHETTO C — Oxidative Stress Link (upstream del UPR)

```
C1  TITLE-ABS-KEY(("oxidative stress") AND (endometriosis) AND ("ovarian cancer" OR "ovarian carcinoma" OR "malignant transformation"))

C2  TITLE-ABS-KEY((ROS OR "reactive oxygen species") AND (endometriosis) AND ("ovarian cancer" OR "malignant transformation"))

C3  TITLE-ABS-KEY(("oxidative stress") AND (endometriosis) AND ("ER stress" OR "UPR" OR "endoplasmic reticulum stress"))

C4  TITLE-ABS-KEY(("oxidative stress") AND ("ovarian cancer") AND ("ER stress" OR "UPR" OR "unfolded protein response"))

C5  TITLE-ABS-KEY((endometriosis) AND ("ovarian cancer") AND ("iron" OR "hemoglobin" OR "heme" OR "ferroptosis") AND ("oxidative stress" OR "ROS"))

C6  TITLE-ABS-KEY(("oxidative stress") AND ("clear cell ovarian" OR OCCC) AND (endometriosis OR "endometriosis-associated"))
```

**Perché questo pacchetto:** Oxidative stress è upstream di ER stress. Se 139 paper esistono su OS ma quasi 0 su UPR, il gap è nel "downstream pathway".

---

## PACCHETTO D — ARID1A Link (mutazione chiave in EAOC)

```
D1  TITLE-ABS-KEY((ARID1A) AND (endometriosis) AND ("ovarian cancer" OR "ovarian carcinoma"))

D2  TITLE-ABS-KEY((ARID1A) AND ("ER stress" OR "UPR" OR "unfolded protein response"))

D3  TITLE-ABS-KEY((ARID1A) AND ("clear cell ovarian" OR "endometrioid ovarian"))

D4  TITLE-ABS-KEY((ARID1A) AND (endometriosis) AND ("malignant transformation" OR "carcinogenesis"))

D5  TITLE-ABS-KEY(("SWI/SNF" OR "BAF complex") AND (endometriosis) AND ("ovarian cancer"))
```

**Perché questo pacchetto:** ARID1A è la mutazione più frequente in EAOC. Se c'è un link ARID1A-UPR, è un meccanismo da esplorare.

---

## PACCHETTO H — Anti-False-Negative (proximity + wildcards)

```
H1  TITLE-ABS-KEY((endometri* W/3 (cancer* OR carcinoma* OR malignan*)) AND ("ER stress" OR "UPR" OR "protein folding stress"))

H2  TITLE-ABS-KEY(("ovarian" PRE/2 (cancer* OR carcinoma*)) AND (endometri*) AND (GRP78 OR BiP OR HSPA5 OR "heat shock protein"))

H3  TITLE-ABS-KEY((endometri* W/5 ovar*) AND ("unfolded protein*" OR "misfolded protein*" OR "protein aggregat*"))

H4  TITLE-ABS-KEY(((clear W/2 cell) OR OCCC) AND (endometri*) AND ("stress response" OR "cellular stress" OR "proteotoxic*"))

H5  TITLE-ABS-KEY((endometri*) AND (ovar* W/3 (cancer* OR neoplas* OR tumor* OR tumour*)) AND ("ER stress" OR "endoplasmic reticulum" W/2 stress))

H6  TITLE-ABS-KEY((EAOC OR "endometriosis associated" OR "endometriosis-associated") AND (stress W/3 (protein* OR endoplasmic OR ER OR cellular)))

H7  TITLE-ABS-KEY((endometri*) AND ("malignant transform*" OR "oncogenic transform*" OR carcino*) AND (chaperon* OR GRP* OR HSP* OR "heat shock"))

H8  TITLE-ABS-KEY((endometri*) AND (ovar*) AND ("proteostasis" OR "protein homeostasis" OR "protein quality control"))
```

**Perché questo pacchetto:** Query con proximity operators per catturare varianti terminologiche e non perdere paper.

---

## PACCHETTO I — Tool/Method Hunting

```
I1  TITLE-ABS-KEY((endometriosis) AND ("ovarian cancer") AND ("UPR" OR "ER stress") AND (tool OR software OR pipeline OR "R package" OR method*))

I2  TITLE-ABS-KEY((endometriosis) AND ("ovarian cancer") AND ("gene signature" OR "expression signature" OR "prognostic signature") AND ("ER stress" OR "UPR"))

I3  TITLE-ABS-KEY(("UPR score" OR "ER stress score" OR "UPR signature") AND ("ovarian cancer" OR "gynecologic* cancer*"))

I4  TITLE-ABS-KEY(("UPR gene*" OR "ER stress gene*") AND (biomarker* OR prognostic*) AND ("ovarian" OR "gynecologic*"))

I5  TITLE-ABS-KEY(("endoplasmic reticulum stress" OR "UPR") AND (signature) AND (cancer OR carcinoma) AND (ovary OR ovarian OR endometri*))
```

**Perché questo pacchetto:** Verificare se esistono signature UPR già applicate a cancro ovarico.

---

## PACCHETTO J — Terminology Traps (sinonimi e varianti)

```
J1  TITLE-ABS-KEY((endometriosis) AND ("ovarian cancer") AND ("integrated stress response" OR "ISR"))

J2  TITLE-ABS-KEY((endometriosis) AND ("ovarian cancer") AND ("autophagy" W/3 ("ER stress" OR "endoplasmic")))

J3  TITLE-ABS-KEY((endometriosis) AND ("ovarian cancer") AND ("apoptosis" OR "cell death") AND (GRP78 OR CHOP OR ATF4))

J4  TITLE-ABS-KEY((endometriosis) AND ("ovarian cancer") AND ("protein folding" OR "protein misfolding" OR "aggregation"))

J5  TITLE-ABS-KEY((endometriosis) AND ("ovarian cancer") AND ("secretory pathway" OR "ER-Golgi" OR "secretion stress"))

J6  TITLE-ABS-KEY((endometriosis) AND (ovar*) AND ("ERAD" OR "ER-associated degradation" OR "retrotranslocation"))

J7  TITLE-ABS-KEY((endometriosis) AND ("ovarian cancer") AND ("calcium homeostasis" OR "calcium signaling") AND ("ER" OR "endoplasmic"))
```

**Perché questo pacchetto:** UPR è collegato a ISR, autophagy, ERAD, calcium signaling. Catturare paper che usano terminologia alternativa.

---

## PACCHETTO K — Mechanism/Pathway

```
K1  TITLE-ABS-KEY((endometriosis) AND ("malignant transformation") AND (pathway* OR mechanism* OR "molecular basis"))

K2  TITLE-ABS-KEY((endometriosis) AND ("ovarian cancer") AND ("cell survival" OR "cell proliferation") AND (stress OR pathway*))

K3  TITLE-ABS-KEY(("clear cell ovarian") AND (mechanism* OR pathway*) AND (stress OR "metabolic reprogram*"))

K4  TITLE-ABS-KEY((endometriosis) AND ("ovarian cancer") AND ("NF-kB" OR "inflammatory pathway") AND stress)

K5  TITLE-ABS-KEY((EAOC OR "endometriosis-associated ovarian") AND ("molecular mechanism*" OR "molecular pathway*" OR "signaling pathway*"))

K6  TITLE-ABS-KEY((endometriosis) AND ("ovarian cancer") AND (PI3K OR AKT OR mTOR) AND ("ER stress" OR "UPR" OR "protein synthesis"))
```

**Perché questo pacchetto:** Pathway studies dove UPR potrebbe essere menzionato come componente.

---

## PACCHETTO M — Exclusion Queries (pulizia)

```
M1  TITLE-ABS-KEY((endometriosis) AND ("ovarian cancer") AND ("ER stress" OR "UPR" OR GRP78 OR XBP1 OR CHOP)) AND NOT (review)

M2  TITLE-ABS-KEY((endometriosis) AND ("ovarian cancer") AND ("unfolded protein response")) AND (LIMIT-TO (DOCTYPE, "ar"))
```

**Perché questo pacchetto:** Separare articoli originali dalle review.

---

## EXECUTION PLAN

1. Eseguire tutte le query A1-A12 (core claim)
2. Per ogni query: registrare N risultati
3. Se N > 0: scaricare titoli, anni, journal, DOI
4. Taggare ogni paper: **UPR-focused** / **UPR-mentioned** / **UPR-absent**
5. Eseguire B1-B8 (competitor transcriptomic)
6. Cross-reference: quali studi transcriptomici hanno analizzato UPR?
7. Eseguire H1-H8 (anti-false-negative)
8. Merge e deduplica
9. Output finale: lista paper con tag

---

## OUTPUT TEMPLATE (per ogni query)

| ID | Query | N | Top papers (anno, titolo, journal) | Tag |
|----|-------|---|-----------------------------------|-----|
| A1 | ... | X | ... | UPR-focused/mentioned/absent |
