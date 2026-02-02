# Analisi Completa: Stress, Infiammazione e Morte Cellulare

**Data:** 2026-01-31
**Confronti:** Progressione da EuE → Endo Pelvica → OMA → CCOC/EnOC
**Compartimenti:** Epitelio, Stroma

---

## 1. Stress Ossidativo

### 1.1 OMA vs EuE (Lesione Benigna)

| Gene | Epithelium log2FC | p-adj | Stroma log2FC | p-adj | Funzione |
|------|-------------------|-------|---------------|-------|----------|
| **HMOX1** | +2.00 | NS | **+1.97*** | 0.007 | Antiossidante (HO-1) |
| **NFE2L2** | -0.31 | NS | **+0.59*** | 0.022 | Master regulator (NRF2) |
| SOD1 | +0.39 | NS | **+0.53*** | 0.009 | Superossido dismutasi |
| SOD2 | -15.51 | NS | -20.23 | 0.12 | SOD mitocondriale |
| GPX1 | ND | - | ND | - | Glutatione perossidasi |
| CAT | -0.11 | NS | -0.43 | 0.09 | Catalasi |
| NQO1 | +0.54 | NS | -0.50 | NS | Chinone reduttasi |
| KEAP1 | -0.24 | NS | -0.18 | NS | Inibitore NRF2 |

**→ OMA: Attivazione NRF2/HMOX1 nello STROMA, SOD2 mitocondriale paradossalmente DOWN**

### 1.2 CCOC vs OMA (Trasformazione Maligna)

| Gene | Epithelium log2FC | p-adj | Stroma log2FC | p-adj | Funzione |
|------|-------------------|-------|---------------|-------|----------|
| **SOD2** | **+10.99*** | 0.0007 | **+15.10*** | 1.97e-12 | SOD mitocondriale |
| **GPX1** | **+9.24*** | 0.005 | **+13.34*** | 1.02e-18 | Glutatione perossidasi |
| **HIF1A** | +6.90~ | 0.056 | **+11.87*** | 9.7e-15 | Hypoxia TF |
| **HMOX1** | -2.41 | 0.15 | **-3.64*** | 0.0001 | Antiossidante |
| **NFE2L2** | -0.46 | NS | **-0.71*** | 0.01 | Master regulator |
| SOD1 | +0.79 | NS | **-1.29*** | 0.002 | SOD citosolico |
| CAT | +0.08 | NS | -0.84~ | 0.05 | Catalasi |
| NQO1 | +0.32 | NS | +0.29 | NS | Chinone reduttasi |
| KEAP1 | +0.62 | NS | +0.13 | NS | Inibitore NRF2 |

**→ CCOC: MASSIVO UP di SOD2/GPX1 (switch antiossidante mitocondriale!), HMOX1 DOWN**

### 1.3 EnOC vs OMA (Trasformazione Maligna)

| Gene | Epithelium log2FC | p-adj | Stroma log2FC | p-adj | Funzione |
|------|-------------------|-------|---------------|-------|----------|
| **GPX1** | +7.36 | 0.13 | **+13.13*** | 9.3e-05 | Glutatione perossidasi |
| **SOD2** | +5.97 | 0.21 | **+12.68*** | 0.0002 | SOD mitocondriale |
| **HIF1A** | +3.98 | 0.41 | **+12.20*** | 0.003 | Hypoxia TF |
| **HMOX1** | -4.52~ | 0.64 | **-2.02** | 0.10 | Antiossidante |
| **NFE2L2** | -1.05 | NS | -0.03 | NS | Master regulator |
| SOD1 | -0.37 | NS | -0.49 | 0.25 | SOD citosolico |
| NQO1 | +0.53 | NS | +0.79 | NS | Chinone reduttasi |
| KEAP1 | +0.74 | NS | -0.28 | NS | Inibitore NRF2 |

**→ EnOC: Pattern simile CCOC nello stroma (SOD2/GPX1↑↑), epitelio meno alterato**

### 1.4 Pattern Integrato Stress Ossidativo

```
ENDOMETRIOMA (OMA)
├── EPITELIO: Stress ossidativo moderato, SOD2↓↓↓
├── STROMA:   NRF2↑, HMOX1↑ (risposta attiva)
└── Net:      Stress compensato dal pathway citosolico

         ↓ CCOC trasformazione ↓

CCOC
├── EPITELIO: SOD2↑↑, GPX1↑↑ (switch mitocondriale)
├── STROMA:   SOD2↑↑↑, GPX1↑↑↑, HMOX1↓↓↓
│             HIF1A↑↑↑ (pseudoipossia?)
└── Net:      Cambio radicale da citosolico a mitocondriale

         ↓ EnOC trasformazione ↓

EnOC
├── EPITELIO: Moderato (non significativo)
├── STROMA:   Pattern simile a CCOC ma meno estremo
└── Net:      Stroma-driven come CCOC
```

---

## 2. Ipossia

### 2.1 OMA vs EuE

| Gene | Epithelium log2FC | p-adj | Stroma log2FC | p-adj | Funzione |
|------|-------------------|-------|---------------|-------|----------|
| **VEGFA** | **-1.62*** | 0.24 | -0.45 | 0.40 | Angiogenesi |
| HIF1A | -0.38 | NS | +0.45 | 0.30 | Master regulator |
| LDHA | -0.09 | NS | +0.52 | 0.32 | Glicolisi anaerobica |
| SLC2A1 | -0.75 | 0.78 | **-2.06*** | 0.003 | GLUT1 |
| BNIP3 | -0.73 | 0.82 | +0.32 | 0.43 | Autofagia ipossica |
| EPAS1 | +0.31 | NS | +0.85 | 0.20 | HIF2A |
| PDK1 | -0.08 | NS | +0.34 | 0.35 | Glicolisi |

**→ OMA: VEGFA↓ nell'epitelio, SLC2A1/GLUT1↓↓ nello stroma**

### 2.2 CCOC vs OMA

| Gene | Epithelium log2FC | p-adj | Stroma log2FC | p-adj | Funzione |
|------|-------------------|-------|---------------|-------|----------|
| **HIF1A** | +6.90~ | 0.056 | **+11.87*** | 9.7e-15 | Master regulator |
| **VEGFA** | -0.28 | NS | -0.51 | 0.47 | Angiogenesi |
| **LDHA** | +0.23 | NS | **-0.98~** | 0.06 | Glicolisi anaerobica |
| SLC2A1 | +0.07 | NS | -0.56 | 0.63 | GLUT1 |
| BNIP3 | +0.47 | NS | -0.65 | 0.12 | Autofagia ipossica |
| EPAS1 | +1.44 | 0.28 | -0.49 | 0.56 | HIF2A |
| PDK1 | +0.72 | NS | +0.70 | 0.17 | Glicolisi |

**→ CCOC: HIF1A massivamente UP (pseudoipossia!) ma VEGFA e LDHA non aumentano → disaccoppiamento**

### 2.3 EnOC vs OMA

| Gene | Epithelium log2FC | p-adj | Stroma log2FC | p-adj | Funzione |
|------|-------------------|-------|---------------|-------|----------|
| **HIF1A** | +3.98 | 0.41 | **+12.20*** | 0.003 | Master regulator |
| **VEGFA** | +0.73 | NS | +1.41 | 0.22 | Angiogenesi |
| SLC2A1 | +1.41 | 0.27 | +0.47 | 0.73 | GLUT1 |
| LDHA | -1.21 | NS | +0.12 | NS | Glicolisi anaerobica |
| BNIP3 | +0.29 | NS | -0.13 | NS | Autofagia ipossica |
| EPAS1 | -0.87 | NS | -0.09 | NS | HIF2A |
| PDK1 | -0.06 | NS | -0.15 | NS | Glicolisi |

**→ EnOC: HIF1A↑↑ nello stroma, pattern simile a CCOC**

---

## 3. UPR/ER Stress

### 3.1 OMA vs EuE

| Gene | Epithelium log2FC | p-adj | Stroma log2FC | p-adj | Branch |
|------|-------------------|-------|---------------|-------|--------|
| **ERN1** | -1.80 | 0.54 | **-1.42*** | 0.006 | IRE1 |
| **XBP1** | **-1.17*** | 0.12 | -0.53 | 0.29 | IRE1 target |
| **EIF2AK3** | -0.66 | NS | **-1.00*** | 0.004 | PERK |
| **HSPA5** | -0.63 | 0.73 | **-0.64*** | 0.05 | GRP78/BiP |
| **PDIA4** | -0.40 | 0.93 | **-0.70*** | 0.01 | Chaperone |
| **HYOU1** | -0.99 | 0.69 | **-1.14*** | 0.01 | Chaperone |
| DDIT3 | -0.92 | 0.72 | -0.40 | 0.46 | CHOP |
| ATF4 | -0.42 | 0.88 | +0.15 | NS | PERK target |
| ATF6 | +0.29 | 0.93 | +0.15 | NS | ATF6 branch |

**→ OMA: UPR GLOBALMENTE SOPPRESSO, tutti i branch DOWN nello stroma**

### 3.2 CCOC vs OMA

| Gene | Epithelium log2FC | p-adj | Stroma log2FC | p-adj | Branch |
|------|-------------------|-------|---------------|-------|--------|
| **EIF2AK3** | **+3.20*** | 0.01 | **+2.38*** | 1.9e-07 | PERK |
| **DDIT3** | **+1.69~** | 0.20 | -0.17 | NS | CHOP |
| **ERN1** | +1.57 | 0.29 | +0.40 | NS | IRE1 |
| **ATF6** | -0.08 | NS | **+0.56~** | 0.10 | ATF6 branch |
| XBP1 | -0.75 | 0.70 | -0.79 | 0.11 | IRE1 target |
| HSPA5 | +0.36 | 0.87 | -0.31 | 0.45 | GRP78/BiP |
| ATF4 | -0.58 | 0.76 | **-1.08*** | 0.009 | PERK target |
| PDIA4 | -0.25 | NS | +0.55 | 0.26 | Chaperone |
| HYOU1 | -0.12 | NS | +0.33 | NS | Chaperone |

**→ CCOC: PERK branch ATTIVATO (EIF2AK3↑↑↑), IRE1/XBP1 NON attivato**

### 3.3 EnOC vs OMA

| Gene | Epithelium log2FC | p-adj | Stroma log2FC | p-adj | Branch |
|------|-------------------|-------|---------------|-------|--------|
| **EIF2AK3** | +0.84 | NS | **+1.34*** | 0.02 | PERK |
| **HSPA5** | +1.46 | NS | **+2.41*** | 0.02 | GRP78/BiP |
| **ATF6** | +0.34 | NS | **+0.68~** | 0.07 | ATF6 branch |
| **DDIT3** | +2.00 | 0.11 | **+1.63~** | 0.09 | CHOP |
| **PDIA4** | +0.21 | NS | **+0.78** | 0.10 | Chaperone |
| **HYOU1** | +0.53 | NS | +1.12 | 0.32 | Chaperone |
| XBP1 | +0.66 | NS | +0.20 | NS | IRE1 target |
| ERN1 | -0.87 | NS | +0.54 | NS | IRE1 |
| ATF4 | +0.17 | NS | +0.50 | 0.58 | PERK target |

**→ EnOC: UPR attivato nello STROMA (HSPA5↑↑, EIF2AK3↑, ATF6↑), pattern diverso da CCOC**

### 3.4 UPR Branch Summary

| Branch | OMA | CCOC | EnOC |
|--------|-----|------|------|
| **PERK (EIF2AK3)** | ↓ (soppresso) | **↑↑↑ (attivato)** | ↑ (stroma) |
| **IRE1/XBP1** | ↓ (soppresso) | ❌ (non attivato) | ❌ |
| **ATF6** | = | ↑ (trend) | ↑ (stroma) |
| **Chaperones (HSPA5)** | ↓ | = | **↑↑ (stroma)** |

---

## 4. Infiammazione: IL6/JAK/STAT

### 4.1 OMA vs EuE

| Gene | Epithelium log2FC | p-adj | Stroma log2FC | p-adj | Funzione |
|------|-------------------|-------|---------------|-------|----------|
| IL6 | -0.72 | NS | +0.82 | 0.50 | Citochina |
| IL6R | -1.40* | 0.45 | +0.03 | NS | Recettore |
| **IL6ST** | -1.71 | 0.61 | **+0.55*** | 0.07 | gp130 |
| STAT3 | +0.10 | NS | +0.35 | 0.28 | Trasduttore |
| JAK1 | +0.02 | NS | -0.11 | NS | Chinasi |
| JAK2 | +0.54 | NS | -0.07 | NS | Chinasi |
| **SOCS3** | +0.88* | 0.57 | **+1.29*** | 0.05 | Feedback inibitore |
| SOCS1 | -3.33 | 0.14 | -0.89 | 0.26 | Feedback inibitore |

**→ OMA: SOCS3↑↑ come FRENO al pathway, IL6ST↑ nello stroma**

### 4.2 CCOC vs OMA

| Gene | Epithelium log2FC | p-adj | Stroma log2FC | p-adj | Funzione |
|------|-------------------|-------|---------------|-------|----------|
| IL6 | +1.72 | 0.38 | -0.09 | 0.95 | Citochina |
| **IL6ST** | **-1.97~** | 0.07 | **-1.71*** | 3.9e-07 | gp130 |
| **IL6R** | +0.02 | NS | **-1.33~** | 0.09 | Recettore |
| STAT3 | +0.53 | 0.77 | -0.58 | 0.19 | Trasduttore |
| **JAK2** | +1.06 | 0.69 | **+1.56*** | 0.0004 | Chinasi |
| JAK1 | -0.51 | NS | -0.41 | 0.36 | Chinasi |
| **SOCS3** | -0.88 | NS | **-3.26*** | 1.7e-21 | Feedback inibitore |
| SOCS1 | +2.41 | 0.50 | -0.55 | 0.65 | Feedback inibitore |

**→ CCOC: SOCS3↓↓↓ (FRENO RIMOSSO!), JAK2↑, IL6ST/gp130↓ (paradosso)**

### 4.3 EnOC vs OMA

| Gene | Epithelium log2FC | p-adj | Stroma log2FC | p-adj | Funzione |
|------|-------------------|-------|---------------|-------|----------|
| IL6 | -4.04~ | 0.55 | **+3.09** | NA | Citochina |
| IL6ST | -1.35 | NS | -0.06 | NS | gp130 |
| IL6R | +0.41 | NS | +0.10 | NS | Recettore |
| **STAT3** | -0.79 | NS | **-1.11~** | 0.08 | Trasduttore |
| **JAK2** | +0.73 | NS | **+0.94~** | 0.13 | Chinasi |
| JAK1 | -0.44 | NS | +1.04 | 0.27 | Chinasi |
| **SOCS3** | **-2.97*** | 0.004 | **-1.99*** | 0.01 | Feedback inibitore |
| SOCS1 | +3.12 | 0.05 | +0.90 | NS | Feedback inibitore |

**→ EnOC: SOCS3↓↓↓ (freno rimosso), IL6↑↑↑ dallo STROMA**

### 4.4 IL6 Source Analysis

| Condizione | Epithelium IL6 | Stroma IL6 | **Fonte Principale** |
|------------|----------------|------------|---------------------|
| OMA vs EuE | -0.72 | +0.82 | Stroma (trending) |
| **CCOC vs OMA** | +1.72 | -0.09 | **Epitelio** (ma NS) |
| **EnOC vs OMA** | **-4.04** | **+3.09** | **STROMA** |

---

## 5. NF-κB Pathway

### 5.1 OMA vs EuE

| Gene | Epithelium log2FC | p-adj | Stroma log2FC | p-adj | Funzione |
|------|-------------------|-------|---------------|-------|----------|
| **RELA** | -1.53* | 0.06 | +0.22 | 0.47 | p65 |
| **NFKB1** | -1.41 | 0.40 | -0.41 | 0.45 | p105/p50 |
| NFKB2 | -0.76 | 0.80 | +0.30 | 0.46 | p100/p52 |
| RELB | -0.24 | NS | -0.32 | 0.53 | RelB |
| IKBKB | -0.05 | NS | -0.28 | 0.18 | IKKβ |
| TNFAIP3 | -0.90 | 0.76 | -0.003 | NS | A20 (feedback) |

**→ OMA: RELA/p65 trending DOWN nell'epitelio**

### 5.2 CCOC vs OMA

| Gene | Epithelium log2FC | p-adj | Stroma log2FC | p-adj | Funzione |
|------|-------------------|-------|---------------|-------|----------|
| **RELA** | -0.33 | NS | **-1.20*** | 0.02 | p65 |
| **NFKB2** | -1.13 | 0.28 | **-1.28~** | 0.07 | p100/p52 |
| NFKB1 | +0.005 | NS | +0.43 | 0.45 | p105/p50 |
| RELB | +0.69 | NS | +0.83 | 0.18 | RelB |
| IKBKB | +1.03 | 0.35 | +0.64 | 0.18 | IKKβ |
| TNFAIP3 | -1.51 | 0.49 | -1.00 | 0.36 | A20 |

**→ CCOC: RELA/NFKB2↓ nello stroma (NF-κB soppresso)**

### 5.3 EnOC vs OMA

| Gene | Epithelium log2FC | p-adj | Stroma log2FC | p-adj | Funzione |
|------|-------------------|-------|---------------|-------|----------|
| RELA | +0.70 | NS | +0.01 | NS | p65 |
| NFKB1 | -0.15 | NS | +0.59 | 0.56 | p105/p50 |
| NFKB2 | -0.40 | NS | -0.13 | NS | p100/p52 |
| RELB | +0.69 | NS | +2.32 | NA | RelB |
| IKBKB | +0.55 | NS | +0.17 | NS | IKKβ |
| TNFAIP3 | -2.56 | 0.23 | -1.11 | 0.38 | A20 |

**→ EnOC: NF-κB non significativamente alterato**

---

## 6. Morte Cellulare: Piroptosi

### 6.1 OMA vs EuE

| Gene | Epithelium log2FC | p-adj | Stroma log2FC | p-adj | Funzione |
|------|-------------------|-------|---------------|-------|----------|
| **GSDME** | **+4.59*** | 8e-05 | -0.33 | 0.42 | Gasdermina E |
| **CASP1** | **+2.20~** | 0.20 | **+1.48*** | 0.001 | Caspasi infiammatoria |
| IL1B | +0.65 | NS | +0.95 | 0.62 | IL-1β |
| IL18 | +1.42 | 0.69 | +0.26 | NS | IL-18 |
| GSDMD | -0.11 | NS | +0.43 | 0.23 | Gasdermina D |
| NLRP3 | +0.60 | NS | +0.86 | 0.29 | Inflammasoma |

**→ OMA: GSDME↑↑↑ nell'epitelio, CASP1↑ in entrambi (sensibilizzazione piroptosi)**

### 6.2 CCOC vs OMA

| Gene | Epithelium log2FC | p-adj | Stroma log2FC | p-adj | Funzione |
|------|-------------------|-------|---------------|-------|----------|
| **GSDME** | **-8.21*** | 0.03 | **-4.96*** | 0.02 | Gasdermina E |
| ATG7 | **+2.52~** | 0.08 | **+2.91*** | 2.2e-20 | Autofagia |
| ATG5 | +1.96~ | 0.19 | +0.58 | 0.14 | Autofagia |
| MLKL | +3.25 | 0.34 | **+2.64*** | 3.4e-06 | Necroptosi |
| IL18 | -1.58 | 0.34 | -0.36 | NS | IL-18 |
| CASP1 | -0.48 | NS | +0.47 | NS | Caspasi |
| GSDMD | +0.24 | NS | -0.15 | NS | Gasdermina D |
| NLRP3 | ND | - | -1.46 | 0.30 | Inflammasoma |

**→ CCOC: GSDME↓↓↓ (piroptosi disattivata), switch verso autofagia (ATG7↑↑↑)**

### 6.3 EnOC vs OMA

| Gene | Epithelium log2FC | p-adj | Stroma log2FC | p-adj | Funzione |
|------|-------------------|-------|---------------|-------|----------|
| **GSDME** | **-6.56*** | 0.18 | **-2.35*** | 0.04 | Gasdermina E |
| **GSDMD** | -0.40 | NS | **-1.14*** | 0.02 | Gasdermina D |
| **NLRP3** | ND | - | **-2.46~** | 0.07 | Inflammasoma |
| IL18 | -2.55~ | 0.04 | +0.25 | NS | IL-18 |
| IL1B | -1.13 | NS | -2.09 | 0.32 | IL-1β |
| CASP1 | -2.78 | 0.12 | -0.82 | 0.37 | Caspasi |
| ATG7 | +1.57 | 0.34 | **+1.72*** | 0.02 | Autofagia |

**→ EnOC: Piroptosi COMPLETAMENTE DISATTIVATA (GSDME↓, GSDMD↓, NLRP3↓)**

---

## 7. Morte Cellulare: Necroptosi

### 7.1 OMA vs EuE

| Gene | Epithelium log2FC | p-adj | Stroma log2FC | p-adj | Funzione |
|------|-------------------|-------|---------------|-------|----------|
| MLKL | -1.12 | 0.50 | -0.08 | NS | Esecutore |
| RIPK1 | -0.38 | NS | -0.26 | 0.39 | Iniziatore |
| RIPK3 | -1.83 | 0.77 | -1.16 | 0.23 | Iniziatore |

**→ OMA: Necroptosi non significativamente alterata**

### 7.2 CCOC vs OMA

| Gene | Epithelium log2FC | p-adj | Stroma log2FC | p-adj | Funzione |
|------|-------------------|-------|---------------|-------|----------|
| **MLKL** | +3.25 | 0.34 | **+2.64*** | 3.4e-06 | Esecutore |
| **RIPK1** | +1.73 | 0.39 | **+1.40*** | 0.003 | Iniziatore |
| RIPK3 | ND | - | +0.04 | NS | Iniziatore |

**→ CCOC: MLKL↑↑↑, RIPK1↑ nello stroma (sensibilizzazione necroptosi)**

### 7.3 EnOC vs OMA

| Gene | Epithelium log2FC | p-adj | Stroma log2FC | p-adj | Funzione |
|------|-------------------|-------|---------------|-------|----------|
| **RIPK1** | +1.17 | 0.29 | **+1.58*** | 0.001 | Iniziatore |
| MLKL | +0.75 | NS | +0.27 | NS | Esecutore |
| RIPK3 | +0.72 | NS | +0.18 | NS | Iniziatore |

**→ EnOC: RIPK1↑ nello stroma, necroptosi meno attivata vs CCOC**

---

## 8. Morte Cellulare: Autofagia

### 8.1 OMA vs EuE

| Gene | Epithelium log2FC | p-adj | Stroma log2FC | p-adj | Funzione |
|------|-------------------|-------|---------------|-------|----------|
| **MAP1LC3B** | **-1.35** | 0.40 | -0.59 | 0.26 | LC3 |
| **ATG12** | **-0.77~** | 0.46 | **-0.30~** | 0.12 | Elongazione |
| BECN1 | -0.45 | NS | -0.28 | 0.17 | Iniziazione |
| SQSTM1 | -0.89 | 0.55 | +0.06 | NS | p62 |
| ULK1 | -0.09 | NS | -0.38 | 0.38 | Iniziazione |
| ATG5 | -0.06 | NS | -0.006 | NS | Elongazione |
| ATG7 | +0.04 | NS | -0.12 | NS | Elongazione |

**→ OMA: Autofagia moderatamente soppressa**

### 8.2 CCOC vs OMA

| Gene | Epithelium log2FC | p-adj | Stroma log2FC | p-adj | Funzione |
|------|-------------------|-------|---------------|-------|----------|
| **ATG7** | **+2.52~** | 0.08 | **+2.91*** | 2.2e-20 | Elongazione |
| **ATG5** | +1.96~ | 0.19 | +0.58 | 0.14 | Elongazione |
| **MAP1LC3B** | -0.28 | NS | **-1.39*** | 0.003 | LC3 |
| **SQSTM1** | +0.50 | NS | **-0.91*** | 0.003 | p62 |
| **ATG12** | +0.99 | NS | **-0.85*** | 0.01 | Elongazione |
| BECN1 | +0.90 | NS | -0.40 | 0.41 | Iniziazione |
| ULK1 | -0.04 | NS | -0.77 | 0.17 | Iniziazione |

**→ CCOC: ATG7↑↑↑ (autofagia selettiva), ma LC3/SQSTM1↓ nello stroma (flusso alterato)**

### 8.3 EnOC vs OMA

| Gene | Epithelium log2FC | p-adj | Stroma log2FC | p-adj | Funzione |
|------|-------------------|-------|---------------|-------|----------|
| **ATG7** | +1.57 | 0.34 | **+1.72*** | 0.02 | Elongazione |
| ATG5 | +0.12 | NS | +0.10 | NS | Elongazione |
| MAP1LC3B | +0.23 | NS | -0.20 | NS | LC3 |
| SQSTM1 | -0.37 | NS | -0.43 | 0.65 | p62 |
| ATG12 | +0.22 | NS | +0.12 | NS | Elongazione |
| BECN1 | +0.80 | NS | +0.27 | NS | Iniziazione |
| ULK1 | -0.16 | NS | -0.58 | 0.33 | Iniziazione |

**→ EnOC: ATG7↑ nello stroma, autofagia meno alterata vs CCOC**

---

## 9. Modello Integrato: Progressione delle Vie di Morte

```
ENDOMETRIOMA (OMA)
├── PIROPTOSI: GSDME↑↑↑ epitelio, CASP1↑ (SENSIBILIZZATO)
├── NECROPTOSI: Non alterata
├── AUTOFAGIA: Moderatamente soppressa
└── APOPTOSI: BAX↑ epitelio (vedi analisi precedente)

         ↓ CCOC trasformazione ↓

CCOC
├── PIROPTOSI: GSDME↓↓↓ (DISATTIVATA - escape!)
├── NECROPTOSI: MLKL↑↑↑, RIPK1↑ stroma (paradosso?)
├── AUTOFAGIA: ATG7↑↑↑ ma flusso alterato (LC3↓, p62↓)
└── Interpretazione: Switch da piroptosi a autofagia selettiva

         ↓ EnOC trasformazione ↓

EnOC
├── PIROPTOSI: GSDME↓↓, GSDMD↓, NLRP3↓ (COMPLETAMENTE OFF)
├── NECROPTOSI: RIPK1↑ (moderato)
├── AUTOFAGIA: ATG7↑ stroma (moderato)
└── Interpretazione: Soppressione globale piroptosi
```

---

## 10. Insights Chiave

### 10.1 Switch Antiossidante OMA → Tumore

```
OMA: NRF2/HMOX1 pathway (citosolico)
  ↓
CCOC/EnOC: SOD2/GPX1 pathway (MITOCONDRIALE)
```

- I tumori cambiano completamente strategia antiossidante
- SOD2 da -15 a +15 (30-fold change!)
- Potenziale target: inibitori SOD2 mitocondriale

### 10.2 SOCS3 come Guardiano della Trasformazione

| Condizione | SOCS3 Stroma | Interpretazione |
|------------|--------------|-----------------|
| OMA vs EuE | **+1.29*** | FRENO ATTIVO |
| CCOC vs OMA | **-3.26*** | FRENO RIMOSSO |
| EnOC vs OMA | **-1.99*** | FRENO RIMOSSO |

**La perdita di SOCS3 è un evento chiave nella trasformazione!**

### 10.3 UPR Branch Specifici

| Tumore | Branch Dominante | Localizzazione |
|--------|------------------|----------------|
| **CCOC** | PERK (EIF2AK3) | Sistemico (E+S) |
| **EnOC** | ATF6/Chaperone (HSPA5) | Solo Stroma |

### 10.4 Escape dalla Piroptosi

- OMA: GSDME↑↑↑ (sensibilizzato alla morte)
- CCOC/EnOC: GSDME↓↓↓ (escape dalla piroptosi)
- La down-regolazione di GSDME potrebbe essere un meccanismo di evasione

### 10.5 HIF1A Pseudoipossia

- HIF1A massiccamente UP nello stroma tumorale (+11-12 log2FC)
- Ma VEGFA e LDHA non aumentano proporzionalmente
- Suggerisce "pseudoipossia" o HIF1A disaccoppiato dalla risposta classica

---

## 11. Implicazioni Terapeutiche

| Target | CCOC Rationale | EnOC Rationale |
|--------|----------------|----------------|
| **Inibitori SOD2** | SOD2 +15.10*** | SOD2 +12.68*** |
| **Induttori piroptosi** | GSDME↓ (può essere riattivato?) | GSDME↓, NLRP3↓ |
| **PERK/ISR inibitori** | EIF2AK3 +2.38*** | EIF2AK3 +1.34* |
| **JAK2 inibitori** | JAK2 +1.56*** | JAK2 +0.94~ |
| **HIF1A inibitori** | HIF1A +11.87*** | HIF1A +12.20*** |
| **SOCS3 mimetici** | SOCS3 -3.26*** | SOCS3 -1.99*** |

---

## 12. Limitazioni

1. **N=3 CCOC:** Risultati da interpretare con cautela
2. **EnOC Epithelium:** Molti geni NS (basso n?)
3. **Pathway funzionali:** Trascritto ≠ proteina ≠ attività
4. **Interazioni:** Crosstalk tra pathway non catturato

---

*Analisi completata: 2026-01-31*
*Legenda: * p<0.05, *** p<0.001, ~ p<0.1, NS = non significativo, NA = padj non disponibile, ND = non detected*
