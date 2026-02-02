# Analisi Vie Proliferative e di Sopravvivenza

**Data:** 2026-01-31
**Progressione:** EuE → Endo Pelvica → OMA → CCOC/EnOC
**Compartimenti:** Epitelio e Stroma

---

## 1. Marcatori di Proliferazione

### 1.1 Progressione Epiteliale

| Gene | Pelvica vs EuE | OMA vs EuE | CCOC vs OMA | EnOC vs OMA | Funzione |
|------|----------------|------------|-------------|-------------|----------|
| **MKI67** | +0.00 | -1.60 | +1.75 | +3.18 (NA) | Proliferazione |
| **TOP2A** | +1.26 | -1.57 | **+2.34~** | +2.56 | DNA replicazione |
| **CCNE1** | -2.09 | -2.87 | **+6.47*** | +2.99 | Ciclo G1/S |
| CCND1 | -0.78 | -0.86 | -1.11 | +2.15 | Ciclo G1 |
| PCNA | +0.23 | -0.45 | -0.18 | +0.33 | DNA repair |
| CDK4 | +0.09 | -0.10 | -0.17 | +0.65 | Ciclo G1 |
| CDK6 | +0.63 | **+2.51*** | -1.72 | -0.23 | Ciclo G1 |
| CDK2 | -0.74 | -0.36 | -0.04 | +0.07 | Ciclo G1/S |
| RB1 | -0.07 | -0.01 | +0.57 | +0.41 | TSG |

### 1.2 Progressione Stromale

| Gene | Pelvica vs EuE | OMA vs EuE | CCOC vs OMA | EnOC vs OMA |
|------|----------------|------------|-------------|-------------|
| **MKI67** | **-2.05~** | **-3.06*** | **+4.24*** | **+4.41*** |
| **TOP2A** | **-1.94~** | **-2.77*** | **+4.25*** | **+3.58*** |
| **CCNE1** | **-1.75*** | -1.36*** | **+3.07*** | +1.73~ |
| CCND1 | **-1.09*** | -0.54 | +1.34~ | **+1.71*** |
| PCNA | +0.19 | +0.05 | +0.32 | +0.67 |
| CDK4 | -0.23 | -0.27 | **-1.27*** | -0.41 |
| CDK6 | **-1.02*** | -0.08 | +0.64 | +0.85 |
| RB1 | **-0.54*** | **-0.45~** | +1.09~ | +0.67 |

### 1.3 Pattern Chiave: Proliferazione

```
ENDOMETRIO EUTOPICO (Baseline)
         ↓
ENDOMETRIOSI PELVICA (vs EuE)
├── EPITHELIUM: Proliferazione stabile
└── STROMA: **Proliferazione DOWN** (MKI67↓, TOP2A↓, CCNE1↓)
         ↓
ENDOMETRIOMA (vs EuE)
├── EPITHELIUM: CDK6 UP (+2.51***), altri stabili/down
└── STROMA: **Proliferazione molto DOWN** (MKI67↓↓, TOP2A↓↓)
         ↓ TRASFORMAZIONE ↓
CCOC (vs OMA)
├── EPITHELIUM: **CCNE1 +6.47***↑**, TOP2A +2.34~
└── STROMA: **Riattivazione massiva** (MKI67↑↑, TOP2A↑↑, CCNE1↑)

EnOC (vs OMA)
├── EPITHELIUM: Trend UP (MKI67, TOP2A, CCND1)
└── STROMA: **Riattivazione forte** (MKI67↑↑, TOP2A↑↑)
```

**Insight:** Lo stroma delle lesioni benigne è "quiescente" (MKI67/TOP2A DOWN). La trasformazione maligna riattiva massivamente la proliferazione.

---

## 2. Vie di Sopravvivenza/Apoptosi

### 2.1 Progressione Epiteliale

| Gene | Pelvica vs EuE | OMA vs EuE | CCOC vs OMA | EnOC vs OMA | Funzione |
|------|----------------|------------|-------------|-------------|----------|
| **BAX** | +0.53~ | +0.51~ | -1.47 | -0.45 | Pro-apoptotico |
| **BAK1** | -0.66 | -4.23 | +3.64 (NA) | +2.92 (NA) | Pro-apoptotico |
| BCL2 | -0.56 | -0.92 | +1.24 | +1.81 | Anti-apoptotico |
| BCL2L1 | -0.31 | -1.07 | +1.56 | +1.17 | Anti-apoptotico |
| MCL1 | -0.73~ | -0.45 | +0.14 | -0.23 | Anti-apoptotico |
| **BIRC5** | -0.07 | -3.28 | +2.96 | +2.94 | Survivin |
| XIAP | -0.33 | -0.83 | +1.44 | +0.71 | Anti-apoptotico |

### 2.2 Progressione Stromale

| Gene | Pelvica vs EuE | OMA vs EuE | CCOC vs OMA | EnOC vs OMA |
|------|----------------|------------|-------------|-------------|
| **BAX** | **-0.44~** | -0.37 | **-1.32~** | -0.77 |
| **BAK1** | -0.51~ | **-1.00*** | -1.03 | +0.63 |
| BCL2 | +0.39 | +0.17 | +0.18 | -0.04 |
| **BCL2L1** | **-0.61*** | **-0.60*** | **+1.63*** | +0.57 |
| MCL1 | -0.47 | +0.14 | **-0.96*** | -0.36 |
| **BIRC5** | -1.25 | **-2.73*** | **+4.33*** | +1.54 |
| XIAP | -0.18 | -0.04 | +0.26 | +0.41 |

### 2.3 Pattern Chiave: Apoptosi

```
LESIONI BENIGNE (Endo/OMA vs EuE)
├── EPITHELIUM: BAX trending UP (+0.51~) → sensibilità apoptosi
├── STROMA: BAX/BAK1 DOWN → resistenza apoptosi
└── STROMA: BCL2L1 DOWN → meno anti-apoptosi (paradosso)

TRASFORMAZIONE MALIGNA (Tumori vs OMA)
├── CCOC Epithelium: BAX DOWN (-1.47), BIRC5/Survivin UP (+2.96)
├── CCOC Stroma: **BIRC5 +4.33***↑↑↑**, BCL2L1 UP
├── EnOC: Pattern simile ma meno marcato
└── → RESISTENZA ALL'APOPTOSI
```

---

## 3. Pathway PI3K/AKT/mTOR

### 3.1 Progressione Epiteliale

| Gene | Pelvica vs EuE | OMA vs EuE | CCOC vs OMA | EnOC vs OMA | Funzione |
|------|----------------|------------|-------------|-------------|----------|
| PIK3CA | +0.04 | +0.62 | +0.31 | +0.58 | Oncogene |
| AKT1 | +0.21 | -0.58 | +0.00 | +0.96 | Sopravvivenza |
| AKT2 | +0.01 | -0.20 | +0.72 | -0.03 | Metabolismo |
| **AKT3** | -0.30 | **+2.79*** | -1.38 | -0.83 | Sopravvivenza |
| **MTOR** | -0.84 | -0.32 | +1.06 | +1.07 | Crescita |
| TSC1 | -0.55 | -0.51 | +0.94 | +0.94 | Inibitore mTOR |
| TSC2 | -0.58 | -1.12 | +0.57 | +0.84 | Inibitore mTOR |
| RPTOR | -0.58 | +0.82 | +1.42 | +0.60 | mTORC1 |
| RICTOR | +0.09 | +0.33 | -0.14 | +0.40 | mTORC2 |

### 3.2 Progressione Stromale

| Gene | Pelvica vs EuE | OMA vs EuE | CCOC vs OMA | EnOC vs OMA |
|------|----------------|------------|-------------|-------------|
| PIK3CA | **+0.61*** | **+0.52~** | +0.07 | -0.19 |
| AKT1 | -0.01 | -0.19 | -0.08 | +0.05 |
| AKT2 | **+0.29~** | +0.02 | **-1.63*** | -0.57 |
| AKT3 | +0.26 | **+0.66*** | +0.62 | +0.37 |
| MTOR | **-0.76*** | **-0.83*** | **+1.70*** | **+0.97~** |
| TSC2 | -0.09 | -0.26 | **-0.88~** | **-0.75*** |
| RPTOR | -0.11 | -0.09 | **+2.06*** | **+1.43~** |

### 3.3 Pattern Chiave: PI3K/AKT/mTOR

```
LESIONI BENIGNE
├── EPITHELIUM: AKT3 UP in OMA (+2.79***)
├── STROMA: PIK3CA UP, MTOR DOWN
└── → Attivazione parziale PI3K, mTOR soppresso

TRASFORMAZIONE MALIGNA
├── CCOC Stroma: **mTOR +1.70***↑**, RPTOR +2.06***↑**, TSC2↓
├── EnOC Stroma: mTOR +0.97~, RPTOR +1.43~
└── → RIATTIVAZIONE mTORC1 (RPTOR UP, TSC2 DOWN)
```

---

## 4. Pathway MAPK/ERK

### 4.1 Confronti Chiave

| Gene | OMA vs EuE (Epi) | OMA vs EuE (Stroma) | CCOC vs OMA (Stroma) | EnOC vs OMA (Stroma) |
|------|------------------|---------------------|----------------------|----------------------|
| KRAS | -0.27 | -0.08 | -0.53 | +0.13 |
| BRAF | -0.77 | -0.04 | +0.77 | +0.15 |
| RAF1 | -0.87~ | -0.11 | +0.22 | +0.02 |
| MAP2K1 | -0.69 | -0.01 | +0.24 | **+0.85~** |
| **MAPK3** | -0.04 | -0.12 | **-1.19~** | -0.41 |
| MAPK1 | +0.06 | -0.09 | -0.57 | +0.14 |

**Pattern:** MAPK pathway non mostra alterazioni significative consistenti. MAPK3 DOWN in CCOC stroma.

---

## 5. Pathway WNT/β-catenina

### 5.1 Progressione

| Gene | Pelvica vs EuE (Epi) | OMA vs EuE (Epi) | CCOC vs OMA (Epi) | EnOC vs OMA (Stroma) | Funzione |
|------|----------------------|------------------|-------------------|----------------------|----------|
| CTNNB1 | -0.19 | +0.78*** | -1.24 | -0.54 | β-catenina |
| **MYC** | **-1.73*** | +0.83 | **-2.36~** | **-1.46~** | Target WNT |
| **LEF1** | -0.33 | +0.24 | **-3.34~** | +1.32 | Fattore trascriz. |
| TCF7 | **+1.18*** | -0.25 | **-2.34~** | **-1.84*** | Fattore trascriz. |
| APC | -0.19 | +0.21 | -1.28 | -0.38 | Inibitore |
| AXIN2 | -0.65 | -1.11 | - | -0.96 | Target/Feedback |

### 5.2 Pattern WNT nello Stroma

| Gene | Pelvica vs EuE | OMA vs EuE | CCOC vs OMA | EnOC vs OMA |
|------|----------------|------------|-------------|-------------|
| CTNNB1 | -0.38 | +0.44 | -0.38 | -0.54 |
| **MYC** | -0.12 | +0.58 | **-3.36*** | **-1.46~** |
| **LEF1** | **-2.46*** | **-1.84*** | **+2.95***↑** | +1.32 |
| **TCF7** | **-1.12*** | -0.46 | **-2.11*** | **-1.84*** |

### 5.3 Pattern Chiave: WNT

```
LESIONI BENIGNE
├── LEF1 DOWN nello stroma (endometriosi e OMA)
├── TCF7 DOWN (endo pelvica stroma)
└── → WNT pathway soppresso nelle lesioni benigne

CCOC vs OMA
├── Epithelium: MYC↓, LEF1↓, TCF7↓ (WNT soppresso)
├── Stroma: **LEF1 +2.95***↑** ma MYC/TCF7 DOWN
└── → Pattern paradossale, LEF1 up isolato

EnOC vs OMA
├── TCF7 DOWN in epitelio e stroma
├── MYC DOWN
└── → WNT soppresso
```

---

## 6. Sintesi: Mappa della Trasformazione

### 6.1 EPITELIO

| Pathway | Endo Pelvica | OMA | CCOC | EnOC |
|---------|--------------|-----|------|------|
| Proliferazione | Stabile | CDK6↑ | **CCNE1↑↑↑** | UP trend |
| Apoptosi | BAX↑ (sensibile) | BAX↑ | BAX↓, BIRC5↑ | BAX↓ |
| PI3K/AKT | Stabile | AKT3↑↑ | mTOR trend | mTOR trend |
| WNT | TCF7↑ | CTNNB1↑ | WNT↓ | WNT↓ |

### 6.2 STROMA

| Pathway | Endo Pelvica | OMA | CCOC | EnOC |
|---------|--------------|-----|------|------|
| Proliferazione | **DOWN** | **DOWN↓↓** | **UP↑↑↑** | **UP↑↑** |
| Apoptosi | BCL2L1↓ | BIRC5↓↓ | **BIRC5↑↑↑** | BIRC5↑ |
| PI3K/AKT | mTOR↓ | mTOR↓ | **mTOR↑↑** | mTOR↑ |
| WNT | LEF1↓ | LEF1↓ | LEF1↑ | LEF1 trend |

---

## 7. Modello Integrato

```
ENDOMETRIO EUTOPICO
│
├── ENDOMETRIOSI PELVICA
│   ├── Stroma: Proliferazione ↓↓, WNT ↓
│   └── Epitelio: Stabile, BAX↑ (sensibile)
│
├── ENDOMETRIOMA (OMA)
│   ├── Stroma: Proliferazione ↓↓↓, mTOR ↓, BIRC5↓
│   ├── Epitelio: CDK6↑, CTNNB1↑, AKT3↑
│   └── → STROMA QUIESCENTE, EPITELIO PRIMED
│
└── TRASFORMAZIONE MALIGNA
    │
    ├── CCOC
    │   ├── Stroma: MKI67↑↑, TOP2A↑↑, BIRC5↑↑↑, mTOR↑↑, LEF1↑
    │   ├── Epitelio: CCNE1↑↑↑, BAX↓, WNT↓
    │   └── → STROMA DRIVER, CICLO G1/S ATTIVO
    │
    └── EnOC
        ├── Stroma: MKI67↑↑, TOP2A↑, mTOR↑, TCF7↓
        ├── Epitelio: CCND1 trend↑
        └── → SIMILE MA MENO ESTREMO
```

---

## 8. Implicazioni Terapeutiche

| Target | CCOC Rationale | EnOC Rationale |
|--------|----------------|----------------|
| **Inibitori CDK4/6** | CCNE1↑ (bypass CDK4), non ideale | Potenzialmente efficace |
| **Inibitori mTOR** | mTOR↑, RPTOR↑ nello stroma | mTOR↑ |
| **BH3 mimetics** | BIRC5↑↑↑, BCL2L1↑ | BIRC5↑ |
| **Inibitori survivina** | BIRC5 +4.33*** nello stroma | Meno marcato |

---

## 9. Limitazioni

1. N=3 per CCOC, N=4 per EnOC
2. Trascritti ≠ proteine (es. fosforilazione AKT)
3. Pattern stromali potrebbero riflettere CAF recruitment, non trasformazione

---

*Analisi completata: 2026-01-31*
*Legenda: * p<0.05, *** p<0.001, ~ p<0.1*
