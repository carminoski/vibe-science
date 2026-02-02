# Analisi Completa dei Pathway: Oncogeni, Ormonali, Ferroptosi

**Data:** 2026-01-31
**Confronti:** Progressione da lesioni benigne (Endo Pelvica, OMA) a tumori (CCOC, EnOC)
**Compartimenti:** Epitelio, Stroma, Immune

---

## 1. Oncogeni e Tumor Suppressori

### 1.1 CCOC vs OMA (Trasformazione a Cellule Chiare)

| Gene | Epithelium log2FC | Stroma log2FC | Immune log2FC | Funzione |
|------|-------------------|---------------|---------------|----------|
| **CDKN2A** | +1.54 | **+1.81~** | +0.24 | TSG (p16/p14ARF) |
| **NF1** | +1.30 | **+1.14*** | **+1.79*** | TSG (neurofibromin) |
| TP53 | -0.75 | -0.15 | +0.01 | TSG |
| KRAS | -0.42 | -0.53 | -0.02 | Oncogene |
| ARID1A | -0.37 | -0.51 | -0.17 | TSG (SWI/SNF) |
| PIK3CA | +0.31 | +0.07 | -0.20 | Oncogene |
| PTEN | -0.09 | +0.27 | -0.08 | TSG |

**→ CCOC: NF1 paradossalmente UP (tutti compartimenti); CDKN2A UP in stroma**

### 1.2 EnOC vs OMA (Trasformazione Endometrioide)

| Gene | Epithelium log2FC | Stroma log2FC | Immune log2FC | Funzione |
|------|-------------------|---------------|---------------|----------|
| **CDKN2A** | +1.81 | +2.38 (NA) | -0.49 | TSG (p16/p14ARF) |
| **NF1** | -0.40 | +0.41 | **+2.44*** | TSG |
| TP53 | -0.30 | -0.56 | +0.59 | TSG |
| PTEN | +0.45 | +0.16 | +0.88~ | TSG |
| PIK3CA | +0.58 | -0.19 | +0.40 | Oncogene |
| KRAS | +0.23 | +0.13 | +0.21 | Oncogene |
| ARID1A | -0.27 | -0.42 | +0.32 | TSG |

**→ EnOC: NF1 fortemente UP nel compartimento immune (+2.44***)**

### 1.3 Lesioni Benigne: OMA vs EuE

| Gene | Epithelium | Stroma | Immune | Trend |
|------|------------|--------|--------|-------|
| PTEN | +0.81~ | -0.31 | +0.05 | Stabile |
| NF1 | +0.74~ | +0.27 | -0.21 | Trending UP epitelio |
| PIK3CA | +0.62 | **+0.52*** | +0.24 | UP stroma |
| CDKN2A | +0.84 | -0.33 | -0.70 | Variabile |
| TP53 | -0.21 | **-0.72*** | -0.39 | DOWN stroma |
| ARID1A | -0.07 | -0.18 | -0.07 | Stabile |
| KRAS | -0.27 | -0.08 | +0.28 | Stabile |

**→ OMA: TP53 DOWN nello stroma (permissivo); PIK3CA UP nello stroma**

---

## 2. Pathway Ormonale

### 2.1 CCOC vs OMA

| Gene | Epithelium | Stroma | Immune | Recettore |
|------|------------|--------|--------|-----------|
| **AR** | **-5.10*** | **-1.21*** | -2.13 | Androgeno |
| **ESR2** | **+5.89*** | **+2.96*** | **+1.89*** | Estrogeno β |
| ESR1 | -2.02 | +0.94 | **+3.14*** | Estrogeno α |
| PGR | -1.83 | -2.42 (NA) | -0.15 | Progesterone |

**→ CCOC: Switch ormonale! AR fortemente DOWN, ESR2 fortemente UP**

### 2.2 EnOC vs OMA

| Gene | Epithelium | Stroma | Immune | Recettore |
|------|------------|--------|--------|-----------|
| **PGR** | **+2.99~** | -0.46 | +2.35 | Progesterone |
| **ESR1** | **+2.40~** | +0.69 | **+3.90*** | Estrogeno α |
| ESR2 | +3.15 | **+2.50*** | -0.36 | Estrogeno β |
| AR | +0.70 | -0.24 | +1.72 | Androgeno |

**→ EnOC: ESR1 e PGR UP, pattern estrogeno-dipendente classico**

### 2.3 Lesioni Benigne: OMA vs EuE

| Gene | Epithelium | Stroma | Immune |
|------|------------|--------|--------|
| **ESR1** | -2.29~ | **-2.82*** | **-3.78*** |
| **PGR** | **-3.65*** | **-2.17*** | **-3.70*** |
| ESR2 | -2.77 | **+3.13*** | **+1.88*** |
| AR | +0.88 | -0.59 | **-1.88~** |

**→ OMA: ESR1/PGR DOWN (resistenza ormonale), ESR2 UP nello stroma**

### 2.4 Pattern Ormonale Integrato

```
ENDOMETRIO EUTOPICO
├── ESR1/PGR: Alti (risposta estrogeno/progesterone)
├── ESR2: Basso
└── AR: Presente

         ↓ Progressione a OMA ↓

ENDOMETRIOMA (OMA)
├── ESR1/PGR: **DOWN** (resistenza ormonale)
├── ESR2: **UP stroma** (switch recettoriale)
└── AR: Variabile

         ↓ Trasformazione tumorale ↓

CCOC                           EnOC
├── AR: **-5.10*** (silenziato) ├── ESR1: **+2.40~** (riattivato)
├── ESR2: **+5.89***            ├── PGR: **+2.99~** (riattivato)
└── ESR1: variabile             └── ESR2: UP stroma
```

---

## 3. Ferroptosi e Metabolismo del Ferro

### 3.1 CCOC vs OMA

| Gene | Epithelium | Stroma | Immune | Funzione |
|------|------------|--------|--------|----------|
| **GPX4** | -0.79 | **-2.06*** | **-1.53*** | Anti-ferroptosi |
| **FTH1** | -1.53 | **-1.46*** | **-2.94*** | Storage ferro |
| **FTL** | +0.54 | **-2.18*** | **-2.62*** | Storage ferro |
| **SLC7A11** | +1.50 | **+2.92*** | -2.11 | Anti-ferroptosi (xCT) |
| ACSL4 | **-2.07~** | **+0.78*** | +0.43 | Pro-ferroptosi |
| TFRC | +0.16 | **-1.07~** | - | Uptake ferro |
| NCOA4 | -0.84 | -0.43 | - | Ferritinofagia |

**→ CCOC Stroma: Paradosso ferroptosi! GPX4↓ + FTH1/FTL↓ ma SLC7A11↑ (compensazione)**

### 3.2 EnOC vs OMA

| Gene | Epithelium | Stroma | Immune | Funzione |
|------|------------|--------|--------|----------|
| **ACSL4** | **-1.55*** | **+1.25*** | +0.47 | Pro-ferroptosi |
| **SLC7A11** | -1.41 | **+2.91*** | -1.06 | Anti-ferroptosi |
| FTL | -0.68 | -1.20 | -1.08 | Storage ferro |
| GPX4 | -1.01 | -0.90 | +0.37 | Anti-ferroptosi |
| FTH1 | -1.08 | +0.07 | +0.35 | Storage ferro |
| TFRC | -0.78 | -0.11 | - | Uptake ferro |
| NCOA4 | -0.74 | +0.08 | - | Ferritinofagia |

**→ EnOC: ACSL4 UP stroma (pro-ferroptosi), ma SLC7A11 UP compensa**

### 3.3 Lesioni Benigne: OMA vs EuE

| Gene | Epithelium | Stroma | Immune |
|------|------------|--------|--------|
| GPX4 | -0.56 | -0.11 | -0.21 |
| FTH1 | -0.37 | +0.24 | +0.84 |
| FTL | +0.63 | **+0.70~** | **+1.45~** |
| SLC7A11 | -0.60 | +0.15 | +1.50 |
| ACSL4 | -0.61 | **+0.59~** | +0.26 |

**→ OMA: FTL UP (accumulo ferro), ACSL4 UP stroma (sensibilizzazione)**

### 3.4 Modello Ferroptosi

```
ENDOMETRIOMA (OMA)
├── FTL↑: Accumulo ferro (chocolate cyst)
├── ACSL4↑ stroma: Sensibilizzazione
└── GPX4: Stabile (protezione residua)

         ↓ CCOC pathway ↓

CCOC
├── GPX4↓ + FTH1↓ + FTL↓: Vulnerabilità teorica
├── SLC7A11↑↑: Compensazione (cisteina uptake)
└── Net: Resistenza acquisita

         ↓ EnOC pathway ↓

EnOC
├── ACSL4↑ stroma: Sensibilizzazione persistente
├── SLC7A11↑ stroma: Compensazione parziale
└── Net: Equilibrio instabile
```

---

## 4. Tabella Comparativa: Tumor vs Precursore

### Oncogeni/TSG più Alterati

| Gene | CCOC (max) | EnOC (max) | Compartimento Principale |
|------|------------|------------|--------------------------|
| NF1 | +1.79*** (immune) | +2.44*** (immune) | **IMMUNE** |
| CDKN2A | +1.81~ (stroma) | +2.38 (stroma) | Stroma |
| TP53 | -0.75 | -0.56 | Non significativo |
| PIK3CA | +0.31 | +0.58 | Non significativo |
| ARID1A | -0.51 | -0.42 | Non significativo |

**Nota Critica:** ARID1A e PIK3CA, mutati frequentemente in EAOC, mostrano scarsi cambiamenti trascrizionali. Le mutazioni potrebbero essere loss-of-function senza riduzione espressione.

### Recettori Ormonali

| Gene | CCOC | EnOC | Pattern |
|------|------|------|---------|
| **AR** | **-5.10***↓ | +0.70 | CCOC-specifico |
| **ESR2** | **+5.89***↑ | +3.15 | UP in entrambi |
| **ESR1** | -2.02/+3.14 | **+2.40~**↑ | EnOC riattiva |
| **PGR** | -1.83 | **+2.99~**↑ | EnOC riattiva |

### Ferroptosi

| Gene | CCOC Stroma | EnOC Stroma | Significato |
|------|-------------|-------------|-------------|
| GPX4 | **-2.06***↓ | -0.90 | CCOC vulnerabile |
| SLC7A11 | **+2.92***↑ | **+2.91***↑ | Compensazione |
| ACSL4 | +0.78*↑ | **+1.25***↑ | Pro-ferroptosi |

---

## 5. Insights Chiave

### 5.1 NF1: Paradosso Tumor Suppressor

NF1 è un noto tumor suppressor (neurofibromatosi), ma nel nostro dataset è **UP** in entrambi i tumori, specialmente nel compartimento immune.

**Possibili spiegazioni:**
1. Risposta compensatoria allo stress
2. Up-regulation nei linfociti infiltranti (non nel tumore)
3. Artefatto da proporzioni cellulari

### 5.2 Switch Ormonale CCOC vs EnOC

```
CCOC: AR↓↓↓ + ESR2↑↑↑ → "De-mascolinizzazione"
EnOC: ESR1↑ + PGR↑ → "Ristabilimento sensibilità estrogeno/progesterone"
```

Questo potrebbe spiegare:
- CCOC: Resistenza a terapie ormonali classiche
- EnOC: Potenziale sensibilità a modulatori ormonali

### 5.3 Ferroptosi come Target Terapeutico

**CCOC:**
- GPX4 DOWN nello stroma/immune → vulnerabilità
- SLC7A11 UP → compensazione tramite import cisteina
- **Target:** Inibitori SLC7A11 (erastin, sulfasalazina)

**EnOC:**
- ACSL4 UP nello stroma → sensibilizzazione lipidica
- SLC7A11 UP → protezione
- **Target:** Combinazione ACSL4 agonisti + inibitori SLC7A11

### 5.4 Progressione OMA → Tumore

| Evento | OMA | CCOC | EnOC |
|--------|-----|------|------|
| ESR1/PGR | DOWN | Variabile | **Riattivato** |
| ESR2 stroma | UP | **UP++** | UP |
| AR | Stabile | **DOWN***| Stabile |
| TP53 stroma | DOWN | Stabile | Stabile |
| PIK3CA stroma | UP | Stabile | Stabile |
| FTL | UP | DOWN | DOWN |
| SLC7A11 | Basso | **UP** | **UP** |

---

## 6. Implicazioni Terapeutiche

| Target | CCOC Rationale | EnOC Rationale |
|--------|----------------|----------------|
| **Inibitori SLC7A11** | SLC7A11 +2.92*** compensa GPX4↓ | SLC7A11 +2.91*** |
| **Modulatori ESR2** | ESR2 +5.89*** | ESR2 +3.15 |
| **Anti-estrogeni (ESR1)** | Probabilmente inefficace | ESR1 riattivato (+2.40~) |
| **Progestinici** | PGR basso | PGR riattivato (+2.99~) |
| **Inibitori PERK/ISR** | EIF2AK3 UP (da analisi UPR) | Meno rilevante |
| **Targeting stromal IL6** | IL6 trending | **IL6 +3.09*** stroma |

---

## 7. Limitazioni

1. **N=3 CCOC:** Risultati da interpretare con cautela
2. **NF1 paradosso:** Richiede validazione (IHC, mutazioni)
3. **Ferroptosi funzionale:** Trascritto ≠ attività proteica
4. **ARID1A/PIK3CA:** Mutazioni vs espressione non correlate

---

*Analisi completata: 2026-01-31*
*Legenda: * p<0.05, *** p<0.001, ~ p<0.1, NA = padj non disponibile*
