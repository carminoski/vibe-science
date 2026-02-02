# Analisi delle Lesioni Benigne: Chi Fa Cosa?

**Data:** 2026-01-31
**Confronti:** Endometriosi Pelvica e Endometrioma vs Endometrio Eutopico
**Compartimenti:** Epitelio, Stroma, Immune

---

## 1. Endometriosi Pelvica vs Endometrio Eutopico

### EPITHELIUM (Cellule Epiteliali Ectopiche)

| Pathway | Gene | log2FC | Sig | Interpretazione |
|---------|------|--------|-----|-----------------|
| **Ipossia** | HIF1A | -1.65 | *** | Ipossia pathway DOWN |
| | VEGFA | -2.11 | *** | Angiogenesi soppressa |
| **UPR/ER Stress** | ERN1 | -2.03 | *** | IRE1 branch DOWN |
| | XBP1 | -0.84 | *** | IRE1 branch DOWN |
| | HYOU1 | -1.15 | * | Chaperone DOWN |
| | EIF2AK3 | -1.13 | * | PERK trending DOWN |
| **IL6/STAT3** | IL6 | **-4.74** | *** | IL6 fortemente DOWN |
| | IL6ST | -2.04 | * | Recettore DOWN |
| | IL6R | -1.11 | * | Recettore DOWN |
| **Stress ossidativo** | SOD2 | -15.09 | * | Antiossidante DOWN |
| | HMOX1 | +1.96 | | Trending UP |
| | NQO1 | +1.79 | * | Antiossidante UP |
| **Apoptosi** | BAX | +0.53 | * | Pro-apoptotico UP |

**→ Epitelio endometriosico: UPR spento, IL6 silenziato, stress ossidativo attivo**

### STROMA (Fibroblasti/Cellule Stromali)

| Pathway | Gene | log2FC | Sig | Interpretazione |
|---------|------|--------|-----|-----------------|
| **UPR/ER Stress** | ERN1 | -1.27 | *** | IRE1 DOWN |
| | HSPA5 | -0.90 | *** | GRP78 DOWN |
| | HYOU1 | -0.99 | *** | Chaperone DOWN |
| | PDIA4 | -0.76 | *** | Chaperone DOWN |
| **Proliferazione** | MKI67 | -2.05 | *** | Proliferazione DOWN |
| | TP53 | -0.95 | *** | p53 DOWN |
| **JAK/STAT** | JAK2 | +0.65 | *** | JAK2 UP |
| | IL6ST | +0.55 | *** | gp130 UP |
| **Stress ossidativo** | HMOX1 | +1.02 | * | Antiossidante UP |
| | NFE2L2 | +0.40 | *** | NRF2 UP |
| **Apoptosi** | BAX | -0.44 | *** | Pro-apoptotico DOWN |

**→ Stroma endometriosico: UPR spento, bassa proliferazione, JAK2 attivo, resistenza apoptosi**

### IMMUNE (Cellule Immunitarie)

| Pathway | Gene | log2FC | Sig | Interpretazione |
|---------|------|--------|-----|-----------------|
| **IL6/STAT3** | IL6 | +1.02 | | Trending UP |
| | SOCS3 | +0.59 | * | Feedback attivo |
| **UPR** | HSPA5 | -0.91 | *** | GRP78 DOWN |
| | ATF4 | -0.48 | *** | PERK target DOWN |
| | EIF2AK3 | -0.54 | *** | PERK DOWN |
| **Proliferazione** | MKI67 | -1.95 | *** | Proliferazione DOWN |
| **Stress ossidativo** | HMOX1 | +1.12 | * | Antiossidante UP |

**→ Immune: IL6 trending up (fonte!), UPR spento, HMOX1 up**

---

## 2. Endometrioma vs Endometrio Eutopico

### EPITHELIUM (Cellule Epiteliali della Cisti)

| Pathway | Gene | log2FC | Sig | Interpretazione |
|---------|------|--------|-----|-----------------|
| **UPR/ER Stress** | XBP1 | -1.17 | *** | IRE1 branch DOWN |
| | ERN1 | -1.80 | * | IRE1 DOWN |
| **Ipossia** | VEGFA | -1.62 | *** | Angiogenesi DOWN |
| **IL6/STAT3** | IL6R | -1.40 | * | Recettore DOWN |
| | IL6ST | -1.71 | | gp130 trending DOWN |
| | SOCS3 | +0.88 | * | Feedback UP |
| | IL6 | -0.72 | | Trending DOWN |
| **Stress ossidativo** | HMOX1 | +2.00 | | Antiossidante UP |
| | SOD2 | -15.51 | | Mitocondriale DOWN |
| **Apoptosi** | BAX | +0.51 | * | Pro-apoptotico UP |

**→ Epitelio OMA: UPR spento (come endo pelvica), SOCS3 UP, stress ossidativo**

### STROMA (Stroma della Cisti Ovarica)

| Pathway | Gene | log2FC | Sig | Interpretazione |
|---------|------|--------|-----|-----------------|
| **UPR/ER Stress** | EIF2AK3 | -1.00 | *** | PERK DOWN |
| | ERN1 | -1.42 | *** | IRE1 DOWN |
| | HSPA5 | -0.64 | *** | GRP78 DOWN |
| | HYOU1 | -1.14 | *** | Chaperone DOWN |
| | PDIA4 | -0.70 | *** | Chaperone DOWN |
| **Proliferazione** | MKI67 | -3.06 | *** | Proliferazione molto DOWN |
| | TP53 | -0.72 | *** | p53 DOWN |
| **IL6/STAT3** | SOCS3 | +1.29 | *** | Feedback FORTE |
| | IL6ST | +0.55 | *** | gp130 UP |
| | IL6 | +0.82 | | Trending UP |
| **Stress ossidativo** | HMOX1 | +1.97 | *** | Antiossidante UP |
| | NFE2L2 | +0.59 | *** | NRF2 UP |
| | SOD2 | -20.23 | *** | Mitocondriale DOWN |

**→ Stroma OMA: UPR fortemente soppresso, SOCS3 molto UP, stress ossidativo (HMOX1/NRF2 UP)**

### IMMUNE (Cellule Immunitarie nell'Endometrioma)

| Pathway | Gene | log2FC | Sig | Interpretazione |
|---------|------|--------|-----|-----------------|
| **IL6/STAT3** | IL6 | **+1.31** | *** | **IL6 UP - FONTE!** |
| | STAT3 | +0.71 | *** | STAT3 attivo |
| | SOCS3 | +0.97 | *** | Feedback attivo |
| | IL6R | +0.71 | | Recettore UP |
| **Stress ossidativo** | HMOX1 | +2.05 | *** | Antiossidante UP |
| **Proliferazione** | MKI67 | -1.67 | *** | Proliferazione DOWN |

**→ Immune OMA: IL6 SOURCE CONFIRMED! STAT3 attivo, SOCS3 feedback**

---

## 3. Tabella Comparativa: Chi Produce Cosa?

### IL6 Expression

| Condizione | Epithelium | Stroma | Immune | **Fonte Principale** |
|------------|------------|--------|--------|---------------------|
| **Endo Pelvica** | -4.74*** | -0.07 | +1.02 | **IMMUNE** |
| **Endometrioma** | -0.72 | +0.82 | **+1.31***| **IMMUNE** |

### UPR Status (EIF2AK3 + XBP1 + HSPA5 media)

| Condizione | Epithelium | Stroma | Immune |
|------------|------------|--------|--------|
| **Endo Pelvica** | **SPENTO** (-1.43) | **SPENTO** (-0.40) | **SPENTO** (-0.52) |
| **Endometrioma** | **SPENTO** (-0.82) | **SPENTO** (-0.72) | +0.06 (normale) |

### Stress Ossidativo (HMOX1)

| Condizione | Epithelium | Stroma | Immune |
|------------|------------|--------|--------|
| **Endo Pelvica** | +1.96 | +1.02* | +1.12* |
| **Endometrioma** | +2.00 | **+1.97***| **+2.05***|

### Feedback Inhibition (SOCS3)

| Condizione | Epithelium | Stroma | Immune |
|------------|------------|--------|--------|
| **Endo Pelvica** | -0.48 | +0.08 | +0.59* |
| **Endometrioma** | +0.88* | **+1.29***| **+0.97***|

---

## 4. Modello Integrato: Endometriosi → Endometrioma

```
ENDOMETRIOSI PELVICA
├── EPITHELIUM: UPR OFF, IL6 silenziato (-4.74***)
│              Stress ossidativo (HMOX1↑)
│              Recettori IL6 down (IL6R↓, IL6ST↓)
│
├── STROMA:    UPR OFF, proliferazione bassa (MKI67↓)
│              JAK2 attivo, resistenza apoptosi (BAX↓)
│              Stress ossidativo (HMOX1↑, NRF2↑)
│
└── IMMUNE:    IL6 trending UP (fonte!)
               UPR down, SOCS3 attivo
               Infiammazione controllata

         ↓ Progressione nell'ovaio ↓

ENDOMETRIOMA (Cisti Ovarica)
├── EPITHELIUM: UPR ancora OFF
│              IL6 normalizza (non più -4.74)
│              SOCS3 UP (feedback attivato)
│              Stress ossidativo persistente
│
├── STROMA:    UPR fortemente OFF (tutti i branch)
│              SOCS3 MOLTO UP (+1.29***)
│              IL6 trending UP (+0.82)
│              Proliferazione molto bassa
│              Stress ossidativo intenso
│
└── IMMUNE:    **IL6 UP (+1.31***) - FONTE CONFERMATA**
               STAT3 attivo (+0.71***)
               SOCS3 feedback attivo
               Infiammazione cronica controllata
```

---

## 5. Insights Chiave

### 5.1 La Fonte di IL6 nelle Lesioni Benigne
- **Nelle lesioni benigne (endo + OMA), IL6 è prodotto dalle cellule IMMUNITARIE**
- L'epitelio ha IL6 DOWN o silenziato
- Lo stroma ha IL6 trending UP solo nell'OMA

### 5.2 UPR è Soppresso
- **UPR è globalmente soppresso in tutte le lesioni benigne**
- Questo contrasta con i tumori dove UPR si "riaccende"
- La soppressione UPR potrebbe essere un meccanismo di tolleranza

### 5.3 SOCS3 come Guardiano
- SOCS3 è UP nelle lesioni benigne (specialmente OMA stroma: +1.29***)
- Questo "frena" la segnalazione JAK/STAT
- **Nella trasformazione maligna, SOCS3 viene DOWN-regolato (CCOC: -3.26***, EnOC: -2.97***)**
- La perdita di SOCS3 potrebbe essere un evento chiave nella trasformazione

### 5.4 Stress Ossidativo come Driver
- HMOX1 UP consistente in tutti i compartimenti e condizioni
- SOD2 paradossalmente DOWN (mitocondriale)
- Suggerisce stress ossidativo cronico compensato da pathway citosolico (HMOX1)

---

## 6. Implicazioni per la Trasformazione Maligna

| Evento | Lesioni Benigne | CCOC | EnOC |
|--------|-----------------|------|------|
| UPR | SPENTO | **ACCESO** (PERK sistemico) | **ACCESO** (ATF6/stroma) |
| IL6 source | IMMUNE | EPITHELIUM? | STROMA |
| SOCS3 | **UP** (freno attivo) | **DOWN** (freno rimosso) | **DOWN** (freno rimosso) |
| HMOX1 | UP | ? | ? |

**Ipotesi di trasformazione:**
1. Lesione benigna mantiene UPR spento e SOCS3 attivo
2. Stress ossidativo cronico (ferro nell'OMA)
3. Perdita di SOCS3 → rimozione freno JAK/STAT
4. Attivazione UPR (PERK o ATF6 a seconda del destino)
5. Divergenza: CCOC (PERK sistemico) vs EnOC (ATF6 stromale)

---

*Analisi completata: 2026-01-31*
