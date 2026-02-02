# Analisi Compartment-Specific: UPR e IL6 Pathway

**Data:** 2026-01-31
**Obiettivo:** Verificare pattern UPR e IL6 separatamente per Epitelio, Stroma, Immune

---

## 1. Precursore: OMA vs EuE (Endometrioma vs Endometrio Eutopico)

| Gene | Epithelium | Stroma | Immune |
|------|------------|--------|--------|
| **PERK Branch** | | | |
| EIF2AK3 | -0.66 | **-1.00*** | +0.12 |
| ATF4 | -0.42 | +0.15 | -0.02 |
| DDIT3 | -0.92 | -0.40 | +0.08 |
| **IRE1 Branch** | | | |
| ERN1 | **-1.80~** | **-1.42*** | +0.44 |
| XBP1 | **-1.17*** | -0.53 | +0.44~ |
| **ATF6 Branch** | | | |
| ATF6 | +0.29 | +0.15 | +0.23 |
| HSPA5 | -0.63 | **-0.64*** | **-0.38*** |
| **IL6/STAT3** | | | |
| IL6 | -0.72 | +0.82 | **+1.31*** |
| STAT3 | +0.10 | +0.35 | **+0.71*** |
| SOCS3 | +0.88~ | **+1.29*** | **+0.97*** |

**Interpretazione OMA:**
- UPR globalmente **SOPPRESSO** in epitelio/stroma (XBP1↓, ERN1↓, HSPA5↓)
- IL6/STAT3 **ATTIVO solo nel compartimento IMMUNE** (+1.31*, +0.71*)
- SOCS3 UP suggerisce feedback negativo attivo

---

## 2. Trasformazione CCOC: CCOC vs OMA

| Gene | Epithelium | Stroma | Immune |
|------|------------|--------|--------|
| **PERK Branch** | | | |
| EIF2AK3 | **+3.20*** | **+2.38*** | **+2.02*** |
| ATF4 | -0.58 | **-1.08*** | **-0.73*** |
| DDIT3 | **+1.69*** | -0.17 | **-1.14*** |
| **IRE1 Branch** | | | |
| ERN1 | **+1.57~** | +0.40 | -0.49 |
| XBP1 | -0.75 | **-0.79*** | +0.37 |
| **ATF6 Branch** | | | |
| ATF6 | -0.08 | **+0.56*** | +0.69 |
| HSPA5 | +0.36 | -0.31 | **-0.61*** |
| **IL6/STAT3** | | | |
| IL6 | +1.72 | -0.09 | -1.74 |
| STAT3 | +0.53 | -0.58~ | -0.43 |
| SOCS3 | -0.88 | **-3.26*** | **-1.93*** |

**Interpretazione CCOC:**
- **EIF2AK3/PERK UP in TUTTI i compartimenti** (pattern "sistemico")
- DDIT3 UP solo in epitelio (pathway pro-apoptotico attenuato in stroma/immune)
- **XBP1 NON attivato** (conferma critica R2)
- IL6 trending UP solo in epitelio (+1.72, NS)
- **SOCS3 DOWN** in stroma/immune → rimozione feedback inibitorio

---

## 3. Trasformazione EnOC: EnOC vs OMA

| Gene | Epithelium | Stroma | Immune |
|------|------------|--------|--------|
| **PERK Branch** | | | |
| EIF2AK3 | +0.84 | **+1.34*** | +0.37 |
| ATF4 | +0.17 | +0.50 | +0.20 |
| DDIT3 | +2.00 | **+1.63*** | +0.69 |
| **IRE1 Branch** | | | |
| ERN1 | -0.87 | +0.54 | -0.19 |
| XBP1 | +0.66 | +0.20 | -0.34 |
| **ATF6 Branch** | | | |
| ATF6 | +0.34 | **+0.68*** | +1.04 |
| HSPA5 | +1.46 | **+2.41*** | +0.46 |
| **IL6/STAT3** | | | |
| IL6 | **-4.04*** | **+3.09*** | -0.37 |
| STAT3 | -0.79 | **-1.11*** | **-1.05*** |
| SOCS3 | **-2.97*** | **-1.99*** | **-2.78*** |

**Interpretazione EnOC:**
- UPR attivo principalmente nello **STROMA** (HSPA5↑, DDIT3↑, ATF6↑)
- **IL6 pattern OPPOSTO per compartimento:**
  - Epitelio: **-4.04*** (DOWN!)
  - Stroma: **+3.09*** (UP!)
- **SOCS3 DOWN in tutti i compartimenti** → rimozione feedback globale
- STAT3 DOWN in stroma/immune

---

## 4. Confronto Diretto: CCOC vs EnOC

### Pattern UPR

| Aspetto | CCOC | EnOC |
|---------|------|------|
| PERK attivazione | **Sistemica** (tutti i compartimenti) | **Solo stroma** |
| Branch dominante | PERK (EIF2AK3, DDIT3 in epi) | ATF6/chaperones (HSPA5 in stroma) |
| XBP1/IRE1 | Non attivo | Non attivo |

### Pattern IL6

| Compartimento | CCOC vs OMA | EnOC vs OMA | Differenza |
|---------------|-------------|-------------|------------|
| **Epithelium** | +1.72 (NS) | **-4.04*** | **OPPOSTO** |
| **Stroma** | -0.09 (NS) | **+3.09*** | **OPPOSTO** |
| Immune | -1.74 (NS) | -0.37 (NS) | Simile |

---

## 5. Insight Chiave

### 5.1 IL6: Sorgente Cellulare Diversa

```
CCOC: IL6 → trending UP in EPITHELIUM (+1.72)
      → source: possibly tumor-intrinsic

EnOC: IL6 → DOWN in epithelium (-4.04***)
      → UP in STROMA (+3.09***)
      → source: STROMA-DERIVED (TME)
```

**Questo risponde alla critica R2 sull'attribuzione cellulare:**
- In EnOC, IL6 NON è tumor-intrinsic, è stroma-derived
- In CCOC, il segnale epiteliale è presente ma debole (N=3, outlier issue)

### 5.2 UPR: Branch Diversi

```
CCOC: PERK → ISR attivato in modo sistemico
      EIF2AK3 UP in tutti i compartimenti
      → stress response coordinato

EnOC: ATF6/chaperone → attivo nello STROMA
      HSPA5/GRP78 UP (+2.41***)
      → ER stress nello stroma, non nel tumore
```

### 5.3 SOCS3: Feedback Removal

```
Entrambi i tumori: SOCS3 DOWN
CCOC: -3.26*** (stroma), -1.93*** (immune)
EnOC: -2.97*** (epi), -1.99*** (stroma), -2.78*** (immune)

→ Rimozione del feedback inibitorio JAK/STAT
→ Ma senza IL6 UP corrispondente, il pathway potrebbe essere
   attivato da altri fattori (IFN, IL-10, etc.)
```

---

## 6. Modello Rivisto

### CCOC Pathway (Rivisto)
```
Endometrioma microenvironment
         ↓
    PERK activation (sistemica)
         ↓
    ├── Epithelium: EIF2AK3↑ + DDIT3↑ + IL6 trending↑
    ├── Stroma: EIF2AK3↑ + SOCS3↓
    └── Immune: EIF2AK3↑ + SOCS3↓
         ↓
    ISR program attivo
         ↓
    Sopravvivenza (ma meccanismo IL6→STAT3 non confermato)
```

### EnOC Pathway (Rivisto)
```
Endometrioma microenvironment
         ↓
    ATF6/chaperone response (stroma-specifico)
         ↓
    ├── Epithelium: IL6↓ + SOCS3↓ + OxPhos/Myc↑ (GSEA)
    ├── Stroma: HSPA5↑ + IL6↑ (source!)
    └── Immune: SOCS3↓
         ↓
    Stroma-derived IL6 (paracrine?)
         ↓
    Metabolic reprogramming (OxPhos, Myc targets)
```

---

## 7. Conclusioni

1. **IL6 source è compartment-specific:**
   - CCOC: epiteliale (ma outlier-sensitive)
   - EnOC: STROMA-derived

2. **UPR branch è tumor-type specific:**
   - CCOC: PERK sistemico
   - EnOC: ATF6/chaperone nello stroma

3. **La critica R2 è validata:** non si può parlare di "tumor-intrinsic IL6" senza specificare il compartimento

4. **Nuovo finding:** EnOC ha pattern "paracrine" (stroma→tumor) mentre CCOC potrebbe avere pattern "autocrine" (tumor-intrinsic, ma N=3 limita la confidenza)

---

## 8. Implicazioni Terapeutiche (Riviste)

| Subtype | Target | Rationale |
|---------|--------|-----------|
| CCOC | PERK/ISR inhibitors | EIF2AK3 UP sistemico |
| CCOC | ISRIB (eIF2α) | ISR program |
| EnOC | TME-targeting | IL6 stroma-derived |
| EnOC | Anti-IL6 (stroma) | Paracrine signaling |
| EnOC | Metabolic (OxPhos) | GSEA: OxPhos enriched |

---

*Legenda: * p<0.05, ~ p<0.1, NS = non significativo*
*Analisi completata: 2026-01-31*
