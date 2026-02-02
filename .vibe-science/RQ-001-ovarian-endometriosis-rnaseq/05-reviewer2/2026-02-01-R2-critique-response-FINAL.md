# Risposta Definitiva alle Critiche R2: compartment-specific-analysis.md

**Data:** 2026-02-01
**Status:** CRITICHE ACCETTATE - DOCUMENTO DA RISCRIVERE

---

## METHODS CAPSULE (Obbligatoria per Ogni Finding)

| Aspetto | Stato | Dettaglio |
|---------|-------|-----------|
| **Dataset** | 52 campioni, 321,232 cellule | 10x Genomics (maggioranza) + Drop-seq (Bertelli) |
| **N donatori per gruppo** | Variabile | EuE: ~10, OMA: ~8, CCOC: **3**, EnOC: ~6, Pelvica: ~5 |
| **Piattaforme** | Eterogenee | 10x 3' v2/v3, Drop-seq |
| **QC/Doublet** | ✅ Applicato | SOLO + Scrublet consensus (soglia 0.95/0.90) |
| **Ambient RNA correction** | ❌ **NON APPLICATO** | SoupX/CellBender non eseguiti |
| **Definizione compartimenti** | Clustering Leiden | Epitelio, Stroma, Immune (res 0.5/1.0) |
| **Integrazione batch** | Harmony su PCA | sample_id come batch variable, 30 PCs |
| **DE pseudobulk** | PyDESeq2 | Aggregazione per sample_id, design ~ batch + condition |
| **Correzione multipla** | Benjamini-Hochberg | FDR < 0.05 per significatività |
| **Soglie** | |log2FC| ≥ 1 + padj < 0.05 | Per ORA; padj < 0.05 per DE |
| **LOOCV / Leave-one-out** | ❌ **NON ESEGUITO** | Sistematico |
| **Compositional analysis** | ❌ **NON ESEGUITO** | scCODA o equivalente |
| **Module scores** | ❌ **NON ESEGUITO** | Per branch UPR |

---

## VERDETTO SULLE CRITICHE

### CRITICA 1: Espressione ≠ Attivazione Pathway
**STATUS: 100% CORRETTA**

**Errore commesso:**
- Ho scritto "PERK activation" basandomi su EIF2AK3 mRNA UP
- PERK è una chinasi attivata post-traduzionalmente (dimerizzazione/fosforilazione)
- mRNA abundance non prova attivazione ISR

**Verifica downstream richiesta:**
La cascata PERK → p-eIF2α → ATF4 (traduzionale) → targets richiede coerenza su:
- ATF4: regolato TRADUZIONALMENTE (uORFs bypass) - mRNA può non riflettere proteina
- Target ISR/ATF4: DDIT3, TRIB3, ATF3, PPP1R15A/GADD34, ASNS, CHAC1

**Nei nostri dati (CCOC vs OMA):**

| Gene | Epi | Stroma | Immune | Coerenza ISR? |
|------|-----|--------|--------|---------------|
| EIF2AK3 | +3.20*** | +2.38*** | +2.02*** | N/A (solo mRNA) |
| DDIT3 | +1.69*** | -0.17 | -1.14*** | **NO** (solo epi) |
| ATF3 | ? | ? | ? | Da verificare |
| TRIB3 | ? | ? | ? | Da verificare |
| ASNS | ? | ? | ? | Da verificare |
| PPP1R15A | ? | ? | ? | Da verificare |

**CONCLUSIONE:** Senza coerenza su gene-set ISR/ATF4, "PERK activation" è **overclaim**.

**CORREZIONE:**
> ~~"PERK activation sistemica"~~
>
> → "EIF2AK3 mRNA elevato in tutti i compartimenti; pattern trascrizionale **compatibile con** stress ER, ma **non dimostrativo** di attività PERK funzionale senza validazione p-eIF2α/ATF4 proteico o coerenza robusta del gene-set ISR."

---

### CRITICA 2: Notazione p vs padj
**STATUS: 100% CORRETTA**

**Errore commesso:**
- Ho usato *** e ~ senza legenda formale
- Non chiaro se p grezzo o padj
- BH applicato su quanti test? Per compartimento? Per gene-set?

**Dalla tesi (3.5):**
> "FDR < 0.05 per il controllo dei falsi positivi"

Quindi i valori sono **padj (BH)**, ma nel file ho scritto "p" per comodità.

**CORREZIONE:**
Ogni tabella deve riportare:
- `padj (BH)` esplicito
- N donatori per gruppo (non N cellule)
- Design formula

**Legenda formale:**
| Simbolo | Significato |
|---------|-------------|
| *** | padj < 0.001 |
| ** | padj < 0.01 |
| * | padj < 0.05 |
| ~ | padj < 0.10 (trend) |
| NS | padj ≥ 0.10 |

---

### CRITICA 3: N donatori mancante
**STATUS: 100% CORRETTA**

**Problema:**
Il file non mostra quanti donatori/samples per gruppo per compartimento.

**Da verificare per ogni contrasto:**

| Contrasto | Compartimento | Gruppo 1 (N) | Gruppo 2 (N) | Cellule min |
|-----------|---------------|--------------|--------------|-------------|
| OMA vs EuE | Epithelium | ? | ? | 50 |
| OMA vs EuE | Stroma | ? | ? | 50 |
| CCOC vs OMA | Epithelium | **3** | ? | 50 |
| CCOC vs OMA | Stroma | **3** | ? | 50 |
| EnOC vs OMA | Epithelium | ? | ? | 50 |

**NOTA CRITICA:** CCOC ha solo **N=3 donatori**. Qualsiasi finding CCOC è underpowered.

---

### CRITICA 4: Harmony non salva DE
**STATUS: 100% CORRETTA**

**Dalla letteratura (Korsunsky et al., 2019):**
> "Harmony modifies cell embeddings... For differential expression, we recommend methods that model batch effects directly in the count space."

**Cosa è stato fatto:**
- Harmony applicato su PCA embedding (corretto per clustering/UMAP)
- DE pseudobulk con `~ batch + condition` quando batch info disponibile

**Problema:**
- Non tutti i campioni hanno batch annotation
- L'eterogeneità piattaforma (10x vs Drop-seq) potrebbe non essere completamente catturata

**CORREZIONE:**
> "L'effetto batch è stato parzialmente modellato nella design formula DE (~ batch + condition) quando l'informazione era disponibile. Tuttavia, l'eterogeneità residua tra piattaforme (10x 3' v2/v3, Drop-seq) rimane un potenziale confondente."

---

### CRITICA 5: IL6/STAT3/SOCS3 narrativa troppo comoda
**STATUS: 100% CORRETTA**

**Errore commesso:**
Ho scritto "SOCS3 DOWN = feedback removal" ma:
- SOCS3 è TARGET di STAT3
- Se STAT3 non è attivo, SOCS3 scende "per inerzia"
- Non significa "freno rimosso"

**Dati EnOC vs OMA (Stroma):**
| Gene | log2FC | padj | Interpretazione |
|------|--------|------|-----------------|
| IL6 | +3.09* | sig | ↑ ligando |
| STAT3 | -1.11* | sig | ↓ trasduttore |
| SOCS3 | -1.99* | sig | ↓ target STAT3 |

**Incoerenza:** Se IL6→JAK→STAT3 fosse attivo, STAT3 targets (incluso SOCS3) dovrebbero essere UP.

**CONCLUSIONE:** Il pattern IL6↑/STAT3↓/SOCS3↓ è **incoerente con attivazione IL6/STAT3** e potrebbe riflettere:
1. Cambiamento composizione cellulare
2. IL6 prodotto ma non segnalante localmente
3. Pathway effettivamente meno attivo

**CORREZIONE:**
> ~~"SOCS3 DOWN → feedback removal → pathway attivato"~~
>
> → "Il pattern IL6↑/STAT3↓/SOCS3↓ nello stroma EnOC è biologicamente **incoerente** con attivazione canonica IL6→JAK→STAT3. Interpretazioni alternative includono shift composizionale, signaling paracrine non-locale, o attenuazione del pathway. Validazione p-STAT3 richiesta."

---

### CRITICA 6: Novità vs cose già note
**STATUS: CORRETTA**

**GIÀ NOTO (non novità):**
- HNF1B in CCOC → firma classica da anni (Tsuchiya et al., 2003)
- UPR/ER stress in tumori ovarici → campo arato
- IL6 nel microambiente ovarico → noto

**POTENZIALMENTE NUOVO (ma non dimostrato):**
- Differenza branch-specific e compartment-specific tra CCOC vs EnOC
- **Requisiti per renderlo nuovo:**
  1. Gene-set scoring per branch (PERK/ATF6/IRE1)
  2. Replicazione su dataset indipendenti
  3. Controllo batch/donatore robusto
  4. Esclusione shift composizionale

---

### CRITICA 7: ATF4 mRNA come "cartina tornasole"
**STATUS: CORRETTA**

**Errore nel mio documento precedente:**
Ho scritto che ATF4 "dovrebbe essere UP" per confermare PERK activation.

**Biologia corretta:**
ATF4 è regolato **traduzionalmente** durante ISR (uORFs bypass). L'mRNA può:
- Non salire
- Persino scendere
- Senza falsificare attivazione a livello proteico

**Cosa posso pretendere dai dati trascrizionali:**
Non "ATF4↑" ma **target program ISR/ATF4 coerente**:
- DDIT3/CHOP ✓ (ma solo in epi)
- TRIB3 - da verificare
- ATF3 - da verificare
- PPP1R15A/GADD34 - da verificare
- ASNS - da verificare
- CHAC1 - da verificare

---

### CRITICA 8: Ambient RNA / SoupX / CellBender
**STATUS: CRITICA DEVASTANTE - NON APPLICATO**

**Verifica nel pipeline (3.2.txt):**
- Doublet removal: ✅ SOLO + Scrublet
- Ambient RNA correction: ❌ **NON MENZIONATO**

**Conseguenza:**
IL6 è gene low-count, zero-inflated, sensibile ad ambient RNA. Senza decontaminazione:
- **IL6 non è interpretabile come "source"**
- Il segnale potrebbe essere contaminazione

**CORREZIONE OBBLIGATORIA:**
> "**Limitazione critica:** Ambient RNA correction (SoupX/CellBender) non è stata applicata. Per geni low-count come IL6, l'interpretazione della 'sorgente cellulare' è fragile e richiede validazione con ISH o sorted populations."

---

### CRITICA 9: Compositional Shift
**STATUS: NON AFFRONTATO**

**Problema:**
"Stroma" include: CAFs, fibroblasti quiescenti, endotelio, periciti, cellule muscolari lisce.
IL6↑ nello stroma potrebbe essere:
- Biologia (attivazione)
- Composizione (più CAFs nel tumore)

**Analisi richiesta:**
- Differential abundance (scCODA o equivalente)
- Controllo per proporzioni cellulari

**STATUS:** Non eseguito.

---

### CRITICA 10: TCGA-OV non è adatto per CCOC/EnOC
**STATUS: CORRETTA**

**Errore nel mio documento:**
Ho scritto "replicazione esterna (TCGA-OV split per istotipo)"

**Problema:**
TCGA-OV è centrato su HGSOC (high-grade serous), non su CCOC/EnOC.

**CORREZIONE:**
> ~~"TCGA-OV split per istotipo"~~
>
> → "Replicazione su coorti multi-istotipo specifiche per CCOC/EnOC (es. AOCS, Japanese Gynecologic Oncology Group) o validazione proteomica/fosfoproteomica (CPTAC) quando disponibile."

---

### CRITICA 11: ISRIB / Terapeutica senza citazioni
**STATUS: CORRETTA**

**Errori corretti:**
1. ISRIB agisce su **eIF2B**, non eIF2α (Sidrauski et al., 2015, eLife)
2. PERK inhibitors hanno tossicità pancreatica (Atkins et al., 2013, PNAS)

**CORREZIONE:**
Rimuovere completamente la sezione "Implicazioni Terapeutiche" o riscriverla con citazioni e caveat su tossicità.

---

## CONTROVERIFICHE RICHIESTE PRIMA DI QUALSIASI FINDING

| # | Verifica | Status | Azione |
|---|----------|--------|--------|
| 1 | Tabella N (donatori per compartimento × condizione) | ❌ Mancante | Generare |
| 2 | DE pseudobulk con design formula esplicita | ✅ Fatto | Documentare |
| 3 | Leave-one-donor-out per geni cardine | ❌ Non fatto | Eseguire per IL6, EIF2AK3, DDIT3, XBP1, HSPA5 |
| 4 | Module score per branch UPR | ❌ Non fatto | Calcolare PERK/IRE1/ATF6 scores |
| 5 | Compositional analysis | ❌ Non fatto | Eseguire scCODA |
| 6 | Ambient RNA check per IL6 | ❌ Non fatto | SoupX post-hoc o dichiarare limitazione |

---

## DOCUMENTO RIVISTO: COSA SI PUÒ AFFERMARE

### Affermazioni Valide (Descrittive)

1. "EIF2AK3 mRNA è differenzialmente espresso tra CCOC ed EnOC vs OMA, con pattern compartmentale diverso."

2. "IL6 mRNA mostra pattern compartmentale opposto tra CCOC (epitelio) ed EnOC (stroma), ma l'interpretazione funzionale è limitata dall'assenza di decontaminazione ambient RNA."

3. "I target a valle del pathway IL6/STAT3 (incluso SOCS3) mostrano pattern incoerente con attivazione canonica."

### Affermazioni NON Valide (Da Rimuovere)

1. ~~"PERK activation sistemica"~~
2. ~~"IL6 source cell-specific"~~
3. ~~"UPR globalmente soppresso"~~
4. ~~"SOCS3 DOWN = feedback rimosso"~~
5. ~~"Implicazioni terapeutiche" specifiche~~
6. ~~Qualsiasi claim meccanicistico senza coerenza ligand→TF→targets~~

---

## NEXT STEPS OBBLIGATORI

1. **Generare tabella N donatori** per ogni compartimento × condizione
2. **Eseguire LOOCV** per IL6, EIF2AK3, DDIT3, XBP1, HSPA5
3. **Calcolare module scores** per branch UPR (PERK/ATF6/IRE1)
4. **Verificare target ISR/ATF4** (TRIB3, ATF3, PPP1R15A, ASNS, CHAC1)
5. **Dichiarare limitazione ambient RNA** esplicitamente
6. **Riscrivere documento** in forma descrittiva, non meccanicistica

---

*Risposta completata: 2026-02-01*
*Tutte le critiche R2 accettate*
