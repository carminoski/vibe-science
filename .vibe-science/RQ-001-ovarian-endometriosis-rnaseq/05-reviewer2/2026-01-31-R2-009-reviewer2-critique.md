# Reviewer #2 Critique (Unfiltered)

**Date:** 2026-01-31
**Finding:** FINDING-009 (Divergent Mechanisms CCOC vs EnOC)
**Verdict:** MAJOR REVISION

---

## Valutazione complessiva

L'idea "**meccanismi divergenti**: **UPR/ER-stress** per Endometrioid ovarian carcinoma vs **asse Interleukin-6–JAK/STAT3** per Ovarian clear cell carcinoma" è *intrigante*, ma allo stato attuale è **sovra-interpretata** rispetto alla solidità dei dati e contiene **incongruenze tecniche** che un revisore serio ti farà a pezzi.

---

## Major issues (da sistemare prima di parlare di "meccanismo")

### 1) Incoerenza tecnica sulla piattaforma: errore consolidato

Scrivi che Gene Expression Omnibus GSE157153 sarebbe "Illumina RNA-seq / GPL17303".

Ma la pagina ufficiale GEO indica chiaramente **GPL17303 = Ion Torrent Proton**. ([NCBI](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE157153))

Questo non è un dettaglio: cambia bias, copertura, quantificazione e comparabilità. **Correggilo ovunque** e metti una nota metodologica (e.g., possibili bias di misurazione su geni low-count come IL6).

### 2) XBP1 splicing: salto logico non giustificato dai tuoi dati

Nel testo proponi una catena causale (IRE1/XBP1 → IL6 → JAK/STAT3).
Problema: **dai dati che descrivi stai misurando espressione**, non **attività** della via:

* XBP1 "attivo" è **XBP1s (spliced)**, che richiede evidenza di splicing/junction reads o saggi dedicati.
* UPR è *branch-specific*: PERK/ATF4, ATF6, IRE1/XBP1 possono muoversi in modo discordante (e in Ciavattini infatti vedi ATF6/GRP78 ↑ ma CHOP/XBP1 ↓ nell'endometrioide). ([PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC6144292/))

**Quindi**: la formulazione corretta oggi è "*pattern compatibile con branch UPR rimodellati*", non "UPR attivata" (né tantomeno "XBP1 splicing guida IL6").

### 3) N troppo piccolo e rischio "storytelling" (specie su IL6)

Indichi che lo scRNA è limitato (es. CCOC n≈3) e lo usi per concludere "meccanismi divergenti".

Con n così basso:

* FDR e ranking TOP50 diventano **instabili**.
* Una singola paziente "outlier" può inventarti IL6 "+10.98 log2FC".

Serve *almeno* una validazione indipendente (bulk o scRNA) prima di fissare la narrativa.

### 4) IL6 "epiteliale": altissimo rischio di contaminazione/mis-annotation

Dici che IL6 è TOP50 **nel comparto epiteliale** CCOC con effetto enorme.
Questa è esattamente la cosa che un revisore sospetta subito:

* IL6 spesso viene da **stromali/immune**; se appare "epiteliale" potrebbe essere **doublet**, **ambient RNA**, o annotazione troppo permissiva.
* Devi mostrare **percentuale di cellule IL6+** per cluster, **UMI distribution**, e controlli anti-doublet/ambient.

Perché, paradossalmente, in letteratura l'asse IL6-STAT3 è sì fortissimo in OCCC, ma non implica automaticamente che IL6 sia *tumor-intrinsic* nel tuo dataset. (In OCCC l'asse IL6-STAT3-HIF è noto e associato a outcome). ([AACR Journals](https://aacrjournals.org/clincancerres/article/17/8/2538/12035/IL6-STAT3-HIF-Signaling-and-Therapeutic-Response))

### 5) La "lacuna" bibliografica è più fragile di come la vendi

Tu sostieni che UPR+EAOC è praticamente ignorato ("2–3 paper").
Sì, l'intersezione *stretta* può essere sottocoperta, ma nel 2025 esistono review su UPR nei tumori ginecologici/ovarici. ([Tandfonline](https://www.tandfonline.com/doi/full/10.1080/15384101.2025.2543091))

Se vuoi che passi, devi rifinire la claim in modo chirurgico:

* non "nessuno parla di UPR", ma "manca un modello comparativo **EnOC vs CCOC** *nel contesto endometriosi-associato* che colleghi UPR/ER-stress a traiettorie divergenti".

---

## Cosa pretendo come controverifica (minimo sindacale)

1. **Correzione piattaforma + sensitivity analysis**

   * aggiorna pipeline assumendo Ion Torrent.
   * rifai DE con filtri su low-expression e shrinkage; verifica se IL6 resta "mostruoso".

2. **Prova/controprova UPR branch-specific**

   * score separati PERK/ATF4, ATF6, IRE1/XBP1 (Reactome/GO curati, non solo Hallmark).
   * se possibile, stima proxy di XBP1s (bulk con junction reads; nello scRNA spesso non ci arrivi).
   * allinea interpretazione a ciò che mostra Ciavattini: ATF6/GRP78 ↑ con CHOP/XBP1 ↓ non è "UPR ON" in senso banale.

3. **Attribuzione cellulare di IL6**

   * per scRNA: (a) IL6 per cluster, (b) doublet score, (c) ambient RNA correction (SoupX/CellBender), (d) conferma con marker epiteliali stringenti.
   * dimostra che STAT3 target genes (SOCS3, IL6R, etc.) sono coerenti nel *medesimo* comparto.

4. **Validazione esterna**

   * trova almeno 1 dataset indipendente (bulk o scRNA) con CCOC/EnOC ± endometriosi adiacente e ripeti: "UPR-skew in EnOC" e "IL6/STAT3-skew in CCOC".
   * senza questo, la tua Finding-009 resta un *bel pattern*, non un risultato difendibile.

5. **Alternative explanation audit (obbligatorio)**

   * confondenti clinici: stadio, grado, trattamento, purezza tumorale.
   * confondenti tecnici: library size, percent mito, stress dissociazione, batch.
   * confondenti biologici: ipossia dell'endometrioma, emorragia, infiltrato macrofagico.

---

## Minor (ma ti possono bocciare lo stesso)

* Definisci in modo univoco "EAOC": include **solo** EnOC+CCOC? includi anche altri? Nel testo è sottinteso ma non blindato.
* Evita frasi tipo "UPR absent in CCOC": troppo assolute. Scrivi "non emerge in modo robusto nei nostri score attuali".
* Quando citi "asse IRE1/XBP1/IL6/STAT3", ancora il link ER-stress→IL6/STAT3 (non solo HCC): esiste letteratura di overview su ER-stress in cancro che include esplicitamente quel wiring. ([Nature](https://www.nature.com/articles/s41420-024-02110-3))

---

## Verdict

**Major revision**. La direzione può diventare Q1-worthy, ma solo se smetti di "chiudere il cerchio" a parole e lo chiudi con: (i) correzione piattaforma, (ii) branch-specific UPR, (iii) attribuzione cellulare IL6 inattaccabile, (iv) almeno una validazione esterna.

---

## References dal Reviewer

1. [GEO GSE157153 - Platform GPL17303 = Ion Torrent Proton](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE157153)
2. [Ciavattini 2018 - UPR in endometrioid](https://pmc.ncbi.nlm.nih.gov/articles/PMC6144292/)
3. [IL6-STAT3-HIF in OCCC](https://aacrjournals.org/clincancerres/article/17/8/2538/12035/IL6-STAT3-HIF-Signaling-and-Therapeutic-Response)
4. [UPR in gynecological tumors - 2025 review](https://www.tandfonline.com/doi/full/10.1080/15384101.2025.2543091)
5. [ER stress in cancer - Nature 2024](https://www.nature.com/articles/s41420-024-02110-3)

---

*Received: 2026-01-31*
*Status: MAJOR REVISION REQUIRED*
