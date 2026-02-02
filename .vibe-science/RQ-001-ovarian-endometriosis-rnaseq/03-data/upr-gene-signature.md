# UPR Gene Signature

Gene lists for Unfolded Protein Response analysis.

## Source 1: MSigDB Hallmark UPR (113 genes)

**Source:** https://www.gsea-msigdb.org/gsea/msigdb/human/geneset/HALLMARK_UNFOLDED_PROTEIN_RESPONSE.html

**Description:** Genes up-regulated during unfolded protein response, a cellular stress response related to the endoplasmic reticulum.

### Gene List

```
ALDH18A1, ARFGAP1, ASNS, ATF3, ATF4, ATF6, ATP6V0D1, BAG3, BANF1, CALR,
CCL2, CEBPB, CEBPG, CHAC1, CKS1B, CNOT2, CNOT4, CNOT6, CXXC1, DCP1A,
DCP2, DCTN1, DDIT4, DDX10, DKC1, DNAJA4, DNAJB9, DNAJC3, EDC4, EDEM1,
EEF2, EIF2AK3, EIF2S1, EIF4A1, EIF4A2, EIF4A3, EIF4E, EIF4EBP1, EIF4G1,
ERN1, ERO1A, EXOC2, EXOSC1, EXOSC10, EXOSC2, EXOSC4, EXOSC5, EXOSC9,
FKBP14, FUS, GEMIN4, GOSR2, H2AX, HERPUD1, HSP90B1, HSPA5, HSPA9, HYOU1,
IARS1, IFIT1, IGFBP1, IMP3, KDELR3, KHSRP, KIF5B, LSM1, LSM4, MTHFD2,
NFYA, NFYB, NHP2, NOLC1, NOP14, NOP56, NPM1, NABP1, PAIP1, PARN, PDIA5,
PDIA6, POP4, PREB, PSAT1, RPS14, RRP9, SDAD1, SEC11A, SEC31A, SERP1,
SHC1, MTREX, SLC1A4, SLC30A5, SLC7A5, SPCS1, SPCS3, SRPRA, SRPRB, SSR1,
STC2, TARS1, TATDN2, TSPYL2, SKIC3, TUBB2A, VEGFA, WFS1, WIPI1, XBP1,
XPOT, YIF1A, YWHAZ, ZBTB17
```

## Source 2: GO BP Response to ER Stress (283 genes)

**Source:** GO:0034976 via MSigDB

**Description:** Any process that results in a change in state or activity of a cell as a result of ER stress.

### Key UPR Pathway Genes (Core)

**Master regulators:**
- HSPA5 (GRP78/BiP) - chaperone, UPR sensor
- ERN1 (IRE1α) - kinase/RNase, XBP1 splicing
- EIF2AK3 (PERK) - kinase, translation attenuation
- ATF6 - transcription factor

**Downstream effectors:**
- XBP1 - transcription factor (spliced form active)
- ATF4 - transcription factor (PERK pathway)
- DDIT3 (CHOP) - pro-apoptotic transcription factor
- CHAC1 - glutathione degradation
- PPP1R15A (GADD34) - feedback dephosphorylation

**ER chaperones:**
- CALR (calreticulin)
- CANX (calnexin)
- PDIA3, PDIA4, PDIA6 - protein disulfide isomerases
- HSP90B1 (GRP94)
- HYOU1 (GRP170)

**ERAD components:**
- EDEM1, EDEM2, EDEM3 - mannose trimming
- SEL1L - adaptor
- SYVN1 (HRD1) - E3 ligase
- DERL1, DERL2, DERL3 - dislocation
- VCP (p97) - AAA+ ATPase

**Apoptosis mediators:**
- BAX, BAK1 - mitochondrial permeabilization
- BBC3 (PUMA), PMAIP1 (NOXA) - BH3-only
- TRIB3 - ATF4 target, AKT inhibitor

## Analysis Strategy

### For GSEA/Pathway Enrichment
Use **Hallmark UPR** (113 genes) - curated, balanced size

### For Targeted Analysis
Use **Core UPR genes** (20-30 genes):
```
HSPA5, ERN1, EIF2AK3, ATF6, XBP1, ATF4, DDIT3, CHAC1, PPP1R15A,
CALR, CANX, PDIA3, PDIA6, HSP90B1, HYOU1, EDEM1, SEL1L, SYVN1,
DERL1, VCP, BAX, BBC3, TRIB3, HERPUD1, DNAJB9, DNAJC3
```

### For Comprehensive Analysis
Use **GO BP ER stress** (283 genes) - complete coverage
