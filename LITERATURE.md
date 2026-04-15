# Literature Review: PROTAC & Macrocycle Chameleonicity

Compiled 2026-04-13. Papers with experimental chameleonicity data for PROTACs,
macrocycles, and bRo5 compounds. Priority-ordered by relevance to expanding
our benchmark.

## Chamelogk Values for Our Benchmark PROTACs

From Garcia Jimenez et al. 2023 (DOI: 10.1021/acs.jmedchem.3c00823).
Threshold: Chamelogk >= 0.6 = chameleonic.

| Benchmark PROTAC | ChameLogK | Label |
|-----------------|-----------|-------|
| MZ1 | 1.15 | CHAMELEONIC |
| dBET1 | 0.80 | CHAMELEONIC |
| dBET6 | 0.86 | CHAMELEONIC |
| ARV-825 | 0.72 | CHAMELEONIC |
| ARV-110 | 0.99 | CHAMELEONIC |
| ARV-771 | — | NOT MEASURED |

Additional PROTACs from the same paper NOT in our benchmark:
| PROTAC | ChameLogK | Label |
|--------|-----------|-------|
| BRD4 degrader AT1 | 1.26 | chameleonic |
| cisMZ1 | 1.27 | chameleonic |
| gefitinib PROTAC 3 | 1.36 | chameleonic |
| MZP-54 | 1.18 | chameleonic |
| PROTAC-1 | 1.07 | chameleonic |
| CRBN-6-5-5-VHL | 1.05 | chameleonic |
| BI-0319 | 0.99 | chameleonic |
| PROTAC Mcl degrader-1 | 0.93 | chameleonic |
| PROTAC BET degrader-10 | 0.83 | chameleonic |
| PROTAC FAK degrader-1 | 0.79 | chameleonic |
| dBET57 | 0.68 | chameleonic |
| ZXH-3-26 | 0.65 | chameleonic |
| BI-4206 | 0.70 | chameleonic |
| **BI-3663** | **0.50** | **NON-CHAMELEONIC** |

From Caron review (ADMET & DMPK 2024, DOI: 10.5599/admet.2334):
- **ARV-110:** Chamelogk = 0.99 (chameleonic, oral)
- **ARV-471:** Chamelogk = 1.35 (chameleonic, oral, not in benchmark)

### DISCREPANCY: Existing Benchmark Labels vs Chamelogk
| Molecule | Our Label | ChameLogK | ChameLogK Label |
|----------|-----------|-----------|-----------------|
| Everolimus | chameleon | 0.45 | non-chameleonic |
| Sirolimus | chameleon | 0.25 | non-chameleonic |
| Cyclosporine A | chameleon | 1.25 | chameleonic (agrees) |

This suggests our benchmark labels may use a more permissive definition
(oral bioavailability of large molecules) than the Chamelogk threshold.
Worth investigating: are these labels from a paper, or were they assigned
based on oral availability alone?

## Anti-Chameleonic Hydrophobic Collapse (Critical for protac_3)

**Poongavanam et al. 2025** (DOI: 10.1021/acsmedchemlett.5c00068)
- Title: Linker-Determined Folding and Hydrophobic Interactions Explain a
  Major Difference in PROTAC Cell Permeability
- Key finding: Alkyl-linker PROTACs show REVERSED chameleonicity —
  hydrophobic collapse in WATER forces folded conformations, while in
  CHCl3 they remain extended. PEG PROTACs adopt similar conformations
  in both solvents.
- **This is exactly what protac_3 (alkyl linker) does.** It's a published,
  named phenomenon — not just our pipeline's failure to detect non-chameleonicity.
- Implication: A descriptor that detects hydrophobic collapse tendency
  (e.g., polar fraction of linker atoms, or ratio of linker SASA in
  water vs CHCl3) could distinguish protac_3.

## Priority 1: Papers with PROTAC Chameleonicity Labels + SMILES

### Garcia Jimenez et al. 2023 — Chamelogk (55 compounds, 18 PROTACs)
- **Journal:** J. Med. Chem. 66, 2023
- **DOI:** 10.1021/acs.jmedchem.3c00823
- **PMC:** PMC10424176 (open access)
- **Data:** 55 molecules (25 Ro5 + 30 bRo5 including 18 PROTACs) with Chamelogk scores
- **STATUS:** NEED SI for SMILES of all 55 compounds

### AstraZeneca 2026 — 68 PROTACs with Chameleonicity Data
- **Title:** Evaluation of Oral PROTAC Guidelines: Efflux Ratio Outweighs
  Chameleonicity Descriptors
- **Journal:** ACS Med. Chem. Lett., 2026
- **DOI:** 10.1021/acsmedchemlett.6c00043
- **Data:** 68 PROTACs with chamelogk curves + physicochemical properties
- **Key finding:** Chameleonicity alone does NOT predict oral bioavailability;
  efflux ratio more important
- **STATUS:** NEED full paper for SMILES + chameleonicity data

### Garcia Jimenez et al. 2024 — cChameCS Index (6 PROTACs)
- **Journal:** J. Med. Chem. 67(13), 2024
- **DOI:** 10.1021/acs.jmedchem.4c01200
- **Figshare:** https://figshare.com/collections/7308421
- **Data:** 6 CRBN-based BET PROTACs with chameleonicity labels + IMHB patterns
- **Defines cChameCS** — IMHB-based chameleonicity predictor; closest to our CI
- **STATUS:** NEED Figshare data + paper for cChameCS formula definition

### Garcia Jimenez et al. 2025 — Linker Methylation (11 PROTACs)
- **Journal:** J. Med. Chem. 68(15), 2025
- **DOI:** 10.1021/acs.jmedchem.5c01497
- **Data:** 11 VHL PROTACs with linker methylation variants + chameleonicity data
- **Key finding:** Efflux ratio (not Caco-2 permeability) predicts oral bioavailability
- **STATUS:** NEED SI for SMILES + data

### Merck 2025 — Metadynamics (32 degraders)
- **Journal:** J. Chem. Inf. Model. 2025
- **DOI:** 10.1021/acs.jcim.5c01600
- **Data:** 32 degraders with experimental permeability + WT-MetaD (Rg, 3D-PSA, IMHB)
- **STATUS:** NEED SI for SMILES + data

### Poongavanam et al. 2025 — Anti-Chameleonic Hydrophobic Collapse (2 PROTACs)
- **Journal:** ACS Med. Chem. Lett. 16(4), 2025
- **DOI:** 10.1021/acsmedchemlett.5c00068
- **Data:** 2 VHL PROTACs (PEG vs alkyl linker) with NMR + MD
- **STATUS:** NEED paper for SMILES

### dBET1 Diastereomers (2026)
- **Journal:** J. Med. Chem. 2026
- **DOI:** 10.1021/acs.jmedchem.5c02791
- **Data:** dBET1 analogs with stereochemistry-dependent chameleonicity

## Priority 2: Papers with PROTAC Conformational Data (Can Derive Labels)

### Poongavanam et al. 2024 — Linker Composition (9 VHL PROTACs)
- **Journal:** J. Med. Chem. 2024
- **DOI:** 10.1021/acs.jmedchem.4c02492
- **Figshare:** https://figshare.com/collections/7590671
- **Data:** 9 VHL PROTACs with PAMPA + NMR-derived conformational ensembles
- **Key data:** Compound 2 (PEG, Rg 5.8 A, chameleonic), Compound 1 (aliphatic,
  PAMPA logPe = -8.29, non-chameleonic). ~1000-fold permeability difference
  from O-for-CH2 substitutions alone.
- **STATUS:** NEED Figshare data

### Poongavanam et al. 2022 — Linker-Dependent Folding (3 PROTACs)
- **Journal:** J. Med. Chem. 65, 2022
- **DOI:** 10.1021/acs.jmedchem.2c00877 (PMC9574858, open access)
- **Data:** 3 CRBN BRD4 PROTACs with different linkers
- **Key:** Compound 1 (PEG, 3D PSA 209 A^2, Caco-2 30 nm/s) = chameleonic;
  Compound 3 (alkyl, 3D PSA 290-330 A^2, Caco-2 6 nm/s) = non-chameleonic
- **NOTE:** These may be protac_1/2/3 from our user_protacs.tsv
- **STATUS:** Open access, NEED to verify SMILES match

### Atilaw et al. 2020 — VHL PROTAC (1 PROTAC)
- **Journal:** ACS Med. Chem. Lett. 12, 2020
- **DOI:** 10.1021/acsmedchemlett.0c00556 (PMC7812666, open access)
- **Data:** 1 VHL PROTAC (MW 1034), NMR-validated chameleonic

### Ermondi et al. 2023 — Computational Validation of PROTAC-1
- **Journal:** Pharmaceutics 15(1), 2023
- **DOI:** 10.3390/pharmaceutics15010272 (open access)
- **Data:** 1 VHL PROTAC with extensive conformational analysis
- **Key:** QM (PM7) with implicit solvation worsened predictions — confirms our
  finding that implicit solvation collapses conformers

### FerroTACs (JACS 2025)
- **DOI:** 10.1021/jacs.4c18354
- **Data:** Ferrocene-linker PROTACs with freely-rotating Cp rings as molecular
  hinges enabling chameleonicity

## Priority 3: Larger Datasets (Permeability as Proxy)

### AbbVie 2025 — ETR ML Model (Figshare)
- **DOI:** 10.1021/acs.jmedchem.5c00536
- **Figshare:** https://figshare.com/collections/7820637
- **Data:** Thousands of compounds with ETR = EPSA/TPSA

### AbbVie 2024 — Polarity Reducers Review
- **DOI:** 10.1021/acs.jmedchem.3c02332
- **Data:** ~10,000 AbbVie compounds with EPSA (mostly proprietary)

### Kihlberg 2023 — PROTAC Permeability ML (~700 PROTACs)
- **DOI:** 10.1021/acsomega.2c07717 (PMC9933238, open access)

### AstraZeneca 2024 — eHBD Rule
- **DOI:** 10.1021/acs.jmedchem.4c01017
- **Key:** Exposed HBD <= 2 in apolar environment is critical for oral PROTACs

### Arvinas 2023 — 1806 PROTACs (proprietary structures)
- **DOI:** 10.1021/acs.jmedchem.3c00740
- **Key thresholds:** unsatisfied HBD <= 2, MW <= 950, TPSA <= 200, RotBonds <= 14

### Scott et al. 2020 — AR PROTAC Permeability (~20 PROTACs)
- **DOI:** 10.1021/acsmedchemlett.0c00194

### SweMacrocycleDB — 4,216 macrocycles
- **URL:** https://swemacrocycledb.com/
- **DOI:** 10.1038/s41597-024-04302-z

### Tang et al. 2023 — aMD of 47 Macrocycles
- **DOI:** 10.1021/acs.jcim.3c01123
- **Data:** 47 peptidic macrocycles with PAMPA + conformational data

### Delta PSA Preprint (2026)
- **URL:** https://www.biorxiv.org/content/10.64898/2026.01.06.697862v1
- **Data:** Macrocyclic peptide permeability benchmark

## Key Methodological Insights from Literature

1. **Implicit-solvent methods systematically fail.** GBSA, GFN-FF, GFN2-xTB all
   collapse conformers to compact states. MD with explicit solvent required for
   accuracy. Our pipeline works AROUND this using ETKDG + MMFF ensemble statistics.

2. **IMHBs are the primary chameleonicity driver** (Caron 2024). The cChameCS index
   is built on IMHB patterns. Our CI's term_imhb is aligned with this.

3. **PEG vs alkyl linkers are the key differentiator.** PEG has gauche preference
   (O-C-C-O dihedral) enabling folding. Alkyl has anti preference preventing it.
   This is well-established physical chemistry.

4. **Anti-chameleonic hydrophobic collapse** (Poongavanam 2025): alkyl-linker
   PROTACs fold in WATER (wrong solvent!) via hydrophobic effect, but extend in
   CHCl3. This means our ETKDG+MMFF ensemble (solvent-agnostic) might be
   capturing the water-folded state for protac_3 rather than the CHCl3-extended
   state, inflating its CI.

5. **Efflux ratio matters more than chameleonicity for oral bioavailability**
   (AstraZeneca 2026, Garcia Jimenez 2025). Chameleonicity enables permeation
   but doesn't prevent efflux. This is a caveat for our pipeline's clinical
   relevance.

6. **ETR = EPSA/TPSA** is a simple, published chameleonicity metric (AbbVie).
   Could be computed from our conformational ensembles for comparison.

## Key References

- Whitty et al. 2016 (DOI: 10.1016/j.drudis.2016.02.005) — foundational framework
- Poongavanam & Kihlberg 2024 (DOI: 10.1038/s41570-023-00563-1) — comprehensive review
- David et al. 2021 (DOI: 10.1002/cmdc.202100306) — chameleonic efficiency
