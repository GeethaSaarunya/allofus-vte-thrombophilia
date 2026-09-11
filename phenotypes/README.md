# Phenotype specifications

This directory contains the count-free phenotype concept specifications used for the manuscript repository.

## VTE endpoint

`vte/` contains the Script 05B phenotype reconstruction used for Supplementary Table S4:

- `15_TableS4_VTE_concept_ids_long.csv` — one row per concept.
- `15_TableS4_VTE_concept_ids_grouped.csv` — manuscript-oriented grouped specification.
- `15_TableS4_VTE_concept_ids_grouped.md` — human-readable grouped version.
- `15_concept_export_qc.csv` — reconstruction QC.
- `15_source_discovery.csv` — source-discovery audit.

The reconstruction uses the locked All of Us CDR and the original Script 05B classifier when a frozen concept-disposition file is not mounted.

## Clinical-risk factors

`clinical/` contains the reconstructed concept specification for the final five-factor clinical-risk framework plus pregnancy/postpartum:

- cancer history
- major surgery
- inpatient hospitalization
- infection/sepsis
- fracture/major trauma
- pregnancy/postpartum (reported separately, not included in the five-factor burden)

The refined infection/sepsis specification follows the original Script 05A1 rules.

### Provenance note

The current Workbench reconstruction did **not** recover the historical Script 05A or 05A1 manual-override CSVs. The clinical export therefore reflects the locked CDR plus the original rule-based Script 05A/05A1 classification logic. The source-discovery and QC files are retained so this boundary is explicit.

No participant-level records, person identifiers, participant counts, or record counts are included in these phenotype specification files.
