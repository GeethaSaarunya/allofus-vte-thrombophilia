# All of Us VTE clinical-genetic burden reproducibility materials

Reproducibility materials for:

**Clinical Risk Burden and Inherited Thrombophilia in Venous Thromboembolism: An EHR-Linked Genomic Cohort Study in the All of Us Research Program**

Authors: Geetha Saarunya, Brian Steffen, Preeti Magesh, Christopher J. Tignanelli, and Todd Costantini.

## Scope

This repository contains the final public analysis code, aggregate-safe source tables, and figure-rendering code for the manuscript. It does **not** contain participant-level All of Us data and cannot be used to reconstruct or redistribute the All of Us Controlled Tier dataset.

The final manuscript analysis used the All of Us Controlled Tier Curated Data Repository v8 (CDR C2024Q3R9). The public code expects authorized users to supply the secure participant-level master inside the All of Us Researcher Workbench.

## Repository structure

- `scripts/01_primary_five_factor_analysis.py` — final five-factor primary manuscript analyses.
- `scripts/02_ancestry_sensitivity.py` — final ancestry-stratified F5/F2 sensitivity analysis.
- `scripts/03_export_vte_concept_specification.py` — internal Workbench utility that reconstructs/exports the count-free VTE concept specification from the locked CDR when the original frozen 05B file is not mounted.
- `scripts/04_reconstruct_clinical_concept_specification.py` — internal Workbench utility that reconstructs the final clinical-risk concept specification from the locked CDR using the original Script 05A and 05A1 rules.
- `scripts/05_render_figures.R` — publication figure rendering from aggregate-safe Script 17 results.
- `data/aggregate/` — aggregate-safe manuscript source tables.
- `figures/` — final Figures 1–3 in PNG/PDF/SVG plus captions.
- `docs/` — reproducibility, data-use, and release notes.
- `environment/` — package/session information available from the final rendering workflow.

## Secure inputs

The Python analysis scripts require a participant-level master that remains inside the authorized All of Us environment. Set:

```bash
export VTE_JTH_MASTER_FILE=/secure/path/to/05_vte_master_analysis_internal.parquet
export WORKSPACE_CDR=<authorized_cdr_dataset>
```

Then run:

```bash
python scripts/01_primary_five_factor_analysis.py
python scripts/02_ancestry_sensitivity.py
```

The scripts write only aggregate outputs.

## Render figures locally

From the repository root:

```bash
Rscript scripts/05_render_figures.R
```

The renderer uses `data/aggregate/VTE_Results_17_20260906.zip` by default. You may override paths with:

```bash
export VTE_REPO_ROOT=/path/to/repository
export VTE_RESULTS17_ZIP=/path/to/VTE_Results_17_20260906.zip
```

## Data access

Participant-level All of Us EHR and genomic data are not redistributed. Authorized researchers can access the Controlled Tier through the All of Us Researcher Workbench subject to applicable program and institutional requirements.

Public aggregate tables in this repository were selected from outputs used in the manuscript and were checked for the All of Us minimum-cell dissemination rule. Users remain responsible for complying with the current All of Us Data and Statistics Dissemination Policy when generating new outputs.

## Final ancestry results

The final ancestry-stratified sensitivity results are in `data/aggregate/ancestry_sensitivity/` and supersede the empty ancestry placeholders created by the earlier primary rebuild.

## Phenotype specifications

The manuscript reports exact VTE and clinical phenotype definitions in the Supplement. `scripts/03_export_vte_concept_specification.py` is included so the frozen VTE concept-disposition files can be exported inside the authorized Workbench before the final public release.

**The count-free VTE and clinical-risk phenotype specifications generated in the authorized Workbench are included under `phenotypes/`. Their QC and source-discovery audits are included alongside them.**

## Citation

See `CITATION.cff`. A Zenodo DOI and journal DOI should be added after they are assigned.

## Funding

Supported by NIH/NIGMS 5R35GM149345 and NIH/NCATS UM1TR004405.

## Status

This repository package is prepared for manuscript submission. Recommended public release tag: `v1.0.0` at acceptance/publication.

## Phenotype reconstruction provenance

The historical local V2 project directory was no longer mounted when the phenotype specifications were prepared for repository deposit. Scripts 03 and 04 therefore reconstruct the vocabulary-level concept specifications from the locked CDR using the original curation logic. The clinical source-discovery audit did not recover the historical optional 05A/05A1 manual-override CSVs; this is documented explicitly in `phenotypes/clinical/04_source_discovery.csv`.
