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

## Phenotype specification provenance

The deposited VTE and clinical-risk concept specifications were reconstructed from All of Us CDR v8 (C2024Q3R9) using the study’s rule-based curation procedures. The repository includes reconstruction scripts, concept-disposition tables, and accompanying source-discovery and quality-control records. Additional methodological details are provided in the [reproducibility notes](docs/REPRODUCIBILITY.md).

**The count-free VTE and clinical-risk phenotype specifications generated in the authorized Workbench are included under `phenotypes/`. Their QC and source-discovery audits are included alongside them.**

## Citation

See `CITATION.cff` for citation metadata. The published **v1.0.0** archive is available at [10.5281/zenodo.22715380](https://doi.org/10.5281/zenodo.22715380); the [all-versions record](https://doi.org/10.5281/zenodo.22715379) links the release series. This working tree prepares **v1.0.1**, a documentation and packaging correction. Cite the version-specific DOI of the release actually used. The v1.0.1 DOI will be assigned when that release is archived.

## Funding

Supported by NIH/NIGMS 5R35GM149345 and NIH/NCATS UM1TR004405.

## Status

**Version 1.0.1 — documentation and packaging correction.**

The v1.0.0 release is already public. This update restores two omitted repository support files, removes the duplicate README provenance section, and refreshes release metadata and checksums. Scientific scripts, aggregate results, figures, and phenotype specifications are unchanged. Publish this correction under **`v1.0.1`**, preserving the existing v1.0.0 tag and archive. The software release does not imply journal acceptance.

See [Reviewer guide](docs/REVIEWER_GUIDE.md), [release notes](docs/RELEASE_NOTES.md), and [release steps](docs/RELEASE_CHECKLIST.md). Any later substantive revision should receive a new release version so the reviewed snapshot remains identifiable.
