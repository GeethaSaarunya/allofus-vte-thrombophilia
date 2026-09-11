# Reproducibility notes

## Authoritative public workflow

The uploaded working directory contained many historical and superseded scripts. The public repository intentionally keeps only the final manuscript-facing workflow:

1. `01_primary_five_factor_analysis.py` — derived from the final Script 17 five-factor rebuild.
2. `02_ancestry_sensitivity.py` — derived from the final Script 18 ancestry analysis and supersedes the ancestry placeholders from Script 17.
3. `03_export_vte_concept_specification.py` — internal Workbench reconstruction/export of the count-free VTE phenotype concept specification.
4. `04_reconstruct_clinical_concept_specification.py` — internal Workbench reconstruction/export of the final clinical-risk concept specification using the original 05A/05A1 rules.
5. `05_render_figures.R` — final warm Okabe-Ito Figure 1–3 renderer.

Historical six-factor, seven-factor, reviewer-development, packaging, and intermediate rendering scripts were not copied into the public repository because they are superseded by this final workflow.

## Analysis population locks

The final WGS analysis cohort contained 358,533 participants, including 7,553 with observed VTE during available follow-up.

## Primary clinical construct

The final primary clinical-risk burden contains the following five factors after the pregnancy-timing revision:

- cancer history
- major surgery
- inpatient hospitalization
- infection/sepsis
- fracture/major trauma

Pregnancy/postpartum is evaluated separately because the all-history EHR indicator is not a coherent time-independent burden component.

## Follow-up

In the final WGS analysis cohort, median available post-enrollment EHR follow-up was 1,343 days (IQR 428–1,640 days).

## Important boundaries

- Primary endpoint: observed VTE occurrence during available follow-up.
- The work is an observational cohort association study, not a deployable prediction model.
- The F5/F2 discrimination analysis is secondary/internal.
- No participant-level output is included in this repository.

## Phenotype reconstruction boundary

The clinical reconstruction source-discovery audit did not recover the historical optional Script 05A or 05A1 manual-override files. The deposited clinical concept specification therefore reflects the locked CDR plus the original rule-based 05A/05A1 classification logic. This provenance boundary is retained in the repository rather than being silently treated as a frozen-file recovery.

## Reviewer access and execution

Version 1.0.0 was publicly released for manuscript review; version 1.0.1 corrects documentation and packaging while preserving the same scientific files. Reviewers can inspect all packaged methods, aggregate results, and phenotype documentation without All of Us access. Local figure rendering uses the included aggregate ZIP. The model scripts start from an existing secure participant master; this package does not by itself rebuild that master from raw CDR tables. Controlled Tier authorization and the required secure inputs are necessary to repeat participant-level analyses. See `REVIEWER_GUIDE.md`.
