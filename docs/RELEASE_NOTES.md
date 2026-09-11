# v1.0.0 — Manuscript submission and peer-review materials

Supporting code, aggregate source results, publication figures, and phenotype specifications for:

**Clinical Risk Burden and Inherited Thrombophilia in Venous Thromboembolism: An EHR-Linked Genomic Cohort Study in the All of Us Research Program**

This version provides a fixed snapshot for manuscript submission and peer review. It does not indicate journal acceptance.

## Included

- Final five-factor primary clinical-risk analysis, based on Script 17.
- Final F5/F2 ancestry sensitivity analysis, based on Script 18.
- Aggregate sources for the main tables, Figure 3, interaction tests, discrimination/calibration, pregnancy timing, and supporting analyses.
- Final Figures 1–3 in PNG, PDF, and SVG, with captions and the v2.3 figure renderer.
- Count-free VTE and clinical concept specifications, reconstruction scripts, and provenance/QC records.
- Citation metadata, environment information, checksums, release validation summary, and reviewer guide.

## Study snapshot

All of Us Controlled Tier CDR v8 (C2024Q3R9); WGS cohort 358,533; observed VTE 7,553. The primary burden includes cancer history, major surgery, inpatient hospitalization, infection/sepsis, and fracture/major trauma. Pregnancy/postpartum is evaluated separately. The primary endpoint is observed VTE occurrence during available follow-up.

## Reproducibility scope

The public aggregate inputs support inspection of reported results and local figure rendering. Re-fitting participant-level models requires authorized Controlled Tier access and the secure analysis master described in the README. Participant-level All of Us data are not distributed.

The clinical concept specification was reconstructed from the locked CDR and original 05A/05A1 rules. Historical optional manual-override CSVs were not recovered; this release does not assert verified identity to those missing historical files.

## Packaging correction

This package revision aligns the citation and documentation with a submission-stage v1.0.0 release and adds review guidance and validation records. Existing scientific scripts, numerical source outputs, and figure bytes are unchanged.
