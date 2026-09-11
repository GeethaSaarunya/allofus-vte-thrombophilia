# Reviewer guide

This submission-stage release contains analysis code, aggregate source outputs, figures, and phenotype documentation. The manuscript and formatted supplementary tables are supplied through the journal submission system; this guide maps their analyses to the underlying public files.

## Where to find the evidence

All paths below are relative to the repository root.

| Manuscript component | Public source |
|---|---|
| Main Table 1: cohort and five-factor burden | `data/aggregate/primary_analysis/17B_Main_Table1_revised_five_factor.csv` |
| Main Table 2: adjusted clinical associations | `data/aggregate/primary_analysis/17C_Main_Table2_adjusted_clinical.csv` |
| Main Table 3 / Figure 2: targeted variants | `data/aggregate/primary_analysis/17D_Main_Table3_targeted_variants.csv` |
| Figure 1: cohort flow | `figures/Figure1_CohortFlow_v2_3.*`, caption, and locked cohort counts in `scripts/05_render_figures.R` |
| Figure 3: burden × F5/F2 observed proportions | `data/aggregate/primary_analysis/17E_Main_Figure3_data.csv` and `17E_Main_Figure3_chisquare.csv` |
| S1: internal discrimination and calibration | `data/aggregate/primary_analysis/17H_S1_*` |
| S2: pregnancy/postpartum timing | `data/aggregate/primary_analysis/17L_pregnancy_postpartum_reviewer_response.csv` |
| S4: VTE phenotype specification | `phenotypes/vte/` and `scripts/03_export_vte_concept_specification.py` |
| S5: clinical phenotype specification | `phenotypes/clinical/` and `scripts/04_reconstruct_clinical_concept_specification.py` |
| S6: formal interaction and model sensitivities | `data/aggregate/primary_analysis/17F_*` and `17G_S6_*` |
| S7: person-time / interval-specific sensitivity | `data/aggregate/primary_analysis/17I_*` and `17J_*` |
| S8: final ancestry sensitivity | `data/aggregate/ancestry_sensitivity/18B_*`, `18C_*`, and `18D_*` |

The nested Script 17 ZIP includes original run provenance and empty ancestry placeholders; the separate Script 18 ancestry files listed above are authoritative. A mapping above does not imply that every formatted supplemental item is reproduced as a separate document in this repository.

## What can be reproduced publicly

With R and the packages listed in `environment/R-packages.txt`, run from the repository root:

```bash
Rscript scripts/05_render_figures.R
```

This renders Figures 1–3 from the included aggregate Script 17 ZIP and writes to `outputs/figures/`. No participant-level model is fitted by the renderer. Existing final figures are available under `figures/` for immediate inspection.

Python syntax can be checked without secure inputs:

```bash
python -m py_compile scripts/*.py
```

On a system with GNU coreutils, verify the release file hashes with:

```bash
sha256sum --check SHA256SUMS.txt
```

## What requires controlled access

Re-fitting the primary and ancestry models requires an authorized All of Us Workbench session, the exact secure participant master expected by the scripts, and the relevant CDR. The master is not included. See the `VTE_JTH_MASTER_FILE` and `WORKSPACE_CDR` settings in the README. CDR access alone does not supply the missing local analysis master automatically.

The VTE and clinical specifications are count-free vocabulary exports. The clinical files represent a reconstruction from original rules, with missing historical optional manual overrides documented in `phenotypes/clinical/04_source_discovery.csv`. The package therefore supports review of the reconstructed specification, not a claim of verified recovery of all historical overrides.

## Interpretation reminders

- Five factors define the primary burden; pregnancy/postpartum is separate.
- The outcome is observed VTE occurrence during available follow-up.
- Figure 3 proportions are unadjusted; formal interaction analyses are in the `17F_*` files.
- Discrimination is internal; no external clinical utility has been established.
- The final ancestry analyses are Script 18 outputs, not Script 17 placeholders.

## Maintainer: update checksums after any file edit

Run this from the repository root before committing/tagging the final snapshot:

```bash
python - <<'PYCODE'
from pathlib import Path
import hashlib
root = Path('.')
excluded = {'.git', '__pycache__', 'outputs'}
files = sorted(p for p in root.rglob('*') if p.is_file()
               and not excluded.intersection(p.parts)
               and p.name != 'SHA256SUMS.txt')
Path('SHA256SUMS.txt').write_text(''.join(
    hashlib.sha256(p.read_bytes()).hexdigest() + '  ' + p.as_posix() + '\n'
    for p in files))
PYCODE
```

Do not regenerate the manifest just to silence a mismatch when verifying a downloaded release; investigate the difference first.
