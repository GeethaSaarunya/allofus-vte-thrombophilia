# v1.0.1 maintenance validation — September 11, 2026

Base public commit: `910a506ba66a4861eb1832bda6155fb6f06d1065`.

- Restored both omitted repository support files from the prepared source package.
- Retained one concise README phenotype-provenance section and its link to detailed reproducibility notes.
- Verified that all 56 scientific scripts, data, figures, phenotype files, and environment records are byte-identical to the public base commit.
- Parsed citation YAML and checked software version 1.0.1, repository URL, and all-versions DOI.
- Regenerated and verified the full release checksum manifest.
- Applied the upload patch to a clean copy of the public base tree and checked that it reconstructs the corrected tree exactly.
- No participant-level analyses or R figure rendering were rerun. No GitHub release or Zenodo record was modified by preparing this patch.

The entries below describe the original prepared package before its initial browser upload; they are retained as historical validation context. The public v1.0.0 upload omitted the two restored support files and later changed README text without refreshing its manifest. Those packaging differences are corrected here.

---

# Submission-release validation — September 11, 2026

| Check | Result |
|---|---|
| Citation metadata | Software 1.0.0; CFF format 1.2.0; YAML parsed, required fields and five authors checked; full external CFF schema validation not run |
| Python syntax | All four scientific Python scripts parsed successfully without execution |
| Scientific asset preservation | 56 files under scripts, data, figures, phenotypes, and environment are byte-identical to the attached GitHub package |
| Recursive ZIP inventory | Nested Script 17 payload inspected; no participant-table identifier columns or restricted participant-data file types found |
| Aggregate count scan | 51 CSVs scanned, including nested entries; 672 integer values in participant-count fields and 232 simple complements checked; no counts 1–20 found |
| Release documentation | Submission-stage release wording; no journal-acceptance timing gate |
| R execution / figure rendering | Not run here because R is unavailable; existing figure bytes preserved |
| Participant-level analytic rerun | Not performed; requires controlled access and secure master |
| Full inferential disclosure audit | Not performed by this bounded count scan |
| External publication | Not performed by this package update; no repository URL, DOI, release date, or license choice invented |

The count scan distinguishes participant counts from model coefficients, percentages, degrees of freedom, category labels, and vocabulary concept IDs. Its detailed output is in `submission_release_aggregate_audit.csv`. It checks straightforward denominator-minus-event complements, not every possible inference across released information. Existing dissemination responsibilities still apply.

After final metadata edits, regenerate the checksum manifest before tagging. A journal DOI is not needed for submission-stage release.
