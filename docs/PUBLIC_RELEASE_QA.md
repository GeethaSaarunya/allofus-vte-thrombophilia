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
