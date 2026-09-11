# Public release at submission: v1.0.0

The purpose of this release is to make the supporting materials accessible for editorial and peer review. Do not wait for journal acceptance or a journal DOI.

## Prepared in these packages

- [x] Set `CITATION.cff` software version to `1.0.0` while retaining CFF format version `1.2.0`.
- [x] Align README, Zenodo metadata, and release notes with submission-stage public release.
- [x] Include the count-free VTE and clinical specifications with QC and source-discovery documentation.
- [x] Preserve the qualification that historical clinical manual-override CSVs were not recovered.
- [x] Include aggregate source results, final ancestry results, and final Figures 1–3.
- [x] Add a reviewer guide mapping findings to files and explaining reproduction requirements.

Local validation results and their scope are recorded in `PUBLIC_RELEASE_QA.md`. Preparing these ZIP files does not itself publish them.

## Publish for review

1. Upload the **contents** of `allofus-vte-thrombophilia-reproducibility/` to the GitHub repository root. `README.md`, `CITATION.cff`, `scripts/`, and `data/` should appear at that root. Include the GitHub validation workflow from the GitHub ZIP.
2. Set the intended repository to **Public**, so reviewers can browse the files without an invitation.
3. Set the intended reuse/license terms before the Zenodo deposit. This packaging update does not assign a new license.
4. Add the actual GitHub URL as `repository-code` in `CITATION.cff`. Keep the software version `1.0.0`. Add `date-released` only when the actual release date is known. Regenerate `SHA256SUMS.txt` after any file edit, using the command in the reviewer guide.
5. Enable the repository in Zenodo **before** creating the GitHub release if using automatic archiving: https://help.zenodo.org/docs/github/enable-repository/
6. Commit the exact reviewed files, create tag **`v1.0.0`**, and publish a GitHub release titled **v1.0.0 — Manuscript submission and peer-review materials**. Use `RELEASE_NOTES.md` for the description. Confirm the repository's validation workflow succeeds.
7. Confirm the Zenodo record is published, downloadable, and has a DOI. If using a manual deposit instead, upload the equivalent Zenodo ZIP and use `ZENODO_METADATA.md`. Use one archiving route for this release rather than creating duplicate records.
8. Put the live repository URL and the DOI for the archived version in the manuscript Data Availability statement and journal submission. Test both links while signed out.

## After release

Record the GitHub tag/commit and Zenodo DOI. Add the DOI to subsequent citation metadata without silently rewriting the released snapshot; make any later file changes as a new release when needed. Add the article citation once available. A later corrected or reviewer-revised analysis should have its own version.
