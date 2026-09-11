# Publish the v1.0.1 correction

The public v1.0.0 release and Zenodo record remain available. This update prepares a new maintenance release.

1. Upload every file inside the supplied `UPLOAD_FILES` folder to the repository root, preserving subdirectories. On macOS, use Command + Shift + period to show hidden items. Confirm `.gitignore` and `.github/workflows/validate.yml` appear in the upload list.
2. Commit the upload to `main` with the message `Correct repository packaging for v1.0.1`.
3. Check the Actions tab: the restored workflow parses Python and R syntax. The local maintenance check does not replace that remote execution.
4. Confirm the checksum manifest matches the final tree. If files are edited again, regenerate it using the command in `REVIEWER_GUIDE.md` before tagging.
5. Confirm the intended Zenodo reuse license. The public v1.0.0 Zenodo record lists CC BY 4.0; this patch does not assign a new repository license.
6. Create a new GitHub release with tag `v1.0.1`, target `main`, and title `v1.0.1 — Documentation and packaging correction`. Use `RELEASE_NOTES.md` for its description. Keep the Zenodo integration enabled.
7. Confirm the new Zenodo version appears and use its version-specific DOI when citing the corrected snapshot. Preserve the v1.0.0 tag and published DOI.

No secure participant-level input, analysis rerun, or figure regeneration is part of this maintenance update.
