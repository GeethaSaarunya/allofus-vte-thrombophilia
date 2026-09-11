# All of Us data-use and dissemination note

This repository must not contain participant-level All of Us data.

The All of Us Data and Statistics Dissemination Policy prohibits publishing or distributing participant counts from 1 through 20 and also prohibits releasing information from which such small counts can be derived. The repository package therefore contains only aggregate-safe outputs selected from the manuscript workflow.

Before every public GitHub or Zenodo release:

1. Confirm that no participant-level files are present.
2. Confirm that no `person_id`, `research_id`, dates tied to individual participants, or participant-level predictions are present.
3. Re-audit aggregate tables for counts from 1 through 20 or derivable small cells.
4. Do not upload secure `.parquet`, `.feather`, pickle, or raw Workbench exports.
5. Follow the current All of Us publication/presentation notification requirements.

Relevant policy pages:
- https://www.researchallofus.org/faq/data-and-statistics-dissemination-policy/
- https://www.researchallofus.org/faq/publication-and-presentation-policy/
