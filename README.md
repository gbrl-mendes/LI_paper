## eDNA Metabarcoding analyses from "Hydrology-Induced Changes Drive Fish eDNA Ecology in a Neotropical Reservoir" paper

### About
This repository hosts reproducible Quarto workflows for the integrative curation, exploration, and analysis of fish eDNA metabarcoding data from Ingleses Lake (Brazil) eDNA metabarcoding project. It standardizes pipeline steps from raw tables to curated taxonomic assignments, exploratory plots, alpha-diversity, multivariate ordinations, and summary tables for publication.

Key deliverables include:

- Curated, metadata-enriched results tables

- Filtering diagnostics and removed-IDs audit tables

- Final species-level summaries across sites, years, and water-level periods

- Exploratory visualizations and ordination analyses

- Ready-to-use figures and tables for manuscripts and supplements

### Repository structure
- scripts/

  - LI_script.qmd — end-to-end curation and data integration, metadata join, filtering, diagnostics, summary tables, modular exploratory plots, alpha-diversity, PCoA, PERMANOVA, accumulation curves, and diversity indices.
  - LI_script.html — HMTL version of LI_script.qmd.
  - LI_tree.qmd — Neighbor-joining tree from ASVs.

- data/

  - curated/ and raw/ inputs as referenced inside the scripts (e.g., curated identifications, BLAST taxonomy, sequencing read counts).

- results/

  - figures/ — exported PDF figures.

  - tables/ — CSV/XLSX exports (summaries, taxonomic levels, alpha-diversity summaries).

- environment/

  - Optional saved environments for quick restoration of a known-good state.

**Paths are defined in-script; adjust to your local environment if folder layout differs.**

### Citation
If this workflow is useful for you, or your academic work, please cite the repository and credit the author of the scripts 😊. 

BLAST references, databases, and R packages should be cited according to their individual guidelines.

### Contact
For more information, contact: [e-mail](mailto:gabrielmendesbrt@gmail.com)
