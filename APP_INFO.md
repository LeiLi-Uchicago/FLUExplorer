# FLUExplorer Methods and Information

## Overview

FLUExplorer is an interactive Shiny application for exploring amino acid variation across Influenza genomes. It is designed to help users examine subtype-specific sequence diversity, compare genetic clades, inspect variation at individual positions, and visualize gene-wide patterns such as conservation and fixed mutations.

The current version of the app is loaded with curated human Influenza A subtype H1N1, H3N2, human Influenza B Yam and Vic lineage amino acid usage data organized by gene, genetic clades, year, and year-month.

FLUExplorer is intended for researchers, bioinformaticians, genomic epidemiologists, and other users who want to investigate FLU evolutionary patterns, mutation dynamics, and lineage-specific amino acid changes through an accessible visual interface.

---

## Data Processing

Data were sourced from GISAID, with all sequences annotated via Nextclade 3. The reference strains utilized for each lineage are listed below: 

| Pathogen/Lineage     | Reference Strain           | GISAID Accession      |
| -------------------- | -------------------------- | --------------------- |
| A(H1N1)pdm09         | A/California/04/2009       | `EPI_ISL_393964`      |
| A(H3N2)              | A/Darwin/9/2021            | `EPI_ISL_2233240`     |
| B/Victoria           | B/Austria/1359417/2021     | `EPI_ISL_1519459`     |
| B/Yamagata           | B/Brisbane/9/2014          | `EPI_ISL_165595`      |


---

## Update Log

### 2026-05-05

- Added the Genetic Clade tab. Users can select a subtype-specific clade annotation, search ranked clade/group choices, and review clade-level summary statistics.

### 2026-04-28

- Added DuckDB-backed usage table loading to reduce memory pressure from large amino acid and nucleotide tables. When DuckDB is available, FLUExplorer now queries only the selected subtype, gene, position, grouping, and time range instead of loading full tables into memory.
- Improved the Single Position Explorer time filter by replacing the draggable Year-Month range slider with explicit Start and End selectors.
- Fixed Year-Month grouping so plots and tables display real `YYYY-MM` groups instead of collapsing into `Unknown`.
- Fixed a Year-Month plotting error that could appear when a selected position had sparse or special time values.

### 2026-04-10

- Updated UI.

### 2026-03-06

- 1st version on-line. Using NextStrain data (Downloaded on 2025).
