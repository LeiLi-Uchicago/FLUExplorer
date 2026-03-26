# Development Log - FLU Divergence Explorer

## Date: March 26, 2026

### 1. Performance & User Experience (UX) Enhancements
- **Global & Modal Loaders:** Implemented full-screen `waiter` loading screens across the app (Global context switches, Single Position Explorer, and Pairwise Comparison) to provide immediate visual feedback and freeze the UI during heavy data processing.
- **Lazy Loading Tabs:** Refactored reactive observers to utilize `session$clientData` hidden states, ensuring that heavy computations only trigger when a specific tab is visible, eliminating background lag.
- **Removed Cross-Tab Syncing:** Decoupled input dependencies between tabs to prevent reactive cascades and double-loading glitches when changing groups or subtypes.

### 2. Memory Management (Posit Server Optimization)
- **Real-Time Memory Monitor:** Added a floating widget to the UI to track RAM usage in real-time, including a manual "Clear Cache" button to gracefully release memory.
- **Aggressive Garbage Collection:** Implemented automatic cache clearing and garbage collection (`gc()`) upon user session termination and during cache evictions.
- **Cache Tuning:** Reduced the `get_lazy_table()` LRU cache size from 5 to 3 tables and added startup memory flushes to keep the application footprint strictly under Posit Server's 1GB memory limit.

### 3. Feature Deprecation
- **MSA Tab Removal:** Commented out the Multiple Sequence Alignment (MSA) tab and its associated heavy Bioconductor dependencies (`msaR`, `Biostrings`, `msa`) to substantially conserve RAM and improve app initialization times.

---

## Date: March 16, 2026

### 1. Migration from RSV to FLU (H1N1/H3N2)
- **Multi-Subtype Data Pipeline:** Overhauled `global.R` to support dynamic loading of multiple influenza subtypes. The app now merges metadata and usage tables from `data/H1N1/` and `data/H3N2/` subdirectories.
- **Neuraminidase (NA) Protein Fix:** Implemented a critical fix in `read_csv` calls (`na = character()`) to prevent the "NA" protein name from being interpreted as a logical missing value.
- **Robust Column Mapping:** Standardized inconsistent naming conventions across source files (e.g., mapping `Protein` to `Gene` and `HA_clade` to `Clade`) using length-stable `rename_with` logic.

### 2. UI/UX Global Overhaul
- **Centralized Header Controls:** Moved the Subtype selector to a high-visibility, fixed-position container in the top-right header, alongside the Data Mode (AA/NT) switch.
- **Streamlined Workflow:** Removed redundant subtype selectors from all individual tabs, enabling a "select once, explore everywhere" workflow.
- **Branding Refresh:** Updated all interface labels, home page content, and documentation to reflect the focus on Influenza Virus diversity.
- **Safe Dropdown Handling:** Applied `na.omit()` and safe `1:nrow()` checks to all selectors to prevent Shiny crashes when handling empty or partially missing data frames (e.g., `important_pos_df`).

### 3. Data Integrity & Visualization
- **Subtype-Specific Coloring:** Implemented explicit color mapping for subtypes (H1N1: Blue, H3N2: Red) in sequencing stats and geographical plots.
- **Enhanced Cache Logic:** Incremented the RDS cache version (`v2`) to force a rebuild, ensuring the "NA" protein fix and merged multi-subtype structure are correctly applied.

---

## Date: March 6, 2026 (Legacy RSV Version)

### 1. Dual Variation Support (AA & NT)
- **Standardized Data Infrastructure:** Implemented a unified data loading pipeline in `global.R` that handles both Amino Acid (AA) and Nucleotide (NT) datasets.
- **Reactive Data Switcher:** Integrated a reactive backend in `server.R` that dynamically swaps data sources based on user selection.

... (rest of legacy logs)
