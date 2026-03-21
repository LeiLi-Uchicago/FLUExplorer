# 🧬 FLU Amino Acid Divergence Explorer

![R Version](https://img.shields.io/badge/R-%3E%3D%204.0.0-blue)
![Shiny](https://img.shields.io/badge/Built_with-R_Shiny-success)
![Bioinformatics](https://img.shields.io/badge/Field-Bioinformatics-purple)
![License](https://img.shields.io/badge/License-MIT-green)

**FLU Amino Acid Divergence Explorer** is a high-resolution, interactive web application designed for the genomic analysis of niun erin
---

## ✨ Key Features

### 🌍 1. Dataset Insights
* **Global Control C:Iab
* **Temporal & Regional Breakdown:** Interactive Plotly charts showing the seasonality of sequencing efforts and geographic hotspots.

### 🔬 2. Single Position Explorer
* Dive deep into the amino acid or nucleotide distribution of any specific position within an Influenza gene (HA, NA, etc.).
* Easily toggle between raw sequence counts and relative frequencies.
* Export publication-ready plots (PNG/PDF) and conditionally-formatted Excel data matrices.

### ⚖️ 3. Pairwise Comparison
* Instantly identify robust, fixed amino acid differences between any two selected viral clades.
* Set custom consensus thresholds (e.g., >90% dominant frequency) to filter out background noise.

### 🏔️ 4. Gene-Wide Landscapes
* **Conservation (Entropy):** Calculates Shannon Entropy to map hypervariable peaks and highly conserved valleys across an entire gene.
* **Mutation Tracker (Lollipop):** Generates staggered lollipop plots to visualize fixed amino acid mutations in a newly emerged target clade against an ancestral reference clade.
* **Interactive Consensus MSA:** Utilizes `msaR` to perform real-time alignments of consensus sequences, highlighting exact regions of structural divergence.

---

## 🚀 Getting Started

### Prerequisites
To run this application locally, you will need **R (>= 4.0.0)** installed on your machine.

```R
# 1. Install CRAN packages
install.packages(c("shiny", "dplyr", "ggplot2", "DT", "readr", "tidyr", 
                   "openxlsx", "plotly", "waiter", "lubridate", "tidyverse", 
                   "leaflet", "leaflet.minicharts", "shinyWidgets", "shinyjs"))

# 2. Install Bioconductor packages
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("Biostrings", "msa", "msaR"))
```

### Installation & Execution
1. Clone the repository:
```Bash
git clone https://github.com/LeiLi-Uchicago/FLUExplorer.git
cd FLUExplorer
```
2. Open the project in RStudio.
3. Run the App: Open `global.R` and click the "Run App" button.

## 📁 Repository Structure
```Plaintext
├── global.R                  # Multi-subtype data pipeline and color mapping
├── ui.R                      # Global control UI and customized styling
├── server.R                  # Backend logic and interactive visualizations
├── data/                     # Organized by Subtype
│   ├── H1N1/                 # Metadata and usage tables for H1N1
│   ├── H3N2/                 # Metadata and usage tables for H3N2
│   └── app_cache_flu_v2.rds  # Optimized binary cache for instant loading
└── www/                      # Static web assets
```

## ✍️ Authors & Citation
Lei Li - Initial work & Development.

## 📜 License
This project is licensed under the MIT License.
