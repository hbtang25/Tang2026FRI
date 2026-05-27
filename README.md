# Bioinformatic Pipeline for Deciphering Cantonese Soy Sauce Microbiome

[![R Version](https://img.shields.io/badge/R-%E2%89%A5%204.3.1-blue.svg)](https://www.r-project.org/)
[![License](https://img.shields.io/badge/License-MIT-green.svg)](https://opensource.org/licenses/MIT)

This repository hosts the standardized, reproducible computational workflow and specialized scripts developed for the manuscript: 
> **"Deciphering the correlation between functional microbiota and flavor metabolite formation in Cantonese soy sauce via microbiome analysis"** (Submitted to *Food Research International*).

The primary objective of this toolkit is to move beyond conventional descriptive sequencing surveys by providing an integrated, mathematically rigorous framework to quantify microbial succession, infer sparse ecological networks, and model microbial-flavor interactions under high-salinity fermentation stress.

---

## 🌟 Key Features & Methodological Rationale

Compared to traditional food microbiology monitoring workflows, this pipeline incorporates specific statistical normalization and modeling approaches to ensure maximum accuracy:

* **Compositionality Correction**: High-throughput amplicon data are inherently compositional. To address compositionality concerns and avoid spurious correlations, we implemented the **SPIEC-EASI** framework for microbial co-occurrence networks instead of traditional Pearson/Spearman coefficients.
* **Data Standardization**: Sample library sizes are normalized using **Total Sum Scaling (TSS)** (rescaling counts to $1 \times 10^6$) followed by **Hellinger transformation** to stabilize variance before multivariate analyses.
* **Predictive Modeling**: Incorporates **Random Forest (RF)** regression and **Redundancy Analysis (RDA)** with permutation tests (999 permutations) to establish statistical associations between key indicator taxa (biomarkers) and volatile organic compounds (VOCs).

---

## 📂 Repository Structure

```text
├── README.md               # Project documentation and reproduction guide
├── RDA.Rmd                 # Complete R Markdown script containing the full workflow
├── scripts/
│   ├── data_preprocessing.R # TSS normalization and Hellinger transformation
│   ├── spiec_easi_network.R # SPIEC-EASI network inference and topological analysis
│   └── random_forest_model.R # Machine learning predictive modeling for flavor profiles
└── data/
    └── anonymized_example/ # Structured toy datasets for pipeline validation
        ├── asv_table.csv    # Example ASV abundance matrix
        ├── metadata.csv     # Sample metadata (Fermentation stages: BFS, FFS, FS)
        └── flavor_vocs.csv  # Matched volatile compound concentrations
