# Bioinformatic Pipeline for Deciphering Cantonese Soy Sauce Microbiome

[![R Version](https://img.shields.io/badge/R-%E2%89%A5%204.0.3-blue.svg)](https://www.r-project.org/)
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
├── 16S_moromi.sh                 # Complete shell script containing the 16S workflow
├── ITS_moromi.sh                 # Complete shell script containing the ITS workflow
├── network.Rmd                 # SPIEC-EASI network inference and topological analysis
├── rf_classification.Rmd                 # Machine learning predictive modeling for microbial feature
├── scripts/                         # Other auxiliary scripts
└── data/
    ├── otutab_bac.txt    # ASV abundance matrix for bacteria
    ├── otutab_fun.txt    # ASV abundance matrix for fungi
    ├── taxonomy_bac.txt    # taxonomy for bacterial ASV
    ├── taxonomy_fun.txt    # taxonomy for fungal ASV
    ├── metadata_bac.tsv     # Sample metadata (Fermentation stages: BFS1-BFS3)
    ├── metadata_fun.tsv     # Sample metadata (Fermentation stages: FFS1-FFS3)
```

<details>
<summary>🛠️ 点击展开查看开发环境与 R 包版本 (sessionInfo)</summary>

```R
R version 4.0.3 
Platform: x86_64-w64-mingw32/x64
Running under: Windows 10 x64 (build 19044)

Matrix products: default
  LAPACK version 3.12.1

locale:
[1] LC_COLLATE=Chinese (Simplified)_China.utf8  LC_CTYPE=Chinese (Simplified)_China.utf8    LC_MONETARY=Chinese (Simplified)_China.utf8
[4] LC_NUMERIC=C                                LC_TIME=Chinese (Simplified)_China.utf8    

time zone: Asia/Shanghai
tzcode source: internal

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] vegan_2.7-3           permute_0.9-10        ggsci_4.2.0           pROC_1.19.0.1         pheatmap_1.0.13      
 [6] ggplot2_3.5.4         randomForest_4.6-14  magrittr_2.0.4        SpiecEasi_1.99.0      Hmisc_5.2-5          
[11] dplyr_1.2.0           igraph_1.6.0          WGCNA_1.74            fastcluster_1.3.0     dynamicTreeCut_1.63-1
[16] reshape2_1.4.5        psych_2.6.3          

loaded via a namespace (and not attached):
 [1] shape_1.4.6.1         gtable_0.3.6          impute_1.84.0         xfun_0.57             htmlwidgets_1.6.4    
 [6] lattice_0.22-9        vctrs_0.7.2           tools_4.5.0           generics_0.1.4        stats4_4.5.0         
[11] parallel_4.5.0        tibble_3.3.1          cluster_2.1.8.2       pkgconfig_2.0.3       Matrix_1.7-5         
[16] huge_1.5              data.table_1.18.2.1   checkmate_2.3.4       RColorBrewer_1.1-3    S7_0.2.1             
[21] lifecycle_1.0.5       compiler_4.5.0        farver_2.1.2          stringr_1.6.0         mnormt_2.1.2         
[26] codetools_0.2-20      htmltools_0.5.9       glmnet_4.1-10         yaml_2.3.12           htmlTable_2.4.3      
[31] preprocessCore_1.72.0 Formula_1.2-5         pillar_1.11.1         MASS_7.3-65           iterators_1.0.14     
[36] rpart_4.1.24          foreach_1.5.2         nlme_3.1-168          tidyselect_1.2.1      digest_0.6.39        
[41] stringi_1.8.7         splines_4.5.0         fastmap_1.2.0         grid_4.5.0            colorspace_2.1-2     
[46] cli_3.6.5             base64enc_0.1-6       dichromat_2.0-0.1     survival_3.8-6        withr_3.0.2          
[51] foreign_0.8-91        scales_1.4.0          backports_1.5.0       rmarkdown_2.31        matrixStats_1.5.0    
[56] otel_0.2.0            nnet_7.3-20           gridExtra_2.3         pulsar_0.3.12         VGAM_1.1-14          
[61] evaluate_1.0.5        knitr_1.51            doParallel_1.0.17     mgcv_1.9-4            rlang_1.1.7          
[66] Rcpp_1.1.1            glue_1.8.0            rstudioapi_0.18.0     R6_2.6.1              plyr_1.8.9  
```
</details>