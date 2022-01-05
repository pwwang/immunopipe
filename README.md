# immunopipe

Integrative analysis for scTCR- and scRNA-seq data

## Requirements & Installation

- `python`: `3.7+`
    - Other python depedencies should be installed via `pip install -U immunopipe`

- `R`
    - `immunarch`(`v0.6.7+`), `Seurat`(`v4.0+`), `scImpute`, `scran`, `scater`
    - `dplyr`, `tidyr`, `tibble`, `ggplot2`, `ggradar`, `ggprism`, `ggrepel`, `reshape2`
    - `ComplexHeatmap`, `RColorBrewer`
    - `future`, `parallel`, `gtools`
    - `enrichR`

- Other
  - VDJtools: https://vdjtools-doc.readthedocs.io/en/master/install.html

## Modules

- Basic TCR data analysis using `immunarch`
- Clone Residency analysis if you have paired samples (i.e. Tumor vs Normal)
- V-J usage, the frequency of various V-J junctions in circos-style plots
- Clustering cells and configurale arguments to separate T and non-T cells
- Clustering T cell, markers for each cluster and enrichment analysis for the markers
- Radar plots to show the composition of cells for clusters
- Markers finder for selected groups of cells
- Expression investigation of genes of interest for selected groups of cells
- UMAPs
- Metabolic landscape analysis (Ref: Xiao, Zhengtao, Ziwei Dai, and Jason W. Locasale. "Metabolic landscape of the tumor microenvironment at single cell resolution." Nature communications 10.1 (2019): 1-12.)

## Documentaion

https://pwwang.github.io/immunopipe
