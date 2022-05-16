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

## Running as a container

### Using docker:

```bash
docker run -w /workdir -v .:/workdir -it justold/immunopipe:dev
```

### Using singularity:

```bash
singularity run -w \  # need it to be writable
  -H /home/immunopipe_user \  # required, used to init conda
  --pwd /workdir -B .:/workdir \  # Could use other directory instead of "."
  # --contain: don't map host filesystem
  # --cleanenv: recommended, to avoid other host's environment variables to be used
  #   For example, $CONDA_PREFIX to affect host's conda environment
  --contain --cleanenv \
  docker://justold/immunopipe:dev

# The mount your data directory to /mnt, which will make startup faster
# For example
# -B .:/workdir,/path/to/data:/mnt
# Where /path/to/data is the data directory containing the data files

# Or you can pull the image first by:
singularity pull --force --dir images/ docker://justold/immunopipe:dev
# Then you can replace "docker://justold/immunopipe:dev" with "images/immunopipe.sif"
```

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
