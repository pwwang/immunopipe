# immunopipe

Integrative analysis for scTCR- and scRNA-seq data

## Requirements & Installation

- `python`: `3.7+`
    - Other python depedencies should be installed via `pip install -U immunopipe`

- `R`
    - A bunch of R packages

- Other
  - VDJtools: https://vdjtools-doc.readthedocs.io/en/master/install.html

- Checking requirements

  ```shell
  pip install -U pipen-cli-require
  pipen require immunopipe.pipeline:pipeline <pipeline arguments>
  ```

- Quick way to install the dependencies using conda
  ```shell
  conda env create --file docker/environment.yml
  # then
  conda activate immunopipe
  ```

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
#   -B .:/workdir,/path/to/data:/mnt
# Where /path/to/data is the data directory containing the data files
# You may also want to bind other directories (i.e. /tmp)
#   -B <other bindings>,/tmp

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
- (Meta-)Markers finder for selected groups/clones of cells
- Psedo-bulk GSEA analysis of two groups of cells
- Seurat cluster statistics, including:
  - Basic statistics of the clusters (e.g. number of cells in each cluster)
  - Gene expressions (e.g. ridge, violin, feature, dot and heatmap plots)
  - Dimensional reduction plots
- TCR clustering using CDR3 sequences and the statistics of the clusters
- Cell group distribution (TCR clones/clusters) in Seurat clusters
- Clone heterogeneity (TCR clone distribution) in Seurat clusters
- Metabolic landscape analysis (Ref: Xiao, Zhengtao, Ziwei Dai, and Jason W. Locasale. "Metabolic landscape of the tumor microenvironment at single cell resolution." Nature communications 10.1 (2019): 1-12.)

## Documentaion

https://pwwang.github.io/immunopipe
