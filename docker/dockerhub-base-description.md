# immunopipe-base

Base Docker image for the Immunopipe single-cell RNA-seq and TCR-seq analysis pipeline.

## Overview

`immunopipe-base` is the foundational layer in a multi-stage Docker build system for [Immunopipe](https://github.com/pwwang/immunopipe). This image contains all core system dependencies, Python packages, R base environment, and Bioconductor packages required by the pipeline, providing a stable foundation that changes infrequently.

## Purpose

This base image serves as:

- **Layer 1** in a 3-tier build strategy (`base` → `rpkgs` → `immunopipe`)
- **Dependency foundation** containing ~100+ conda packages including Python 3.12, R 4.3, and essential bioinformatics tools
- **Cache optimization** layer that rebuilds only when core dependencies change, significantly reducing build times

## Contents

### System Tools
- Node.js 22 (for pipen-report)
- Build essentials (gcc, git, curl, wget)
- System libraries (fontconfig, graphviz, libiconv)
- Utilities (procps-ng, time)

### Python Ecosystem
- Python 3.12
- Poetry (package management)
- Core scientific packages (via pip): liana, celltypist2, tensorflow 2.20, keras, scikit-learn, biopython

### R Ecosystem
- R 4.3 base
- Bioconductor core packages (fgsea, scater, glmgampoi, destiny, DESeq2, edgeR)
- Seurat ecosystem (seurat, seuratobject, seuratdisk, seuratwrappers)
- Single-cell analysis tools (doubletfinder, scdblfinder, scRepertoire)
- Visualization packages (ggplot2, complexheatmap, patchwork, plotroc)
- TCR/BCR analysis (immunarch via scRepertoire dependencies)
- Statistical packages (harmony, presto, glmnet, metap)

### Specialized Tools
- FAISS (CPU version 1.13) for similarity search
- ClusTCR 2025 for TCR clustering
- HDF5R for reading 10x data

## Build Process

This image is built from `Dockerfile.base` using `environment_base.yml` as the conda environment specification. The build includes aggressive cleanup to minimize image size:

- Removal of static libraries, headers, and build tools
- Deletion of Python bytecode and JavaScript source maps
- Stripping of R package documentation and help files
- Removal of compiler compatibility files

**Automated builds**: Triggered automatically via GitHub Actions when `environment_base.yml` or `Dockerfile.base` changes on the `dev` branch.

## Usage

### As a Base for immunopipe-rpkgs

```dockerfile
FROM justold/immunopipe-base:latest

# Install custom R packages built from source
COPY --from=builder /opt/conda/lib/R/library/gglogger /opt/conda/lib/R/library/gglogger
# ...additional R packages...
```

### For Development/Testing

```bash
docker pull justold/immunopipe-base:latest
docker run -it justold/immunopipe-base:latest bash
```

## Image Size

Optimized through multi-stage cleanup, typically ~5-8GB (before R custom packages).

## Maintenance

- **Source**: [immunopipe/docker/environment_base.yml](https://github.com/pwwang/immunopipe/blob/main/docker/environment_base.yml)
- **Dockerfile**: [immunopipe/Dockerfile.base](https://github.com/pwwang/immunopipe/blob/main/Dockerfile.base)
- **Update frequency**: Changes only when core dependencies need version updates or new system-level packages are required

## Related Images

- **justold/immunopipe-rpkgs** (builds on this image): Adds custom R packages (gglogger, enrichit, plotthis, scplotter, biopipen.utils)
- **justold/immunopipe** (final image): Complete pipeline with Python code and configuration

## License

See the [Immunopipe repository](https://github.com/pwwang/immunopipe) for license information.

## Support

For issues or questions:
- GitHub Issues: https://github.com/pwwang/immunopipe/issues
- Documentation: https://pwwang.github.io/immunopipe/
