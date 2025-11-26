# immunopipe-rpkgs

Intermediate Docker image for Immunopipe containing custom R packages built from source.

## Overview

`immunopipe-rpkgs` is the second layer in Immunopipe's multi-stage Docker build system. This image extends `immunopipe-base` by adding custom R packages that are not available in standard conda/CRAN repositories, providing specialized visualization and analysis capabilities for single-cell genomics.

## Purpose

This intermediate image serves as:

- **Layer 2** in a 3-tier build strategy (`base` → `rpkgs` → `immunopipe`)
- **Custom R package layer** containing GitHub-sourced packages built from source
- **Build optimization** layer that separates slow-changing custom packages from frequently updated pipeline code

## Contents

This image includes all dependencies from `immunopipe-base` plus the following custom R packages:

### Custom R Packages

- **gglogger** - Enhanced logging utilities for ggplot2-based workflows
- **enrichit** - Gene set enrichment analysis utilities
- **plotthis** - Comprehensive visualization toolkit for single-cell data
- **scplotter** - Specialized single-cell plotting functions
- **biopipen.utils** - Utility functions for the biopipen ecosystem

All packages are installed from source and compiled during the build process to ensure compatibility with the base environment.

## Build Process

Built using a multi-stage approach:

1. **Builder stage**: Installs packages from `environment_rpkgs.yml` using micromamba
2. **Final stage**: Copies only compiled R libraries from builder to `immunopipe-base`, avoiding metadata bloat
3. **Cleanup**: Removes conda package cache and build artifacts

**Automated builds**: Triggered via GitHub Actions when:
- `environment_rpkgs.yml` changes
- `Dockerfile.rpkgs` changes
- Base dependencies change (`environment_base.yml` or `Dockerfile.base`)

All triggers require changes on the `dev` branch.

## Usage

### As a Base for immunopipe

```dockerfile
FROM justold/immunopipe-rpkgs:latest

COPY --chown=$MAMBA_USER:$MAMBA_USER . /immunopipe

# Install Python pipeline code with Poetry
RUN python -m poetry install --no-cache -v --all-extras
```

### For Development/Testing

```bash
docker pull justold/immunopipe-rpkgs:latest
docker run -it justold/immunopipe-rpkgs:latest R

# Test custom packages in R
> library(plotthis)
> library(scplotter)
```

## Build Strategy Rationale

**Why separate custom R packages?**

1. **Build time optimization**: R packages from source can take 30-60 minutes to compile; separating them from frequently-changed Python code reduces rebuild times
2. **Layer caching**: Custom R packages change infrequently compared to pipeline code
3. **Size efficiency**: Multi-stage build avoids including conda metadata and build tools in final layer

## Image Size

Adds ~1-2GB to the base image size (total ~6-10GB).

## Maintenance

- **Source**: [immunopipe/docker/environment_rpkgs.yml](https://github.com/pwwang/immunopipe/blob/main/docker/environment_rpkgs.yml)
- **Dockerfile**: [immunopipe/Dockerfile.rpkgs](https://github.com/pwwang/immunopipe/blob/main/Dockerfile.rpkgs)
- **Update frequency**: Changes when custom R packages need version updates or new packages are added

## Package Sources

Custom R packages are typically sourced from:
- pwwang's GitHub repositories (gglogger, plotthis, scplotter)
- Bioconductor development versions
- Specialized bioinformatics tools

## Related Images

- **justold/immunopipe-base** (parent): Base system and conda dependencies
- **justold/immunopipe** (child): Complete pipeline with Python code

## License

See the [Immunopipe repository](https://github.com/pwwang/immunopipe) for license information.

## Support

For issues or questions:
- GitHub Issues: https://github.com/pwwang/immunopipe/issues
- Documentation: https://pwwang.github.io/immunopipe/
