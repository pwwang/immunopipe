# immunopipe

Complete Docker image for running the Immunopipe single-cell RNA-seq and TCR-seq analysis pipeline.

## Overview

`immunopipe` is the production-ready Docker image containing the complete Immunopipe pipeline for comprehensive single-cell immune repertoire analysis. This image combines all system dependencies, R packages, and Python pipeline code into a ready-to-use container for analyzing scRNA-seq and TCR-seq data.

## Purpose

This final image serves as:

- **Layer 3** in a 3-tier build strategy (`base` → `rpkgs` → `immunopipe`)
- **Production runtime** for executing complete single-cell analysis workflows
- **Reproducible environment** ensuring consistent results across different systems

## What is Immunopipe?

Immunopipe is a comprehensive pipeline built on the [pipen](https://github.com/pwwang/pipen) workflow framework that orchestrates:

- **Single-cell RNA-seq analysis**: Quality control, normalization, clustering, marker identification
- **TCR/BCR repertoire analysis**: Clonotype clustering, diversity metrics, CDR3 physicochemical properties
- **Cell type annotation**: Automated annotation using reference databases
- **Downstream analyses**: Gene set enrichment, metabolic profiling, cell-cell communication
- **Interactive reports**: HTML reports with visualizations and quality metrics

## Features

- **Dual workflow routing**: Automatic pathway selection based on TCR/BCR presence
- **Seurat-based**: Built on the widely-used Seurat ecosystem for scRNA-seq
- **Comprehensive TCR tools**: Integration of immunarch, scRepertoire, and ClusTCR
- **Batch processing**: Support for multi-sample and multi-condition experiments
- **Customizable**: Extensive configuration options via TOML files
- **Reproducible**: Containerized environment with pinned dependencies

## Contents

This image includes:

- All dependencies from `immunopipe-base` and `immunopipe-rpkgs`
- Python pipeline code installed via Poetry
- Process scripts (R and Python) in `/immunopipe/scripts/`
- Report templates in `/immunopipe/reports/`
- Configuration validation tools
- pipen-board for workflow monitoring

## Quick Start

### Pull the Image

```bash
docker pull justold/immunopipe:<tag>
```

### Run the Pipeline

```bash
docker run -v /path/to/data:/workdir \
  justold/immunopipe:<tag> \
  @config.toml
```

### Interactive Shell

```bash
docker run -it -v /path/to/data:/workdir \
  justold/immunopipe:latest \
  bash
```

## Input Requirements

Immunopipe accepts:

- **Sample metadata**: TSV file with sample information
- **scRNA-seq data**: Seurat objects (.RDS), 10x outputs, or count matrices
- **TCR/BCR data** (optional): Outputs from cellranger vdj or similar tools
- **Configuration**: TOML file specifying pipeline parameters

See the [documentation](https://pwwang.github.io/immunopipe/preparing-input/) for detailed input format specifications.

## Configuration

Configure the pipeline using TOML files:

```toml
[SeuratPreparing]
cell_qc = "nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 20"
use_sct = true

[SeuratClustering]
resolution = [0.4, 0.8, 1.2]

[MarkersFinder]
dbs = ["panglaodb", "cellmarker"]
```

Generate configuration templates:

```bash
docker run justold/immunopipe:latest \
  immunopipe --config-template > config.toml
```

## Output Structure

Pipeline outputs include:

- **Processed Seurat objects**: `.qs` files with clustered data
- **Visualizations**: UMAP plots, heatmaps, QC plots
- **Tables**: Marker genes, clonotype information, enrichment results
- **HTML reports**: Interactive reports for each analysis step
- **Logs**: Detailed execution logs for reproducibility

## Versioning

Images are tagged by branch:

- `justold/immunopipe:<version>` - Latest stable release (from `main` branch)
- `justold/immunopipe:dev` - Development version (from `dev` branch)

## Build Process

Built using a 3-stage strategy:

1. **immunopipe-base**: System dependencies and conda packages
2. **immunopipe-rpkgs**: Custom R packages from source
3. **immunopipe** (this image): Python pipeline code via Poetry

**Automated builds**: Triggered on every push to `main` or `dev` branches via GitHub Actions.

## Image Size

Approximately 8-12GB (complete with all dependencies and pipeline code).

## Advanced Usage

### With pipen-board (Web UI)

```bash
docker run -p 8521:8521 -v /data:/workdir \
  justold/immunopipe:<tag> \
  board immunopipe:Immunopipe
```

Then open http://localhost:8521 in your browser.

## Documentation

- **Full documentation**: https://pwwang.github.io/immunopipe/
- **GitHub repository**: https://github.com/pwwang/immunopipe
- **Getting started guide**: https://pwwang.github.io/immunopipe/getting-started/
- **Process documentation**: https://pwwang.github.io/immunopipe/processes/

## Citation

If you use Immunopipe in your research, please cite:

> Immunopipe: A comprehensive pipeline for single-cell immune repertoire analysis
> https://github.com/pwwang/immunopipe

## System Requirements

- Docker Engine 20.10+
- Minimum 8GB RAM (16GB+ recommended for large datasets)
- 50GB+ free disk space

## Troubleshooting

### Permission Issues

```bash
# Run with your user ID
docker run -u $(id -u):$(id -g) -v /data:/workdir \
  justold/immunopipe:<tag> @config.toml
```

### Memory Issues

```bash
# Increase Docker memory limit
docker run --memory=32g -v /data:/workdir \
  justold/immunopipe:<tag> @config.toml
```

## Related Images

- **justold/immunopipe-base**: Base dependencies layer
- **justold/immunopipe-rpkgs**: Intermediate layer with custom R packages

## License

See the [Immunopipe repository](https://github.com/pwwang/immunopipe) for license information.

## Support

For issues, questions, or feature requests:

- **GitHub Issues**: https://github.com/pwwang/immunopipe/issues
- **Documentation**: https://pwwang.github.io/immunopipe/
- **Discussions**: https://github.com/pwwang/immunopipe/discussions
