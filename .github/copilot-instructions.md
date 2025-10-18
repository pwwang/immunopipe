# Immunopipe AI Coding Agent Instructions

## Project Overview
Immunopipe is a comprehensive single-cell RNA-seq and TCR-seq analysis pipeline built on the **pipen** workflow framework. It orchestrates complex bioinformatics workflows for immune repertoire analysis, cell type annotation, metabolic profiling, and cell-cell communication.

## Architecture

### Core Components
- **`immunopipe/pipeline.py`**: Main pipeline orchestrator defining the workflow DAG
- **`immunopipe/processes.py`**: Individual process definitions (each process = one analysis step)
- **`immunopipe/router.py`**: Workflow routing logic based on input data types (TCR/non-TCR routes)
- **`immunopipe/gbatch.py`**: Batch processing utilities for multiple samples
- **`immunopipe/inhouse.py`**: Institution-specific customizations
- **`immunopipe/scripts/`**: R/Python analysis scripts executed by processes
- **`immunopipe/reports/`**: Report generation templates (likely Rmarkdown/HTML)
- **`immunopipe/mcp/`**: Model Context Protocol server for AI agent integration

### Workflow Framework (pipen)
- Each "process" is a step in the pipeline with inputs, outputs, and a script
- Processes are connected via dependencies forming a DAG
- Configuration via `board.toml` (pipen-board UI config) and TOML config files
- Pipeline execution creates output directories with intermediate results

### Key Data Flow
1. Input: Sample metadata + scRNA-seq data (Seurat objects, counts matrices)
2. Preprocessing: `SeuratPreparing` → QC and normalization
3. Branching based on TCR presence (see `router.py`):
   - **TCR route**: TCR clustering → T-cell selection → specialized analyses
   - **No-TCR route**: Standard scRNA-seq analysis
4. Core analyses: Clustering → marker identification → annotation → downstream (GSEA, metabolism, communication)
5. Outputs: Processed Seurat objects, plots, HTML reports, tables

## Development Workflows

### Setup and Testing
```bash
# Install with Poetry (preferred)
poetry install

# Run tests
poetry run pytest tests/

# Build Docker image
make docker  # or docker build -t immunopipe .

# Run with Docker
docker run -v /data:/data immunopipe [args]
```

### Adding a New Process
1. Define process class in `immunopipe/processes.py` (inherit from `Proc` or similar base)
2. Add corresponding script in `immunopipe/scripts/` (R or Python)
3. Update `immunopipe/pipeline.py` to wire the process into the DAG
4. Add documentation in `docs/processes/ProcessName.md`
5. Update routing logic in `immunopipe/router.py` if the process is route-specific

### Configuration Validation
- `immunopipe/validate_config.py` validates user configurations
- Process-specific configs are defined in process classes (check `envs` attributes)
- See `docs/configurations.md` for user-facing config documentation

## Project-Specific Conventions

### File Organization
- **Process scripts naming**: `immunopipe/scripts/{ProcessName}.R` or `.py` (matches process class name)
- **Process docs**: `docs/processes/{ProcessName}.md` (one doc per process)
- **Tests**: `tests/test_{module}.py` (pytest-based)

### Coding Patterns
- **Process definitions**: Use pipen's process class system with `input`, `output`, `script` attributes
- **R scripts**: Typically use Seurat ecosystem (Seurat, SeuratPlus, immunarch for TCR)
- **Python scripts**: Use scanpy/anndata for scRNA-seq, pandas for data manipulation
- **Config access**: Processes access config via `envs` dictionary passed from pipeline

### Dependencies
- **Python**: pipen framework, pydantic (validation), pandas, scanpy
- **R**: Seurat, tidyverse, immunarch (check `docker/environment.yml` for conda env)
- **Containerization**: Use Dockerfile.full for complete environment with all R/Python deps

## Critical Knowledge

### Router System
The `router.py` determines workflow paths based on:
- Presence of TCR/BCR data in input (`has_tcr` flag)
- Sample groupings and batch structures
- See `docs/routes-tcr.png` and `docs/routes-notcr.png` for visual workflow diagrams

### MCP Server Integration
- `immunopipe/mcp/` implements Model Context Protocol for AI agent interaction
- Allows AI agents to query pipeline state, generate configs, and understand process relationships
- Test MCP functionality with `tests/test_mcp_*.py`

### Testing Strategy
- **Unit tests**: Individual process logic (`test_seuratpreparing.py`, `test_seuratclustering.py`)
- **Integration tests**: Sample data workflows in `tests/running/`
- **Config validation**: `test_mcp_options.py`, `test_mcp_toml_generator.py`

### Documentation Generation
- Uses mkdocs (see `mkdocs.yml`)
- Process documentation auto-extracted from docstrings (see `test_doc_extractor.py`)
- Build docs: `mkdocs build` (requires deps from `docs/requirements.txt`)

## Common Pitfalls
- **Process dependencies**: Ensure upstream processes are properly wired in pipeline.py before downstream processes
- **R/Python script paths**: Scripts must be findable relative to process execution context
- **Configuration schema**: When adding new process options, update `validate_config.py` schemas
- **Docker builds**: Use `Dockerfile.full` for complete reproducibility; base `Dockerfile` is minimal
- **TCR-specific processes**: Check router.py to ensure processes only run when TCR data is present

## Key Files to Reference
- **Pipeline structure**: `immunopipe/pipeline.py`
- **Process catalog**: `immunopipe/processes.py` + `docs/processes/`
- **Workflow routing**: `immunopipe/router.py`, `docs/routes-*.png`
- **Dependencies**: `pyproject.toml`, `docker/environment*.yml`
- **User docs**: `docs/getting-started.md`, `docs/running.md`, `docs/preparing-input.md`
