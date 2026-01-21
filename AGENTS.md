# IMMUNOPIPE PROJECT KNOWLEDGE BASE

**Generated:** 2026-01-20 15:08 MST  
**Commit:** 95a0d60  
**Branch:** dev

## OVERVIEW

Bioinformatics pipeline for integrated scRNA-seq + scTCR/scBCR-seq analysis. Built on pipen workflow framework with multi-language execution (Python orchestration, R analysis scripts). Poetry-managed, Docker-ready, MCP-enabled for AI configuration generation.

## STRUCTURE

```
immunopipe/
├── processes.py          # ALL process definitions (1048 lines, 20+ processes)
├── pipeline.py           # Minimal orchestrator (20 lines) - delegates to processes.py
├── router.py             # Command dispatcher + TCR/non-TCR workflow branching
├── validate_config.py    # TOML config validation (220 lines)
├── gbatch.py             # Google Cloud Batch integration (211 lines)
├── inhouse.py            # Institution-specific processes (312 lines)
├── mcp/                  # Model Context Protocol server (2937 lines) - AI config generation
│   ├── tools.py          # MCP tool implementations (1107 lines)
│   ├── options.py        # Config discovery (486 lines)
│   ├── toml_generator.py # TOML generation utilities
│   └── doc_extractor.py  # Process documentation extraction
├── cli_utils/            # Utility commands (check-genes, check-dim)
├── scripts/              # R analysis scripts (5 files: TCellSelection, CloneHeterogeneity, etc.)
└── reports/              # Svelte report templates (4 files)

tests/running/configs/    # TOML config examples for integration tests
docs/processes/           # Auto-generated process docs (25+ files)
board.toml               # pipen-board web wizard config
```

## WHERE TO LOOK

| Task | Location | Notes |
|------|----------|-------|
| Add new process | `processes.py` + `scripts/` + `docs/processes/` | All processes in ONE file (pipen pattern) |
| Modify workflow routing | `router.py` | TCR vs non-TCR branching, command dispatch |
| Validate config | `validate_config.py` | TOML schema, error/warning system |
| Cloud execution | `gbatch.py` | Google Cloud Batch integration |
| AI config generation | `mcp/` | MCP server for natural language config |
| CLI utilities | `cli_utils/` | Gene/dimension checking tools |
| Entry points | `pyproject.toml` line 56, `__main__.py` | `immunopipe = "immunopipe.router:run"` |
| Build/deploy | `.github/workflows/` + `Dockerfile*` | 3-stage Docker, conditional CI |
| Test configs | `tests/running/configs/*.toml` | Real integration test examples |
| Dependencies | `pyproject.toml` + `docker/environment*.yml` | Poetry (Python) + Conda (R packages) |

## CONVENTIONS

### Workflow Framework (pipen)
- **Process-oriented architecture**: Each process = one analysis step with inputs/outputs/script
- **Single file for all processes**: `processes.py` contains ALL process definitions (not one-per-file)
- **Process ordering**: Use `order` attribute (e.g., `order = -1` runs first, `order = 99` runs last)
- **Conditional execution**: `@when()` decorator controls which processes run based on config
- **Process dependencies**: Defined via `requires` parameter in decorators

### Multi-Language Execution
- **Python**: Orchestration, workflow logic, data validation
- **R scripts**: Seurat/Bioconductor analysis (in `scripts/` directory)
- **Script naming**: Match process class name (e.g., `TCellSelection.R` for `TCellSelection` process)

### Configuration Management
- **TOML configs**: User-facing configuration format (see `tests/running/configs/`)
- **Validation**: `validate_config.py` catches misconfigurations before pipeline runs
- **board.toml**: Web wizard UI configuration for non-programmers
- **Process envs**: Each process accesses config via `envs` dictionary

### Entry Point System
- **Router-based CLI**: `router.py:run()` dispatches to 5 subcommands (gbatch, mcp, utils, help, pipeline)
- **No Click/Typer**: Custom command-line argument parsing
- **Hybrid frameworks**: Mixes argparse, argx, and pipen-args

### Testing
- **Integration tests**: Real TOML configs in `tests/running/configs/` (not just unit tests)
- **pipen-dry**: Use `pipen-dry` for workflow testing without execution
- **Fixtures**: Custom pytest fixtures in `conftest.py` for pipen setup

### Code Style
- **Black formatter**: 88 char line length (configured in pyproject.toml)
- **Flake8**: Ignores E203, W503, E731 (tox.ini)
- **No line length in docstrings**: `# noqa: E501` acceptable for documentation
- **Unused imports OK**: `# noqa: F401` for intentional exports in `__init__.py`

## ANTI-PATTERNS (THIS PROJECT)

**None found** - Codebase is clean of critical anti-patterns:
- ❌ No "DO NOT" / "NEVER" / "HACK" comments
- ❌ No TODO/FIXME/XXX comments
- ❌ No assertions in production code
- ✅ No hardcoded values (configuration-driven throughout)

## UNIQUE STYLES

### Process Definition Pattern
```python
# processes.py - ALL processes in one file
class MyProcess(Proc):
    """Process documentation (auto-extracted to docs/)"""
    input = "infile:file"
    output = "outfile:file"
    script = "Rscript {{proc.dir}}/scripts/MyProcess.R"
    order = 5  # Execution order
    envs = {"param1": "default"}  # Config params

@when("MyProcess" in config, requires=UpstreamProcess)
class ConditionalProcess(Proc):
    # Only runs if MyProcess configured
```

### MCP Server for AI Configuration
- **Unique feature**: Natural language → TOML config generation
- **Tools**: `configure_immunopipe`, `generate_full_config`, `suggest_processes`
- **Access**: `immunopipe mcp` starts server for Claude/GPT integration
- **2937 lines**: Largest module in project (shows commitment to AI usability)

### Hybrid Dependency Management
- **Poetry** (Python packages): pipen, biopipen, scanpy
- **Conda** (R packages + system deps): Seurat, Bioconductor, immunarch
- **Three-stage Docker**: base (conda) → rpkgs (R selective copy) → final (poetry install)

### Workflow Branching
- **Dynamic routes**: Router determines TCR vs non-TCR workflows at runtime
- **Conditional processes**: `@when()` decorator controls process execution
- **Not static**: Workflow adapts to input data structure (unlike traditional pipelines)

## COMMANDS

```bash
# Install
poetry install

# Run pipeline
immunopipe --config config.toml

# Subcommands
immunopipe gbatch       # Google Cloud Batch execution
immunopipe mcp          # Start MCP server for AI config generation
immunopipe utils        # Utility commands (check-genes, check-dim)
immunopipe help         # Process-specific help

# Testing
poetry run pytest tests/                       # All tests
poetry run pytest tests/test_processes.py      # Process tests
poetry run pytest tests/test_mcp_*.py          # MCP server tests

# Docker
docker build -t immunopipe .
docker run -v /data:/data immunopipe --config config.toml

# Documentation
mkdocs build           # Build docs
mkdocs serve           # Local docs server
```

## NOTES

### Numpy Version Constraint
**Critical**: `numpy==2.3` pinned in Dockerfile. scanpy needs numba, which requires numpy < 2.4. Monitor numba compatibility updates.

### Process Discovery
Processes are auto-discovered from `processes.py`. No manual registration needed—pipen framework introspects process classes.

### biopipen Wrapper Pattern
Most processes imported from `biopipen` package. Immunopipe = configuration layer + immunology-specific customizations, not a standalone tool.

### MCP Server Logging
MCP server uses `logging.WARNING` level to avoid stdout interference (stdio transport requirement).

### Validation Before Execution
**Always run validation**: `validate_config.py` catches configuration errors early. Saves hours of compute time on misconfigured runs.

### Router Workflow Visualization
See `docs/routes-tcr.png` and `docs/routes-notcr.png` for visual workflow diagrams. Essential for understanding branching logic.
