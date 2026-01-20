# IMMUNOPIPE/CLI_UTILS/ KNOWLEDGE BASE

**Generated:** 2026-01-20
**Lines:** 234 (4 Python + 2 R files)

## OVERVIEW
Utility commands for post-run inspection of Seurat objects and QC metrics.

## WHERE TO LOOK

| Task | Location | Purpose |
|------|----------|---------|
| check-genes utility | check_genes.py + check_genes.R | Verify gene symbols exist in SeuratPreparing output |
| check-dim utility | check_dim.py + check_dim.R | Display cell/gene counts before/after QC |
| R script runner | utils.py:run_r() | Execute R code via temporary files |
| Entry point | __init__.py:main() | Router "utils" subcommand dispatcher |

## COMMANDS

### check-genes
Verify gene availability before SeuratClusterStats visualization.

```bash
# Check comma-separated genes
immunopipe utils check-genes -w ./pipen/pipeline1 -g CD3D,CD4,CD8A

# Check genes from file
immunopipe utils check-genes -w ./pipen/pipeline1 -g file://genes.txt

# With specific assay
immunopipe utils check-genes -w ./pipen/pipeline1 -g CD3D,CD4 --assay RNA
```

### check-dim
Inspect QC filtering effectiveness by comparing cell/gene counts.

```bash
immunopipe utils check-dim -w ./pipen/pipeline1
```

**Output**: Prints full qc/cell_qc.txt and qc/gene_qc.txt tables showing dimensions before/after filtering.

## CONVENTIONS

### R Script Runner Pattern
```python
from utils import run_r
script = ("source('script.R')", "func(arg1='val1', arg2='val2')")
rc, stdout, stderr = run_r("Rscript", script)
```
- Pass script as tuple/list of lines (auto-joined with newline)
- Returns tuple: (returncode, stdout, stderr)
- Temporary file auto-deleted after execution

### Signature File Resolution
Utilities locate Seurat objects via job.signature.toml in workdir/SeuratPreparing/0/
- Handles Google Batch container path conversion (/mnt/disks/.cwd → local)
- Resolves symlinks automatically
