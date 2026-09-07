# CLI utilities

Besides running the pipeline, `immunopipe` ships a set of utility commands for inspecting the pipeline outputs after a run:

```shell
$ immunopipe utils --help
```

## check-genes

Verify that gene symbols exist in the data before visualizing their expressions in [`SeuratClusterStats`](./processes/SeuratClusterStats.md). The Seurat object is read from the output of the [`SeuratPreparing`](./processes/SeuratPreparing.md) process.

```shell
$ immunopipe utils check-genes -w ./pipen/<pipeline name> -g CD3D,CD4,CD8A
```

### Options

- `-w, --workdir` (required): Working directory of the pipeline, typically `./pipen/<pipeline name>`
- `-g, --genes` (required): Comma-separated gene symbols to check, or a file path prefixed with `file://` containing gene symbols (one per line). Tab-separated files are also accepted, in which case only the first column is used.
- `--assay`: The assay to check the genes against (e.g. `RNA`, `ADT`). The genes are checked against the default assay if not specified.
- `--rscript`: Path to the `Rscript` executable (default: `Rscript`)

### Example

```shell
# Check comma-separated genes
$ immunopipe utils check-genes -w ./pipen/pipeline1 -g CD3D,CD4,CD8A

# Check genes from a file
$ immunopipe utils check-genes -w ./pipen/pipeline1 -g file://genes.txt

# Check genes in a specific assay
$ immunopipe utils check-genes -w ./pipen/pipeline1 -g CD3D,CD4 --assay ADT
```

## check-dim

Print the cell and gene counts before and after QC to inspect the effect of the QC filtering in [`SeuratPreparing`](./processes/SeuratPreparing.md). The output tables are the full contents of `qc/cell_qc.txt` and `qc/gene_qc.txt` in the `SeuratPreparing` output directory.

```shell
$ immunopipe utils check-dim -w ./pipen/<pipeline name>
```

### Options

- `-w, --workdir` (required): Working directory of the pipeline, typically `./pipen/<pipeline name>`
- `--rscript`: Path to the `Rscript` executable (default: `Rscript`)

## select-markers

Select the top marker genes for each cluster from the [`ClusterMarkers`](./processes/ClusterMarkers.md) process, for downstream analysis and visualization. The selected markers are written to stdout as a tab-separated table, so you can redirect it to a file:

```shell
$ immunopipe utils select-markers -o ./pipen/<pipeline name> > selected_markers.tsv
```

### Options

- `-o, --outdir` (required): Output directory of the pipeline (the directory containing the `ClusterMarkers` output)
- `-t, --top-n`: Number of top markers to keep per cluster (default: `10`)
- `--order-by`: An R expression to order the markers by before selecting the top ones. Available variables: `avg_log2FC`, `pct.1`, `pct.2`, `p_val`, `p_val_adj` (default: `desc(avg_log2FC)`)
- `-f, --filter`: An R expression to filter the markers before selection. Available variables: `avg_log2FC`, `pct.1`, `pct.2`, `p_val`, `p_val_adj` (default: `p_val_adj < 0.05`)
- `--rscript`: Path to the `Rscript` executable (default: `Rscript`)

/// Attention
The output directory must contain exactly one `ClusterMarkers/*.markers/` directory with exactly one case directory inside. If multiple are found, the command will list them and exit; specify the correct directory with `-o`.
///
