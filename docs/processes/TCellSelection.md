# TCellSelection

Separate T and non-T cells and select T cells.

If all of your cells are T cells, do not set any configurations for this process.<br />

In such a case, [`SeuratClusteringOfAllCells`](SeuratClusteringOfAllCells.md) will
be clustering all T cells and
[`SeuratClusteringOfTCells`](SeuratClusteringOfTCells.md) will be skipped.<br />

There are two ways to separate T and non-T cells:<br />

1. Use the T cell indicator directly from the metadata.<br />
2. Use the expression values of indicator genes, and the clonotype percentage
of the clusters.<br />

## Environment Variables

- `tcell_selector`:
    The expression passed to `tidyseurat::mutate(is_TCell = ...)`
    to indicate whether a cell is a T cell. For example, `Clonotype_Pct > 0.25`
    to indicate cells with clonotype percentage > 25% are T cells.<br />
    If `indicator_genes` is provided, the expression values can also be used
    in the expression. For example, `Clonotype_Pct > 0.25 & CD3E > 0`.<br />
    The expression is first normalized and scaled to have mean 0 and
    standard deviation 1, and then averaged for each cluster.<br />
    If not provided, a kmeans clustering will be performed on the expression
    values of `indicator_genes` and `Clonotype_Pct`, with K=2, and the cluster
    with higher clonotype percentage will be selected as T cells.<br />
- `indicator_genes` *(`list`)*: *Default: `['CD3E']`*. <br />
    A list of indicator genes whose expression values and
    clonotype percentage will be used to determine T cells.<br />
    The markers could be either positive, such as `CD3E`, `CD3D`, `CD3G`, or
    negative, such as `CD19`, `CD14`, `CD68`.<br />

## Examples


### Use T cell indicator directly

If you have a metadata like this:<br />

| id | Clonotype_Pct | seurat_clusters |
|----|---------------|-----------------|
| 1  | 0.1           | 1               |
| 2  | 0.3           | 2               |
| 3  | 0.5           | 3               |

With the configuration below:<br />

```toml
[TCellSelection.envs]
tcell_selector = "Clonotype_Pct > 0.25"
```

The T cells will be selected as:<br />

| id | Clonotype_Pct | seurat_clusters | is_TCell |
|----|---------------|-----------------|----------|
| 1  | 0.1           | 1               | FALSE    |
| 2  | 0.3           | 2               | TRUE     |
| 3  | 0.5           | 3               | TRUE     |

### Use indicator genes

Let's say we set the indicator genes to `["CD3D", "CD3E", "CD3G"]`.<br />

The mean expression values will be calculated for each cluster:<br />

| id | Clonotype_Pct | seurat_clusters | CD3D | CD3E | CD3G |
|----|---------------|-----------------|------|------|------|
| 1  | 0.1           | 1               | 0.1  | 0.0  | 0.1  |
| 2  | 0.3           | 2               | 1.2  | 1.3  | 0.6  |
| 3  | 0.5           | 3               | 1.5  | 0.8  | 0.9  |

Then a kmeans clustering will be performed on the mean expression values of
the indicator genes, together with `Clonotype_Pct`, with K=2.<br />

| id | Clonotype_Pct | seurat_clusters | CD3D | CD3E | CD3G | is_TCell |
|----|---------------|-----------------|------|------|------|----------|
| 1  | 0.1           | 1               | 0.1  | 0.0  | 0.1  | FALSE    |
| 2  | 0.3           | 2               | 1.2  | 1.3  | 0.6  | TRUE     |
| 3  | 0.5           | 3               | 1.5  | 0.8  | 0.9  | TRUE     |

![kmeans](images/TCellSelection-kmeans.png)

The cluster with higher clonoype percentage will be selected as T cells
(`is_TCell = TRUE`), and sent to
[`SeuratClusteringOfTCells`](SeuratClusteringOfTCells.md) for
further clustering and downstream analysis.<br />

