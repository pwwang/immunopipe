# TOrBCellSelection

Separate T and non-T cells and select T cells; or separate B and non-B cells and select B cells.

If all of your cells are T/B cells, do not set any configurations for this process.<br />

In such a case, [`SeuratClusteringOfAllCells`](SeuratClusteringOfAllCells.md) should
not be used, and [`SeuratClustering`](SeuratClustering.md) will be clustering all
of the cells, which are all T/B cells.<br />

There are two ways to separate T and non-T cells; or B and non-B cells:<br />

1. Use the an expression indicator directly from the metadata.<br />
2. Use the expression values of indicator genes, and the clonotype percentage
of the clusters.<br />

You can also use indicator gene expression values only to select T/B cells by
setting `envs.ignore_vdj` to true.<br />

## Input

- `srtobj`:
    Seurat object file in RDS/qs2
- `immdata`:
    Immune repertoire data file in RDS/qs2

## Output

- `outfile`: *Default: `{{in.srtobj | stem}}.qs`*. <br />
    Seurat object file in qs2 format
- `outdir`: *Default: `details`*. <br />
    Output directory with details

## Environment Variables

- `ignore_vdj` *(`flag`)*: *Default: `False`*. <br />
    Ignore VDJ information for T/B cell selection.<br />
    Use only the expression values of indicator genes if True.<br />
    In this case, the `Clonotype_Pct` column does not exist in the metadata.<br />
    If you want to use `k-means` to select T/B cells, you must have more than
    1 indicator gene, and the first indicator gene in `envs.indicator_genes`
    must be a positive marker, which will be used to select the cluster with
    higher expression values as T/B cells.<br />
- `selector`:
    The expression passed to `tidyseurat::mutate(is_TCell = ...)`
    to indicate whether a cell is a T cell. For example, `Clonotype_Pct > 0.25`
    to indicate cells with clonotype percentage > 25% are T cells.<br />
    If `indicator_genes` is provided, the expression values can also be used
    in the expression. For example, `Clonotype_Pct > 0.25 & CD3E > 0`.<br />
    If `selector` is not provided, a kmeans clustering will be performed
    on the expression values of `indicator_genes` and `Clonotype_Pct`,
    with K=2, and the cluster with higher clonotype percentage will be selected
    as T/B cells.<br />
- `indicator_genes` *(`list`)*: *Default: `['CD3E']`*. <br />
    A list of indicator genes whose expression values and
    clonotype percentage will be used to determine T/B cells.<br />
    The markers could be either positive, such as `CD3E`, `CD3D`, `CD3G`, or
    negative, such as `CD19`, `CD14`, `CD68`, for T cells. For B cells,
    markers such as `CD19`, `MS4A1` (CD20), `CD79A`, `CD79B` could be used.<br />

- `kmeans` *(`type=json`)*: *Default: `{'nstart': 25}`*. <br />
    The parameters for `kmeans` clustering.<br />
    Other arguments for [`stats::kmeans`](https://rdrr.io/r/stats/kmeans.html)
    can be provided here. If there are dots in the argument names, replace them
    with `-`.<br />

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
[TOrBCellSelection.envs]
selector = "Clonotype_Pct > 0.25"
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

The cluster with higher clonoype percentage will be selected as T/B cells
(`is_selected = TRUE`), and sent to
[`SeuratClustering`](SeuratClustering.md) for
further clustering and downstream analysis.<br />

