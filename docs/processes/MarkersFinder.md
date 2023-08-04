# MarkersFinder

This process is extended from [`MarkersFinder`][1] from the [`biopipen`][2] package. `MarkersFinder` is a `pipen` process that wraps the [`Seurat::FindMarkers()`][3] function, and performs enrichment analysis for the markers found.

## Environment variables

- `ncores`: Number of cores to use for parallel computing for some `Seurat` procedures.
    - Used in `future::plan(strategy = "multicore", workers = <ncores>)` to parallelize some Seurat procedures.
    - See also: <https://satijalab.org/seurat/articles/future_vignette.html>
- `mutaters`: The mutaters to mutate the metadata.
    - See also [mutating the metadata](../configurations.md#mutating-the-metadata).
- `ident-1`: The first group of cells to compare
- `ident-2`: The second group of cells to compare
    If not provided, the rest of the cells are used for `ident-2`.
- `group-by`: The column name in metadata to group the cells. If only `group-by` is specified, and `ident-1` and `ident-2` are not specified, markers will be found for all groups in this column in the manner of "group vs rest" comparison. Default: `seurat_clusters`.
    - `NA` group will be ignored.
- `each`: The column name in metadata to separate the cells into different cases.
- `prefix_each` (`flag`): Whether to prefix the `each` column name to the value as the case/section name.
- `dbs` (`list`): The dbs to do enrichment analysis for significant
    markers See below for all libraries.
    <https://maayanlab.cloud/Enrichr/#libraries>
- `sigmarkers`: An expression passed to [`dplyr::filter()`][4] to filter the
    significant markers for enrichment analysis.
    Available variables are `p_val`, `avg_log2FC`, `pct.1`, `pct.2` and
    `p_val_adj`. For example, `"p_val_adj < 0.05 & abs(avg_log2FC) > 1"`
- `section`: The section name for the report.
    Worked only when `each` is not specified and `ident-2` is specified.
    Otherwise, the section name will be constructed from `each` and
    `group-by`.
    If `DEFAULT`, and it's the only section, it not included in the
    case/section names. Default: `DEFAULT`.
- `rest` (`ns`): Rest arguments for [`Seurat::FindMarkers()`][3]. Use `-` to replace `.` in the argument name. For example, use `min-pct` instead of `min.pct`.
    - `<more>`: See https://satijalab.org/seurat/reference/findmarkers
- `cases` (`type=json`): If you have multiple cases, you can specify them
    here. The keys are the names of the cases and the values are the
    above options except `ncores` and `mutaters`. If some options are
    not specified, the default values specified above will be used.
    If no cases are specified, the default case will be added with
    the default values under `envs` with the name `Cluster`.

## Examples

The examples are for more general use of `MarkersFinder`, in order to demonstrate how the final cases are constructed.

Suppose we have a metadata like this:

| id | seurat_clusters | Group |
|----|-----------------|-------|
| 1  | 1               | A     |
| 2  | 1               | A     |
| 3  | 2               | A     |
| 4  | 2               | A     |
| 5  | 3               | B     |
| 6  | 3               | B     |
| 7  | 4               | B     |
| 8  | 4               | B     |

### Default

By default, `group-by` is `seurat_clusters`, and `ident-1` and `ident-2` are not specified. So markers will be found for all clusters in the manner of "cluster vs rest" comparison.

- Cluster
    - 1 (vs 2, 3, 4)
    - 2 (vs 1, 3, 4)
    - 3 (vs 1, 2, 4)
    - 4 (vs 1, 2, 3)

Each case will have the markers and the enrichment analysis for the markers as the results.

### With `each` group

`each` is used to separate the cells into different cases. `group-by` is still `seurat_clusters`.

```toml
[<Proc>.envs]
group-by = "seurat_clusters"
each = "Group"
```

- A:Cluster
    - 1 (vs 2)
    - 2 (vs 1)
- B:Cluster
    - 3 (vs 4)
    - 4 (vs 3)

### With `ident-1` only

`ident-1` is used to specify the first group of cells to compare. Then the rest of the cells in the case are used for `ident-2`.

```toml
[<Proc>.envs]
group-by = "seurat_clusters"
group-1 = "1"
```

- Cluster
    - 1 (vs 2, 3, 4)

### With both `ident-1` and `ident-2`

`ident-1` and `ident-2` are used to specify the two groups of cells to compare.

```toml
[<Proc>.envs]
group-by = "seurat_clusters"
group-1 = "1"
group-2 = "2"
```

- Cluster
    - 1 (vs 2)

### Multie cases

```toml
[<Proc>.envs.cases]
c1_vs_c2 = {group-1 = "1", group-2 = "2"}
c3_vs_c4 = {group-1 = "3", group-2 = "4"}
```

- DEFAULT:c1_vs_c2
    - 1 (vs 2)
- DEFAULT:c3_vs_c4
    - 3 (vs 4)

The `DEFAULT` section name will be ignored in the report. You can specify a section name other than `DEFAULT` for each case to group them in the report.


[1]: https://pwwang.github.io/biopipen/api/biopipen.ns.scrna/#biopipen.ns.scrna.MarkersFinder
[2]: https://pwwang.github.io/biopipen
[3]: https://satijalab.org/seurat/reference/findmarkers
[4]: https://dplyr.tidyverse.org/reference/filter.html