# MarkersFinder

Find markers between different groups of cells

`MarkersFinder` is a process that wraps the
[`Seurat::FindMarkers()`](https://satijalab.org/seurat/reference/findmarkers)
function, and performs enrichment analysis for the markers found.<br />

## Environment Variables

- `ncores` *(`type=int`)*: *Default: `1`*. <br />
    Number of cores to use for parallel computing for some `Seurat` procedures.<br />
    * Used in `future::plan(strategy = "multicore", workers = <ncores>)` to parallelize some Seurat procedures.<br />
    * See also: <https://satijalab.org/seurat/articles/future_vignette.html>
- `mutaters` *(`type=json`)*: *Default: `{}`*. <br />
    The mutaters to mutate the metadata There are also also 4 helper functions, `expanded`, `collapsed`, `emerged` and `vanished`, which can be used to identify the expanded/collpased/emerged/vanished groups (i.e. TCR clones).<br />
    See also <../configurations/#mutater-helpers>.<br />
    For example, you can use `{"Patient1_Tumor_Collapsed_Clones": "expanded(., Source, 'Tumor', subset = Patent == 'Patient1', uniq = FALSE)"}` to create a new column in metadata named `Patient1_Tumor_Collapsed_Clones` with the collapsed clones in the tumor sample (compared to the normal sample) of patient 1.<br />
    The values in this columns for other clones will be `NA`.<br />
    Those functions take following arguments:<br />
    * `df`: The metadata data frame. You can use the `.` to refer to it.<br />
    * `group.by`: The column name in metadata to group the cells.<br />
    * `idents`: The first group or both groups of cells to compare (value in `group.by` column). If only the first group is given, the rest of the cells (with non-NA in `group.by` column) will be used as the second group.<br />
    * `subset`: An expression to subset the cells, will be passed to `dplyr::filter()`. Default is `TRUE` (no filtering).<br />
    * `each`: A column name (without quotes) in metadata to split the cells.<br />
    Each comparison will be done for each value in this column (typically each patient or subject).<br />
    * `id`: The column name in metadata for the group ids (i.e. `CDR3.aa`).<br />
    * `compare`: Either a (numeric) column name (i.e. `Clones`) in metadata to compare between groups, or `.n` to compare the number of cells in each group.<br />
    If numeric column is given, the values should be the same for all cells in the same group.<br />
    This will not be checked (only the first value is used).<br />
    It is helpful to use `Clones` to use the raw clone size from TCR data, in case the cells are not completely mapped to RNA data.<br />
    Also if you have `subset` set or `NA`s in `group.by` column, you should use `.n` to compare the number of cells in each group.<br />
    * `uniq`: Whether to return unique ids or not. Default is `TRUE`. If `FALSE`, you can mutate the meta data frame with the returned ids. For example, `df |> mutate(expanded = expanded(...))`.<br />
    * `debug`: Return the data frame with intermediate columns instead of the ids. Default is `FALSE`.<br />
    * `order`: The expression passed to `dplyr::arrange()` to order intermediate dataframe and get the ids in order accordingly.<br />
    The intermediate dataframe includes the following columns:<br />
    * `<id>`: The ids of clones (i.e. `CDR3.aa`).<br />
    * `<each>`: The values in `each` column.<br />
    * `ident_1`: The size of clones in the first group.<br />
    * `ident_2`: The size of clones in the second group.<br />
    * `.diff`: The difference between the sizes of clones in the first and second groups.<br />
    * `.sum`: The sum of the sizes of clones in the first and second groups.<br />
    * `.predicate`: Showing whether the clone is expanded/collapsed/emerged/vanished.<br />
    * `include_emerged`: Whether to include the emerged group for `expanded` (only works for `expanded`). Default is `FALSE`.<br />
    * `include_vanished`: Whether to include the vanished group for `collapsed` (only works for `collapsed`). Default is `FALSE`..<br />
    See also
    [mutating the metadata](../configurations.md#mutating-the-metadata).<br />
- `ident-1`:
    The first group of cells to compare
- `ident-2`:
    The second group of cells to compare
    If not provided, the rest of the cells are used for `ident-2`.<br />
- `group-by`: *Default: `seurat_clusters`*. <br />
    The column name in metadata to group the cells.<br />
    If only `group-by` is specified, and `ident-1` and `ident-2` are
    not specified, markers will be found for all groups in this column
    in the manner of "group vs rest" comparison.<br />
    `NA` group will be ignored.<br />
- `each`:
    The column name in metadata to separate the cells into different
    cases.<br />
- `prefix_each` *(`flag`)*: *Default: `True`*. <br />
    Whether to prefix the `each` column name to the
    value as the case/section name.<br />
- `dbs` *(`list`)*: *Default: `['KEGG_2021_Human', 'MSigDB_Hallmark_2020']`*. <br />
    The dbs to do enrichment analysis for significant
    markers See below for all libraries.<br />
    <https://maayanlab.cloud/Enrichr/#libraries>
- `sigmarkers`: *Default: `p_val_adj < 0.05`*. <br />
    An expression passed to `dplyr::filter()` to filter the
    significant markers for enrichment analysis.<br />
    Available variables are `p_val`, `avg_log2FC`, `pct.1`, `pct.2` and
    `p_val_adj`. For example, `"p_val_adj < 0.05 & abs(avg_log2FC) > 1"`
    to select markers with adjusted p-value < 0.05 and absolute log2
    fold change > 1.<br />
- `assay`:
    The assay to use.<br />
- `volcano_genes` *(`type=auto`)*: *Default: `True`*. <br />
    The genes to label in the volcano plot if they are
    significant markers.<br />
    If `True`, all significant markers will be labeled. If `False`, no
    genes will be labeled. Otherwise, specify the genes to label.<br />
    It could be either a string with comma separated genes, or a list
    of genes.<br />
- `section`: *Default: `DEFAULT`*. <br />
    The section name for the report. It must not contain colon (`:`).<br />
    Ignored when `each` is not specified and `ident-1` is specified.<br />
    When neither `each` nor `ident-1` is specified, case name will be used
    as section name.<br />
    If `each` is specified, the section name will be constructed from
    `each` and case name.<br />
- `subset`:
    An expression to subset the cells for each case.<br />
- `use_presto`: *Default: `False`*. <br />
    Whether to use [`presto::wilcoxauc`](https://rdrr.io/github/immunogenomics/presto/man/wilcoxauc.html)
    to find markers.<br />
    [`presto`](https://github.com/immunogenomics/presto) is a package performs
    fast Wilcoxon rank sum test and auROC analysis.<br />
- `rest` *(`ns`)*:
    Rest arguments for `Seurat::FindMarkers()`.<br />
    Use `-` to replace `.` in the argument name. For example,
    use `min-pct` instead of `min.pct`.<br />
    This only works when `use_presto` is `False`.<br />
    - `<more>`:
        See <https://satijalab.org/seurat/reference/findmarkers>
- `dotplot` *(`ns`)*:
    Arguments for `Seurat::DotPlot()`.<br />
    Use `-` to replace `.` in the argument name. For example,
    use `group-bar` instead of `group.bar`.<br />
    Note that `object`, `features`, and `group-by` are already specified
    by this process. So you don't need to specify them here.<br />
    - `devpars` *(`ns`)*:
        The device parameters for the plots.<br />
        - `res` *(`type=int`)*:
            The resolution of the plots.<br />
        - `height` *(`type=int`)*:
            The height of the plots.<br />
        - `width` *(`type=int`)*:
            The width of the plots.<br />
    - `<more>`:
        See <https://satijalab.org/seurat/reference/doheatmap>
- `cases` *(`type=json`)*: *Default: `{}`*. <br />
    If you have multiple cases, you can specify them
    here. The keys are the names of the cases and the values are the
    above options except `ncores` and `mutaters`. If some options are
    not specified, the default values specified above will be used.<br />
    If no cases are specified, the default case will be added with
    the default values under `envs` with the name `DEFAULT`.<br />
- `overlap` *(`list`)*: *Default: `[]`*. <br />
    The sections to do overlap analysis.<br />

## Examples

The examples are for more general use of `MarkersFinder`, in order to
demonstrate how the final cases are constructed.<br />

Suppose we have a metadata like this:<br />

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

By default, `group-by` is `seurat_clusters`, and `ident-1` and `ident-2`
are not specified. So markers will be found for all clusters in the manner
of "cluster vs rest" comparison.<br />

- Cluster
    - 1 (vs 2, 3, 4)
    - 2 (vs 1, 3, 4)
    - 3 (vs 1, 2, 4)
    - 4 (vs 1, 2, 3)

Each case will have the markers and the enrichment analysis for the
markers as the results.<br />

### With `each` group

`each` is used to separate the cells into different cases. `group-by`
is still `seurat_clusters`.<br />

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

`ident-1` is used to specify the first group of cells to compare.<br />
Then the rest of the cells in the case are used for `ident-2`.<br />

```toml
[<Proc>.envs]
group-by = "seurat_clusters"
ident-1 = "1"
```

- Cluster
    - 1 (vs 2, 3, 4)

### With both `ident-1` and `ident-2`

`ident-1` and `ident-2` are used to specify the two groups of cells to
compare.<br />

```toml
[<Proc>.envs]
group-by = "seurat_clusters"
ident-1 = "1"
ident-2 = "2"
```

- Cluster
    - 1 (vs 2)

### Multiple cases

```toml
[<Proc>.envs.cases]
c1_vs_c2 = {ident-1 = "1", ident-2 = "2"}
c3_vs_c4 = {ident-1 = "3", ident-2 = "4"}
```

- DEFAULT:c1_vs_c2
    - 1 (vs 2)
- DEFAULT:c3_vs_c4
    - 3 (vs 4)

The `DEFAULT` section name will be ignored in the report. You can specify
a section name other than `DEFAULT` for each case to group them
in the report.<br />

