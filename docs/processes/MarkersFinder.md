# MarkersFinder

Find markers between different groups of cells

`MarkersFinder` is a process that wraps the
[`Seurat::FindMarkers()`](https://satijalab.org/seurat/reference/findmarkers)
function, and performs enrichment analysis for the markers found.<br />

## Input

- `srtobj`:
    The seurat object loaded by `SeuratPreparing`
    If you have your `Seurat` object prepared by yourself, you can also
    use it here, but you should make sure that the object has been processed
    by `PrepSCTFindMarkers` if data is not normalized using `SCTransform`.<br />

## Output

- `outdir`: *Default: `{{in.srtobj | stem0}}.markers`*. <br />
    The output directory for the markers and plots

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
    * `include_vanished`: Whether to include the vanished group for `collapsed` (only works for `collapsed`). Default is `FALSE`.<br />
    You can also use `top()` to get the top clones (i.e. the clones with the largest size) in each group.<br />
    For example, you can use `{"Patient1_Top10_Clones": "top(subset = Patent == 'Patient1', uniq = FALSE)"}` to create a new column in metadata named `Patient1_Top10_Clones`.<br />
    The values in this columns for other clones will be `NA`.<br />
    This function takes following arguments:<br />
    * `df`: The metadata data frame. You can use the `.` to refer to it.<br />
    * `id`: The column name in metadata for the group ids (i.e. `CDR3.aa`).<br />
    * `n`: The number of top clones to return. Default is `10`.<br />
    If n < 1, it will be treated as the percentage of the size of the group.<br />
    Specify `0` to get all clones.<br />
    * `compare`: Either a (numeric) column name (i.e. `Clones`) in metadata to compare between groups, or `.n` to compare the number of cells in each group.<br />
    If numeric column is given, the values should be the same for all cells in the same group.<br />
    This will not be checked (only the first value is used).<br />
    It is helpful to use `Clones` to use the raw clone size from TCR data, in case the cells are not completely mapped to RNA data.<br />
    Also if you have `subset` set or `NA`s in `group.by` column, you should use `.n` to compare the number of cells in each group.<br />
    * `subset`: An expression to subset the cells, will be passed to `dplyr::filter()`. Default is `TRUE` (no filtering).<br />
    * `each`: A column name (without quotes) in metadata to split the cells.<br />
    Each comparison will be done for each value in this column (typically each patient or subject).<br />
    * `uniq`: Whether to return unique ids or not. Default is `TRUE`. If `FALSE`, you can mutate the meta data frame with the returned ids. For example, `df |> mutate(expanded = expanded(...))`.<br />
    * `debug`: Return the data frame with intermediate columns instead of the ids. Default is `FALSE`.<br />
    * `with_ties`: Whether to include ties (i.e. clones with the same size as the last clone) or not. Default is `FALSE`..<br />
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
- `prefix_group` *(`flag`)*: *Default: `True`*. <br />
    When neither `ident-1` nor `ident-2` is specified,
    should we prefix the group name to the section name?<br />
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
- `error` *(`flag`)*: *Default: `True`*. <br />
    Error out if no/not enough markers are found or no pathways are enriched.<br />
    If `False`, empty results will be returned.<br />
- `site`: *Default: `Enrichr`*. <br />
    The site to use for the `enrichR` enrichment analysis.<br />
- `subset`:
    An expression to subset the cells for each case.<br />
- `cache` *(`type=auto`)*: *Default: `/tmp/user`*. <br />
    Where to cache to `FindAllMarkers` results.<br />
    If `True`, cache to `outdir` of the job. If `False`, don't cache.<br />
    Otherwise, specify the directory to cache to.<br />
- `rest` *(`ns`)*:
    Rest arguments for `Seurat::FindMarkers()`.<br />
    Use `-` to replace `.` in the argument name. For example,
    use `min-pct` instead of `min.pct`.<br />
    This only works when `use_presto` is `False`.<br />
    - `<more>`:
        See <https://satijalab.org/seurat/reference/findmarkers>
- `allmarker_plots_defaults` *(`ns`)*:
    Default options for the plots for all markers when `ident-1` is not specified.<br />
    - `plot_type`:
        The type of the plot.<br />
        See <https://pwwang.github.io/scplotter/reference/FeatureStatPlot.html>.<br />
        Available types are `violin`, `box`, `bar`, `ridge`, `dim`, `heatmap` and `dot`.<br />
    - `more_formats` *(`list`)*: *Default: `[]`*. <br />
        The extra formats to save the plot in.<br />
    - `save_code` *(`flag`)*: *Default: `False`*. <br />
        Whether to save the code to generate the plot.<br />
    - `devpars` *(`ns`)*:
        The device parameters for the plots.<br />
        - `res` *(`type=int`)*: *Default: `100`*. <br />
            The resolution of the plots.<br />
        - `height` *(`type=int`)*:
            The height of the plots.<br />
        - `width` *(`type=int`)*:
            The width of the plots.<br />
    - `order_by`: *Default: `desc(abs(avg_log2FC))`*. <br />
        an expression to order the markers, passed by `dplyr::arrange()`.<br />
    - `genes`: *Default: `10`*. <br />
        The number of top genes to show or an expression passed to `dplyr::filter()` to filter the genes.<br />
    - `<more>`:
        Other arguments passed to [`scplotter::FeatureStatPlot()`](https://pwwang.github.io/scplotter/reference/FeatureStatPlot.html).<br />
- `allmarker_plots` *(`type=json`)*: *Default: `{}`*. <br />
    All marker plot cases.<br />
    The keys are the names of the cases and the values are the dicts inherited from `allmarker_plots_defaults`.<br />
- `marker_plots_defaults` *(`ns`)*:
    Default options for the plots to generate for the markers.<br />
    - `plot_type`:
        The type of the plot.<br />
        See <https://pwwang.github.io/scplotter/reference/FeatureStatPlot.html>.<br />
        Available types are `violin`, `box`, `bar`, `ridge`, `dim`, `heatmap` and `dot`.<br />
        There are two additional types available - `volcano_pct` and `volcano_log2fc`.<br />
    - `more_formats` *(`list`)*: *Default: `[]`*. <br />
        The extra formats to save the plot in.<br />
    - `save_code` *(`flag`)*: *Default: `False`*. <br />
        Whether to save the code to generate the plot.<br />
    - `devpars` *(`ns`)*:
        The device parameters for the plots.<br />
        - `res` *(`type=int`)*: *Default: `100`*. <br />
            The resolution of the plots.<br />
        - `height` *(`type=int`)*:
            The height of the plots.<br />
        - `width` *(`type=int`)*:
            The width of the plots.<br />
    - `order_by`: *Default: `desc(abs(avg_log2FC))`*. <br />
        an expression to order the markers, passed by `dplyr::arrange()`.<br />
    - `genes`: *Default: `10`*. <br />
        The number of top genes to show or an expression passed to `dplyr::filter()` to filter the genes.<br />
    - `<more>`:
        Other arguments passed to [`scplotter::FeatureStatPlot()`](https://pwwang.github.io/scplotter/reference/FeatureStatPlot.html).<br />
        If `plot_type` is `volcano_pct` or `volcano_log2fc`, they will be passed to
        [`scplotter::VolcanoPlot()`](https://pwwang.github.io/plotthis/reference/VolcanoPlot.html).<br />
- `marker_plots` *(`type=json`)*: *Default: `{'Volcano Plot (diff_pct)': Diot({'plot_type': 'volcano_pct'}), 'Volcano Plot (log2FC)': Diot({'plot_type': 'volcano_log2fc'}), 'Dot Plot': Diot({'plot_type': 'dot'})}`*. <br />
    Cases of the plots to generate for the markers.<br />
    Plot cases. The keys are the names of the cases and the values are the dicts inherited from `marker_plots_defaults`.<br />
- `enrich_plots_defaults` *(`ns`)*:
    Default options for the plots to generate for the enrichment analysis.<br />
    - `plot_type`:
        The type of the plot.<br />
        See <https://pwwang.github.io/scplotter/reference/EnrichmentPlot.html>.<br />
        Available types are `bar`, `dot`, `lollipop`, `network`, `enrichmap` and `wordcloud`.<br />
    - `more_formats` *(`list`)*: *Default: `[]`*. <br />
        The extra formats to save the plot in.<br />
    - `save_code` *(`flag`)*: *Default: `False`*. <br />
        Whether to save the code to generate the plot.<br />
    - `devpars` *(`ns`)*:
        The device parameters for the plots.<br />
        - `res` *(`type=int`)*: *Default: `100`*. <br />
            The resolution of the plots.<br />
        - `height` *(`type=int`)*:
            The height of the plots.<br />
        - `width` *(`type=int`)*:
            The width of the plots.<br />
    - `<more>`:
        See <https://pwwang.github.io/scplotter/reference/EnrichmentPlot.htmll>.<br />
- `enrich_plots` *(`type=json`)*: *Default: `{'Bar Plot': Diot({'plot_type': 'bar', 'ncol': 1, 'top_term': 10})}`*. <br />
    Cases of the plots to generate for the enrichment analysis.<br />
    The keys are the names of the cases and the values are the dicts inherited from `enrich_plots_defaults`.<br />
- `cases` *(`type=json`)*: *Default: `{}`*. <br />
    If you have multiple cases for marker discovery, you can specify them
    here. The keys are the names of the cases and the values are the above options. If some options are
    not specified, the default values specified above (under `envs`) will be used.<br />
    If no cases are specified, the default case will be added with the default values under `envs` with the name `DEFAULT`.<br />
    If you want to put some cases under the same section in the report, you can specify the section name in the case name
    as a prefix separated by `::`. For example, `section1::case1` and `section1::case2` will be put `case1` and `case2`
    under the section `section1`.<br />
- `overlaps_defaults` *(`ns`)*:
    Default options for investigating the overlapping of significant markers between different cases.<br />
    - `cases` *(`list`)*: *Default: `[]`*. <br />
        The cases to do the overlapping analysis, including the prefix section name.<br />
        The case must have `ident-1` specified. When `each` is specified, the case will be expanded.<br />
        For example, `case1` with `each = "group"`, where `group` has `g1` and `g2`, will be expanded to
        `case1::g1` and `case1::g2`, or `case1::group - g1` and `case1::group - g2` if `prefix_each` is `True`.<br />
        There must be at least 2 cases to do the overlapping analysis.<br />
    - `sigmarkers`:
        The expression to filter the significant markers for each case.<br />
        If not provided, `envs.sigmarkers` will be used.<br />
    - `venn` *(`ns`)*:
        The options for the Venn diagram.<br />
        - `enabled` *(`flag`)*: *Default: `auto`*. <br />
            Whether to enable the Venn diagram.<br />
            Default is "auto", which means enabled when there are no more than 5 cases.<br />
        - `more_formats` *(`list`)*: *Default: `[]`*. <br />
            The extra formats to save the plot in.<br />
        - `save_code` *(`flag`)*: *Default: `False`*. <br />
            Whether to save the code to generate the plot.<br />
        - `devpars` *(`ns`)*:
            The device parameters for the plots.<br />
            - `res` *(`type=int`)*: *Default: `100`*. <br />
                The resolution of the plots.<br />
            - `height` *(`type=int`)*:
                The height of the plots.<br />
            - `width` *(`type=int`)*:
                The width of the plots.<br />
        - `<more>`:
            More arguments pased to `plotthis::VennDiagram()`.<br />
            https://pwwang.github.io/plotthis/reference/venndiagram1.html
    - `upset` *(`ns`)*:
        The options for the UpSet plot.<br />
        - `enabled` *(`flag`)*: *Default: `True`*. <br />
            Whether to enable the UpSet plot.<br />
        - `more_formats` *(`list`)*: *Default: `[]`*. <br />
            The extra formats to save the plot in.<br />
        - `save_code` *(`flag`)*: *Default: `False`*. <br />
            Whether to save the code to generate the plot.<br />
        - `devpars` *(`ns`)*:
            The device parameters for the plots.<br />
            - `res` *(`type=int`)*: *Default: `100`*. <br />
                The resolution of the plots.<br />
            - `height` *(`type=int`)*:
                The height of the plots.<br />
            - `width` *(`type=int`)*:
                The width of the plots.<br />
        - `<more>`:
            More arguments pased to `plotthis::UpsetPlot()`.<br />
            https://pwwang.github.io/plotthis/reference/upsetplot1.html
- `overlaps` *(`type=json`)*: *Default: `{}`*. <br />
    Cases for investigating the overlapping of significant markers between different cases.<br />
    The keys are the names of the cases and the values are the dicts inherited from `overlaps_defaults`.<br />

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

