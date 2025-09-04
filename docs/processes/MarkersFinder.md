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
    The mutaters to mutate the metadata.<br />
    You can also use the clone selectors to select the TCR clones/clusters.<br />
    See <https://pwwang.github.io/scplotter/reference/clone_selectors.html>..<br />
    See also
    [mutating the metadata](../configurations.md#mutating-the-metadata).<br />
- `group_by`:
    The column name in metadata to group the cells.<br />
    If only `group_by` is specified, and `ident-1` and `ident-2` are
    not specified, markers will be found for all groups in this column
    in the manner of "group vs rest" comparison.<br />
    `NA` group will be ignored.<br />
    If `None`, `Seurat::Idents(srtobj)` will be used, which is usually
    `"seurat_clusters"` after unsupervised clustering.<br />
- `ident_1`:
    The first group of cells to compare
    When this is empty, the comparisons will be expanded to each group v.s. the rest of the cells in `group_by`.<br />
- `ident_2`:
    The second group of cells to compare
    If not provided, the rest of the cells are used for `ident-2`.<br />
- `each`:
    The column name in metadata to separate the cells into different
    cases.<br />
    When this is specified, the case will be expanded for each value of
    the column in metadata. For example, when you have `envs.cases."Cluster Markers".each = "Sample"`,
    then the case will be expanded as `envs.cases."Cluster Markers - Sample1"`, `envs.cases."Cluster Markers - Sample2"`, etc.<br />
    You can specify `allmarker_plots` and `overlaps` to plot the markers for all cases in the same plot and plot the overlaps of the markers
    between different cases by values in this column.<br />
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
- `enrich_style` *(`choice`)*: *Default: `enrichr`*. <br />
    The style of the enrichment analysis.<br />
    The enrichment analysis will be done by `EnrichIt()` from [`enrichit`](https://pwwang.github.io/enrichit/).<br />
    Two styles are available:<br />
    - `enrichr`:
        `enrichr` style enrichment analysis (fisher's exact test will be used).<br />
    - `clusterprofiler`:
        `clusterProfiler` style enrichment analysis (hypergeometric test will be used).<br />
    - `clusterProfiler`:
        alias for `clusterprofiler`
- `assay`:
    The assay to use.<br />
- `error` *(`flag`)*: *Default: `False`*. <br />
    Error out if no/not enough markers are found or no pathways are enriched.<br />
    If `False`, empty results will be returned.<br />
- `subset`:
    An expression to subset the cells for each case.<br />
- `cache` *(`type=auto`)*: *Default: `/tmp/m161047`*. <br />
    Where to cache the results.<br />
    If `True`, cache to `outdir` of the job. If `False`, don't cache.<br />
    Otherwise, specify the directory to cache to.<br />
- `rest` *(`ns`)*:
    Rest arguments for `Seurat::FindMarkers()`.<br />
    Use `-` to replace `.` in the argument name. For example,
    use `min-pct` instead of `min.pct`.<br />
    - `<more>`:
        See <https://satijalab.org/seurat/reference/findmarkers>
- `allmarker_plots_defaults` *(`ns`)*:
    Default options for the plots for all markers when `ident-1` is not specified.<br />
    - `plot_type`:
        The type of the plot.<br />
        See <https://pwwang.github.io/scplotter/reference/FeatureStatPlot.html>.<br />
        Available types are `violin`, `box`, `bar`, `ridge`, `dim`, `heatmap` and `dot`.<br />
    - `more_formats` *(`type=list`)*: *Default: `[]`*. <br />
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
- `allenrich_plots_defaults` *(`ns`)*:
    Default options for the plots to generate for the enrichment analysis.<br />
    - `plot_type`: *Default: `heatmap`*. <br />
        The type of the plot.<br />
    - `devpars` *(`ns`)*:
        The device parameters for the plots.<br />
        - `res` *(`type=int`)*: *Default: `100`*. <br />
            The resolution of the plots.<br />
        - `height` *(`type=int`)*:
            The height of the plots.<br />
        - `width` *(`type=int`)*:
            The width of the plots.<br />
    - `<more>`:
        See <https://pwwang.github.io/scplotter/reference/EnrichmentPlot.html>.<br />
- `allenrich_plots` *(`type=json`)*: *Default: `{}`*. <br />
    Cases of the plots to generate for the enrichment analysis.<br />
    The keys are the names of the cases and the values are the dicts inherited from `allenrich_plots_defaults`.<br />
    The cases under `envs.cases` can inherit this options.<br />
- `marker_plots_defaults` *(`ns`)*:
    Default options for the plots to generate for the markers.<br />
    - `plot_type`:
        The type of the plot.<br />
        See <https://pwwang.github.io/scplotter/reference/FeatureStatPlot.html>.<br />
        Available types are `violin`, `box`, `bar`, `ridge`, `dim`, `heatmap` and `dot`.<br />
        There are two additional types available - `volcano_pct` and `volcano_log2fc`.<br />
    - `more_formats` *(`type=list`)*: *Default: `[]`*. <br />
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
    The cases under `envs.cases` can inherit this options.<br />
- `enrich_plots_defaults` *(`ns`)*:
    Default options for the plots to generate for the enrichment analysis.<br />
    - `plot_type`:
        The type of the plot.<br />
        See <https://pwwang.github.io/scplotter/reference/EnrichmentPlot.html>.<br />
        Available types are `bar`, `dot`, `lollipop`, `network`, `enrichmap` and `wordcloud`.<br />
    - `more_formats` *(`type=list`)*: *Default: `[]`*. <br />
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
    The cases under `envs.cases` can inherit this options.<br />
- `overlaps_defaults` *(`ns`)*:
    Default options for investigating the overlapping of significant markers between different cases or comparisons.<br />
    This means either `ident-1` should be empty, so that they can be expanded to multiple comparisons.<br />
    - `sigmarkers`:
        The expression to filter the significant markers for each case.<br />
        If not provided, `envs.sigmarkers` will be used.<br />
    - `plot_type` *(`choice`)*: *Default: `venn`*. <br />
        The type of the plot to generate for the overlaps.<br />
        - `venn`:
            Use `plotthis::VennDiagram()`.<br />
        - `upset`:
            Use `plotthis::UpsetPlot()`.<br />
    - `more_formats` *(`type=list`)*: *Default: `[]`*. <br />
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
        More arguments pased to `plotthis::VennDiagram()`
        (<https://pwwang.github.io/plotthis/reference/venndiagram1.html>)
        or `plotthis::UpsetPlot()`
        (<https://pwwang.github.io/plotthis/reference/upsetplot1.html>)
- `overlaps` *(`type=json`)*: *Default: `{}`*. <br />
    Cases for investigating the overlapping of significant markers between different cases or comparisons.<br />
    The keys are the names of the cases and the values are the dicts inherited from `overlaps_defaults`.<br />
    There are two situations that we can perform overlaps:<br />
    1. If `ident-1` is not specified, the overlaps can be performed between different comparisons.<br />
    2. If `each` is specified, the overlaps can be performed between different cases, where in each case, `ident-1` must be specified.<br />
- `cases` *(`type=json`)*: *Default: `{}`*. <br />
    If you have multiple cases for marker discovery, you can specify them
    here. The keys are the names of the cases and the values are the above options. If some options are
    not specified, the default values specified above (under `envs`) will be used.<br />
    If no cases are specified, the default case will be added with the default values under `envs` with the name `Marker Discovery`.<br />

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

By default, `group_by` is `seurat_clusters`, and `ident_1` and `ident_2`
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

`each` is used to separate the cells into different cases. `group_by`
is still `seurat_clusters`.<br />

```toml
[<Proc>.envs]
group_by = "seurat_clusters"
each = "Group"
```

- A:Cluster
    - 1 (vs 2)
    - 2 (vs 1)
- B:Cluster
    - 3 (vs 4)
    - 4 (vs 3)

### With `ident_1` only

`ident_1` is used to specify the first group of cells to compare.<br />
Then the rest of the cells in the case are used for `ident_2`.<br />

```toml
[<Proc>.envs]
group_by = "seurat_clusters"
ident_1 = "1"
```

- Cluster
    - 1 (vs 2, 3, 4)

### With both `ident_1` and `ident_2`

`ident_1` and `ident_2` are used to specify the two groups of cells to
compare.<br />

```toml
[<Proc>.envs]
group_by = "seurat_clusters"
ident_1 = "1"
ident_2 = "2"
```

- Cluster
    - 1 (vs 2)

### Multiple cases

```toml
[<Proc>.envs.cases]
c1_vs_c2 = {ident_1 = "1", ident_2 = "2"}
c3_vs_c4 = {ident_1 = "3", ident_2 = "4"}
```

- DEFAULT:c1_vs_c2
    - 1 (vs 2)
- DEFAULT:c3_vs_c4
    - 3 (vs 4)

The `DEFAULT` section name will be ignored in the report. You can specify
a section name other than `DEFAULT` for each case to group them
in the report.<br />

