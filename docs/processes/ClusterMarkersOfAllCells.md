# ClusterMarkersOfAllCells

Markers for clusters of all cells.

See also [ClusterMarkers](./ClusterMarkers.md).<br />

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
- `group_by`:
    The column name in metadata to group the cells.<br />
    If only `group_by` is specified, and `ident-1` and `ident-2` are
    not specified, markers will be found for all groups in this column
    in the manner of "group vs rest" comparison.<br />
    `NA` group will be ignored.<br />
    If `None`, `Seurat::Idents(srtobj)` will be used, which is usually
    `"seurat_clusters"` after unsupervised clustering.<br />
- `dbs` *(`list`)*: *Default: `['KEGG_2021_Human', 'MSigDB_Hallmark_2020']`*. <br />
    The dbs to do enrichment analysis for significant
    markers See below for all libraries.<br />
    <https://maayanlab.cloud/Enrichr/#libraries>
- `sigmarkers`: *Default: `p_val_adj < 0.05 & avg_log2FC > 0`*. <br />
    An expression passed to `dplyr::filter()` to filter the
    significant markers for enrichment analysis.<br />
    Available variables are `p_val`, `avg_log2FC`, `pct.1`, `pct.2` and
    `p_val_adj`. For example, `"p_val_adj < 0.05 & abs(avg_log2FC) > 1"`
    to select markers with adjusted p-value < 0.05 and absolute log2
    fold change > 1.<br />
- `enrich_style` *(`choice`)*: *Default: `enrichr`*. <br />
    The style of the enrichment analysis.<br />
    The enrichment analysis will be done by `EnrichIt()` from [`enrichit`](https://user.github.io/enrichit/).<br />
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
- `cache` *(`type=auto`)*: *Default: `/tmp`*. <br />
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
        See <https://user.github.io/biopipen.utils.R/reference/VizDEGs.html>.<br />
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
    - `<more>`:
        Other arguments passed to [`biopipen.utils::VizDEGs()`](https://user.github.io/biopipen.utils.R/reference/VizDEGs.html).<br />
- `allmarker_plots` *(`type=json`)*: *Default: `{'Top 10 markers of all clusters': Diot({'plot_type': 'heatmap'})}`*. <br />
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
        See <https://user.github.io/scplotter/reference/EnrichmentPlot.html>.<br />
- `allenrich_plots` *(`type=json`)*: *Default: `{}`*. <br />
    Cases of the plots to generate for the enrichment analysis.<br />
    The keys are the names of the cases and the values are the dicts inherited from `allenrich_plots_defaults`.<br />
    The cases under `envs.cases` can inherit this options.<br />
- `marker_plots_defaults` *(`ns`)*:
    Default options for the plots to generate for the markers.<br />
    - `plot_type`:
        The type of the plot.<br />
        See <https://user.github.io/biopipen.utils.R/reference/VizDEGs.html>.<br />
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
    - `<more>`:
        Other arguments passed to [`biopipen.utils::VizDEGs()`](https://user.github.io/biopipen.utils.R/reference/VizDEGs.html).<br />
        If `plot_type` is `volcano_pct` or `volcano_log2fc`, they will be passed to
        [`scplotter::VolcanoPlot()`](https://user.github.io/plotthis/reference/VolcanoPlot.html).<br />
- `marker_plots` *(`type=json`)*: *Default: `{'Volcano Plot (diff_pct)': Diot({'plot_type': 'volcano_pct'}), 'Volcano Plot (log2FC)': Diot({'plot_type': 'volcano_log2fc'}), 'Dot Plot': Diot({'plot_type': 'dot'})}`*. <br />
    Cases of the plots to generate for the markers.<br />
    Plot cases. The keys are the names of the cases and the values are the dicts inherited from `marker_plots_defaults`.<br />
    The cases under `envs.cases` can inherit this options.<br />
- `enrich_plots_defaults` *(`ns`)*:
    Default options for the plots to generate for the enrichment analysis.<br />
    - `plot_type`:
        The type of the plot.<br />
        See <https://user.github.io/scplotter/reference/EnrichmentPlot.html>.<br />
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
        See <https://user.github.io/scplotter/reference/EnrichmentPlot.html>.<br />
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
        (<https://user.github.io/plotthis/reference/venndiagram1.html>)
        or `plotthis::UpsetPlot()`
        (<https://user.github.io/plotthis/reference/upsetplot1.html>)
- `overlaps` *(`type=json`)*: *Default: `{}`*. <br />
    Cases for investigating the overlapping of significant markers between different cases or comparisons.<br />
    The keys are the names of the cases and the values are the dicts inherited from `overlaps_defaults`.<br />
    There are two situations that we can perform overlaps:<br />
    1. If `ident-1` is not specified, the overlaps can be performed between different comparisons.<br />
    2. If `each` is specified, the overlaps can be performed between different cases, where in each case, `ident-1` must be specified.<br />

