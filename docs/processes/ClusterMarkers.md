# ClusterMarkers

Markers for clusters of all or selected T cells.

This process is extended from [`MarkersFinder`](https://pwwang.github.io/biopipen/api/biopipen.ns.scrna/#biopipen.ns.scrna.MarkersFinder)
from the [`biopipen`](https://pwwang.github.io/biopipen) package.<br />
`MarkersFinder` is a `pipen` process that wraps the
[`Seurat::FindMarkers()`](https://satijalab.org/seurat/reference/findmarkers)
function, and performs enrichment analysis for the markers found.<br />

The enrichment analysis is done by [`enrichr`](https://maayanlab.cloud/Enrichr/).<br />

/// Note
Since this process is extended from `MarkersFinder`, other environment variables from `MarkersFinder` are also available.<br />
However, they should not be used in this process. Other environment variables are used for more complicated cases for marker finding
(See [`MarkersFinder`](https://pwwang.github.io/biopipen/api/biopipen.ns.scrna/#biopipen.ns.scrna.MarkersFinder) for more details).<br />

If you are using `pipen-board` to run the pipeline
(see [here](../running.md#run-the-pipeline-via-pipen-board) and
[here](../running.md#run-the-pipeline-via-pipen-board-using-docker-image)),
you may see the other environment variables of this process are hidden and readonly.<br />
///

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

