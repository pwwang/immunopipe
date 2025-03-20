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
- `group-by`: *Default: `seurat_clusters`*. <br />
    The column name in metadata to group the cells.<br />
    If only `group-by` is specified, and `ident-1` and `ident-2` are
    not specified, markers will be found for all groups in this column
    in the manner of "group vs rest" comparison.<br />
    `NA` group will be ignored.<br />
- `prefix_group` *(`flag`)*: *Default: `True`*. <br />
    When neither `ident-1` nor `ident-2` is specified,
    should we prefix the group name to the section name?<br />
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

