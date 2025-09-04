# PseudoBulkDEG

Pseduo-bulk differential gene expression analysis

This process performs differential gene expression analysis, instead of
on single-cell level, on the pseudo-bulk data, aggregated from the single-cell data.<br />

## Input

- `sobjfile`:
    The seurat object file in RDS or qs/qs2 format.<br />

## Output

- `outdir`: *Default: `{{in.sobjfile | stem}}.pseudobulk_deg`*. <br />
    The output containing the results of the differential gene expression
    analysis.<br />

## Environment Variables

- `mutaters` *(`type=json`)*: *Default: `{}`*. <br />
    Mutaters to mutate the metadata of the
    seurat object. Keys are the new column names and values are the
    expressions to mutate the columns. These new columns can be
    used to define your cases.<br />
    You can also use the clone selectors to select the TCR clones/clusters.<br />
    See <https://pwwang.github.io/scplotter/reference/clone_selectors.html>.<br />
- `each`:
    The column name in metadata to separate the cells into different cases.<br />
    When specified, the case will be expanded to multiple cases for
    each value in the column.<br />
- `subset`:
    An expression in string to subset the cells.<br />
- `aggregate_by`:
    The column names in metadata to aggregate the cells.<br />
- `layer`: *Default: `counts`*. <br />
    The layer to pull and aggregate the data.<br />
- `assay`: *Default: `RNA`*. <br />
    The assay to pull and aggregate the data.<br />
- `error` *(`flag`)*: *Default: `False`*. <br />
    Error out if no/not enough markers are found or no pathways are enriched.<br />
    If `False`, empty results will be returned.<br />
- `group_by`:
    The column name in metadata to group the cells.<br />
- `ident_1`:
    The first identity to compare.<br />
- `ident_2`:
    The second identity to compare.<br />
    If not specified, the rest of the identities will be compared with `ident_1`.<br />
- `paired_by`:
    The column name in metadata to mark the paired samples.<br />
    For example, subject. If specified, the paired test will be performed.<br />
- `dbs` *(`list`)*: *Default: `['KEGG_2021_Human', 'MSigDB_Hallmark_2020']`*. <br />
    The databases to use for enrichment analysis.<br />
    The databases are passed to `biopipen.utils::Enrichr()` to do the
    enrichment analysis. The default databases are `KEGG_2021_Human` and
    `MSigDB_Hallmark_2020`.<br />
    See <https://maayanlab.cloud/Enrichr/#libraries> for the available
    libraries.<br />
- `sigmarkers`: *Default: `p_val_adj < 0.05`*. <br />
    An expression passed to `dplyr::filter()` to filter the
    significant markers for enrichment analysis.<br />
    The default is `p_val_adj < 0.05`.<br />
    If `tool = 'DESeq2'`, the variables that can be used for filtering
    are: `baseMean`, `log2FC`, `lfcSE`, `stat`, `p_val`, `p_val_adj`.<br />
    If `tool = 'edgeR'`, the variables that can be used for filtering
    are: `logCPM`, `log2FC`, `LR`, `p_val`, `p_val_adj`.<br />
- `enrich_style` *(`choice`)*: *Default: `enrichr`*. <br />
    The style of the enrichment analysis.<br />
    - `enrichr`:
        Use `enrichr`-style for the enrichment analysis.<br />
    - `clusterProfiler`:
        Use `clusterProfiler`-style for the enrichment analysis.<br />
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
    - `order_by`: *Default: `desc(abs(log2FC))`*. <br />
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
    - `order_by`: *Default: `desc(abs(log2FC))`*. <br />
        an expression to order the markers, passed by `dplyr::arrange()`.<br />
    - `genes`: *Default: `10`*. <br />
        The number of top genes to show or an expression passed to `dplyr::filter()` to filter the genes.<br />
    - `<more>`:
        Other arguments passed to [`scplotter::FeatureStatPlot()`](https://pwwang.github.io/scplotter/reference/FeatureStatPlot.html).<br />
        If `plot_type` is `volcano_pct` or `volcano_log2fc`, they will be passed to
        [`scplotter::VolcanoPlot()`](https://pwwang.github.io/plotthis/reference/VolcanoPlot.html).<br />
- `marker_plots` *(`type=json`)*: *Default: `{'Volcano Plot': Diot({'plot_type': 'volcano'})}`*. <br />
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
- `tool` *(`choice`)*: *Default: `DESeq2`*. <br />
    The method to use for the differential expression analysis.<br />
    - `DESeq2`:
        Use DESeq2 for the analysis.<br />
    - `edgeR`:
        Use edgeR for the analysis.<br />
- `plots_defaults` *(`ns`)*:
    The default parameters for the plots.<br />
    - `<more>`:
        Parameters passed to `biopipen.utils::VizBulkDEGs()`.<br />
        See: <https://pwwang.github.io/biopipen.utils.R/reference/VizBulkDEGs.html>
- `plots` *(`type=json`)*:
    The parameters for the plots.<br />
    The keys are the names of the plots and the values are the parameters
    for the plots. The parameters will override the defaults in `plots_defaults`.<br />
    If not specified, no plots will be generated.<br />
- `cases` *(`type=json`)*: *Default: `{}`*. <br />
    The cases for the analysis.<br />
    The keys are the names of the cases and the values are the arguments for
    the analysis. The arguments include the ones inherited from `envs`.<br />
    If no cases are specified, a default case will be added with
    the name `DEG Analysis` and the default values specified above.<br />

