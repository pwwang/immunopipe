# TopExpressingGenes

Top expressing genes for clusters of all or selected T/B cells.



This process finds the top expressing genes of clusters of T/B cells, and also
performs the enrichment analysis against the genes.<br />

The enrichment analysis is done by
[`enrichr`](https://maayanlab.cloud/Enrichr/).<br />

/// Note
There are other environment variables also available. However, they should not
be used in this process. Other environment variables are used for more
complicated cases for investigating top genes
(See [`biopipen.ns.scrna.TopExpressingGenes`](https://pwwang.github.io/biopipen/api/biopipen.ns.scrna/#biopipen.ns.scrna.TopExpressingGenes) for more details).<br />

If you are using `pipen-board` to run the pipeline
(see [here](../running.md#run-the-pipeline-via-pipen-board) and
[here](../running.md#run-the-pipeline-via-pipen-board-using-docker-image)),
you may see the other environment variables of this process are hidden and
readonly.<br />
///

## Input

- `srtobj`:
    The seurat object in RDS or qs/qs2 format

## Output

- `outdir`: *Default: `{{in.srtobj | stem}}.top_expressing_genes`*. <br />
    The output directory for the tables and plots

## Environment Variables

- `dbs` *(`list`)*: *Default: `['KEGG_2021_Human', 'MSigDB_Hallmark_2020']`*. <br />
    The dbs to do enrichment analysis for significant
    markers See below for all libraries.<br />
    <https://maayanlab.cloud/Enrichr/#libraries>
- `n` *(`type=int`)*: *Default: `250`*. <br />
    The number of top expressing genes to find.<br />
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
- `subset`:
    An expression to subset the cells for each case.<br />

## SeeAlso

- [TopExpressingGenesOfAllCells](./TopExpressingGenesOfAllCells.md)
- [ClusterMarkers](./ClusterMarkers.md) for examples of enrichment plots

