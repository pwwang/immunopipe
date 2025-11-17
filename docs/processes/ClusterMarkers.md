# ClusterMarkers

Markers for clusters of all or selected T/B cells.

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
    The dbs to do enrichment analysis for significant markers.<br />
    You can use built-in dbs in `enrichit`, or provide your own gmt files.<br />
    See also <https://pwwang.github.io/enrichit/reference/FetchGMT.html>.<br />
    The built-in dbs include:<br />
    * "BioCarta" or "BioCarta_2016"
    * "GO_Biological_Process" or "GO_Biological_Process_2025"
    * "GO_Cellular_Component" or "GO_Cellular_Component_2025"
    * "GO_Molecular_Function" or "GO_Molecular_Function_2025"
    * "KEGG", "KEGG_Human", "KEGG_2021", or "KEGG_2021_Human"
    * "Hallmark", "MSigDB_Hallmark", or "MSigDB_Hallmark_2020"
    * "Reactome", "Reactome_Pathways", or "Reactome_Pathways_2024"
    * "WikiPathways", "WikiPathways_2024", "WikiPathways_Human", or "WikiPathways_2024_Human"
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
    Default options for the plots for all markers when `ident_1` is not specified.<br />
    - `plot_type`:
        The type of the plot.<br />
        See <https://pwwang.github.io/biopipen.utils.R/reference/VizDEGs.html>.<br />
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
        Other arguments passed to [`biopipen.utils::VizDEGs()`](https://pwwang.github.io/biopipen.utils.R/reference/VizDEGs.html).<br />
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
        See <https://pwwang.github.io/scplotter/reference/EnrichmentPlot.html>.<br />
- `allenrich_plots` *(`type=json`)*: *Default: `{}`*. <br />
    Cases of the plots to generate for the enrichment analysis.<br />
    The keys are the names of the cases and the values are the dicts inherited from `allenrich_plots_defaults`.<br />
    The cases under `envs.cases` can inherit this options.<br />
- `marker_plots_defaults` *(`ns`)*:
    Default options for the plots to generate for the markers.<br />
    - `plot_type`:
        The type of the plot.<br />
        See <https://pwwang.github.io/biopipen.utils.R/reference/VizDEGs.html>.<br />
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
        Other arguments passed to [`biopipen.utils::VizDEGs()`](https://pwwang.github.io/biopipen.utils.R/reference/VizDEGs.html).<br />
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
        See <https://pwwang.github.io/scplotter/reference/EnrichmentPlot.html>.<br />
- `enrich_plots` *(`type=json`)*: *Default: `{'Bar Plot': Diot({'plot_type': 'bar', 'ncol': 1, 'top_term': 10})}`*. <br />
    Cases of the plots to generate for the enrichment analysis.<br />
    The keys are the names of the cases and the values are the dicts inherited from `enrich_plots_defaults`.<br />
    The cases under `envs.cases` can inherit this options.<br />
- `overlaps_defaults` *(`ns`)*:
    Default options for investigating the overlapping of significant markers between different cases or comparisons.<br />
    This means either `ident_1` should be empty, so that they can be expanded to multiple comparisons.<br />
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
    1. If `ident_1` is not specified, the overlaps can be performed between different comparisons.<br />
    2. If `each` is specified, the overlaps can be performed between different cases, where in each case, `ident_1` must be specified.<br />

## SeeAlso

- [MarkersFinder](./MarkersFinder.md)
- [ClusterMarkersOfAllCells](./ClusterMarkersOfAllCells.md)
- [biopipen.ns.scrna.MarkersFinder](https://pwwang.github.io/biopipen/api/biopipen.ns.scrna/#biopipen.ns.scrna.MarkersFinder)

## Examples

### Visualize Log2 Fold Change of Markers

```toml
[ClusterMarkers.envs.marker_plots."Volcano Plot (log2FC)"]
plot_type = "volcano_log2fc"
```

![Volcano Plot (log2FC)](https://raw.githubusercontent.com/pwwang/immunopipe/tests-output/clustermarkers/ClusterMarkers/sampleinfo.markers/Cluster/seurat_clusters-c1/markers.Volcano-Plot-log2FC-.png)

### Visualize differential percentage of expression of Markers

```toml
[ClusterMarkers.envs.marker_plots."Volcano Plot (pct_diff)"]
plot_type = "volcano_pct"
```

![Volcano Plot (pct_diff)](https://raw.githubusercontent.com/pwwang/immunopipe/tests-output/clustermarkers/ClusterMarkers/sampleinfo.markers/Cluster/seurat_clusters-c1/markers.Volcano-Plot-diff_pct-.png)

### Visualize Average Expression of Markers with Dot Plot

```toml
[ClusterMarkers.envs.marker_plots."Dot Plot (AvgExp)"]
plot_type = "dotplot"
order_by = "desc(avg_log2FC)"
```

![Dot Plot (AvgExp)](https://raw.githubusercontent.com/pwwang/immunopipe/tests-output/clustermarkers/ClusterMarkers/sampleinfo.markers/Cluster/seurat_clusters-c1/markers.Dot-Plot.png)

### Visualize Average Expression of Markers with Heatmap

```toml
[ClusterMarkers.envs.marker_plots."Heatmap (AvgExp)"]
plot_type = "heatmap"
order_by = "desc(avg_log2FC)"
```

![Heatmap (AvgExp)](https://raw.githubusercontent.com/pwwang/immunopipe/tests-output/clustermarkers/ClusterMarkers/sampleinfo.markers/Cluster/seurat_clusters-c1/markers.Heatmap-of-Expressions-of-Top-Markers.png)

### Visualize Expression of Markers with Violin Plots

```toml
[ClusterMarkers.envs.marker_plots."Violin Plots"]
plot_type = "violin"
```

![Violin Plots](https://raw.githubusercontent.com/pwwang/immunopipe/tests-output/clustermarkers/ClusterMarkers/sampleinfo.markers/Cluster/seurat_clusters-c1/markers.Violin-Plots-for-Top-Markers.png)

### Visualize enrichment analysis results with Bar/EnrichMap/Network/WordCloud Plots

```toml
# Visualize enrichment of markers
[ClusterMarkers.envs.enrich_plots."Bar Plot"]  # Default
plot_type = "bar"

[ClusterMarkers.envs.enrich_plots."Network"]
plot_type = "network"

[ClusterMarkers.envs.enrich_plots."Enrichmap"]
plot_type = "enrichmap"

[ClusterMarkers.envs.enrich_plots."Word Cloud"]
plot_type = "wordcloud"
```

![Bar Plot](https://raw.githubusercontent.com/pwwang/immunopipe/tests-output/clustermarkers/ClusterMarkers/sampleinfo.markers/Cluster/seurat_clusters-c1/enrich.MSigDB_Hallmark_2020.Bar-Plot.png)
![Network](https://raw.githubusercontent.com/pwwang/immunopipe/tests-output/clustermarkers/ClusterMarkers/sampleinfo.markers/Cluster/seurat_clusters-c1/enrich.MSigDB_Hallmark_2020.Network.png)
![Enrichmap](https://raw.githubusercontent.com/pwwang/immunopipe/tests-output/clustermarkers/ClusterMarkers/sampleinfo.markers/Cluster/seurat_clusters-c1/enrich.MSigDB_Hallmark_2020.Enrichmap.png)
![Word Cloud](https://raw.githubusercontent.com/pwwang/immunopipe/tests-output/clustermarkers/ClusterMarkers/sampleinfo.markers/Cluster/seurat_clusters-c1/enrich.MSigDB_Hallmark_2020.Word-Cloud.png)

### Visualize top markers of all clusters with Heatmap

```toml
[ClusterMarkers.envs.allmarker_plots."Top 10 markers of all clusters"]
plot_type = "heatmap"
```

![Top 10 markers of all clusters](https://raw.githubusercontent.com/pwwang/immunopipe/tests-output/clustermarkers/ClusterMarkers/sampleinfo.markers/Cluster/seurat_clusters-All-Markers-/Top-10-markers-of-all-clusters.png)

### Visualize Log2 Fold Change of all markers

```toml
[ClusterMarkers.envs.allmarker_plots."Log2 Fold Change of all markers"]
plot_type = "heatmap_log2fc"
subset_by = "seurat_clusters"
```

![Log2 Fold Change of all markers](https://raw.githubusercontent.com/pwwang/immunopipe/tests-output/clustermarkers/ClusterMarkers/sampleinfo.markers/Cluster/seurat_clusters-All-Markers-/Log2FC-of-all-clusters.png)

### Visualize all markers in all clusters with Jitter Plots

```toml
[ClusterMarkers.envs.allmarker_plots."Jitter Plots of all markers"]
plot_type = "jitter"
subset_by = "seurat_clusters"
```

![Jitter Plots of all markers](https://raw.githubusercontent.com/pwwang/immunopipe/tests-output/clustermarkers/ClusterMarkers/sampleinfo.markers/Cluster/seurat_clusters-All-Markers-/Jitter-Plots-for-all-clusters.png)

### Visualize all enrichment analysis results of all clusters

```toml
[ClusterMarkers.envs.allenrich_plots."Heatmap of enriched terms of all clusters"]
plot_type = "heatmap"
```

![Heatmap of enriched terms of all clusters](https://raw.githubusercontent.com/pwwang/immunopipe/tests-output/clustermarkers/ClusterMarkers/sampleinfo.markers/Cluster/seurat_clusters-All-Enrichments-/allenrich.MSigDB_Hallmark_2020.Heatmap-of-enriched-terms-of-all-clusters.png)

### Overlapping markers

```toml
[ClusterMarkers.envs.overlaps."Overlapping Markers"]
plot_type = "venn"
```

![Overlapping Markers](https://raw.githubusercontent.com/pwwang/immunopipe/tests-output/clustermarkers/ClusterMarkers/sampleinfo.markers/Cluster/seurat_clusters-Overlaps-/Overlapping-Markers.png)

