# ClusterMarkers

This process is extended from [`MarkersFinder`][1] from the [`biopipen`][2] package. `MarkersFinder` is a `pipen` process that wraps the [`Seurat::FindMarkers()`][3] function, and performs enrichment analysis for the markers found.

The enrichment analysis is done by [`enrichr`][5].

## Environment variables

- `ncores`: Number of cores to use for parallel computing for some `Seurat` procedures.
    - Used in `future::plan(strategy = "multicore", workers = <ncores>)` to parallelize some Seurat procedures.
    - See also: <https://satijalab.org/seurat/articles/future_vignette.html>
- `group-by`: The column name in metadata with the cluster information. Default: `seurat_clusters`.
- `dbs` (`list`): The dbs to do enrichment analysis for significant
    markers See below for all libraries.
    <https://maayanlab.cloud/Enrichr/#libraries>
    Default:
    - `GO_Biological_Process_2021`
    - `GO_Cellular_Component_2021`
    - `GO_Molecular_Function_2021`
    - `KEGG_2021_Human`
- `sigmarkers`: An expression passed to [`dplyr::filter()`][4] to filter the
    significant markers for enrichment analysis.
    Available variables are `p_val`, `avg_log2FC`, `pct.1`, `pct.2` and
    `p_val_adj`.
    For example, `"p_val_adj < 0.05 & abs(avg_log2FC) > 1"` to select markers with adjusted p-value < 0.05 and absolute log2 fold change > 1.
- `rest` (`ns`): Rest arguments for [`Seurat::FindMarkers()`][3]. Use `-` to replace `.` in the argument name. For example, use `min-pct` instead of `min.pct`.
    - `<more>`: See <https://satijalab.org/seurat/reference/findmarkers>

!!! note
    Since this process is extended from `MarkersFinder`, other environment variables from `MarkersFinder` are also available. However, they should not be used in this process. Other environment variables are used for more complicated cases for marker finding (See [`MarkersFinder`][1] for more details).

    If you are using `pipen-board` to run the pipeline (see [here](../running.md#run-the-pipeline-via-pipen-board) and [here](../running.md#run-the-pipeline-via-pipen-board-using-docker-image)), you may see the other environment variables of this process are hidden and readonly.


[1]: https://pwwang.github.io/biopipen/api/biopipen.ns.scrna/#biopipen.ns.scrna.MarkersFinder
[2]: https://pwwang.github.io/biopipen
[3]: https://satijalab.org/seurat/reference/findmarkers
[4]: https://dplyr.tidyverse.org/reference/filter.html
[5]: https://maayanlab.cloud/Enrichr/
