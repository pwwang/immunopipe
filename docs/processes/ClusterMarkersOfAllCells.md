# ClusterMarkersOfAllCells

Markers for clusters of all cells.

/// Tip | Added in 0.9.0
`ClusterMarkersOfAllCells` is added in `0.9.0` and is optional by default.<br />
///

See also [ClusterMarkers](./ClusterMarkers.md).<br />

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
- `dbs` *(`list`)*: *Default: `['GO_Biological_Process_2021', 'GO_Cellular_Component_2021', 'GO_Molecular_Function_2021', 'KEGG_2021_Human']`*. <br />
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
- `volcano_genes` *(`type=auto`)*: *Default: `True`*. <br />
    The genes to label in the volcano plot if they are
    significant markers.<br />
    If `True`, all significant markers will be labeled. If `False`, no
    genes will be labeled. Otherwise, specify the genes to label.<br />
    It could be either a string with comma separated genes, or a list
    of genes.<br />
- `subset`:
    An expression to subset the cells for each case.<br />
- `rest` *(`ns`)*:
    Rest arguments for `Seurat::FindMarkers()`.<br />
    Use `-` to replace `.` in the argument name. For example,
    use `min-pct` instead of `min.pct`.<br />
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
- `overlap` *(`list`)*: *Default: `[]`*. <br />
    The sections to do overlap analysis.<br />

