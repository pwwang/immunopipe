# ClusterMarkersOfAllCells

Markers for clusters of all cells.

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
    This only works when `use_presto` is `False`.<br />
    - `<more>`:
        See <https://satijalab.org/seurat/reference/findmarkers>
- `dotplot` *(`ns`)*:
    Arguments for `Seurat::DotPlot()`.<br />
    Use `-` to replace `.` in the argument name. For example,
    use `group-bar` instead of `group.bar`.<br />
    Note that `object`, `features`, and `group-by` are already specified
    by this process. So you don't need to specify them here.<br />
    - `maxgenes` *(`type=int`)*: *Default: `20`*. <br />
        The maximum number of genes to plot.<br />
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
- `overlap_defaults` *(`ns`)*:
    The default options for overlapping analysis.<br />
    - `venn` *(`ns`)*:
        The options for the Venn diagram.<br />
        Venn diagram can only be plotted for sections with no more than 4 cases.<br />
        - `devpars` *(`ns`)*:
            The device parameters for the plots.<br />
            - `res` *(`type=int`)*: *Default: `100`*. <br />
                The resolution of the plots.<br />
            - `height` *(`type=int`)*: *Default: `600`*. <br />
                The height of the plots.<br />
            - `width` *(`type=int`)*: *Default: `1000`*. <br />
                The width of the plots.<br />
    - `upset` *(`ns`)*:
        The options for the UpSet plot.<br />
        - `devpars` *(`ns`)*:
            The device parameters for the plots.<br />
            - `res` *(`type=int`)*: *Default: `100`*. <br />
                The resolution of the plots.<br />
            - `height` *(`type=int`)*: *Default: `600`*. <br />
                The height of the plots.<br />
            - `width` *(`type=int`)*: *Default: `800`*. <br />
                The width of the plots.<br />
- `overlap` *(`json`)*: *Default: `{}`*. <br />
    The sections to do overlaping analysis, including
    Venn diagram and UpSet plot. The Venn diagram and UpSet plot
    will be plotted for the overlapping of significant markers between
    different cases.<br />
    The keys of this option are the names of the sections. The values are
    a dict of options with keys `venn` and `upset`, values will
    be inherited from `envs.overlap_defaults`, recursively.<br />
    You can set `envs.overlap.<section>.venn` to `False`/`None` to disable
    the Venn diagram for the section.<br />
    It works when `each` is specified. In such a case, the sections will be
    the case names.<br />
    This does not work for the cases where `ident-1` is not specified. In case
    you want to do such analysis for those cases, you should enumerate the
    idents in different cases and specify them here.<br />
- `cache` *(`type=auto`)*: *Default: `/tmp/user`*. <br />
    Where to cache to `FindAllMarkers` results.<br />
    If `True`, cache to `outdir` of the job. If `False`, don't cache.<br />
    Otherwise, specify the directory to cache to.<br />
    Only works when `use_presto` is `False` (presto works fast enough).<br />

