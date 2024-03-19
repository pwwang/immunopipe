# SeuratClusterStats

Statistics of the clustering.

Including the number/fraction of cells in each cluster, the gene expression values
and dimension reduction plots. It's also possible to perform stats on
TCR clones/clusters or other metadata for each T-cell cluster.<br />

## Environment Variables

- `mutaters` *(`type=json`)*: *Default: `{}`*. <br />
    The mutaters to mutate the metadata to subset the cells.<br />
    The mutaters will be applied in the order specified.<br />
- `hists_defaults` *(`ns`)*:
    The default parameters for histograms.<br />
    This will plot histograms for the number of cells along `x`.<br />
    For example, you can plot the number of cells along cell activity score.<br />
    - `x`:
        The column name in metadata to plot as the x-axis.<br />
        The NA values will be removed.<br />
        It could be either numeric or factor/character.<br />
    - `x_order` *(`list`)*: *Default: `[]`*. <br />
        The order of the x-axis, only works for factor/character `x`.<br />
        You can also use it to subset `x` (showing only a subset values of `x`).<br />
    - `cells_by`:
        A column name in metadata to group the cells.<br />
        The NA values will be removed. It should be a factor/character.<br />
        if not specified, all cells will be used.<br />
    - `cells_order` *(`list`)*: *Default: `[]`*. <br />
        The order of the cell groups for the plots.<br />
        It should be a list of strings. You can also use `cells_orderby` and `cells_n`
        to determine the order.<br />
    - `cells_orderby`:
        An expression passed to `dplyr::arrange()` to order the cell groups.<br />
    - `cells_n`: *Default: `10`*. <br />
        The number of cell groups to show.<br />
        Ignored if `cells_order` is specified.<br />
    - `ncol` *(`type=int`)*: *Default: `2`*. <br />
        The number of columns for the plots, split by `cells_by`.<br />
    - `subset`:
        An expression to subset the cells, will be passed to `dplyr::filter()`.<br />
    - `each`:
        Whether to plot each group separately.<br />
    - `bins`: *Default: `30`*. <br />
        The number of bins to use, only works for numeric `x`.<br />
    - `plus` *(`list`)*: *Default: `[]`*. <br />
        The extra elements to add to the `ggplot` object.<br />
    - `devpars` *(`ns`)*:
        The device parameters for the plots.<br />
        - `res` *(`type=int`)*: *Default: `100`*. <br />
            The resolution of the plots.<br />
        - `height` *(`type=int`)*:
            The height of the plots.<br />
        - `width` *(`type=int`)*:
            The width of the plots.<br />
- `hists` *(`type=json`)*: *Default: `{}`*. <br />
    The cases for histograms.<br />
    Keys are the names of the plots and values are the dicts inherited from `env.hists_defaults`.<br />
    There is no default case.<br />
- `stats_defaults` *(`ns`)*:
    The default parameters for `stats`.<br />
    This is to do some basic statistics on the clusters. For more comprehensive analysis,
    see `RadarPlots` and `CellsDistribution`.<br />
    The parameters from the cases can overwrite the default parameters.<br />
    - `frac` *(`flag`)*: *Default: `False`*. <br />
        Whether to output the fraction of cells instead of number.<br />
    - `pie` *(`flag`)*: *Default: `False`*. <br />
        Also output a pie chart?<br />
    - `circos` *(`flag`)*: *Default: `False`*. <br />
        Also output a circos plot?<br />
    - `table` *(`flag`)*: *Default: `False`*. <br />
        Whether to output a table (in tab-delimited format) and in the report.<br />
    - `frac_ofall(flag)`:
        Whether to output the fraction against all cells,
        instead of the fraction in each group.<br />
        Does not work for circos plot.<br />
        Only works when `frac` is `True` and `group-by` is specified.<br />
    - `transpose` *(`flag`)*: *Default: `False`*. <br />
        Whether to transpose the cluster and group, that is,
        using group as the x-axis and cluster to fill the plot.<br />
        For circos plot, when transposed, the arrows will be drawn from the idents (by `ident`) to the
        the groups (by `group-by`).<br />
        Only works when `group-by` is specified.<br />
    - `position` *(`choice`)*: *Default: `auto`*. <br />
        The position of the bars. Does not work for pie and circos plots.<br />
        - `stack`:
            Use `position_stack()`.<br />
        - `fill`:
            Use `position_fill()`.<br />
        - `dodge`:
            Use `position_dodge()`.<br />
        - `auto`:
            Use `stack` when there are more than 5 groups, otherwise use `dodge`.<br />
    - `ident`: *Default: `seurat_clusters`*. <br />
        The column name in metadata to use as the identity.<br />
    - `group-by`:
        The column name in metadata to group the cells.<br />
        Does NOT support for pie charts.<br />
    - `split-by`:
        The column name in metadata to split the cells into different plots.<br />
        Does NOT support for circos plots.<br />
    - `subset`:
        An expression to subset the cells, will be passed to
        `dplyr::filter()` on metadata.<br />
    - `circos_labels_rot` *(`flag`)*: *Default: `False`*. <br />
        Whether to rotate the labels in the circos plot.<br />
        In case the labels are too long.<br />
    - `circos_devpars` *(`ns`)*:
        The device parameters for the circos plots.<br />
        - `res` *(`type=int`)*: *Default: `100`*. <br />
            The resolution of the plots.<br />
        - `height` *(`type=int`)*: *Default: `600`*. <br />
            The height of the plots.<br />
        - `width` *(`type=int`)*: *Default: `600`*. <br />
            The width of the plots.<br />
    - `pie_devpars` *(`ns`)*:
        The device parameters for the pie charts.<br />
        - `res` *(`type=int`)*: *Default: `100`*. <br />
            The resolution of the plots.<br />
        - `height` *(`type=int`)*: *Default: `600`*. <br />
            The height of the plots.<br />
        - `width` *(`type=int`)*: *Default: `800`*. <br />
            The width of the plots.<br />
    - `devpars` *(`ns`)*:
        The device parameters for the plots.<br />
        - `res` *(`type=int`)*: *Default: `100`*. <br />
            The resolution of the plots.<br />
        - `height` *(`type=int`)*: *Default: `600`*. <br />
            The height of the plots.<br />
        - `width` *(`type=int`)*: *Default: `800`*. <br />
            The width of the plots.<br />
    - `frac_ofall`: *Default: `False`*. <br />
- `stats` *(`type=json`)*: *Default: `{'Number of cells in each cluster': Diot({'pie': True}), 'Number of cells in each cluster by Sample': Diot({'group-by': 'Sample', 'table': True, 'frac': True})}`*. <br />
    The number/fraction of cells to plot.<br />
    Keys are the names of the plots and values are the dicts inherited from `env.stats_defaults`.<br />
    Here are some examples -

    ```python
    {
        "nCells_All": {},
        "nCells_Sample": {"group-by": "Sample"},
        "fracCells_Sample": {"frac": True, "group-by": "Sample"},
    }
    ```

- `ngenes_defaults` *(`ns`)*:
    The default parameters for `ngenes`.<br />
    The default parameters to plot the number of genes expressed in each cell.<br />
    - `ident`: *Default: `seurat_clusters`*. <br />
        The column name in metadata to use as the identity.<br />
    - `group-by`:
        The column name in metadata to group the cells.<br />
        Dodge position will be used to separate the groups.<br />
    - `split-by`:
        The column name in metadata to split the cells into different plots.<br />
    - `subset`:
        An expression to subset the cells, will be passed to `tidyrseurat::filter()`.<br />
    - `devpars` *(`ns`)*:
        The device parameters for the plots.<br />
        - `res` *(`type=int`)*: *Default: `100`*. <br />
            The resolution of the plots.<br />
        - `height` *(`type=int`)*: *Default: `800`*. <br />
            The height of the plots.<br />
        - `width` *(`type=int`)*: *Default: `1000`*. <br />
            The width of the plots.<br />
- `ngenes` *(`type=json`)*: *Default: `{'Number of genes expressed in each cluster': Diot({})}`*. <br />
    The number of genes expressed in each cell.<br />
    Keys are the names of the plots and values are the dicts inherited from `env.ngenes_defaults`.<br />
- `features_defaults` *(`ns`)*:
    The default parameters for `features`.<br />
    - `features`:
        The features to plot.<br />
        It can be either a string with comma separated features, a list of features, a file path with `file://` prefix with features
        (one per line), or an integer to use the top N features from `VariantFeatures(srtobj)`.<br />
    - `ident`: *Default: `seurat_clusters`*. <br />
        The column name in metadata to use as the identity.<br />
        If it is from subclustering (reduction `sub_umap_<ident>` exists), the reduction will be used.<br />
    - `cluster_orderby` *(`type=auto`)*:
        The order of the clusters to show on the plot.<br />
        An expression passed to `dplyr::summarise()` on the grouped data frame (by `seurat_clusters`).<br />
        The summary stat will be passed to `dplyr::arrange()` to order the clusters. It's applied on the whole meta.data before grouping and subsetting.<br />
        For example, you can order the clusters by the activation score of
        the cluster: `desc(mean(ActivationScore, na.rm = TRUE))`, suppose you have a column
        `ActivationScore` in the metadata.<br />
        You may also specify the literal order of the clusters by a list of strings.<br />
    - `subset`:
        An expression to subset the cells, will be passed to `tidyrseurat::filter()`.<br />
    - `devpars` *(`ns`)*:
        The device parameters for the plots. Does not work for `table`.<br />
        - `res` *(`type=int`)*: *Default: `100`*. <br />
            The resolution of the plots.<br />
        - `height` *(`type=int`)*:
            The height of the plots.<br />
        - `width` *(`type=int`)*:
            The width of the plots.<br />
    - `plus`:
        The extra elements to add to the `ggplot` object. Does not work for `table`.<br />
    - `group-by`:
        Group cells in different ways (for example, orig.ident). Works for `ridge`, `vln`, and `dot`.<br />
        It also works for `feature` as `shape.by` being passed to [`Seurat::FeaturePlot`](https://satijalab.org/seurat/reference/featureplot).<br />
    - `split-by`:
        The column name in metadata to split the cells into different plots.<br />
        It works for `vln`, `feature`, and `dot`.<br />
    - `assay`:
        The assay to use.<br />
    - `layer`:
        The layer to use.<br />
    - `reduction`:
        The reduction to use. Only works for `feature`.<br />
    - `section`:
        The section to put the plot in the report.<br />
        If not specified, the case title will be used.<br />
    - `ncol` *(`type=int`)*: *Default: `2`*. <br />
        The number of columns for the plots.<br />
    - `kind` *(`choice`)*:
        The kind of the plot or table.<br />
        - `ridge`:
            Use `Seurat::RidgePlot`.<br />
        - `ridgeplot`:
            Same as `ridge`.<br />
        - `vln`:
            Use `Seurat::VlnPlot`.<br />
        - `vlnplot`:
            Same as `vln`.<br />
        - `violin`:
            Same as `vln`.<br />
        - `violinplot`:
            Same as `vln`.<br />
        - `feature`:
            Use `Seurat::FeaturePlot`.<br />
        - `featureplot`:
            Same as `feature`.<br />
        - `dot`:
            Use `Seurat::DotPlot`.<br />
        - `dotplot`:
            Same as `dot`.<br />
        - `bar`:
            Bar plot on an aggregated feature.<br />
            The features must be a single feature, which will be either an  existing feature or an expression
            passed to `dplyr::summarise()` (grouped by `ident`) on the existing features to create a new feature.<br />
        - `barplot`:
            Same as `bar`.<br />
        - `heatmap`:
            Use `Seurat::DoHeatmap`.<br />
        - `avgheatmap`:
            Plot the average expression of the features in each cluster as a heatmap.<br />
        - `table`:
            The table for the features, only gene expressions are supported.<br />
            (supported keys: ident, subset, and features).<br />
- `features` *(`type=json`)*: *Default: `{}`*. <br />
    The plots for features, include gene expressions, and columns from metadata.<br />
    Keys are the titles of the cases and values are the dicts inherited from `env.features_defaults`. It can also have other parameters from
    each Seurat function used by `kind`. Note that for argument name with `.`, you should use `-` instead.<br />
- `dimplots_defaults` *(`ns`)*:
    The default parameters for `dimplots`.<br />
    - `ident`: *Default: `seurat_clusters`*. <br />
        The identity to use.<br />
        If it is from subclustering (reduction `sub_umap_<ident>` exists), this reduction will be used if `reduction`
        is set to `dim` or `auto`.<br />
    - `group-by`:
        Same as `ident` if not specified, to define how the points are colored.<br />
    - `na_group`:
        The group name for NA values, use `None` to ignore NA values.<br />
    - `split-by`:
        The column name in metadata to split the cells into different plots.<br />
    - `shape-by`:
        The column name in metadata to use as the shape.<br />
    - `subset`:
        An expression to subset the cells, will be passed to `tidyrseurat::filter()`.<br />
    - `devpars` *(`ns`)*:
        The device parameters for the plots.<br />
        - `res` *(`type=int`)*: *Default: `100`*. <br />
            The resolution of the plots.<br />
        - `height` *(`type=int`)*: *Default: `800`*. <br />
            The height of the plots.<br />
        - `width` *(`type=int`)*: *Default: `1000`*. <br />
            The width of the plots.<br />
    - `reduction` *(`choice`)*: *Default: `dim`*. <br />
        Which dimensionality reduction to use.<br />
        - `dim`:
            Use `Seurat::DimPlot`.<br />
            First searches for `umap`, then `tsne`, then `pca`.<br />
            If `ident` is from subclustering, `sub_umap_<ident>` will be used.<br />
        - `auto`:
            Same as `dim`
        - `umap`:
            Use `Seurat::UMAPPlot`.<br />
        - `tsne`:
            Use `Seurat::TSNEPlot`.<br />
        - `pca`:
            Use `Seurat::PCAPlot`.<br />
    - `<more>`:
        See <https://satijalab.org/seurat/reference/dimplot>
- `dimplots` *(`type=json`)*: *Default: `{'Dimensional reduction plot': Diot({'label': True, 'label-box': True, 'repel': True}), 'TCR presence': Diot({'ident': 'TCR_Presence', 'order': 'TCR_absent', 'cols': ['#FF000066', 'gray']})}`*. <br />
    The dimensional reduction plots.<br />
    Keys are the titles of the plots and values are the dicts inherited from `env.dimplots_defaults`. It can also have other parameters from
    [`Seurat::DimPlot`](https://satijalab.org/seurat/reference/dimplot).<br />

## Examples

### Number of cells in each cluster

```toml
[SeuratClusterStats.envs.stats]
# suppose you have nothing set in `envs.stats_defaults`
# otherwise, the settings will be inherited here
nCells_All = { }
```

![nCells_All](../processes/images/SeuratClusterStats_nCells_All.png){: width="80%" }

### Number of cells in each cluster by groups

```toml
[SeuratClusterStats.envs.stats]
nCells_Sample = { group-by = "Sample" }
```

![nCells_Sample](../processes/images/SeuratClusterStats_nCells_Sample.png){: width="80%" }

### Violin plots for the gene expressions

```toml
[SeuratClusterStats.envs.features]
features = "CD4,CD8A"
# Remove the dots in the violin plots
vlnplots = { pt-size = 0, kind = "vln" }
# Don't use the default genes
vlnplots_1 = { features = ["FOXP3", "IL2RA"], pt-size = 0, kind = "vln" }
```

![vlnplots](../processes/images/SeuratClusterStats_vlnplots.png){: width="80%" }
![vlnplots_1](../processes/images/SeuratClusterStats_vlnplots_1.png){: width="80%" }

### Dimension reduction plot with labels

```toml
[SeuratClusterStats.envs.dimplots.Idents]
label = true
label-box = true
repel = true
```

![dimplots](../processes/images/SeuratClusterStats_dimplots.png){: width="80%" }

