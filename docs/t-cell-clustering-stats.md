We have 3 types of statistics for T-cell clustering using Seurat.

## stats

Basic statistics of the clusters, including:

- `nCells` - Number of cells for each cluster
- `nCellsPerSample` - Number of cells per sample for each cluster
- `percCellsPerSample` - Percentage of cells per sample for each cluster

## exprs

The expression values to plot

- `genes` - The set of genes for the plots, unless features for those plots is specified. Could also specify a file with genes (one per line)
- `ridgeplots` - The ridge plots for the gene expressions. See ?Seurat::RidgePlot.
- `vlnplots` - Violin plots for the gene expressions. See ?Seurat::VlnPlot. You can have boxplot key to add geom_boxplot() to the violin plots
- `featureplots` - The feature plots for the gene expressions. See ?Seurat::FeaturePlot.
- `dotplot` - Dot plots for the gene expressions. See ?Seurat::DotPlot.
- `heatmap` - Heatmap for the gene expressions. See ?Seurat::DoHeatmap. You can specify average=True to plot on the average of the expressions.

All the above with devpars to define the output figures and plus to add elements to the ggplot object. You can have subset to subset the data. Multiple cases can be distinguished by `ridgeplots` and `ridgeplots.1`

## dimplots

The dimensional reduction plots

- `<case>` - The case to plot. Keys are the arguments for `Seurat::Dimplot()`, and `devpars` to control the plots.

## Example configurations

```toml

[SeuratClusterStats.envs]
# # Statistics about Seurat clusters
[SeuratClusterStats.envs.stats]
# # The statistics to plot
# # - nCells - Number of cells for each cluster
# # - nCellsPerSample - Number of cells per sample for each cluster
# # - percCellsPerSample - Percentage of cells per sample for each cluster
nCells = {res = 100, height = 1000, width = 1000}
nCellsPerSample = {res = 100, height = 1000, width = 1000}
percCellsPerSample = {res = 100, height = 1000, width = 1000}
# # The expression values to plot
[SeuratClusterStats.envs.exprs]
# # ridgeplots - The ridge plots for the gene expressions.
# # See `?Seurat::RidgePlot`.
# # vlnplots - Violin plots for the gene expressions.
# # See `?Seurat::VlnPlot`. You can have `boxplot` key to add
# # `geom_boxplot()` to the violin plots
# # featureplots - The feature plots for the gene expressions.
# # See `?Seurat::FeaturePlot`.
# # dotplot - Dot plots for the gene expressions.
# # See `?Seurat::DotPlot`.
# # heatmap - Heatmap for the gene expressions.
# # See `?Seurat::DoHeatmap`. You can specify `average=True` to plot on
# # the average of the expressions.
# # All the above with `devpars` to define the output figures
# # and `plus` to add elements to the `ggplot` object.
# # You can have `subset` to subset the data. Multiple cases can be
# # distinguished by `ridgeplots` and `ridgeplots.1`
"ridgeplots.1"  = {title = "Gene expressions in Tumor", subset = "region == 'Tumor'"}
"ridgeplots.2"  = {title = "Gene expressions in Normal", subset = "region == 'Normal adjacent tissue'"}
"vlnplots.1"  = {boxplot = {}, "pt.size" = 0, title = "Gene expressions in Tumor", subset = "region == 'Tumor'"}
"vlnplots.2"  = {boxplot = {}, "pt.size" = 0, title = "Gene expressions in Normal", subset = "region == 'Normal adjacent tissue'"}
"featureplots.1"  = {title = "Gene expressions in Tumor", subset = "region == 'Tumor'"}
"featureplots.2"  = {title = "Gene expressions in Normal", subset = "region == 'Normal adjacent tissue'"}
"dotplot.1"  = {plus = "RotatedAxis()", title = "Gene expressions in Tumor", subset = "region == 'Tumor'"}
"dotplot.2"  = {plus = "RotatedAxis()", title = "Gene expressions in Normal", subset = "region == 'Normal adjacent tissue'"}
"heatmap.1"  = {downsample = "average", title = "Gene expressions in Tumor", subset = "region == 'Tumor'"}
"heatmap.2"  = {downsample = "average", title = "Gene expressions in Normal", subset = "region == 'Normal adjacent tissue'"}
# # The dimensional reduction plots
[SeuratClusterStats.envs.dimplots]
# # `<case>` - The case to plot. Keys are the arguments for `Seurat::Dimplot()`, add `devpars`.
Ident = {"group.by" = "ident", devpars = {res = 100, width = 1000, height = 1000}}
```