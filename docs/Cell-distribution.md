This module generates a set of pie charts showing how given groups of cells (e.g. a TCR clone/cluster) are distributed in the Seurat clusters (clustering by scRNA-seq data).

An example configuration:

```toml

[CellsDistribution.envs]
# # Show cells/clone/TCR cluster distribution in Seurat clusters
# # Cases (subsets/specific clones/clusters/etc)
[CellsDistribution.envs.cases.Tumor_vs_Normal]
# # Focus on a subset of the cells
# filter = "region != 'other'"
# # You can also create new columns in the data
[CellsDistribution.envs.cases.Tumor_vs_Normal.mutaters]
Region = "if_else (region == 'Tumor', 'Tumor', 'Normal')"
# # Columns of the pie charts
[CellsDistribution.envs.cases.Tumor_vs_Normal.group]
by = "Region"
order = ["Tumor", "Normal"]
# # Rows of the pie charts
[CellsDistribution.envs.cases.Tumor_vs_Normal.cells]
# # Select the cells/clone/clusters by?
by = "TCR_Cluster"
# # How many?
n = 10
# # order by? Will be passed to `dplyr::arrange()`
# # Available size columns to use: .CloneSize, .CloneGroupSize, .CloneGroupClusterSize
orderby = "desc(.CloneSize)"
```