Dimention reduction plots.

They are produced by process `DimPlots`. Configurations could specified to
`Dimplots.envs.cases`, where keys are names for the plot file, values are the
arguments that are passed to `Seurat::DimPlots()`

One of the important argument is `group.by`, which uses the meta data to group the cells.
You can also add new meta data columns based off of existing columns. See examples below
with argument `mutate`:

```toml
[DimPlots.envs.cases.Ident_UMAP]
"group.by" = "ident"
reduction = "umap"

[DimPlots.envs.cases.Source_TimePoint_TSNE]
"group.by" = "Source_TimePoint"
# Create a new column (Source_TimePoint) in metadata by concatenate those 2 original columns
mutate = "paste0(Source, '/', TimePoint)"
reduction = "tsne"
```
