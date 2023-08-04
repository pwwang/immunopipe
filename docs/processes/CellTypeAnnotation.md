# CellTypeAnnotation

Annotate the T cell clusters. Currently, three ways are supported:

1. Pass the cell type annotation directly
2. Use [`ScType`][1]
3. Use [`scCATCH`][2]

The annotated cell types will replace the original `seurat_clusters` column in the metadata, so that the downstream processes will use the annotated cell types.

## Environment variables

- `tool` (`choice`): The tool to use for cell type annotation.
    - `sctype`: Use `ScType` to annotate cell types.
       See <https://github.com/IanevskiAleksandr/sc-type>
    - `sccatch`: Use `scCATCH` to annotate cell types.
       See <https://github.com/ZJUFanLab/scCATCH>
    - `direct`: Directly assign cell types
- `sctype_tissue`: The tissue to use for `sctype`.
   Avaiable tissues should be the first column (`tissueType`) of `sctype_db`.
   Examples are `Immune system`, `Pancreas`, `Liver`, `Eye`, `Kidney`,
   `Brain`, `Lung`, `Adrenal`, `Heart`, `Intestine`, `Muscle`,
   `Placenta`, `Spleen`, `Stomach` and `Thymus` for the example `sctype_db`.
- `sctype_db`: The database to use for sctype.
   Check examples at <https://github.com/IanevskiAleksandr/sc-type/blob/master/ScTypeDB_full.xlsx>
- `cell_types` (`list`): The cell types to use for direct annotation. You can use `"-"` or `""` as the placeholder for the clusters that you want to keep the original cell types (`seurat_clusters`). If the length of `cell_types` is shorter than the number of clusters, the remaining clusters will be kept as the original cell types.

    !!! note
        If `tool` is `direct` and `cell_types` is not specified or
        an empty list, the original cell types will be kept and nothing
        will be changed.

- `sccatch_args` (`ns`): The arguments for [`scCATCH::findmarkergene()`][3] if `tool` is `sccatch`.
    - `species` (`choice`): The specie of cells.
        - `Human`: The sample is from human.
        - `Mouse`: The sample is from mouse.
    - `cancer`: If the sample is from cancer tissue, then the cancer type may be defined.
    - `tissue`: Tissue origin of cells must be defined.
    - <more>: Other arguments for `scCATCH::findmarkergene()`
       See <https://www.rdocumentation.org/packages/scCATCH/versions/3.2.2/topics/findmarkergene>.
       You can pass an RDS file to `sccatch_args.marker` to work as custom marker. If so,
       `if_use_custom_marker` will be set to `TRUE` automatically.

## Examples

### Directly assign cell types

```toml
[CellTypeAnnotation.envs]
tool = "direct"
cell_types = ["CellType1", "CellType2", "-", "CellType4"]
```

The cell types will be assigned as:

```
0 -> CellType1
1 -> CellType2
2 -> 2
3 -> CellType4
```

[1]: https://github.com/IanevskiAleksandr/sc-type
[2]: https://github.com/ZJUFanLab/scCATCH
[3]: https://rdrr.io/github/ZJUFanLab/scCATCH/man/findmarkergene.html