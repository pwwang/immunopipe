# CellTypeAnnotation

Annotate the T cell clusters.

Annotate the cell clusters. Currently, four ways are supported:<br />

1. Pass the cell type annotation directly
2. Use [`ScType`](https://github.com/IanevskiAleksandr/sc-type)
3. Use [`scCATCH`](https://github.com/ZJUFanLab/scCATCH)
4. Use [`hitype`](https://github.com/pwwang/hitype)

The annotated cell types will replace the original `seurat_clusters` column in the metadata,
so that the downstream processes will use the annotated cell types.<br />

The old `seurat_clusters` column will be renamed to `seurat_clusters_id`.<br />

If you are using `ScType`, `scCATCH`, or `hitype`, a text file containing the mapping from
the old `seurat_clusters` to the new cell types will be generated and saved to
`cluster2celltype.tsv` under `<workdir>/<pipline_name>/CellTypeAnnotation/0/output/`.<br />

The `<workdir>` is typically `./.pipen` and the `<pipline_name>` is `Immunopipe`
by default.<br />

/// Note
When supervised clustering [`SeuratMap2Ref`](./SeuratMap2Ref.md) is used, this
process will be ignored.<br />
///

/// Note
When cell types are annotated, the old `seurat_clusters` column will be renamed
to `seurat_clusters_id`, and the new `seurat_clusters` column will be added.<br />
///

## Environment Variables

- `tool` *(`choice`)*: *Default: `direct`*. <br />
    The tool to use for cell type annotation.<br />
    - `sctype`:
        Use `scType` to annotate cell types.<br />
        See <https://github.com/IanevskiAleksandr/sc-type>
    - `hitype`:
        Use `hitype` to annotate cell types.<br />
        See <https://github.com/pwwang/hitype>
    - `sccatch`:
        Use `scCATCH` to annotate cell types.<br />
        See <https://github.com/ZJUFanLab/scCATCH>
    - `direct`:
        Directly assign cell types
- `sctype_tissue`:
    The tissue to use for `sctype`.<br />
    Avaiable tissues should be the first column (`tissueType`) of `sctype_db`.<br />
    If not specified, all rows in `sctype_db` will be used.<br />
- `sctype_db`:
    The database to use for sctype.<br />
    Check examples at <https://github.com/IanevskiAleksandr/sc-type/blob/master/ScTypeDB_full.xlsx>
- `hitype_tissue`:
    The tissue to use for `hitype`.<br />
    Avaiable tissues should be the first column (`tissueType`) of `hitype_db`.<br />
    If not specified, all rows in `hitype_db` will be used.<br />
- `hitype_db`:
    The database to use for hitype.<br />
    Compatible with `sctype_db`.<br />
    See also <https://pwwang.github.io/hitype/articles/prepare-gene-sets.html>
    You can also use built-in databases, including `hitypedb_short`, `hitypedb_full`, and `hitypedb_pbmc3k`.<br />
- `cell_types` *(`list`)*: *Default: `[]`*. <br />
    The cell types to use for direct annotation.<br />
    You can use `"-"` or `""` as the placeholder for the clusters that
    you want to keep the original cell types (`seurat_clusters`).<br />
    If the length of `cell_types` is shorter than the number of
    clusters, the remaining clusters will be kept as the original cell
    types.<br />
    You can also use `NA` to remove the clusters from downstream analysis. This
    only works when `envs.newcol` is not specified.<br />

    /// Note
    If `tool` is `direct` and `cell_types` is not specified or an empty list,
    the original cell types will be kept and nothing will be changed.<br />
    ///

- `sccatch_args` *(`ns`)*:
    The arguments for `scCATCH::findmarkergene()` if `tool` is `sccatch`.<br />
    - `species` *(`choice`)*:
        The specie of cells.<br />
        - `Human`:
            Human cells.<br />
        - `Mouse`:
            Mouse cells.<br />
    - `cancer`:
        If the sample is from cancer tissue, then the cancer type may be defined.<br />
    - `tissue`:
        Tissue origin of cells must be defined.<br />
    - `<more>`:
        Other arguments for [`scCATCH::findmarkergene()`](https://rdrr.io/cran/scCATCH/man/findmarkergene.html).<br />
        You can pass an RDS file to `sccatch_args.marker` to work as custom marker. If so,
        `if_use_custom_marker` will be set to `TRUE` automatically.<br />
- `newcol`:
    The new column name to store the cell types.<br />
    If not specified, the `seurat_clusters` column will be overwritten.<br />
    If specified, the original `seurat_clusters` column will be kept and `Idents` will be kept as the original `seurat_clusters`.<br />

## Examples

```toml
[CellTypeAnnotation.envs]
tool = "direct"
cell_types = ["CellType1", "CellType2", "-", "CellType4"]
```

The cell types will be assigned as:<br />

```
0 -> CellType1
1 -> CellType2
2 -> 2
3 -> CellType4
```

## Metadata

When `envs.tool` is `direct` and `envs.cell_types` is empty, the metadata of
the `Seurat` object will be kept as is.<br />

When `envs.newcol` is specified, the original `seurat_clusters` column will
be kept is, and the annotated cell types will be saved in the new column.<br />
Otherwise, the original `seurat_clusters` column will be replaced by the
annotated cell types and the original `seurat_clusters` column will be
saved at `seurat_clusters_id`.<br />

![CellTypeAnnotation-metadata](../processes/images/CellTypeAnnotation-metadata.png)

