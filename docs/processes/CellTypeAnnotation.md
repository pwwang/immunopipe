# CellTypeAnnotation

Annotate all or selected T/B cell clusters.

Annotate the cell clusters. Currently, four ways are supported:<br />

1. Pass the cell type annotation directly
2. Use [`ScType`](https://github.com/IanevskiAleksandr/sc-type)
3. Use [`scCATCH`](https://github.com/ZJUFanLab/scCATCH)
4. Use [`hitype`](https://github.com/pwwang/hitype)
5. Use [`celltypist`](https://github.com/Teichlab/celltypist)

The annotated cell types will replace the original identity column in the metadata,
so that the downstream processes will use the annotated cell types.<br />

/// Note

When cell types are annotated, the original identity column (e.g. `seurat_clusters`) will be renamed
to `envs.backup_col` (e.g. `seurat_clusters_id`), and the new identity column will be added.<br />

///

If you are using `ScType`, `scCATCH`, or `hitype`, a text file containing the mapping from
the original identity to the new cell types will be generated and saved to
`cluster2celltype.tsv` under `<workdir>/<pipline_name>/CellTypeAnnotation/0/output/`.<br />

The `<workdir>` is typically `./.pipen` and the `<pipline_name>` is `Immunopipe`
by default.<br />

/// Note

If you have other annotation processes, including [`SeuratClustering`](./SeuratClustering.md)
process or [`SeuratMap2Ref`](./SeuratMap2Ref.md) process enabled in the same run,
you may want to specify a different name for the column to store the annotated cell types
using `envs.newcol`, so that the results from different annotation processes won't overwrite each other.<br />

///

## Input

- `sobjfile`:
    The single-cell object in RDS/qs/qs2/h5ad format.<br />

## Output

- `outfile`: *Default: `{{in.sobjfile | stem}}.annotated.{{- ext0(in.sobjfile) if envs.outtype == 'input' else envs.outtype -}}`*. <br />
    The rds/qs/qs2/h5ad file of seurat object with cell type annotated.<br />
    A text file containing the mapping from the old identity to the new cell types
    will be generated and saved to `cluster2celltype.tsv` under the job output directory.<br />
    Note that if `envs.ident` is specified, the output Seurat object will have
    the identity set to the specified column in metadata.<br />

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
    - `celltypist`:
        Use `celltypist` to annotate cell types.<br />
        See <https://github.com/Teichlab/celltypist>
    - `direct`:
        Directly assign cell types
- `sctype_tissue`:
    The tissue to use for `sctype`.<br />
    Avaiable tissues should be the first column (`tissueType`) of `sctype_db`.<br />
    If not specified, all rows in `sctype_db` will be used.<br />
- `sctype_db`:
    The database to use for sctype.<br />
    Check examples at <https://github.com/IanevskiAleksandr/sc-type/blob/master/ScTypeDB_full.xlsx>
- `ident`:
    The column name in metadata to use as the clusters.<br />
    If not specified, the identity column will be used when input is rds/qs/qs2 (supposing we have a Seurat object).<br />
    If input data is h5ad, this is required to run cluster-based annotation tools.<br />
    For `celltypist`, this is a shortcut to set `over_clustering` in `celltypist_args`.<br />
- `backup_col`: *Default: `seurat_clusters_id`*. <br />
    The backup column name to store the original identities.<br />
    If not specified, the original identity column will not be stored.<br />
    If `envs.newcol` is specified, this will be ignored.<br />
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
    you want to keep the original cell types.<br />
    If the length of `cell_types` is shorter than the number of
    clusters, the remaining clusters will be kept as the original cell
    types.<br />
    You can also use `NA` to remove the clusters from downstream analysis. This
    only works when `envs.newcol` is not specified.<br />

    /// Note
    If `tool` is `direct` and `cell_types` is not specified or an empty list,
    the original cell types will be kept and nothing will be changed.<br />
    ///

- `more_cell_types` *(`type=json`)*:
    The additional cell type annotations to add to the metadata.<br />
    The keys are the new column names and the values are the cell types lists.<br />
    The cell type lists work the same as `cell_types` above.<br />
    This is useful when you want to keep multiple annotations of cell types.<br />

- `sccatch_args` *(`ns`)*:
    The arguments for `scCATCH::findmarkergene()` if `tool` is `sccatch`.<br />
    - `species`:
        The specie of cells.<br />
    - `cancer`: *Default: `Normal`*. <br />
        If the sample is from cancer tissue, then the cancer type may be defined.<br />
    - `tissue`:
        Tissue origin of cells must be defined.<br />
    - `marker`:
        The marker genes for cell type identification.<br />
    - `if_use_custom_marker` *(`flag`)*: *Default: `False`*. <br />
        Whether to use custom marker genes. If `True`, no `species`, `cancer`, and `tissue` are needed.<br />
    - `<more>`:
        Other arguments for [`scCATCH::findmarkergene()`](https://rdrr.io/cran/scCATCH/man/findmarkergene.html).<br />
        You can pass an RDS file to `sccatch_args.marker` to work as custom marker. If so,
        `if_use_custom_marker` will be set to `TRUE` automatically.<br />
- `celltypist_args` *(`ns`)*:
    The arguments for `celltypist::celltypist()` if `tool` is `celltypist`.<br />
    - `model`:
        The path to model file.<br />
    - `python`: *Default: `python`*. <br />
        The python path where celltypist is installed.<br />
    - `majority_voting`: *Default: `True`*. <br />
        When true, it refines cell identities within local subclusters after an over-clustering approach
        at the cost of increased runtime.<br />
    - `over_clustering` *(`type=auto`)*:
        The column name in metadata to use as clusters for majority voting.<br />
        Set to `False` to disable over-clustering.<br />
        When `in.sobjfile` is rds/qs/qs2 (supposing we have a Seurat object), the default ident is used by default.<br />
        Otherwise, it is False by default.<br />
    - `assay`:
        When converting a Seurat object to AnnData, the assay to use.<br />
        If input is h5seurat, this defaults to RNA.<br />
        If input is Seurat object in RDS, this defaults to the default assay.<br />
- `merge` *(`flag`)*: *Default: `False`*. <br />
    Whether to merge the clusters with the same cell types.<br />
    Otherwise, a suffix will be added to the cell types (ie. `.1`, `.2`, etc).<br />
- `newcol`:
    The new column name to store the cell types.<br />
    If not specified, the identity column will be overwritten.<br />
    If specified, the original identity column will be kept and `Idents` will be kept as the original identity.<br />
- `outtype` *(`choice`)*: *Default: `input`*. <br />
    The output file type. Currently only works for `celltypist`.<br />
    An RDS file will be generated for other tools.<br />
    - `input`:
        Use the same file type as the input.<br />
    - `rds`:
        Use RDS file.<br />
    - `qs`:
        Use qs2 file.<br />
    - `qs2`:
        Use qs2 file.<br />
    - `h5ad`:
        Use AnnData file.<br />

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

When `envs.newcol` is specified, the original identity column (e.g. `seurat_clusters`) will
be kept is, and the annotated cell types will be saved in the new column.<br />
Otherwise, the original identity column will be replaced by the
annotated cell types and the original identity column will be
saved at `envs.backup_col` (e.g. `seurat_clusters_id`).<br />

![CellTypeAnnotation-metadata](images/CellTypeAnnotation-metadata.png)

