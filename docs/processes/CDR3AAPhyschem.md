# CDR3AAPhyschem

CDR3 AA physicochemical feature analysis

The idea is to perform a regression between two groups of cells
(e.g. Treg vs Tconv) at different length of CDR3 AA sequences.
The regression will be performed for each physicochemical feature of the
AA (hydrophobicity, volume and isolectric point).

## Reference:

- [Stadinski, Brian D., et al. "Hydrophobic CDR3 residues promote the development of self-reactive T cells." Nature immunology 17.8 (2016): 946-955.](https://www.nature.com/articles/ni.3491)
- [Lagattuta, Kaitlyn A., et al. "Repertoire analyses reveal T cell antigen receptor sequence features that influence T cell fate." Nature immunology 23.3 (2022): 446-457.](https://www.nature.com/articles/s41590-022-01129-x)
- [Wimley, W. C. & White, S. H. Experimentally determined hydrophobicity scale for proteins at membrane - interfaces. Nat. Struct. Biol. 3, 842-848 (1996).](https://www.nature.com/articles/nsb1096-842)
- [Handbook of chemistry & physics 72nd edition. (CRC Press, 1991).](https://books.google.com/books?hl=en&lr=&id=bNDMBQAAQBAJ&oi=fnd&pg=PP1&dq=Hdbk+of+chemistry+%26+physics&ots=H9fzwhwz-C&sig=EXHI9N3q4OW9TYEBWlldqkvADfM#v=onepage&q=Hdbk%20of%20chemistry%20%26%20physics&f=false)
- [Zamyatnin, A. A. Protein volume in solution. Prog. Biophys. Mol. Biol. 24, 107-123 (1972).](https://www.sciencedirect.com/science/article/pii/0079610772900053)

## Environment variables

- `group`: The key of group in metadata to define the groups to
    compare. For example, `CellType`, which has cell types annotated
    for each cell in the combined object (immdata + Seurat metadata)
- `comparison` (`type=json`): A dict of two groups, with keys as the
    group names and values as the group labels. For example,
    ```toml
    Treg = ["CD4 CTL", "CD4 Naive", "CD4 TCM", "CD4 TEM"]
    Tconv = "Tconv"
    ```

- `prefix`: The prefix of the cell names (rownames) in the metadata.
    The prefix is usually not needed in immdata, as the data is stored
    in the `immdata` object separately for each sample. However, the
    `Seurat` object has a combined `meta.data` for all the samples,
    so the prefix is needed. Usually, the prefix is the sample name.
    For example, `Sample1-AACGTTGAGGCTACGT-1`.
    We need this prefix to add the sample name to the cell names in
    immdata, so that we can match the cells in `immdata` and
    `Seurat` object. Set it to `None` or an empty string if the
    `Seurat` object has the same cell names as `immdata`. You can use
    placeholders to specify the prefix, e.g., `{Sample}_`. In such a
    case, the `Sample` column must exist in the `Seurat` object.
- `target`: Which group to use as the target group. The target
    group will be labeled as 1, and the other group will be labeled as
    0 in the regression.
- `subset`: A column, or a list of columns separated by comma,
    in the merged object to subset the cells to perform the regression,
    for each group in the columns.
    If not provided, all the cells will be used.
