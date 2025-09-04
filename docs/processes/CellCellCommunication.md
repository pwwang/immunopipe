# CellCellCommunication

Cell-cell communication inference

This is implemented based on [LIANA](https://liana-py.readthedocs.io/en/latest/index.html),
which is a Python package for cell-cell communication inference and provides a list of existing
methods including [CellPhoneDB](https://github.com/ventolab/CellphoneDB),
[Connectome](https://github.com/msraredon/Connectome/), log2FC,
[NATMI](https://github.com/forrest-lab/NATMI),
[SingleCellSignalR](https://github.com/SCA-IRCM/SingleCellSignalR), Rank_Aggregate, Geometric Mean,
[scSeqComm](https://gitlab.com/sysbiobig/scseqcomm), and [CellChat](https://github.com/jinworks/CellChat).<br />

You can also try `python -c 'import liana; liana.mt.show_methods()'` to see the methods available.<br />

Note that this process does not do any visualization. You can use `CellCellCommunicationPlots`
to visualize the results.<br />

## Input

- `sobjfile`:
    The seurat object file in RDS or h5seurat format or AnnData file.<br />

## Output

- `outfile`: *Default: `{{in.sobjfile | stem}}-ccc.txt`*. <br />
    The output file with the 'liana_res' data frame.<br />
    Stats are provided for both ligand and receptor entities, more specifically: ligand and receptor are
    the two entities that potentially interact. As a reminder, CCC events are not limited to secreted signalling,
    but we refer to them as ligand and receptor for simplicity.<br />
    Also, in the case of heteromeric complexes, the ligand and receptor columns represent the subunit with minimum
    expression, while *_complex corresponds to the actual complex, with subunits being separated by _.<br />
    source and target columns represent the source/sender and target/receiver cell identity for each interaction, respectively
    * `*_props`: represents the proportion of cells that express the entity.<br />
    By default, any interactions in which either entity is not expressed in above 10%% of cells per cell type
    is considered as a false positive, under the assumption that since CCC occurs between cell types, a sufficient
    proportion of cells within should express the genes.<br />
    * `*_means`: entity expression mean per cell type.<br />
    * `lr_means`: mean ligand-receptor expression, as a measure of ligand-receptor interaction magnitude.<br />
    * `cellphone_pvals`: permutation-based p-values, as a measure of interaction specificity.<br />

## Environment Variables

- `method` *(`choice`)*: *Default: `cellchat`*. <br />
    The method to use for cell-cell communication inference.<br />
    - `CellPhoneDB`:
        Use CellPhoneDB method.<br />
        Magnitude Score: lr_means; Specificity Score: cellphone_pvals.<br />
    - `Connectome`:
        Use Connectome method.<br />
    - `log2FC`:
        Use log2FC method.<br />
    - `NATMI`:
        Use NATMI method.<br />
    - `SingleCellSignalR`:
        Use SingleCellSignalR method.<br />
    - `Rank_Aggregate`:
        Use Rank_Aggregate method.<br />
    - `Geometric_Mean`:
        Use Geometric Mean method.<br />
    - `scSeqComm`:
        Use scSeqComm method.<br />
    - `CellChat`:
        Use CellChat method.<br />
    - `cellphonedb`:
        alias for `CellPhoneDB`
    - `connectome`:
        alias for `Connectome`
    - `log2fc`:
        alias for `log2FC`
    - `natmi`:
        alias for `NATMI`
    - `singlesignaler`:
        alias for `SingleCellSignalR`
    - `rank_aggregate`:
        alias for `Rank_Aggregate`
    - `geometric_mean`:
        alias for `Geometric_Mean`
    - `scseqcomm`:
        alias for `scSeqComm`
    - `cellchat`:
        alias for `CellChat`
- `subset`:
    An expression in string to subset the cells.<br />
    When a `.rds` or `.h5seurat` file is provided for `in.sobjfile`, you can provide an expression in `R`,
    which will be passed to `base::subset()` in `R` to subset the cells.<br />
    But you can always pass an expression in `python` to subset the cells.<br />
    See <https://anndata.readthedocs.io/en/latest/tutorials/notebooks/getting-started.html#subsetting-using-metadata>.<br />
    You should use `adata` to refer to the AnnData object. For example, `adata.obs.groups == "g1"` will subset the cells
    with `groups` equal to `g1`.<br />
- `subset_using`: *Default: `auto`*. <br />
    The method to subset the cells.<br />
    - `auto`:
        Automatically detect the method to use.<br />
        Note that this is not always accurate. We simply check if `[` is in the expression.<br />
        If so, we use `python` to subset the cells; otherwise, we use `R`.<br />
    - `python`:
        Use python to subset the cells.<br />
    - `r`:
        Use R to subset the cells.<br />
- `split_by`:
    The column name in metadata to split the cells to run the method separately.<br />
    The results will be combined together with this column in the final output.<br />
- `assay`:
    The assay to use for the analysis.<br />
    Only works for Seurat object.<br />
- `seed` *(`type=int`)*: *Default: `1337`*. <br />
    The seed for the random number generator.<br />
- `ncores` *(`type=int`)*: *Default: `1`*. <br />
    The number of cores to use.<br />
- `groupby`: *Default: `seurat_clusters`*. <br />
    The column name in metadata to group the cells.<br />
    Typically, this column should be the cluster id.<br />
- `species` *(`choice`)*: *Default: `human`*. <br />
    The species of the cells.<br />
    - `human`:
        Human cells, the 'consensus' resource will be used.<br />
    - `mouse`:
        Mouse cells, the 'mouseconsensus' resource will be used.<br />
- `expr_prop` *(`type=float`)*: *Default: `0.1`*. <br />
    Minimum expression proportion for the ligands and
    receptors (+ their subunits) in the corresponding cell identities. Set to 0
    to return unfiltered results.<br />
- `min_cells` *(`type=int`)*: *Default: `5`*. <br />
    Minimum cells (per cell identity if grouped by `groupby`)
    to be considered for downstream analysis.<br />
- `n_perms` *(`type=int`)*: *Default: `1000`*. <br />
    Number of permutations for the permutation test.<br />
    Relevant only for permutation-based methods (e.g., `CellPhoneDB`).<br />
    If `0` is passed, no permutation testing is performed.<br />
- `rscript`: *Default: `Rscript`*. <br />
    The path to the Rscript executable used to convert RDS file to AnnData.<br />
    if `in.sobjfile` is an RDS file, it will be converted to AnnData file (h5ad).<br />
    You need `Seurat`, `SeuratDisk` and `digest` installed.<br />
- `<more>`:
    Other arguments for the method.<br />
    The arguments are passed to the method directly.<br />
    See the method documentation for more details and also
    `help(liana.mt.<method>.__call__)` in Python.<br />

## Reference

- [Review](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9184522/).<br />
- [LIANA](https://www.biorxiv.org/content/10.1101/2023.08.19.553863v1).<br />

