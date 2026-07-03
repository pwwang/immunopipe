# ScVelo

Velocity analysis for single-cell RNA-seq data

This process is implemented based on the Python package `scvelo` (v0.3.3).<br />
Note that it doesn't work with `numpy>=2`.<br />

## Input

- `sobjfile`:
    The seurat object file in RDS or h5seurat format or AnnData file.<br />

## Output

- `outfile`: *Default: `{{in.sobjfile | stem}}-scvelo.{{ext0(in.sobjfile) if envs.outtype == '<input>' else envs.outtype}}`*. <br />
    The output object with the velocity embeddings and information.<br />
    In either RDS, h5seurat or h5ad format, depending on the `envs.outtype`.<br />
    There will be also plots generated in the output directory
    (parent directory of `outfile`).<br />
    Note that these plots will not be used in the report, but can be used as
    supplementary information for the velocity analysis.<br />
    To visualize the velocity embeddings, you can use the `SeuratClusterStats`
    process with `v_reduction` provided to one of the `envs.dimplots`.<br />

## Environment Variables

- `ncores` *(`type=int`)*: *Default: `1`*. <br />
    Number of cores to use.<br />
- `group_by`:
    The column name in metadata to group the cells.<br />
    Typically, this column should be the cluster id.<br />
    If provided input is a Seurat object, the default identity will be used by
    default. Otherwise, it is recommended to provide this parameter.<br />
    "seurat_clusters" will be used with a warning if the input is in AnnData
    format and this parameter is not provided.<br />
- `mode` *(`type=list`)*: *Default: `['deterministic', 'stochastic', 'dynamical']`*. <br />
    The mode to use for the velocity analysis.<br />
    It should be a subset of `['deterministic', 'stochastic', 'dynamical']`,
    meaning that we can perform the velocity analysis in multiple modes.<br />
- `fitting_by` *(`choice`)*: *Default: `stochastic`*. <br />
    The mode to use for fitting the velocities.<br />
    - `stochastic`:
        Stochastic mode
    - `deterministic`:
        Deterministic mode
- `min_shared_counts` *(`type=int`)*: *Default: `30`*. <br />
    Minimum number of counts
    (both unspliced and spliced) required for a gene.<br />
- `n_neighbors` *(`type=int`)*: *Default: `30`*. <br />
    The number of neighbors to use for the velocity graph.<br />
- `n_pcs` *(`type=int`)*: *Default: `30`*. <br />
    The number of PCs to use for the velocity graph.<br />
- `denoise` *(`flag`)*: *Default: `False`*. <br />
    Whether to denoise the data.<br />
- `denoise_topn` *(`type=int`)*: *Default: `3`*. <br />
    Number of genes with highest likelihood selected to
    infer velocity directions.<br />
- `kinetics` *(`flag`)*: *Default: `False`*. <br />
    Whether to compute the RNA velocity kinetics.<br />
- `kinetics_topn` *(`type=int`)*: *Default: `100`*. <br />
    Number of genes with highest likelihood selected to
    infer velocity directions.<br />
- `calculate_velocity_genes` *(`flag`)*: *Default: `False`*. <br />
    Whether to calculate the velocity genes.<br />
- `top_n` *(`type=int`)*: *Default: `6`*. <br />
    The number of top features to plot.<br />
- `rscript`: *Default: `Rscript`*. <br />
    The path to the Rscript executable used to convert RDS file to AnnData.<br />
    if `in.sobjfile` is an RDS file, it will be converted to AnnData file
    (h5ad). You need `Seurat`, `SeuratDisk` and `digest` installed.<br />
- `outtype` *(`choice`)*: *Default: `<input>`*. <br />
    The output file type.<br />
    - `<input>`:
        The same as the input file type.<br />
    - `h5seurat`:
        h5seurat file
    - `h5ad`:
        h5ad file
    - `qs`:
        qs/qs2 file
    - `qs2`:
        qs2 file
    - `rds`:
        RDS file

