# MetabolicExprImputation

This process imputes the dropout values in scRNA-seq data.

It takes the Seurat object as input and outputs the Seurat object with
imputed expression data.<br />

You can turn off the imputation by setting the `noimpute` option
of the process group to `True`.<br />

## Input

- `infile`:
    The input file in RDS/qs format of Seurat object

## Output

- `outfile`: *Default: `{{in.infile | stem}}.imputed.qs`*. <br />
    The output file in RDS format of Seurat object
    Note that with rmagic and alra, the original default assay will be
    renamed to `RAW` and the imputed RNA assay will be
    renamed to `RNA` and set as default assay.<br />

## Environment Variables

- `tool` *(`choice`)*: *Default: `alra`*. <br />
    Either alra, scimpute or rmagic
    - `alra`:
        Use RunALRA() from Seurat
    - `scimpute`:
        Use scImpute() from scimpute
    - `rmagic`:
        Use magic() from Rmagic
- `scimpute_args` *(`ns`)*:
    The arguments for scimpute
    - `drop_thre` *(`type=float`)*: *Default: `0.5`*. <br />
        The dropout threshold
    - `kcluster` *(`type=int`)*:
        Number of clusters to use
    - `ncores` *(`type=int`)*: *Default: `1`*. <br />
        Number of cores to use
    - `refgene`: *Default: `""`*. <br />
        The reference gene file
- `rmagic_args` *(`ns`)*:
    The arguments for rmagic
    - `python`: *Default: `python`*. <br />
        The python path where magic-impute is installed.<br />
    - `threshold` *(`type=float`)*: *Default: `0.5`*. <br />
        The threshold for magic imputation.<br />
        Only the genes with dropout rates greater than this threshold (No. of
        cells with non-zero expression / total number of cells) will be imputed.<br />
- `alra_args` *(`type=json`)*: *Default: `{}`*. <br />
    The arguments for `RunALRA()`

## Reference

- [Linderman, George C., Jun Zhao, and Yuval Kluger. "Zero-preserving imputation of scRNA-seq data using low-rank approximation." BioRxiv (2018): 397588.](https://www.nature.com/articles/s41467-021-27729-z)
- [Li, Wei Vivian, and Jingyi Jessica Li. "An accurate and robust imputation method scImpute for single-cell RNA-seq data." Nature communications 9.1 (2018): 997.](https://www.nature.com/articles/s41467-018-03405-7)
- [Dijk, David van, et al. "MAGIC: A diffusion-based imputation method reveals gene-gene interactions in single-cell RNA-sequencing data." BioRxiv (2017): 111591.](https://www.cell.com/cell/abstract/S0092-8674(18)30724-4)

