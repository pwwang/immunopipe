# MetabolicExprImpution

This process imputes the dropout values in scRNA-seq data. It takes the Seurat object as input and outputs the Seurat object with imputed expression data. The Seurat object will be saved as `RDS` file and passed to the next processes in the [`ScrnaMetabolicLandscape`](./ScrnaMetabolicLandscape.md) group for further analysis.

You can also use the `noimpute` argument from [`ScrnaMetabolicLandscape`](./ScrnaMetabolicLandscape.md#group-arguments) to skip this process.

## Environment variables

- `tool` (`choice`): Either alra, scimpute or rmagic
    - `alra`: Use [RunALRA()][1] from Seurat
    - `scimpute`: Use [scImpute()][2] from scimpute
    - `rmagic`: Use [magic()][3] from Rmagic
- `scimpute_args` (`ns`): The arguments for [`scimpute`][4]
    - `drop_thre` (`type=float`): The dropout threshold
    - `kcluster` (`type=int`): Number of clusters to use
    - `ncores` (`type=int`): Number of cores to use
    - `refgene`: The reference gene file
- `rmagic_args` (`ns`): The arguments for `rmagic`
    - `python`: The python path where [`magic-impute`][6] is installed.
- `alra_args` (`type=json`): The arguments for [`RunALRA()`][5]

## Reference

- [Linderman, George C., Jun Zhao, and Yuval Kluger. "Zero-preserving imputation of scRNA-seq data using low-rank approximation." BioRxiv (2018): 397588.][7]
- [Li, Wei Vivian, and Jingyi Jessica Li. "An accurate and robust imputation method scImpute for single-cell RNA-seq data." Nature communications 9.1 (2018): 997.][8]
- [Dijk, David van, et al. "MAGIC: A diffusion-based imputation method reveals gene-gene interactions in single-cell RNA-sequencing data." BioRxiv (2017): 111591.][9]

[1]: https://github.com/satijalab/seurat-wrappers/blob/master/docs/alra.md
[2]: https://www.nature.com/articles/s41467-018-03405-7
[3]: https://github.com/cran/Rmagic
[4]: https://rdrr.io/github/Vivianstats/scImpute/man/scimpute.html
[5]: https://rdrr.io/github/satijalab/seurat-wrappers/man/RunALRA.html
[6]: https://pypi.org/project/magic-impute/
[7]: https://www.nature.com/articles/s41467-021-27729-z
[8]: https://www.nature.com/articles/s41467-018-03405-7
[9]: https://www.cell.com/cell/abstract/S0092-8674(18)30724-4
