# TCRClustering

This process is used to cluster TCR clones based on their CDR3 sequences.

It uses either:

[GIANA](https://github.com/s175573/GIANA)

> Zhang, Hongyi, Xiaowei Zhan, and Bo Li.
> "GIANA allows computationally-efficient TCR clustering and multi-disease
> repertoire classification by isometric transformation."
> Nature communications 12.1 (2021): 1-11.

Or [ClusTCR](https://github.com/svalkiers/clusTCR)

> Sebastiaan Valkiers, Max Van Houcke, Kris Laukens, Pieter Meysman,
> ClusTCR: a Python interface for rapid clustering of large sets of CDR3
> sequences with unknown antigen specificity,
> Bioinformatics, 2021.

Both methods are based on the [Faiss Clustering Library](https://github.com/facebookresearch/faiss), for efficient similarity search and clustering of dense vectors, so both methods yield similar results.

A text file will be generated with the cluster assignments for each cell, together with the `immunarch` object (in `R`) with the cluster assignments at `TCR_Clsuter` column. This information will then be merged to the `Seurat` object by [TCRClusters2Seurat](./TCRClusters2Seurat.md). Futher downstream analysis can be performed using the cluster assignments.

The cluster assignments are prefixed with `S_` or `M_` to indicate whether a cluster has only one unique CDR3 sequence or multiple CDR3 sequences. Note that a cluster with `S_` prefix may still have multiple cells, as the same CDR3 sequence may be shared by multiple cells.

## Environment variables

- `tool` (`choice`): The tool used to do the clustering, either
    [GIANA](https://github.com/s175573/GIANA) or
    [ClusTCR](https://github.com/svalkiers/clusTCR).
    For GIANA, using TRBV mutations is not supported
    - `GIANA`: by Li lab at UT Southwestern Medical Center
    - `ClusTCR`: by Sebastiaan Valkiers, etc
- `python`: The path of python with `GIANA`'s dependencies installed
    or with `clusTCR` installed. Depending on the `tool` you choose.
- `args` (`type=json`): The arguments for the clustering tool
    For GIANA, they will be passed to `python GIAna.py`
    See <https://github.com/s175573/GIANA#usage>.
    For ClusTCR, they will be passed to `clustcr.Clustering(...)`
    See <https://svalkiers.github.io/clusTCR/docs/clustering/how-to-use.html#clustering>.
- `on_multi` (`flag`): Whether to run clustering on multi-chain seq or the seq read and processed by `immunarch`. Currently, `immunarch` is only able to process single-chain seq, and the data with multi-chain seq is saved in `immdata$multi` in `R`. This option is kept for future use.
