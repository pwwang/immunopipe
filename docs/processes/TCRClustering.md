# TCRClustering

Cluster the TCR clones by their CDR3 sequences

You can disable this by remving the whole sections of
TCRClustering and TCRClusterStats in the config file.<br />

This process is used to cluster TCR clones based on their CDR3 sequences.<br />

It uses either

[GIANA](https://github.com/s175573/GIANA)

> Zhang, Hongyi, Xiaowei Zhan, and Bo Li.<br />
> "GIANA allows computationally-efficient TCR clustering and multi-disease
> repertoire classification by isometric transformation."
> Nature communications 12.1 (2021): 1-11.<br />

Or [ClusTCR](https://github.com/svalkiers/clusTCR)

> Sebastiaan Valkiers, Max Van Houcke, Kris Laukens, Pieter Meysman,
> ClusTCR: a Python interface for rapid clustering of large sets of CDR3
> sequences with unknown antigen specificity,
> Bioinformatics, 2021.<br />

Both methods are based on the
[Faiss Clustering Library](https://github.com/facebookresearch/faiss),
for efficient similarity search and clustering of dense vectors, so both methods
yield similar results.<br />

A text file will be generated with the cluster assignments for each cell, together
with the `immunarch` object (in `R`) with the cluster assignments at `TCR_Clsuter`
column. This information will then be merged to a `Seurat` object for further
downstream analysis.<br />

The cluster assignments are prefixed with `S_` or `M_` to indicate whether a
cluster has only one unique CDR3 sequence or multiple CDR3 sequences.<br />
Note that a cluster with `S_` prefix may still have multiple cells, as the same
CDR3 sequence may be shared by multiple cells.<br />

## Environment Variables

- `tool` *(`choice`)*: *Default: `GIANA`*. <br />
    The tool used to do the clustering, either
    [GIANA](https://github.com/s175573/GIANA) or
    [ClusTCR](https://github.com/svalkiers/clusTCR).<br />
    For GIANA, using TRBV mutations is not supported
    - `GIANA`:
        by Li lab at UT Southwestern Medical Center
    - `ClusTCR`:
        by Sebastiaan Valkiers, etc
- `prefix`:
    The prefix to the barcodes. You can use placeholder like `{Sample}_`
    The prefixed barcodes will be used to match the barcodes in `in.metafile`.<br />
    Not used if `in.metafile` is not specified.<br />
    If `None` (default), `immdata$prefix` will be used.<br />
- `python`: *Default: `python`*. <br />
    The path of python with `GIANA`'s dependencies installed
    or with `clusTCR` installed. Depending on the `tool` you choose.<br />
- `args` *(`type=json`)*: *Default: `{}`*. <br />
    The arguments for the clustering tool
    For GIANA, they will be passed to `python GIAna.py`
    See <https://github.com/s175573/GIANA#usage>.<br />
    For ClusTCR, they will be passed to `clustcr.Clustering(...)`
    See <https://svalkiers.github.io/clusTCR/docs/clustering/how-to-use.html#clustering>.<br />

