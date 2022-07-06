## Clustering

Cluster the TCR clones by their CDR3 sequences

With GIANA

https://github.com/s175573/GIANA

> Zhang, Hongyi, Xiaowei Zhan, and Bo Li. "GIANA allows computationally-efficient TCR clustering and multi-disease repertoire classification by isometric transformation." Nature communications 12.1 (2021): 1-11.

Or ClusTCR

https://github.com/svalkiers/clusTCR

> Sebastiaan Valkiers, Max Van Houcke, Kris Laukens, Pieter Meysman, ClusTCR: a Python interface for rapid clustering of large sets of CDR3 sequences with unknown antigen specificity, Bioinformatics, 2021.


## Statistics

- `shared_clusters` — Stats about shared TCR clusters
    - `numbers_on_heatmap`: Whether to show the numbers on the heatmap
    - `heatmap_meta`: The metadata to show on the heatmap
    - `grouping`: The groups to investigate the shared clusters

- `sample_diversity` — Sample diversity using TCR clusters instead of clones keys are the methods and values, currently, by to plot the diversities by groups

## Configurations

To select the tool for clustering:

```toml
[TCRClustering.envs]
# # The tool used to do the clustering, either GIANA or ClusTCR
tool = "ClusTCR"
# # Arguments for GIANA/ClusTCR
# # For GIANA, they will be passed to `python GIAna.py <args>`
# # For ClusTCR, they will be passed to `clustcr.Clustering(<args>)`
# args = {}
```

To control the statistics:

```toml
[TCRClusteringStats.envs]
# # numbers_on_heatmap: Whether to show the numbers on the heatmap
# # heatmap_meta: The metadata to show on the heatmap
# # grouping: The groups to investigate the shared clusters
shared_clusters = { numbers_on_heatmap = true, heatmap_meta = [], grouping = "null" }
# # Sample diversity using TCR clusters instead of clones
# # keys are the methods and values, currently, `by` to plot
# # the diversities by groups
sample_diversity = { gini = { by = ["region"] } }
```
