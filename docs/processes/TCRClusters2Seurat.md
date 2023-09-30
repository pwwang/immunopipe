# TCRClusters2Seurat

This process is used to merge the cluster assignments from [TCRClustering](./TCRClustering.md) to the `Seurat` object. The cluster assignments are prefixed with `S_` or `M_` to indicate whether a cluster has only one unique CDR3 sequence or multiple CDR3 sequences. Note that a cluster with `S_` prefix may still have multiple cells, as the same CDR3 sequence may be shared by multiple cells. The cluster assignments are saved in the `Seurat` object at `TCR_Cluster` column in `seurat_object@meta.data` in `R`.

Other two columns are also added to the `Seurat` object: `TCR_Cluster_Size` and `TCR_Cluster_Size1`. The `TCR_Cluster_Size` column contains the number of cells in each cluster, while the `TCR_Cluster_Size1` column contains the number of unique CDR3 sequences in each cluster.

/// Tip | New in `0.7.0`
`TCR_Cluster_Size` and `TCR_Cluster_Size1` are added in `0.7.0`.
///

Those columns can be then used for further downstream analysis. For example, you can find the markers for the TCR cluster (i.e. `S_1` vs `S_2`) in each seurat cluster by:

```toml
[MarkersFinder.envs]
group-by = "TCR_Cluster"
ident-1 = "S_1"
ident-2 = "S_2"
each = "seurat_clusters"
```

There is no environment variables for this process.
