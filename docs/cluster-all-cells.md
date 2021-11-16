Cell clustering is done by `biopipen.namespaces.scrna.SeuratClustering` using `Seurat::FindClusters()`. See the full configurations by:

```shell
pipen run scrna SeuratClustering
```

You can specify some arguments to `Seurat::FindClusters()` by `envs.FindClusters`. For example:

```toml
[SeuratClustering.envs]
FindClusters = { resolution = 0.8 }
```

After the clustering, we need to separate T cells and non-T cells, since not all cells from scRNA seq data belong to a T-cell clone. To do that, we use two metrics to separate the cells. Other than clonotype percentage, we also examine the mean expression of CD3E [Ref: Wu, Thomas D., et al. "Peripheral T cell expansion predicts tumour infiltration and clinical response." Nature 579.7798 (2020): 274-278.].

To setup the cutoff to separate T and non-T cells, you can pass some filters:

```toml
[SelectTCells.envs]
tcell_filter = "Clonotype_pct > 0.25 & CD3E_means > mean(CD3E_means)"
```

To futher help determine the cell types, the biomarkers for each cluster are found by `Seurat::FindMarkers()` and enrichment analysis is also done for those markers. To specify the libraries for enrichment analysis:

```toml
[MarkersForClustersOfAllCells.envs]
dbs = ["KEGG_2021_Human", "GO_Biological_Process_2021"]
```

See all available libraries here:
[https://maayanlab.cloud/Enrichr/#libraries](https://maayanlab.cloud/Enrichr/#libraries)
