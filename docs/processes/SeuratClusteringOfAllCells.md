# SeuratClusteringOfAllCells

This process will perform clustering on all cells using [`Seurat`][1] package. The clusters will then be used to select T cells by [`TCellSelection`](TCellSelection.md) process. The clustering results will also be used to find markers for each cluster by [`MarkersForClustersOfAllCells`](MarkersForClustersOfAllCells.md) process, which is useful to examine the cell type composition of the samples.

!!! note
    If all you cells are T cells, this process will perform clustering on all cells and the results will be used on downstream T-cell analyses and the integrative analyses. At the same time, you should leave the [`TCellSelection`](TCellSelection.md) process out of the pipeline, by not setting anything for the process in the configuration file.

To perform the clustering, you have two routes to choose from:

1. Performing integration on datasets normalized with `SCTransform`
    - See: [https://satijalab.org/seurat/articles/integration_rpca.html#performing-integration-on-datasets-normalized-with-sctransform-1](https://satijalab.org/seurat/articles/integration_rpca.html#performing-integration-on-datasets-normalized-with-sctransform-1)
2. Fast integration using reciprocal PCA (`RPCA`)
    - See: [https://satijalab.org/seurat/articles/integration_rpca.html](https://satijalab.org/seurat/articles/integration_rpca.html)

For the first route, following procedures will be performed in the order:

- [`SplitObject`][2]: Split the `Seurat` object into a list of `Seurat` objects, one for each sample.
- [`SCTransform*`][3]: Normalize the data using `Seurat::SCTransform()`.
- [`SelectIntegrationFeatures`][4]: Select integration features.
- [`PrepSCTIntegration`][5]: Prepare the integration.
- [`RunPCA*`][6]: Perform PCA on each sample individually.
- [`FindIntegrationAnchors`][7]: Find integration anchors.
- [`IntegrateData`][8]: Integrate the data.
- [`RunPCA`][6]: Perform PCA on the integrated data.
- [`RunUMAP`][9]: Perform UMAP on the integrated data.
- [`FindNeighbors`][10]: Find neighbors.
- [`FindClusters`][11]: Perform clustering.
- `*`: These steps are performed on each sample individually.

For the second route, following procedures will be performed in the order:

- [`SplitObject`][2]: Split the `Seurat` object into a list of `Seurat` objects, one for each sample.
- [`NormalizeData*`][12]: Normalize the data using `Seurat::NormalizeData()`.
- [`FindVariableFeatures*`][13]: Find variable features.
- [`SelectIntegrationFeatures`][4]: Select integration features.
- [`ScaleData*`][14]: Scale the data.
- [`RunPCA*`][6]: Perform PCA on each sample individually.
- [`FindIntegrationAnchors`][7]: Find integration anchors.
- [`IntegrateData`][8]: Integrate the data.
- [`ScaleData`][14]: Scale the integrated data.
- [`RunPCA`][6]: Perform PCA on the integrated data.
- [`RunUMAP`][9]: Perform UMAP on the integrated data.
- [`FindNeighbors`][10]: Find neighbors.
- [`FindClusters`][11]: Perform clustering.
- `*`: These steps are performed on each sample individually.

## Environment variables

The variables marked with (`ns`) are `namespace` variables. See also: [Namespace environment variables](../configurations.md/#namespace-environment-variables).

- `ncores`: Number of cores to use for parallel computing for some `Seurat` procedures.
    - Used in `future::plan(strategy = "multicore", workers = <ncores>)` to parallelize some Seurat procedures.
    - See also: <https://satijalab.org/seurat/articles/future_vignette.html>
- `use_sct`: Whether use `SCTransform` route (the first route mentioned above) or not (default: `False`, using `RPCA`).
- `SCTransform` (`ns`): Arguments for [`SCTransform()`](https://satijalab.org/seurat/reference/sctransform).
    `object` is specified internally, and `-` in the key will be replaced with `.`.
    - `<more>`: See more arguments: <https://satijalab.org/seurat/reference/sctransform>
- `SelectIntegrationFeatures` (`ns`): Arguments for [`SelectIntegrationFeatures()`](https://satijalab.org/seurat/reference/selectintegrationfeatures).
    `object.list` is specified internally, and `-` in the key will be replaced with `.`.
    - `nfeatures` (`type=int`): The number of features to select
    - `<more>`: See more arguments: <https://satijalab.org/seurat/reference/selectintegrationfeatures>
- `PrepSCTIntegration` (`ns`): Arguments for [`PrepSCTIntegration()`](https://satijalab.org/seurat/reference/prepsctintegration).
    `object.list` and `anchor.features` is specified internally, and `-` in the key will be replaced with `.`.
    - `<more>`: See more arguments: <https://satijalab.org/seurat/reference/prepsctintegration>
- `NormalizeData` (`ns`): Arguments for [`NormalizeData()`](https://satijalab.org/seurat/reference/normalizedata).
    `object` is specified internally, and `-` in the key will be replaced with `.`.
    - `<more>`: See more arguments: <https://satijalab.org/seurat/reference/normalizedata>
- `FindVariableFeatures` (`ns`): Arguments for [`FindVariableFeatures()`](https://satijalab.org/seurat/reference/findvariablefeatures).
    `object` is specified internally, and `-` in the key will be replaced with `.`.
    - `<more>`: See more arguments: <https://satijalab.org/seurat/reference/findvariablefeatures>
- `FindIntegrationAnchors` (`ns`): Arguments for [`FindIntegrationAnchors()`](https://satijalab.org/seurat/reference/findintegrationanchors).
    `object.list` and `anchor.features` is specified internally, and `-` in the key will be replaced with `.`.
    `dims=N` will be expanded to `dims=1:N`; The maximal value of `N` will be the minimum of `N` and the number of columns for each sample.
    Sample names can also be specified in `reference` instead of indices only.
    `reduction` defaults to `rpca`.
    `normalization.method` defaults to `SCT` if `use_sct` is `True`.
    **If you want to use reference-based integration, you can also set `reference` to a list of sample names, instead of a list of indices.**
    - `<more>`: See more arguments: <https://satijalab.org/seurat/reference/findintegrationanchors>
- `IntegrateData` (`ns`): Arguments for [`IntegrateData()`](https://satijalab.org/seurat/reference/integratedata).
    `anchorset` is specified internally, and `-` in the key will be replaced with `.`.
    `dims=N` will be expanded to `dims=1:N`; The maximal value of `N` will be the minimum of `N` and the number of columns for each sample.
    `normalization.method` defaults to `SCT` if `use_sct` is `True`.
    - <more>: See <https://satijalab.org/seurat/reference/integratedata>
- `ScaleData` (`ns`): Arguments for [`ScaleData()`](https://satijalab.org/seurat/reference/scaledata).
    `object` and `features` is specified internally, and `-` in the key will be replaced with `.`.
    - `verbose` (flag): Whether to print the progress
    - `<more>`: See more arguments: <https://satijalab.org/seurat/reference/scaledata>
- `ScaleData1` (`ns`): Arguments for [`ScaleData()`](https://satijalab.org/seurat/reference/scaledata) on each sample.
    `object` and `features` is specified internally, and `-` in the key will be replaced with `.`.
    - `verbose` (flag): Whether to print the progress
    - `<more>`: See more arguments: <https://satijalab.org/seurat/reference/scaledata>
- `RunPCA` (`ns`): Arguments for [`RunPCA()`](https://satijalab.org/seurat/reference/runpca).
    `object` and `features` is specified internally, and `-` in the key will be replaced with `.`.
    - `npcs` (`type=int`): The number of PCs to compute.
        For each sample, `npcs` will be no larger than the number of columns - 1.
    - `verbose` (flag): Whether to print the progress
    - `<more>`: See more arguments: <https://satijalab.org/seurat/reference/runpca>
- `RunPCA1` (`ns`): Arguments for [`RunPCA()`](https://satijalab.org/seurat/reference/runpca) on each sample.
    `object` and `features` is specified internally, and `-` in the key will be replaced with `.`.
    - `npcs` (`type=int`): The number of PCs to compute.
        For each sample, `npcs` will be no larger than the number of columns - 1.
    - `verbose` (flag): Whether to print the progress
    - `<more>`: See more arguments: <https://satijalab.org/seurat/reference/runpca>
- `RunUMAP` (`ns`): Arguments for [`RunUMAP()`](https://satijalab.org/seurat/reference/runumap).
    `object` is specified internally, and `-` in the key will be replaced with `.`.
    `dims=N` will be expanded to `dims=1:N`; The maximal value of `N` will be the minimum of `N` and the number of columns - 1 for each sample.
    - `dims` (`type=int`): The number of PCs to use
    - `reduction`: The reduction to use for UMAP
    - `<more>`: See more arguments: <https://satijalab.org/seurat/reference/runumap>
- `FindNeighbors` (`ns`): Arguments for [`FindNeighbors()`](https://satijalab.org/seurat/reference/findneighbors).
    `object` is specified internally, and `-` in the key will be replaced with `.`.
    - `<more>`: See more arguments: <https://satijalab.org/seurat/reference/findneighbors>
- `FindClusters` (`ns`): Arguments for [`FindClusters()`](https://satijalab.org/seurat/reference/findclusters).
    `object` is specified internally, and `-` in the key will be replaced with `.`.
    - `resolution` (`type=float`): The resolution of the clustering
    - `<more>`: See more arguments: <https://satijalab.org/seurat/reference/findclusters>

[1]: https://satijalab.org/seurat/
[2]: https://satijalab.org/seurat/reference/splitobject
[3]: https://satijalab.org/seurat/reference/sctransform
[4]: https://satijalab.org/seurat/reference/selectintegrationfeatures
[5]: https://satijalab.org/seurat/reference/prepsctintegration
[6]: https://satijalab.org/seurat/reference/runpca
[7]: https://satijalab.org/seurat/reference/findintegrationanchors
[8]: https://satijalab.org/seurat/reference/integratedata
[9]: https://satijalab.org/seurat/reference/runumap
[10]: https://satijalab.org/seurat/reference/findneighbors
[11]: https://satijalab.org/seurat/reference/findclusters
[12]: https://satijalab.org/seurat/reference/normalizedata
[13]: https://satijalab.org/seurat/reference/findvariablefeatures
[14]: https://satijalab.org/seurat/reference/scaledata
