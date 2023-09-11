# FAQ

/// details | Why I am getting "Error writing to connection:  No space left on device" while running [`ImmunarchLoading`](processes/ImmunarchLoading.md)?

If you are running the pipeline and it complains about "No space left on device" at `ImmunarchLoading`, and you are pretty sure that your working directory is way from full, then  it is likely that your temporary directory does not have enough space. This is because that the `ImmunarchLoading` process will create a temporary directory to store the intermediate files, and the default temporary directory is `/tmp`. Make sure that you have enough space in `/tmp` or you can change the temporary directory by setting the environment variable of the process: [`envs.tmpdir`](processes/ImmunarchLoading.md#environment-variables).

It is also likely that you are running the pipeline in a docker container and the docker container does not have enough space in `/tmp`. In such case, you can try to run the pipeline with the `-v` option of docker to local directory to `/tmp` in the container. For example:

```shell
docker run --rm -w /workdir -v .:/workdir -v path/to/tmp:/tmp \
#                                         ^^^^^^^^^^^^^^^^^^^
    justold/immunopipe:<tag> @config.toml
```

If you are using singularity, you can try to use the `-B` option to bind the local directory to `/tmp` in the container.

///

/// details | Why does the pipeline stop at [`SeuratClusteringOfAllCells`](processes/SeuratClusteringOfAllCells.md) and family without a clear error message?

This is likely because that the pipeline is running out of memory. The [`SeuratClusteringOfAllCells`](processes/SeuratClusteringOfAllCells.md) and family processes (e.g. [`SeuratClusteringOfTCells`](processes/SeuratClusteringOfTCells.md)) will run the a series of `Seurat` functions to perform the clustering, especially the [`IntegrateData`](https://satijalab.org/seurat/reference/integratedata) and [`FindIntegrationAnchors`](https://satijalab.org/seurat/reference/findintegrationanchors) functions.

Please see the following issues for more details:

- <https://github.com/satijalab/seurat/issues/3326>
- <https://github.com/satijalab/seurat/issues/1720>
- <https://github.com/satijalab/seurat/issues/2828>
- <https://github.com/satijalab/seurat/issues/1254>
- <https://github.com/satijalab/seurat/issues/7027>

Also check out the tips by the Seurat team:

<https://satijalab.org/seurat/articles/integration_large_datasets>

Two possible solutions are:

- Use `reduction = "rpca"` for [`FindIntegrationAnchors`](https://satijalab.org/seurat/reference/findintegrationanchors) under `[SeuratClusteringOfAllCells.envs.FindIntegrationAnchors]`.
- Use Reference-based integration `reference = [1, 2]` for [`FindIntegrationAnchors`](https://satijalab.org/seurat/reference/findintegrationanchors) under `[SeuratClusteringOfAllCells.envs.FindIntegrationAnchors]`.

!!! Tip
    You can also pass a list of sample names instead of the sample indices. For example, `reference = ["sample1", "sample2"]` under `[SeuratClusteringOfAllCells.envs.FindIntegrationAnchors]` to use `sample1` and `sample2` as the reference samples.

    See also description about `FindIntegrationAnchors` [here](processes/SeuratClusteringOfAllCells.md#environment-variables).

///

<p> </p>
