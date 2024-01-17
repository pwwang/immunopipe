# FAQ

/// details | `immunopipe` command not found?

Please make sure if you have installed `immunopipe` on the right `python`. If you have used `pip` to install `immunopipe`, make sure the `pip` is associated with the right `python`.

You may try `/path/to/python -m pip install -U immunopipe` to ensure `immunopipe` is installed with the `python` you wanted.

If `immunopipe` still can't be found from command line, try `/path/to/python -m immunopipe`.

///

//// details | Why I am getting "Error writing to connection:  No space left on device" while running [`ImmunarchLoading`](processes/ImmunarchLoading.md)?

If you are running the pipeline and it complains about "No space left on device" at `ImmunarchLoading`, and you are pretty sure that your working directory is way from full,  it is likely that your temporary directory does not have enough space. This is because that the `ImmunarchLoading` process will create a temporary directory to store the intermediate files, and the default temporary directory is `/tmp`. Make sure that you have enough space in `/tmp` or you can change the temporary directory by setting the environment variable of the process: [`envs.tmpdir`](processes/ImmunarchLoading.md#environment-variables).

It is also likely that you are running the pipeline in a docker container and the docker container does not have enough space in `/tmp`. In such case, you can try to run the pipeline with the `-v` option of docker to local directory to `/tmp` in the container. For example:

```shell
docker run --rm -w /workdir -v .:/workdir -v path/to/tmp:/tmp \
#                                         ^^^^^^^^^^^^^^^^^^^
    justold/immunopipe:<tag> @config.toml
```

If you are using `singularity`/`apptainer`, you can try to use the `-B` option to bind the local directory to `/tmp` in the container.

/// tab | Singularity
```
singularity run \
    --pwd /workdir -B .:/workdir -c -e -w \
    -B path/to/tmp:/tmp \
  # ^^^^^^^^^^^^^^^^^^^
    docker://justold/immunopipe:<tag> \
    @config.toml
```
///
/// tab | Apptainer
```
apptainer run \
    --pwd /workdir -B .:/workdir -c -e -w --unsquash \
    -B path/to/tmp:/tmp \
  # ^^^^^^^^^^^^^^^^^^^
    docker://justold/immunopipe:<tag> \
    @config.toml
```
///
////

//// details | Why does the pipeline stop at [`SeuratClusteringOfAllCells`](processes/SeuratClusteringOfAllCells.md) and family without a clear error message?

This is likely because that the pipeline is running out of memory. The [`SeuratClusteringOfAllCells`](processes/SeuratClusteringOfAllCells.md) and family processes (e.g. [`SeuratClustering`](processes/SeuratClustering.md)) will run the a series of `Seurat` functions to perform the clustering, especially the [`IntegrateData`](https://satijalab.org/seurat/reference/integratedata) and [`FindIntegrationAnchors`](https://satijalab.org/seurat/reference/findintegrationanchors) functions, and [`IntegrateLayers`](https://satijalab.org/seurat/reference/integratelayers) with `Seurat` v5.

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
- For `Seurat` v5, use corresponding parameters for [`IntegrateLayers`](https://satijalab.org/seurat/reference/integratelayers).

/// Tip
You can also pass a list of sample names instead of the sample indices. For example, `reference = ["sample1", "sample2"]` under `[SeuratPreparing.envs.IntegrateLayers]` to use `sample1` and `sample2` as the reference samples.

See also description about `IntegrateLayers` [here](processes/SeuratPreparing.md#environment-variables).
///

////

/// details | Can I run one of the processes from the pipeline separately if I have the input files prepared?

Only for some of the processes. `immunopipe` depends on [`biopipen`][1]. Most of the processes in `immunopipe` are subclasses of processes in `biopipen`. You can run the processes in `biopipen` separately by:

```shell
pipen run scrna SeuratClustering [options]
```

Note that only the processes from [`biopipen`][1] can be run separately. The processes in `immunopipe` are not designed to be run separately. For example, the `SeuratClusteringOfAllCells` process in `immunopipe` is a subclass of the `SeuratClustering` process in `biopipen`. It's specialized for the `immunopipe` pipeline. If you want to run a similar process separately, you should use the `SeuratClustering` process in `biopipen` instead.

Like `immunopipe`, you can also either provide a configuration file:

```shell
pipen run scrna SeuratClustering @config.toml
```

or specify the options in the command line:

```shell
pipen run scrna SeuratClustering --in.srtobj path/to/srtobj.RDS ...
```

You can also use the `-h`/`--help` option to see the brief options of the process, or use `-h+`/`--help+` to see the full options of the process.

///

/// details | How to run the pipeline on a cluster?
    attrs: {id: how-to-run-the-pipeline-on-a-cluster}

To run the pipeline on a cluster, it's recommended to install the pipeline locally so that the cluster nodes can access the pipeline.

The idea is that you have to make sure the file system can be accessed within the docker container on the nodes. One trick to use docker on a cluster is to mount the first level of the file system to the container. For example, if your data and environment are under `/research`, you can mount `/research` to `/research` in the container. Then you can use `/research` in the container to access the data and environment under `/research` on the host.

`immunopipe` is built on top of [`pipen`][2] and [`xqute`][3]. A set of schedulers are supported by default. These schedulers are:

- `local`: Run the pipeline locally.
- `slurm`: Run the pipeline on a slurm cluster.
- `sge`: Run the pipeline on a sge cluster.
- `ssh`: Run the pipeline on a remote host via ssh.

The scheduler can be specified via `scheduler_opts` for the whole pipeline or for a specific process. For example, to run the whole pipeline on a slurm cluster, you can use the following configuration file:

```toml
scheduler = "slurm"

[scheduler_opts]
sbatch_partition = "1-day"
```

To run a specific process on a slurm cluster, you can use the following configuration file:

```toml
[<Process>]
scheduler = "slurm"

[<Process>.scheduler_opts]
sbatch_partition = "1-day"
```

You can also use profiles to switch between different schedulers. See also <https://pwwang.github.io/pipen/configurations/#profiles>

///

/// details | Do I have to re-run the entire pipeline if I want to change some parameters?

If you want to change some parameters for a specific process, you just modify the configuration file and re-run the pipeline. The pipeline will detect the changes and re-run the necessary processes. For example, if you are changing some environment variables for [`ScFGSEA`](processes/ScFGSEA.md), the prior processes, such as the ones for clustering and differential expression analysis, will be cached and will not be re-run.

///

/// details | Why I am getting this error when running with [`apptainer`][4]: `FATAL:   no SIF writable overlay partition found in /tmp/apptainer_cache_xxx/...`?

You may need to add `--unsquash` option for `apptainer run`.

///

/// details | How can I use data with soft links while using docker image to run the pipeline?

The container does not have access to the host filesystem directly. You need to mount the directory containing the data to the container.

For example, if your real data is under `/path/to/data`, you can mount it to `/data` in the container (using `-v /path/to/data:/data` option for `docker` or `-B /path/to/data:/data` option for `singularity` or `apptainer`).

Then you can use `/data` in the container to access the data under `/path/to/data` on the host. Also remember to change the path of `RNAData` and `TCRData` in the file (e.g. `samples.txt`) that is passed to `SampleInfo` process.

Other than `/data`, there are other directories that you can use for mounting inside the container, including `/mnt` and `/tmp`, in case your want to mount multiple directories.

See also [`The directory structure in the container`](./installation.md#the-directory-structure-in-the-container).

///

/// details | Why I am getting `disk quota exceeded` error while pulling the docker image using `apptainer` with still plenty of space on the disk?

It's probably because that the cache directory of `apptainer` is full. You can try to use a different cache directory by setting the environment variable `APPTAINER_CACHEDIR` to a different directory. For example:

```shell
export APPTAINER_CACHEDIR=/path/to/cache
apptainer pull justold/immunopipe:<tag>
```

See also: <https://apptainer.org/docs/user/main/build_env.html#cache-folders>

///

/// details | `Unable to fork: Cannot allocate memory` or `long vectors not supported yet` during clustering using Seurat?

This is likely because that the pipeline is running out of memory. The [`SeuratClusteringOfAllCells`](processes/SeuratClusteringOfAllCells.md) and family processes (e.g. [`SeuratClustering`](processes/SeuratClustering.md)) will run the a series of `Seurat` functions to perform the clustering, especially the [`IntegrateData`](https://satijalab.org/seurat/reference/integratedata) and [`FindIntegrationAnchors`](https://satijalab.org/seurat/reference/findintegrationanchors) functions, and [`IntegrateLayers`](https://satijalab.org/seurat/reference/integratelayers) with `Seurat` v5.

You can try to set `envs.ncores` to a smaller number to reduce the memory usage. For example:

```toml
[SeuratClusteringOfAllCells.envs]
ncores = 4  # instead of 16
```

The other strategy is to use Reference-based integration `reference = [1, 2]` for [`IntegrateLayers`](https://satijalab.org/seurat/reference/findintegrationanchors) with method `rpca` or `cca`. See also description about `IntegrateLayers` [here](processes/SeuratPreparing.md#environment-variables). For example:

```toml
[SeuratPreparing.envs.IntegrateLayers]
method = "rpca"
reference = [1, 2]  # You can also use sample names instead of indices
```

See also these issues for more details:

- <https://github.com/satijalab/seurat/issues/1029>
- <https://github.com/satijalab/seurat/issues/7419>

///

/// details | Got error `Not all stats values are finite numbers` while running [`ScFGSEA`](processes/ScFGSEA.md)?

It's probably because that there are too many missing values in the expression matrix and `signal_to_noise` is not able to detect the rank of the genes. You can try a different method for gene preranking. For example:

```toml
[ScFGSEA.envs]
method = "diff_of_classes"
```

///

<p> </p>

[1]: https://github.com/pwwang/biopipen
[2]: https://github.com/pwwang/pipen
[3]: https://github.com/pwwang/xqute
[4]: https://apptainer.org/
