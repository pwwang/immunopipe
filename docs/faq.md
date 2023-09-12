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

<p> </p>

[1]: https://github.com/pwwang/biopipen
[2]: https://github.com/pwwang/pipen
[3]: https://github.com/pwwang/xqute