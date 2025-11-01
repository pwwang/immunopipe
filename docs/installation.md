# Installation

## Install the pipline and the dependencies using conda

/// Tip
If you plan to use the docker image to run the pipeline locally, you can skip this section.
///

`immunopipe` is built upon [`pipen`](https://github.com/pwwang/pipen) framework, and a number of packages written in `R` and `python`. It's not recommended to install the packages manually. Instead, you can use the provided `environment.yml` to create a conda environment.

```shell
$ conda env create \
    -n immunopipe \
    -f https://raw.githubusercontent.com/pwwang/immunopipe/master/docker/environment.yml
```

/// Attention
The `environment.yml` includes only a subset of the packages required by the pipeline. The dependencies of the [`TESSA`](processes/TESSA.md) process are not included.

If you want to enable [`TESSA`](processes/TESSA.md) process , please use the `environment_full.yml` (<https://raw.githubusercontent.com/pwwang/immunopipe/master/docker/environment_full.yml>) instead.

See also [Choose the right tag of the docker image](running.md#choose-the-right-tag-of-the-docker-image).
///

If the URL doesn't work, you can download the file and create the environment locally.

For more detailed instructions of `conda env create`, please refer to [conda docs](https://docs.conda.io/projects/conda/en/latest/commands/env/create.html).

/// Attention
The pipeline itself is NOT included in the conda environment. You need to install it separately.

```shell
$ conda activate immunopipe
$ pip install -U immunopipe
$ # If you want to create diagram and generate running information
$ pip install -U immunopipe[diagram,runinfo]
$ # You also need to install the frontend dependencies to generate reports
$ pipen report update
```
///

## Use the docker image

You can also use the docker image to run the pipeline. The image is built upon `miniconda3` and `micromamba` is used as the package manager. The image is available at [Docker Hub](https://hub.docker.com/r/justold/immunopipe).

To pull the image:

/// tab | Using docker
```shell
$ docker pull justold/immunopipe:<tag>
```
///

/// tab | Using singularity
If you are using `singularity`, you can pull and convert the image to `sif` format:

```shell
$ singularity pull docker://justold/immunopipe:<tag>
```
///

/// tab | Using apptainer
```shell
$ apptainer pull docker://justold/immunopipe:<tag>
```
///

To run the pipeline use the image, please refer to [Running the pipeline](./running.md).


### The directory structure in the container

The docker image is build upon [`mambaorg/micromamba:1.4.3`][1]. The OS is linux/amd64. Other than the default directories, the following directories are also created or should be mapped during the run:

- `/immunopipe`: The directory where the source code of the pipeline is. It is general a clone of the [repository][2]. The pipeline is also installed from this directory.
- `/workdir`: The working directory. It is the directory where the pipeline is run. It is recommended to map the current directory (`.`) to this directory.
- `/data`: An empty directory. You can map your data directory to this directory.

## Prepare to run the pipeline via Google Batch Jobs

There are two ways of running the pipeline via Google Batch Jobs: using the `gbatch` scheduler (provided by [`xqute`][3]) or using [`pipen-cli-gbatch`][4]. See more details in [Running the pipeline via Google Batch Jobs](./running.md#run-the-pipeline-using-google-cloud-batch-jobs).

In addition to prepare the docker image in the artifact registry (or docker hub if your google cloud project allows pulling from docker hub), you also need to install some dependencies locally.

If you choose to use the `gbatch` scheduler, in addition to installing the pipeline:

```shell
$ pip install -U immunopipe

# install cloud dependencies
$ pip install -U cloudpathlib[gs]
```

You still have to install the following dependencies to generate reports:

- [`nodejs`](https://nodejs.org/): Follow the instructions at <https://nodejs.org/en/download/package-manager/> to install `nodejs` for your system (v20+ is required); or
- [`bunjs`](https://bun.sh/): Follow the instructions at <https://bun.sh/docs/install> to install `bunjs` for your system.

Then you need to install frontend dependencies for report generation:

```shell
$ pipen report update
```

If you choose to use [`pipen-cli-gbatch`][4] (running the pipeline via `immunopipe gbatch`), you just need to install the pipeline with the `cli-gbatch` extra:

```shell
$ pip install -U immunopipe[cli-gbatch]
```

[1]: https://hub.docker.com/layers/mambaorg/micromamba/1.4.3/images/sha256-0251b94151c021c85d3e4f4ffe1fc81c436f18e01337d3b367d0f7c76ee716ac?context=explore
[2]: https://github.com/pwwang/immunopipe
[3]: https://github.com/pwwang/xqute
[4]: https://github.com/pwwang/pipen-cli-gbatch
