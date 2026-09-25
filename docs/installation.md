# Installation

## Install the pipline and the dependencies using conda/mamba

/// Attention
**The conda/mamba route below is Linux-only.** The R packages immunopipe depends on
(`bioconductor-screpertoire`, `r-biopipen.utils`, `r-hitype`, `r-plotthis`,
`r-scplotter`, `r-seuratdisk`, `r-seuratwrappers`) are published on the `pwwang`
channel as `linux-64` builds only. On macOS the solve therefore fails with a list of
packages that "does not exist (perhaps a typo or a missing channel)", regardless of
whether the machine is Intel or Apple silicon. On macOS, install with `pip` and provide
`R` separately, or use the [docker image](#use-the-docker-image).

This is verified continuously: `.github/workflows/install-matrix.yml` performs the
documented installation and then runs the pipeline on a minimal dataset on both
`ubuntu-latest` and `macos-latest`.
///

/// Tip
If you plan to use the docker image to run the pipeline locally, you can skip this section.
///

`immunopipe` is built upon [`pipen`](https://github.com/pwwang/pipen) framework, and a number of packages written in `R` and `python`. It's not recommended to install the packages manually. Instead, you can use the provided `environment_base.yml` to create a conda environment. Download the environment files first:

```shell
$ curl -fLO https://raw.githubusercontent.com/pwwang/immunopipe/master/docker/environment_base.yml
$ curl -fLO https://raw.githubusercontent.com/pwwang/immunopipe/master/docker/environment_rpkgs.yml
```

Then create the environment with the base file:

```shell
$ conda env create \
    -n immunopipe \
    -f ./environment_base.yml
```

And update it with the essential `R` packages:

```shell
$ conda env update \
    -n immunopipe \
    -f ./environment_rpkgs.yml
```

/// Attention
Do not pass the environment file as a URL (that is, do not use `-f https://...`). In `conda env create`, `-f` accepts a *list* of paths; when the argument is a URL, conda fails in its own pip step with `AttributeError: 'list' object has no attribute 'split'` (see `conda/env/pip_util.py`). The dependency solve succeeds and every package is installed before that failure, so the environment is left **silently incomplete** — the `pip:` entries (`tensorflow`, `scvelo`, `liana`, `keras`, ...) are missing, and the error only appears at the very end of a long installation. Downloading the files first, as above, avoids this. Verified with conda 26.7.2.
///

For more detailed instructions of `conda env create`, please refer to [conda docs](https://docs.conda.io/projects/conda/en/latest/commands/env/create.html).

/// Attention
The pipeline itself is NOT included in the conda environment. You need to install it separately.

```shell
$ conda activate immunopipe
$ pip install -U immunopipe
$ # If you want to create diagram and generate running information
$ # or use the dry run scheduler, install with the extras
$ pip install -U immunopipe[diagram,runinfo,dry]
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

The docker image is build upon [`mambaorg/micromamba:2.5.0`][1]. The OS is linux/amd64. Other than the default directories, the following directories are also created or should be mapped during the run:

- `/immunopipe`: The directory where the source code of the pipeline is. It is general a clone of the [repository][2]. The pipeline is also installed from this directory.
- `/workdir`: The working directory. It is the directory where the pipeline is run. It is recommended to map the current directory (`.`) to this directory.

## Prepare to run the pipeline via Google Batch Jobs

There are two ways of running the pipeline via Google Batch Jobs: using the `gbatch` scheduler (provided by [`xqute`][3]) or using [`pipen-cli-gbatch`][4]. See more details in [Running the pipeline via Google Batch Jobs](./running.md#run-the-pipeline-using-google-cloud-batch-jobs).

In addition to prepare the docker image in the artifact registry (or docker hub if your google cloud project allows pulling from docker hub), you also need to install some dependencies locally.

If you choose to use the `gbatch` scheduler, in addition to installing the pipeline:

```shell
$ pip install -U immunopipe

# install cloud dependencies
$ pip install -U panpath[gs,async-gs]
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

[1]: https://hub.docker.com/r/mambaorg/micromamba
[2]: https://github.com/pwwang/immunopipe
[3]: https://github.com/pwwang/xqute
[4]: https://github.com/pwwang/pipen-cli-gbatch
