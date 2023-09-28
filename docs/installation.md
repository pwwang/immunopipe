# Installation

## Install the pipline and the dependencies using conda

/// Tip
If you plan to use the docker image, you can skip this section.
///

`immunopipe` is built upon [`pipen`](https://github.com/pwwang/pipen) framework, and a number of packages in `R` and `python`. It's not recommended to install the packages manually. Instead, you can use the provided `environment.yml` to create a conda environment.

```shell
$ conda env create \
    -n immunopipe \
    -f https://raw.githubusercontent.com/pwwang/immunopipe/dev/docker/environment.yml
```

If the URL doesn't work, you can download the file and create the environment locally.

For more detailed instructions of `conda env create`, please refer to [conda docs](https://docs.conda.io/projects/conda/en/latest/commands/env/create.html).

/// Reminder
The pipeline itself is NOT included in the conda environment. You need to install it separately.

```shell
$ conda activate immunopipe
$ pip install -U immunopipe
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

/// Note
Note the GPU support is not available for TCR clustering using `faiss-gpu`. If you want to use the GPU support, please install the pipeline and the dependencies using conda, and then install the packages with GPU support manually.
///

### The directory structure in the container

The docker image is build upon [`mambaorg/micromamba:1.4.3`][1]. The OS is linux/amd64. Other than the default directories, the following directories are also created or should be mapped during the run:

- `/immunopipe`: The directory where the source code of the pipeline is. It is general a clone of the [repository][2]. The pipeline is also installed from this directory.
- `/workdir`: The working directory. It is the directory where the pipeline is run. It is recommended to map the current directory (`.`) to this directory.

[1]: https://hub.docker.com/layers/mambaorg/micromamba/1.4.3/images/sha256-0251b94151c021c85d3e4f4ffe1fc81c436f18e01337d3b367d0f7c76ee716ac?context=explore
[2]: https://github.com/pwwang/immunopipe
