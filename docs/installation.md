# Installation

## Install the pipeline itself

```shell
$ pip install -U immunopipe
```

## Install the dependencies

`immunopipe` is built upon [`pipen`](https://github.com/pwwang/pipen) framework, and a number of packages in `R` and `python`. It's not recommended to install the packages manually. Instead, you can use the provided `environment.yml` to create a conda environment.

```shell
$ conda env create \
    -n immunopipe \
    -f https://raw.githubusercontent.com/pwwang/immunopipe/dev/docker/environment.yml
```

If the URL doesn't work, you can download the file and create the environment locally.

For more detailed instructions of `conda env create`, please refer to [conda docs](https://docs.conda.io/projects/conda/en/latest/commands/env/create.html).

!!! note
    The pipeline itself is not included in the conda environment. You need to install it separately.

    ```shell
    $ conda activate immunopipe
    $ pip install -U immunopipe
    ```

## Use the docker image

You can also use the docker image to run the pipeline. The image is built upon `miniconda3` and `micromamba` is used as the package manager. The image is available at [Docker Hub](https://hub.docker.com/r/justold/immunopipe).

To pull the image:

```shell
$ docker pull justold/immunopipe:<tag>
```

If you are using `singularity`, you can pull and convert the image to `sif` format:

```shell
$ singularity pull docker://justold/immunopipe:<tag>
```

To run the pipeline use the image, please refer to [Running the pipeline](./running.md).
