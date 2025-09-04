# Getting started

You can find the nessary files and source code of this tutorial in the [example repository](https://github.com/pwwang/immunopipe-example).

In this tutorial we will show you how to run the immunopipe pipeline on a small dataset of 6 patients from 3 groups: colitis (n=2), non-colitis(n=2) and control(n=2). The dataset is part of the data used in the publication below:

- [Luoma, Adrienne M., et al. "Molecular pathways of colon inflammation induced by cancer immunotherapy." Cell 182.3 (2020): 655-671.][1]

We are using a small subset of the data to make the tutorial run faster. The full dataset can be downloaded from Gene Expression Omnibus (GEO) [GSE144469](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE144469).

## Download and prepare the data

The data can be downloaded and prepared by running the following commands:

```bash
# Clone the example repository
git clone https://github.com/pwwang/immunopipe-example.git

# Enter the example directory
cd immunopipe-example

# Download and prepare the data
bash prepare-data.sh
# The data from GSE144469 (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE144469)
# will be downloaded and extracted into:
#
#   ./prepared-data/C1
#   ./prepared-data/C2
#   ...
#
```

You may also check other files in the `data/` directory, especially the `samples.txt` file, which contains the sample information for the dataset we prepared above.

## Prepare the configuration file

To run the pipeline, we need to prepare a configuration file (recommended) or pass the arguments directly via command line. Here we will use the configuration file. See also [Configurations](./configurations.md) for more details.

As explained in the [Configurations](./configurations.md) page, we can provide a configuration file with [a minimal set of configuration items](./configurations.md#minimal-configurations) to get the pipeline running. The only required configuration item is the input file for the [`SampleInfo`](./processes/SampleInfo.md) process. However, here we want to give the pipeline a different name and output directory to distinguish it from other runs with a different set of configurations.

The configuration file shall be in the [TOML](https://toml.io/en/) format. We can create a file named `ImmunopipeMinimal.config.toml` with the following content:

```toml
name = "ImmunopipeMinimal"
outdir = "minimal"

[SampleInfo.in]
infile = [ "data/samples.txt" ]
```

## Run the pipeline

The easiest way to run the pipeline is to run it within the docker container. We can use the following command to run the pipeline with the configuration file we just created:

/// tab | Using docker

```bash
docker run \
    --rm -w /workdir -v .:/workdir \
    justold/immunopipe:master \
    @ImmunopipeMinimal.config.toml
```

///

/// tab | Using singularity

```bash
singularity run \
    --pwd /workdir -B .:/workdir,/tmp -c -e --writable-tmpfs \
    docker://justold/immunopipe:master \
    @ImmunopipeMinimal.config.toml
```

///

/// tab | Using apptainer

```bash
apptainer run \
    --pwd /workdir -B .:/workdir,/tmp -c -e --unsquash --writable-tmpfs \
    docker://justold/immunopipe:master \
    @ImmunopipeMinimal.config.toml
```

///

/// Tip
`docker`, `singularity` and `apptainer` commands map the current directory (`.`) to the `/workdir` directory in the container. To get the detailed directory structure in the container, please refer to the [The directory structure in the container](https://pwwang.github.io/immunopipe/installation/#the-directory-structure-in-the-container).
///

/// Tip
If you want to install and run the pipeline without docker, please refer to the [Installation](https://pwwang.github.io/immunopipe/installation/) and [Running the pipeline](https://pwwang.github.io/immunopipe/running/) pages for more details.
///

/// Note
You need at least 16GB of memory to run the pipeline with the example dataset and minimal configuration.

You may also need to decrease `ncores` of some processes to avoid running out of memory. For example:

```diff
[SeuratClusteringOfAllCells.envs]
- ncores = 8
+ ncores = 4
```

///

## Check the results

With that "minimal" configuration file, only a subset of the processes will be run. See also [Enabling/Disabling processes](./configurations.md#enablingdisabling-processes). The results will be saved in the `minimal` directory. You can also check the reports at `minimal/REPORTS/index.html` with a web browser.

You can also visit the following link to see the reports of the pipeline we just ran:

<http://imp.pwwang.com/minimal/REPORTS/index.html>

## Next steps

You may read through this documentation to learn more about the pipeline and how to configure it. There is also a configuration file, named [`Immunopipe.config.toml`][2] in the example repository, with more processes enabled. You can use it to run the pipeline with the dataset prepared above. Check out the following link for the reports:

<http://imp.pwwang.com/output/REPORTS/index.html>

/// nOTE
The results provided by this example configuration files are for demonstration purpose only. They are not intended to be used for any scientific analysis.
///

You may also want to try other routes of the pipeline with the prepared data. These routes are defined in:

- `ImmunopipeMinimalNoTCR.config.toml`: The configuration for minimal analyses without scTCR-/scBCR-seq data.
- `ImmunopipeMinimalSupervised.config.toml`: The configuration for minimal analyses with supervised clustering of T cells.
- `ImmunopipeNoTCR.config.toml`: The configuration for full analyses without scTCR-/scBCR-seq data.
- `ImmunopipeWSNoTCR.config.toml`: The configuration for full analyses without scTCR-/scBCR-seq data, but with selection of T cells.
- `ImmunopipeSupervised.config.toml`: The configuration for full analyses with supervised clustering of T cells.

Also check out the [gallery](./gallery.md) for more real-world examples.

[1]: https://www.cell.com/cell/fulltext/S0092-8674(20)30688-7
[2]: https://github.com/pwwang/immunopipe-example/blob/master/Immunopipe.config.toml
