# Getting started

In this tutorial we will show you how to run the immunopipe pipeline on a small dataset of 3 patients with paired tumor and normal samples. The dataset is part of the data used in the publications below:

- [Wu TD, Madireddi S, de Almeida PE, Banchereau R et al. Peripheral T cell expansion predicts tumour infiltration and clinical response. Nature 2020 Mar;579(7798):274-278.][1]
- [Banta KL, Xu X, Chitre AS, Au-Yeung A et al. Mechanistic convergence of the TIGIT and PD-1 inhibitory pathways necessitates co-blockade to optimize anti-tumor CD8+ T cell responses. Immunity 2022 Mar 8;55(3):512-526.e9.][2]

We are using a small subset of the data to make the tutorial run faster. The full dataset can be downloaded from Gene Expression Omnibus (GEO) [GSE139555](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE139555).

## Download and prepare the data

The data can be downloaded and prepared by running the following commands:

```bash
# Clone the example repository
git clone https://github.com/pwwang/immunopipe-example.git

# Enter the example directory
cd immunopipe-example

# Download and prepare the data
bash prepare-data.sh
# The data from GSE139555 (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE139555)
# will be downloaded and extracted into:
#
#   ./prepared-data/LT1
#   ./prepared-data/LN2
#   ...
#
# Only first 6 samples (LN, LT) are used.
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

The easiest way to run the pipeline is to run it with docker. We can use the following command to run the pipeline with the configuration file we just created:

```bash
docker run \
    --rm -w /workdir -v .:/workdir \
    justold/immunopipe:0.6.0 \
    @ImmunopipeMinimal.config.toml
```

or with singularity:

```bash
singularity run \
    --pwd /workdir -B .:/workdir -c -e -w \
    docker://justold/immunopipe:0.6.0 \
    @ImmunopipeMinimal.config.toml
```

!!! tip

    Both the `docker` and `singularity` commands map the current directory (`.`) to the `/workdir` directory in the container. To get the detailed directory structure in the container, please refer to the [The directory structure in the container](./installation.md#the-directory-structure-in-the-container).

!!! tip

    If you want to install and run the pipeline without docker, please refer to the [Installation](./installation.md) and [Running the pipeline](./running.md) pages for more details.

!!! note

    You need at least 8GB of memory to run the pipeline with the example dataset and minimal configuration. 16GB or more is recommended.

## Check the results

With that "minimal" configuration file, only a subset of the processes will be run. See also [Enabling/Disabling processes](./configurations.md#enablingdisabling-processes). The results will be saved in the `minimal` directory. You can also check the reports at `minimal/REPORTS/index.html` with a web browser.

You can also visit the following link to see the reports of the pipeline we just ran:

<http://imp.pwwang.com/minimal/REPORTS/index.html>

## Next steps

You may read through this documentation to learn more about the pipeline and how to configure it. There is also a configuration file, named [`Immunopipe.config.toml`][3] in the example repository, with more processes enabled. You can use it to run the pipeline with the dataset prepared above. Check out the following link for the reports:

<http://imp.pwwang.com/output/REPORTS/index.html>


[1]: https://www.ncbi.nlm.nih.gov/pubmed/32103181
[2]: https://www.ncbi.nlm.nih.gov/pubmed/35263569
[3]: https://github.com/pwwang/immunopipe-example/blob/master/Immunopipe.config.toml
