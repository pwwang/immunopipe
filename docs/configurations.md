
## Arguments from command line

`immunopipe` is written using `pipen`. The pipeline is composed of a number processes. Each process' input (start processes), envs and other configurations can be passed from command line.

For example, to run V-J usage analysis on 8 cores:
```shell
immunopipe --VJUsage.envs.ncores 8 ...
```

## Arguments from configuration file

As the argument parsing for pipen is supported by `pipen-args`, which allows us to pass arguments from both command line and a configuration file.

To specify the configuration file:
```shell
immunopipe --config config.toml ...
```

The configuration file is in [TOML](https://github.com/toml-lang/toml) format.

You can specify the same argument from both command line and the configuration file, but the one from the command line has higher priority.

To specify the arguments, for example, `--VJUsage.envs.cores=8` from the configuration file:

```toml
[VJUsage.envs]
cores = 8
```

You can also specify some arguments for the entire pipeline, for example, the scheduler of the pipeline:

```toml
scheduler = "sge"

[VJUsage.envs]
cores = 8
```

## Exploration of all configurable items for the pipeline or for a specific process

To see all available configuration items for the pipeline, run:

```shell
immunopipe --full
```

For a process, for example, for `VJUsage`:

```shell
pipen run tcr VJUsage --full
```

Or you can check the docstring of the process in `biopipen`'s source code.

## Configurations for other modules

### `MARKERS_FINDER`

- See: [https://pwwang.github.io/immunopipe/markers-finder/](https://pwwang.github.io/immunopipe/markers-finder/)

### `GENE_EXPR_INVESTIGATION_CLUSTERS`

- See: [https://pwwang.github.io/immunopipe/gene-expr-investigation-for-each-cluster/](https://pwwang.github.io/immunopipe/gene-expr-investigation-for-each-cluster/)

### `DIM_PLOTS`

- See: [https://pwwang.github.io/immunopipe/dimplots/](https://pwwang.github.io/immunopipe/dimplots/)]

### `RADAR_PLOTS`

- See: [https://pwwang.github.io/immunopipe/radar/](https://pwwang.github.io/immunopipe/radar/)]
