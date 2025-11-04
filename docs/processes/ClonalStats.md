# ClonalStats

Visualize the clonal information.

Using [`scplotter`](https://github.com/pwwang/scplotter) to visualize the clonal
information.<br />

## Input

- `screpfile`:
    The `scRepertoire` object in RDS/qs format

## Output

- `outdir`: *Default: `{{in.screpfile | stem}}.clonalstats`*. <br />
    The output directory containing the plots

## Environment Variables

- `mutaters` *(`type=json;order=-9`)*: *Default: `{}`*. <br />
    The mutaters passed to `dplyr::mutate()` to add new variables.<br />
    When the object loaded form `in.screpfile` is a list, the mutaters will be applied to each element.<br />
    The keys are the names of the new variables, and the values are the expressions.<br />
    When it is a `Seurat` object, typically an output of `scRepertoire::combineExpression()`,
    the mutaters will be applied to the `meta.data`.<br />
- `viz_type` *(`choice`)*:
    The type of visualization to generate.<br />
    - `volume`:
        The volume of the clones using [`ClonalVolumePlot`](https://pwwang.github.io/scplotter/reference/ClonalVolumePlot.html)
    - `abundance`:
        The abundance of the clones using [`ClonalAbundancePlot`](https://pwwang.github.io/scplotter/reference/ClonalAbundancePlot.html)
    - `length`:
        The length of the CDR3 sequences using [`ClonalLengthPlot`](https://pwwang.github.io/scplotter/reference/ClonalLengthPlot.html)
    - `residency`:
        The residency of the clones using [`ClonalResidencyPlot`](https://pwwang.github.io/scplotter/reference/ClonalResidencyPlot.html)
    - `dynamics`:
        The dynamics of the clones using [`ClonalDynamicsPlot`](https://pwwang.github.io/scplotter/reference/ClonalDynamicsPlot.html)
    - `composition`:
        The composition of the clones using [`ClonalCompositionPlot`](https://pwwang.github.io/scplotter/reference/ClonalCompositionPlot.html)
    - `overlap`:
        The overlap of the clones using [`ClonalOverlapPlot`](https://pwwang.github.io/scplotter/reference/ClonalOverlapPlot.html)
    - `diversity`:
        The diversity of the clones using [`ClonalDiversityPlot`](https://pwwang.github.io/scplotter/reference/ClonalDiversityPlot.html)
    - `geneusage`:
        The gene usage of the clones using [`ClonalGeneUsagePlot`](https://pwwang.github.io/scplotter/reference/ClonalGeneUsagePlot.html)
    - `positional`:
        The positional information of the clones using [`ClonalPositionalPlot`](https://pwwang.github.io/scplotter/reference/ClonalPositionalPlot.html)
    - `kmer`:
        The kmer information of the clones using [`ClonalKmerPlot`](https://pwwang.github.io/scplotter/reference/ClonalKmerPlot.html)
    - `rarefaction`:
        The rarefaction curve of the clones using [`ClonalRarefactionPlot`](https://pwwang.github.io/scplotter/reference/ClonalRarefactionPlot.html)
- `subset`:
    An expression to subset the data before plotting.<br />
    Similar to `mutaters`, it will be applied to each element by `dplyr::filter()` if the object
    loaded form `in.screpfile` is a list; otherwise, it will be applied to
    `subset(sobj, subset = <expr>)` if the object is a `Seurat` object.<br />
- `devpars` *(`ns`)*:
    The parameters for the plotting device.<br />
    - `width` *(`type=int`)*:
        The width of the device
    - `height` *(`type=int`)*:
        The height of the device
    - `res` *(`type=int`)*: *Default: `100`*. <br />
        The resolution of the device
- `more_formats` *(`list`)*: *Default: `[]`*. <br />
    The extra formats to save the plots in, other than PNG.<br />
- `save_code` *(`flag`)*: *Default: `False`*. <br />
    Whether to save the code used to generate the plots
    Note that the data directly used to generate the plots will also be saved in an `rda` file.<br />
    Be careful if the data is large as it may take a lot of disk space.<br />
- `descr`:
    The description of the plot, used to show in the report.<br />
- `<more>`:
    The arguments for the plot function
    See the documentation of the corresponding plot function for the details
- `cases` *(`type=json`)*: *Default: `{'Clonal Volume': Diot({'viz_type': 'volume'}), 'Clonal Abundance': Diot({'viz_type': 'abundance'}), 'CDR3 Length': Diot({'viz_type': 'length'}), 'Clonal Diversity': Diot({'viz_type': 'diversity'})}`*. <br />
    The cases to generate the plots if we have multiple cases.<br />
    The keys are the names of the cases, and the values are the arguments for the plot function.<br />
    The arguments in `envs` will be used if not specified in `cases`, except for `mutaters`.<br />
    Sections can be specified as the prefix of the case name, separated by `::`.<br />
    For example, if you have a case named `Clonal Volume::Case1`, the plot will be put in the
    section `Clonal Volume`. By default, when there are multiple cases for the same 'viz_type', the name of the 'viz_type' will be used
    as the default section name (for example, when 'viz_type' is 'volume', the section name will be 'Clonal Volume').<br />
    When there is only a single case, the section name will default to 'DEFAULT', which will not be shown
    in the report.<br />

