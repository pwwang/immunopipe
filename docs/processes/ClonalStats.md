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
    - `stats`:
        The stats of the clones using [`ClonalStatsPlot`](https://pwwang.github.io/scplotter/reference/ClonalStatsPlot.html)
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

## Examples

### Clonal Volume

```toml
[ClonalStats.envs.cases."Clonal Volume"]
viz_type = "volume"
x_text_angle = 45
```

![Clonal_Volume](https://raw.githubusercontent.com/pwwang/immunopipe/tests-output/clonalstats/ClonalStats/sampleinfo.scRep.clonalstats/Number-of-Clones/Clonal-Volume.png){: width="80%"}

### Clonal Volume by Diagnosis

```toml
[ClonalStats.envs.cases."Clonal Volume by Diagnosis"]
viz_type = "volume"
x = "seurat_clusters"
group_by = "Diagnosis"
comparisons = true
```

![Clonal_Volume_by_Diagnosis](https://raw.githubusercontent.com/pwwang/immunopipe/tests-output/clonalstats/ClonalStats/sampleinfo.scRep.clonalstats/Number-of-Clones/Clonal-Volume-by-Diagnosis.png){: width="80%"}

### Clonal Abundance

```toml
[ClonalStats.envs.cases."Clonal Abundance"]
viz_type = "abundance"
```

![Clonal_Abundance](https://raw.githubusercontent.com/pwwang/immunopipe/tests-output/clonalstats/ClonalStats/sampleinfo.scRep.clonalstats/Clonal-Abundance/Clonal-Abundance.png){: width="80%"}

### Clonal Abundance Density

```toml
[ClonalStats.envs.cases."Clonal Abundance Density"]
viz_type = "abundance"
plot_type = "density"
```

![Clonal_Abundance_Density](https://raw.githubusercontent.com/pwwang/immunopipe/tests-output/clonalstats/ClonalStats/sampleinfo.scRep.clonalstats/Clonal-Abundance/Clonal-Abundance-Density.png){: width="80%"}

### CDR3 Length

```toml
[ClonalStats.envs.cases."CDR3 Length"]
viz_type = "length"
```

![CDR3_Length](https://raw.githubusercontent.com/pwwang/immunopipe/tests-output/clonalstats/ClonalStats/sampleinfo.scRep.clonalstats/Clonal-Sequence-Length/CDR3-Length.png){: width="80%"}

### CDR3 Length (Beta Chain)

```toml
[ClonalStats.envs.cases."CDR3 Length (Beta Chain)"]
viz_type = "length"
chain = "TRB"
```

![CDR3_Length_Beta_Chain](https://raw.githubusercontent.com/pwwang/immunopipe/tests-output/clonalstats/ClonalStats/sampleinfo.scRep.clonalstats/Clonal-Sequence-Length/CDR3-Length-Beta-Chain-.png){: width="80%"}

### Clonal Residency

```toml
[ClonalStats.envs.cases."Clonal Residency"]
viz_type = "residency"
group_by = "Diagnosis"
chain = "TRB"
clone_call = "gene"
groups = ["Colitis", "NoColitis"]
```

![Clonal_Residency](https://raw.githubusercontent.com/pwwang/immunopipe/tests-output/clonalstats/ClonalStats/sampleinfo.scRep.clonalstats/Clonal-Residency/Clonal-Residency.png){: width="80%"}

### Clonal Residency (UpSet Plot)

```toml
[ClonalStats.envs.cases."Clonal Residency (UpSet Plot)"]
viz_type = "residency"
plot_type = "upset"
group_by = "Diagnosis"
chain = "TRB"
clone_call = "gene"
groups = ["Colitis", "NoColitis"]
devpars = {width = 800}
```

![Clonal_Residency_UpSet_Plot](https://raw.githubusercontent.com/pwwang/immunopipe/tests-output/clonalstats/ClonalStats/sampleinfo.scRep.clonalstats/Clonal-Residency/Clonal-Residency-UpSet-Plot-.png){: width="80%"}

### Clonal Statistics with Expanded Clones

```toml
[ClonalStats.envs.cases."Clonal Statistics with Expanded Clones"]
viz_type = "stat"
plot_type = "pies"
group_by = "Diagnosis"
groups = ["Colitis", "NoColitis"]
clones = {"Expanded Clones In Colitis" = "sel(Colitis > 2)", "Expanded Clones In NoColitis" = "sel(NoColitis > 2)"}
subgroup_by = "seurat_clusters"
pie_size = "sqrt"
show_row_names = true
show_column_names = true
devpars = {width = 720}
```

![Clonal_Statistics_with_Expanded_Clones](https://raw.githubusercontent.com/pwwang/immunopipe/tests-output/clonalstats/ClonalStats/sampleinfo.scRep.clonalstats/Clonal-Statistics/Clonal-Statistics-with-Expanded-Clones.png){: width="80%"}

### Hyperexpanded Clonal Dynamics

```toml
[ClonalStats.envs.cases."Hyperexpanded Clonal Dynamics"]
viz_type = "stat"
plot_type = "sankey"
group_by = "Diagnosis"
chain = "TRB"
groups = ["Colitis", "NoColitis"]
clones = {"Hyper-Expanded Clones In Colitis" = "sel(Colitis > 5)", "Hyper-Expanded Clones In NoColitis" = "sel(NoColitis > 5)"}
devpars = {width = 800}
```

![Hyperexpanded_Clonal_Dynamics](https://raw.githubusercontent.com/pwwang/immunopipe/tests-output/clonalstats/ClonalStats/sampleinfo.scRep.clonalstats/Clonal-Statistics/Hyperexpanded-Clonal-Dynamics.png){: width="80%"}

### Clonal Composition

```toml
[ClonalStats.envs.cases."Clonal Composition"]
viz_type = "composition"
x_text_angle = 45
```

![Clonal_Composition](https://raw.githubusercontent.com/pwwang/immunopipe/tests-output/clonalstats/ClonalStats/sampleinfo.scRep.clonalstats/Clonal-Composition.png){: width="80%"}

### Clonal Overlapping

```toml
viz_type = "overlap"
chain = "TRB"
clone_call = "gene"
```

![Clonal_Overlapping](https://raw.githubusercontent.com/pwwang/immunopipe/tests-output/clonalstats/ClonalStats/sampleinfo.scRep.clonalstats/Clonal-Overlapping.png){: width="80%"}

### Clonal Diversity

```toml
[ClonalStats.envs.cases."Clonal Diversity"]
# method = "shannon"  # default
viz_type = "diversity"
x_text_angle = 45
```

![Clonal_Diversity](https://raw.githubusercontent.com/pwwang/immunopipe/tests-output/clonalstats/ClonalStats/sampleinfo.scRep.clonalstats/Clonal-Diversity/Clonal-Diversity.png){: width="80%"}

### Clonal Diversity (gini.coeff, by Diagnosis)

```toml
[ClonalStats.envs.cases."Clonal Diversity (gini.coeff, by Diagnosis)"]
method = "gini.coeff"
viz_type = "diversity"
plot_type = "box"
group_by = "Diagnosis"
comparisons = true
devpars = {height = 600, width = 600}
```

![Clonal_Diversity_gini_coeff_by_Diagnosis](https://raw.githubusercontent.com/pwwang/immunopipe/tests-output/clonalstats/ClonalStats/sampleinfo.scRep.clonalstats/Clonal-Diversity/Clonal-Diversity-gini-coeff-by-Diagnosis-.png){: width="80%"}

### Gene Usage Frequency

```toml
[ClonalStats.envs.cases."Gene Usage Frequency"]
viz_type = "geneusage"
devpars = {width = 1200}
```

![Gene_Usage_Frequency](https://raw.githubusercontent.com/pwwang/immunopipe/tests-output/clonalstats/ClonalStats/sampleinfo.scRep.clonalstats/Gene-Usage-Frequency.png){: width="80%"}

### Positional amino acid frequency

```toml
[ClonalStats.envs.cases."Positional amino acid frequency"]
viz_type = "positional"
# method = "AA"  # default
devpars = {width = 1600}
```

![Positional_amino_acid_frequency](https://raw.githubusercontent.com/pwwang/immunopipe/tests-output/clonalstats/ClonalStats/sampleinfo.scRep.clonalstats/Positional-Properties/Positional-amino-acid-frequency.png){: width="80%"}

### Positional shannon entropy

```toml
[ClonalStats.envs.cases."Positional shannon entropy"]
viz_type = "positional"
method = "shannon"
devpars = {width = 1200}
```

![Positional_shannon_entropy](https://raw.githubusercontent.com/pwwang/immunopipe/tests-output/clonalstats/ClonalStats/sampleinfo.scRep.clonalstats/Positional-Properties/Positional-shannon-entropy.png){: width="80%"}

### 3-Mer Frequency

```toml
[ClonalStats.envs.cases."3-Mer Frequency"]
viz_type = "kmer"
k = 3  # default is 3
devpars = {width = 800}
```

![3_Mer_Frequency](https://raw.githubusercontent.com/pwwang/immunopipe/tests-output/clonalstats/ClonalStats/sampleinfo.scRep.clonalstats/3-Mer-Frequency.png){: width="80%"}

### Rarefaction Curve

```toml
[ClonalStats.envs.cases."Rarefaction Curve"]
viz_type = "rarefaction"
```

![Rarefaction_Curve](https://raw.githubusercontent.com/pwwang/immunopipe/tests-output/clonalstats/ClonalStats/sampleinfo.scRep.clonalstats/Rarefaction-Curve.png){: width="80%"}

