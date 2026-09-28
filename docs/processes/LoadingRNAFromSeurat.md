# LoadingRNAFromSeurat

Load RNA data from a Seurat object, instead of RNAData from SampleInfo



## Input

- `infile`:
    An [RDS](https://rdrr.io/r/base/readRDS.html) or [qs/qs2](https://github.com/qsbase/qs2)
    format file containing a Seurat object.<br />

## Output

- `outfile`: *Default: `{{in.infile | basename}}`*. <br />

## Environment Variables

- `prepared` *(`flag`)*: *Default: `False`*. <br />
    Whether the Seurat object is well-prepared for the
    pipeline (so that SeuratPreparing process is not needed).<br />
- `clustered` *(`flag`)*: *Default: `False`*. <br />
    Whether the Seurat object is clustered, so that
    `SeuratClustering` (`SeuratClusteringOfAllCells`) process or
    `SeuratMap2Ref` is not needed.<br />
    Force `prepared` to be `True` if this is `True`.<br />
- `sample`: *Default: `Sample`*. <br />
    The column name in the metadata of the Seurat object that
    indicates the sample name.<br />
    Multiple columns will be concatenated with `_` to form the sample name.<br />
- `mutaters` *(`type=json`)*: *Default: `{}`*. <br />
    The mutaters to mutate the metadata
    Keys are the names of the mutaters and values are the R expressions
    passed by `dplyr::mutate()` to mutate the metadata.<br />
- `subset`:
    An expression to subset the cells, will be passed to `dplyr::filter()`.<br />
    This will be applied after mutating the metadata.<br />
- `ncores` *(`type=int`)*: *Default: `1`*. <br />
    The number of threads used to load/save the Seurat object.<br />

## SeeAlso

- [Preparing the input](../preparing-input.md#single-cell-rna-seq-scrna-seq-data).<br />
- [Routes of the pipeline](../introduction.md#routes-of-the-pipeline).<br />

## Description

Loads the RNA data from a pre-existing `Seurat` object (an RDS or
qs/qs2 file) instead of the `RNAData` directories listed by
`SampleInfo`. This is not a wrapper of an upstream analysis tool: the
process is immunopipe's own and runs an R script that ships with
immunopipe (`immunopipe/scripts/LoadingRNAFromSeurat.R`), which reads
and writes the object with `tidyseurat` and `qs2`.<br />

## Base class

`biopipen.core.proc.Proc` - biopipen's bare process class, which
declares no parameters of its own.<br />

## Deviations

All six options of this process are added by immunopipe, since the base
class declares none: `prepared` (default `False`), `clustered`
(default `False`), `sample` (default `Sample`), `mutaters`
(default `{}`), `subset` (default `None`) and `ncores` (default `1`).<br />

