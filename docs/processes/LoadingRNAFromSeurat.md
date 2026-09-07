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

