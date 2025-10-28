# LoadingRNAFromSeurat

Load RNA data from a Seurat object, instead of RNAData from SampleInfo



## Input

- `infile`:
    An RDS or qs2 format file containing a Seurat object.<br />

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

