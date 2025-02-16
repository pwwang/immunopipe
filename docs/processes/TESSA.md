# TESSA

Tessa is a Bayesian model to integrate T cell receptor (TCR) sequence profiling with transcriptomes of T cells.

Enabled by the recently developed single cell sequencing techniques, which provide
both TCR sequences and RNA sequences of each T cell concurrently, Tessa maps the
functional landscape of the TCR repertoire, and generates insights into
understanding human immune response to diseases. As the first part of tessa,
BriseisEncoder is employed prior to the Bayesian algorithm to capture the TCR
sequence features and create numerical embeddings. We showed that the reconstructed
Atchley Factor matrices and CDR3 sequences, generated through the numerical
embeddings, are highly similar to their original counterparts. The CDR3 peptide
sequences are constructed via a RandomForest model applied on the reconstructed
Atchley Factor matrices.<br />

See <https://github.com/jcao89757/TESSA>

When finished, two columns will be added to the `meta.data` of the `Seurat` object:<br />

- `TESSA_Cluster`: The cluster assignments from TESSA.<br />
- `TESSA_Cluster_Size`: The number of cells in each cluster.<br />

These columns can be then used for further downstream analysis to explore the
functional landscape of the TCR repertoire.<br />

/// Note
The dependencies of TESSA are not included in the docker image of immunopipe
with tag without `-full` suffix. If you want to use TESSA, please use the
docker image with tag with `-full` suffix, or install the dependencies manually.<br />
///

## Input

- `immdata`:
    The immunarch object in RDS file or text file of TCR data loaded by
    [`ImmunarchLoading`](!!#biopipennstcrimmunarchloading)
- `srtobj`:
    The `Seurat` object, saved in RDS format, with dimension
    reduction performed if you want to use them to represent the
    transcriptome of T cells.<br />
    This could also be a tab delimited file (can be gzipped) with
    expression matrix or dimension reduction results.<br />

## Output

- `outfile`: *Default: `
        {%- if in.srtobj.lower().endswith(".rds") -%}
        {{in.srtobj | stem}}.tessa.RDS
        {%- else -%}
        {{in.immdata | stem}}.tessa.txt
        {%- endif -%}`*. <br />
    The tab-delimited file with three columns
    (`barcode`, `TESSA_Cluster` and `TESSA_Cluster_Size`) or
    an RDS file if  `in.srtobj` is an RDS file of a Seurat object, with
    `TESSA_Cluster` and `TESSA_Cluster_Size` added to the `meta.data`

## Environment Variables

- `python`: *Default: `python`*. <br />
    The path of python with `TESSA`'s dependencies installed
- `prefix`:
    The prefix of the cell barcodes in the `Seurat` object.<br />
    Once could use a fixed prefix, or a placeholder with the column
    name in meta data. For example, `"{Sample}_"` will replace the
    placeholder with the value of the column `Sample` in meta data.<br />
    If `in.immdata` is text file, the prefix will be ignored and the
    barcode should be already prefixed.<br />
    If `None` and `in.immdata` is RDS file, `immdata$prefix` will be used.<br />
- `within_sample` *(`flag`)*: *Default: `False`*. <br />
    Whether the TCR networks are constructed only
    within TCRs from the same sample/patient (True) or with all the
    TCRs in the meta data matrix (False).<br />
- `assay`:
    Which assay to use to extract the expression matrix.<br />
    Only works if `in.srtobj` is an RDS file of a Seurat object.<br />
    By default, if `SCTransform` is performed, `SCT` will be used.<br />
- `predefined_b` *(`flag`)*: *Default: `False`*. <br />
    Whether use the predefined `b` or not.<br />
    Please check the paper of tessa for more details about the b vector.<br />
    If True, the tessa will not update b in the MCMC iterations.<br />
- `max_iter` *(`type=int`)*: *Default: `1000`*. <br />
    The maximum number of iterations for MCMC.<br />
- `save_tessa` *(`flag`)*: *Default: `False`*. <br />
    Save tessa detailed results to seurat object?<br />
    Only works if `in.srtobj` is an RDS file of a Seurat object.<br />
    It will be saved to `sobj@misc$tessa`.<br />

## Reference

- 'Mapping the Functional Landscape of TCR Repertoire.',
    Zhang, Z., Xiong, D., Wang, X. et al. 2021.<br />
    [link](https://www.nature.com/articles/s41592-020-01020-3)
- 'Deep learning-based prediction of the T cell receptor-antigen
    binding specificity.', Lu, T., Zhang, Z., Zhu, J. et al. 2021.<br />
    [link](https://www.nature.com/articles/s42256-021-00383-2)

## Metadata

The metadata of the `Seurat` object will be updated with the TESSA clusters
and the cluster sizes:<br />

![TESSA-metadata](../..//processes/images/TESSA-metadata.png)

