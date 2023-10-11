# TESSA

[Tessa][1] is a Bayesian model to integrate T cell receptor (TCR) sequence profiling with transcriptomes of T cells. Enabled by the recently developed single cell sequencing techniques, which provide both TCR sequences and RNA sequences of each T cell concurrently, Tessa maps the functional landscape of the TCR repertoire, and generates insights into understanding human immune response to diseases. As the first part of tessa, BriseisEncoder is employed prior to the Bayesian algorithm to capture the TCR sequence features and create numerical embeddings. We showed that the reconstructed Atchley Factor matrices and CDR3 sequences, generated through the numerical embeddings, are highly similar to their original counterparts. The CDR3 peptide sequences are constructed via a RandomForest model applied on the reconstructed Atchley Factor matrices.

For more information, please refer to the following papers:

- [Mapping the Functional Landscape of TCR Repertoire, Zhang, Z., Xiong, D., Wang, X. et al. 2021.][2]
- [Deep learning-based prediction of the T cell receptor–antigen binding specificity, Lu, T., Zhang, Z., Zhu, J. et al. 2021.][3]

/// Note
The dependencies of TESSA are not included in the docker image of immunopipe with tag without `-full` suffix. If you want to use TESSA, please use the docker image with tag with `-full` suffix, or install the dependencies manually.
///

When finished, two columns will be added to the `meta.data` of the `Seurat` object:

- `TESSA_Cluster`: The cluster assignments from TESSA.
- `TESSA_Cluster_Size`: The number of cells in each cluster.

These columns can be then used for further downstream analysis to explore the functional landscape of the TCR repertoire.

## Environment variables

- `python`: The path of python with `TESSA`'s dependencies installed
- `prefix`: The prefix to the barcodes of TCR data. You can use placeholder
    like `{Sample}_` to use the meta data from the immunarch object.
- `within_sample` (`flag`): Whether the TCR networks are constructed only
    within TCRs from the same sample/patient (True) or with all the
    TCRs in the meta data matrix (False).
- `predefined_b` (`flag`): Whether use the predefined `b` or not.
    Please check the paper of tessa for more details about the b vector.
    If True, the tessa will not update b in the MCMC iterations.
- `max_iter` (`type=int`): The maximum number of iterations for MCMC.
- `assay`: Which assay to use to extract the expression matrix.
- `save_tessa` (`flag`): Save tessa detailed results to seurat object?
    They will be saved to `sobj@misc$tessa`.

[1]: https://github.com/jcao89757/TESSA
[2]: https://www.nature.com/articles/s41592-020-01020-3
[3]: https://www.nature.com/articles/s42256-021-00383-2
