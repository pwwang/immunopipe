# TopExpressingGenes

This process finds the top expressing genes of clusters of all cells, and also performs the enrichment analysis against the genes.

The enrichment analysis is done by [`enrichr`][1].

## Environment variables

- `group-by`: The column name in metadata with the cluster information. Default: `seurat_clusters`.
- `dbs` (`list`): The dbs to do enrichment analysis for significant
    markers See below for all libraries.
    <https://maayanlab.cloud/Enrichr/#libraries>.
    Default: `["GO_Biological_Process_2021", "GO_Cellular_Component_2021", "GO_Molecular_Function_2021", "KEGG_2021_Human"]`
- `n` (`type=int`): The number of genes to find for each cluster. Default: `250`.


/// Note
There are other environment variables also available. However, they should not be used in this process. Other environment variables are used for more complicated cases for investigating top genes (See [`biopipen.ns.scrna.TopExpressingGenes`][2] for more details).

If you are using `pipen-board` to run the pipeline (see [here](../running.md#run-the-pipeline-via-pipen-board) and [here](../running.md#run-the-pipeline-via-pipen-board-using-docker-image)), you may see the other environment variables of this process are hidden and readonly.
///

[1]: https://maayanlab.cloud/Enrichr/
[2]: https://pwwang.github.io/biopipen/api/biopipen.ns.scrna/#biopipen.ns.scrna.TopExpressingGenes
