# TopExpressingGenes

Top expressing genes for clusters of all or selected T cells.



This process finds the top expressing genes of clusters of T cells, and also
performs the enrichment analysis against the genes.<br />

The enrichment analysis is done by
[`enrichr`](https://maayanlab.cloud/Enrichr/).<br />

/// Note
There are other environment variables also available. However, they should not
be used in this process. Other environment variables are used for more
complicated cases for investigating top genes
(See [`biopipen.ns.scrna.TopExpressingGenes`](https://pwwang.github.io/biopipen/api/biopipen.ns.scrna/#biopipen.ns.scrna.TopExpressingGenes) for more details).<br />

If you are using `pipen-board` to run the pipeline
(see [here](../running.md#run-the-pipeline-via-pipen-board) and
[here](../running.md#run-the-pipeline-via-pipen-board-using-docker-image)),
you may see the other environment variables of this process are hidden and
readonly.<br />
///

## Environment Variables

- `dbs` *(`list`)*: *Default: `['KEGG_2021_Human', 'MSigDB_Hallmark_2020']`*. <br />
    The dbs to do enrichment analysis for significant
    markers See below for all libraries.<br />
    <https://maayanlab.cloud/Enrichr/#libraries>
- `n` *(`type=int`)*: *Default: `250`*. <br />
    The number of top expressing genes to find.<br />
- `subset`:
    An expression to subset the cells for each case.<br />

