# TopExpressingGenesOfAllCells

Top expressing genes for clusters of all cells.



See also [TopExpressingGenes](./TopExpressingGenes.md).<br />

## Input

- `srtobj`:
    The seurat object in RDS format

## Output

- `outdir`: *Default: `{{in.srtobj | stem}}.top_expressing_genes`*. <br />
    The output directory for the tables and plots

## Environment Variables

- `dbs` *(`list`)*: *Default: `['KEGG_2021_Human', 'MSigDB_Hallmark_2020']`*. <br />
    The dbs to do enrichment analysis for significant
    markers See below for all libraries.<br />
    <https://maayanlab.cloud/Enrichr/#libraries>
- `n` *(`type=int`)*: *Default: `250`*. <br />
    The number of top expressing genes to find.<br />
- `subset`:
    An expression to subset the cells for each case.<br />

