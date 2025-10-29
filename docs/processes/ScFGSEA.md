# ScFGSEA

Gene set enrichment analysis for cells in different groups using `fgsea`

This process allows us to do Gene Set Enrichment Analysis (GSEA) on the expression data,
but based on variaties of grouping, including the from the meta data and the
scTCR-seq data as well.<br />

The GSEA is done using the
[fgsea](https://bioconductor.org/packages/release/bioc/html/fgsea.html) package,
which allows to quickly and accurately calculate arbitrarily low GSEA P-values
for a collection of gene sets.<br />
The fgsea package is based on the fast algorithm for preranked GSEA described in
[Subramanian et al. 2005](https://www.pnas.org/content/102/43/15545).<br />

For each case, the process will generate a table with the enrichment scores for
each gene set, and GSEA plots for the top gene sets.<br />

## Input

- `srtobj`:
    The seurat object in RDS format

## Output

- `outdir`: *Default: `{{(in.casefile or in.srtobj) | stem0}}.fgsea`*. <br />
    The output directory for the results and plots

## Environment Variables

- `ncores` *(`type=int`)*: *Default: `1`*. <br />
    Number of cores for parallelization
    Passed to `nproc` of `fgseaMultilevel()`.<br />
- `mutaters` *(`type=json`)*: *Default: `{}`*. <br />
    The mutaters to mutate the metadata.<br />
    The key-value pairs will be passed the `dplyr::mutate()` to mutate the metadata.<br />
    You can also use the clone selectors to select the TCR clones/clusters.<br />
    See <https://user.github.io/scplotter/reference/clone_selectors.html>.<br />

- `group_by`:
    The column name in metadata to group the cells.<br />
- `ident_1`:
    The first group of cells to compare
- `ident_2`:
    The second group of cells to compare, if not provided, the rest of the cells that are not `NA`s in `group_by` column are used for `ident-2`.<br />
- `each`:
    The column name in metadata to separate the cells into different subsets to do the analysis.<br />
- `subset`:
    An expression to subset the cells.<br />
- `gmtfile`: *Default: `KEGG_2021_Human`*. <br />
    The pathways in GMT format, with the gene names/ids in the same format as the seurat object.<br />
    One could also use a URL to a GMT file. For example, from <https://download.baderlab.org/EM_Genesets/current_release/Human/symbol/Pathways/>.<br />
- `method` *(`choice`)*: *Default: `s2n`*. <br />
    The method to do the preranking.<br />
    - `signal_to_noise`:
        Signal to noise.<br />
        The larger the differences of the means (scaled by the standard deviations);
        that is, the more distinct the gene expression is in each phenotype and the more the gene
        acts as a "class marker".<br />
    - `s2n`:
        Alias of signal_to_noise.<br />
    - `abs_signal_to_noise`:
        The absolute value of signal_to_noise.<br />
    - `abs_s2n`:
        Alias of abs_signal_to_noise.<br />
    - `t_test`:
        T test.<br />
        Uses the difference of means scaled by the standard deviation and number of samples.<br />
    - `ratio_of_classes`:
        Also referred to as fold change.<br />
        Uses the ratio of class means to calculate fold change for natural scale data.<br />
    - `diff_of_classes`:
        Difference of class means.<br />
        Uses the difference of class means to calculate fold change for nature scale data
    - `log2_ratio_of_classes`:
        Log2 ratio of class means.<br />
        Uses the log2 ratio of class means to calculate fold change for natural scale data.<br />
        This is the recommended statistic for calculating fold change for log scale data.<br />
- `top` *(`type=auto`)*: *Default: `20`*. <br />
    Do gsea table and enrich plot for top N pathways.<br />
    If it is < 1, will apply it to `padj`, selecting pathways with `padj` < `top`.<br />
- `eps` *(`type=float`)*: *Default: `0`*. <br />
    This parameter sets the boundary for calculating the p value.<br />
    See <https://rdrr.io/bioc/fgsea/man/fgseaMultilevel.html>
- `alleach_plots_defaults` *(`ns`)*:
    Default options for the plots to generate for all pathways.<br />
    - `plot_type`: *Default: `heatmap`*. <br />
        The type of the plot, currently either dot or heatmap (default)
    - `devpars` *(`ns`)*:
        The device parameters for the plots.<br />
        - `res` *(`type=int`)*: *Default: `100`*. <br />
            The resolution of the plots.<br />
        - `height` *(`type=int`)*:
            The height of the plots.<br />
        - `width` *(`type=int`)*:
            The width of the plots.<br />
    - `<more>`:
        See <https://user.github.io/biopipen.utils.R/reference/VizGSEA.html>.<br />
- `alleach_plots` *(`type=json`)*: *Default: `{}`*. <br />
    Cases of the plots to generate for all pathways.<br />
    The keys are the names of the cases and the values are the dicts inherited from `alleach_plots_defaults`.<br />
- `minsize` *(`type=int`)*: *Default: `10`*. <br />
    Minimal size of a gene set to test. All pathways below the threshold are excluded.<br />
- `maxsize` *(`type=int`)*: *Default: `100`*. <br />
    Maximal size of a gene set to test. All pathways above the threshold are excluded.<br />
- `rest` *(`type=json;order=98`)*: *Default: `{}`*. <br />
    Rest arguments for [`fgsea()`](https://rdrr.io/bioc/fgsea/man/fgsea.html)
    See also <https://rdrr.io/bioc/fgsea/man/fgseaMultilevel.html>
- `cases` *(`type=json;order=99`)*: *Default: `{}`*. <br />
    If you have multiple cases, you can specify them here.<br />
    The keys are the names of the cases and the values are the above options except `mutaters`.<br />
    If some options are not specified, the default values specified above will be used.<br />
    If no cases are specified, the default case will be added with the name `GSEA`.<br />

