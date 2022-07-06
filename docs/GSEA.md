This module performs a "pseudo-bulk" GSEA using the [FGSEA][1] package.

An example configuration:

```toml
[ScFGSEA.envs]
# # Do GSEA on two groups of cells
# # Like seudo-bulk GSEA
# # gmtfile with pathways or gene sets
gmtfile = "MSigDB_Hallmark_v7.5.1.gmt"
# # One could also use placeholders for the cases.
# # Currently only cluster is supported. One could use `{cluster}` or
# # `{ident}` to denote the clusters.
[ScFGSEA.envs.cases.Tumor_vs_Normal]
# # Add columns to meta data
mutaters = {Group = "if_else (region == 'Tumor', 'Tumor', 'Normal')"}
# # The first group of cells
"ident.1" = "Tumor"
# # The second group of cells
"ident.2" = "Normal"
# # From which column?
"group.by" = "Group"
# # Could also focus on a subset of the cells
# filter = "region != 'other'"
```


[1]: https://github.com/ctlab/fgsea
