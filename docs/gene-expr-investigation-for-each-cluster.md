Other than finding markers for each cluster, you can also investigate the expressions of a set of marker genes. You can also limit the expressions to a subset of samples:

```toml
# [[]]: we can do multiple cases
[[GENE_EXPR_INVESTIGATION_CLUSTERS]]
name = "Gene expressions of T-Cell clusters for pre-CART samples"
target = [
    # the pre-CART samples
]
# ncol: split the boxplots in 3 columns
# res: The resolution of the figure
# width: The width of the plot
# height: The height of the plot
plots = { boxplot = {ncol = 3, res = 100, width = 1200, height = 3200} }
```
