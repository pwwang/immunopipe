To investigate the cell composition for each cluster, a radar chart can be drawn to show the composition cells from different groups.

For example, you can check the cell proportion of cluster#0 from pre-CART and post-CART samples. This is the `inter-cluster` way.

There is also another way: `intra-cluster`, that is, to check the proportion of cells from pre-CART samples in cluster#0 out of all cells from pre-CART samples (all clusters)

```toml
# [[ ]] indicates that we can do multiple cases
[[RADAR_PLOTS]]
# name of the case
name = "Cell proportion of control samples"
# set the limits
breaks = [0, 15, 30]
direction = "intra-cluster"
[RADAR_PLOTS.filters.BM]
# Select cells of the first group
"by.meta" = { TimePoint = "pre-CART" }
[RADAR_PLOTS.filters.PB]
# Select cells of the second group
"by.meta" = { TimePoint = "post-CART" }
```

See https://pwwang.github.io/immunopipe/markers-finder/ for configurations to select cells