To investigate the cell composition for each cluster, a radar chart can be drawn to show the composition cells from different groups.

For example, you can check the cell proportion of cluster#0 from pre-CART and post-CART samples. This is the `inter-cluster` way.

There is also another way: `intra-cluster`, that is, to check the proportion of cells from pre-CART samples in cluster#0 out of all cells from pre-CART samples (all clusters)

```toml
[RadarPlots.envs.cases.Tumor_vs_Normal]
# set the limits
breaks = [0, 15, 30]
# # Direction to calculate the percentages
# # inter-cluster: the percentage of the cells in all groups
# # intra-cluster: the percentage of the cells in all clusters
direction = "intra-cluster"
# # Which column to use to separate the cells in different groups
by = "Group"
# # The order of the values in `by`.
order = ["Tumor", "Normal"]
# # Whether the percentages are calculated before or after filtering out the NAs
perc_with_na = false
```
