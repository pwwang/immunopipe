## Find markers through pair-wised comparisons

`immunopipe` is able to find markers for selected group of cells, not only the cells from clusters, using `Seruat`'s `FindMarkers()`

You can use the metadata to select the cells. For example, to find markers for cells from post-CART treatment samples, compared to the cells from pre-CART samples, one can do:

```toml
# [[]] indicates that we can have multiple cases
[[MARKERS_FINDER]]
# name of this case
name = "Post-CART vs pre-CART treatment"
[MARKERS_FINDER.design]
# The design of comparison
Post_vs_Pre = ["Post", "Pre"]
# Define the subsetings (how clones are selected)
[MARKERS_FINDER.subsetting.Post]
"by.meta" = { TimePoint = "Post" }
[MARKERS_FINDER.subsetting.Pre]
"by.meta" = { Source = "Pre" }
```

To check the filters to select cells, check out:
```shell
pipen run tcr ImmunarchFilter
```

## Investigate markers across comparisons

If you have more than 1 comparison (say `Post_vs_Pre` and `Post_vs_Early` (samples sequenced at an eariler time point)), and you want to investigate how the markers are overlapping across the comparisons, for example, if there are any shared markers between the two comparisons, your configurations could look like:

```toml
[[MARKERS_FINDER]]
name = "Post-CART vs pre-CART vs early"
[MARKERS_FINDER.design]
Post_vs_Pre = ["Post", "Pre"]
Post_vs_Early = ["Post", "Early"]
[MARKERS_FINDER.overlap]
# Find the overlapping markers between the two comparisons
PostPre_vs_PostEarly = ["Post_vs_Pre", "Post_vs_Early"]
[MARKERS_FINDER.subsetting.Post]
"by.meta" = { TimePoint = "Post" }
[MARKERS_FINDER.subsetting.Pre]
"by.meta" = { Source = "Pre" }
[MARKERS_FINDER.subsetting.Early]
"by.meta" = { Source = "Early" }
```

## Find meta-markers from more than 2 groups

For 2-group or pair-wised comparisons, `Seurat::FindMarkers()` can do the job. However, if we want to look for some markers that are determinants across 3 groups, we need to perform additional analysis. `immunopipe` uses `ANOVA` to do it. The configurations is like:

```toml
[[MARKERS_FINDER]]
name = "Post-CART vs pre-CART vs early"
[MARKERS_FINDER.design]
Post_vs_Pre = ["Post", "Pre"]
Post_vs_Early = ["Post", "Early"]
[MARKERS_FINDER.meta]
Post_vs_Pre_vs_Early = ["Post", "Pre", "Early"]
[MARKERS_FINDER.subsetting.Post]
"by.meta" = { TimePoint = "Post" }
[MARKERS_FINDER.subsetting.Pre]
"by.meta" = { Source = "Pre" }
[MARKERS_FINDER.subsetting.Early]
"by.meta" = { Source = "Early" }
```
