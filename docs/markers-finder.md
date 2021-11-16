`immunopipe` is able to find markers for selected group of cells, not only the cells from clusters, using `Seruat`'s `FindMarkers()`

You can use the metadata to select the cells. For example, to find markers for cells from post-CART treatment samples, compared to the cells from pre-CART samples, one can do:

```toml
# [[]] indicates that we can have multiple cases
[[MARKERS_FINDER.filters]]
# name of this case
name = "[pre-CART BM DR] Markers of top 20 TCR clones"
[MARKERS_FINDER.filters."ident.1"]
# Select cells for the first group
"by.meta" = { TimePoint = "pre-CART" }
[MARKERS_FINDER.filters."ident.2"]
# Select cells for the second(control) group
"by.meta" = { TimePoint = "post-CART" }
```

To check the filters to select cells, check out:
```shell
pipen run tcr ImmunarchFilter
```
