## Find markers through pair-wised comparisons

`immunopipe` is able to find markers for selected group of cells, not only the cells from clusters, using `Seruat`'s `FindMarkers()`

You can use the metadata to select the cells, or mutate the metadata before selection.

```toml
# # Do marker finding for clones
# # Similar as `ScFGSEA.envs.cases`. "ident.2" is optional. In such a case,
# # The rest of the cells are used as the second group.
[MarkersFinderForClones.envs.cases.Top10_Clones_Tumor]
"ident.1" = "Top10_Clones"
"ident.2" = "Rest_Clones"
"group.by" = "Top10_Clones"
filter = "region == 'Tumor' & !is.na(Clones)"
[MarkersFinderForClones.envs.cases.Top10_Clones_Tumor.mutaters]
".CloneSize" = "group_by(meta, CDR3.aa) |> mutate(n=n()) |> pull('n')"
Top10_Clones = """if_else (
  .CloneSize %in% sort(unique(.CloneSize), decreasing=TRUE)[1:10],
  'Top10_Clones',
  'Rest_Clones'
)"""
```

## Investigate markers across comparisons

If you have more than 1 comparison (say `Post_vs_Pre` and `Post_vs_Early` (samples sequenced at an eariler time point)), and you want to investigate how the markers are overlapping across the comparisons, for example, if there are any shared markers between the two comparisons, your configurations could look like:

```toml
[MarkersOverlapping.envs]
# # You can also inspect how the markers are overlapping (requiring >=2 cases)
[MarkersOverlapping.envs.cases.Top_10Clones_Makers_Tumor_Normal]
# # Markers finder cases
overlaps = ["Top10_Clones_Tumor", "Top10_Clones_Normal"]
# # Devpars for venn plot
devpars = {res = 100, height = 1000, width = 1200}
```

## Find meta-markers from more than 2 groups

For 2-group or pair-wised comparisons, `Seurat::FindMarkers()` can do the job. However, if we want to look for some markers that are determinants across 3 groups, we need to perform additional analysis. `immunopipe` uses `ANOVA` to do it. The configurations is like:

```toml
[MetaMarkersForClones.envs]
# # Find markers for multiple groups (> 2)
# # Number of cores to use
ncores = 4
# # cases
[MetaMarkersForClones.envs.cases.Top5_10_Rest]
"group.by" = "CloneGroup"
filter = "!is.na(Clones)"
# # Also, add helper columns
[MetaMarkersForClones.envs.cases.Top5_10_Rest.mutaters]
".CloneSize" = "group_by(meta, CDR3.aa) |> mutate(n=n()) |> pull('n')"
CloneGroup = """case_when (
  .CloneSize %in% sort(unique(.CloneSize), decreasing=TRUE)[1:5] ~ 'Top5_Clones',
  .CloneSize %in% sort(unique(.CloneSize), decreasing=TRUE)[6:10] ~ 'Top6_10_Clones',
  TRUE ~ 'Rest'
)"""
```
