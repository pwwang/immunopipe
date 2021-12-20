Metabolic landscape of single cells in the tumor microenvironment

See:

*Xiao, Zhengtao, Ziwei Dai, and Jason W. Locasale. "Metabolic landscape of the tumor microenvironment at single cell resolution." Nature communications 10.1 (2019): 1-12.*

and

https://github.com/LocasaleLab/Single-Cell-Metabolic-Landscape

## Configurations

Example:

```toml

[METABOLIC]
# The GMT file containing the metabolic pathways
gmtfile = "KEGG_metabolism.gmt"

[[METABOLIC.cases]]
# The name of the case
name = "BM samples"
# How should we do the groupings
# Idents: Use the Seurat clustering
grouping = "Idents"
# What's the name to show for the grouping on plots/reports
grouping_name = "Cluster"
[METABOLIC.cases.design]
# Defines how we should do the comparisons and genes are ranked for each design
Post_DR-vs-CNTRL = ["Post_DR", "CNTRL"]
Post_PD-vs-CNTRL = ["Post_PD", "CNTRL"]
Post_CR-vs-CNTRL = ["Post_CR", "CNTRL"]
[METABOLIC.cases.subsetting]
# Define the subsets.
# The conditions will be passed to `subset(<seruat_obj>, subset = ...)`
Post_DR = 'Source == "BM" & TimePoint == "Post-CART" & Response == "DR"'
Post_PD = 'Source == "BM" & TimePoint == "Post-CART" & Response == "PD"'
Post_CR = 'Source == "BM" & TimePoint == "Post-CART" & Response == "CR"'
CNTRL = 'Source == "BM" & Response == "CNTRL"'

[[METABOLIC.cases]]
# more cases
```

If you want to subset the cells by counts, you can use `ImmunarchFilter` by specifying `subset_using = "immunarch"` for each cases. Then for each subset, you can use `by.count` or `by.meta` to select the cells.

See details by running this command:

```shell
pipen run tcr ImmunarchFilter
```
