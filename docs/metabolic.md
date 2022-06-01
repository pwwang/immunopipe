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
# The name of the case, used in reports, optional
name = "BM samples"

[METABOLIC.cases.grouping]
# How should we do the groupings
# Idents: Use the Seurat clustering
groupby = "Idents"

[METABOLIC.cases.subsetting]
# Define the subsets.
groupby = "Group"

[METABOLIC.cases.subsetting.mutaters]
Group = """
    case_when(
        Source == "BM" & TimePoint == "Post-CART" & Response == "DR" ~ "Post_DR",
        Source == "BM" & TimePoint == "Post-CART" & Response == "PD" ~ "Post_PD",
        Source == "BM" & TimePoint == "Post-CART" & Response == "CR" ~ "Post_CR",
        Source == "BM" & Response == "CNTRL" ~ "CNTRL",
        TRUE ~ NA_character_
    )
"""

[METABOLIC.cases.design]
# Defines how we should do the comparisons and genes are ranked for each design
Post_DR-vs-CNTRL = ["Post_DR", "CNTRL"]
Post_PD-vs-CNTRL = ["Post_PD", "CNTRL"]
Post_CR-vs-CNTRL = ["Post_CR", "CNTRL"]

[[METABOLIC.cases]]
# more cases
```
