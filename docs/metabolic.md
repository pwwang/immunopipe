Metabolic landscape of single cells in the tumor microenvironment

See:

*Xiao, Zhengtao, Ziwei Dai, and Jason W. Locasale. "Metabolic landscape of the tumor microenvironment at single cell resolution." Nature communications 10.1 (2019): 1-12.*

and

https://github.com/LocasaleLab/Single-Cell-Metabolic-Landscape

## Configurations

Example:

```toml

[METABOLIC]
# # The metabolic pathway in gmt file
gmtfile = "KEGG_metabolism.gmt"
[[METABOLIC.cases]]
# # The name of the case
name = "Case1"
# # How do we group the cells?
[METABOLIC.cases.grouping]
# # Add new columns to the meta.data
# mutaters = {}
# # The columns to group the cells
groupby = "seurat_clusters"
# # How do we subset the data. The imputation will be done in each subset separately
[METABOLIC.cases.subsetting]
# # Add new columns to the meta.data
# mutaters = {}
# # The columns to subset the data
groupby = "region"
# # The alias of the subset working as a prefix to subset names
alias = "SampleType"
# # What kind of comparisons are we doing?
# # It should be the values of subsetting `groupby`s
[METABOLIC.cases.design]
Tumor_vs_Normal = ["Tumor", "Normal adjacent tissue"]

# # More cases
# [[METABOLIC.cases]]
# name = "Case2"
# ...
"""

See also: https://pwwang.github.io/biopipen/pipelines/scrna_metabolic/
