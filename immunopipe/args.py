
from pipen_args import Args, config  # noqa: F401

args = Args(prog="immunopipe")

# Add more arguments if you want
args.add_param(
    "METABOLIC",
    type="ns",
    desc="Configurations for metabolic landscape analysis",
)

args.add_param(
    "METABOLIC.gmtfile",
    type="file",
    desc="The metabolic pathway in GMT file",
)

args.add_param(
    "METABOLIC.cases",
    type="list:json",
    desc=(
        'The list of cases for metabolic landscape analysis. '
        'Each case is a dictionary with the following keys: '
        '"name", "grouping", "subsetting" and "design". ',
        '"name" is the name of the case.',
        '"grouping" is the grouping of the case, with "mutaters" to add '
        'new columns to the metadata and "groupby" the column to group '
        'the cells.',
        '"subsetting" is the subsetting of the case, with "mutaters" to '
        'add new columns to the metadata, "groupby" the column to '
        'subset the cells and "alias" as prefix to add to subset names.',
        '"design" is the design of the case to define the comparisons. '
    )
)
