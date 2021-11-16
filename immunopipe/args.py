import sys
from simpleconf import Config
from pipen_args import Args

args = Args(prog="immunopipe")

# Add more arguments if you want
args.add_param(
    "markers_finder",
    default={},
    type="json",
    desc=[
        "Sets of configurations for markers finding. Keys are process names, "
        "and values are:",
        "- `desc`: The description of the process",
        "- `filters`: The filters for the clones/cells "
        "(ident.1 for the cases, ident.2 for the controls)",
        "- `dbs`: The databases to do enrichment analysis against"
    ]
)

config = Config(with_profile=False)

try:
    confidx = sys.argv.index("--config")
    conffile = sys.argv[confidx + 1]
except (ValueError, IndexError):
    pass
else:
    config._load(conffile)
