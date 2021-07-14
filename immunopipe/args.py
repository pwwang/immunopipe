"""A pipeline for integrative analysis for scTCR- and scRNA-seq data"""
import sys

import toml
from diot import Diot
from pipen_args import Args

from .defaults import PIPELINE_OPTION_TITLE, PIPELINE_DESCRIPTION
from .commands import full_opts

params = Args(pipen_opt_group=PIPELINE_OPTION_TITLE,
              desc=[PIPELINE_DESCRIPTION,
                    'See full options by `{prog} full-opts`'])
opts_to_hide = (
    'profile',
    'loglevel',
    'dirsig',
    'cache',
    'forks',
    'workdir',
    'error_strategy',
    'num_retries',
    'submission_batch',
    'envs',
    'scheduler',
    'scheduler_opts',
    'plugins',
    'plugin_opts'
)
for opt in opts_to_hide:
    params.get_param(opt).show = False

params.add_command(
    'full-opts',
    'Show full options',
    help_on_void=False
)
params.add_param(
    'c,config',
    type='path',
    desc='The configuration file, in toml.'
)
params.add_param(
    'samples',
    type='path',
    desc='The file defines sample information.',
    required=True
)
params.add_param(
    'meta',
    type='path',
    desc=(
        'The meta information about the patients, '
        'using `Patient` from `samples` as key.'
    ),
    required=True
)
params.add_param(
    'perl',
    type=str,
    default='perl',
    desc='Path to perl'
)
params.add_param(
    'meta_excl',
    type=str,
    default='',
    desc='Values to be excluded in metadata grouping analysis (i.e. controls)'
)
params.add_param(
    'datadir',
    type='path',
    desc=(
        'The directory of the samples listed in samples file. If the samples '
        'are given with absolute paths, this options is ignored'
    )
)
params.add_param(
    'ncores',
    type=int,
    default=1,
    desc=('Number of cores to use for wherever parallelization applies.')
)
params.add_param(
    'rscript',
    show=False,
    default='Rscript',
    desc='The path to Rscript to run processes with scripts in R.'
)
params.add_param(
    'extra_config',
    type='file',
    show=True,
    required=False,
    desc='Extra configurations in TOML for individual processes.'
)

params.param_groups[PIPELINE_OPTION_TITLE] = params.param_groups.pop(
    PIPELINE_OPTION_TITLE
)

args = params.parse(ignore_errors=True)
if args.__command__ == 'full-opts':
    full_opts(params)
    sys.exit(0)

else:
    if args.config and args.config.is_file():
        params.from_file(args.config, force=True)
    args = params.parse()

# parse extra config
extra_config = Diot(DE=None, PatientSamples=None)
if args.extra_config and args.extra_config.is_file():
    with open(args.extra_config):
        xconfig = toml.load(args.extra_config)
    if 'DE' in xconfig:
        extra_config['DE'] = toml.dumps(xconfig['DE'])
    if 'PatientSamples' in xconfig:
        extra_config['PatientSamples'] = toml.dumps(xconfig['PatientSamples'])
    if 'TCellClusterGeneExprs' in xconfig:
        extra_config['TCellClusterGeneExprs'] = toml.dumps(
            xconfig['TCellClusterGeneExprs']
        )
    if 'TCellClusters' in xconfig:
        extra_config['TCellClusters'] = toml.dumps(
            xconfig['TCellClusters']
        )

args.extra_config = extra_config
