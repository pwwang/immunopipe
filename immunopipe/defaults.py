"""Default settings"""
from pathlib import Path
from pyparam import defaults as pyparam_defaults

# adjust cli help page
pyparam_defaults.CONSOLE_WIDTH = 100
pyparam_defaults.HELP_OPTION_WIDTH = 36

PIPELINE_OPTION_TITLE = 'PIPELINE_OPTIONS'

PIPELINE_DESCRIPTION = (
    'A pipeline for integrative analysis for scTCR- and scRNA-seq data'
)

# other constants
SCRIPT_DIR = Path(__file__).parent.joinpath('scripts').resolve()
REPORT_DIR = Path(__file__).parent.joinpath('reports').resolve()
