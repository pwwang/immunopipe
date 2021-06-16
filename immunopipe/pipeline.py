from pipen import Pipen

from .defaults import PIPELINE_DESCRIPTION
from .args import args
from .processes import SampleInfo

def pipeline():
    """Get the pipeline object"""

    return Pipen(
        'immunopipe',
        desc=PIPELINE_DESCRIPTION,
        # plugin_opts={'report_debug': True}
    ).starts(SampleInfo)
