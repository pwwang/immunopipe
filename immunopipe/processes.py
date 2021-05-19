from pipen import Proc

from .args import args
from .defaults import SCRIPT_DIR, REPORT_DIR

class LoadSamples(Proc):
    """List sample information"""
    input_keys = 'samplefile:file'
    input = [args.samples]
    output = 'outfile:file:samples.RData'
    script = f'file://{SCRIPT_DIR}/LoadSamples.R'
    lang = args.rscript
    args = dict(
        datadir=args.datadir
    )
    plugin_opts = {
        'report': f'file://{REPORT_DIR}/LoadSamples.svx'
    }

class LoadTCR(Proc):
    """Load TCR data into immunarch object"""
    requires = LoadSamples
    input_keys = 'samples:file'
    output = 'outfile:file:immunarch.RData'
    script = f'file://{SCRIPT_DIR}/LoadTCR.R'
    lang = args.rscript

class CDR3LengthDistribution(Proc):
    """CDR3 length distribution"""
    requires = LoadTCR
    input_keys = 'immdata:file'
    output = 'outdir:file:CDR3Lengths'
    lang = args.rscript
    script = f'file://{SCRIPT_DIR}/CDR3LengthDistribution.R'
    plugin_opts = {
        'report': f'file://{REPORT_DIR}/CDR3LengthDistribution.svx'
    }
