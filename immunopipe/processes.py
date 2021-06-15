import sys
from pipen import Proc

from .args import args
from .defaults import SCRIPT_DIR, REPORT_DIR

class LoadSamples(Proc):
    """List sample information

    Output file has 4 variables:
    - samples: The sample information
    - metadata: The metadata for the samples
    - tcrdir: The directory with TCR files, easier for immunarch to load
    - metagroups: How the meta variables grouped
    """
    input_keys = 'samplefile:file, metafile:file'
    input = [(args.samples, args.meta)]
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

class VJUsage(Proc):
    """V-J usage in circular plot"""
    requires = LoadSamples, LoadTCR
    input_keys = 'samples:file, immdata:file'
    output = 'outdir:file:VJUsage'
    lang = args.rscript
    script = f'file://{SCRIPT_DIR}/VJUsage.R'
    plugin_opts = {
        'report': f'file://{REPORT_DIR}/VJUsage.svx'
    }
    args = {
        'vdjtools_patch': f'{SCRIPT_DIR}/vdjtools-patch.sh',
        'vdjtools': 'vdjtools',
        'ncores': args.ncores
    }

class RepertoireOverlap(Proc):
    """Repertorie overlaps between samples"""
    requires = LoadTCR
    input_keys = 'immdata:file'
    output = 'outdir:file:RepertoireOverlap'
    lang = args.rscript
    script = f'file://{SCRIPT_DIR}/RepertoireOverlap.R'
    plugin_opts = {
        'report': f'file://{REPORT_DIR}/RepertoireOverlap.svx'
    }

class BasicStatistics(Proc):
    """Basic statistics and clonality"""
    requires = LoadSamples, LoadTCR
    input_keys = 'samples:file, immdata:file'
    output = 'outdir:file:BasicStats'
    lang = args.rscript
    script = f'file://{SCRIPT_DIR}/BasicStatistics.R'
    args = {'exclude': args.meta_excl}
    plugin_opts = {
        'report': f'file://{REPORT_DIR}/BasicStatistics.svx'
    }

class LoadTCRForIntegration(Proc):
    """Load TCR data into R object for futher RNA data integration"""
    requires = LoadSamples
    input_keys = 'samples:file'
    output = 'outdir:file:TCR-counts'
    lang = sys.executable
    script = f'file://{SCRIPT_DIR}/TCR-counts/LoadTCRForIntegration.py'
    args = {
        'perl': args.perl,
        'rscript': args.rscript,
        'master_loader': f'{SCRIPT_DIR}/TCR-counts/assign-TCR-clonotypes.pl',
        'count_loader': f'{SCRIPT_DIR}/TCR-counts/TCR-counts.R',
    }

class ResidencyColors(Proc):
    """Define residency colors"""
    requires = LoadSamples, LoadTCRForIntegration
    input_keys = 'samples:file, tcr_counts:file'
    output = 'outdir:file:ResidencyColors'
    lang = args.rscript
    script = f'file://{SCRIPT_DIR}/ResidencyColors.R'
    args = {'colpat': f'{SCRIPT_DIR}/utils/color-palettes.R'}

class ResidencyPlots(Proc):
    """Clonotype residency plots"""
    requires = LoadSamples, LoadTCRForIntegration, ResidencyColors
    input_keys = 'samples:file, tcr_counts:file, rc_colors:file'
    output = 'outdir:file:ResidencyPlots'
    lang = args.rscript
    script = f'file://{SCRIPT_DIR}/ResidencyPlots.R'
    plugin_opts = {
        'report': f'file://{REPORT_DIR}/ResidencyPlots.svx'
    }
