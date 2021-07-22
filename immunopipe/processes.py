import sys
from pipen import Proc

from .args import args
from .defaults import SCRIPT_DIR, REPORT_DIR

class SampleInfo(Proc):
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
    script = f'file://{SCRIPT_DIR}/SampleInfo.R'
    lang = args.rscript
    args = dict(
        datadir=args.datadir
    )
    plugin_opts = {
        'report': f'file://{REPORT_DIR}/SampleInfo.svx'
    }

class LoadTCR(Proc):
    """Load TCR data into immunarch object"""
    requires = SampleInfo
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
    requires = SampleInfo, LoadTCR
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
    requires = SampleInfo, LoadTCR
    input_keys = 'samples:file, immdata:file'
    output = 'outdir:file:BasicStats'
    lang = args.rscript
    script = f'file://{SCRIPT_DIR}/BasicStatistics.R'
    args = {
        'exclude': args.meta_excl,
        'config': args.extra_config.BasicStatistics
    }
    plugin_opts = {
        'report': f'file://{REPORT_DIR}/BasicStatistics.svx'
    }

class LoadTCRForIntegration(Proc):
    """Load TCR data into R object for futher RNA data integration"""
    requires = SampleInfo
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
    requires = SampleInfo, LoadTCRForIntegration
    input_keys = 'samples:file, tcr_counts:file'
    output = 'outdir:file:ResidencyColors'
    lang = args.rscript
    script = f'file://{SCRIPT_DIR}/ResidencyColors.R'
    args = {'colpat': f'{SCRIPT_DIR}/utils/color-palettes.R'}

class ResidencyPlots(Proc):
    """Clonotype residency plots"""
    requires = SampleInfo, LoadTCRForIntegration, ResidencyColors
    input_keys = 'samples:file, tcr_counts:file, rc_colors:file'
    output = 'outdir:file:ResidencyPlots'
    lang = args.rscript
    script = f'file://{SCRIPT_DIR}/ResidencyPlots.R'
    plugin_opts = {
        'report': f'file://{REPORT_DIR}/ResidencyPlots.svx'
    }

class LoadExprData(Proc):
    """Load expression data"""
    # Matrix-0 directory
    requires = SampleInfo
    input_keys = 'samples:file'
    output = 'outdir:file:LoadExprData'
    lang = args.rscript
    script = f'file://{SCRIPT_DIR}/LoadExprData.R'
    args = {'ncores': args.ncores}

class DEAnalysis(Proc):
    """Differential gene expression analysis between different groups"""
    if args.extra_config.DE:
        requires = SampleInfo, LoadExprData

    input_keys = 'samples:file, exprdir:file'
    output = 'outdir:file:DEAnalysis'
    lang = args.rscript
    script = f'file://{SCRIPT_DIR}/DEAnalysis.R'
    args = {'ncores': args.ncores, 'config': args.extra_config.DE}
    plugin_opts = {
        'report': f'file://{REPORT_DIR}/DEAnalysis.svx'
    }

class CrossSampleClonotypeComparison(Proc):
    """Clonetype changes across samples of the same patient"""
    if args.extra_config.PatientSamples:
        requires = LoadTCR

    input_keys = 'immdata:file'
    output = 'outdir:file:CrossSampleClonotypeComparison'
    lang = args.rscript
    script = f'file://{SCRIPT_DIR}/CrossSampleClonotypeComparison.R'
    args = {'ncores': args.ncores, 'config': args.extra_config.PatientSamples}
    plugin_opts = {
        'report': f'file://{REPORT_DIR}/CrossSampleClonotypeComparison.svx'
    }

class DEAnalysisChangedClonotypes(Proc):
    """Differential expression analysis for changed clonotypes"""
    if args.extra_config.PatientSamples:
        requires = (
            SampleInfo, LoadTCR, LoadExprData, CrossSampleClonotypeComparison
        )

    input_keys = 'samples:file, immdata:file, exprdir:file, ccdir:file'
    output = 'outdir:file:DEAnalysisChangedClonotypes'
    lang = args.rscript
    script = f'file://{SCRIPT_DIR}/DEAnalysisChangedClonotypes.R'
    args = {
        'ncores': args.ncores,
        'config': args.extra_config.PatientSamples
    }
    plugin_opts = {
        'report': f'file://{REPORT_DIR}/DEAnalysisChangedClonotypes.svx'
    }

class IntegrateTCRExprData(Proc):
    """Integrate TCR and Expression data"""
    # Directories: ,TCR-counts, Matrices-0
    requires = SampleInfo, LoadTCRForIntegration, LoadExprData
    input_keys = 'samples:file, tcr_counts:file, exprdir:file'
    output = 'outdir:file:IntegrateTCRExprData'
    cache = 'force'
    lang = args.rscript
    script = f'file://{SCRIPT_DIR}/Seurat-0/IntegrateTCRExprData.R'
    args = {
        'ncores': args.ncores,
        'seurate_source': f'{SCRIPT_DIR}/Seurat-0/seurate-source.R',
    }

class SeparateTnonTCells(Proc):
    """Separate T and non-T cells"""
    requires = IntegrateTCRExprData
    input_keys = 'itgdir:file'
    output = 'outdir:file:SeparateTnonTCells'
    lang = args.rscript
    script = f'file://{SCRIPT_DIR}/SeparateTnonTCells.R'
    args = { 'ncores': args.ncores }
    plugin_opts = {
        'report': f'file://{REPORT_DIR}/SeparateTnonTCells.svx'
    }

class ClusterTCells(Proc):
    """Cluster T cells"""
    requires = SampleInfo, LoadTCRForIntegration, LoadExprData, SeparateTnonTCells
    input_keys = 'samples:file, tcr_counts:file, exprdir:file, septdir:file'
    output = 'outdir:file:ClusterTCells'
    cache = 'force'
    lang = args.rscript
    script = f'file://{SCRIPT_DIR}/ClusterTCells.R'
    args = {
        'ncores': args.ncores,
        'seurate_source': f'{SCRIPT_DIR}/Seurat-0/seurate-source.R',
    }
    plugin_opts = {
        'report': f'file://{REPORT_DIR}/ClusterTCells.svx'
    }

class TCellClusterGeneExprs(Proc):
    """Expressions in different T-cell clusters for genes of interest"""
    if args.extra_config.TCellClusterGeneExprs:
        requires = ClusterTCells

    input_keys = 'cldir:file'
    output = 'outdir:file:TCellClusterGeneExprs'
    lang = args.rscript
    script = f'file://{SCRIPT_DIR}/TCellClusterGeneExprs.R'
    args = {
        'genes': args.extra_config.TCellClusterGeneExprs,
        'tclusters': args.extra_config.TCellClusters
    }
    plugin_opts = {
        'report': f'file://{REPORT_DIR}/TCellClusterGeneExprs.svx'
    }
