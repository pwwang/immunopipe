import pytest  # noqa: F401

from .conftest import dry_run


@pytest.mark.forked
def test_route_sampleinfo_minimal(tmp_path):
    """Test route when only SampleInfo is provided."""
    config = tmp_path / "route_sampleinfo_minimal.config.toml"
    config.write_text(
        f"""
        workdir = "{tmp_path}/workdir"
        outdir = "{tmp_path}/outdir"

        [SampleInfo]
        in.infile = "{__file__}"
    """
    )

    output = dry_run(f"@{config}")
    assert "SampleInfo: <<< [START]" in output
    assert "SampleInfo: >>> ['SeuratPreparing']" in output
    assert "SeuratPreparing: <<< ['SampleInfo']" in output
    assert "SeuratPreparing: >>> ['SeuratClustering']" in output
    assert "SeuratClustering: <<< ['SeuratPreparing']" in output
    assert "SeuratClustering: >>> ['ClusterMarkers', 'SeuratClusterStats']" in output
    assert "SeuratClusterStats: <<< ['SeuratClustering']" in output
    assert "SeuratClusterStats: >>> [END]" in output
    assert "ClusterMarkers: <<< ['SeuratClustering']" in output
    assert "ClusterMarkers: >>> [END]" in output


@pytest.mark.forked
def test_route_sampleinfo_torbselection(tmp_path):
    """Test route when SampleInfo with TCR/BCR data is provided."""
    config = tmp_path / "route_sampleinfo_torbselection.config.toml"
    sampleinfo = tmp_path / "sampleinfo_torbselection.txt"
    sampleinfo.write_text(
        "SampleID\tTCRData\tRNAData\n" "S1\tdata/S1_tcr.fastq\tdata/S1_rna.fastq\n"
    )
    config.write_text(
        f"""
        workdir = "{tmp_path}/workdir"
        outdir = "{tmp_path}/outdir"

        [SampleInfo]
        in.infile = "{sampleinfo}"

        [TOrBCellSelection]
    """
    )

    output = dry_run(f"@{config}")
    assert "SampleInfo: <<< [START]" in output
    assert "SampleInfo: >>> ['ScRepLoading', 'SeuratPreparing']" in output
    assert "ScRepLoading: <<< ['SampleInfo']" in output
    assert (
        "ScRepLoading: >>> ['TOrBCellSelection', 'ScRepCombiningExpression']" in output
    )
    assert (
        "ScRepCombiningExpression: <<< ['ScRepLoading', 'SeuratClustering']" in output
    )
    assert (
        "ScRepCombiningExpression: >>> ['SeuratClusterStats', 'ClonalStats']" in output
    )
    assert "SeuratClusterStats: <<< ['ScRepCombiningExpression']" in output
    assert "SeuratClusterStats: >>> [END]" in output
    assert "SeuratPreparing: <<< ['SampleInfo']" in output
    assert "SeuratPreparing: >>> ['SeuratClusteringOfAllCells']" in output
    assert "SeuratClusteringOfAllCells: <<< ['SeuratPreparing']" in output
    assert (
        "SeuratClusteringOfAllCells: >>> "
        "['ClusterMarkersOfAllCells', 'TOrBCellSelection']" in output
    )
    assert "ClusterMarkersOfAllCells: <<< ['SeuratClusteringOfAllCells']" in output
    assert "ClusterMarkersOfAllCells: >>> [END]" in output
    assert (
        "TOrBCellSelection: <<< ['SeuratClusteringOfAllCells', 'ScRepLoading']"
        in output
    )
    assert "TOrBCellSelection: >>> ['SeuratClustering']" in output
    assert "SeuratClustering: <<< ['TOrBCellSelection']" in output
    assert (
        "SeuratClustering: >>> ['ClusterMarkers', 'ScRepCombiningExpression']" in output
    )
    assert "ClusterMarkers: <<< ['SeuratClustering']" in output
    assert "ClusterMarkers: >>> [END]" in output
    assert "ClonalStats: <<< ['ScRepCombiningExpression']" in output
    assert "ClonalStats: >>> [END]" in output


@pytest.mark.forked
def test_route_loadingrnafromseurat_only(tmp_path):
    """Test route when only LoadingRNAFromSeurat is provided."""
    config = tmp_path / "route_loadingrnafromseurat_only.config.toml"
    seuratfile = tmp_path / "seurat_object.rds"
    seuratfile.write_text("MOCK SEURAT OBJECT CONTENT")
    config.write_text(
        f"""
        workdir = "{tmp_path}/workdir"
        outdir = "{tmp_path}/outdir"

        [LoadingRNAFromSeurat]
        in.infile = "{seuratfile}"
    """
    )

    output = dry_run(f"@{config}")
    assert "LoadingRNAFromSeurat: <<< [START]" in output
    assert "LoadingRNAFromSeurat: >>> ['SeuratPreparing']" in output
    assert "SeuratPreparing: <<< ['LoadingRNAFromSeurat']" in output
    assert "SeuratPreparing: >>> ['SeuratClustering']" in output
    assert "SeuratClustering: <<< ['SeuratPreparing']" in output
    assert "SeuratClustering: >>> ['ClusterMarkers', 'SeuratClusterStats']" in output
    assert "SeuratClusterStats: <<< ['SeuratClustering']" in output
    assert "SeuratClusterStats: >>> [END]" in output
    assert "ClusterMarkers: <<< ['SeuratClustering']" in output
    assert "ClusterMarkers: >>> [END]" in output


@pytest.mark.forked
def test_route_loadingrnafromseurat_only_prepared(tmp_path):
    """Test route when only LoadingRNAFromSeurat is provided."""
    config = tmp_path / "route_loadingrnafromseurat_only_prepared.config.toml"
    seuratfile = tmp_path / "seurat_object.rds"
    seuratfile.write_text("MOCK SEURAT OBJECT CONTENT")
    config.write_text(
        f"""
        workdir = "{tmp_path}/workdir"
        outdir = "{tmp_path}/outdir"

        [LoadingRNAFromSeurat]
        in.infile = "{seuratfile}"

        [LoadingRNAFromSeurat.envs]
        prepared = true
    """
    )

    output = dry_run(f"@{config}")
    assert "LoadingRNAFromSeurat: <<< [START]" in output
    assert "LoadingRNAFromSeurat: >>> ['SeuratClustering']" in output
    assert "SeuratClustering: <<< ['LoadingRNAFromSeurat']" in output
    assert "SeuratClustering: >>> ['ClusterMarkers', 'SeuratClusterStats']" in output
    assert "SeuratClusterStats: <<< ['SeuratClustering']" in output
    assert "SeuratClusterStats: >>> [END]" in output
    assert "ClusterMarkers: <<< ['SeuratClustering']" in output
    assert "ClusterMarkers: >>> [END]" in output


@pytest.mark.forked
def test_route_loadingrnafromseurat_only_clustered(tmp_path):
    """Test route when only LoadingRNAFromSeurat is provided."""
    config = tmp_path / "route_loadingrnafromseurat_only_clustered.config.toml"
    seuratfile = tmp_path / "seurat_object.rds"
    seuratfile.write_text("MOCK SEURAT OBJECT CONTENT")
    config.write_text(
        f"""
        workdir = "{tmp_path}/workdir"
        outdir = "{tmp_path}/outdir"

        [LoadingRNAFromSeurat]
        in.infile = "{seuratfile}"

        [LoadingRNAFromSeurat.envs]
        clustered = true
    """
    )

    output = dry_run(f"@{config}")
    assert "LoadingRNAFromSeurat: <<< [START]" in output
    assert (
        "LoadingRNAFromSeurat: >>> ['ClusterMarkers', 'SeuratClusterStats']" in output
    )
    assert "SeuratClusterStats: <<< ['LoadingRNAFromSeurat']" in output
    assert "SeuratClusterStats: >>> [END]" in output
    assert "ClusterMarkers: <<< ['LoadingRNAFromSeurat']" in output
    assert "ClusterMarkers: >>> [END]" in output


@pytest.mark.forked
def test_route_loadingrnafromseurat_torbselection(tmp_path):
    """Test route when LoadingRNAFromSeurat with TCR/BCR data is provided."""
    config = tmp_path / "route_loadingrnafromseurat_torbselection.config.toml"
    seuratfile = tmp_path / "seurat_object.rds"
    seuratfile.write_text("MOCK SEURAT OBJECT CONTENT")
    config.write_text(
        f"""
        workdir = "{tmp_path}/workdir"
        outdir = "{tmp_path}/outdir"

        [LoadingRNAFromSeurat]
        in.infile = "{seuratfile}"

        [TOrBCellSelection]
    """
    )

    output = dry_run(f"@{config}")
    assert "LoadingRNAFromSeurat: <<< [START]" in output
    assert "LoadingRNAFromSeurat: >>> ['SeuratPreparing']" in output
    assert "SeuratPreparing: <<< ['LoadingRNAFromSeurat']" in output
    assert "SeuratPreparing: >>> ['SeuratClusteringOfAllCells']" in output
    assert "SeuratClusteringOfAllCells: <<< ['SeuratPreparing']" in output
    assert (
        "SeuratClusteringOfAllCells: >>> "
        "['ClusterMarkersOfAllCells', 'TOrBCellSelection']" in output
    )
    assert "ClusterMarkersOfAllCells: <<< ['SeuratClusteringOfAllCells']" in output
    assert "ClusterMarkersOfAllCells: >>> [END]" in output
    assert "TOrBCellSelection: <<< ['SeuratClusteringOfAllCells']" in output
    assert "TOrBCellSelection: >>> ['SeuratClustering']" in output
    assert "SeuratClustering: <<< ['TOrBCellSelection']" in output
    assert "SeuratClustering: >>> ['ClusterMarkers', 'SeuratClusterStats']" in output
    assert "SeuratClusterStats: <<< ['SeuratClustering']" in output
    assert "SeuratClusterStats: >>> [END]" in output
    assert "ClusterMarkers: <<< ['SeuratClustering']" in output
    assert "ClusterMarkers: >>> [END]" in output


@pytest.mark.forked
def test_route_loadingrnafromseurat_and_sampleinfo(tmp_path):
    """Test route when both LoadingRNAFromSeurat and SampleInfo are provided."""
    config = tmp_path / "route_loadingrnafromseurat_and_sampleinfo.config.toml"
    seuratfile = tmp_path / "seurat_object.rds"
    seuratfile.write_text("MOCK SEURAT OBJECT CONTENT")
    sampleinfo = tmp_path / "sampleinfo.txt"
    sampleinfo.write_text(
        "SampleID\tTCRData\tRNAData\n" "S1\tdata/S1_tcr.fastq\tdata/S1_rna.fastq\n"
    )
    config.write_text(
        f"""
        workdir = "{tmp_path}/workdir"
        outdir = "{tmp_path}/outdir"

        [LoadingRNAFromSeurat]
        in.infile = "{seuratfile}"

        [SampleInfo]
        in.infile = "{sampleinfo}"
    """
    )

    output = dry_run(f"@{config}")
    assert "SampleInfo: <<< [START]" in output
    assert "SampleInfo: >>> ['ScRepLoading']" in output
    assert "LoadingRNAFromSeurat: <<< [START]" in output
    assert "LoadingRNAFromSeurat: >>> ['SeuratPreparing']" in output
    assert "ScRepLoading: <<< ['SampleInfo']" in output
    assert "ScRepLoading: >>> ['ScRepCombiningExpression']" in output
    assert "SeuratPreparing: <<< ['LoadingRNAFromSeurat']" in output
    assert "SeuratPreparing: >>> ['SeuratClustering']" in output
    assert "SeuratClustering: <<< ['SeuratPreparing']" in output
    assert (
        "SeuratClustering: >>> ['ClusterMarkers', 'ScRepCombiningExpression']" in output
    )
    assert (
        "ScRepCombiningExpression: <<< ['ScRepLoading', 'SeuratClustering']" in output
    )
    assert (
        "ScRepCombiningExpression: >>> ['SeuratClusterStats', 'ClonalStats']" in output
    )
    assert "SeuratClusterStats: <<< ['ScRepCombiningExpression']" in output
    assert "SeuratClusterStats: >>> [END]" in output
    assert "ClusterMarkers: <<< ['SeuratClustering']" in output
    assert "ClusterMarkers: >>> [END]" in output
    assert "ClonalStats: <<< ['ScRepCombiningExpression']" in output
    assert "ClonalStats: >>> [END]" in output


@pytest.mark.forked
def test_route_loadingrnafromseurat_and_sampleinfo_torbselection(tmp_path):
    """Test route when both LoadingRNAFromSeurat and SampleInfo with TCR/BCR data are
    provided."""
    config = (
        tmp_path / "route_loadingrnafromseurat_and_sampleinfo_torbselection.config.toml"
    )
    seuratfile = tmp_path / "seurat_object.rds"
    seuratfile.write_text("MOCK SEURAT OBJECT CONTENT")
    sampleinfo = tmp_path / "sampleinfo_torbselection.txt"
    sampleinfo.write_text(
        "SampleID\tTCRData\tRNAData\n" "S1\tdata/S1_tcr.fastq\tdata/S1_rna.fastq\n"
    )
    config.write_text(
        f"""
        workdir = "{tmp_path}/workdir"
        outdir = "{tmp_path}/outdir"

        [LoadingRNAFromSeurat]
        in.infile = "{seuratfile}"

        [SampleInfo]
        in.infile = "{sampleinfo}"

        [TOrBCellSelection]
    """
    )

    output = dry_run(f"@{config}")
    assert "SampleInfo: <<< [START]" in output
    assert "SampleInfo: >>> ['ScRepLoading']" in output
    assert "LoadingRNAFromSeurat: <<< [START]" in output
    assert "LoadingRNAFromSeurat: >>> ['SeuratPreparing']" in output
    assert "ScRepLoading: <<< ['SampleInfo']" in output
    assert (
        "ScRepLoading: >>> ['TOrBCellSelection', 'ScRepCombiningExpression']" in output
    )
    assert "SeuratPreparing: <<< ['LoadingRNAFromSeurat']" in output
    assert "SeuratPreparing: >>> ['SeuratClusteringOfAllCells']" in output
    assert "SeuratClusteringOfAllCells: <<< ['SeuratPreparing']" in output
    assert (
        "SeuratClusteringOfAllCells: >>> "
        "['ClusterMarkersOfAllCells', 'TOrBCellSelection']" in output
    )
    assert "ClusterMarkersOfAllCells: <<< ['SeuratClusteringOfAllCells']" in output
    assert "ClusterMarkersOfAllCells: >>> [END]" in output
    assert (
        "TOrBCellSelection: <<< ['SeuratClusteringOfAllCells', 'ScRepLoading']"
        in output
    )
    assert "TOrBCellSelection: >>> ['SeuratClustering']" in output
    assert "SeuratClustering: <<< ['TOrBCellSelection']" in output
    assert (
        "SeuratClustering: >>> ['ClusterMarkers', 'ScRepCombiningExpression']" in output
    )
    assert (
        "ScRepCombiningExpression: <<< ['ScRepLoading', 'SeuratClustering']" in output
    )
    assert (
        "ScRepCombiningExpression: >>> ['SeuratClusterStats', 'ClonalStats']" in output
    )
    assert "SeuratClusterStats: <<< ['ScRepCombiningExpression']" in output
    assert "SeuratClusterStats: >>> [END]" in output
    assert "ClusterMarkers: <<< ['SeuratClustering']" in output
    assert "ClusterMarkers: >>> [END]" in output
    assert "ClonalStats: <<< ['ScRepCombiningExpression']" in output
    assert "ClonalStats: >>> [END]" in output


@pytest.mark.forked
def test_route_sampleinfo_full(tmp_path):
    """Test route when full SampleInfo is provided."""
    config = tmp_path / "route_sampleinfo_full.config.toml"
    sampleinfo = tmp_path / "sampleinfo_full.txt"
    sampleinfo.write_text(
        "SampleID\tTCRData\tRNAData\n" "S1\tdata/S1_tcr.fastq\tdata/S1_rna.fastq\n"
    )
    config.write_text(
        f"""
        workdir = "{tmp_path}/workdir"
        outdir = "{tmp_path}/outdir"

        [SampleInfo]
        in.infile = "{sampleinfo}"

        [TOrBCellSelection]
        [TopExpressingGenesOfAllCells]
        [ModuleScoreCalculator]
        [SeuratClustering]
        [SeuratMap2Ref]
        [CellTypeAnnotation]
        [TopExpressingGenes]
        [SeuratSubClustering]
        [CDR3Clustering]
        [TESSA]
        [MarkersFinder]
        [ScFGSEA]
        [CellCellCommunication]
        [PseudoBulkDEG]
        [CDR3AAPhyschem]
        [ScrnaMetabolicLandscape]
    """
    )
    output = dry_run(f"@{config}")

    assert "SampleInfo: <<< [START]" in output
    assert "SampleInfo: >>> ['ScRepLoading', 'SeuratPreparing']" in output
    assert "ScRepLoading: <<< ['SampleInfo']" in output
    assert (
        "ScRepLoading: >>> ['TOrBCellSelection', 'ScRepCombiningExpression']" in output
    )
    assert "SeuratPreparing: <<< ['SampleInfo']" in output
    assert "SeuratPreparing: >>> ['SeuratClusteringOfAllCells']" in output
    assert "SeuratClusteringOfAllCells: <<< ['SeuratPreparing']" in output
    assert (
        "SeuratClusteringOfAllCells: >>> ['ClusterMarkersOfAllCells', "
        "'TopExpressingGenesOfAllCells', 'TOrBCellSelection']" in output
    )
    assert "ClusterMarkersOfAllCells: <<< ['SeuratClusteringOfAllCells']" in output
    assert "ClusterMarkersOfAllCells: >>> [END]" in output
    assert "TopExpressingGenesOfAllCells: <<< ['SeuratClusteringOfAllCells']" in output
    assert "TopExpressingGenesOfAllCells: >>> [END]" in output
    assert (
        "TOrBCellSelection: <<< ['SeuratClusteringOfAllCells', 'ScRepLoading']"
        in output
    )
    assert "TOrBCellSelection: >>> ['ModuleScoreCalculator']" in output
    assert "ModuleScoreCalculator: <<< ['TOrBCellSelection']" in output
    assert "ModuleScoreCalculator: >>> ['SeuratClustering']" in output
    assert "SeuratMap2Ref: <<< ['CellTypeAnnotation']" in output
    assert "SeuratMap2Ref: >>> ['SeuratSubClustering']" in output
    assert "SeuratClustering: <<< ['ModuleScoreCalculator']" in output
    assert "SeuratClustering: >>> ['CellTypeAnnotation']" in output
    assert "CellTypeAnnotation: <<< ['SeuratClustering']" in output
    assert "CellTypeAnnotation: >>> ['SeuratMap2Ref']" in output
    assert "SeuratSubClustering: <<< ['SeuratMap2Ref']" in output
    assert (
        "SeuratSubClustering: >>> ['ClusterMarkers', 'TopExpressingGenes', "
        "'ScRepCombiningExpression']" in output
    )
    assert (
        "ScRepCombiningExpression: <<< ['ScRepLoading', 'SeuratSubClustering']"
        in output
    )
    assert "ScRepCombiningExpression: >>> ['CDR3Clustering']" in output
    assert "ClusterMarkers: <<< ['SeuratSubClustering']" in output
    assert "ClusterMarkers: >>> [END]" in output
    assert "TopExpressingGenes: <<< ['SeuratSubClustering']" in output
    assert "TopExpressingGenes: >>> [END]" in output
    assert "CDR3Clustering: <<< ['ScRepCombiningExpression']" in output
    assert "CDR3Clustering: >>> ['TESSA']" in output
    assert "TESSA: <<< ['CDR3Clustering']" in output
    assert (
        "TESSA: >>> "
        "['CellCellCommunication', 'SeuratClusterStats', 'ClonalStats', 'ScFGSEA', "
        "'PseudoBulkDEG', 'MarkersFinder', 'CDR3AAPhyschem', 'MetabolicInput']"
        in output
    )
    assert "SeuratClusterStats: <<< ['TESSA']" in output
    assert "SeuratClusterStats: >>> [END]" in output
    assert "CellCellCommunication: <<< ['TESSA']" in output
    assert "CellCellCommunication: >>> ['CellCellCommunicationPlots']" in output
    assert "CellCellCommunicationPlots: <<< ['CellCellCommunication']" in output
    assert "CellCellCommunicationPlots: >>> [END]" in output
    assert "ClonalStats: <<< ['TESSA']" in output
    assert "ClonalStats: >>> [END]" in output
    assert "ScFGSEA: <<< ['TESSA']" in output
    assert "ScFGSEA: >>> [END]" in output
    assert "PseudoBulkDEG: <<< ['TESSA']" in output
    assert "PseudoBulkDEG: >>> [END]" in output
    assert "MarkersFinder: <<< ['TESSA']" in output
    assert "MarkersFinder: >>> [END]" in output
    assert "CDR3AAPhyschem: <<< ['TESSA']" in output
    assert "CDR3AAPhyschem: >>> [END]" in output
    assert "MetabolicInput: <<< ['TESSA']" in output
    assert (
        "MetabolicInput: >>> "
        "['MetabolicPathwayActivity', 'MetabolicPathwayHeterogeneity', "
        "'MetabolicFeatures']" in output
    )
    assert "MetabolicPathwayActivity: <<< ['MetabolicInput']" in output
    assert "MetabolicPathwayActivity: >>> [END]" in output
    assert "MetabolicFeatures: <<< ['MetabolicInput']" in output
    assert "MetabolicFeatures: >>> [END]" in output
    assert "MetabolicPathwayHeterogeneity: <<< ['MetabolicInput']" in output
    assert "MetabolicPathwayHeterogeneity: >>> [END]" in output


@pytest.mark.forked
def test_route_no_input(tmp_path):
    """Test route when neither LoadingRNAFromSeurat nor SampleInfo is provided."""
    config = tmp_path / "route_no_input.config.toml"
    config.write_text(
        f"""
        workdir = "{tmp_path}/workdir"
        outdir = "{tmp_path}/outdir"
    """
    )
    with pytest.raises(RuntimeError):
        output = dry_run(f"@{config}")

        assert "SampleInfo: <<< [START]" in output
        assert (
            "SampleInfo: This is a start process, "
            "but no 'input_data' specified." in output
        )
