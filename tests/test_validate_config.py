import pytest
import re
import sys
from pathlib import Path
from subprocess import Popen, PIPE
from immunopipe.validate_config import WARNINGS, _get_arg_from_cli, _log_error

configs_dir = Path(__file__).parent / "running" / "configs"


def _run_with_argv(args):
    cmd = [sys.executable, "-m", "immunopipe"] + args
    process = Popen(cmd, stdout=PIPE, stderr=PIPE, text=True)
    stdout, _ = process.communicate()
    return stdout


def test_get_arg_from_cli_eq():
    args = ["script.py", f"@{configs_dir}/sampleinfo.config.toml", "--param=value"]
    result = _get_arg_from_cli("param", default="default", args=args)
    assert result == "value"


def test_get_arg_from_cli_eq_flag():
    args = ["script.py", f"@{configs_dir}/sampleinfo.config.toml", "--flag=true"]
    result = _get_arg_from_cli("flag", default=False, is_flag=True, args=args)
    assert result is True


def test_get_arg_from_cli_sp():
    args = ["script.py", f"@{configs_dir}/sampleinfo.config.toml", "--param", "value"]
    result = _get_arg_from_cli("param", default="default", args=args)
    assert result == "value"


def test_get_arg_from_cli_sp_flag():
    args = ["script.py", f"@{configs_dir}/sampleinfo.config.toml", "--flag"]
    result = _get_arg_from_cli("flag", default=False, is_flag=True, args=args)
    assert result is True


def test_log_error(caplog):
    WARNINGS.extend(["Warning 1", "Warning 2"])
    with pytest.raises(SystemExit):
        _log_error("This is an error message.")

    assert "ERROR" in caplog.text
    assert "Warning 1" in caplog.text
    assert "Warning 2" in caplog.text
    assert "This is an error message." in caplog.text


def test_log_error_no_message(caplog):
    WARNINGS.extend(["Warning A", "Warning B"])
    _log_error()

    assert "WARNING" in caplog.text
    assert "Warning A" in caplog.text
    assert "Warning B" in caplog.text


def test_validate_config_malformed_config():
    args = [f"@{configs_dir}/malformed.config.toml"]
    stdout = _run_with_argv(args)
    stdout = re.sub(r"[\s\n]+", " ", stdout)
    assert "Failed to load configuration." in stdout


def test_validate_config_gbatch():
    args = ["gbatch"]
    stdout = _run_with_argv(args)
    stdout = re.sub(r"[\s\n]+", " ", stdout)
    assert "options to run pipeline on Google Cloud Batch" in stdout


def test_validate_config_torbcellselection_seuratclusteringofallcells():
    args = [f"@{configs_dir}/torbcellselection_seuratclusteringofallcells.config.toml"]
    stdout = _run_with_argv(args)
    stdout = re.sub(r"[\s\n]+", " ", stdout)
    assert "All cells are T cells" in stdout


def test_validate_config_mixed():
    args = [f"@{configs_dir}/mixed.config.toml"]
    stdout = _run_with_argv(args)
    stdout = re.sub(r"[\s\n]+", " ", stdout)
    assert "Miscofigurations detected:" in stdout
    assert (
        "[ClusterMarkersOfAllCells] should not be used and will be ignored."
        in stdout
    )
    assert (
        "[TopExpressingGenesOfAllCells] should not be used and will be ignored."
        in stdout
    )
    assert (
        "Cannot do both supervised [SeuratMap2Ref]"
        in stdout
    )


def test_validate_config_celltypeannotation_seuratmap2ref():
    args = [
        f"@{configs_dir}/celltypeannotation_seuratmap2ref.config.toml"
    ]
    stdout = _run_with_argv(args)
    stdout = re.sub(r"[\s\n]+", " ", stdout)
    assert "Miscofigurations detected:" in stdout
    assert (
        "[CellTypeAnnotation] is ignored when [SeuratMap2Ref] is used."
        in stdout
    )


def test_validate_config_loadingrnafromseurat():
    args = [
        f"@{configs_dir}/loadingrnafromseurat.config.toml"
    ]
    # Runs normally
    stdout = _run_with_argv(args)
    stdout = re.sub(r"[\s\n]+", " ", stdout)
    assert "Integrative analysis for" in stdout


def test_validate_config_loadingrnafromseurat_clustered():
    args = [
        f"@{configs_dir}/loadingrnafromseurat.config.toml",
        "--LoadingRNAFromSeurat.envs.clustered",
    ]
    # Runs normally
    stdout = _run_with_argv(args)
    stdout = re.sub(r"[\s\n]+", " ", stdout)
    assert "Integrative analysis for" in stdout


def test_validate_config_sampleinfo_multiple_infiles():
    args = [
        f"@{configs_dir}/sampleinfo_multiple_infiles.config.toml"
    ]

    stdout = _run_with_argv(args)
    stdout = re.sub(r"[\s\n]+", " ", stdout)
    assert (
        "More than one input file specified in configuration file"
        in stdout
    )
