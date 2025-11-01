from __future__ import annotations

import re
import sys
from typing import Any, Dict

from yunpath import AnyPath
from pipen.utils import logger

WARNINGS = []


def _get_arg_from_cli(
    argname: str,
    default: Any | None = None,
    is_flag: bool = False,
    args: list[str] | None = None,
) -> Any:
    """Get argument from command line.

    Args:
        argname: The argument name.
        default: The default value if the argument is not found.
        is_flag: Whether the argument is a flag.
        args: The argument list to search. If None, use sys.argv.

    Returns:
        The argument value.
    """
    if args is None:
        args = sys.argv

    rex_eq = re.compile(f"--{re.escape(argname)}=(.+)")
    match_eq = [re.match(rex_eq, item) for item in args]

    rex_sp = re.compile(f"--{re.escape(argname)}")
    match_sp = [re.match(rex_sp, item) for item in args]

    for m1, m2 in zip(match_eq, match_sp):
        if not m1 and not m2:
            continue
        if m1:
            index = match_eq.index(m1)
            value = args[index].split("=", 1)[1]
            if is_flag:
                value = value.lower() in ("1", "true", "yes", "on")
            return value
        elif is_flag:
            return True
        else:
            index = match_sp.index(m2)

            if index + 1 < len(args):
                value = args[index + 1]
                return value

    return default


def _log_error(message: str | None = None) -> None:
    """Log an error message and exit.

    Args:
        message: The error message.
    """
    logger.error("Miscofigurations detected:")
    for warning in WARNINGS:
        logger.warning(f"- {warning}")

    if message:
        for line in message.splitlines():
            logger.error(f"- {line}")
        logger.error("")
        sys.exit(1)

    logger.warning("")


def validate_config(args: list[str] | None = None) -> Dict[str, Any]:
    """Validate the configuration.

    Args:
        config: The configuration.
    """
    if args is None:
        args = sys.argv

    try:
        from pipen_args import config
    except Exception as e:
        _log_error(f"Failed to load configuration.\n{e}")

    config.has_vdj = True  # Default to True, will be updated later
    if len(args) > 1 and args[1] == "gbatch":
        # Let immunopipe in the VM handle the validation
        # But we need to enable essential processes if we do have VDJ data
        config.has_vdj = (
            ("SampleInfo" in config and "LoadingRNAFromSeurat" not in config)
            or "ScRepLoading" in config
            or "TOrBCellSelection" in config
            or "ScRepCombiningExpression" in config
            or "TCRClustering" in config
            or "TESSA" in config
            or "ClonalStats" in config
            or "CDR3AAPhyschem" in config
        )
        return config

    if "TOrBCellSelection" not in config and "SeuratClusteringOfAllCells" in config:
        _log_error(
            "All cells are T cells ([TOrBCellSelection] is not set), "
            "so [SeuratClusteringOfAllCells] should not be used, "
            "use [SeuratClustering] instead."
        )

    if "TOrBCellSelection" not in config and "ClusterMarkersOfAllCells" in config:
        WARNINGS.append(
            "All cells are T cells ([TOrBCellSelection] is not set), "
            "so [ClusterMarkersOfAllCells] should not be used and will be ignored."
        )

    if "TOrBCellSelection" not in config and "TopExpressingGenesOfAllCells" in config:
        WARNINGS.append(
            "All cells are T cells ([TOrBCellSelection] is not set), "
            "so [TopExpressingGenesOfAllCells] should not be used and will be ignored."
        )

    # Input from Seurat object
    if "LoadingRNAFromSeurat" in config:
        LoadingRNAFromSeurat_prepared = _get_arg_from_cli(
            "LoadingRNAFromSeurat.envs.prepared",
            config.get("LoadingRNAFromSeurat", {}).get("envs", {}).get("prepared"),
            is_flag=True,
        )
        if LoadingRNAFromSeurat_prepared is None:
            LoadingRNAFromSeurat_prepared = False

        LoadingRNAFromSeurat_clustered = _get_arg_from_cli(
            "LoadingRNAFromSeurat.envs.clustered",
            config.get("LoadingRNAFromSeurat", {}).get("envs", {}).get("clustered"),
            is_flag=True,
        )
        if LoadingRNAFromSeurat_clustered is None:
            LoadingRNAFromSeurat_clustered = False

        if LoadingRNAFromSeurat_clustered:
            LoadingRNAFromSeurat_prepared = True

        config.setdefault("LoadingRNAFromSeurat", {}).setdefault("envs", {})
        config.LoadingRNAFromSeurat.envs.prepared = LoadingRNAFromSeurat_prepared
        config.LoadingRNAFromSeurat.envs.clustered = LoadingRNAFromSeurat_clustered
        config.has_vdj = "SampleInfo" in config

    # Input from sample info file
    else:
        infiles = _get_arg_from_cli(
            "SampleInfo.in.infile",
            config.get("SampleInfo", {}).get("in", {}).get("infile"),
        )

        if not isinstance(infiles, list):
            infiles = [infiles]

        if len(infiles) > 1:
            _log_error(
                "More than one input file specified in configuration file "
                "[SampleInfo.in.infile]."
            )

        if len(infiles) == 1 and infiles[0] is not None:
            infile = AnyPath(infiles[0])
            if infile.is_file():
                header = infile.read_text().splitlines()[0]
                config.has_vdj = "TCRData" in header or "BCRData" in header
            else:
                mount = config.get("scheduler_opts", {}).get("mount", [])
                if not isinstance(mount, list):
                    mount = [mount]

                # Let's check if infile a mounted path
                # Say infile is /mnt/disks/datadir/data/samples.txt
                # and we have fast_mout
                # gs://bucket/path:/mnt/disks/datadir
                # Then we can restore the cloud path for it:
                # gs://bucket/path/data/samples.txt
                if mount:
                    for mount in mount:
                        p1, p2 = mount.rsplit(":", 1)
                        p2 = AnyPath(p2)
                        if not infile.is_relative_to(p2):
                            continue
                        p1 = p1.rstrip("/")
                        cloud_path = f"{p1}/{infile.relative_to(p2)}"
                        if AnyPath(cloud_path).is_file():
                            header = AnyPath(cloud_path).read_text().splitlines()[0]
                            config.has_vdj = "TCRData" in header or "BCRData" in header
                        else:
                            _log_error(
                                f"Input file {infile} does not exist, "
                                "and the restored cloud path does not exist either: "
                                f"{cloud_path}."
                            )
                        break
                    else:
                        _log_error(
                            f"Input file {infile} does not exist, "
                            "and no mount can restore it as a cloud path."
                        )
                else:
                    _log_error(
                        f"Input file {infile} does not exist, "
                        "and no mount is specified to restore it as a cloud path."
                    )

    if WARNINGS:
        _log_error()

    return config
