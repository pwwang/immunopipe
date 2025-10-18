from __future__ import annotations

import re
import sys
from typing import Any, Dict

from diot import Diot
from yunpath import AnyPath
from pipen.utils import logger, update_dict

WARNINGS = []


def _get_arg_from_cli(
    argname: str,
    ns: bool = False,
    default: Any | None = None,
) -> Any:
    """Get argument first from configration file then from command line.

    Args:
        config: The configuration from the configuration file.
        argname: The argument name.
        ns: Whether to match namespace arguments.
            If True, will find `--argname.something` for `argname`.
        default: The default value if the argument is not found.

    Returns:
        The argument value.
    """
    rex_eq = re.compile(
        f"--{re.escape(argname)}\\.(.+)=(.+)"
        if ns
        else f"--{re.escape(argname)}=(.+)"
    )
    match_eq = [re.match(rex_eq, item) for item in sys.argv]

    rex_sp = re.compile(
        f"--{re.escape(argname)}\\.(.+)"
        if ns
        else f"--{re.escape(argname)}"
    )
    match_sp = [re.match(rex_sp, item) for item in sys.argv]

    if not ns:
        if any(match_eq):
            index = match_eq.index(True)
            value = sys.argv[index].split("=", 1)[1]
            return value
        elif any(match_sp):
            index = match_sp.index(True)
            if index + 1 < len(sys.argv):
                value = sys.argv[index + 1]
                return value

        return default

    out = Diot()
    for i, match in enumerate(match_eq):
        if match:
            key = match.group(1)
            value = match.group(2)
            out[key] = value
    for i, match in enumerate(match_sp):
        if match:
            key = match.group(1)
            if i + 1 < len(sys.argv):
                value = sys.argv[i + 1]
                out[key] = value

    return out or default


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


def validate_config() -> Dict[str, Any]:
    """Validate the configuration.

    Args:
        config: The configuration.
    """

    try:
        from pipen_args import config
    except Exception as e:
        config = Diot()

    if len(sys.argv) > 1 and sys.argv[1] == "gbatch":
        # Let immunopipe in the VM handle the validation
        return config

    torbcellselection = _get_arg_from_cli(
        "TOrBCellSelection", True, config.get("TOrBCellSelection")
    )
    seuratclusteringofallcells = _get_arg_from_cli(
        "SeuratClusteringOfAllCells", True, config.get("SeuratClusteringOfAllCells")
    )
    if torbcellselection is not None and seuratclusteringofallcells is not None:
        _log_error(
            "All cells are T cells ([TOrBCellSelection] is not set), "
            "so [SeuratClusteringOfAllCells] should not be used, "
            "use [SeuratClustering] instead."
        )

    clustermarkersofallcells = _get_arg_from_cli(
        "ClusterMarkersOfAllCells", True,  config.get("ClusterMarkersOfAllCells")
    )
    if torbcellselection is not None and clustermarkersofallcells is not None:
        WARNINGS.append(
            "All cells are T cells ([TOrBCellSelection] is not set), "
            "so [ClusterMarkersOfAllCells] should not be used and will be ignored."
        )

    topexpressinggenesofallcells = _get_arg_from_cli(
        "TopExpressingGenesOfAllCells", True, config.get("TopExpressingGenesOfAllCells")
    )
    seuratclustering = _get_arg_from_cli(
        "SeuratClustering", True, config.get("SeuratClustering")
    )
    if torbcellselection is not None and topexpressinggenesofallcells is not None:
        WARNINGS.append(
            "All cells are T cells ([TOrBCellSelection] is not set), "
            "so [TopExpressingGenesOfAllCells] should not be used and will be ignored."
        )

    seuratmap2ref = _get_arg_from_cli(
        "SeuratMap2Ref", True, config.get("SeuratMap2Ref")
    )
    if seuratmap2ref is not None and seuratclustering is not None:
        _log_error(
            "Cannot do both supervised [SeuratMap2Ref] and "
            "unsupervised [SeuratClustering] clustering."
        )

    celltypeannotation = _get_arg_from_cli(
        "CellTypeAnnotation", True, config.get("CellTypeAnnotation")
    )
    if seuratmap2ref is not None and celltypeannotation is not None:
        WARNINGS.append(
            "[CellTypeAnnotation] is ignored when [SeuratMap2Ref] is used."
        )

    loadrnafromseurat = _get_arg_from_cli(
        "LoadRNAFromSeurat", True, config.get("LoadRNAFromSeurat")
    )
    sampleinfo = _get_arg_from_cli("SampleInfo", True, config.get("SampleInfo"))
    # Input from Seurat object
    if loadrnafromseurat is not None:
        loadrnafromseurat_prepared = _get_arg_from_cli(
            "LoadRNAFromSeurat.envs.prepared",
            False,
            config.get("LoadRNAFromSeurat", {}).get("envs", {}).get("prepared"),
        )
        if loadrnafromseurat_prepared is None:
            loadrnafromseurat_prepared = False

        loadrnafromseurat_clustered = _get_arg_from_cli(
            "LoadRNAFromSeurat.envs.clustered",
            False,
            config.get("LoadRNAFromSeurat", {}).get("envs", {}).get("clustered"),
        )
        if loadrnafromseurat_clustered is None:
            loadrnafromseurat_clustered = False

        if loadrnafromseurat_clustered:
            loadrnafromseurat_prepared = True

        config.setdefault("LoadRNAFromSeurat", {}).setdefault("envs", {})
        config.LoadRNAFromSeurat.envs.prepared = loadrnafromseurat_prepared
        config.LoadRNAFromSeurat.envs.clustered = loadrnafromseurat_clustered
        config.has_vdj = sampleinfo is not None

    # Input from sample info file
    else:
        infiles = _get_arg_from_cli(
            "SampleInfo.in.infile",
            False,
            config.get("SampleInfo", {}).get("in", {}).get("infile"),
        )
        if infiles is None:
            _log_error(
                "No input file specified in configuration file "
                "[SampleInfo.in.infile], assuming passing from CLI."
            )

        if not isinstance(infiles, list):
            infiles = [infiles]

        if len(infiles) > 1:
            _log_error(
                "More than one input file specified in configuration file "
                "[SampleInfo.in.infile]."
            )

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
