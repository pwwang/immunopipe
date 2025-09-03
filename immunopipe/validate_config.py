import sys
from typing import Any, Dict
from yunpath import AnyPath
from pipen.utils import logger


def validate_config() -> Dict[str, Any]:
    """Validate the configuration.

    Args:
        config: The configuration.
    """
    WARNINGS = []
    ERRORS = []

    try:
        from pipen_args import config
    except Exception as e:
        ERRORS.append(repr(e))

    if not ERRORS and not config:  # no arguments
        config.has_vdj = True
        return config

    if not ERRORS and (
        "TOrBCellSelection" not in config and "SeuratClusteringOfAllCells" in config
    ):
        ERRORS.append(
            "All cells are T cells ([TOrBCellSelection] is not set), "
            "so [SeuratClusteringOfAllCells] should not be used, "
            "use [SeuratClustering] instead."
        )

    if not ERRORS and (
        "TOrBCellSelection" not in config and "ClusterMarkersOfAllCells" in config
    ):
        WARNINGS.append(
            "All cells are T cells ([TOrBCellSelection] is not set), "
            "so [ClusterMarkersOfAllCells] should not be used and will be ignored."
        )

    if not ERRORS and (
        "TOrBCellSelection" not in config and "TopExpressingGenesOfAllCells" in config
    ):
        WARNINGS.append(
            "All cells are T cells ([TOrBCellSelection] is not set), "
            "so [TopExpressingGenesOfAllCells] should not be used and will be ignored."
        )

    if not ERRORS and "SeuratMap2Ref" in config and "SeuratClustering" in config:
        ERRORS.append(
            "Cannot do both supervised [SeuratMap2Ref] and "
            "unsupervised [SeuratClustering] clustering."
        )

    if not ERRORS and "SeuratMap2Ref" in config and "CellTypeAnnotation" in config:
        WARNINGS.append(
            "[CellTypeAnnotation] is ignored when [SeuratMap2Ref] is used."
        )

    if not ERRORS:
        if "LoadRNAFromSeurat" in config:
            config.LoadRNAFromSeurat.setdefault("envs", {})
            config.LoadRNAFromSeurat.envs.setdefault("prepared", False)
            config.LoadRNAFromSeurat.envs.setdefault("clustered", False)
            if config.LoadRNAFromSeurat.envs.clustered:
                config.LoadRNAFromSeurat.envs.prepared = True

            config.has_vdj = "SampleInfo" in config

        else:
            infiles = config.get("SampleInfo", {}).get("in", {}).get("infile", [])
            if not isinstance(infiles, list):
                infiles = [infiles]

            if not infiles:
                WARNINGS.append(
                    "No input file specified in configuration file "
                    "[SampleInfo.in.infile], assuming passing from CLI."
                )
                WARNINGS.append("Assuming scTCR-seq/scBCR-seq data is present")

            elif len(infiles) > 1:
                WARNINGS.append(
                    "More than one input file specified in configuration file "
                    "[SampleInfo.in.infile], only the first one will be used."
                )
                config["SampleInfo"]["in"]["infile"] = [infiles[0]]

            infile = AnyPath(infiles[0])
            if infile.is_file():
                header = infile.read_text().splitlines()[0]
                config.has_vdj = "TCRData" in header or "BCRData" in header
            else:
                mount = (
                    config.get("scheduler_opts", {}).get("mount", [])
                    or config.get("cli-gbatch", {}).get("mount", [])
                )
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
                            ERRORS.append(
                                f"Input file {infile} does not exist, "
                                "and the restored cloud path does not exist either: "
                                f"{cloud_path}."
                            )
                        break
                    else:
                        ERRORS.append(
                            f"Input file {infile} does not exist, "
                            "and no mount can restore it as a cloud path."
                        )
                else:
                    ERRORS.append(
                        f"Input file {infile} does not exist, "
                        "and no mount is specified to restore it as a cloud path."
                    )

    if ERRORS:
        logger.error("Miscofigurations detected:")
        for error in ERRORS:
            logger.error(f"- {error}")
        logger.error("")
        sys.exit(1)

    if WARNINGS:
        logger.warning("Miscofigurations detected:")
        for warning in WARNINGS:
            logger.warning(f"- {warning}")
        logger.warning("")

    return config
