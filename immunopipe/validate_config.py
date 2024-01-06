import sys
from typing import Any, Dict
from pathlib import Path
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
        config.has_tcr = True
        return config

    if not ERRORS and (
        "TCellSelection" not in config and "SeuratClusteringOfAllCells" in config
    ):
        ERRORS.append(
            "All cells are T cells ([TCellSelection] is not set), "
            "so [SeuratClusteringOfAllCells] should not be used, "
            "use [SeuratClustering] instead."
        )

    if not ERRORS and (
        "TCellSelection" not in config and "ClusterMarkersOfAllCells" in config
    ):
        WARNINGS.append(
            "All cells are T cells ([TCellSelection] is not set), "
            "so [ClusterMarkersOfAllCells] should not be used and will be ignored."
        )

    if not ERRORS and (
        "TCellSelection" not in config and "TopExpressingGenesOfAllCells" in config
    ):
        WARNINGS.append(
            "All cells are T cells ([TCellSelection] is not set), "
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
        config.has_tcr = True
        infiles = config.get("SampleInfo", {}).get("in", {}).get("infile", [])
        if not isinstance(infiles, list):
            infiles = [infiles]
        if not infiles:
            WARNINGS.append(
                "No input file specified in configuration file [SampleInfo.in.infile], "
                "assuming passing from CLI."
            )
            WARNINGS.append("Assuming scTCR-seq data is present")
        else:
            if len(infiles) > 1:
                WARNINGS.append(
                    "More than one input file specified in configuration file "
                    "[SampleInfo.in.infile], only the first one will be used."
                )
                config["SampleInfo"]["in"]["infile"] = [infiles[0]]

            infile = Path(infiles[0])
            if not infile.is_file():
                ERRORS.append(f"Input file {infile} does not exist.")
            else:
                header = infile.read_text().splitlines()[0]
                config.has_tcr = "TCRData" in header

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
