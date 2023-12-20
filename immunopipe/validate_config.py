import sys
from typing import Any, Dict
from pipen.utils import logger


def validate_config(config: Dict[str, Any]) -> Dict[str, Any]:
    """Validate the configuration.

    Args:
        config: The configuration.
    """
    WARNINGS = []
    ERRORS = []

    if "TCellSelection" not in config and "SeuratClusteringOfAllCells" in config:
        ERRORS.append(
            "All cells are T cells ([TCellSelection] is not set), "
            "so [SeuratClusteringOfAllCells] should not be used."
        )

    if "TCellSelection" not in config and "ClusterMarkersOfAllCells" in config:
        WARNINGS.append(
            "All cells are T cells ([TCellSelection] is not set), "
            "so [ClusterMarkersOfAllCells] should not be used and will be ignored."
        )

    if "TCellSelection" not in config and "TopExpressingGenesOfAllCells" in config:
        WARNINGS.append(
            "All cells are T cells ([TCellSelection] is not set), "
            "so [TopExpressingGenesOfAllCells] should not be used and will be ignored."
        )

    if "SeuratMap2Ref" in config and "SeuratClustering" in config:
        ERRORS.append(
            "Cannot do both supervised [SeuratMap2Ref] and "
            "unsupervised [SeuratClustering] clustering."
        )

    if "SeuratMap2Ref" in config and "CellTypeAnnotation" in config:
        WARNINGS.append("[CellTypeAnnotation] is ignored when [SeuratMap2Ref] is used.")

    if WARNINGS or ERRORS:
        logger.warning("Miscofigurations detected:")
        for warning in WARNINGS:
            logger.warning(f"- {warning}")
        for error in ERRORS:
            logger.error(f"- {error}")
        logger.warning("")

    if ERRORS:
        sys.exit(1)

    return config
