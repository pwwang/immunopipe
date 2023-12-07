"""Assemable the pipeline"""
from pipen import Pipen
from pipen_args import parser

from .version import __version__  # noqa: F401

# Import your processeses
from .processes import SampleInfo


class Immunopipe(Pipen):
    """The pipeline class"""
    starts = [SampleInfo]
    desc = (
        f"Immunopipe (v{__version__}): "
        "Integrative analysis for scRNA-seq and scTCR-seq data"
    )


def main():
    """Run the pipeline"""
    parser.description = Immunopipe.desc
    Immunopipe().run()
