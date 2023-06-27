"""Assemable the pipeline"""
from pipen import Pipen
from .version import __version__

# Import your processeses
from .processes import STARTS


class Immunopipe(Pipen):
    """The pipeline class"""
    starts = STARTS
    desc="A pipeline for integrative analysis for scTCR- and scRNA-seq data"


def main():
    """Run the pipeline"""
    Immunopipe().run()
