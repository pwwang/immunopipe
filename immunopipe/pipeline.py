"""Assemable the pipeline"""
from pipen import Pipen
from pipen_args import parser

from .version import __version__, desc  # noqa: F401

# Import your processeses
from .processes import start_processes


class Immunopipe(Pipen):
    """The pipeline class"""
    starts = start_processes
    desc = desc


def main():
    """Run the pipeline"""
    parser.description = desc
    Immunopipe().run()
