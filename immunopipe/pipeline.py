"""Assemable the pipeline"""
from pipen import Pipen
from .version import __version__
# Import your processeses
from .processes import starts

# Set the name and description of your pipeline here

# Make sure plugins and configurations get setup
Pipen.SETUP = False

pipeline = Pipen(
    name="immunopipe",
    desc=(
        "A pipeline for integrative analysis for scTCR- and scRNA-seq data\n"
        f"version: {__version__}"
    )
)

def main():
    """Run the pipeline"""
    pipeline.set_starts(starts).run()
