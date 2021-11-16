"""Assemable the pipeline"""
from pipen import Pipen
# Import your processeses
from .processes import starts

# Set the name and description of your pipeline here
pipeline = Pipen(
    name="immunopipe",
    desc="A pipeline for integrative analysis for scTCR- and scRNA-seq data"
)

def main():
    """Run the pipeline"""
    pipeline.set_starts(starts).run()
