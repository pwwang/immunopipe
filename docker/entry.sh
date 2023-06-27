#!/usr/bin/env bash

# Check if we have any arguments
if [ $# -eq 0 ]; then
    echo "There are two ways to run the pipeline:"
    echo "1. Run the pipeline from the command line"
    echo "   docker run <options> justold/immunopipe:master @<configfile> [options]"
    echo "   singularity run <options> docker://justold/immunopipe:master @<configfile> [options]"
    echo "2. Run the pipeline using pipen-board:"
    echo "   docker run <options> justold/immunopipe:master board immunopipe:Immunopipe -a /immunopipe/board.toml [options]"
    echo "   singularity run <options> docker://biopipen/scrna-basic:master board immunopipe:Immunopipe -a /workdir/board.toml [options]"
    exit 1
else if [ "$1" == "board" ]; then
    # Run the pipeline using pipen-board
    pipen $@
else
    # Run the pipeline from the command line
    immunopipe $@
fi
