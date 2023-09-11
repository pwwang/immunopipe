#!/usr/bin/env bash

# Check if we have any arguments
if [ $# -eq 0 ]; then
    echo "There are two ways to run the pipeline:"
    echo "1. Run the pipeline from the command line"
    echo "   $ docker run <options> justold/immunopipe:<tag> \\"
    echo "       @<configfile> [options]"
    echo "   $ singularity run <options> docker://justold/immunopipe:<tag> \\"
    echo "       @<configfile> [options]"
    echo "   Or:"
    echo "   $ docker run <options> justold/immunopipe:<tag> \\"
    echo "       immunopipe @<configfile> [options]"
    echo "   $ singularity run <options> docker://justold/immunopipe:<tag> \\"
    echo "       immunopipe @<configfile> [options]"
    echo ""
    echo "2. Run the pipeline using \`pipen-board\`:"
    echo "   $ docker run <options> justold/immunopipe:<tag> \\"
    echo "       board immunopipe:Immunopipe -a /immunopipe/board.toml [options]"
    echo "   $ singularity run <options> docker://biopipen/scrna-basic:<tag> \\"
    echo "       board immunopipe:Immunopipe -a /workdir/board.toml [options]"
    echo "   Or:"
    echo "   $ docker run <options> justold/immunopipe:<tag> \\"
    echo "       pipen board immunopipe:Immunopipe -a /immunopipe/board.toml [options]"
    echo "   $ singularity run <options> docker://biopipen/scrna-basic:<tag> \\"
    echo "       pipen board immunopipe:Immunopipe -a /workdir/board.toml [options]"
    exit 1
elif [[ "$1" == "board" ]]; then
    # Run the pipeline using pipen-board
    pipen $@
elif [[ "$1" == "@"* ]]; then
    # Run the pipeline using immunopipe
    immunopipe $@
else
    # Run the pipeline from the command line
    eval "$@"
fi
