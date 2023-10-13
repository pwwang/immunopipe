#!/usr/bin/env bash

# Check if we have any arguments
if [ $# -eq 0 ]; then
    Color_Off='\033[0m'       # Text Reset
    UPurple='\033[4;35m'      # Purple
    # Background
    On_Cyan='\033[46m'        # Cyan
    # Bold High Intensity
    BIGreen='\033[1;92m'      # Green
    BIBlue='\033[1;94m'       # Blue
    # High Intensity backgrounds
    On_IWhite='\033[0;107m'   # White

    echo -e "${On_IWhite}\n"
    echo "  There are two options for you to run the pipeline:"
    echo ""
    echo -e "  ${BIBlue}1. Running the pipeline from the command line directly${Color_Off}${On_IWhite}"
    echo ""
    echo -e "     ${On_Cyan}    Docker   ${Color_Off}${On_IWhite}${BIGreen} $ docker run \\"
    echo -e "                     --rm -w /workdir -v .:/workdir \\"
    echo "                     justold/immunopipe:<tag> \\"
    echo "                     @config.toml"
    echo -e "     ${On_Cyan} Singularity ${Color_Off}${On_IWhite}${BIGreen} $ singularity run \\"
    echo -e "                     --pwd /workdir -B .:/workdir -c -e -w \\"
    echo "                     docker://justold/immunopipe:<tag> \\"
    echo "                     @config.toml"
    echo -e "     ${On_Cyan}  Apptainer  ${Color_Off}${On_IWhite}${BIGreen} $ apptainer run \\"
    echo -e "                     --pwd /workdir -B .:/workdir -c -e -w --unsquash \\"
    echo "                     docker://justold/immunopipe:<tag> \\"
    echo "                     @config.toml"
    echo ""
    echo -e "  ${BIBlue}2. Run the pipeline using ${UPurple}pipen-board${Color_Off}${On_IWhite}"
    echo ""
    echo -e "     ${On_Cyan}    Docker   ${Color_Off}${On_IWhite}${BIGreen} $ docker run -p 18521:18521 \\"
    echo -e "                     --rm -w /workdir -v .:/workdir \\"
    echo "                     justold/immunopipe:<tag> board \\"
    echo "                     immunopipe:Immunopipe -a /immunopipe/board.toml"
    echo -e "     ${On_Cyan} Singularity ${Color_Off}${On_IWhite}${BIGreen} $ singularity run \\"
    echo -e "                     --pwd /workdir -B .:/workdir -c -e -w \\"
    echo "                     docker://justold/immunopipe:<tag> board \\"
    echo "                     immunopipe:Immunopipe -a /immunopipe/board.toml"
    echo -e "     ${On_Cyan}  Apptainer  ${Color_Off}${On_IWhite}${BIGreen} $ apptainer run \\"
    echo -e "                     --pwd /workdir -B .:/workdir -c -e -w --unsquash \\"
    echo "                     docker://justold/immunopipe:<tag> board \\"
    echo "                     immunopipe:Immunopipe -a /immunopipe/board.toml"
    echo ""
    echo -e "${Color_Off}"
    exit 1
elif [[ "$1" == "board" ]]; then
    # Run the pipeline using pipen-board
    # shellcheck disable=SC2068
    pipen $@
elif [[ "$1" == "@"* ]]; then
    # Run the pipeline using immunopipe
    # shellcheck disable=SC2068
    immunopipe $@
else
    # Run the pipeline from the command line
    "$@" || exit
fi
