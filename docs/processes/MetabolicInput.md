# MetabolicInput

This process takes Seurat object as input and pass it to the next processes in the `ScrnaMetabolicLandscape` group.

There is no configuration for this process.<br />

## Input

- `infile`:
    The input file

## Output

- `outfile`: *Default: `{{in.infile | basename}}`*. <br />
    The output symbolic link to the input file

## Environment Variables

- `copy` *(`flag`)*: *Default: `False`*. <br />
    Whether to copy the input file to the output file instead of
    creating a symbolic link.<br />

