# ScRepLoading

Load the single cell TCR/BCR data into a `scRepertoire` compatible object

This process loads the single cell TCR/BCR data into a `scRepertoire`
(>= v2.0.8, < v2.3.2) compatible object. Later, `scRepertoire::combineExpression`
can be used to combine the expression data with the TCR/BCR data.<br />

For the data path specified at `TCRData`/`BCRData` in the input file
(`in.metafile`), will be used to find the TCR/BCR data files and
`scRepertoire::loadContigs()` will be used to load the data.<br />

A directory can be specified in `TCRData`/`BCRData`, then
`scRepertoire::loadContigs()` will be used directly to load the data from the
directory. Otherwise if a file is specified, it will be symbolically linked to
a directory for `scRepertoire::loadContigs()` to load.<br />
Note that when the file name can not be recognized by `scRepertoire::loadContigs()`,
`envs.format` must be set for the correct format of the data.<br />

## Input

- `metafile`:
    The meta data of the samples
    A tab-delimited file
    Two columns are required:<br />
    * `Sample` to specify the sample names.<br />
    * `TCRData`/`BCRData` to assign the path of the data to the samples,
    and this column will be excluded as metadata.<br />

## Output

- `outfile`: *Default: `{{in.metafile | stem}}.scRep.qs`*. <br />
    The `scRepertoire` compatible object in qs/qs2 format

## Environment Variables

- `type` *(`choice`)*: *Default: `auto`*. <br />
    The type of the data to load.<br />
    - `TCR`:
        T cell receptor data
    - `BCR`:
        B cell receptor data
    - `auto`:
        Automatically detect the type from the metadata.<br />
        If `auto` is selected, the type will be determined by the presence of
        `TCRData` or `BCRData` columns in the metadata. If both columns are
        present, `TCR` will be selected by default.<br />
- `combineTCR` *(`type=json`)*: *Default: `{'samples': True}`*. <br />
    The extra arguments for `scRepertoire::combineTCR`
    function.<br />
    See also <https://www.borch.dev/uploads/screpertoire/reference/combinetcr>
- `combineBCR` *(`type=json`)*: *Default: `{'samples': True}`*. <br />
    The extra arguments for `scRepertoire::combineBCR`
    function.<br />
    See also <https://www.borch.dev/uploads/screpertoire/reference/combinebcr>
- `exclude` *(`auto`)*: *Default: `['BCRData', 'TCRData', 'RNAData']`*. <br />
    The columns to exclude from the metadata to add to the object.<br />
    A list of column names to exclude or a string with column names separated
    by `,`. By default, `BCRData`, `TCRData` and `RNAData` will be excluded.<br />
- `tmpdir`: *Default: `/tmp`*. <br />
    The temporary directory to store the symbolic links to the
    TCR/BCR data files.<br />
- `format` *(`choice`)*:
    The format of the TCR/BCR data files.<br />
    - `10X`:
        10X Genomics data, which is usually in a directory with
        `filtered_contig_annotations.csv` file.<br />
    - `AIRR`:
        AIRR format, which is usually in a file with
        `airr_rearrangement.tsv` file.<br />
    - `BD`:
        Becton Dickinson data, which is usually in a file with
        `Contigs_AIRR.tsv` file.<br />
    - `Dandelion`:
        Dandelion data, which is usually in a file with
        `all_contig_dandelion.tsv` file.<br />
    - `Immcantation`:
        Immcantation data, which is usually in a file with
        `data.tsv` file.<br />
    - `JSON`:
        JSON format, which is usually in a file with `.json` extension.<br />
    - `ParseBio`:
        ParseBio data, which is usually in a file with
        `barcode_report.tsv` file.<br />
    - `MiXCR`:
        MiXCR data, which is usually in a file with `clones.tsv` file.<br />
    - `Omniscope`:
        Omniscope data, which is usually in a file with `.csv`
        extension.<br />
    - `TRUST4`:
        TRUST4 data, which is usually in a file with
        `barcode_report.tsv` file.<br />
    - `WAT3R`:
        WAT3R data, which is usually in a file with
        `barcode_results.csv` file.<br />
        See also: <https://rdrr.io/github/ncborcherding/scRepertoire/man/loadContigs.html>
        If not provided, the format will be guessed from the file name by `scRepertoire::loadContigs()`.<br />

