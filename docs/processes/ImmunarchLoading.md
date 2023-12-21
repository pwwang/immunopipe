# ImmunarchLoading

Immuarch - Loading data

Load the raw data into [`immunarch`](https://immunarch.com) object,
using [`immunarch::repLoad()`](https://immunarch.com/reference/repLoad.html).<br />

For the data path specified at `TCRData` in the input file, we will first find
`filtered_contig_annotations.csv` and `filtered_config_annotations.csv.gz` in the
path. If neighter of them exists, we will find `all_contig_annotations.csv` and
`all_contig_annotations.csv.gz` in the path and a warning will be raised
(You can find it at `./.pipen/<pipeline-name>/ImmunarchLoading/0/job.stderr`).<br />

If none of the files exists, an error will be raised.<br />

This process will also generate a text file with the information for each cell.<br />
The file will be saved at
`./.pipen/<pipeline-name>/ImmunarchLoading/0/output/<prefix>.tcr.txt`.<br />
The file can be used by the [`IntegratingTCR`](./IntegratingTCR.md) process to integrate the
TCR-seq data into the `Seurat` object for further integrative analysis.<br />
`envs.metacols` can be used to specify the columns to be exported to the text file.<br />

## Environment Variables

- `prefix`: *Default: `{Sample}_`*. <br />
    The prefix to the barcodes. You can use placeholder like `{Sample}_`
    to use the meta data from the `immunarch` object. The prefixed barcodes will
    be saved in `out.metatxt`. The `immunarch` object keeps the original barcodes, but
    the prefix is saved at `immdata$prefix`.<br />

    /// Note
    This option is useful because the barcodes for the cells from scRNA-seq
    data are usually prefixed with the sample name, for example,
    `Sample1_AAACCTGAGAAGGCTA-1`. However, the barcodes for the cells from
    scTCR-seq data are usually not prefixed with the sample name, for example,
    `AAACCTGAGAAGGCTA-1`. So we need to add the prefix to the barcodes for
    the scTCR-seq data, and it is easier for us to integrate the data from
    different sources later.<br />
    ///

- `tmpdir`: *Default: `/tmp/user`*. <br />
    The temporary directory to link all data files.<br />
    `Immunarch` scans a directory to find the data files. If the data files
    are not in the same directory, we can link them to a temporary directory
    and pass the temporary directory to `Immunarch`.<br />
    This option is useful when the data files are in different directories.<br />
- `mode`: *Default: `single`*. <br />
    Either "single" for single chain data or "paired" for
    paired chain data. For `single`, only TRB chain will be kept
    at `immdata$data`, information for other chains will be
    saved at `immdata$tra` and `immdata$multi`.<br />
- `extracols` *(`list`)*: *Default: `[]`*. <br />
    The extra columns to be exported to the text file.<br />
    You can refer to the
    [immunarch documentation](https://immunarch.com/articles/v2_data.html#immunarch-data-format)
    to get a sense for the full list of the columns.<br />
    The columns may vary depending on the data source.<br />
    The columns from `immdata$meta` and some core columns, including
    `Barcode`, `CDR3.aa`, `Clones`, `Proportion`, `V.name`, `J.name`, and
    `D.name` will be exported by default. You can use this option to
    specify the extra columns to be exported.<br />

