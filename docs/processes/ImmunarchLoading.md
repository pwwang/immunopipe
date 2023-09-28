# ImmunarchLoading

Load the raw data into [`immunarch`][1] object, using [`immunarch::repLoad()`][2].

For the data path specified at `TCRData` in the input file, we will first find `filtered_contig_annotations.csv` and `filtered_config_annotations.csv.gz` in the path. If neighter of them exists, we will find `all_contig_annotations.csv` and `all_contig_annotations.csv.gz` in the path and a warning will be raised (You can find it at `./.pipen/<pipeline-name>/ImmunarchLoading/0/job.stderr`).

If none of the files exists, an error will be raised.

This process will also generate a text file with the information for each cell. The file will be saved at `./.pipen/<pipeline-name>/ImmunarchLoading/0/output/<prefix>.tcr.txt`. The file will be used by the [SeuratMetadataMutater](./SeuratMetadataMutater.md) process to integrate the TCR-seq data into the `Seurat` object for further integrative analysis. `envs.metacols` can be used to specify the columns to be exported to the text file.

## Environment variables

- `prefix`: The prefix to the barcodes. You can use placeholder like `{Sample}_` to use the meta data from the `immunarch` object.

    /// Note
    This option is useful because the barcodes for the cells from scRNA-seq data are usually prefixed with the sample name, for example, `Sample1_AAACCTGAGAAGGCTA-1`. However, the barcodes for the cells from scTCR-seq data are usually not prefixed with the sample name, for example, `AAACCTGAGAAGGCTA-1`. So we need to add the prefix to the barcodes for the scTCR-seq data, and it is easier for us to integrate the data from different sources later.
    ///

- `tmpdir`: The temporary directory to link all data files.
    `Immunarch` scans a directory to find the data files. If the data files are not in the same directory, we can link them to a temporary directory and pass the temporary directory to `Immunarch`. This option is useful when the data files are in different directories.
- `mode`: Either "single" for single chain data or "paired" for
      paired chain data. For `single`, only TRB chain will be kept
      at `immdata$data`, information for other chains will be
      saved at `immdata$tra` and `immdata$multi`.
- `metacols` (`list`): The columns to be exported to the text file. The default value is `["Clones", "Proportion", "CDR3.aa"]`. You can refer to the [immunarch documentation][3] for the full list of the columns.

[1]: https://immunarch.com
[2]: https://immunarch.com/reference/repLoad.html
[3]: https://immunarch.com/articles/v2_data.html#immunarch-data-format
