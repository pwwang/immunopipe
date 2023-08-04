# ImmunarchLoading

Load the raw data into [`immunarch`][1] object, using [`immunarch::repLoad()`][2].

For the data path specified at `TCRData` in the input file, we will first find `filtered_contig_annotations.csv` and `filtered_config_annotations.csv.gz` in the path. If neighter of them exists, we will find `all_contig_annotations.csv` and `all_contig_annotations.csv.gz` in the path and a warning will be raised (You can find it at `./.pipen/<pipeline-name>/ImmunarchLoading/0/job.stderr`).

If none of the files exists, an error will be raised.

## Environment variables

- `prefix`: The prefix to the barcodes. You can use placeholder like `{Sample}_` to use the meta data from the `immunarch` object.

    !!! note
        This option is useful because the barcodes for the cells from scRNA-seq data are usually prefixed with the sample name, for example, `Sample1_AAACCTGAGAAGGCTA-1`. However, the barcodes for the cells from scTCR-seq data are usually not prefixed with the sample name, for example, `AAACCTGAGAAGGCTA-1`. So we need to add the prefix to the barcodes for the scTCR-seq data, and it is easier for us to integrate the data from different sources later.

[1]: https://immunarch.com
[2]: https://immunarch.com/reference/repLoad.html
