# SeuratMetadataMutater

This process is used to integrate scTCR-seq data into the `Seurat` object. The scTCR-seq data is loaded by [ImmunarchLoading](./ImmunarchLoading.md) process. The integration is done by matching the barcodes from the `Seurat` object and the scTCR-seq data. The barcodes from the scTCR-seq data are prefixed with the sample name, for example, `Sample1_AAACCTGAGAAGGCTA-1`. The prefix is specified by the `prefix` environment variable in the [ImmunarchLoading](./ImmunarchLoading.md) process.

[ImmunarchLoading](./ImmunarchLoading.md) process will generate a text file with the information for each cell. `ImmunarchLoading.envs.metacols` can be used to specify the columns to be exported to the text file, which will then be integrated into the `Seurat` object by this process.

There is no environment variables for this process.
