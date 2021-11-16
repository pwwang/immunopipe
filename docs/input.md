Currently only 10X data are supported.

To start the pipeline, you need 3 files:

- A metadata file, delimited by `TAB`, with 3 required columns:
  - `Sample`: A unique id for each sample
  - `TCRDir`: The directory for single-cell TCR data for this sample. Specifically, it should contain `filtered_contig_annotations.csv` from cellranger
  - `RNADir`: The directory for single-cell RNA data for this sample. Specifically, it should be able to be read by `Seurat::Read_10X()`
  - Other columns are meta data for the samples

- Two gene list files with first column the gene name in features and second the name to show in the plots
  - One is for us to investigate the gene expressions (a pile of boxplots show expressions of these genes for certain groups of cells)
  - The other is for us to plot heatmaps for these genes for certain groups of genes.
