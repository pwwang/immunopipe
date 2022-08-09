Currently only 10X data are supported.

To start the pipeline, you need 3 files:

- A metadata file, delimited by `TAB`, with 3 required columns:
  - `Sample`: A unique id for each sample
  - `TCRDir`: The directory for single-cell TCR data for this sample. Specifically, it should contain `filtered_contig_annotations.csv` from cellranger
  - `RNADir`: The directory for single-cell RNA data for this sample. Specifically, it should be able to be read by `Seurat::Read_10X()`
  - Other columns are meta data for the samples

Other optional files:

- If `METABOLIC` is enabled, a gmt file with metabolic pathways (can be found [here][1]) and a GTF file with exons are required for the analysis. The GTF file can be obtained by filter the type with `exon` from http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.refGene.gtf.gz or http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/genes/hg38.refGene.gtf.gz


[1]: https://github.com/pwwang/biopipen/blob/master/tests/data/scrna_metabolic/KEGG_metabolism.gmt
