# Change Log

## 2.3.1

- chore: update r-plotthis dependency version to 0.9.4
- chore: update r-scplotter dependency version to 0.6.5
  - feat: clone selectors gain a `within` option to select clones within a specific subset
  - feat: clone selectors `and`/`or` support multiple selectors (more than 2)
- chore: update biopipen dependency version to 1.1.4
  - fix(scrna.MarkersFinder): update database handling in enrichment plots in case it is from gmt files
  - fix(scrna.PseudoBulkDEG): update database handling in enrichment plots to use unique databases in case gmt files were used

## 2.3.0

- chore: bump pipen to v1.1 so asynchronous path operations are supported
- refactor: convert synchronous methods to asynchronous in gbatch.py
- feat: add example configuration file
- fix: update header reading method in validate_config to use read_text
- docs: update process documentation with additional references and image links
- docs: update installation and running documentation to reflect dependency changes
- docs: add FAQ entry for mouse data support in the pipeline
- chore: update r-biopipen.utils version to 0.3.5 in environment_rpkgs.yml
  - fix(RunSeuratDEAnalysis): handle +Inf/-Inf values in avg_log2FC to make visualization works
- chore: bump biopipen to 1.1.1
  - docs(tcr.ClonalStats): correct typo in documentation for ClonalStatsPlot reference

## 2.2.2

- chore: force conda to rebuild with biopipen.utils.R v0.3.4-7
  - make sure 'direct' method of CellTypeAnnotation keeps the order
- chore: update r-plotthis dependency version to 0.9.0=r43_3
- chore: update biopipen dependency to version 1.0.0
  - fix(scrna.CellTypeAnnotation): handle NA values in cell type processing to prevent errors
  - chore(scrna.SeuratPreparing): update min_cells parameter to require at least 3 cells for SCTransform
- fix: enhance seurat path resolution and error handling in check_genes function
- fix: improve seurat path resolution for containerized environments in check_dim function

## 2.2.1

- docs: add FAQ entries for process reruns after updating and skipping essential processes
- docs: add installation note for pipen-dry package in FAQ
- chore: add pipen-dry as an optional dependency and include it in extras
- ci(docker): use multi-staging building and make the packages relying numpy v1 to be merged to the base environment, which largely reduce the size of the image
- chore: bump pipen-cli-gbatch to 0.1.5 (make default running profile work)
- chore: update pipen-report version to 0.23.15
- chore: update clustcr version in environment_base.yml to be compatible with numpy v2
- chore: add .dockerignore file to exclude unnecessary files from Docker builds
- chore: update Python version to 3.12 and adjust biopipen.utils version in environment files
- chore: add Dockerhub descriptions for immunopipe, immunopipe-base, and immunopipe-rpkgs
- chore: update r-scplotter version to 0.6.4=r43_2 in environment_rpkgs.yml
  - fix(ClonalDiversityPlot): fix calculating Gini coefficient
  - fix(MarkersPlot): handle errors when subsetting object features
  - fix(MarkersPlot): exclude groups that failed DE analysis from plotting
- chore: update r-plotthis version to 0.9.0=r43_1 in environment_rpkgs.yml
  - fix(boxviolinplot): skip processing for data frames with less than 2 rows
  - fix(network): update ggplot2 version check for link_type_by support
  - fix(velocityplot): update ggplot2 version check for arrow length handling
  - chore(heatmap): add warning for unknown arguments in HeatmapAtomic function
  - feat: add list_fonts and import_font utilities
- chore: update r-biopipen.utils version in environment_rpkgs.yml to 0.3.4=r43_7
  - fix(RunSeuratDEAnalysis): improve error handling by returning empty data frame on logic error and adding traceback on stop
  - feat: set default font to "LiberationSans" for plotthis plots
  - fix: update RunUMAPArgs to use the correct reduction and improve dims handling
  - fix(RunSeuratUMAP): avoid changing of identity when using features for RunUMAP
  - chore(RunSeuratTransformation, RunSeuratMap2Ref): ensure return.only.var.genes defaults to FALSE when using SCTransform
- chore: update biopipen version to ^0.34.28 in pyproject.toml
  - feat(tcr.ClonalStats): add save_data parameter to ClonalStats for saving plot data
  - fix(scrna.PseudoBulkDEG): change default assay from "RNA" to None so that default assay can be used by default
  - feat(scrna.ScFGSEA): add assay parameter to allow specification of assay in analysis
  - fix(scrna.CellTypeAnnotation): correct assignment of identities in rename_idents function
  - feat(tcr.CDR3Clustering): add verbose output option for GIANA command
  - fix(scrna.CellTypeAnnotation): update package requirement from celltypist to celltypist2 (a version adopts numpy v2)
  - fix(tcr.TESSA): add Keras model migration support for v2 and v3
  - fix(scrna.CellTypeAnnotation): improve handling of over_clustering assignment from Seurat object for celltypist

## 2.2.0

- BREAKING: restructure Docker workflows and update environment configurations using multiple-image/stage building, no "-full" suffix for image tags needed anymore.
- BREAKING: rename TCRClustering to CDR3Clustering across documentation and codebase to adopt for BCR data
- BREAKING: prevent export of RDS/qs files in SeuratPreparing and SeuratMap2Ref classes.
  > NOTE that this may break reruning with cached SeuratPreparing and SeuratMap2Ref processes.
  > The RDS/qs output files are no longer exported to output directory to reduce the report size. Instead, they are staying in the working directory for downstream processes to use.
  > But the the figures will be still displayed in the report as before.
- fix: swap SeuratMap2Ref and SeuratClustering to make sure graph doesn't get overwritten
- fix: refactor VDJ configuration detection in validate_config when running with gbatch
- fix: default group_by to None for ClusterMarkers/ClusterMarkersOfAllCells so default ident will be used
- docs: update configuration options in getting started guide
- docs: update notes in Seurat processes to prevent overwriting cluster information
- docs: update column naming in SeuratClustering and CellTypeAnnotation documentation
- docs: update figures after CDR3Clustering renaming
- ci: add cache deletion step to save space in Docker workflow
- test: update configuration files and tests
- test: update assertions in test_route_sampleinfo_full for correct flow after refactoring
- test: reorder assertions in test_route_sampleinfo_full for correct flow due to the swapping of SeuratMap2Ref and SeuratClustering
- chore: bump biopipen version to 0.34.23 (related changes below)
  - docs(scrna.SeuratPreparing): improve documentation for SeuratPreparing to make sure new lines in code blocks
  - docs(scrna.SeuratClusterStats): enhance documentation for cluster statistics plots
  - docs(tcr.ClonalStats): add exapmle configurations and plots
  - docs(scrna.ScFGSEA): add examples and usage details for GSEA plots
  - docs(CellCellCommunication): add example output for ligand-receptor interactions
  - docs(CellCellCommunicationPlots): add examples for various plot types
  - docs(scrna.MetabolicPathwayActivity): add merged heatmap example
  - docs: improve documentation for SeuratPreparing to make sure new line… (#180)
  - fix(scrna): update group_by assignment to use GetIdentityColumn
  - feat(scrna.SeuratClustering): add envs.ident as shortcut for custom cluster name
  - chore(scrna.SeuratClusterStats): set default group_by for dimplots to None
  - docs(scrna): correct parameter names in docstring (ident-1 to ident_1, and ident-2 to ident-2)
  - feat(scrna): add return option for identity column in convert_seurat_to_anndata function
  - feat(scrna.CellCellCommunication): enhance groupby handling for Seurat objects
  - feat(scrna.CellTypeAnnotation): add support for specifying identity column in CellTypeAnnotation for celltypist
  - feat(scrna.CellTypeAnnotation): enhance identity column handling and add backup column support
  - feat(scrna.CellTypeAnnotation): enhance classifier initialization and input handling for celltypist
  - enh(scrna.MarkersFinder): add logging for plot processing in marker and enrichment functions
  - feat(scrna.CellTypeAnnotation): require celltypist package with specified version and python interpreter
  - docs(scrna): update dbs/gmtfile parameter description for enrichment analysis
  - chore(scrna): update future.globals.maxSize to Inf for improved memory handling
  - feat(tcr.CDR3Clustering): rename TCRClustering to CDR3Clustering and adopt BCR data
  - docs(scrna.CellTypeAnnotation): clarify renaming of original identity column during cell type annotation
  - see full changelog at <https://pwwang.github.io/biopipen/CHANGELOG/>.

## 2.1.3

- docs: update docs for ClusterMarkers with examples
- ci: update biopipen dependency for tests
- feat: add new tests and configuration for ScRepLoading, ScRepCombiningExpression and SeuratClusterStats
- test: add test for clustermarkers
- test: add test for TOrBCellSelection
- docs: enhance documentation for PseudoBulkDEG and MarkersFinder
- test: add ClonalStats tests and configuration
- test: add ScFGSEA tests and configuration
- test: add CellCellCommunication and CellCellCommunicationPlots tests
- docs: update ScrnaMetabolicLandscape doc images
- chore: bump biopipen to 0.34.17
  - fix(scrna.MarkersFinder): enhance error handling in enrichment plotting and fix all-enrich plots when ident_1 is NULL
  - fix(tcr.ClonalStats): replace deprecated ClonalDynamicsPlot with ClonalStatPlot

## 2.1.2

- test: update order of ScRepLoading and SeuratClustering in assertions in route tests
- docs: update configuration section name from cli-gbatch to gbatch
- docs: enhance installation and introduction documentation
- fix: better infer has_vdj when running by gbatch

## 2.1.1

- fix: correct order of input requirements for ScRepCombiningExpression

## 2.1.0

- BREAKING: rename process LoadRNAFromSeurat to LoadingRNAFromSeurat
- feat: enhance configuration validation logic
- feat: implement LoadRNAFromSeurat process to allow specifying sample column
- feat: allow SeuratClustering, SeuratMap2Ref and CellTypeAnnotation at the same time in a run
- feat: add "help" subcommand to show help for specific processes
- fix: make sure the routes of the pipeline work as expected (tests added)
- fix: ensure correct handling of kp_mount in main function
- docs: documentation and update processes for T/B cell selection and clustering
- ci: improve immunopipe installation check
- docs: update configurations and introduction for LoadingRNAFromSeurat process
- chore: bump xqute to 0.10.17
- chore: bump pipen to 0.17.24
- chore: bump pipen-cli-gbatch to v0.1.3
- chore: bump pipen-args to v0.17.7
- chore: bump pipen-report to 0.23.13
- chore: update biopipen to version 0.34.16
  - feat(scrna.MarkersFinder): allow using other metadata columns from object for enrichment plot of all subcases

## 2.0.3

- feat: adopt pipen-cli-gbatch, allowing --mount-as-cwd to infer workdir and outdir
- fix: fix when running with a configuration file solely
- fix: fix gbatch deamon workdir not created inside the pipeline working folder for gbatch cli
- fix: fix outdir not following the pipeline name with gbatch cli
- docs: update volume mounting syntax in docker and singularity commands
- docs: update `immunopipe gbatch` options for running the pipeline on Google Cloud Batch Jobs
- chore: bump pipen-cli-gbatch to 0.0.7 and pipen-report to 0.23.12
- chore: bump biopipen to 0.34.14
  - fix(scrna_metabolic_landscape): fix report paging issue
  - docs(scrna.MarkersFinder): fix links in docs
  - fix(scrna.SeuratClusterStats): improve error handling in feature plotting when save_code (due to upgrade to ggplot2 v4)
  - feat(MarkersFinder): use `scplotter::MarkersPlot` (wrapped by `biopipen.utils::VizDEGs` to visualize marker
  - fix(scrna.CellTypeAnnotation): update logging for celltypist command execution

## 2.0.2

- fix: fix immunopipe -h/--help not working as expected
- feat(gbatch): enhance validation and default handling in main function
- feat(gbatch): add --mount-as-cwd option to mount a cloudpath as working directory
- chore: bump biopipen to 0.34.10
  - docs(scrna.SeuratPreparing): enhance cell_qc parameter description in SeuratPreparing
  - docs(scrna.ModuleScoreCalculator): update link format in ModuleScoreCalculator docstring

## 2.0.1

- feat: add utility cli for gene checking and dimension verification
- chore: make tests order work as expected
- chore: bump biopipen to 0.34.9
  - fix(scrna.CellCellCommunication): handle numpy product attribute error
  - feat(scrna.ModuleScoreCalculator): add post mutaters functionality to allow compound modules based on added modules
  - docs(scrna.MarkersFinder): correct URL in documentation
  - feat(scrna.CellTypeAnnotation): add support for additional direct cell type annotations
  - feat(scrna.MarkersFinder): enhance enrichment plot descriptions
  - chore(scrna.CellCellCommunicationPlots): set default case to "Cell-Cell Communication"
  - feat(scrna.CellCellCommunicationPlots): add table output option for ccc data
- chore: bump scplotter to 0.6.0 in docker image
  - chore: optimize clonal_size_data
  - feat(featurestatplot): add pos_only parameter to filter positive feature values
  - chore(enrichmentplot) remove preset label_nudge parameter in BarPlot call
  - feat: add downsample parameter to feature statistic plots
  - feat(featurestatplot) default options to TRUE to show row and column names when plot_type is heatmap-alike

## 2.0.0 (Cloud support, enhanced visualization, new analyses/features, and more ...)

### Cloud support

- feat: immunopipe can now be run on Google Cloud Batch Jobs, allowing for scalable and efficient processing of larger datasets.
  - You can either run the pipeline using the [gbatch scheduler](https://pwwang.github.io/pipen/scheduler/#gbatch); or
  - run the entire pipeline on Gooble Batch Jobs using `immunopipe gbatch` command.
  - See [Run the pipeline using Google Cloud Batch Jobs¶](https://pwwang.github.io/immunopipe/latest/running/#run-the-pipeline-using-google-cloud-batch-jobs) for more details.

### Enhanced visualization

- feat: [scplotter](https://github.com/pwwang/scplotter) and [plotthis](https://github.com/pwwang/plotthis) are now used for plotting, providing enhanced visualization capabilities and uniformity across different processes.
- feat: default descriptions/captions are now added to plots, making them more informative.

### Enhanced performance

- feat: the pipeline now uses [`qs2`](https://github.com/qsbase/qs2) for store the R objects, which speeds up the loading and saving of Seurat objects.
- feat: step-wise caching (in addition to process-wise) is now supported, especially for Seurat processes, allowing for faster re-running of the pipeline by caching intermediate results and improving results reproducibility.
- refactor: [MetabolicLandscapeAnalysis](https://pwwang.github.io/immunopipe/dev/processes/ScrnaMetabolicLandscape/) is refactored for flexibility and performance improvement.
- BREAKING: [enrichR](https://cran.r-project.org/web/packages/enrichR/vignettes/enrichR.html) is retired and replaced by [enrichit](https://github.com/pwwang/enrichit) for enrichment analysis, making it offline and more flexible. This enables the entire pipline to run without internet connection.

### New analyses/features

- feat: the pipeline now supports [cell-cell communication analysis](https://pwwang.github.io/immunopipe/latest/processes/CellCellCommunication/).
- feat: plots are supported for all cases for [MarkersFinder](https://pwwang.github.io/immunopipe/latest/processes/MarkersFinder/) and [ScFGSEA](https://pwwang.github.io/immunopipe/latest/processes/ScFGSEA/), allowing plotting the markers (DEGs) and enriched pathways for all cases (e.g. all seurat clusters) in a single plot.
- BREAKING: the [immunarch](https://immunarch.com) package is now replaced by [scRepertoire](https://github.com/BorchLab/scRepertoire) for more features and allowing customized clonotype definition.
- feat: `envs.mutaters` is now supported for [SeuratPreparing](https://pwwang.github.io/immunopipe/latest/processes/SeuratPreparing/) to allow create factor (categorical) columns in the metadata.
- feat: [PseudoBulkDEG](https://pwwang.github.io/immunopipe/latest/processes/PseudoBulkDEG/) is added to perform pseudo-bulk differential expression analysis.
- feat: BCR-seq data is now supported, allowing users to analyze BCR-seq data paired with scRNA-seq data.
- feat: Now `Seurat` object (in RDS or qs2 format) is supported as input for scRNA-seq data.
- feat: Now loom format is supported for scRNA-seq data, allowing users to use loom files as input for the pipeline.
- feat: add mcp server functionality to launch mcp server to help compose the configuration file.

### House keeping

- build: docker images are now built based on the [`biopipen`](https://github.com/pwwang/biopipen) base image.
- ci: the test workflow now caches the running intermediate files to speed up the tests.
- docs: the citation information is now added to the documentation, allowing users to easily cite the pipeline in their publications.
- chore(deps): biopipen is bumped to 0.34.8, which includes various bug fixes and enhancements. See the [biopipen releases](https://github.com/pwwang/biopipen/releases) for more details.

## 1.4.4

- chore(deps): add gcc_linux-64 to Docker environment dependencies
- docs: add input and output sections to multiple process documentation files
- feat: add PDF output for K-means and T cell plots, enhancing report generation
- deps: bump biopipen to 0.32.3
  - 0.32.1:
    - fix(scrna.ScFGSEA): fix case gmtfile not working
    - fix(TopExpressingGenes): add InlineNotification component to TopExpressingGenes.svelte
    - fix(scrna.SeuratClusterStats): fix kind not being added to the figure file name for plots of features
    - feat(scrna.SeuratPreparing): support percent.mt, percent.ribo, percent.hb and percent.plat for mouse
  - 0.32.2:
    - feat: add PDF output option for SampleInfo plots
    - feat: add PDF output options for violin and scatter plots in Seurat preparation scripts
    - feat: add PDF output options for volcano, dotplot, venn, and upset plots
    - feat: add PDF output option for Enrichr plots in TopExpressingGenes script
    - feat: add PDF output options for UMAP plots in SeuratMap2Ref script; update image handling in misc.liq
    - feat: add PDF output options for cluster size distribution, shared clusters, and sample diversity plots; update plotting functions to handle multiple output formats
    - feat: add PDF output options for various Immunarch scripts; enhance reporting with downloadable PDF files
    - feat: add PDF output options for cluster size distribution, dimension plots, and feature plots; enhance reporting with downloadable PDF files
    - feat: add PDF output options for radar and bar plots; enhance reporting with downloadable PDF files
    - feat: add PDF output options for CloneResidency script; enhance reporting with downloadable PDF files
    - feat: add PDF output options for GSEA table and enrichment plots; enhance reporting with downloadable PDF files
    - feat: add PDF output options for pie charts, heatmaps, Venn plots, and UpSet plots; enhance reporting with downloadable PDF files
    - feat: add PDF output options for Enrichr plots; enhance reporting with downloadable PDF files
    - feat: add PDF output options for estimated coefficients and distribution plots; enhance reporting with downloadable PDF files
  - 0.32.3:
    - chore: add descriptive summaries for fgsea and enrichr results

## 1.4.3

- deps: update pipen-runinfo dependency to version 0.8.0 (pipen to 0.15.2)
- deps: update biopipen dependency to version 0.31.4
  - fix(scrna.SeuratMap2Ref): fix refnorm not detected for NormlizeData'ed reference

## 1.4.2

- deps: add bioconductor-destiny dependency in docker environment files for ModuleScoreCalculator

## 1.4.1

- docs: update tutorial dataset information and links

## 1.4.0

- docs: update Singularity and Apptainer commands to include --writable-tmpfs flag
- docs: allow collapsing ns/choice items in the docs for processes
- docker: update Dockerfile to include npm cache configuration to allow the pipeline to run on read-only file system
- tests: update SeuratPreparing config to use DoubletFinder for doublet detection
- ci: use latest actions
- deps: add r-clustree as a dependency in docker environment files
- deps: update dependencies in docker environment files (python3.10, R4.3)
- deps: bump biopipen to 0.29.0
  - fix(tcr.TCRClusterStats): fix envs.shared_clusters.heatmap_meta being broken by envs.shared_clusters.sample_order (@<li.ying@mayo.edu>)
  - choir(scrna.SeuratMap2Ref): present better error message when envs.use or values of envs.MapQuery.refdata not in reference (@<li.ying@mayo.edu>)
  - fix(scrna.MarkersFinder): run PrepSCTFindMarkers when needed
  - choir(scrna.SeuratClustering): use FindClusters to run for multiple resolutions
  - choir(scrna.SeuratSubClustering): use FindClusters to run for multiple resolutions
  - feat(scrna.SeuratClustering): add clustree plot (@<li.ying@mayo.edu>)
  - feat(scrna.SeuratSubClustering): add clustree plot
  - tests(scrna.SeuratClusterStats): add assertion for clustree plot generation
- deps: bump biopipen to 0.29.1
  - fix(delim.SampleInfo): fix numbers not split up when each is specified.
  - enh(delim.SampleInfo): make sizes of pie charts proportional to number of samples when each is specified
  - enh(scrna.MarkersFinder): run PrepSCTFindMarkers when necessary before calling FindMarkers
  - feat(scrna.SeuratPreparing): add option to cache Seurat object at different steps
  - feat(scrna.SeuratPreparing): allow doubletfinder to run with a different number of cores
  - chore(scrna.SeuratClustering): record PrepSCTFindMarkers command in sobj@commands
  - tests(scrna.SeuratClusterStats): use less stringent p-value cutoff for DEG/MarkersFinder
  - tests(scrna.SeuratPreparing): add doubletfinder in tests
- deps: bump biopipen to 0.29.2
  - chore(scrna.SeuratClusterStats): use ident label length to adjust default height for feature plots
  - fix(scrna.MetaMarkers): fix seurat object not updated when expanding cases and run `PrepSCTFindMarkers` when necessary before calling meta-markers
  - fix(scrna.MarkersFinder): fix fetching command when composing the PrepSCTFindMarkers command
  - fix(scrna_metabolic_landscape): handle null values in for loop in MetabolicFeatures and MetabolicFeaturesIntraSubset for report generation
- deps: bump biopipen to 0.30.0
  - BREAKING(scrna): move clustree plots from SeuratClustering/SeuratSubClustering to SeuratClusterStats
  - feat(scrna.CellTypeAnnotation): allow to merge/not to merge (envs.merge) the clusters with the same labels predicted
  - feat(scrna.SeuratPreparing): add scDblFinder to detect doublets
  - feat(scrna.SeuratMap2Ref): add envs.skip_if_normalized option to skip normalization if query is already normalized using the same - method as the reference
  - refactor(tcr.Immunarch): source the files for Immunarch scripts for better debugging
  - refactor(scnra.SeuratClustering): refactor the script for better debugging
  - refactor(scnra.SeuratPreparing): refactor the script for better debugging
  - fix(scrna): fix resolution expansion for SeuratClustering and SeuratSubClustering
  - fix(scrna): Fix generating PrepSCTFindMarkers command when no previous commands present
  - tests(scrna.ScFGSEA): fix unavailable urls to GMT files
  - chore(scrna.SeuratMap2Ref): optimize memory usage
  - chore(scrna.MetaMarkers): remove plugin_opts.poplog_max
  - chore(tcr.CloneResidency): improve logging when handling subjects
- deps: bump biopipen to 0.31.3
  - enh(scrna.SeuratMap2Ref): check if reference has SCTModel if SCTransform'ed (likely prepared by old Seurat)
  - fix(tcr.CDR3AAPhyschem): use sequence from TRB chain only
  - fix(tcr.CDR3AAPhyschem): fix when chain is not available
  - fix(tcr.TCRClustering): fix for multi-chain TCRs, use TRB only if on_multi is false
  - fix(tcr.TCRClustering): fix when chain is not available

## 1.3.9

- docs: update docs for TCellSelection to avoid confusion
- deps: bump biopipen to 0.27.9
  - feat(tcr.TCRClusterStats): add sample_order to set sample order on heatmap and cluster_rows to switch row clustering on/off

## 1.3.8

- docs: remove -w option for apptainer/singularity as no writing is necessary since pipen-board 0.15.1
- deps: update biopipen to version 0.27.8
  - fix(scrna.SeuratClusterStats): fix selected columns not unique for stats
  - feat(scrna.SeuratMap2Ref): allow non-SCTransform'ed reference
  - feat(scrna.SeuratMap2Ref): allow splitting query object for mapping (pwwang/immunopipe#61)
  - deps: update pipen-board to version 0.15.1 (allow configuration file path in the URL box on Web UI)

## 1.3.7

- ci: fix docker images building when no essential changes made

## 1.3.6

- ci: fix deploy workflow (#59)
- ci: add README.md to tests-output branch
- ci: fix test/test workflow
- tests: add `make test`
- tests: init test data preparation
- tests: add test for ImmunarchLoading
- tests: add tests for SeuratPreparing
- tests: Update configs for SeuratPreparing test to subset cells so tests can run on CI
- tests: update SeuratPreparing test to disable export
- tests: add tests for SeuratClusteringOfAllCells/SeuratClustering
- docs: update installation instructions (@<stein.mariam@mayo.edu>)
- deps: bump biopipen to version 0.27.7 (0.27.5-0.27.7)
  - fix(scrna.SeuratClusterStats): fix color palette for ridge plots (@<stein.mariam@mayo.edu>)
  - feat(scrna.SeuratPreparing): add envs.cell_qc_per_sample to filter cells before merging instead after
  - fix(scrna_metabolic_landscape.MetabolicFeatures): fix return value of groups with less than 5 cells in `do_one_group`
  - fix(scrna_metabolic_landscape): fix mutaters not working.
  - fix(scrna_metabolic_landscape.MetabolicFeatures/MetabolicFeaturesIntraSubset): skip groups with less than 5 cells in `do_one_group` and save a warning file under the case
  - chore: fix typo in class name `ExprImpution` to `ExprImputation`

## 1.3.5

- ci/test: add tests in CI and deploy output in a different branch
- deps: bump biopipen to 0.27.4
  - choir(delim.SampleInfo): add alpha to the colors of the plots using biopipen color pallete
  - docs(tcr/scrna/scrna_metabolic_landscape): update links of images in docs

## 1.3.4-post

- ci/test: init ci for tests
- docs: introduce versioning for docs

## 1.3.4

- deps: bump biopipen to 0.27.3
  - deps: bump pipen-poplog to 0.1.2 (quick fix for populating logs when job fails)
  - deps: temporary fix copier breaks with pyyaml-include v2 (copier-org/copier#1568)
  - choir(scrna.ScFGSEA): Skip cases when no cells found (#50)
  - choir(scrna.MarkersFinder): Skip cases when no cells found (#50)
  - choir(scrna.MetaMarkers): Skip cases when no cells found (#50)
  - feat(scrna.SeuratPreparing): support `DoubletFinder` (#52)

## 1.3.3

- deps: temporary fix copier breaks with pyyaml-include v2 (copier-org/copier#1568)
- docs: update FAQ.md with instructions for running pipeline on a cluster
- deps: bump biopipen to 0.27.2
  - fix(scrna.RadarPlots): fix mutaters not working
  - feat(tcr.CloneResidency): support `envs.upset_ymax` to set the max value of y axis in upset bar plot.

## 1.3.2

- deps: bump pipen to 0.14.5
- deps: add r-complexupset package to environment.yml and environment_full.yml for CloneResidency
- deps: pin tensorflow to 2.15 for TESSA
- deps: bump biopipen to 0.27.1
  - depr(scrna.MarkersFinder): remove `envs.use_presto` as it's used by Seurat v5 by default
  - enh(tcr.CloneResidency): support log scale for y axis of upset bar plots
  - enh(scrna.SeuratClusterStats): allow to rotate labels in circos plot (pwwang/immunopipe#48) @<li.ying@mayo.edu>
  - enh(scrna.SeuratClusterStats): use `pal_biopipen` for ident colors in circos plot
  - fix(scrna.CellsDistribution): fix the row order of the heatmaps
  - fix(scrna.SeuratClusterStats): fix when `envs.split-by` is specified
  - feat(scrna.CellsDistribution): support `envs.prefix_each`
  - feat(scrna.MarkersFinder): allow to set max number of genes to plot in dotplots
  - feat(scrna.MarkersFinder): support setting detailed arguments for overlapping plots
  - feat(scrna.MarkersFinder): support `envs.prefix_group`
  - feat(scrna.ScFGSEA): support `envs.prefix_each`
  - feat(scrna.RadarPlots): support `envs.prefix_each` and `envs.subset`
  - choir(scrna.SeuratClusterStats): use logger instead of print for log messages
  - choir(tcr.TCRClustering): print session info for `clustcr` script
  - choir(scrna.MarkersFinder): flatten toc when no `section` and no `ident-1` specified
  - docs: add more detailed docs for `envs.section` for multiple processes
  - BREAKING(scrna.SeuratMap2Ref): rename envs.name to envs.ident so envs.MapQuery.refdata is not - required anymore. It will be inferred from envs.ident and envs.use. @<li.ying@mayo.edu>

## 1.3.1

- deps: bump pipen to 0.14.3
- deps: pin ggplot2 to 3.4 for docker due to breaking changes of 3.5
- deps: bump biopipen to 0.26.2
  - deps: bump datar-pandas to 0.5.5 to dismiss deprecated warnings
  - fix(utils.misc.R): replace latin and greek characters with closest ascii chars for slugify()
  - feat(scrna.TopExpressingGenes): support `subset`
  - fix(scrna.CellsDistribution): fix the row order of the heatmaps.
  - enh(tcr.CloneResidency): add legend for multiplets in upset plots.
  - feat(scrna.SeuratClusterStats): add circos plot for cell composition stats (#46).

## 1.3.0

- deps: bump pipen to 0.14.1
- deps: bump pipen-report to 0.18.2
- deps: bump biopipen to 0.26.0
  - fix(scrna.CellTypeAnnotation): keep factor meta data when input and output are RDS for celltypist
  - deps: bump datar to 0.15.4 (support pandas 2.2)
  - fix(utils.single_cell.R): fix `immdata_from_expanded` missing other data columns
  - fix(tcr.Immunarch): fix `mutaters` not working when no subset is set
  - fix(scrna.CellsDistribution): fix `hm_devpars` not working
  - fix(scrna.CellsDistribution): fix multiple `cells_by` columns and speed up plotting
  - choir(tcr.CloneResidency): mark singletons in Venn diagrams more clear
  - fix(scrna.RadarPlots): fix the order of groups on radar plots
  - choir(scrna.RadarPlots): transpose the count/percentage table to save to files
  - fix(scrna.MarkersFinder): fix generating report json file when no significant genes found
  - choir(scrna.MarkersFinder): Plot maximum 20 genes in dotplots
  - choir(scrna.MarkersFinder): Do not convert dashes in case names to dots
  - see more at <https://github.com/pwwang/biopipen/releases/tag/0.26.0>

## 1.2.0

- docs: update FAQs to align with Seurat v5
- docs: add image from manuscript to README.md
- docs: center the flowchart image in README.md
- docs: mention `celltypist` model prep in preparing input data
- deps: bump pipen to 0.13.2
- deps: bump biopipen to 0.25.2:
  - scrna.MarkersFinder: allow to cache `FindAllMarkers` results
  - scrna.CellTypeAnnotation: support `celltypist` (pwwang/biopipen#111)
  - scrna.SeuratSubClustering: add `envs_depth = 1` to replace whole `envs.cases` when new case assigned
  - scrna_metabolic_landscape.MetabolicPathwayHeterogeneit): fix output directory path is not slugified
  - tcr.Immunarch: change case filling log to debug level

## 1.1.1

- deps: Bump biopipen to 0.24.2
  - chore: use internal slugify instead of slugify library
  - tcr.Immunarch: fix spectratyping output file extension is not png
  - scrna.SeuratPreparing: fix displaying filters in report
  - scrna.SeuratPreparing: fix logging Seurat procedure arguments

## 1.1.0

- docs: update table in gallery
- deps: use [pipen-poplog](https://github.com/pwwang/pipen-poplog) to populate job logs to pipeline running log
- deps: bump biopipen to 0.24. Hights:
  - scrna.ScFGSEA: add subset to filter cells (pwwang/biopipen#112) @yuey11
  - scrna.SeuratClustering/SeuratSubClustering: cache Seurat procedures step by step (#40) @xyfqwlzoe
  - tcr.Immunarch: add `plot_type` to support boxplots for diversity metrics
  - see more at <https://github.com/pwwang/biopipen/releases/tag/0.24.0>

## 1.0.5

- change: do not rescale gene expression in `TCellSelection` any more
- fix: fix column names of indicators not aligned with `indicator_genes`
- feat: add feature plots in `TCellSelection`
- deps: bump biopipen to 0.23.8
  - scrna.SeuratPreparing: log Seurat procedure arguments
  - scrna.ScFGSEA: add `subset` to filter cells (pwwang/biopipen#112)

## 1.0.4

- deps: bump biopipen to 0.23.7
  - scrna.SeuratPreparing: update log message for transformation/scaling step
  - scrna_metabolic_landscape.MetabolicPathwayHeterogeneity: add utils.gsea script source to support `localizeGmtfile`

## 1.0.3

- deps: add r-seuratdisk dependency to conda env files. @yuey11
- deps: pin r-matrixstats to 1.1.0 in conda env files to fix `useNames = NA` error. @yuey11
- refactor: optimize configuration file validation
- deps: bump biopipen to 0.23.6
  - feat: support url for gmtfile wherever GSEA is performed (pwwang/biopipen#113)
  - tcr.Immunarch: add error message for empty filtered/subset data in diversity
  - scrna.SeuratPreparing: correct description of default assay in docstr
  - scrna.SeuratPreparing: run also the normal normalization procedures when `SCTransform` is used (useful for visualization purposes on RNA assay)
  - scrna.ModuleScoreCalculator: document the names added by cell cycle score (#34)
  - scrna.SeuratPreparing: support sample names as `reference` for `IntegrateLayers`

## 1.0.2

- deps: add bioconductor-glmgampoi to conda env files (#33)
- docs: correct the Seurat object assay description
- deps: bump biopipen to 0.23.5
  - fix: fix when no enriched items found for `scrna.MarkersFinder`, `scrna.MetaMarkers` and `scrna.TopExpressingGenes`
  - scrna.SeuratClusterStats: fix when `frac` or `frac_ofall` is true and no `group-by` nor `split-by` is specified for `stats`
  - utils.gsea.R: fix when no enriched items found for `runEnrichr`
  - scrna_metabolic_landscript: fix adding report when `ncores` > 1

## 1.0.1

- docs: add gallery section to README.md
- change: set default `nstart` of kmeans to 25 in `TCellSelection`
- deps: add `r-hdf5r` in conda env files to support `Read_10x_h5` from Seurat. @yuey11
- deps: bump biopipen to 0.23.4
  - scrna.TopExpressingGenes: fix colnames while pulling average expression
  - scrna.CellsDistribution: fix when `cells_by` has multiple column names
  - scrna.CellTypeAnnotation: fix the order of the clusters for `direct` method
  - scrna.SeuratClusterStats: add `position` options for bar plots for stats
  - scrna.RadarPlots: add `colors` to set the colors of the loops in radar and bar plots
  - tcr.Immunarch: add `split_by` and `split_order` to put subplots together in one single plots

## 1.0.0

### Highlights

- feat: support `Seurat` v5 (integration is now down by `Seurat::IntegrateLayers`)
- feat: support supervised clustering (mapping cells to reference by `Seurat`)
- feat: support dataset with scRNA-seq data only (no scTCR-seq data)
- feat: support diffusion map calculation (by `ModuleScoreCalculator`)
- feat: support subclassing to cluster subsets of cells (by `SeuratSubClustering`)
- feat: allow to ignore TCR data in `TCellSelection` and pass kmeans arguments
- feat: allow to set multiple resolutions (`envs.FindClusters.resolution`) in `SeuratClustering`/`SeuratClusteringOfTCells`
- change: change unsuperved cluster labels to `c1`, `c2`, ... in `SeuratClustering` by default
- docs: add gallery, which contains real-world examples of datasets from publications

### Breaking changes

- change: rename `SeuratMetadataMutater` to `IntegratingTCR`
- change: rename `SeuratClusteringOfTCells` to `SeuratClustering`
- change: rename `TCRClusters2Seurat` to `IntegratingTCRClusters`
- refactor: make `SeuratClustering` (instead of `SeuratClusteringOfAllCells`) work for all cells when all are T cells
- change: move data preparation and integration from `SeuratClustering` to `SeuratPreparing`
- change: default `mode` of `ImmunarchLoading` to `paired` (instead of `single`), which requires both alpha and beta chains (instead of beta chain only) to define a clonotype
- change: default `dbs` for enrichment analysis wherever applies to `KEGG_2021_Human` and `MSigDB_Hallmark_2020`

### Changes

- feat: make `TopExpressingGenes` optional
- feat: add `validate_config` to validate configuration schematically

### Features

- feat(SeuratPreparing): allow to filter genes directly (by specifying `envs.gene_qc.excludes`)
- feat(SeuratClusterStats): add `ngenes` to plot the number of genes expressed in each cluster
- feat(SeuratClusterStats): add `barplot` for features and allow aggregation of features
- feat(SeuratClusterStats): add `envs.mutaters` to mutate meta data
- feat(SeuratClusterStats): add histograms to plot number of cells against another variable
- feat(SeuratClusterStats): Add `frac_ofall` and `transpose` for stats to calculate fraction within group or against all cells, and transpose ident and group, respectively

### Dependencies

- deps: add `r-presto` to conda environment files to support using presto to fastly find markers
- deps: add `bioconductor-destiny` to conda environment file to support add diffusion map components in ModuleScoreCalculator
- deps: add `r-harmony` to support harmony integration by Seurat v5 in conda env file
- deps: add `r-sf` to conda env file
- deps: remove `vdjtools` from conda env files
- deps: bump `pipen-report` to [0.16.3](https://github.com/pwwang/pipen-report/releases/tag/0.16.3)
- deps: bump `biopipen` to [0.23.3](https://github.com/pwwang/biopipen/releases). Hightlight changes:
  - scrna.MarkersFinder: Add `envs.use_presto` to use presto to speed up finding markers
  - scrna.SeuratPreparing: Set envs.gene_qc.min_cells to 0 by default (instead of 3)
  - scrna.ScFGSEA: Allow to ignore small group when fgsea fails due to all NAs for pre-ranks
  - scrna.CellsDistribution: Allow to order clusters by envs.cluster_orderby
  - scrna.CellsDistribution: Add heatmaps
  - tcr.CloneResidency: Make section works in report
  - tcr.Immunarch: Support paired chain data for VJ conjuction plots
  - tcr.TESSA: Change envs.assay to None to use default assay of Seurat object
  - scrna.SeuratClusterStats: Add `avgheatmap` to plot more elegant heatmap for average gene expressions
  - scrna.SeuratClusterStats: Fix ident not working for dimplots
  - scrna.SeuratClusterStats: Add `cluster_orderby` to order clusters for features
  - scrna.SeuratClusterStats: Add `na_group` to keep NA values in `group-by`
  - utils.mutate_helpers: Change arguments `id_col` and `compare_col` of `paired` to `id` and `compare`, respectively
  - utils.mutate_helpers: Fix that subset can't be an expression for expanded family
  - utils.mutate_helpers: Add `top` to select top entities (e.g clones)
  - scrna.RadarPlots: Add breakdown and test to break down the cell distribution and run statistic test on the fractions

## 0.11.2

- docs: move `Immunarch` to the later position in process list
- docs: Use `master` tag in getting-started

## 0.11.1

- chore: change line length to 88 for flake8
- chore: dismissing warning about wasting columns for `SeuratClusteringOfTCells`
- docs: update CHANGELOG.md with missing changes of last version
- docs: add version of renaming `envs.tcell_indicator` to `envs.tcell_selector`
- docs: remove unused doc files
- docs: add metadata illustration
- deps: bump biopipen to 0.22.8. Highlights:
  - deps: bump pipen-board to 0.13.10 (pipen-report to 0.16.2)
  - CellsDistribution: Don't add rownames to the output table file
  - MarkersFinder (ClusterMarkers/ClusterMarkersOfAllCells): Optimize to use `FindAllMarkers` if `ident.1` is not specified
  - SeuratClusterStats: Fix path of expression table file
  - CellTypeAnnotation: Allow using `NA` to exclude clusters from output `Seurat` object
  - utils.mutate_helpers: Return ids only when subset is true and group is not `NA` for `uniq = TRUE` in `expanded`, `collapsed`, `emerged` and `vanished`

## 0.11.0

- deps: update biopipen to 0.22.1, highlights:
  - add V-J junction circos plots to `Immunarch` process
  - add cache option to cache the clustering results if nothing changed except ncores, to `SeuratClustering` process
  - add dot plots to `MarkersFinder` (`ClusterMarkersOfAllCells`, `ClusterMarkers`) process
  - save exported table with only necessary columns for `CellsDistribution` process
  - add `descr` to describe cases cases in report for `CellsDistribution` process
  - add `subset` for dimplots in `SeuratClusterStats` process
  - use a new palette (`biopipen`) for related processes
  - optimize report rendering (using `render_job()` filter from `pipen-report`)
  - change metacols to extracols so essential columns get exported for `ImmunarchLoading` process
  - add cache option to cache the clustering results if nothing changed except ncores for `SeuratClustering` (`SeuratClusteringOfAllCells`) process
  - see more at <https://github.com/pwwang/biopipen/releases/tag/0.22.0> and <https://github.com/pwwang/biopipen/releases/tag/0.22.1>
- deps: update pipen-report to 0.16, highlights:
  - scroll anchor into view on the page
  - build report page when each process is done, instead of the whole pipeline
  - see more at <https://github.com/pwwang/pipen-report/releases/tag/0.16.0>
- change: remove `Immunarch2VDJtools` and `VJUsage` processes (vj usage analysis can be done in `Immunarch` process)
- change: change `tcell_indicator` to `tcell_selector` in `TCellSelection` process
- enhance: provide better error message when none barcode matches from RNA and TCR data for `TCRClustering` process
- docs: add memory usage reduction tips in FAQ
- chore: dismiss warnings of wasted input columns for multiple processes

## 0.10.1

- chore: update pipeline description to include version in the logs
- fix: add fc-cache command to Dockerfile to solve `Fontconfig error`
- docker: optimize building full image based off the base image

## 0.10.0

- docker: lock r-matrix version to 1.6_1 for compatibility
- docs: adopt mkdocs-rtd 0.0.10 (add scrollbar to the table of contents)
- deps: bump biopipen to 0.21.1
  - use `r-logger` for logging in R scripts
  - docs: fix internal references in API docs
  - deps: bump pipen-board to 0.13.6
  - SampleInfo: refactor data subset logic using `subset` instead of `distinct`
  - Immunarch: add `in.metafile` to allow other meta info (i.e. seurat clusters) for future subsetting (#22)
  - Immunarch: fix empty groups in diversity plot after subsetting
  - Immunarch: allow `subset` to subset cells for analyses
  - Immunarch: allow `separate_by` also works on other diversity plots
  - Immunarch: add `ymin` and `ymax` to align diversity plots by `separate_by`
  - Immunarch: add `ncol` to specify # columns in the combined plots
  - RadarPlots: fix `envs.order` not working
  - MarkersFinder: add `overlap` to find overlapping markers between cases (#24)
  - MarkersFinder: allow `subset` to subset cells for analyses
  - MarkersFinder: add dot plots for significant markers
  - CellsDistribution: allow multiple columns for `cells_by`
  - CellsDistribution: allow `subset` to subset cells for analyses
  - utils.mutate_helpers.R: add `include_emerged` for `expanded()` and `include_vanished` for `collapsed()`

## 0.9.3

- deps: Bump biopipen to 0.20.7
- deps: Bump pipen-board to 0.13.4
- ClusterMarkers/ClusterMarkersOfAllCells: Choose avg_log2FC > 0 markers by default
- MarkersFinder: Allow to set assay and set assay to `RNA` by default
- CellsDistribution: Add venn/upset plot for overlapping cell groups in different cases
- SampleInfo: Add `distinct` to case to perform stats on distinct records

## 0.9.2

- ➕ Add `r-ggnewscale` as dependency for `CDR3AAPhyschem` in docker image
- ⬆️ Bump biopipen to 0.20.5
  - 🧱 CloneResidency: Integrate RNA data to allow more flexible analysis (i.e. within specific seurat clusters)
  - 🏗️ CloneResidency: Rename envs.sample_groups to envs.section to be consistent with other processes
  - 📝 ScFGSEA: Remove the link in the summary of the docstring (since they are not transformed in the report)
  - 🎨 CDR3AAPhyschem: Give better error message when wrong group items are given
- ⬆️ Bump pipen-board to 0.13.3
  - Add items automatically when blurred for list options
  - Add other sections to description on the UI for processes

## 0.9.1

- 🐛 Fix docstring for `RadarPlots`
- ➕ Add `pipen-diagram` as dependency
- ➕ Set `pipen-runinfo` as optional
- ⬆️ Bump biopipen to 0.20.4
- 📝 Update version in docs

## 0.9.0

### Housekeeping and docs

- Bump biopipen to 0.20.3 (pipen to 0.12)
- Use [`pipen-cli-ref`](https://github.com/pwwang/pipen-cli-ref) to generate API for processes (it uses docstring of the process class so that we don't need to maintain two copies of docs)

### Fixed/Enhanced

- Make `/data` directory in container, so it can be mounted
- Fix a bug when a single gene provided to `indicator_genes` in `TCellSelection`
- Move `ModuleScoreCalculator` before clustering so that the scores can be used in `vars.to.regress` of `SCTransform` while clustering
- Set default assay to RNA in case module scores only caculated using integrated features in `ModuleScoreCalculator`
- Improve QC plots in `SeuratPreparing` by marking the cells that are removed in the plots instead of doing before/after plots
- Fix type annotation for envs.features_defaults.ncol in docstring for `SeuratPreparing` (causing `pipen-board` not converting to int)
- Fix the cluster order in pie charts for `CellsDistribution`
- Fix the cluster order in pie charts for `SeuratClusterStats`
- Fix order in pie charts for `SampleInfo`
- Fix docstring for `envs.div.args` of `Immunarch` (more clear description of method)
- Allow mutiple columns in the file for `envs.features_defaults.features` in `SeuratClusterStats`
- Allow order to be optional for `CloneResidency` (errored when not provided)
- Add number of clusters at the end of log for `SeuratClusteringOfAllCells`/`SeuratClusteringOfTCells`
- Add stricter checker for input file (#13)
- Indicate the case name in logs when pie is enabled for `group-by` in `SeuratClusterStats`
- Allow to skip overlap and gene usage analyses by setting method to `none` for `Immunoarch` (#11, #12)
- Don't cluster on heatmap when there are only 2 samples for `TCRClusterStats` (#11)
- Import Seurat explictly to avoid satijalab/seurat#2853 in `MetabolicFeatures`
- Fix when NA values in data for heatmap in `MetabolicPathwayActivity`
- Fix error when no significant pathways selected in `MetabolicPathwayHeterogeneity`
- Give better error message in CellsDistribution if group value not found for `CellsDistribution` (#16)
- Try including more genes (even though insignificant) in volcano plot for `MarkersFinder`/`ClusterMarkers`/`ClusterMarkersOfAllCells` (#17)
- Add margins to volcano plot for `MarkersFinder`/`ClusterMarkers`/`ClusterMarkersOfAllCells`
- Fix when `envs.cell_qc` is `None` (not provided) for `SeuratPreparing`
- Fix `ident` in cases of `envs.dimplots` not working for `SeuratClusterStats`

### Added

- Add `ClusterMarkersOfAllCells` and `TopExpressingGenesOfAllCells` and set them as optional
- Add dim plots in `SeuratClusterStats` to overlay TCR presence/absence of cells (#14)

### Breaking changes-0.9.0

- Rename `TCRClusteringStats` to `TCRClusterStats` (#15)

## 0.8.3

- 📝 Fix typos in docs
- 📝 Add links to some optional input files (#9, 5)
- 🔨 Add apptainer to docker entry.sh (#9, 6)
- 💄 Adjust process order in reports (#9, 1)
- ⬆️ Bump pipen-report to 0.13.1 (#9, 2)

## 0.8.2

- Bump biopipen to 0.18.3 to fix when either ident is empty for `MarkersFinder`

## 0.8.1

- Bump biopipen to 0.18.2 to fix a bug when the min length of CDR3 seqs > 12 for `CDR3AAphyschem`

## 0.8.0

### Housekeeping and docs updates

- Bump biopipen to 0.18.1
- Mention function changes with versions in docs
- Add apptainer in board.toml so the command can be generated in pipen-board
- Make logo shorter in docs
- Add docker image with `-full` tags to include all dependencies
- Print command help message if run test failed in CI
- Add singularity/apptainer in FAQ for "no space left" question
- Add -w fro apptainer in docs (as we need to save pipen-board file in home directory)

### Added-0.8.0

- Add `TESSA` process for [tessa analysis](https://pwwang.github.io/immunopipe/processes/TESSA/)
- Add volcano plot for `MarkersFinder` and `ClusterMarkers`

### Fixed

- Fix when `Sample` is the only column in meta for `ImmunarchLoading`
- Add clear message when `k.weight` is too large for `IntegrateData` in `SeuratClustering`
- Allow `unique:` prefix for `on` in `SampleInfo`
- Fix sample order in plots for `SampleInfo`
- Remove `tidyseurat::` prefix for `filter` in scripts of `MetaMarkers`, `ScFGSEA` and `SeuratClusterStats` in case `tidyseurat::filter` is not exported when installed from `conda` (but it will make `dplyr::filter` work anyway on seurat object)

### Breaking changes-0.8.0

- Redesign envs for `SeuratClusteringStats` to allow setting defaults for cases and switch identities for plots

## 0.7.0

### Housekeeping and docs updates-0.7.0

- Fix typos in docs/configurations
  - `TCRClustering` should be `TCRClusteringStats` in Multi-case variable design section
  - `infile` of `[SampleInfo.in]` should be `samples.txt` rather than `sample.txt`
- Remove unused scripts by deprecated processes
- Bump `pipen-report` to [0.12.8](https://github.com/pwwang/pipen-report/releases/tag/0.12.8)
- Add `master` branch and `master` tag as stable tag for docker image
- Add pdf version of the flowchart (#4)
- Add warning for the results in getting started tutorial
- Bump `pipen-board` to [0.11.5](https://github.com/pwwang/pipen-board/releases/tag/0.11.5)
- Add apptainer to the docs

### Added-0.7.0

- Add `ModuleScoreCalculator` to calculate module scores or cell cycle scores
  - See: <https://pwwang.github.io/immunopipe/processes/ModuleScoreCalculator/>
- Allow `SampleInfo` to perform statistics on the sample information
  - See: <https://pwwang.github.io/immunopipe/processes/SampleInfo/>
- Add `TCR_Cluster_Size` and `TCR_Cluster_Size1` from `TCRClustering` to metadata for further integrative analysis
  - See: <https://pwwang.github.io/immunopipe/processes/TCRClusters2Seurat/>

### Fixed-0.7.0

- Fix default height and width for plots in `SeuratClusterStats`
- Fix cluster order not kept after annotation using `hitype` in `CellTypeAnnotation`

### Breaking changes-0.7.0

- Change `seurat_clusters_old` to `seurat_clusters_id` to save old `seurat_clusters` in `CellTypeAnnotation`
- Remove `MarkersForClustersOfAllCells` and `TopExpressingGenesOfAllCells` processes
- Rename `MarkersForClustersOfTCells` to `ClusterMarkers`
- Rename `TopExpressingGenesOfTCells` to `TopExpressingGenes`
- Rename `envs.exprs` to `envs.features` for `SeuratClusterStats`
  - `envs.exprs.genes` is also renamed to `envs.features.features`

## 0.6.0

- ⬆️ Bump biopipen to 0.16
- 📝 Add documentation
- 💚 Fix docs building in CI
- 📝 Update README with flowchart

## 0.5.1

- ✨ Add `TopExpressingGenes`
- 🎨 Move `RadarPlots` to `biopipen`
- ⬆️ Bump biopipen to 0.15.2

## 0.5.0

- ⬆️ Upgrade biopipen to 0.15.0
- 💚 Use better strategy docker image building

## 0.4.0

- ⬆️ Bump biopipen to 0.6
- ⬆️ Upgrade other dependencies
- 💚 Use micromamba for docker image building
- ⬆️ Add procps-ng for vdjtools for docker  building

## 0.3.0

- 💚 Use build 2 for genomeinfodbdata from bioconda (0.2.4)
- 👽️ Use config from pipen_args
- ⬆️ Pump biopipen to 0.5.3, pipen-args to 0.3.2
- ⬆️ Upgrade deps for docker
- 📝 Add flowchart in README.md
- 🐛 Fix error when --config not passed

## 0.2.4

- 💚 Use lastest miniconda3 for docker build
- 💚 Use conda channel pwwang for bioconductor-genomeinfodbdata for fix (bioconda/bioconda-recipes#31349)
- ⬆️ Upgrade biopipen to 0.4.9
- 📝 Add URL to example in README

## 0.2.3

- ⬆️ Upgrade biopipen to 0.4.8

## 0.2.2

- ⬆️ Upgrade biopipen to 0.4.7 to fix SeuratPreparing

## 0.2.1

- 🔥 Fix the bug of the wrong arguments in help page
- ⬆️ Upgrade clustcr to 1.0.2
- 📝 Fix docs for metabolic analysis

## 0.2.0

- ♻️ Move in-house processes out of processes.py
- ♻️ Split up MARKERS_FINDER
- ♻️ Refactor RadarPlots
- ✨ Add an example config file
- ⚡️ Add `filter` for RadarPlots
- 📝 Update docs
- ⬆️ Upgrade deps
- 🔧 Update docker/environment.yml
- 🐛 Fix CloneHeterogeneity when only 1 row in continency table

## 0.1.1

- 💚 Try fix pip in environment.yml
- 📝 Update readme for requirement checking
- 📝 Update docs to fix #1
- 📝 Update CHANGELOG
- ⬆️ Adopt biopipen 0.4.0

## 0.1.0

- 🩹 Disable force-caching for some procs
- ⬆️ Upgrade datar to 0.8.*
- ✨ Add dockerfile
- ⬆️ Upgrade pipen to 0.3
- 💥 Remove gene lists from start processes
- ⬆️ Upgrade biopipen to 0.3
- ⬆️ Upgrade pipen to 0.3.5

## 0.0.7

- Add CloneHeterogeneity
- Allow setting `indicator_gene` for `TCellSelection`
- Adopt latest datar and biopipen

## 0.0.6

- ✨ Allow dimplots with clonal information

## 0.0.5

- ✨ Allow more flexible dim plots

## 0.0.4

- ✨ Refactor markers finder module and add meta-marker analysis

## 0.0.3

-✨ Add metabolic pathway analysis

## 0.0.2

- Adopt biopipen 0.1.3

## 0.0.1

- First release
