# Change Log

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

### Breaking changes

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

### Added

- Add `TESSA` process for [tessa analysis](https://pwwang.github.io/immunopipe/processes/TESSA/)
- Add volcano plot for `MarkersFinder` and `ClusterMarkers`

### Fixed

- Fix when `Sample` is the only column in meta for `ImmunarchLoading`
- Add clear message when `k.weight` is too large for `IntegrateData` in `SeuratClustering`
- Allow `unique:` prefix for `on` in `SampleInfo`
- Fix sample order in plots for `SampleInfo`
- Remove `tidyseurat::` prefix for `filter` in scripts of `MetaMarkers`, `ScFGSEA` and `SeuratClusterStats` in case `tidyseurat::filter` is not exported when installed from `conda` (but it will make `dplyr::filter` work anyway on seurat object)

### Breaking changes

- Redesign envs for `SeuratClusteringStats` to allow setting defaults for cases and switch identities for plots

## 0.7.0

### Housekeeping and docs updates

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

### Added

- Add `ModuleScoreCalculator` to calculate module scores or cell cycle scores
    - See: <https://pwwang.github.io/immunopipe/processes/ModuleScoreCalculator/>
- Allow `SampleInfo` to perform statistics on the sample information
    - See: <https://pwwang.github.io/immunopipe/processes/SampleInfo/>
- Add `TCR_Cluster_Size` and `TCR_Cluster_Size1` from `TCRClustering` to metadata for further integrative analysis
    - See: <https://pwwang.github.io/immunopipe/processes/TCRClusters2Seurat/>

### Fixed

- Fix default height and width for plots in `SeuratClusterStats`
- Fix cluster order not kept after annotation using `hitype` in `CellTypeAnnotation`

### Breaking changes

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
