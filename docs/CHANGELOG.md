# Change Log

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
