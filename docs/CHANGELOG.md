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
