library(dplyr)
# library(multidplyr)
library(tibble)
library(foreach)
library(doParallel)

samples = "{{ in.samples }}"
tcrcount_dir = "{{ in.tcr_counts }}"
exprdir = "{{ in.exprdir }}"
outdir = "{{ out.outdir }}"
seurat_source = "{{ args.seurate_source }}"
ncores = {{ args.ncores }}

set.seed(8525)
registerDoParallel(ncores)
# cluster = new_cluster(ncores)
source(seurat_source)
dir.create(outdir, showWarnings = FALSE)

load(samples) #samples

# get paired tumor-normal patients
samples = samples %>%
    filter(Type=='scRNA') %>%
    group_by(Patient) %>%
    filter(n() == 2)

# Sample   Type  Patient  Source Path        Genes         Matrix        Prefix
# <fct>    <fct> <fct>    <fct>  <chr>       <chr>         <chr>         <chr>
# 1 MM003BM… scRNA MM003-e… BM     .pipen/sam… /research/la… /research/la… BM-9
# 2 MM003WB… scRNA MM003-e… WBC    .pipen/sam… /research/la… /research/la… WBC-9
# 3 MM005BM… scRNA MM005-e… BM     .pipen/sam… /research/la… /research/la… BM-12
# 4 MM005WB… scRNA MM005-e… WBC    .pipen/sam… /research/la… /research/la… WBC-12

print("Processing expression matrices ...")
foreach (i=seq_len(nrow(samples))) %dopar% {
    row = samples[i, ]
    process.mat(
        row$Prefix,
        row$Source,
        row$Patient,
        outdir,
        exprdir,
        tcrcount_dir
    )
}

print("Combining samples from the same patient ...")
foreach (group=group_split(samples)) %dopar% {
    combine2(
        unlist(group$Patient)[1],
        unlist(group$Prefix),
        outdir
    )
}

print("Merging all patients ...")
global.obj = combine_all_patients(
    samples %>% pull(Patient) %>% unique() %>% as.character(),
    outdir
)

print("Finding markers for each cluster ...")
markers_dir = file.path(outdir, "Markers")
dir.create(markers_dir, showWarnings = FALSE)

idents = names(table(Idents(global.obj)))
foreach (ident = idents) %dopar% {
    print(paste("* For ident", ident, "..."))
    markers <- FindMarkers(global.obj, ident.1=ident)
    save(
        markers,
        file=file.path(markers_dir, paste0("Ident.",ident,".markers.RData"))
    )
}
