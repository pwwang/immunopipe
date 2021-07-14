library(dplyr)
library(tibble)
library(foreach)
library(doParallel)
library(Seurat)


samples = "{{ in.samples }}"
tcrcount_dir = "{{ in.tcr_counts }}"
exprdir = "{{ in.exprdir }}"
septdir = "{{ in.septdir }}"
outdir = "{{ out.outdir }}"
seurat_source = "{{ args.seurate_source }}"
ncores = {{ args.ncores }}

registerDoParallel(ncores)
source(seurat_source)
load(samples)

dir.create(outdir, showWarnings = FALSE)

# Matrices-1
matrices_1 = file.path(outdir, "Matrices-1")
dir.create(matrices_1, showWarnings = FALSE)
foreach (prefmat = Sys.glob(file.path(exprdir, "*.mat.rds"))) %dopar% {
    prefix = sub(".mat.rds", "", basename(prefmat), fixed=T)
    divide.mat(prefix, septdir, exprdir, matrices_1)
}

# get paired tumor-normal patients
samples = samples %>%
    filter(Type=='scRNA') %>%
    group_by(Patient) %>%
    filter(n() == 2)

# Seurat-1
seurat_1 = file.path(outdir, "Seurat-1")
dir.create(seurat_1, showWarnings = FALSE)
foreach (i=seq_len(nrow(samples))) %dopar% {
    row = samples[i, ]
    process.mat(
        row$Prefix,
        row$Source,
        row$Patient,
        seurat_1,
        matrices_1,
        tcrcount_dir,
        cell="tcell"
    )
}

foreach (group=group_split(samples)) %dopar% {
    combine2(
        unlist(group$Patient)[1],
        unlist(group$Prefix),
        seurat_1,
        cell="tcell"
    )
}

# Combine-1
# combine_1 = file.path(outdir, "Combine-1")
# dir.create(combine_1, showWarnings = FALSE)
global.obj = combine_all_patients(
    samples %>% pull(Patient) %>% unique() %>% as.character(),
    seurat_1,
    cell="tcell",
    outpref="global.tnb.tcell"
)

# Markers
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
    write.table(
        markers,
        file.path(markers_dir, paste0("Ident.",ident,".markers.txt")),
        col.names = TRUE,
        row.names = TRUE,
        sep="\t",
        quote = FALSE
    )
}
