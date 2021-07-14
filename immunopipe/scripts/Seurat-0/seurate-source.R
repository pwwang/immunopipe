library(Seurat)

# Use sctransform.  Replaces NormalizeData, ScaleData, and FindVariableFeatures
process.mat <- function (
    prefix, # BM1
    source, # BM
    patient, # MM003-eariler
    outdir,
    exprdir,  # Matrices-0
    tcrcount_dir, # TCR-counts
    cell="",
    find.clusters=FALSE
) {
    if (nzchar(cell)) {
        suf = paste0(".", cell)
    } else {
        suf = ""
    }
    matfile = paste0(prefix, suf, ".mat.rds")
    mat <- readRDS(file.path(exprdir, matfile))
    clonotypes.file = paste0(patient, ".pheno.", source)
    pheno <- read.table(file.path(tcrcount_dir, clonotypes.file), header=TRUE)
    idx <- match(colnames(mat), pheno$barcode)

    # Now done in Matrices-0
    #colnames(mat) <- paste(prefix, colnames(mat), sep="_")
    obj <- CreateSeuratObject(counts=mat, project=prefix)
    obj$patient <- patient
    obj$sample <- prefix
    obj$source <- source
    obj$clonotype <- pheno$clonotype[idx]
    obj$gamma <- pheno$gamma[idx]
    obj$delta <- pheno$delta[idx]
    obj$igh <- pheno$igh[idx]
    obj$igk <- pheno$igk[idx]
    obj$igl <- pheno$igl[idx]

    # Puts results into SCT assay.  Trying to regress nUMI fails
    # vars.to.regress="source" -- we are processing each source separately
    obj <- SCTransform(object=obj, return.only.var.genes=FALSE, verbose=FALSE)
    outfile = file.path(outdir, paste0(prefix, suf, ".obj.rds"))
    saveRDS(obj, file=outfile)
}


#sample <- "lung1.tn"
#tumorfile <- "lt1.obj.RData"
#normalfile <- "ln1.obj.RData"

combine2 <- function(
    sample,
    prefices,
    outdir,
    cell="",
    find.clusters=FALSE
) {
    if (tolower(prefices[1]) %in% c('wbc', 'blood', 'normal', 'control', 'ctrl')) {
        normal_pref = prefices[1]
        tumor_pref = prefices[2]
    } else {
        normal_pref = prefices[2]
        tumor_pref = prefices[1]
    }
    if (nzchar(cell)) {
        suf = paste0(".", cell)
    } else {
        suf = ""
    }
    tumorfile = file.path(outdir, paste0(tumor_pref, suf, ".obj.rds"))
    normalfile = file.path(outdir, paste0(normal_pref, suf, ".obj.rds"))
    tumor.obj <- readRDS(tumorfile)
    normal.obj <- readRDS(normalfile)

    set.seed(31415926)
    reference.list <- list()
    reference.list[[1]] <- tumor.obj
    reference.list[[2]] <- normal.obj

    k.filter <- 200
    k.filter <- min(k.filter, ncol(tumor.obj))
    k.filter <- min(k.filter, ncol(normal.obj))
    combined.anchors <- FindIntegrationAnchors(
        reference.list,
        dims=1:30,
        k.filter=k.filter
    )

    # Maximize data that is put into the combined object.
    # Otherwise only var.features are stored
    joint.features <- intersect(
        rownames(tumor.obj@assays$SCT@data),
        rownames(normal.obj@assays$SCT@data)
    )
    combined.obj <- IntegrateData(
        combined.anchors,
        features.to.integrate=joint.features
    )

    saveRDS(combined.obj, file=file.path(outdir, paste0(sample, suf, ".rds")))
}

combine_all_patients = function(samples, outdir, cell="", outpref="global.1.6") {
    if (nzchar(cell)) {
        suf = paste0(".", cell)
    } else {
        suf = ""
    }

    reference.list <- list()
    i <- 1
    for (sample in samples) {
        file <- file.path(outdir, paste0(sample, suf, ".rds"))
        reference.list[[i]] <- readRDS(file)
        i <- i + 1
    }

    combined.anchors <- FindIntegrationAnchors(reference.list, dims=1:30)
    global.obj <- IntegrateData(combined.anchors, dims=1:30)
    rm(combined.anchors)

    # source("./samples.R")
    patient <- factor(global.obj$sample, levels=samples)
    # levels(patient) <- patient.levels
    global.obj@meta.data$patient <- patient
    #saveRDS(global.obj, file="global.integrated.rds")

    global.obj <- ScaleData(global.obj, verbose=FALSE)
    global.obj <- RunPCA(global.obj, npcs=30, verbose=FALSE)

    global.obj <- RunUMAP(global.obj, reduction="pca", dims=1:30)
    #global.obj <- RunTSNE(global.obj, reduction="pca", dims=1:30)

    global.obj <- FindNeighbors(global.obj)
    #saveRDS(global.obj, file="global.precluster.rds")

    global.obj <- FindClusters(global.obj, resolution=.8)
    save(global.obj, file=file.path(outdir, paste0(outpref, ".RData")))
    #saveRDS(global.obj@assays$integrated@scale.data, "all_integrated.rds")

    p2 <- DimPlot(global.obj, reduction="umap", label=FALSE)
    write.table(
        cbind(
            p2$data,
            patient=global.obj$patient,
            sample=global.obj$sample,
            source=global.obj$source,
            clonotype=global.obj$clonotype
        ),
        file=file.path(outdir, paste0(outpref, ".umap")),
        quote=FALSE,
        sep="\t",
        col.names=NA,
        row.names=TRUE
    )

    global.obj
}

divide.mat <- function (prefix, septdir, exprdir, outdir) {
    prefix.uc <- paste0(toupper(prefix),"_")
    tt <- c()
    file <- file.path(septdir, "tcell.main.xls")

    xls <- read.csv(file, row.names=1, sep="\t")
    tt <- rbind(tt, xls[grep(prefix.uc,rownames(xls)),])
    tcell.barcodes <- rownames(tt)

    # nn <- c()
    # file <- file.path(septdir, "nont.main.xls")

    # xls <- read.csv(file, row.names=1, sep="\t")
    # nn <- rbind(nn, xls[grep(prefix.uc,rownames(xls)),])
    # nont.barcodes <- rownames(nn)

    mat <- readRDS(file.path(exprdir, paste0(prefix,".mat.rds")))
    tcell.mat <- mat[,which(colnames(mat) %in% tcell.barcodes)]
    # nont.mat <- mat[,which(colnames(mat) %in% nont.barcodes)]

    saveRDS(tcell.mat, file=file.path(outdir, paste0(prefix,".tcell.mat.rds")))
    # saveRDS(nont.mat,file=paste0(prefix,".nont.mat.rds"))
}
