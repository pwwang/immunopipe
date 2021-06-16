library(dplyr)
library(foreach)
library(doParallel)

samples = "{{ in.samples }}"
outdir = "{{ out.outdir }}"
ncores = {{ args.ncores }}

registerDoParallel(ncores)

get.cellranger <- function (sample, prefix, path, genes, matfile) {

    fp <- gzfile(path, "r")
    cell.info <- read.delim2(fp,header=FALSE)
    close(fp)
    cell.names <- paste(prefix, cell.info[,1], sep="_")

    # Reads columns in as factors
    fp <- gzfile(genes,"r")
    gene.info <- read.delim2(fp,header=FALSE)
    close(fp)

    # Handle problem with GeneID:1454, a readthrough also annotated as CSNK1E
    gene.names <- as.character(gene.info[,2])
    gene.names <- gsub("_","-",gene.names)  # Seurat v3 cannot handle underscores
    gene.names[which(gene.info$V1 == "GeneID:102800317")] <- "LOC400927-CSNK1E"

    fp <- gzfile(matfile,"r")
    line <- readLines(fp,n=2)
    line <- readLines(fp,n=1)
    fields = unlist(strsplit(line," "))
    ngenes <- as.numeric(fields[1])
    ncells <- as.numeric(fields[2])

    # Put genes by row, as needed by Seurat
    mat <- matrix(0,nrow=ngenes,ncol=ncells)

    while (TRUE) {
        line <- readLines(fp,n=1)
        if (length(line) == 0) {
            break
        }
        fields = unlist(strsplit(line," "))

        mat[as.numeric(fields[1]),as.numeric(fields[2])] <- as.numeric(fields[3])
    }
    close(fp)

    rownames(mat) <- gene.names
    colnames(mat) <- cell.names
    mat
}



dir.create(outdir, showWarnings = FALSE)
load(samples)

samples = samples %>% filter(Type == 'scRNA')

foreach(patient=unique(samples$Patient)) %dopar% {
    sdata = filter(samples, Patient == patient)
    print(paste('Doing patient:', patient, '...'))
    for (j in seq_len(nrow(sdata))) {
        sample = as.character(sdata[j, 'Sample'])
        prefix = as.character(sdata[j, 'Prefix'])
        path = as.character(sdata[j, 'Path'])
        genes = as.character(sdata[j, 'Genes'])
        matfile = as.character(sdata[j, 'Matrix'])

        mat = get.cellranger(sample, prefix, path, genes, matfile)
        saveRDS(mat, file=file.path(outdir, paste0(prefix,".mat.rds")))
    }
}
