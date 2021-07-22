
# clonotypes is a vector of length(ncells) -> clonotype
# clonotype.subtypedist is a matrix of length(nclonotypes) x nsubtypes -> count of cells
# clonotype.subtypedist.tumor is a matrix of length(nclonotypes) x nsubtypes -> count of cells in tumor
# clonotype.subtypedist.nat is a matrix of length(nclonotypes) x nsubtypes -> count of cells in nat
# clonotype.subtypedist.blood is a matrix of length(nclonotypes) x nsubtypes -> count of cells in blood
# clonotype.subtype is a vector of length(nclonotypes) -> subtype
# clonotype.tcounts is a vector of length(nclonotypes) -> count of cells
# clonotype.ncounts is a vector of length(nclonotypes) -> count of cells


compute.counts.1 <- function (sample, source, outdir) {

    pheno <- read.table(file.path(outdir, paste0(sample,".pheno.",source)),header=TRUE)
    cell.source <- rep(source,nrow(pheno))

    clonotypes <- as.character(pheno$clonotype)

    pheno.all <- clonotypes[!is.na(clonotypes)]
    names.presort <- unique(pheno.all)
    length(names.presort)

    oo <- order(as.numeric(sub(".*[.]C","",names.presort)))
    names.all <- names.presort[oo]

    table.1 <- table(pheno$clonotype)
    names.1 <- names(table.1)[table.1 > 0]
    clonotype.counts <- numeric(length(names.all))
    names(clonotype.counts) <- names.all
    clonotype.counts[match(names.1,names.all)] <- table.1[table.1 > 0]

    # table.i <- table(pheno.2$clonotype)
    # names.i <- names(table.i)[table.i > 0]
    # clonotype.counts2 <- numeric(length(names.all))
    # names(clonotype.counts2) <- names.all
    # clonotype.counts2[match(names.i,names.all)] <- table.i[table.i > 0]

    # clonotype.bcounts <- rep(0, length(clonotype.tcounts))

    # clonotype.residency <- character(length(names.all))
    # names(clonotype.residency) <- names.all
    # clonotype.residency[clonotype.counts2 == 1 & clonotype.counts1 == 0] <- paste(source2, "singletons")
    # clonotype.residency[clonotype.counts1 == 1 & clonotype.counts2 == 0] <- paste(source1, "singletons")
    # clonotype.residency[clonotype.counts2 > 1 & clonotype.counts1 == 0] <- paste(source2, "multiplets")
    # clonotype.residency[clonotype.counts1 > 1 & clonotype.counts2 == 0] <- paste(source1, "multiplets")
    # clonotype.residency[clonotype.counts1 >= 1 & clonotype.counts2 >= 1] <- "Dual resident"
    # clonotype.residency <- factor(clonotype.residency, levels=c(
    #     paste(source2, "singletons"),
    #     paste(source1, "singletons"),
    #     paste(source2, "multiplets"),
    #     paste(source1, "multiplets"),
    #     "Dual resident"
    # ), labels=c(tolower(source2), tolower(source1), toupper(source2), toupper(source1), "Dual"))

    # cell.residency <- clonotype.residency[clonotypes]
    # sources = c(source1, source2)
    clonotype.residency = NULL
    cell.residency = NULL
    clonotype.counts1 = clonotype.counts
    clonotype.counts2 = NULL
    sources = source
    save(clonotypes, clonotype.counts1, clonotype.counts2, sources,
         clonotype.residency, cell.residency, cell.source,
         file=file.path(outdir, paste0(sample,".clonotype.counts.RData")))
}

compute.counts.2 <- function (sample, source1, source2, outdir) {

    pheno.1 <- read.table(file.path(outdir, paste0(sample,".pheno.",source1)),header=TRUE)
    pheno.2 <- read.table(file.path(outdir, paste0(sample,".pheno.",source2)),header=TRUE)
    cell.source <- c(rep(source1,nrow(pheno.1)), rep(source2,nrow(pheno.2)))

    clonotypes <- c(as.character(pheno.1$clonotype),as.character(pheno.2$clonotype)) # Includes NA

    pheno.all <- clonotypes[!is.na(clonotypes)]
    names.presort <- unique(pheno.all)
    length(names.presort)

    oo <- order(as.numeric(sub(".*[.]C","",names.presort)))
    names.all <- names.presort[oo]

    table.i <- table(pheno.1$clonotype)
    names.i <- names(table.i)[table.i > 0]
    clonotype.counts1 <- numeric(length(names.all))
    names(clonotype.counts1) <- names.all
    clonotype.counts1[match(names.i,names.all)] <- table.i[table.i > 0]

    table.i <- table(pheno.2$clonotype)
    names.i <- names(table.i)[table.i > 0]
    clonotype.counts2 <- numeric(length(names.all))
    names(clonotype.counts2) <- names.all
    clonotype.counts2[match(names.i,names.all)] <- table.i[table.i > 0]

    # clonotype.bcounts <- rep(0, length(clonotype.tcounts))

    clonotype.residency <- character(length(names.all))
    names(clonotype.residency) <- names.all
    clonotype.residency[clonotype.counts2 == 1 & clonotype.counts1 == 0] <- paste(source2, "singletons")
    clonotype.residency[clonotype.counts1 == 1 & clonotype.counts2 == 0] <- paste(source1, "singletons")
    clonotype.residency[clonotype.counts2 > 1 & clonotype.counts1 == 0] <- paste(source2, "multiplets")
    clonotype.residency[clonotype.counts1 > 1 & clonotype.counts2 == 0] <- paste(source1, "multiplets")
    clonotype.residency[clonotype.counts1 >= 1 & clonotype.counts2 >= 1] <- "Dual resident"
    clonotype.residency <- factor(clonotype.residency, levels=c(
        paste(source2, "singletons"),
        paste(source1, "singletons"),
        paste(source2, "multiplets"),
        paste(source1, "multiplets"),
        "Dual resident"
    ), labels=c(tolower(source2), tolower(source1), toupper(source2), toupper(source1), "Dual"))

    cell.residency <- clonotype.residency[clonotypes]
    sources = c(source1, source2)
    save(clonotypes, clonotype.counts1, clonotype.counts2, sources,
         clonotype.residency, cell.residency, cell.source,
         file=file.path(outdir, paste0(sample,".clonotype.counts.RData")))
}


#sample <- "renal1.tnb"
#tumor.phenofile <- "../TCR-counts-0/renal1.tnb.pheno1"
#normal.phenofile <- "../TCR-counts-0/renal1.tnb.pheno2"
#blood.phenofile <- "../TCR-counts-0/renal1.tnb.pheno3"

# compute.counts.3 <- function (sample) {
#     pheno.1 <- read.table(paste0(sample,".pheno1"),header=TRUE)
#     pheno.2 <- read.table(paste0(sample,".pheno2"),header=TRUE)
#     pheno.3 <- read.table(paste0(sample,".pheno3"),header=TRUE)
#     cell.source <- c(rep("BM",nrow(pheno.1)), rep("WBC",nrow(pheno.2)), rep("Blood",nrow(pheno.3)))

#     clonotypes <- c(as.character(pheno.1$clonotype),as.character(pheno.2$clonotype),as.character(pheno.3$clonotype)) # Includes NA

#     pheno.all <- clonotypes[!is.na(clonotypes)]
#     names.presort <- unique(pheno.all)
#     length(names.presort)

#     oo <- order(as.numeric(sub(".*[.]C","",names.presort)))
#     names.all <- names.presort[oo]

#     table.i <- table(pheno.1$clonotype)
#     names.i <- names(table.i)[table.i > 0]
#     clonotype.tcounts <- numeric(length(names.all))
#     names(clonotype.tcounts) <- names.all
#     clonotype.tcounts[match(names.i,names.all)] <- table.i[table.i > 0]

#     table.i <- table(pheno.2$clonotype)
#     names.i <- names(table.i)[table.i > 0]
#     clonotype.ncounts <- numeric(length(names.all))
#     names(clonotype.ncounts) <- names.all
#     clonotype.ncounts[match(names.i,names.all)] <- table.i[table.i > 0]

#     table.i <- table(pheno.3$clonotype)
#     names.i <- names(table.i)[table.i > 0]
#     clonotype.bcounts <- numeric(length(names.all))
#     names(clonotype.bcounts) <- names.all
#     clonotype.bcounts[match(names.i,names.all)] <- table.i[table.i > 0]

#     clonotype.residency <- character(length(names.all))
#     names(clonotype.residency) <- names.all
#     clonotype.residency[clonotype.ncounts == 0 & clonotype.tcounts == 0 & clonotype.bcounts == 1] <- "Blood singletons"
#     clonotype.residency[clonotype.ncounts == 0 & clonotype.tcounts == 0 & clonotype.bcounts > 1] <- "Blood multiplets"
#     clonotype.residency[clonotype.ncounts == 1 & clonotype.tcounts == 0] <- "WBC singletons"
#     clonotype.residency[clonotype.tcounts == 1 & clonotype.ncounts == 0] <- "BM singletons"
#     clonotype.residency[clonotype.ncounts > 1 & clonotype.tcounts == 0] <- "WBC multiplets"
#     clonotype.residency[clonotype.tcounts > 1 & clonotype.ncounts == 0] <- "BM multiplets"
#     clonotype.residency[clonotype.tcounts >= 1 & clonotype.ncounts >= 1] <- "Dual resident"
#     clonotype.residency <- factor(clonotype.residency, levels=c("WBC singletons","WBC multiplets","Dual resident","BM multiplets","BM singletons","Blood singletons", "Blood multiplets"), labels=c("n","N","D","T","t","b","B"))

#     cell.residency <- clonotype.residency[clonotypes]

#     save(clonotypes, clonotype.tcounts, clonotype.ncounts, clonotype.bcounts,
#          clonotype.residency, cell.residency, cell.source,
#          file=paste0(sample,".clonotype.counts.RData"))
# }

args <- commandArgs(trailingOnly=TRUE)
print(args)
patient = args[1]
outdir = args[2]
sources = args[3:length(args)]

if (length(sources) == 2) {
    compute.counts.2(
        patient,
        sources[1],
        sources[2],
        outdir
        # args[length(args)-3], # patient
        # args[length(args)-2], # source1
        # args[length(args)-1], # source2
        # args[length(args)] # outdir
    )
} else {
    compute.counts.1(
        patient,
        sources[1],
        outdir
    )
}
