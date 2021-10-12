
library(RcppTOML)
library(tidyr)
library(dplyr)
library(tibble)
library(Seurat)
library(enrichR)

pccdir = "{{ in.pccdir }}"
exprdir = "{{ in.exprdir }}"
samples = "{{ in.samples }}"
outdir = "{{ out.outdir }}"
tclusters = '{{ args.tclusters }}'
cd4cd8clusters = '{{ args.cd4cd8clusters }}'
commoncfg = '{{ args.commoncfg }}'
multipt_samples = '{{ args.multipt_samples }}'

tclusters = parseTOML(tclusters, fromFile=FALSE)
cd4cd8clusters = parseTOML(cd4cd8clusters, fromFile=FALSE)
multipt_samples = parseTOML(multipt_samples, fromFile=FALSE)

idents.tcell.pal = unlist(tclusters$colors)
names(idents.tcell.pal) = unlist(tclusters$names)

dir.create(outdir, showWarnings = FALSE)

load(samples)
tcrdir = "{{ in.tcrdir }}"
load(file.path(pccdir, "global.clonotype.ident.RData"))
# load(file.path(pccdir, "global.counts.RData"))
setEnrichrSite("Enrichr")

gsea_dbs = parseTOML(commoncfg, fromFile=FALSE)$GSEA_DBs
patients = setdiff(names(multipt_samples), "COMPARING")
normal_ident = clonotype.normal.ident[!is.na(clonotype.normal.ident)]
tumor_ident = clonotype.tumor.ident[!is.na(clonotype.tumor.ident)]

is.normal = function(src, controls = c("wbc", "blood", "normal", "control")) {
    tolower(src) %in% controls
}

find_shared_clones = function(presample, postsample) {
    premaster = read.table(
        file.path(tcrdir, paste0(presample, ".master")),
        header=FALSE,
        row.names=NULL,
        sep="\t"
    )
    postmaster = read.table(
        file.path(tcrdir, paste0(postsample, ".master")),
        header=FALSE,
        row.names=NULL,
        sep="\t"
    )
    shared_clones = inner_join(premaster, postmaster, by=c("V2", "V4")) %>%
        select(
            pre_clone=V1.x,
            post_clone=V1.y,
            alpha_consensus=V2,
            beta_consensus=V4
        ) %>%
        mutate(
            ID=paste0("SC_", row_number()),
            .before=1
        )
    shared_clones
}

merge_exprs = function(prefices, pheno) {
    # prefices: BM-1, WBC-1
    # pheno: barcode, ... , clonotype (cols)
    # Merge the exprs of columns from barcode (WBC-13_AAACCTGCAGCTCGCA-1)
    #   to clonotype (MM005-postr.C140)
    # It is a many-to-one mapping
    pheno1 = pheno %>% select(barcode, clonotype) %>% filter(!is.na(clonotype))

    exprs = list()
    i = 1
    for (pref in prefices) {
        # before t():
        # cols: WBC-13_AAACCTGCAGCTCGCA-1, ...
        # rows: Genes
        exprs[[i]] = readRDS(file.path(exprdir, paste0(pref, ".mat.rds"))) %>%
            t() %>%
            as.data.frame() %>%
            rownames_to_column("barcode") %>%
            right_join(pheno1) %>%
            select(-"barcode") %>%
            group_by(clonotype) %>%
            summarise(across(everything(), sum)) %>%
            column_to_rownames("clonotype") %>%
            t()
        i = i + 1
    }

    if (length(exprs) == 1) {
        return (exprs[[1]])
    }

    common_clones = intersect(colnames(exprs[[1]]), colnames(exprs[[2]]))
    cbind(
        exprs[[1]][, setdiff(colnames(exprs[[1]]), common_clones)],
        exprs[[2]][, setdiff(colnames(exprs[[2]]), common_clones)],
        exprs[[1]][, common_clones] + exprs[[2]][, common_clones]
    )
}

find_markers = function(idents, type, pre_or_post, name, cluster, cldir, exprs) {
    # idents: the clone idents
    # type: Expanded or Contracted
    # pre_or_post: Whether find markers in pre sample or post
    # name: Tumor_and_Normal/Tumor/Normal
    # cluster: The cluster name
    # cldir: <outdir>/<patient>/<name>/<cluster>
    # exprs: The expr matrix (rows: Gene, cols: MM005-earlier.C3, ...)
    if (nrow(idents) < 3) { return(NULL) }

    marksdir = file.path(cldir, paste0(type, "-", pre_or_post))
    if (name == "Tumor_and_Normal") {
        # MM005-earlier.C3, ...
        normal_cl_clones = names(normal_ident[normal_ident == cluster])
        tumor_cl_clones = names(tumor_ident[tumor_ident == cluster])
        clones = c(normal_cl_clones, tumor_cl_clones)
    } else if (name == "Tumor") {
        clones = names(tumor_ident[tumor_ident == cluster])
    } else {
        clones = names(normal_ident[normal_ident == cluster])
    }
    clones = intersect(colnames(exprs), unique(clones))

    if (pre_or_post == "Pre") {
        target_clones = intersect(clones, idents$pre_clone)
    } else {
        target_clones = intersect(clones, idents$post_clone)
    }

    if (length(target_clones) == 0) { return (NULL) }

    exprs0 = exprs[, clones, drop=F]
    exprs1 = exprs0[, target_clones, drop=F]
    exprs2 = exprs0[, setdiff(clones, target_clones), drop=F]
    de = tryCatch({
        s1 = CreateSeuratObject(counts=exprs1, min.cells=3, min.features=200)
        s1$group = "ChangedClones"
        s2 = CreateSeuratObject(counts=exprs2, min.cells=3, min.features=200)
        s2$group = 'RestClones'
        merge(s1, s2)
    }, error=function(e) {
        NULL
    })
    if (is.null(de)) { return(NULL) }
    dir.create(marksdir, showWarnings=FALSE)

    de$percent.mt = PercentageFeatureSet(de, pattern = '^MT-')
    de$percent.mt = de$percent.mt[de$percent.mt < 7.5]
    de <- NormalizeData(de, normalization.method = "LogNormalize")
    Idents(de) = "group"

    markers = FindMarkers(object = de, ident.1 = 'ChangedClones')
    sig_markers = markers %>% filter(p_val_adj < 0.05)

    write.table(
        sig_markers,
        file.path(marksdir, 'markers.txt'),
        row.names=T,
        col.names=T,
        sep="\t",
        quote=F
    )

    genes = rownames(sig_markers)
    enriched = enrichr(genes, gsea_dbs)
    for (db in gsea_dbs) {
        outtable = file.path(marksdir, paste0('enrichr_', db, '.txt'))
        outfig = file.path(marksdir, paste0('enrichr_', db, '.png'))

        write.table(enriched[[db]], outtable, col.names=T, row.names=F, sep="\t", quote=F)

        png(outfig, width=1000, height=1000, res=100)
        print(plotEnrich(enriched[[db]], title=db))
        dev.off()
    }
}

handle_kind = function(patient, name, idents) {
    kind_dir = file.path(outdir, patient, name)
    dir.create(kind_dir, showWarnings = FALSE, recursive = TRUE)
    idents = idents %>% as.data.frame() %>% rownames_to_column('Clone')
    clusters = setdiff(colnames(idents), "Clone")
    presample = multipt_samples[[patient]][1]
    postsample = multipt_samples[[patient]][2]
    # convert sample to prefix
    pre_prefices = samples %>%
        filter(Patient == presample) %>%
        distinct(Prefix, .keep_all=T) %>%
        pull(Prefix)
    post_prefices = samples %>%
        filter(Patient == postsample) %>%
        distinct(Prefix, .keep_all=T) %>%
        pull(Prefix)

    pre_pheno_names = unique(unlist(lapply(strsplit(pre_prefices, "-"), function(s) s[1])))
    post_pheno_names = unique(unlist(lapply(strsplit(post_prefices, "-"), function(s) s[1])))

    pre_pheno = NULL
    post_pheno = NULL
    for (pre_pheno_name in pre_pheno_names) {
        phenofile = file.path(tcrdir, paste0(presample, ".pheno.", pre_pheno_name))
        pre_pheno = bind_rows(
            pre_pheno,
            read.table(phenofile, header=TRUE, row.names=NULL, sep="\t", check.names=FALSE)
        )
    }
    for (post_pheno_name in post_pheno_names) {
        phenofile = file.path(tcrdir, paste0(postsample, ".pheno.", post_pheno_name))
        post_pheno = bind_rows(
            post_pheno,
            read.table(phenofile, header=TRUE, row.names=NULL, sep="\t", check.names=FALSE)
        )
    }
    pre_exprs = merge_exprs(pre_prefices, pre_pheno)
    post_exprs = merge_exprs(post_prefices, post_pheno)

    shared_clones = find_shared_clones(presample, postsample)

    idents_pre = idents %>%
        filter(startsWith(Clone, paste0(presample, "."))) %>%
        left_join(shared_clones, by=c("Clone" = "pre_clone"), keep=TRUE)
    idents_pre$ID[is.na(idents_pre$ID)] = paste0("PREUC_", seq_len(sum(is.na(idents_pre$ID))))

    idents_post = idents %>%
        filter(startsWith(Clone, paste0(postsample, "."))) %>%
        left_join(shared_clones, by=c("Clone" = "post_clone"), keep=TRUE)
    idents_post$ID[is.na(idents_post$ID)] = paste0("POSTUC_", seq_len(sum(is.na(idents_post$ID))))

    pt_idents = full_join(
            idents_pre,
            idents_post,
            by=c("ID", "pre_clone", "post_clone"),
            suffix=c("_pre", "_post")
        ) %>%
        select(!starts_with("alpha_consensus_") & !starts_with("beta_consensus_")) %>%
        mutate(
            pre_clone=as.character(pre_clone),
            post_clone=as.character(post_clone),
            pre_clone=if_else(is.na(pre_clone), Clone_pre, pre_clone),
            post_clone=if_else(is.na(post_clone), Clone_post, post_clone)
        ) %>%
        select(!starts_with("Clone_"))

    for (cluster in clusters) {
        cldir = file.path(kind_dir, cluster)
        dir.create(cldir, showWarnings=FALSE)
        pre_cluster = paste0(cluster, "_pre")
        post_cluster = paste0(cluster, "_post")
        cols = c("ID", "pre_clone", "post_clone", pre_cluster, post_cluster)
        cl_idents = pt_idents %>%
            select(cols) %>%
            mutate(diff=pt_idents[[post_cluster]] - pt_idents[[pre_cluster]]) %>%
            filter(!is.na(diff))
        cl_expanded_idents = cl_idents %>%
            filter(diff > 3) %>%
            arrange(desc(diff))
        cl_contracted_idents = cl_idents %>%
            filter(diff < -3) %>%
            arrange(diff)

        write.table(cl_expanded_idents, file.path(cldir, "Expanded.txt"), row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)
        write.table(cl_contracted_idents, file.path(cldir, "Contracted.txt"), row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)

        find_markers(cl_expanded_idents, "Expanded", "Pre", name, cluster, cldir, pre_exprs)
        find_markers(cl_expanded_idents, "Expanded", "Post", name, cluster, cldir, post_exprs)
        find_markers(cl_contracted_idents, "Contracted", "Pre", name, cluster, cldir, pre_exprs)
        find_markers(cl_contracted_idents, "Contracted", "Post", name, cluster, cldir, post_exprs)

    }

    for (cd4cd8 in names(cd4cd8clusters)) {
        cldir = file.path(kind_dir, paste0(cd4cd8, "-groups"))
        dir.create(cldir, showWarnings=FALSE)
        groups = unname(unlist(tclusters$names[cd4cd8clusters[[cd4cd8]]]))

        pre_cluster = paste0(groups, "_pre")
        post_cluster = paste0(groups, "_post")

        cols = c("ID", "pre_clone", "post_clone", pre_cluster, post_cluster)
        cl_idents = pt_idents %>%
            select(cols) %>%
            mutate(
                pre_sum=rowSums(pt_idents[, pre_cluster]),
                post_sum=rowSums(pt_idents[, post_cluster]),
            ) %>%
            mutate(diff=post_sum - pre_sum) %>%
            filter(!is.na(diff))
        cl_expanded_idents = cl_idents %>%
            filter(diff > 3) %>%
            arrange(desc(diff))
        cl_contracted_idents = cl_idents %>%
            filter(diff < -3) %>%
            arrange(diff)

        write.table(cl_expanded_idents, file.path(cldir, "Expanded.txt"), row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)
        write.table(cl_contracted_idents, file.path(cldir, "Contracted.txt"), row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)

        find_markers(cl_expanded_idents, "Expanded", "Pre", name, cd4cd8, cldir, pre_exprs)
        find_markers(cl_expanded_idents, "Expanded", "Post", name, cd4cd8, cldir, post_exprs)
        find_markers(cl_contracted_idents, "Contracted", "Pre", name, cd4cd8, cldir, pre_exprs)
        find_markers(cl_contracted_idents, "Contracted", "Post", name, cd4cd8, cldir, post_exprs)
    }

}


for (patient in patients) {
    handle_kind(patient, "Tumor_and_Normal", clonotype.idents)
    handle_kind(patient, "Tumor", clonotype.tumor.idents)
    handle_kind(patient, "Normal", clonotype.normal.idents)
}
