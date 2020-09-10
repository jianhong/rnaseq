#!/usr/bin/env Rscript

# Command line argument processing
args = commandArgs(trailingOnly=TRUE)
if (length(args) < 3) {
  stop("Usage: DESeq2FromFeatureCounts.r <species> <design table> <featureCounts table file>", call.=FALSE)
}
species <- args[1]
design <- args[2]
counts <- args[3]

# Debug messages (stderr)
message("Input species  (Arg 1): ", species)
message("Input design   (Arg 2): ", design)
message("Input counts   (Arg 3): ", counts)

if(file.exists(design) && file.exists(counts)){
  # Load / install packages
  if(!require("BiocManager")){
    install.packages("BiocManager", suppressUpdates=TRUE)
  }
  
  pkgs <- BiocManager::available("^TxDb")
  TxDb <- pkgs[grepl(paste0("^TxDb.*?", species), pkgs)]
  if(any(grepl("knownGene", TxDb))){
    TxDb <- TxDb[grepl("knownGene", TxDb)][1]
  }
  TxDb
  askPkg <- function(pkg){
    if (!require(pkg, character.only = TRUE)){
      install(pkg, suppressUpdates=TRUE)
      library(pkg, character.only = TRUE)
    }
  }
  shortStrs <- function(strs, len=60){
    if(length(strs)==0) return(strs)
    strs <- as.character(strs)
    shortStr <- function(str, len=60){
      stopifnot(length(str)==1)
      stopifnot(is.character(str))
      if(nchar(str)<=len) return(str)
      strs <- strsplit(str, " ")[[1]]
      nc <- nchar(strs)
      nclast <- nc[length(nc)] + 3
      paste0(substring(str, first = 1, last = len-nclast), "...", strs[length(strs)])
    }
    strs <- sapply(strs, shortStr, len=len)
    make.unique(strs)
  }
  orgShort <- function(organism){
    x <- strsplit(tolower(organism), " ")[[1]]
    paste0(substr(x[1], start = 1, stop = 1), x[2])
  }
  for(pkg in c("ChIPpeakAnno", "DESeq2", "WriteXLS", "clusterProfiler", "scales",
               "EnhancedVolcano", TxDb)){
    askPkg(pkg)
  }
  
  txdb <- get(TxDb)
  organism <- organism(txdb)
  (org <- ChIPpeakAnno::egOrgMap(organism))
  askPkg(org)
  
  cts <- read.delim(counts, comment.char = "#", stringsAsFactors = FALSE)
  anno <- cts[, !grepl("bam$", colnames(cts))]
  cts <- cts[, grepl("bam$", colnames(cts))]
  rownames(cts) <- anno$Geneid
  
  samples <- read.delim(design)
  stopifnot(rownames(samples)==colnames(cts))
  cn <- colnames(samples)
  uniqueCn <- apply(samples, 2, FUN=function(.ele){
    length(unique(.ele))>1
  })
  cn <- cn[uniqueCn]
  cn <- cn[!cn %in% c("sampleName")]
  fl <- as.formula(paste("~", cn, collapse = " "))
  dds <- DESeqDataSetFromMatrix(countData = cts,
                                colData = samples,
                                design = fl)
  if("techRep" %in% cn){
    dds <- collapseReplicates(dds, groupby = dds$techRep)
  }
  dds <- DESeq(dds)
  conditions <- unique(as.character(samples$condition))
  contrasts <- combn(conditions, 2, simplify = FALSE)
  names(contrasts) <- sapply(contrasts, paste, collapse=".vs.")
  
  for(.id in seq_along(contrasts)){
    pf <- names(contrasts)[.id]
    contr <- contrasts[[.id]]
    contr <- contr[order(grepl("control|contr|ctl|day0|minus|neg|sham|WT",
                               contr, ignore.case = TRUE))]
    dds1 <- dds
    res.raw <- results(dds1, contrast = c("condition", contr))
    res <- lfcShrink(dds1, contrast = c("condition", contr), res = res.raw, type="normal")
    
    counts <- counts(dds1)
    rld <- counts(dds1, normalized=TRUE)
    
    genes <- exonsBy(txdb, by="gene")
    rr <- rowRanges(dds1)
    dbExt <- ifelse(grepl("^ENS", names(rr)), "ENSEMBL2EG", "ALIAS2EG")
    dbExt <- sort(table(dbExt), decreasing=TRUE)
    dbExt <- names(dbExt)[1]
    eg <- ChIPpeakAnno::xget(names(rr), get(sub(".db", dbExt, org)), output = "first")
    eg[is.na(eg)] <- "NA"
    rr.x <- genes[eg[eg %in% names(genes)]]
    names(rr.x) <- names(rr[eg %in% names(genes)])
    rr[eg %in% names(genes)] <- rr.x
    rowRanges(dds1) <- rr
    fpkm <- fpkm(dds1)
    
    colnames(counts) <- paste0("counts.", colnames(counts))
    colnames(fpkm) <- paste0("fpkm.", colnames(fpkm))
    colnames(rld) <- paste0("nromlized.counts.", colnames(rld))
    
    stopifnot(identical(rownames(counts), rownames(res)))
    stopifnot(identical(rownames(fpkm), rownames(res)))
    stopifnot(identical(rownames(rld), rownames(res)))
    stopifnot(identical(rownames(res.raw), rownames(res)))
    
    data <- cbind(gene=rownames(res),
                  log2FoldChangeWithoutShrink=res.raw$log2FoldChange,
                  stat=res.raw$stat,
                  as.data.frame(res), counts, fpkm, rld)
    WriteXLS(data, paste0(pf, ".DESeq2.featureCounts.diff.xls"))
    metadata <- as.data.frame(res@elementMetadata)
    WriteXLS(metadata, paste0(pf, ".DESeq2.featureCounts.metadata.xls"))
    pdf(paste0(pf, ".pvalue.dist.pdf"))
    hist(data$pvalue, breaks = 50)
    dev.off()
    pdf(paste0(pf, ".MAplot.pdf"))
    plotMA(res)
    dev.off()
    data.s <- data[!is.na(data$padj), ]
    data.s <- data.s[data.s$padj < 0.05, ]
    WriteXLS(data.s, paste0(pf, ".DESeq2.featureCounts.diff.fdr.0.05.xls"))
    data.s <- data.s[data.s$padj < 0.05 & abs(data.s$log2FoldChangeWithoutShrink)>1, ]
    WriteXLS(data.s, paste0(pf, ".DESeq2.featureCounts.diff.fdr.0.05.lfc.1.xls"))
    
    pdf(paste0(pf, ".DESeq2.featureCounts.volcanonPlot.pdf"), width=9, height=6)
    EnhancedVolcano(data, lab = data$gene, x='log2FoldChange', y='padj',
                    title = metadata[4, 2],
                    titleLabSize = 12,
                    subtitle = NULL,
                    pCutoff = 0.05,
                    ylab = bquote(~-Log[10]~italic(adjP)),
                    FCcutoff = 1,
                    pointSize = 1.0,
                    labSize = 1.0,
                    legendPosition = 'right')
    dev.off()
    
    gene.df <- bitr(data.s$gene,
                    fromType=ifelse(grepl("^ENS", as.character(data.s$gene)[1]),
                                    "ENSEMBL", "SYMBOL"),
                    toType = "ENTREZID", OrgDb = org)
    if(nrow(gene.df)>1){
      ego <- sapply(c("BP", "MF", "CC"), function(.onto){
        enrichGO(  gene = gene.df$ENTREZID,
                   OrgDb = org,
                   ont = .onto,
                   readable = TRUE
        )
      })
      null <- mapply(ego, names(ego), FUN=function(.ele, .name){
        write.csv(.ele, paste0(pf, paste0("GO.", .name, ".enrichment.csv")))
        
        .ele <- as.data.frame(.ele)
        if(nrow(.ele)>1){
          .ele$qvalue <- -log10(.ele$p.adjust)
          plotdata <- .ele[!is.na(.ele$qvalue), c("Description", "qvalue", "Count")]
          if(nrow(plotdata)>20) plotdata <- plotdata[1:20, ]
          plotdata$Description <- shortStrs(plotdata$Description)
          ggplot(plotdata, aes(x=reorder(Description, -qvalue), y=qvalue, fill=Count, label=Count)) +
            scale_fill_gradient2(low = muted("blue"), high = muted("red"), oob = scales::squish) +
            geom_bar(stat="identity") + scale_y_continuous(expand = expand_scale(mult = c(0, .1))) +
            geom_text(vjust=-.1) +
            xlab("") + ylab("-log10(p-value)") +
            theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))
          ggsave(paste0(pf, paste(".GO.", .name, ".enrichment.top.pdf", sep = ".")), width = 6, height = 6)
        }
      })
      
      kk <- enrichKEGG(gene = gene.df$ENTREZID, 
                       organism = substr(orgShort(organism), start = 1, stop = 3))
      kk <- as.data.frame(kk)
      eid <- strsplit(kk$geneID, "\\/")
      symbol <- lapply(eid, function(.ele) gene.df[match(.ele, gene.df$ENTREZID), "SYMBOL"])
      symbol <- sapply(symbol, paste, collapse="/")
      kk$geneSYM <- symbol
      write.csv(kk, file.path(pf, "KEGGenrichment.csv"))
    }
    tryCatch({
      if(orgShort(organism)!="hsapiens"){
        askPkg("biomaRt")
        mart <- useMart('ensembl', dataset = paste0(orgShort(organism), "_gene_ensembl"))
        homologs <- getBM(attributes = c("ensembl_gene_id", 
                                         "hsapiens_homolog_ensembl_gene", 
                                         "hsapiens_homolog_associated_gene_name"),
                          filters = "ensembl_gene_id",
                          values = data$gene,
                          mart = mart)
        homologs <- homologs[match(data$gene, homologs$ensembl_gene_id), 
                             c("hsapiens_homolog_ensembl_gene", "hsapiens_homolog_associated_gene_name")]
        data <- cbind(data, homologs)
        
        rnk <- data[order(data$stat, decreasing = TRUE), ]
        rnk <- rnk[!is.na(rnk$pvalue), ]
        rnk <- rnk[!is.na(rnk$hsapiens_homolog_associated_gene_name), ]
        rnk <- rnk[rnk$hsapiens_homolog_associated_gene_name!="", ]
        rnk <- rnk[, c("hsapiens_homolog_associated_gene_name", "stat")]
        write.table(rnk, paste0(pf, ".gsea.rnk"), row.names = FALSE, col.names = FALSE, quote = FALSE, sep="\t")
      }
    }, error = function(.e) message(.e))
  }
  
}


