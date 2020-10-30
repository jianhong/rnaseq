#!/usr/bin/env Rscript

# Command line argument processing
args = commandArgs(trailingOnly=TRUE)
if (length(args) < 3) {
  stop("Usage: DESeq2FromFeatureCounts.r <species> <design table> <featureCounts table file>", call.=FALSE)
}
gtf <- args[1]
design <- args[2]
counts <- args[3]

# Debug messages (stderr)
message("Input gtf  (Arg 1): ", gtf)
message("Input design   (Arg 2): ", design)
message("Input counts   (Arg 3): ", counts)

if(file.exists(design) && file.exists(counts) && file.exists(gtf)){
  # Load / install packages
  if(!require("BiocManager")){
    install.packages("BiocManager", suppressUpdates=TRUE)
  }
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
  for(pkg in c("ChIPpeakAnno", "DESeq2", "WriteXLS", "clusterProfiler", "scales",
               "EnhancedVolcano", "GenomicFeatures", "rtracklayer")){
    askPkg(pkg)
  }
  
  txdb <- makeTxDbFromGFF(gtf)
  gtf <- import(gtf)
  
  cts <- read.delim(counts, comment.char = "#", stringsAsFactors = FALSE)
  anno <- cts[, !grepl("bam$", colnames(cts))]
  cts <- cts[, grepl("bam$", colnames(cts))]
  rownames(cts) <- anno$Geneid
  
  samples <- read.delim(design)
  stopifnot(all(c("condition", "R1") %in% colnames(samples)))
  samples <- samples[match(sub("Aligned.sortedByCoord.*.bam", "", colnames(cts)),
                           sub(".(fastq|fq).gz", "", basename(as.character(samples$R1)))), ]
  rownames(samples) <- colnames(cts)
  dds <- DESeqDataSetFromMatrix(countData = cts,
                                colData = samples,
                                design = ~ condition)
  if("techRep" %in% colnames(samples)){
    dds <- collapseReplicates(dds, groupby = paste(dds$condition, dds$techRep))
  }
  dds <- DESeq(dds)
  conditions <- unique(as.character(samples$condition))
  contrasts <- combn(conditions, 2, simplify = FALSE)
  names(contrasts) <- sapply(contrasts, paste, collapse=".vs.")
  
  for(.id in seq_along(contrasts)){
    pf <- file.path("DESeq2", names(contrasts)[.id])
    dir.create(pf, recursive = TRUE)
    contr <- contrasts[[.id]]
    contr <- contr[order(grepl("control|contr|ctl|day0|minus|neg|sham|WT",
                               contr, ignore.case = TRUE))]
    dds1 <- dds
    res.raw <- results(dds1, contrast = c("condition", contr))
    res <- lfcShrink(dds1, contrast = c("condition", contr), res = res.raw, type="normal")
    
    counts <- counts(dds1)
    rld <- counts(dds1, normalized=TRUE)
    
    fpkm <- NULL
    tryCatch({
      genes <- exonsBy(txdb, by="gene")
      rr <- rowRanges(dds1)
      rowRanges(dds1) <- genes[names(rr)]
      fpkm <- fpkm(dds1)
      }, error=function(.e) message(.e))
    
    colnames(counts) <- paste0("counts.", colnames(counts))
    if(length(fpkm)>0) colnames(fpkm) <- paste0("fpkm.", colnames(fpkm))
    colnames(rld) <- paste0("nromlized.counts.", colnames(rld))
    
    stopifnot(identical(rownames(counts), rownames(res)))
    if(length(fpkm)>0) stopifnot(identical(rownames(fpkm), rownames(res)))
    stopifnot(identical(rownames(rld), rownames(res)))
    stopifnot(identical(rownames(res.raw), rownames(res)))
    
    data <- cbind(gene=rownames(res),
                  log2FoldChangeWithoutShrink=res.raw$log2FoldChange,
                  as.data.frame(res), counts, rld)
    if(length(fpkm)>0) data <- cbind(data, fpkm)
    WriteXLS(data, file.path(pf, paste0(names(contrasts)[.id], ".DESeq2.featureCounts.diff.xls")))
    metadata <- as.data.frame(res@elementMetadata)
    WriteXLS(metadata, file.path(pf, paste0(names(contrasts)[.id], ".DESeq2.featureCounts.metadata.xls")))
    pdf(file.path(pf, paste0(names(contrasts)[.id], ".pvalue.dist.pdf")))
    hist(data$pvalue, breaks = 50)
    dev.off()
    pdf(file.path(pf, paste0(names(contrasts)[.id], ".MAplot.pdf")))
    plotMA(res)
    dev.off()
    data.s <- data[!is.na(data$padj), ]
    data.s <- data.s[data.s$padj < 0.05, ]
    WriteXLS(data.s, file.path(pf, paste0(names(contrasts)[.id], ".DESeq2.featureCounts.diff.fdr.0.05.xls")))
    data.s <- data.s[data.s$padj < 0.05 & abs(data.s$log2FoldChangeWithoutShrink)>1, ]
    WriteXLS(data.s, file.path(pf, paste0(names(contrasts)[.id], ".DESeq2.featureCounts.diff.fdr.0.05.lfc.1.xls")))
    
    pdf(file.path(pf, paste0(names(contrasts)[.id], ".DESeq2.featureCounts.volcanonPlot.pdf")), width=9, height=6)
    EnhancedVolcano(data[, c("gene", "log2FoldChangeWithoutShrink", "padj")], lab = data$gene, x='log2FoldChangeWithoutShrink', y='padj',
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
  }
  
}


