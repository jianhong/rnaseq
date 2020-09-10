#!/usr/bin/env Rscript

# Command line argument processing
args = commandArgs(trailingOnly=TRUE)
rmdpath <- args[1]
results <- args[2]

rmd <- readLines(rmdpath)
id <- which(rmd=="##insertBlockHere##")
blocks <- paste0("resultsFolder <- '", results, "'")
rmd <- c(rmd[seq.int(id-1)], blocks, rmd[-seq.int(id)])
writeLines(rmd, "index.rmd")
rmarkdown::render(input = "index.rmd")

