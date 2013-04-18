#! /usr/bin/env Rscript

library(tools)

rdfiles <- list.files("igraph/man", pattern=".*\\.Rd$", full.names=TRUE)
out <- file("igraph-Ex.R", open="w")
cat("### Load the package\nlibrary(igraph)\n\n", file=out)
sapply(rdfiles, Rd2ex, out=out)
close(out)

