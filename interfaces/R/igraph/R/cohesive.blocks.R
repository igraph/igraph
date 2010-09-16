
#   IGraph R package
#   Copyright (C) 2010  Gabor Csardi <csardi.gabor@gmail.com>
#   Rue de l'Industrie 5, 1005 Lausanne, Switzerland
#   
#   This program is free software; you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation; either version 2 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#   
#   You should have received a copy of the GNU General Public License
#   along with this program; if not, write to the Free Software
#   Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA
#   02110-1301 USA
#
###################################################################

cohesive.blocks <- function(graph, labels=TRUE) {
  # Argument checks
  if (!is.igraph(graph)) { stop("Not a graph object") }

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  # Function call
  res <- .Call("R_igraph_cohesive_blocks", graph,
        PACKAGE="igraph")
  class(res) <- "cohesiveBlocks"
  if (labels && "name" %in% list.vertex.attributes(graph)) {
    res$labels <- V(graph)$name
  }
  res$vcount <- vcount(graph)
  res
}

length.cohesiveBlocks <- function(x) {
  length(x$blocks)
}

blocks <- function(blocks) {
  blocks$blocks
}

blockGraphs <- function(blocks, graph) {
  lapply(blocks(blocks), induced.subgraph, graph=graph)
}

cohesion <- function(blocks) {
  blocks$cohesion
}

hierarchy <- function(blocks) {
  blocks$blockTree
}

parent <- function(blocks) {
  blocks$parent
}

print.cohesiveBlocks <- function(x, ...) {
  cat("Structurally cohesive block structure, with",
      length(blocks(x)), "blocks.\n\n")
  cat("Hierarchy:\n")
  print(hierarchy(x))
  cat("\nBlocks and their cohesion:\n")
  myb <- blocks(x)
  lapply(seq_along(myb), function(b) {
    cat(sep="", "[[", b-1, "]], cohesion: ", cohesion(x)[b],
        ", parent: ", parent(x)[b], "\n")
    if (!is.null(x$labels)) {
      print(x$labels[ myb[[b]]+1 ])
    } else {
      print(myb[[b]])
    }
  })
  invisible(x)
}

summary.cohesiveBlocks <- function(x, ...) {
  cat("Structurally cohesive block structure, with",
      length(blocks(x)), "blocks.\n")
  invisible(x)
}

plot.cohesiveBlocks <- function(x, y,
                                colbar=rainbow(max(cohesion(x))+1),
                                col=colbar[maxcohesion(x)+1],
                                mark.groups=blocks(x)[-1],
                                ...) {
  plot(y, mark.groups=mark.groups,
       vertex.color=col, ...)
}

plotHierarchy <- function(blocks,
                          layout=layout.reingold.tilford(hierarchy(blocks),
                            root=0), ...) {
  plot(hierarchy(blocks), layout=layout, ...)
}

exportPajek.cohesiveblocks.pf <- function(blocks, graph, file) {

  closeit <- FALSE
  if (is.character(file)) {
    file <- file(file, open = "w+b")
    closeit <- TRUE
  }
  if (!isOpen(file)) {
    file <- open(file)
    closeit <- TRUE
  }

  ## The original graph
  cat(file=file, sep="", "*Network cohesive_blocks_input.net\r\n")
  write.graph(graph, file=file, format="pajek")

  ## The hierarchy graph
  cat(file=file, sep="", "\r\n*Network hierarchy.net\r\n")
  write.graph(hierarchy(blocks), file=file, format="pajek")

  ## The blocks
  myb <- blocks(blocks)
  for (b in seq_along(myb)) {
    thisb <- rep(0, vcount(graph))
    thisb[ myb[[b]]+1 ] <- 1
    cat(file=file, sep="", "\r\n*Partition block_", b, ".clu\r\n",
        "*Vertices ", vcount(graph), "\r\n   ")    
    cat(thisb, sep="\r\n   ", file=file)
  }
  
  if (closeit) {
    close(file)
  }
  invisible(NULL)
}

exportPajek.cohesiveblocks.nopf <- function(blocks, graph, file) {

  ## The original graph
  write.graph(graph, file=paste(sep="", file, ".net"), format="pajek")

  ## The hierarchy graph
  write.graph(hierarchy(blocks), file=paste(sep="", file, "_hierarchy.net"),
              format="pajek")

  ## The blocks
  myb <- blocks(blocks)
  for (b in seq_along(myb)) {
    thisb <- rep(0, vcount(graph))
    thisb[ myb[[b]]+1 ] <- 1
    cat(file=paste(sep="", file, "_block_", b, ".clu"), sep="\r\n",
        paste("*Vertices", vcount(graph)), thisb)
  }
  
  invisible(NULL)
}

exportPajek <- function(blocks, graph, file,
                        project.file=TRUE) {
  
  if (!project.file && !is.character(file)) {
    stop(paste("`file' must be a filename (without extension) when writing",
               "to separate files"))
  }

  if (project.file) {
    return(exportPajek.cohesiveblocks.pf(blocks, graph, file))
  } else {
    return(exportPajek.cohesiveblocks.nopf(blocks, graph, file))
  }
}

maxcohesion <- function(blocks) {
  res <- numeric(blocks$vcount)
  myb <- blocks(blocks)
  coh <- cohesion(blocks)
  oo <- order(coh)
  myb <- myb[oo]
  coh <- coh[oo]
  for (b in seq_along(myb)) {
    res[ myb[[b]]+1 ] <- coh[b]
  }
  res
}
