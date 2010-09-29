
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
  cat("Cohesive block structure:\n")
  myb <- blocks(x)
  ch <- cohesion(x)
  pp <- parent(x)
  si <- sapply(myb, length)

  cs <- 3 + 2 + nchar(length(x)) +
    max(shortest.paths(hierarchy(x), mode="out", v=1)) * 3
  
  .plot <- function(b, ind="") {
    if (b!=1) {
      he <- format(paste(sep="", ind, "'- B-", b), width=cs)
      ind <- paste("  ", ind)
    } else {
      he <- format(paste(sep="", "B-", b), width=cs)
    }
    cat(sep="", he,
        "c ", format(ch[b], width=nchar(max(ch)), justify="right"),
        ", n ", format(si[b], width=nchar(x$vcount), justify="right"))

    if (x$vcount <= options("width")$width-40 && b != 1) {
      o <- rep(".", x$vcount)
      o[ myb[[b]] ] <- "o"
      oo <- character()
      for (i in 1:floor(x$vcount/10)) {
        oo <- c(oo, o[((i-1)*10+1):(i*10)], " ")
      }
      if (x$vcount %% 10) { oo <- c(oo, o[(i*10+1):length(o)]) }
      cat("  ", paste(oo, collapse=""), "\n")
    } else {
      cat("\n")
    }

    wc <- which(pp==b)
    sapply(wc, .plot, ind=ind)
  }
  if (length(x) >0) .plot(1) else cat("No cohesive blocks found.")
  
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
                            root=1), ...) {
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
    thisb[ myb[[b]] ] <- 1
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
    thisb[ myb[[b]] ] <- 1
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
    res[ myb[[b]] ] <- coh[b]
  }
  res
}

#########################################################
## Various designs to print the cohesive blocks

## Cohesive block structure:
## B-1          c. 1, n. 34
## '- B-2       c. 2, n. 28    1,2,3,4,8,9,10,13,14,15,16,18,19,20,21,22,
##    |                        23,24,25,26,27,28,29,30,31,32,33,34
##    '- B-4    c. 4, n.  5    1,2,3,4,8
##    '- B-5    c. 3, n.  7    1,2,3,9,31,33,34
##    '- B-7    c. 4, n.  5    1,2,3,4,14
##    '- B-8    c. 3, n. 10    3,24,25,26,28,29,30,32,33,34
## '- B-3       c. 2, n.  6    1,5,6,7,11,17
##    '- B-6    c. 3, n.  5    1,5,6,7,11

## Cohesive block structure:
## B-1          c. 1, n. 23
## '- B-2       c. 2, n. 14    1,2,3,4,5,6,7,8,12,19,20,21,22,23
##    '- B-4    c. 5, n.  7    1,2,3,4,5,6,7
## '- B-3       c. 2, n. 10    7,9,10,11,13,14,15,16,17,18
##    '- B-5    c. 3, n.  4    7,9,10,11

## #########################################################

## Cohesive block structure:
## B-1        c 1, n 34    
## '- B-2     c 2, n 28    oooo...ooo ..oooo.ooo oooooooooo oooo
##    '- B-4  c 4, n  5    oooo...o.. .......... .......... ....
##    '- B-5  c 3, n  7    ooo.....o. .......... .......... o.oo
##    '- B-7  c 4, n  5    oooo...... ...o...... .......... ....
##    '- B-8  c 3, n 10    ..o....... .......... ...ooo.ooo .ooo
## '- B-3     c 2, n  6    o...ooo... o.....o... .......... ....
##    '- B-6  c 3, n  5    o...ooo... o......... .......... ....

## Cohesive block structure:
## B-1        c 1, n 23    oooooooooo oooooooooo ooo
## '- B-2     c 2, n 14    oooooooo.. .o......oo ooo
##    '- B-4  c 5, n  7    ooooooo... .......... ...
## '- B-3     c 2, n 10    ......o.oo o.oooooo.. ...
##    '- B-5  c 3, n  4    ......o.oo o......... ...

## #########################################################

## Cohesive block structure:
## B-1          c. 1, n. 34
## '- B-2       c. 2, n. 28     1, 2, 3, 4, 8, 9,10,13,14,15,16,18,19,20,21,
##    |                        22,23,24,25,26,27,28,29,30,31,32,33,34
##    '- B-4    c. 4, n.  5     1, 2, 3, 4, 8
##    '- B-5    c. 3, n.  7     1, 2, 3, 9,31,33,34
##    '- B-7    c. 4, n.  5     1, 2, 3, 4,14
##    '- B-8    c. 3, n. 10     3,24,25,26,28,29,30,32,33,34
## '- B-3       c. 2, n.  6     1, 5, 6, 7,11,17
##    '- B-6    c. 3, n.  5     1, 5, 6, 7,11

## Cohesive block structure:
## B-1          c. 1, n. 23
## '- B-2       c. 2, n. 14     1, 2, 3, 4, 5, 6, 7, 8,12,19,20,21,22,23
##    '- B-4    c. 5, n.  7     1, 2, 3, 4, 5, 6, 7
## '- B-3       c. 2, n. 10     7, 9,10,11,13,14,15,16,17,18
##    '- B-5    c. 3, n.  4     7, 9,10,11

## #########################################################

## Cohesive block structure:
## B-1          c. 1, n. 34
## '- B-2       c. 2, n. 28    1-4, 8-10, 13-16, 18-34
##    '- B-4    c. 4, n.  5    1-4, 8
##    '- B-5    c. 3, n.  7    1-3, 9, 31, 33-34
##    '- B-7    c. 4, n.  5    1-4, 14
##    '- B-8    c. 3, n. 10    3, 24-26, 28-30, 32-34
## '- B-3       c. 2, n.  6    1, 5-7, 11, 17
##    '- B-6    c. 3, n.  5    1, 5-7, 11

## Cohesive block structure:
## B-1          c. 1, n. 23
## '- B-2       c. 2, n. 14    1-8, 12, 19-23
##    '- B-4    c. 5, n.  7    1-7
## '- B-3       c. 2, n. 10    7, 9-11, 13-18
##    '- B-5    c. 3, n.  4    7, 9-11

## ##########################################################

## Cohesive block structure:
## B-1        c. 1, n. 34    
## |- B-2     c. 2, n. 28  [ 1] oooo...ooo ..oooo.ooo 
## |  |                    [21] oooooooooo oooo
## |  |- B-4  c. 4, n.  5  [ 1] oooo...o.. .......... 
## |  |                    [21] .......... ....
## |  |- B-5  c. 3, n.  7  [ 1] ooo.....o. .......... 
## |  |                    [21] .......... o.oo
## |  |- B-7  c. 4, n.  5  [ 1] oooo...... ...o...... 
## |  |                    [21] .......... ....
## |  |- B-8  c. 3, n. 10  [ 1] ..o....... .......... 
## |                       [21] ...ooo.ooo .ooo
## '- B-3     c. 2, n.  6  [ 1] o...ooo... o.....o... 
##    |                    [21] .......... ....
##    '- B-6  c. 3, n.  5  [ 1] o...ooo... o......... 
##                         [21] .......... ....

## Cohesive block structure:
## B-1          c. 1, n. 23  [ 1] oooooooooo oooooooooo 
## |                         [21] ooo
## |- B-2       c. 2, n. 14  [ 1] oooooooo.. .o......oo 
## |  |                      [21] ooo
## |  '- B-4    c. 5, n.  7  [ 1] ooooooo... .......... 
## |                         [21] ...
## '- B-3       c. 2, n. 10  [ 1] ......o.oo o.oooooo.. 
##    |                      [21] ...
##    '- B-5    c. 3, n.  4  [ 1] ......o.oo o.........
##                           [21] ...
