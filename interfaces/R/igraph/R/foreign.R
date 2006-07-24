
#   IGraph R package
#   Copyright (C) 2005  Gabor Csardi <csardi@rmki.kfki.hu>
#   MTA RMKI, Konkoly-Thege Miklos st. 29-33, Budapest 1121, Hungary
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

###################################################################
# Reading foreign file formats
###################################################################

read.graph.toraw <- function(filename) {
  if (is.character(filename)) {
    filename <- file(filename)
  }
  if (!isOpen(filename)) {
    open(filename)
  }

  tmpbufsize <- 20000
  buffer <- tmpbuffer <- readChar(filename, tmpbufsize)
  while (nchar(tmpbuffer) == tmpbufsize) {
    tmpbuffer <- readChar(filename, tmpbufsize)
    buffer <- paste(sep="", buffer, tmpbuffer)
  }
  close(filename)
  rm(tmpbuffer)
  
  charToRaw(buffer)  
}

write.graph.fromraw <- function(buffer, file) {

  closeit <- FALSE
  if (is.character(file)) {
    file <- file(file, open="w+b")
    closeit <- TRUE
  }
  
  if (!isOpen(file)) {
    file <- open(file)
    closeit <- TRUE
  }

  writeBin(buffer, file)

  if (closeit) {
    close(file)
  }

  invisible(NULL)
}

read.graph <- function(file, format="edgelist", ...) {

  if (igraph.i.have.fmemopen) {
    file <- read.graph.toraw(file)
  } else {
    if (!is.character(file) || length(grep("://", file, fixed=TRUE))>0) {
      buffer <- read.graph.toraw(file)
      file <- tempfile()
      write.graph.fromraw(buffer, file)
    }
  }
  
  res <- switch(format,
                "pajek"=read.graph.pajek(file, ...),
                "ncol"=read.graph.ncol(file, ...),
                "edgelist"=read.graph.edgelist(file, ...),
                "lgl"=read.graph.lgl(file, ...),
                stop(paste("Unknown file format:",format))
                )
  res
}

write.graph <- function(graph, file, format="edgelist", ...) {

  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }
  if (!igraph.i.have.open.memstream) {
    if (!is.character(file) || length(grep("://", file, fixed=TRUE))>0) {
      tmpfile <- TRUE
      origfile <- file
      file <- tempfile()
    } else {
      tmpfile <- FALSE
    }
  }
  
  res <- switch(format,
                "pajek"=write.graph.pajek(graph, file, ...),
                "edgelist"=write.graph.edgelist(graph, file, ...),
                "ncol"=write.graph.ncol(graph, file, ...),
                "lgl"=write.graph.lgl(graph, file, ...),
                stop(paste("Unknown file format:",format))
                )

  if (igraph.i.have.open.memstream) {
    write.graph.fromraw(res, file)
  } else {
    if (tmpfile) {
      buffer <- read.graph.toraw(file)
      write.graph.fromraw(buffer, origfile)
    }
  }
  
  invisible(res)
}

################################################################
# Plain edge list format, not sorted
################################################################

read.graph.edgelist <- function(file, n=0,
                                directed=TRUE, ...) {

  .Call("R_igraph_read_graph_edgelist", file,
        as.numeric(n), as.logical(directed),
        PACKAGE="igraph")
}

write.graph.edgelist <- function(graph, file, 
                                 ...) {
  
  .Call("R_igraph_write_graph_edgelist", graph, file,
        PACKAGE="igraph")
}

################################################################
# NCOL and LGL formats, quite simple
################################################################

read.graph.ncol <- function(file, predef=character(0), names=TRUE,
                           weights=TRUE, directed=FALSE, ...) {

  .Call("R_igraph_read_graph_ncol", file, as.character(predef),
        as.logical(names), as.logical(weights), as.logical(directed),
        PACKAGE="igraph")
}

write.graph.ncol <- function(graph, file, 
                             names="name", weights="weight", ...) {
  names <- as.character(names)
  weights <- as.character(weights)
  if (length(names)==0 || ! names %in% list.vertex.attributes(graph)) { names <- NULL }
  if (length(weights)==0 || ! weights %in% list.edge.attributes(graph)) { weights <- NULL }
  
  .Call("R_igraph_write_graph_ncol", graph, file,
        names, weights,
        PACKAGE="igraph")
}  

read.graph.lgl <- function(file, names=TRUE,
                           weights=TRUE, ...) {

  .Call("R_igraph_read_graph_lgl", file,
        as.logical(names), as.logical(weights),
        PACKAGE="igraph")
}

write.graph.lgl <- function(graph, file, 
                            names="name", weights="weight",
                            isolates=FALSE, ...) {
  names <- as.character(names)
  weights <- as.character(weights)
  if (length(names)==0 || ! names %in% list.vertex.attributes(graph)) { names <- NULL }
  if (length(weights)==0 || ! weights %in% list.edge.attributes(graph)) { weights <- NULL }
  
  .Call("R_igraph_write_graph_lgl", graph, file,
        names, weights, as.logical(isolates),
        PACKAGE="igraph")
}  

read.graph.pajek <- function(file, ...) {

  .Call("R_igraph_read_graph_pajek", file,
        PACKAGE="igraph")
}

write.graph.pajek <- function(graph, file, ...) {

  .Call("R_igraph_write_graph_pajek", graph, file,
        PACKAGE="igraph")
}
