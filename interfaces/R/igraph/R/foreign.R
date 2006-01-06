
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
#   Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
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

  res <- switch(format,
#                "pajek"=read.graph.pajek(file, ...),
                "ncol"=read.graph.ncol(file, ...),
                "edgelist"=read.graph.edgelist(file, ...),
                "lgl"=read.graph.lgl(file, ...),
                stop(paste("Unknown file format:",format))
                )
  res
}

write.graph <- function(graph, file, format="edgelist", ...) {
  
  res <- switch(format,
#                "pajek"=write.graph.pajek(graph, file, ...),
                "edgelist"=write.graph.edgelist(graph, file, ...),
                "ncol"=write.graph.ncol(graph, file, ...),
                stop(paste("Unknown file format:",format))
                )
  invisible(res)
}

################################################################
# Plain edge list format, not sorted
################################################################

read.graph.edgelist <- function(filename, n=0,
                                directed=TRUE, ...) {

  buffer <- read.graph.toraw(filename)
  .Call("R_igraph_read_graph_edgelist", buffer,
        as.numeric(n), as.logical(directed),
        PACKAGE="igraph")
}

write.graph.edgelist <- function(graph, file, 
                                 ...) {

  buffer <- .Call("R_igraph_write_graph_edgelist", graph,
                  PACKAGE="igraph")
  write.graph.fromraw(buffer, file)

  invisible(NULL)
}

################################################################
# NCOL and LGL formats, quite simple
################################################################

read.graph.ncol <- function(filename, names=TRUE,
                           weights=TRUE, ...) {

  buffer <- read.graph.toraw(filename)
  .Call("R_igraph_read_graph_ncol", buffer,
        as.logical(names), as.logical(weights),
        PACKAGE="igraph")
}

write.graph.ncol <- function(graph, file, 
                             names="name", weights="weight", ...) {
  names <- as.character(names)
  weights <- as.character(weights)
  if (length(names)==0 || ! names %in% v.a(graph)) { names <- NULL }
  if (length(weights)==0 || ! weights %in% e.a(graph)) { weights <- NULL }
  
  buffer <- .Call("R_igraph_write_graph_ncol", graph,
                  names, weights,
                  PACKAGE="igraph")
  write.graph.fromraw(buffer, file)
  
  invisible(NULL)
}  

read.graph.lgl <- function(filename, names=TRUE,
                           weights=TRUE, ...) {

  buffer <- read.graph.toraw(filename)
  .Call("R_igraph_read_graph_lgl", buffer,
        as.logical(names), as.logical(weights),
        PACKAGE="igraph")
}

write.graph.lgl <- function(graph, file, 
                            names="name", weights="weight",
                            isolates=FALSE, ...) {
  names <- as.character(names)
  weights <- as.character(weights)
  if (length(names)==0 || ! names %in% v.a(graph)) { names <- NULL }
  if (length(weights)==0 || ! weights %in% e.a(graph)) { weights <- NULL }
  
  buffer <- .Call("R_igraph_write_graph_lgl", graph,
                  names, weights, as.logical(isolates),
                  PACKAGE="igraph")
  write.graph.fromraw(buffer, file)
  
  invisible(NULL)
}  
