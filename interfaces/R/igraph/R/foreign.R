
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
    filename <- file(filename, open="rb")
  }
  if (!isOpen(filename)) {
    open(filename, open="rb")
  }

  tmpbufsize <- 20000
  buffer <- tmpbuffer <- readBin(filename, what=raw(0), n=tmpbufsize)
  while (length(tmpbuffer) == tmpbufsize) {
    tmpbuffer <- readBin(filename, what=raw(0), n=tmpbufsize)
    buffer <- c(buffer, tmpbuffer)
  }
  close(filename)
  rm(tmpbuffer)
  
  buffer
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

read.graph <- function(file, format=c("edgelist", "pajek", "ncol", "lgl",
                               "graphml", "dimacs", "graphdb", "gml"), ...) {

  if (igraph.i.have.fmemopen) {
    file <- read.graph.toraw(file)
  } else {
    if (!is.character(file) || length(grep("://", file, fixed=TRUE))>0) {
      buffer <- read.graph.toraw(file)
      file <- tempfile()
      write.graph.fromraw(buffer, file)
    }
  }

  format <- igraph.match.arg(format)
  res <- switch(format,
                "pajek"=read.graph.pajek(file, ...),
                "ncol"=read.graph.ncol(file, ...),
                "edgelist"=read.graph.edgelist(file, ...),
                "lgl"=read.graph.lgl(file, ...),
                "graphml"=read.graph.graphml(file, ...),
                "dimacs"=read.graph.dimacs(file, ...),
                "graphdb"=read.graph.graphdb(file, ...),
                "gml"=read.graph.gml(file, ...),
                stop(paste("Unknown file format:",format))
                )
  res
}

write.graph <- function(graph, file, format=c("edgelist", "pajek", "ncol", "lgl",
                                       "graphml", "dimacs", "gml", "dot"), ...) {

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
  
  format <- igraph.match.arg(format)
  res <- switch(format,
                "pajek"=write.graph.pajek(graph, file, ...),
                "edgelist"=write.graph.edgelist(graph, file, ...),
                "ncol"=write.graph.ncol(graph, file, ...),
                "lgl"=write.graph.lgl(graph, file, ...),
                "graphml"=write.graph.graphml(graph, file, ...),
                "dimacs"=write.graph.dimacs(graph, file, ...),
                "gml"=write.graph.gml(graph, file, ...),
                "dot"=write.graph.dot(graph, file, ...),
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

  if (length(list(...))>0) {
    stop("Unknown arguments to read.graph (edgelist format)")
  }
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph0") )
  .Call("R_igraph_read_graph_edgelist", file,
        as.numeric(n), as.logical(directed),
        PACKAGE="igraph0")
}

write.graph.edgelist <- function(graph, file, ...) {
  
  if (length(list(...))>0) {
    stop("Unknown arguments to write.graph (edgelist format)")
  }
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph0") )
  .Call("R_igraph_write_graph_edgelist", graph, file,
        PACKAGE="igraph0")
}

################################################################
# NCOL and LGL formats, quite simple
################################################################

read.graph.ncol <- function(file, predef=character(0), names=TRUE,
                           weights=TRUE, directed=FALSE, ...) {

  if (length(list(...))>0) {
    stop("Unknown arguments to read.graph (NCOL format)")
  }
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph0") )
  .Call("R_igraph_read_graph_ncol", file, as.character(predef),
        as.logical(names), as.logical(weights), as.logical(directed),
        PACKAGE="igraph0")
}

write.graph.ncol <- function(graph, file, 
                             names="name", weights="weight", ...) {
  if (length(list(...))>0) {
    stop("Unknown arguments to write.graph (NCOL format)")
  }
  names <- as.character(names)
  weights <- as.character(weights)
  if (length(names)==0 || ! names %in% list.vertex.attributes(graph)) { names <- NULL }
  if (length(weights)==0 || ! weights %in% list.edge.attributes(graph)) { weights <- NULL }
  
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph0") )
  .Call("R_igraph_write_graph_ncol", graph, file,
        names, weights,
        PACKAGE="igraph0")
}  

read.graph.lgl <- function(file, names=TRUE,
                           weights=TRUE, ...) {

  if (length(list(...))>0) {
    stop("Unknown arguments to read.graph (LGL format)")
  }
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph0") )
  .Call("R_igraph_read_graph_lgl", file,
        as.logical(names), as.logical(weights),
        PACKAGE="igraph0")
}

write.graph.lgl <- function(graph, file, 
                            names="name", weights="weight",
                            isolates=FALSE, ...) {
  if (length(list(...))>0) {
    stop("Unknown arguments to write.graph (LGL format)")
  }
  names <- as.character(names)
  weights <- as.character(weights)
  if (length(names)==0 || ! names %in% list.vertex.attributes(graph)) { names <- NULL }
  if (length(weights)==0 || ! weights %in% list.edge.attributes(graph)) { weights <- NULL }
  
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph0") )
  .Call("R_igraph_write_graph_lgl", graph, file,
        names, weights, as.logical(isolates),
        PACKAGE="igraph0")
}  

read.graph.pajek <- function(file, ...) {

  if (length(list(...))>0) {
    stop("Unknown arguments to read.graph (Pajek format)")
  }
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph0") )
  .Call("R_igraph_read_graph_pajek", file,
        PACKAGE="igraph0")
}

write.graph.pajek <- function(graph, file, ...) {

  if (length(list(...))>0) {
    stop("Unknown arguments to write.graph (Pajek format)")
  }
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph0") )
  .Call("R_igraph_write_graph_pajek", graph, file,
        PACKAGE="igraph0")
}

read.graph.dimacs <- function(file, directed=TRUE, ...) {

  if (length(list(...))>0) {
    stop("Unknown arguments to read.graph (DIMACS format)")
  }
  res <- .Call("R_igraph_read_graph_dimacs", file, as.logical(directed),
               PACKAGE="igraph0")
  if (res[[1]][1] == "max") {
    graph <- res[[2]]
    graph <- set.graph.attribute(graph, "problem", res[[1]])
    graph <- set.graph.attribute(graph, "source", res[[3]])
    graph <- set.graph.attribute(graph, "target", res[[4]])
    E(graph)$capacity <- res[[5]]
    graph
  } else if (res[[1]][1] == "edge") {
    graph <- res[[2]]
    graph <- set.graph.attribute(graph, "problem", res[[1]])
    V(graph)$label <- res[[3]]
    graph
  }
}

write.graph.dimacs <- function(graph, file,
                               source=NULL, target=NULL, capacity=NULL, ...) {

  if (length(list(...))>0) {
    stop("Unknown arguments to write.graph (DIMACS format)")
  }
  if (is.null(source)) {
    source <- get.graph.attribute(graph, "source")
  }
  if (is.null(target)) {
    target <- get.graph.attribute(graph, "target")
  }
  if (is.null(capacity)) {
    capacity <- E(graph)$capacity
  }
  
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph0") )
  .Call("R_igraph_write_graph_dimacs", graph, file, as.numeric(source),
        as.numeric(target), as.numeric(capacity),
        PACKAGE="igraph0")
}

################################################################
# GraphML
################################################################

read.graph.graphml <- function(file, index=0, ...) {

  if (length(list(...))>0) {
    stop("Unknown arguments to read.graph (GraphML format)")
  }
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph0") )
  .Call("R_igraph_read_graph_graphml", file, as.numeric(index),
        PACKAGE="igraph0")
}

write.graph.graphml <- function(graph, file, ...) {

  if (length(list(...))>0) {
    stop("Unknown arguments to write.graph (GraphML format)")
  }
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph0") )
  .Call("R_igraph_write_graph_graphml", graph, file,
        PACKAGE="igraph0")
}

################################################################
# GML
################################################################

read.graph.gml <- function(file, ...) {
  if (length(list(...))>0) {
    stop("Unknown arguments to read.graph (GML format)")
  }
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph0") )
  .Call("R_igraph_read_graph_gml", file,
        PACKAGE="igraph0")
}

write.graph.gml <- function(graph, file, id=NULL, creator=NULL, ...) {
  if (length(list(...))>0) {
    stop("Unknown arguments to write.graph (GML format)")
  }
  if (!is.null(id)) {
    id <- as.numeric(id)
  }
  if (!is.null(creator)) {
    creator <- as.character(creator)
  }
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph0") )
  .Call("R_igraph_write_graph_gml", graph, file, id, creator,
        PACKAGE="igraph0")
}

################################################################
# Dot
################################################################

write.graph.dot <- function(graph, file, ...) {
  if (length(list(...))>0) {
    stop("Unknown arguments to write.graph (DOT format)")
  }
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph0") )
  .Call("R_igraph_write_graph_dot", graph, file,
        PACKAGE="igraph0")
}

################################################################
# Download a file from the graph database for
# isomorphic problems
################################################################

graph.graphdb <- function(url=NULL,
                          prefix="iso", type="r001", nodes=NULL, pair="A", which=0,
                          base="http://cneurocvs.rmki.kfki.hu/graphdb/gzip",
                          compressed=TRUE, directed=TRUE) {
  
  if (is.null(nodes) && is.null(url)) {
    stop("The `nodes' or the `url' argument must be non-null")
  }

  if (is.null(url)) {
  
    prefixes <- c("iso", "si6", "mcs10", "mcs30", "mcs50", "mcs70", "mcs90")
    types <- c("r001", "r005", "r01", "r02", "m2D", "m2Dr2", "m2Dr4", "m2Dr6",
               "m3D", "m3Dr2", "m3Dr4", "m3Dr6", "m4D", "m4Dr2", "m4Dr4",
               "m4Dr6", "b03", "b03m", "b06", "b06m", "b09", "b09m")
    sizecode <- if (nodes<=100) "s" else if (nodes<2000) "m" else "l" # "l" ????
    typegroups <- c("rand", "rand", "rand", "rand",
                    "m2D", "m2D", "m2D", "m2D",
                    "m2D", "m3D", "m3D", "m3D",
                    "m4D", "m4D", "m4D", "m4D",
                    "bvg", "bvg", "bvg", "bvg", "bvg", "bvg")
    typegroup <- typegroups[which(types==type)]
    
    if (!prefix %in% prefixes) {
      stop("Invalid prefix!")
    }
    if (!type %in% types) {
      stop("Invalid graph type!")
    }
    suff <- if (compressed) ".gz" else ""
    filename <- paste(sep="", base, "/", prefix, "/", typegroup, "/", type, "/",
                      prefix, "_", type, "_", sizecode, nodes,
                      ".", pair, formatC(which, width=2, flag="0"), suff)
  } else {
    filename <- url
  }
  
  ## ok, we have the filename

  f <- try(gzcon(file(filename, open="rb")))
  if (inherits(f, "try-error")) {
    stop(paste("Cannot open URL:", filename));
  }

  if (igraph.i.have.fmemopen) {
    f <- read.graph.toraw(f)
  } else {
    buffer <- read.graph.toraw(f)
    f <- tempfile()
    write.graph.fromraw(buffer, f)
  }

  .Call("R_igraph_read_graph_graphdb", f, as.logical(directed),
        PACKAGE="igraph0")
  
}

read.graph.graphdb <- function(file, directed=TRUE, ...) {
  if (length(list(...))>0) {
    stop("Unknown arguments to read.graph (GraphDB format)")
  }
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph0") )
  .Call("R_igraph_read_graph_graphdb", file, as.logical(directed),
        PACKAGE="igraph0")
}
