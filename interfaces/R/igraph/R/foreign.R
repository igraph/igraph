#   IGraph R package
#   Copyright (C) 2005-2012  Gabor Csardi <csardi.gabor@gmail.com>
#   334 Harvard street, Cambridge, MA 02139 USA
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



#' Reading foreign file formats
#' 
#' The \code{read_graph} function is able to read graphs in various
#' representations from a file, or from a http connection. Currently some
#' simple formats are supported.
#' 
#' The \code{read_graph} function may have additional arguments depending on
#' the file format (the \code{format} argument). See the details separately for
#' each file format, below.
#' 
#' @aliases read.graph LGL Pajek GraphML GML DL UCINET
#' @param file The connection to read from. This can be a local file, or a
#' \code{http} or \code{ftp} connection. It can also be a character string with
#' the file name or URI.
#' @param format Character constant giving the file format. Right now
#' \code{as_edgelist}, \code{pajek}, \code{graphml}, \code{gml}, \code{ncol},
#' \code{lgl}, \code{dimacs} and \code{graphdb} are supported, the default is
#' \code{edgelist}. As of igraph 0.4 this argument is case insensitive.
#' @param \dots Additional arguments, see below.
#' @return A graph object.
#' @section Edge list format: This format is a simple text file with numeric
#' vertex ids defining the edges. There is no need to have newline characters
#' between the edges, a simple space will also do.
#' 
#' Additional arguments: \describe{ \item{n}{The number of vertices in the
#' graph. If it is smaller than or equal to the largest integer in the file,
#' then it is ignored; so it is safe to set it to zero (the default).}
#' \item{directed}{Logical scalar, whether to create a directed graph. The
#' default value is \code{TRUE}.} }
#' @author Gabor Csardi \email{csardi.gabor@@gmail.com}
#' @seealso \code{\link{write_graph}}
#' @keywords graphs
#' @export

read_graph <- function(file, format=c("edgelist", "pajek", "ncol", "lgl",
                               "graphml", "dimacs", "graphdb", "gml", "dl"),
                       ...) {

  if (!is.character(file) || length(grep("://", file, fixed=TRUE)) > 0 ||
      length(grep("~", file, fixed=TRUE)) > 0) {
    buffer <- read.graph.toraw(file)
    file <- tempfile()
    write.graph.fromraw(buffer, file)
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
                "dl"=read.graph.dl(file, ...),
                stop(paste("Unknown file format:",format))
                )
  res
}



#' Writing the graph to a file in some format
#' 
#' \code{write_graph} is a general function for exporting graphs to foreign
#' file formats, however not many formats are implemented right now.
#' 
#' @aliases write.graph
#' @param graph The graph to export.
#' @param file A connection or a string giving the file name to write the graph
#' to.
#' @param format Character string giving the file format. Right now
#' \code{pajek}, \code{graphml}, \code{dot}, \code{gml}, \code{edgelist},
#' \code{lgl}, \code{ncol} and \code{dimacs} are implemented. As of igraph 0.4
#' this argument is case insensitive.
#' @param \dots Other, format specific arguments, see below.
#' @return A NULL, invisibly.
#' @section Edge list format: The \code{edgelist} format is a simple text file,
#' with one edge in a line, the two vertex ids separated by a space character.
#' The file is sorted by the first and the second column. This format has no
#' additional arguments.
#' @author Gabor Csardi \email{csardi.gabor@@gmail.com}
#' @seealso \code{\link{read_graph}}
#' @references Adai AT, Date SV, Wieland S, Marcotte EM. LGL: creating a map of
#' protein function with an algorithm for visualizing very large biological
#' networks. \emph{J Mol Biol.} 2004 Jun 25;340(1):179-90.
#' @export
#' @keywords graphs
#' @examples
#' 
#' g <- make_ring(10)
#' \dontrun{write_graph(g, "/tmp/g.txt", "edgelist")}
#' 
write_graph <- function(graph, file, format=c("edgelist", "pajek", "ncol", "lgl",
                                       "graphml", "dimacs", "gml", "dot", "leda"), ...) {

  if (!is_igraph(graph)) {
    stop("Not a graph object")
  }
  if (!is.character(file) || length(grep("://", file, fixed=TRUE)) > 0 ||
      length(grep("~", file, fixed=TRUE)) > 0) {
    tmpfile <- TRUE
    origfile <- file
    file <- tempfile()
  } else {
    tmpfile <- FALSE
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
                "leda"=write.graph.leda(graph, file, ...),
                stop(paste("Unknown file format:",format))
                )

  if (tmpfile) {
    buffer <- read.graph.toraw(file)
    write.graph.fromraw(buffer, origfile)
  }
  
  invisible(res)
}

################################################################
# Plain edge list format, not sorted
################################################################

read.graph.edgelist <- function(file, n=0,
                                directed=TRUE, ...) {

  if (length(list(...))>0) {
    stop("Unknown arguments to read_graph (edgelist format)")
  }
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  .Call("R_igraph_read_graph_edgelist", file,
        as.numeric(n), as.logical(directed),
        PACKAGE="igraph")
}

write.graph.edgelist <- function(graph, file, ...) {
  
  if (length(list(...))>0) {
    stop("Unknown arguments to write_graph (edgelist format)")
  }
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  .Call("R_igraph_write_graph_edgelist", graph, file,
        PACKAGE="igraph")
}

################################################################
# NCOL and LGL formats, quite simple
################################################################

read.graph.ncol <- function(file, predef=character(0), names=TRUE,
                            weights=c("auto", "yes", "no"),
                            directed=FALSE, ...) {

  if (length(list(...))>0) {
    stop("Unknown arguments to read_graph (NCOL format)")
  }
  weights <- switch(igraph.match.arg(weights), "no"=0, "yes"=1, "auto"=2)
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  .Call("R_igraph_read_graph_ncol", file, as.character(predef),
        as.logical(names), as.numeric(weights), as.logical(directed),
        PACKAGE="igraph")
}

write.graph.ncol <- function(graph, file, 
                             names="name", weights="weight", ...) {
  if (length(list(...))>0) {
    stop("Unknown arguments to write_graph (NCOL format)")
  }
  names <- as.character(names)
  weights <- as.character(weights)
  if (length(names)==0 || ! names %in% vertex_attr_names(graph)) { names <- NULL }
  if (length(weights)==0 || ! weights %in% edge_attr_names(graph)) { weights <- NULL }
  
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  .Call("R_igraph_write_graph_ncol", graph, file,
        names, weights,
        PACKAGE="igraph")
}  

read.graph.lgl <- function(file, names=TRUE,
                           weights=c("auto", "yes", "no"),
                           directed=FALSE,
                            ...) {

  if (length(list(...))>0) {
    stop("Unknown arguments to read_graph (LGL format)")
  }
  weights <- switch(igraph.match.arg(weights), "no"=0, "yes"=1, "auto"=2)
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  .Call("R_igraph_read_graph_lgl", file,
        as.logical(names), as.numeric(weights), as.logical(directed),
        PACKAGE="igraph")
}

write.graph.lgl <- function(graph, file, 
                            names="name", weights="weight",
                            isolates=FALSE, ...) {
  if (length(list(...))>0) {
    stop("Unknown arguments to write_graph (LGL format)")
  }
  names <- as.character(names)
  weights <- as.character(weights)
  if (length(names)==0 || ! names %in% vertex_attr_names(graph)) { names <- NULL }
  if (length(weights)==0 || ! weights %in% edge_attr_names(graph)) { weights <- NULL }
  
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  .Call("R_igraph_write_graph_lgl", graph, file,
        names, weights, as.logical(isolates),
        PACKAGE="igraph")
}  

read.graph.pajek <- function(file, ...) {

  if (length(list(...))>0) {
    stop("Unknown arguments to read_graph (Pajek format)")
  }
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  res <- .Call("R_igraph_read_graph_pajek", file,
               PACKAGE="igraph")
  if ("type" %in% vertex_attr_names(res)) {
    type <- as.logical(V(res)$type)
    res <- delete_vertex_attr(res, "type")
    res <- set_vertex_attr(res, "type", value=type)
  }
  res
}

write.graph.pajek <- function(graph, file, ...) {

  if (length(list(...))>0) {
    stop("Unknown arguments to write_graph (Pajek format)")
  }
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  .Call("R_igraph_write_graph_pajek", graph, file,
        PACKAGE="igraph")
}

read.graph.dimacs <- function(file, directed=TRUE, ...) {

  if (length(list(...))>0) {
    stop("Unknown arguments to read_graph (DIMACS format)")
  }
  res <- .Call("R_igraph_read_graph_dimacs", file, as.logical(directed),
               PACKAGE="igraph")
  if (res[[1]][1] == "max") {
    graph <- res[[2]]
    graph <- set_graph_attr(graph, "problem", res[[1]])
    graph <- set_graph_attr(graph, "source", res[[3]])
    graph <- set_graph_attr(graph, "target", res[[4]])
    E(graph)$capacity <- res[[5]]
    graph
  } else if (res[[1]][1] == "edge") {
    graph <- res[[2]]
    graph <- set_graph_attr(graph, "problem", res[[1]])
    V(graph)$label <- res[[3]]
    graph
  }
}

write.graph.dimacs <- function(graph, file,
                               source=NULL, target=NULL, capacity=NULL, ...) {

  if (length(list(...))>0) {
    stop("Unknown arguments to write_graph (DIMACS format)")
  }
  if (is.null(source)) {
    source <- graph_attr(graph, "source")
  }
  if (is.null(target)) {
    target <- graph_attr(graph, "target")
  }
  if (is.null(capacity)) {
    capacity <- E(graph)$capacity
  }
  
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  .Call("R_igraph_write_graph_dimacs", graph, file, as.numeric(source),
        as.numeric(target), as.numeric(capacity),
        PACKAGE="igraph")
}

################################################################
# GraphML
################################################################

read.graph.graphml <- function(file, index=0, ...) {

  if (length(list(...))>0) {
    stop("Unknown arguments to read_graph (GraphML format)")
  }
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  .Call("R_igraph_read_graph_graphml", file, as.numeric(index),
        PACKAGE="igraph")
}

write.graph.graphml <- function(graph, file, prefixAttr=TRUE, ...) {

  if (length(list(...))>0) {
    stop("Unknown arguments to write_graph (GraphML format)")
  }
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  .Call("R_igraph_write_graph_graphml", graph, file, as.logical(prefixAttr),
        PACKAGE="igraph")
}

################################################################
# GML
################################################################

read.graph.gml <- function(file, ...) {
  if (length(list(...))>0) {
    stop("Unknown arguments to read_graph (GML format)")
  }
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  .Call("R_igraph_read_graph_gml", file,
        PACKAGE="igraph")
}

write.graph.gml <- function(graph, file, id=NULL, creator=NULL, ...) {
  if (length(list(...))>0) {
    stop("Unknown arguments to write_graph (GML format)")
  }
  if (!is.null(id)) {
    id <- as.numeric(id)
  }
  if (!is.null(creator)) {
    creator <- as.character(creator)
  }
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  .Call("R_igraph_write_graph_gml", graph, file, id, creator,
        PACKAGE="igraph")
}

################################################################
# UCINET DL
################################################################

read.graph.dl <- function(file, directed=TRUE, ...) {
  if (length(list(...))>0) {
    stop("Unknown arguments to read_graph (DL format)")
  }
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  .Call("R_igraph_read_graph_dl", file, as.logical(directed),
        PACKAGE="igraph")
}  

################################################################
# Dot
################################################################

write.graph.dot <- function(graph, file, ...) {
  if (length(list(...))>0) {
    stop("Unknown arguments to write_graph (DOT format)")
  }
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  .Call("R_igraph_write_graph_dot", graph, file,
        PACKAGE="igraph")
}

################################################################
# Download a file from the graph database for
# isomorphic problems
################################################################



#' Load a graph from the graph database for testing graph isomorphism.
#' 
#' This function downloads a graph from a database created for the evaluation
#' of graph isomorphism testing algothitms.
#' 
#' \code{graph_from_graphdb} reads a graph from the graph database from an FTP or
#' HTTP server or from a local copy. It has two modes of operation:
#' 
#' If the \code{url} argument is specified then it should the complete path to
#' a local or remote graph database file. In this case we simply call
#' \code{\link{read_graph}} with the proper arguments to read the file.
#' 
#' If \code{url} is \code{NULL}, and this is the default, then the filename is
#' assembled from the \code{base}, \code{prefix}, \code{type}, \code{nodes},
#' \code{pair} and \code{which} arguments.
#' 
#' Unfortunately the original graph database homepage is now defunct, but see
#' its old version at
#' \url{http://web.archive.org/web/20090215182331/http://amalfi.dis.unina.it/graph/db/doc/graphdbat.html}
#' for the actual format of a graph database file and other information.
#'
#' @aliases graph.graphdb
#' @param url If not \code{NULL} it is a complete URL with the file to import.
#' @param prefix Gives the prefix. See details below. Possible values:
#' \code{iso}, \code{i2}, \code{si4}, \code{si6}, \code{mcs10}, \code{mcs30},
#' \code{mcs50}, \code{mcs70}, \code{mcs90}.
#' @param type Gives the graph type identifier. See details below. Possible
#' values: \code{r001}, \code{r005}, \code{r01}, \code{r02}, \code{m2D},
#' \code{m2Dr2}, \code{m2Dr4}, \code{m2Dr6} \code{m3D}, \code{m3Dr2},
#' \code{m3Dr4}, \code{m3Dr6}, \code{m4D}, \code{m4Dr2}, \code{m4Dr4},
#' \code{m4Dr6}, \code{b03}, \code{b03m}, \code{b06}, \code{b06m}, \code{b09},
#' \code{b09m}.
#' @param nodes The number of vertices in the graph.
#' @param pair Specifies which graph of the pair to read. Possible values:
#' \code{A} and \code{B}.
#' @param which Gives the number of the graph to read. For every graph type
#' there are a number of actual graphs in the database. This argument specifies
#' which one to read.
#' @param base The base address of the database. See details below.
#' @param compressed Logical constant, if TRUE than the file is expected to be
#' compressed by gzip. If \code{url} is \code{NULL} then a \sQuote{\code{.gz}}
#' suffix is added to the filename.
#' @param directed Logical constant, whether to create a directed graph.
#' @return A new graph object.
#' @author Gabor Csardi \email{csardi.gabor@@gmail.com}
#' @seealso \code{\link{read_graph}}, \code{\link{graph.isomorphic.vf2}}
#' @references M. De Santo, P. Foggia, C. Sansone, M. Vento: A large database
#' of graphs and its use for benchmarking graph isomorphism algorithms,
#' \emph{Pattern Recognition Letters}, Volume 24, Issue 8 (May 2003)
#' @export
#' @keywords graphs
#' @examples
#' 
#' \dontrun{
#' g <- graph_from_graphdb(prefix="iso", type="r001", nodes=20, pair="A",
#'   which=10, compressed=TRUE)
#' g2 <- graph_from_graphdb(prefix="iso", type="r001", nodes=20, pair="B",
#'   which=10, compressed=TRUE)
#' graph.isomorphic.vf2(g, g2)	% should be TRUE
#' g3 <- graph_from_graphdb(url=paste(sep="/",
#'                               "http://cneurocvs.rmki.kfki.hu",
#'                               "graphdb/gzip/iso/bvg/b06m",
#'                               "iso_b06m_m200.A09.gz"))
#' }
graph_from_graphdb <- function(url=NULL,
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

  buffer <- read.graph.toraw(f)
  f <- tempfile()
  write.graph.fromraw(buffer, f)

  .Call("R_igraph_read_graph_graphdb", f, as.logical(directed),
        PACKAGE="igraph")
  
}

read.graph.graphdb <- function(file, directed=TRUE, ...) {
  if (length(list(...))>0) {
    stop("Unknown arguments to read_graph (GraphDB format)")
  }
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  .Call("R_igraph_read_graph_graphdb", file, as.logical(directed),
        PACKAGE="igraph")
}

write.graph.leda <- function(graph, file, vertex.attr=NULL, edge.attr=NULL,
                             ...) {
  if (length(list(...))>0) {
    stop("Unknown arguments to write_graph (LEDA format)")
  }
  if (!is.null(vertex.attr)) { vertex.attr <- as.character(vertex.attr) }
  if (!is.null(edge.attr))   { edge.attr   <- as.character(edge.attr)   }
  on.exit(.Call("R_igraph_finalizer", PACKAGE="igraph") )
  .Call("R_igraph_write_graph_leda", graph, file, vertex.attr, edge.attr,
        PACKAGE="igraph")
}
