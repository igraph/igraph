
## ----------------------------------------------------------------
##
##   IGraph R package
##   Copyright (C) 2005-2014  Gabor Csardi <csardi.gabor@gmail.com>
##   334 Harvard street, Cambridge, MA 02139 USA
##
##   This program is free software; you can redistribute it and/or modify
##   it under the terms of the GNU General Public License as published by
##   the Free Software Foundation; either version 2 of the License, or
##   (at your option) any later version.
##
##   This program is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY; without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##   GNU General Public License for more details.
##
##   You should have received a copy of the GNU General Public License
##   along with this program; if not, write to the Free Software
##   Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA
##   02110-1301 USA
##
## -----------------------------------------------------------------

#' Creating igraph graphs from data frames or vice-versa
#'
#' This function creates an igraph graph from one or two data frames containing
#' the (symbolic) edge list and edge/vertex attributes.
#'
#' \code{graph_from_data_frame} creates igraph graphs from one or two data frames.
#' It has two modes of operatation, depending whether the \code{vertices}
#' argument is \code{NULL} or not.
#'
#' If \code{vertices} is \code{NULL}, then the first two columns of \code{d}
#' are used as a symbolic edge list and additional columns as edge attributes.
#' The names of the attributes are taken from the names of the columns.
#'
#' If \code{vertices} is not \code{NULL}, then it must be a data frame giving
#' vertex metadata. The first column of \code{vertices} is assumed to contain
#' symbolic vertex names, this will be added to the graphs as the
#' \sQuote{\code{name}} vertex attribute. Other columns will be added as
#' additional vertex attributes. If \code{vertices} is not \code{NULL} then the
#' symbolic edge list given in \code{d} is checked to contain only vertex names
#' listed in \code{vertices}.
#'
#' Typically, the data frames are exported from some speadsheat software like
#' Excel and are imported into R via \code{\link{read.table}},
#' \code{\link{read.delim}} or \code{\link{read.csv}}.
#'
#' \code{as_data_frame} converts the igraph graph into one or more data
#' frames, depending on the \code{what} argument.
#'
#' If the \code{what} argument is \code{edges} (the default), then the edges of
#' the graph and also the edge attributes are returned. The edges will be in
#' the first two columns, named \code{from} and \code{to}. (This also denotes
#' edge direction for directed graphs.)  For named graphs, the vertex names
#' will be included in these columns, for other graphs, the numeric vertex ids.
#' The edge attributes will be in the other columns. It is not a good idea to
#' have an edge attribute named \code{from} or \code{to}, because then the
#' column named in the data frame will not be unique. The edges are listed in
#' the order of their numeric ids.
#'
#' If the \code{what} argument is \code{vertices}, then vertex attributes are
#' returned. Vertices are listed in the order of their numeric vertex ids.
#'
#' If the \code{what} argument is \code{both}, then both vertex and edge data
#' is returned, in a list with named entries \code{vertices} and \code{edges}.
#'
#' @aliases graph_from_data_frame graph.data.frame as_data_frame get.data.frame
#' @param d A data frame containing a symbolic edge list in the first two
#' columns. Additional columns are considered as edge attributes.  Since
#' version 0.7 this argument is coerced to a data frame with
#' \code{as.data.frame}.
#' @param directed Logical scalar, whether or not to create a directed graph.
#' @param vertices A data frame with vertex metadata, or \code{NULL}. See
#' details below. Since version 0.7 this argument is coerced to a data frame
#' with \code{as.data.frame}, if not \code{NULL}.
#' @return An igraph graph object for \code{graph_from_data_frame}, and either a
#' data frame or a list of two data frames named \code{edges} and
#' \code{vertices} for \code{as.data.frame}.
#' @note For \code{graph_from_data_frame} \code{NA} elements in the first two
#' columns \sQuote{d} are replaced by the string \dQuote{NA} before creating
#' the graph. This means that all \code{NA}s will correspond to a single
#' vertex.
#'
#' \code{NA} elements in the first column of \sQuote{vertices} are also
#' replaced by the string \dQuote{NA}, but the rest of \sQuote{vertices} is not
#' touched. In other words, vertex names (=the first column) cannot be
#' \code{NA}, but other vertex attributes can.
#' @author Gabor Csardi \email{csardi.gabor@@gmail.com}
#' @seealso \code{\link{graph_from_literal}}
#' for another way to create graphs, \code{\link{read.table}} to read in tables
#' from files.
#' @keywords graphs
#' @examples
#'
#' ## A simple example with a couple of actors
#' ## The typical case is that these tables are read in from files....
#' actors <- data.frame(name=c("Alice", "Bob", "Cecil", "David",
#'                             "Esmeralda"),
#'                      age=c(48,33,45,34,21),
#'                      gender=c("F","M","F","M","F"))
#' relations <- data.frame(from=c("Bob", "Cecil", "Cecil", "David",
#'                                "David", "Esmeralda"),
#'                         to=c("Alice", "Bob", "Alice", "Alice", "Bob", "Alice"),
#'                         same.dept=c(FALSE,FALSE,TRUE,FALSE,FALSE,TRUE),
#'                         friendship=c(4,5,5,2,1,1), advice=c(4,5,5,4,2,3))
#' g <- graph_from_data_frame(relations, directed=TRUE, vertices=actors)
#' print(g, e=TRUE, v=TRUE)
#'
#' ## The opposite operation
#' as_data_frame(g, what="vertices")
#' as_data_frame(g, what="edges")
#'
graph_from_data_frame <- function(d, directed=TRUE, vertices=NULL) {

  d <- as.data.frame(d)
  if (!is.null(vertices)) { vertices <- as.data.frame(vertices) }

  if (ncol(d) < 2) {
    stop("the data frame should contain at least two columns")
  }

  ## Handle if some elements are 'NA'
  if (any(is.na(d[,1:2]))) {
    warning("In `d' `NA' elements were replaced with string \"NA\"")
    d[,1:2][ is.na(d[,1:2]) ] <- 'NA'
  }
  if (!is.null(vertices) && any(is.na(vertices[,1]))) {
    warning("In `vertices[,1]' `NA' elements were replaced with string \"NA\"")
    vertices[,1][is.na(vertices[,1])] <- 'NA'
  }

  names <- unique( c(as.character(d[,1]), as.character(d[,2])) )
  if (!is.null(vertices)) {
    names2 <- names
    vertices <- as.data.frame(vertices)
    if (ncol(vertices) < 1) {
      stop("Vertex data frame contains no rows")
    }
    names <- as.character(vertices[,1])
    if (any(duplicated(names))) {
      stop("Duplicate vertex names")
    }
    if (any(! names2 %in% names)) {
      stop("Some vertex names in edge list are not listed in vertex data frame")
    }
  }

  # create graph
  g <- make_empty_graph(n=0, directed=directed)

  # vertex attributes
  attrs <- list(name=names)
  if (!is.null(vertices)) {
    if (ncol(vertices) > 1) {
      for (i in 2:ncol(vertices)) {
        newval <- vertices[,i]
        if (class(newval) == "factor") {
          newval <- as.character(newval)
        }
        attrs[[ names(vertices)[i] ]] <- newval
      }
    }
  }

  # add vertices
  g <- add_vertices(g, length(names), attr=attrs)

  # create edge list
  from <- as.character(d[,1])
  to <- as.character(d[,2])
  edges <- rbind(match(from, names), match(to,names))

  # edge attributes
  attrs <- list()
  if (ncol(d) > 2) {
    for (i in 3:ncol(d)) {
      newval <- d[,i]
      if (class(newval) == "factor") {
        newval <- as.character(newval)
      }
      attrs[[ names(d)[i] ]] <- newval
    }
  }

  # add the edges
  g <- add_edges(g, edges, attr=attrs)
  g
}

#' @rdname graph_from_data_frame
#' @param ... Passed to \code{graph_from_data_frame}.
#' @export

from_data_frame <- function(...) constructor_spec(graph_from_data_frame, ...)

## -----------------------------------------------------------------

#' Create a graph from an edge list matrix
#'
#' \code{graph_from_edgelist} creates a graph from an edge list. Its argument
#' is a two-column matrix, each row defines one edge. If it is
#' a numeric matrix then its elements are interpreted as vertex ids. If
#' it is a character matrix then it is interpreted as symbolic vertex
#' names and a vertex id will be assigned to each name, and also a
#' \code{name} vertex attribute will be added.
#'
#' @aliases graph.edgelist
#' @concept Edge list
#' @param el The edge list, a two column matrix, character or numeric.
#' @param directed Whether to create a directed graph.
#' @return An igraph graph.
#'
#' @family determimistic constructors
#' @export
#' @examples
#' el <- matrix( c("foo", "bar", "bar", "foobar"), nc = 2, byrow = TRUE)
#' graph_from_edgelist(el)
#'
#' # Create a ring by hand
#' graph_from_edgelist(cbind(1:10, c(2:10, 1)))

graph_from_edgelist <- function(el, directed=TRUE) {

  if (!is.matrix(el) || ncol(el) != 2) {
    stop("graph_from_edgelist expects a matrix with two columns")
  }

  if (nrow(el) == 0) {
    res <- make_empty_graph(directed=directed)
  } else {
    if (is.character(el)) {
      ## symbolic edge list
      names <- unique(as.character(t(el)))
      ids <- seq(names)
      names(ids) <- names
      res <- graph( unname(ids[t(el)]), directed=directed)
      rm(ids)
      V(res)$name <- names
    } else {
      ## normal edge list
      res <- graph( t(el), directed=directed )
    }
  }

  res
}

#' @rdname graph_from_edgelist
#' @param ... Passed to \code{graph_from_edgelist}.
#' @export

from_edgelist <- function(...) constructor_spec(graph_from_edgelist, ...)
