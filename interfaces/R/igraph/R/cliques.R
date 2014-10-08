#   IGraph R package
#   Copyright (C) 2006-2012  Gabor Csardi <csardi.gabor@gmail.com>
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



#' The functions find cliques, ie. complete subgraphs in a graph
#' 
#' These functions find all, the largest or all the maximal cliques in an
#' undirected graph. The size of the largest clique can also be calculated.
#' 
#' \code{cliques} find all complete subgraphs in the input graph, obeying the
#' size limitations given in the \code{min} and \code{max} arguments.
#' 
#' \code{largest_cliques} finds all largest cliques in the input graph. A
#' clique is largest if there is no other clique including more vertices.
#' 
#' \code{max_cliques} finds all maximal cliques in the input graph.  A
#' clique in maximal if it cannot be extended to a larger clique. The largest
#' cliques are always maximal, but a maximal clique is not neccessarily the
#' largest.
#' 
#' \code{count_max_cliques} counts the maximal cliques.
#' 
#' \code{clique_num} calculates the size of the largest clique(s).
#' 
#' The current implementation of these functions searches for maximal
#' independent vertex sets (see \code{\link{ivs}}) in the
#' complementer graph.
#' 
#' @aliases cliques largest_cliques maximal.cliques maximal.cliques.count
#' clique.number clique_num largest.cliques count_max_cliques max_cliques
#' @param graph The input graph, directed graphs will be considered as
#' undirected ones, multiple edges and loops are ignored.
#' @param min Numeric constant, lower limit on the size of the cliques to find.
#' \code{NULL} means no limit, ie. it is the same as 0.
#' @param max Numeric constant, upper limit on the size of the cliques to find.
#' \code{NULL} means no limit.
#' @return \code{cliques}, \code{largest_cliques} and \code{clique_num}
#' return a list containing numeric vectors of vertex ids. Each list element is
#' a clique.
#' 
#' \code{max_cliques} returns \code{NULL}, invisibly, if its \code{file}
#' argument is not \code{NULL}. The output is written to the specified file in
#' this case.
#' 
#' \code{clique_num} and \code{count_max_cliques} return an integer
#' scalar.
#' @author Tamas Nepusz \email{ntamas@@gmail.com} and Gabor Csardi
#' \email{csardi.gabor@@gmail.com}
#' @seealso \code{\link{ivs}}
#' @references For maximal cliques the following algorithm is implemented:
#' David Eppstein, Maarten Loffler, Darren Strash: Listing All Maximal Cliques
#' in Sparse Graphs in Near-optimal Time.  \url{http://arxiv.org/abs/1006.5440}
#' @export
#' @keywords graphs
#' @examples
#' 
#' # this usually contains cliques of size six
#' g <- sample_gnp(100, 0.3)
#' clique_num(g)
#' cliques(g, min=6)
#' largest_cliques(g)
#' 
#' # To have a bit less maximal cliques, about 100-200 usually
#' g <- sample_gnp(100, 0.03)
#' max_cliques(g)
#' 
#' 
cliques <- function(graph, min=NULL, max=NULL) {
  if (!is_igraph(graph)) {
    stop("Not a graph object")
  }

  if (is.null(min)) {
    min <- 0
  }
  if (is.null(max)) {
    max <- 0
  }

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  res <- .Call("R_igraph_cliques", graph, as.numeric(min), as.numeric(max),
               PACKAGE="igraph")
  lapply(res, function(x) x+1)
}

#' @export

largest_cliques <- function(graph) {
  if (!is_igraph(graph)) {
    stop("Not a graph object")
  }

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  res <- .Call("R_igraph_largest_cliques", graph,
               PACKAGE="igraph")
  lapply(res, function(x) x+1)  
}

#' @rdname cliques
#' @param subset If not \code{NULL}, then it must be a vector of vertex ids,
#' numeric or symbolic if the graph is named. The algorithm is run from these
#' vertices only, so only a subset of all maximal cliques is returned. See the
#' Eppstein paper for details. This argument makes it possible to easily
#' parallelize the finding of maximal cliques.
#' @param file If not \code{NULL}, then it must be a file name, i.e. a
#' character scalar. The output of the algorithm is written to this file. (If
#' it exists, then it will be overwritten.) Each clique will be a separate line
#' in the file, given with the numeric ids of its vertices, separated by
#' whitespace.
#' @export

max_cliques <- function(graph, min=NULL, max=NULL,
                            subset=NULL, file=NULL) {
  if (!is_igraph(graph)) {
    stop("Not a graph object");
  }

  if (is.null(min)) { min <- 0 }
  if (is.null(max)) { max <- 0 }

  if (!is.null(subset)) {
    subset <- as.integer(as.igraph.vs(graph, subset)-1)
  }
  
  if (!is.null(file)) {
    if (!is.character(file) ||
        length(grep("://", file, fixed=TRUE)) > 0 ||
        length(grep("~", file, fixed=TRUE)) > 0) {
      tmpfile <- TRUE
      origfile <- file
      file <- tempfile()
    } else {
      tmpfile <- FALSE
    }
    on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
    res <- .Call("R_igraph_maximal_cliques_file", graph, subset, file,
                 as.numeric(min), as.numeric(max), PACKAGE="igraph")
    if (tmpfile) {
      buffer <- read.graph.toraw(file)
      write.graph.fromraw(buffer, origfile)
    }
    invisible(NULL)
  } else { 
    on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
    res <- .Call("R_igraph_maximal_cliques", graph, subset,
                 as.numeric(min), as.numeric(max),
                 PACKAGE="igraph")
    lapply(res, function(x) x+1)
  }
}

#' @export

count_max_cliques <- function(graph, min=NULL, max=NULL,
                                  subset=NULL) {
  # Argument checks
  if (!is_igraph(graph)) { stop("Not a graph object") }

  if (is.null(min)) { min <- 0 }
  if (is.null(max)) { max <- 0 }
  min <- as.integer(min)
  max <- as.integer(max)

  if (!is.null(subset)) {
    subset <- as.integer(as.igraph.vs(graph, subset)-1)
  }

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  # Function call
  res <- .Call("R_igraph_maximal_cliques_count", graph, subset, min, max,
               PACKAGE="igraph")

  res
}

#' @export

clique_num <- function(graph) {
  if (!is_igraph(graph)) {
    stop("Not a graph object");
  }
  
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  .Call("R_igraph_clique_number", graph,
        PACKAGE="igraph")
}



#' Independent vertex sets
#' 
#' A vertex set is called independent if there no edges between any two
#' vertices in it. These functions find independent vertex sets in undirected
#' graphs
#' 
#' \code{ivs} finds all independent vertex sets in the
#' network, obeying the size limitations given in the \code{min} and \code{max}
#' arguments.
#' 
#' \code{largest_ivs} finds the largest independent vertex
#' sets in the graph. An independent vertex set is largest if there is no
#' independent vertex set with more vertices.
#' 
#' \code{maximal_ivs} finds the maximal independent vertex
#' sets in the graph. An independent vertex set is maximal if it cannot be
#' extended to a larger independent vertex set. The largest independent vertex
#' sets are maximal, but the opposite is not always true.
#' 
#' \code{independece.number} calculate the size of the largest independent
#' vertex set(s).
#' 
#' These functions use the algorithm described by Tsukiyama et al., see
#' reference below.
#' 
#' @aliases independent.vertex.sets largest.independent.vertex.sets
#' maximal.independent.vertex.sets independence.number ivs_size ivs
#' largest_ivs maximal_ivs
#' @param graph The input graph, directed graphs are considered as undirected,
#' loop edges and multiple edges are ignored.
#' @param min Numeric constant, limit for the minimum size of the independent
#' vertex sets to find. \code{NULL} means no limit.
#' @param max Numeric constant, limit for the maximum size of the independent
#' vertex sets to find. \code{NULL} means no limit.
#' @return \code{ivs},
#' \code{largest_ivs} and
#' \code{maximal_ivs} return a list containing numeric
#' vertex ids, each list element is an independent vertex set.
#' 
#' \code{ivs_size} returns an integer constant.
#' @author Tamas Nepusz \email{ntamas@@gmail.com} ported it from the Very Nauty
#' Graph Library by Keith Briggs (\url{http://keithbriggs.info/}) and Gabor
#' Csardi \email{csardi.gabor@@gmail.com} wrote the R interface and this manual
#' page.
#' @seealso \code{\link{cliques}}
#' @references S. Tsukiyama, M. Ide, H. Ariyoshi and I. Shirawaka. A new
#' algorithm for generating all the maximal independent sets. \emph{SIAM J
#' Computing}, 6:505--517, 1977.
#' @export
#' @keywords graphs
#' @examples
#' 
#' # A quite dense graph
#' set.seed(42)
#' g <- sample_gnp(100, 0.9)
#' ivs_size(g)
#' ivs(g, min=ivs_size(g))
#' largest_ivs(g)
#' # Empty graph
#' induced_subgraph(g, largest_ivs(g)[[1]])
#' 
#' length(maximal_ivs(g))
#' 
ivs <- function(graph, min=NULL, max=NULL) {
  if (!is_igraph(graph)) {
    stop("Not a graph object");
  }

  if (is.null(min)) {
    min <- 0
  }

  if (is.null(max)) {
    max <- 0
  }
  
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  res <- .Call("R_igraph_independent_vertex_sets", graph, as.numeric(min),
               as.numeric(max),
               PACKAGE="igraph")
  lapply(res, function(x) x+1)
}

#' @export

largest_ivs <- function(graph) {
  if (!is_igraph(graph)) {
    stop("Not a graph object");
  }

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  res <- .Call("R_igraph_largest_independent_vertex_sets", graph,
               PACKAGE="igraph")
  lapply(res, function(x) x+1)
}

#' @export

maximal_ivs <- function(graph) {
  if (!is_igraph(graph)) {
    stop("Not a graph object");
  }

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  res <- .Call("R_igraph_maximal_independent_vertex_sets", graph,
               PACKAGE="igraph")
  lapply(res, function(x) x+1)
}

#' @export

ivs_size <- function(graph) {
  if (!is_igraph(graph)) {
    stop("Not a graph object");
  }

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  .Call("R_igraph_independence_number", graph,
        PACKAGE="igraph")
}
