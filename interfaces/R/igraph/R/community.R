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
# Community structure
###################################################################



#' Functions to deal with the result of network community detection
#'
#' igraph community detection functions return their results as an object from
#' the \code{communities} class. This manual page describes the operations of
#' this class.
#' 
#' Community structure detection algorithms try to find dense subgraphs in
#' directed or undirected graphs, by optimizing some criteria, and usually
#' using heuristics.
#' 
#' igraph implements a number of community detection methods (see them below),
#' all of which return an object of the class \code{communities}. Because the
#' community structure detection algorithms are different, \code{communities}
#' objects do not always have the same structure. Nevertheless, they have some
#' common operations, these are documented here.
#' 
#' The \code{print} generic function is defined for \code{communities}, it
#' prints a short summary.
#' 
#' The \code{length} generic function call be called on \code{communities} and
#' returns the number of communities.
#' 
#' The \code{sizes} function returns the community sizes, in the order of their
#' ids.
#' 
#' \code{membership} gives the division of the vertices, into communities. It
#' returns a numeric vector, one value for each vertex, the id of its
#' community. Community ids start from one. Note that some algorithms calculate
#' the complete (or incomplete) hierarchical structure of the communities, and
#' not just a single partitioning. For these algorithms typically the
#' membership for the highest modularity value is returned, but see also the
#' manual pages of the individual algorithms.
#' 
#' \code{communities} is also the name of a function, that returns a list of
#' communities, each identified by their vertices. The vertices will have
#' symbolic names if the \code{add.vertex.names} igraph option is set, and the
#' graph itself was named. Otherwise numeric vertex ids are used.
#' 
#' \code{modularity} gives the modularity score of the partitioning. (See
#' \code{\link{modularity.igraph}} for details. For algorithms that do not
#' result a single partitioning, the highest modularity value is returned.
#' 
#' \code{algorithm} gives the name of the algorithm that was used to calculate
#' the community structure.
#' 
#' \code{crossing} returns a logical vector, with one value for each edge,
#' ordered according to the edge ids. The value is \code{TRUE} iff the edge
#' connects two different communities, according to the (best) membership
#' vector, as returned by \code{membership()}.
#' 
#' \code{is_hierarchical} checks whether a hierarchical algorithm was used to
#' find the community structure. Some functions only make sense for
#' hierarchical methods (e.g. \code{merges}, \code{cut_at} and
#' \code{as.dendrogram}).
#' 
#' \code{merges} returns the merge matrix for hierarchical methods. An error
#' message is given, if a non-hierarchical method was used to find the
#' community structure. You can check this by calling \code{is_hierarchical} on
#' the \code{communities} object.
#' 
#' \code{cut_at} cuts the merge tree of a hierarchical community finding method,
#' at the desired place and returns a membership vector. The desired place can
#' be expressed as the desired number of communities or as the number of merge
#' steps to make. The function gives an error message, if called with a
#' non-hierarchical method.
#' 
#' \code{as.dendrogram} converts a hierarchical community structure to a
#' \code{dendrogram} object. It only works for hierarchical methods, and gives
#' an error message to others. See \code{\link[stats]{dendrogram}} for details.
#' 
#' \code{as.hclust} is similar to \code{as.dendrogram}, but converts a
#' hierarchical community structure to a \code{hclust} object.
#' 
#' \code{as_phylo} converts a hierarchical community structure to a \code{phylo}
#' object, you will need the \code{ape} package for this.
#' 
#' \code{show_trace} works (currently) only for communities found by the leading
#' eigenvector method (\code{\link{cluster_leading_eigen}}), and
#' returns a character vector that gives the steps performed by the algorithm
#' while finding the communities.
#' 
#' \code{code_len} is defined for the InfoMAP method
#' (\code{\link{cluster_infomap}} and returns the code length of the
#' partition.
#' 
#' It is possibly to call the \code{plot} function on \code{communities}
#' objects. This will plot the graph (and uses \code{\link{plot.igraph}}
#' internally), with the communities shown. By default it colores the vertices
#' according to their communities, and also marks the vertex groups
#' corresponding to the communities. It passes additional arguments to
#' \code{\link{plot.igraph}}, please see that and also
#' \code{\link{igraph.plotting}} on how to change the plot.
#' 
#' @rdname communities
#' @aliases communities membership algorithm crossing cutat merges sizes cut_at
#' is.hierarchical print.communities plot.communities length.communities
#' as.dendrogram.communities as.hclust.communities code_len
#' asPhylo asPhylo.communities showtrace code.length
#' as_phylo as_phylo.communities show_trace is_hierarchical
#' @param communities,x,object A \code{communities} object, the result of an
#' igraph community detection function.
#' @param graph An igraph graph object, corresponding to \code{communities}.
#' @param y An igraph graph object, corresponding to the communities in
#' \code{x}.
#' @param no Integer scalar, the desired number of communities. If too low or
#' two high, then an error message is given. Exactly one of \code{no} and
#' \code{steps} must be supplied.
#' @param steps The number of merge operations to perform to produce the
#' communities. Exactly one of \code{no} and \code{steps} must be supplied.
#' @param colbar A vector of colors, in any format that is accepted by the
#' regular R plotting methods. E.g. it may be an integer vector, a character
#' vector of color names, a character vector of RGB colors. This vector gives
#' the color bar for the vertices. The length of the vector should be the same
#' as the number of communities.
#' @param col A vector of colors, in any format that is accepted by the regular
#' R plotting methods. This vector gives the colors of the vertices explicitly.
#' @param mark.groups A list of numeric vectors. The communities can be
#' highlighted using colored polygons. The groups for which the polygons are
#' drawn are given here. The default is to use the groups given by the
#' communities. Supply \code{NULL} here if you do not want to highlight any
#' groups.
#' @param edge.color The colors of the edges. By default the edges within
#' communities are colored green and other edges are red.
#' @param hang Numeric scalar indicating how the height of leaves should be
#' computed from the heights of their parents; see \code{\link{plot.hclust}}.
#' @param use.modularity Logical scalar, whether to use the modularity values
#' to define the height of the branches.
#' @param \dots Additional arguments. \code{plot.communities} passes these to
#' \code{\link{plot.igraph}}. The other functions silently ignore
#' them.
#' @param membership Numeric vector, one value for each vertex, the membership
#' vector of the community structure. Might also be \code{NULL} if the
#' community structure is given in another way, e.g. by a merge matrix.
#' @param algorithm If not \code{NULL} (meaning an unknown algorithm), then a
#' character scalar, the name of the algorithm that produced the community
#' structure.
#' @param merges If not \code{NULL}, then the merge matrix of the hierarchical
#' community structure. See \code{merges} below for more information on its
#' format.
#' @param modularity Numeric scalar or vector, the modularity value of the
#' community structure. It can also be \code{NULL}, if the modularity of the
#' (best) split is not available.
#' @return \code{print} returns the \code{communities} object itself,
#' invisibly.
#' 
#' \code{length} returns an integer scalar.
#' 
#' \code{sizes} returns a numeric vector.
#' 
#' \code{membership} returns a numeric vector, one number for each vertex in
#' the graph that was the input of the community detection.
#' 
#' \code{modularity} returns a numeric scalar.
#' 
#' \code{algorithm} returns a character scalar.
#' 
#' \code{crossing} returns a logical vector.
#' 
#' \code{is_hierarchical} returns a logical scalar.
#' 
#' \code{merges} returns a two-column numeric matrix.
#' 
#' \code{cut_at} returns a numeric vector, the membership vector of the
#' vertices.
#' 
#' \code{as.dendrogram} returns a \code{\link[stats]{dendrogram}} object.
#' 
#' \code{show_trace} returns a character vector.
#' 
#' \code{code_len} returns a numeric scalar for communities found with the
#' InfoMAP method and \code{NULL} for other methods.
#' 
#' \code{plot} for \code{communities} objects returns \code{NULL}, invisibly.
#' 
#' #' @author Gabor Csardi \email{csardi.gabor@@gmail.com}
#' @seealso See \code{\link{plot_dendrogram}} for plotting community structure
#' dendrograms.
#' 
#' See \code{\link{compare}} for comparing two community structures
#' on the same graph.
#' 
#' The different methods for finding communities, they all return a
#' \code{communities} object: \code{\link{cluster_edge_betweenness}},
#' \code{\link{cluster_fast_greedy}},
#' \code{\link{cluster_label_prop}},
#' \code{\link{cluster_leading_eigen}},
#' \code{\link{cluster_louvain}}, \code{\link{cluster_optimal}},
#' \code{\link{cluster_spinglass}}, \code{\link{cluster_walktrap}}.
#' @keywords graphs
#' @export
#' @examples
#' 
#' karate <- make_graph("Zachary")
#' wc <- cluster_walktrap(karate)
#' modularity(wc)
#' membership(wc)
#' plot(wc, karate)
#' 

membership <- function(communities) {
  if (!is.null(communities$membership)) {
    res <- communities$membership
  } else if (!is.null(communities$merges) &&
             !is.null(communities$modularity)) {
    res <- community.to.membership2(communities$merges, communities$vcount,
                                    which.max(communities$modularity))
  } else {
    stop("Cannot calculate community membership")
  }
  if (igraph_opt("add.vertex.names") && !is.null(communities$names)) {
    names(res) <- communities$names
  }
  class(res) <- "membership"
  res
}

#' @method print membership
#' @export

print.membership <- function(x, ...) print(unclass(x), ...)

#' Declare a numeric vector as a membership vector
#'
#' This is useful if you want to use functions defined on
#' membership vectors, but your membership vector does not
#' come from an igraph clustering method.
#'
#' @param x The input vector.
#' @return The input vector, with the \code{membership} class added.
#' @export
#' @examples
#' ## Compare to the correct clustering
#' g <- (make_full_graph(10) + make_full_graph(10)) %>%
#'   rewire(each_edge(p = 0.2))
#' correct <- rep(1:2, each = 10) %>% as_membership
#' fc <- cluster_fast_greedy(g)
#' compare(correct, fc)
#' compare(correct, membership(fc))

as_membership <- function(x) add_class(x, "membership")


#' @rdname communities
#' @method print communities
#' @export
#' @importFrom printr head_print

print.communities <- function(x, ...) {

  noc <- if (!is.null(x$membership)) max(membership(x)) else NA
  mod <- if (!is.null(x$modularity)) {
    modularity(x) %>% format(digits = 2)
  } else {
    NA_real_
  }
  alg <- x$algorithm %||% "unknown"

  cat("IGRAPH clustering ", alg, ", groups: ", noc, ", mod: ", mod, "\n", sep="")

  if (!is.null(x$membership)) {
    grp <- groups(x)
    cat("+ groups:\n")
    hp <- function(o) {
      head_print(o, max_lines = igraph_opt("auto.print.lines"),
                 omitted_footer = "+ ... omitted several groups/vertices\n",)
    }
    indent_print(grp, .printer = hp, .indent = "  ")

  } else {
    cat(" + groups not available\n")
  }

  invisible(x)
}

#' Creates a communities object.
#'
#' This is useful to integrate the results of community finding algorithms
#' that are not included in igraph.
#'
#' @param graph The graph of the community structure.
#' @param membership The membership vector of the community structure, a
#'   numeric vector denoting the id of the community for each vertex. It
#'   might be \code{NULL} for hierarchical community structures.
#' @param algorithm Character string, the algorithm that generated
#'   the community structure, it can be arbitrary.
#' @param merges A merge matrix, for hierarchical community structures (or
#'   \code{NULL} otherwise.
#' @param modularity Modularity value of the community structure. If this
#'   is \code{TRUE} and the membership vector is available, then it the
#'   modularity values is calculated automatically.
#' @return A \code{communities} object.
#'
#' @export

make_clusters <- function(graph, membership = NULL, algorithm = NULL,
                          merges = NULL, modularity = TRUE) {

  stopifnot(is.null(membership) || is.numeric(membership))
  stopifnot(is.null(algorithm) ||
            (is.character(algorithm) && length(algorithm)==1))
  stopifnot(is.null(merges) ||
              (is.matrix(merges) && is.numeric(merges) && ncol(merges)==2))
  stopifnot(is.null(modularity) ||
            (is.logical(modularity) && length(modularity) == 1) ||
            (is.numeric(modularity) &&
             length(modularity) %in% c(1, length(membership))))

  if (is.logical(modularity)) {
    if (modularity && !is.null(membership)) {
      modularity <- modularity(graph, membership)
    } else {
      modularity <- NULL
    }
  }

  res <- list(membership=membership,
              algorithm=if (is.null(algorithm)) "unknown" else algorithm,
              modularity=modularity)
  if (!is.null(merges)) {
    res$merges <- merges
  }
  if (!is.null(membership)) {
    res$vcount <- length(membership)
  } else if (!is.null(merges)) {
    res$vcount <- nrow(merges) + 1
  }
  class(res) <- "communities"
  res
}

#' @export

modularity <- function(x, ...)
  UseMethod("modularity")

#' Modularity of a community structure of a graph
#' 
#' This function calculates how modular is a given division of a graph into
#' subgraphs.
#' 
#' \code{modularity} calculates the modularity of a graph with respect to the
#' given \code{membership} vector.
#' 
#' The modularity of a graph with respect to some division (or vertex types)
#' measures how good the division is, or how separated are the different vertex
#' types from each other. It defined as \deqn{Q=\frac{1}{2m} \sum_{i,j}
#' (A_{ij}-\frac{k_ik_j}{2m})\delta(c_i,c_j),}{Q=1/(2m) * sum( (Aij-ki*kj/(2m)
#' ) delta(ci,cj),i,j),} here \eqn{m} is the number of edges, \eqn{A_{ij}}{Aij}
#' is the element of the \eqn{A} adjacency matrix in row \eqn{i} and column
#' \eqn{j}, \eqn{k_i}{ki} is the degree of \eqn{i}, \eqn{k_j}{kj} is the degree
#' of \eqn{j}, \eqn{c_i}{ci} is the type (or component) of \eqn{i},
#' \eqn{c_j}{cj} that of \eqn{j}, the sum goes over all \eqn{i} and \eqn{j}
#' pairs of vertices, and \eqn{\delta(x,y)}{delta(x,y)} is 1 if \eqn{x=y} and 0
#' otherwise.
#' 
#' If edge weights are given, then these are considered as the element of the
#' \eqn{A} adjacency matrix, and \eqn{k_i}{ki} is the sum of weights of
#' adjacent edges for vertex \eqn{i}.
#' 
#' \code{modularity_matrix} calculates the modularity matrix. This is a dense matrix,
#' and it is defined as the difference of the adjacency matrix and the
#' configuration model null model matrix. In other words element
#' \eqn{M_{ij}}{M[i,j]} is given as \eqn{A_{ij}-d_i
#' d_j/(2m)}{A[i,j]-d[i]d[j]/(2m)}, where \eqn{A_{ij}}{A[i,j]} is the (possibly
#' weighted) adjacency matrix, \eqn{d_i}{d[i]} is the degree of vertex \eqn{i},
#' and \eqn{m} is the number of edges (or the total weights in the graph, if it
#' is weighed).
#'
#' @aliases modularity
#' @param x,graph The input graph.
#' @param membership Numeric vector, for each vertex it gives its community.
#' The communities are numbered from one.
#' @param weights If not \code{NULL} then a numeric vector giving edge weights.
#' @param \dots Additional arguments, none currently.
#' @return For \code{modularity} a numeric scalar, the modularity score of the
#' given configuration.
#' 
#' For \code{modularity_matrix} a numeic square matrix, its order is the number of
#' vertices in the graph.
#' @author Gabor Csardi \email{csardi.gabor@@gmail.com}
#' @seealso \code{\link{cluster_walktrap}},
#' \code{\link{cluster_edge_betweenness}},
#' \code{\link{cluster_fast_greedy}}, \code{\link{cluster_spinglass}} for
#' various community detection methods.
#' @references Clauset, A.; Newman, M. E. J. & Moore, C. Finding community
#' structure in very large networks, \emph{Phyisical Review E} 2004, 70, 066111
#' @method modularity igraph
#' @export
#' @keywords graphs
#' @examples
#' 
#' g <- make_full_graph(5) %du% make_full_graph(5) %du% make_full_graph(5)
#' g <- add_edges(g, c(1,6, 1,11, 6, 11))
#' wtc <- cluster_walktrap(g)
#' modularity(wtc)
#' modularity(g, membership(wtc))
#' 

modularity.igraph <- function(x, membership, weights=NULL, ...) {
  # Argument checks
  if (!is_igraph(x)) { stop("Not a graph object") }
  membership <- as.numeric(membership)
  if (!is.null(weights)) weights <- as.numeric(weights)

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  # Function call
  res <- .Call("R_igraph_modularity", x, membership-1, weights,
        PACKAGE="igraph")
  res
}

#' @rdname communities
#' @method modularity communities
#' @export

modularity.communities <- function(x, ...) {
  if (!is.null(x$modularity)) {
    max(x$modularity)
  } else {
    stop("Modularity was not calculated")
  }
}

#' @rdname modularity.igraph
#' @aliases mod.matrix
#' @export

modularity_matrix <- function(graph, membership, weights=NULL) {
  # Argument checks
  if (!is_igraph(graph)) { stop("Not a graph object") }
  membership <- as.numeric(membership)-1
  if (is.null(weights) && "weight" %in% edge_attr_names(graph)) { 
  weights <- E(graph)$weight 
  } 
  if (!is.null(weights) && any(!is.na(weights))) { 
  weights <- as.numeric(weights) 
  } else { 
  weights <- NULL 
  }

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  # Function call
  res <- .Call("R_igraph_modularity_matrix", graph, membership, weights,
        PACKAGE="igraph")

  res
}

#' @rdname communities
#' @method length communities
#' @export

length.communities <- function(x) {
  m <- membership(x)
  max(m)
}

#' @rdname communities
#' @export

sizes <- function(communities) {
  m <- membership(communities)
  table(`Community sizes`=m)
}

#' @rdname communities
#' @export

algorithm <- function(communities) {
  communities$algorithm
}

#' @rdname communities
#' @export

merges <- function(communities) {
  if (!is.null(communities$merges)) {
    communities$merges
  } else {
    stop("Not a hierarchical community structure")
  }
}

#' @rdname communities
#' @export

crossing <- function(communities, graph) {
  m <- membership(communities)
  el <- as_edgelist(graph, names=FALSE)
  m1 <- m[el[,1]]
  m2 <- m[el[,2]]
  res <- m1 != m2
  if (!is.null(names(m1))) {
    names(res) <- paste(names(m1), names(m2), sep="|")
  }
  res
}

#' @rdname communities
#' @export

code_len <- function(communities) {
  communities$codelength
}

#' @rdname communities
#' @export

is_hierarchical <- function(communities) {
  ! is.null(communities$merges)
}

complete.dend <- function(comm, use.modularity) {
  merges <- comm$merges
  if (nrow(merges) < comm$vcount-1) {
    if (use.modularity) {
      stop(paste("`use.modularity' requires a full dendrogram,",
                 "i.e. a connected graph"))
    }
    miss <- seq_len(comm$vcount + nrow(merges))[-as.vector(merges)]
    miss <- c(miss, seq_len(length(miss)-2) + comm$vcount+nrow(merges))
    miss <- matrix(miss, byrow=TRUE, ncol=2)
    merges <- rbind(merges, miss)
  }
  storage.mode(merges) <- "integer"

  merges
}

# The following functions were adapted from the stats R package

#' @rdname communities
#' @importFrom stats as.dendrogram
#' @method as.dendrogram communities
#' @export
 
as.dendrogram.communities <- function(object, hang=-1, use.modularity=FALSE,
                                      ...) {
  if (!is_hierarchical(object)) {
    stop("Not a hierarchical community structure")
  }

  .memberDend <- function(x) {
    r <- attr(x,"x.member")
    if(is.null(r)) {
      r <- attr(x,"members")
    if(is.null(r)) r <- 1:1
    }
    r
  }

  ## If multiple components, then we merge them in arbitrary order
  merges <- complete.dend(object, use.modularity)
  
  storage.mode(merges) <- "integer"
  
  if (is.null(object$names)) {
    object$names <- 1:(nrow(merges)+1)
  }
  z <- list()
  if (!use.modularity || is.null(object$modularity)) {
    object$height <- 1:nrow(merges)
  } else {
    object$height <- object$modularity[-1]
    object$height <- cumsum(object$height - min(object$height))
  }
  nMerge <- length(oHgt <- object$height)
  if (nMerge != nrow(merges))
    stop("'merge' and 'height' do not fit!")
  hMax <- oHgt[nMerge]
  one <- 1L
  two <- 2L
  leafs <- nrow(merges)+1
  for (k in 1:nMerge) {
    x <- merges[k, ]# no sort() anymore!
    if (any(neg <- x < leafs+1))
      h0 <- if (hang < 0) 0 else max(0, oHgt[k] - hang * hMax)
    if (all(neg)) {                  # two leaves
      zk <- as.list(x)
      attr(zk, "members") <- two
      attr(zk, "midpoint") <- 0.5 # mean( c(0,1) )
      objlabels <- object$names[x]
      attr(zk[[1]], "label") <- objlabels[1]
      attr(zk[[2]], "label") <- objlabels[2]
      attr(zk[[1]], "members") <- attr(zk[[2]], "members") <- one
      attr(zk[[1]], "height") <- attr(zk[[2]], "height") <- h0
      attr(zk[[1]], "leaf") <- attr(zk[[2]], "leaf") <- TRUE
    }
    else if (any(neg)) {            # one leaf, one node
      X <- as.character(x)
      ## Originally had "x <- sort(..) above => leaf always left, x[1];
      ## don't want to assume this
      isL <- x[1] < leafs+1 ## is leaf left?
      zk <-
        if(isL) list(x[1], z[[X[2]]])
        else    list(z[[X[1]]], x[2])
      attr(zk, "members") <- attr(z[[X[1 + isL]]], "members") + one
      attr(zk, "midpoint") <-
        (.memberDend(zk[[1]]) + attr(z[[X[1 + isL]]], "midpoint"))/2
      attr(zk[[2 - isL]], "members") <- one
      attr(zk[[2 - isL]], "height") <- h0
      attr(zk[[2 - isL]], "label") <- object$names[x[2 - isL]]
      attr(zk[[2 - isL]], "leaf") <- TRUE
      }
    else {                        # two nodes
      x <- as.character(x)
      zk <- list(z[[x[1]]], z[[x[2]]])
      attr(zk, "members") <- attr(z[[x[1]]], "members") +
        attr(z[[x[2]]], "members")
      attr(zk, "midpoint") <- (attr(z[[x[1]]], "members") +
                               attr(z[[x[1]]], "midpoint") +
                               attr(z[[x[2]]], "midpoint"))/2
    }
    attr(zk, "height") <- oHgt[k]
    z[[k <- as.character(k+leafs)]] <- zk
  }
  z <- z[[k]]
  class(z) <- "dendrogram"
  z
}

#' @rdname communities
#' @importFrom stats as.hclust
#' @method as.hclust communities
#' @export
 
as.hclust.communities <- function(x, hang=-1, use.modularity=FALSE,
                                  ...) {
  as.hclust(as.dendrogram(x, hang=hang, use.modularity=use.modularity))
}

#' @rdname communities
#' @export

as_phylo <- function(x, ...)
  UseMethod("as_phylo")

#' @rdname communities
#' @method as_phylo communities
#' @export

as_phylo.communities <- function(x, use.modularity=FALSE, ...) {

  if (!is_hierarchical(x, full=TRUE)) {
    stop("Not a fully hierarchical community structure")
  }

  require(ape, quietly = TRUE)
  
  ## If multiple components, then we merge them in arbitrary order
  merges <- complete.dend(x, use.modularity)

  if (!use.modularity || is.null(x$modularity)) {
    height <- 1:nrow(merges)
  } else {
    height <- x$modularity[-1]
    height <- cumsum(height - min(height))
  }

  if (is.null(x$names)) {
    labels <- 1:(nrow(merges)+1)
  } else {
    labels <- x$names
  }

  N <- nrow(merges)
  edge <- matrix(0L, 2*N, 2)
  edge.length <- numeric(2*N)
  node <- integer(N)
  node[N] <- N + 2L
  cur.nod <- N + 3L
  j <- 1L
  for (i in N:1) {
    edge[j:(j+1), 1] <- node[i]
    for (l in 1:2) {
      k <- j + l -1L
      y <- merges[i, l]
      if (y > N+1) {
        edge[k, 2] <- node[y-N-1] <- cur.nod
        cur.nod <- cur.nod + 1L
        edge.length[k] <- height[i] - height[y-N-1]
      } else {
        edge[k, 2] <- y
        edge.length[k] <- height[i]
      }
    }
    j <- j + 2L    
  }

  obj <- list(edge=edge, edge.length=edge.length/2, tip.label=labels,
              Nnode=N)
  class(obj) <- "phylo"
  reorder(obj)
}

#' @rdname communities
#' @export

cut_at <- function(communities, no, steps) {

  if (!inherits(communities, "communities")) {
    stop("Not a community structure")
  }
  if (!is_hierarchical(communities)) {
    stop("Not a hierarchical communitity structure")
  }

  if ((!missing(no) && !missing(steps)) ||
      ( missing(no) &&  missing(steps))) {
    stop("Please give either `no' or `steps' (but not both)")
  }

  if (!missing(steps)) {
    mm <- merges(communities)
    if (steps > nrow(mm)) {
      warning("Cannot make that many steps")
      steps <- nrow(mm)
    }
    community.to.membership2(mm, communities$vcount, steps)
  } else {
    mm <- merges(communities)
    noc <- communities$vcount - nrow(mm) # final number of communities
    if (no<noc) {
      warning("Cannot have that few communities")
      no=noc
    }
    steps <- communities$vcount-no
    community.to.membership2(mm, communities$vcount, steps)    
  }
}

#' @rdname communities
#' @export

show_trace <- function(communities) {

  if (!inherits(communities, "communities")) {
    stop("Not a community structure")
  }
  if (is.null(communities$history)) {
    stop("History was not recorded")
  }

  res <- character()
  i <- 1
  while (i <= length(communities$history)) {
    if (communities$history[i] == 2) {  # IGRAPH_LEVC_HIST_SPLIT
      resnew <- paste("Splitting community", communities$history[i+1],
                      "into two.")
      i <- i + 2
    } else if (communities$history[i]==3) { # IGRAPH_LEVC_HIST_FAILED
      resnew <- paste("Failed splitting community",
                      communities$history[i+1], "into two.")
      i <- i + 2
    } else if (communities$history[i]==4) { # IGRAPH_LEVC_START_FULL
      resnew <- "Starting with the whole graph as a community."
      i <- i + 1
    } else if (communities$history[i]==5) { # IGRAPH_LEVC_START_GIVEN
      resnew <- paste("Starting from the", communities$history[i+1],
                      "given communities.")
      i <- i + 2
    }

    res <- c(res, resnew)
  }
  res
}

#####################################################################

community.to.membership2 <- function(merges, vcount, steps) {
  mode(merges) <- "numeric"
  mode(vcount) <- "numeric"
  mode(steps)  <- "numeric"
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  res <- .Call("R_igraph_community_to_membership2", merges-1, vcount, steps,
               PACKAGE="igraph")
  res+1
}

#####################################################################



#' Finding communities in graphs based on statistical meachanics
#' 
#' This function tries to find communities in graphs via a spin-glass model and
#' simulated annealing.
#' 
#' This function tries to find communities in a graph. A community is a set of
#' nodes with many edges inside the community and few edges between outside it
#' (i.e. between the community itself and the rest of the graph.)
#' 
#' This idea is reversed for edges having a negative weight, ie. few negative
#' edges inside a community and many negative edges between communities. Note
#' that only the \sQuote{neg} implementation supports negative edge weights.
#' 
#' The \code{spinglass.cummunity} function can solve two problems related to
#' community detection. If the \code{vertex} argument is not given (or it is
#' \code{NULL}), then the regular community detection problem is solved
#' (approximately), i.e. partitioning the vertices into communities, by
#' optimizing the an energy function.
#' 
#' If the \code{vertex} argument is given and it is not \code{NULL}, then it
#' must be a vertex id, and the same energy function is used to find the
#' community of the the given vertex. See also the examples below.
#'
#' @aliases spinglass.community
#' @param graph The input graph, can be directed but the direction of the edges
#' is neglected.
#' @param weights The weights of the edges. Either a numeric vector or
#' \code{NULL}. If it is null and the input graph has a \sQuote{weight} edge
#' attribute then that will be used. If \code{NULL} and no such attribute is
#' present then the edges will have equal weights. Set this to \code{NA} if the
#' graph was a \sQuote{weight} edge attribute, but you don't want to use it for
#' community detection.
#' @param vertex This parameter can be used to calculate the community of a
#' given vertex without calculating all communities. Note that if this argument
#' is present then some other arguments are ignored.
#' @param spins Integer constant, the number of spins to use. This is the upper
#' limit for the number of communities. It is not a problem to supply a
#' (reasonably) big number here, in which case some spin states will be
#' unpopulated.
#' @param parupdate Logical constant, whether to update the spins of the
#' vertices in parallel (synchronously) or not. This argument is ignored if the
#' second form of the function is used (ie. the \sQuote{\code{vertex}} argument
#' is present). It is also not implemented in the \dQuote{neg} implementation.
#' @param start.temp Real constant, the start temperature.  This argument is
#' ignored if the second form of the function is used (ie. the
#' \sQuote{\code{vertex}} argument is present).
#' @param stop.temp Real constant, the stop temperature. The simulation
#' terminates if the temperature lowers below this level.  This argument is
#' ignored if the second form of the function is used (ie. the
#' \sQuote{\code{vertex}} argument is present).
#' @param cool.fact Cooling factor for the simulated annealing.  This argument
#' is ignored if the second form of the function is used (ie. the
#' \sQuote{\code{vertex}} argument is present).
#' @param update.rule Character constant giving the \sQuote{null-model} of the
#' simulation. Possible values: \dQuote{simple} and \dQuote{config}.
#' \dQuote{simple} uses a random graph with the same number of edges as the
#' baseline probability and \dQuote{config} uses a random graph with the same
#' vertex degrees as the input graph.
#' @param gamma Real constant, the gamma argument of the algorithm. This
#' specifies the balance between the importance of present and non-present
#' edges in a community. Roughly, a comunity is a set of vertices having many
#' edges inside the community and few edges outside the community. The default
#' 1.0 value makes existing and non-existing links equally important. Smaller
#' values make the existing links, greater values the missing links more
#' important.
#' @param implementation Character scalar. Currently igraph contains two
#' implementations for the Spin-glass community finding algorithm. The faster
#' original implementation is the default. The other implementation, that takes
#' into account negative weights, can be chosen by supplying \sQuote{neg} here.
#' @param gamma.minus Real constant, the gamma.minus parameter of the
#' algorithm. This specifies the balance between the importance of present and
#' non-present negative weighted edges in a community. Smaller values of
#' gamma.minus, leads to communities with lesser negative intra-connectivity.
#' If this argument is set to zero, the algorithm reduces to a graph coloring
#' algorithm, using the number of spins as the number of colors. This argument
#' is ignored if the \sQuote{orig} implementation is chosen.
#' @return If the \code{vertex} argument is not given, ie. the first form is
#' used then a \code{\link{cluster_spinglass}} returns a
#' \code{\link{communities}} object.
#' 
#' If the \code{vertex} argument is present, ie. the second form is used then a
#' named list is returned with the following components:
#' \item{community}{Numeric vector giving the ids of the vertices in the same
#' community as \code{vertex}.} \item{cohesion}{The cohesion score of the
#' result, see references.} \item{adhesion}{The adhesion score of the result,
#' see references.} \item{inner.links}{The number of edges within the community
#' of \code{vertex}.} \item{outer.links}{The number of edges between the
#' community of \code{vertex} and the rest of the graph. }
#' @author Jorg Reichardt
#' (\url{http://theorie.physik.uni-wuerzburg.de/~reichardt/}) for the original
#' code and Gabor Csardi \email{csardi.gabor@@gmail.com} for the igraph glue
#' code.
#' 
#' Changes to the original function for including the possibility of negative
#' ties were implemented by Vincent Traag (\url{http://www.traag.net/}).
#' @seealso \code{\link{communities}}, \code{\link{components}}
#' @references J. Reichardt and S. Bornholdt: Statistical Mechanics of
#' Community Detection, \emph{Phys. Rev. E}, 74, 016110 (2006),
#' \url{http://arxiv.org/abs/cond-mat/0603718}
#' 
#' M. E. J. Newman and M. Girvan: Finding and evaluating community structure in
#' networks, \emph{Phys. Rev. E} 69, 026113 (2004)
#' 
#' V.A. Traag and Jeroen Bruggeman: Community detection in networks with
#' positive and negative links, \url{http://arxiv.org/abs/0811.2329} (2008).
#' @export
#' @keywords graphs
#' @examples
#' 
#'   g <- sample_gnp(10, 5/10) %du% sample_gnp(9, 5/9)
#'   g <- add_edges(g, c(1, 12))
#'   g <- induced_subgraph(g, subcomponent(g, 1))
#'   cluster_spinglass(g, spins=2)
#'   cluster_spinglass(g, vertex=1)
#' 
cluster_spinglass <- function(graph, weights=NULL, vertex=NULL, spins=25,
                                parupdate=FALSE, start.temp=1,
                                stop.temp=0.01, cool.fact=0.99,
                                update.rule=c("config", "random", "simple"),
                                gamma=1.0, implementation=c("orig", "neg"),
                                gamma.minus=1.0) {

  if (!is_igraph(graph)) {
    stop("Not a graph object")
  }

  if (is.null(weights) && "weight" %in% edge_attr_names(graph)) {
    weights <- E(graph)$weight
  }
  if (!is.null(weights) && any(!is.na(weights))) {
    weights <- as.numeric(weights)
  } else {
    weights <- NULL
  }

  update.rule <- igraph.match.arg(update.rule)
  update.rule <- switch(update.rule, "simple"=0, "random"=0, "config"=1)
  implementation <- switch(igraph.match.arg(implementation),
                                            "orig"=0, "neg"=1)

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  if (is.null(vertex)) {    
    res <- .Call("R_igraph_spinglass_community", graph, weights,
                 as.numeric(spins), as.logical(parupdate),
                 as.numeric(start.temp),
                 as.numeric(stop.temp), as.numeric(cool.fact),
                 as.numeric(update.rule), as.numeric(gamma),
                 as.numeric(implementation), as.numeric(gamma.minus),
                 PACKAGE="igraph")
    res$algorithm  <- "spinglass"
    res$vcount     <- vcount(graph)
    res$membership <- res$membership + 1
    if (igraph_opt("add.vertex.names") && is_named(graph)) {
      res$names <- vertex_attr(graph, "name")
    }
    class(res) <- "communities"
  } else {
    res <- .Call("R_igraph_spinglass_my_community", graph, weights,
                 as.igraph.vs(graph, vertex)-1, as.numeric(spins), 
                 as.numeric(update.rule), as.numeric(gamma),
                 PACKAGE="igraph")
    res$community <- res$community + 1
  }
  res
}



#' Community strucure via short random walks
#' 
#' This function tries to find densely connected subgraphs, also called
#' communities in a graph via random walks. The idea is that short random walks
#' tend to stay in the same community.
#' 
#' This function is the implementation of the Walktrap community finding
#' algorithm, see Pascal Pons, Matthieu Latapy: Computing communities in large
#' networks using random walks, http://arxiv.org/abs/physics/0512106
#'
#' @aliases walktrap.community
#' @param graph The input graph, edge directions are ignored in directed
#' graphs.
#' @param weights The edge weights.
#' @param steps The length of the random walks to perform.
#' @param merges Logical scalar, whether to include the merge matrix in the
#' result.
#' @param modularity Logical scalar, whether to include the vector of the
#' modularity scores in the result. If the \code{membership} argument is true,
#' then it will be always calculated.
#' @param membership Logical scalar, whether to calculate the membership vector
#' for the split corresponding to the highest modularity value.
#' @return \code{cluster_walktrap} returns a \code{\link{communities}}
#' object, please see the \code{\link{communities}} manual page for details.
#' @author Pascal Pons (\url{http://psl.pons.free.fr/}) and Gabor Csardi
#' \email{csardi.gabor@@gmail.com} for the R and igraph interface
#' @seealso See \code{\link{communities}} on getting the actual membership
#' vector, merge matrix, modularity score, etc.
#' 
#' \code{\link{modularity}} and \code{\link{cluster_fast_greedy}},
#' \code{\link{cluster_spinglass}},
#' \code{\link{cluster_leading_eigen}},
#' \code{\link{cluster_edge_betweenness}} for other community detection
#' methods.
#' @references Pascal Pons, Matthieu Latapy: Computing communities in large
#' networks using random walks, http://arxiv.org/abs/physics/0512106
#' @export
#' @keywords graphs
#' @examples
#' 
#' g <- make_full_graph(5) %du% make_full_graph(5) %du% make_full_graph(5)
#' g <- add_edges(g, c(1,6, 1,11, 6, 11))
#' cluster_walktrap(g)
#' 
cluster_walktrap <- function(graph, weights=E(graph)$weight, steps=4,
                               merges=TRUE, modularity=TRUE,
                               membership=TRUE) {
  if (!is_igraph(graph)) {
    stop("Not a graph object!")
  }

  if (membership && !modularity) {
    modularity <- TRUE
  }
  
  if (!is.null(weights)) {
    weights <- as.numeric(weights)
  }

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  res <- .Call("R_igraph_walktrap_community", graph, weights, as.numeric(steps),
        as.logical(merges), as.logical(modularity), as.logical(membership),
        PACKAGE="igraph")
  if (igraph_opt("add.vertex.names") && is_named(graph)) {
    res$names <- V(graph)$name
  }

  res$vcount <- vcount(graph)
  res$algorithm <- "walktrap"
  res$membership <- res$membership + 1
  res$merges <- res$merges + 1
  class(res) <- "communities"
  res
}



#' Community structure detection based on edge betweenness
#' 
#' Many networks consist of modules which are densely connected themselves but
#' sparsely connected to other modules.
#' 
#' The edge betweenness score of an edge measures the number of shortest paths
#' through it, see \code{\link{edge_betweenness}} for details. The idea of the
#' edge betweenness based community structure detection is that it is likely
#' that edges connecting separate modules have high edge betweenness as all the
#' shortest paths from one module to another must traverse through them. So if
#' we gradually remove the edge with the highest edge betweenness score we will
#' get a hierarchical map, a rooted tree, called a dendrogram of the graph. The
#' leafs of the tree are the individual vertices and the root of the tree
#' represents the whole graph.
#' 
#' \code{cluster_edge_betweenness} performs this algorithm by calculating the
#' edge betweenness of the graph, removing the edge with the highest edge
#' betweenness score, then recalculating edge betweenness of the edges and
#' again removing the one with the highest score, etc.
#' 
#' \code{edge.betweeness.community} returns various information collected
#' throught the run of the algorithm. See the return value down here.
#' 
#' @aliases edge.betweenness.community cluster_edge_betweenness
#' @param graph The graph to analyze.
#' @param weights The edge weights. Supply \code{NULL} to omit edge weights. By
#' default the \sQuote{\code{weight}} edge attribute is used, if it is present.
#' @param directed Logical constant, whether to calculate directed edge
#' betweenness for directed graphs. It is ignored for undirected graphs.
#' @param edge.betweenness Logical constant, whether to return the edge
#' betweenness of the edges at the time of their removal.
#' @param merges Logical constant, whether to return the merge matrix
#' representing the hierarchical community structure of the network.  This
#' argument is called \code{merges}, even if the community structure algorithm
#' itself is divisive and not agglomerative: it builds the tree from top to
#' bottom. There is one line for each merge (i.e. split) in matrix, the first
#' line is the first merge (last split). The communities are identified by
#' integer number starting from one. Community ids smaller than or equal to
#' \eqn{N}, the number of vertices in the graph, belong to singleton
#' communities, ie. individual vertices. Before the first merge we have \eqn{N}
#' communities numbered from one to \eqn{N}. The first merge, the first line of
#' the matrix creates community \eqn{N+1}, the second merge creates community
#' \eqn{N+2}, etc.
#' @param bridges Logical constant, whether to return a list the edge removals
#' which actually splitted a component of the graph.
#' @param modularity Logical constant, whether to calculate the maximum
#' modularity score, considering all possibly community structures along the
#' edge-betweenness based edge removals.
#' @param membership Logical constant, whether to calculate the membership
#' vector corresponding to the highest possible modularity score.
#' @return \code{cluster_edge_betweenness} returns a
#' \code{\link{communities}} object, please see the \code{\link{communities}}
#' manual page for details.
#' @author Gabor Csardi \email{csardi.gabor@@gmail.com}
#' @seealso \code{\link{edge_betweenness}} for the definition and calculation
#' of the edge betweenness, \code{\link{cluster_walktrap}},
#' \code{\link{cluster_fast_greedy}},
#' \code{\link{cluster_leading_eigen}} for other community detection
#' methods.
#' 
#' See \code{\link{communities}} for extracting the results of the community
#' detection.
#' @references M Newman and M Girvan: Finding and evaluating community
#' structure in networks, \emph{Physical Review E} 69, 026113 (2004)
#' @export
#' @keywords graphs
#' @examples
#' 
#' g <- barabasi.game(100,m=2)
#' eb <- cluster_edge_betweenness(g)
#' 
#' g <- make_full_graph(10) %du% make_full_graph(10)
#' g <- add_edges(g, c(1,11))
#' eb <- cluster_edge_betweenness(g)
#' eb
#' 
cluster_edge_betweenness <- function(graph, weights=E(graph)$weight,
                                       directed=TRUE,
                                       edge.betweenness=TRUE,
                                       merges=TRUE, bridges=TRUE,
                                       modularity=TRUE,
                                       membership=TRUE) {
  if (!is_igraph(graph)) {
    stop("Not a graph object!")
  }

  if (!is.null(weights)) {
    weights <- as.numeric(weights)
  }

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  res <- .Call("R_igraph_community_edge_betweenness", graph, weights,
               as.logical(directed),
               as.logical(edge.betweenness),
               as.logical(merges), as.logical(bridges),
               as.logical(modularity), as.logical(membership),
               PACKAGE="igraph")
  if (igraph_opt("add.vertex.names") && is_named(graph)) {
    res$names <- V(graph)$name
  }
  res$vcount <- vcount(graph)
  res$algorithm <- "edge betweenness"
  res$membership <- res$membership + 1
  res$merges <- res$merges + 1
  res$removed.edges <- res$removed.edges + 1
  res$bridges <- res$bridges + 1
  class(res) <- "communities"
  res
}

#' Community structure via greedy optimization of modularity
#' 
#' This function tries to find dense subgraph, also called communities in
#' graphs via directly optimizing a modularity score.
#' 
#' This function implements the fast greedy modularity optimization algorithm
#' for finding community structure, see A Clauset, MEJ Newman, C Moore: Finding
#' community structure in very large networks,
#' http://www.arxiv.org/abs/cond-mat/0408187 for the details.
#'
#' @aliases fastgreedy.community
#' @param graph The input graph
#' @param merges Logical scalar, whether to return the merge matrix.
#' @param modularity Logical scalar, whether to return a vector containing the
#' modularity after each merge.
#' @param membership Logical scalar, whether to calculate the membership vector
#' corresponding to the maximum modularity score, considering all possible
#' community structures along the merges.
#' @param weights If not \code{NULL}, then a numeric vector of edge weights.
#' The length must match the number of edges in the graph.  By default the
#' \sQuote{\code{weight}} edge attribute is used as weights. If it is not
#' present, then all edges are considered to have the same weight.
#' @return \code{cluster_fast_greedy} returns a \code{\link{communities}}
#' object, please see the \code{\link{communities}} manual page for details.
#' @author Tamas Nepusz \email{ntamas@@gmail.com} and Gabor Csardi
#' \email{csardi.gabor@@gmail.com} for the R interface.
#' @seealso \code{\link{communities}} for extracting the results.
#' 
#' See also \code{\link{cluster_walktrap}},
#' \code{\link{cluster_spinglass}},
#' \code{\link{cluster_leading_eigen}} and
#' \code{\link{cluster_edge_betweenness}} for other methods.
#' @references A Clauset, MEJ Newman, C Moore: Finding community structure in
#' very large networks, http://www.arxiv.org/abs/cond-mat/0408187
#' @export
#' @keywords graphs
#' @examples
#' 
#' g <- make_full_graph(5) %du% make_full_graph(5) %du% make_full_graph(5)
#' g <- add_edges(g, c(1,6, 1,11, 6, 11))
#' fc <- cluster_fast_greedy(g)
#' membership(fc)
#' sizes(fc)
#' 
cluster_fast_greedy <- function(graph, merges=TRUE, modularity=TRUE,
                                 membership=TRUE,
                                 weights=E(graph)$weight) {
  if (!is_igraph(graph)) {
    stop("Not a graph object")
  }

  if (!is.null(weights)) {
    weights <- as.numeric(weights)
  }

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  res <- .Call("R_igraph_community_fastgreedy", graph, as.logical(merges),
               as.logical(modularity), as.logical(membership), weights,
               PACKAGE="igraph")
  if (igraph_opt("add.vertex.names") && is_named(graph)) {
    res$names <- V(graph)$name
  }
  res$algorithm <- "fast greedy"
  res$vcount <- vcount(graph)
  res$membership <- res$membership + 1
  res$merges <- res$merges + 1
  class(res) <- "communities"
  res
}

igraph.i.levc.arp <- function(externalP, externalE) {
  f <- function(v) {
    v <- as.numeric(v)
    base::.Call("R_igraph_i_levc_arp", externalP, externalE, v,
                PACKAGE="igraph");
  }
  f
}



#' Community structure detecting based on the leading eigenvector of the
#' community matrix
#' 
#' This function tries to find densely connected subgraphs in a graph by
#' calculating the leading non-negative eigenvector of the modularity matrix of
#' the graph.
#' 
#' The function documented in these section implements the \sQuote{leading
#' eigenvector} method developed by Mark Newman, see the reference below.
#' 
#' The heart of the method is the definition of the modularity matrix,
#' \code{B}, which is \code{B=A-P}, \code{A} being the adjacency matrix of the
#' (undirected) network, and \code{P} contains the probability that certain
#' edges are present according to the \sQuote{configuration model}. In other
#' words, a \code{P[i,j]} element of \code{P} is the probability that there is
#' an edge between vertices \code{i} and \code{j} in a random network in which
#' the degrees of all vertices are the same as in the input graph.
#' 
#' The leading eigenvector method works by calculating the eigenvector of the
#' modularity matrix for the largest positive eigenvalue and then separating
#' vertices into two community based on the sign of the corresponding element
#' in the eigenvector. If all elements in the eigenvector are of the same sign
#' that means that the network has no underlying comuunity structure.  Check
#' Newman's paper to understand why this is a good method for detecting
#' community structure.
#' 
#' @aliases leading.eigenvector.community
#' @param graph The input graph. Should be undirected as the method needs a
#' symmetric matrix.
#' @param steps The number of steps to take, this is actually the number of
#' tries to make a step. It is not a particularly useful parameter.
#' #' @param weights An optional weight vector. The \sQuote{weight} edge attribute
#' is used if present. Supply \sQuote{\code{NA}} here if you want to ignore the
#' \sQuote{weight} edge attribute.
#' @param start \code{NULL}, or a numeric membership vector, giving the start
#' configuration of the algorithm.
#' @param options A named list to override some ARPACK options.
#' @param callback If not \code{NULL}, then it must be callback function. This
#' is called after each iteration, after calculating the leading eigenvector of
#' the modularity matrix. See details below.
#' @param extra Additional argument to supply to the callback function.
#' @param env The environment in which the callback function is evaluated.
#' @return \code{cluster_leading_eigen} returns a named list with the
#' following members: \item{membership}{The membership vector at the end of the
#' algorithm, when no more splits are possible.} \item{merges}{The merges
#' matrix starting from the state described by the \code{membership} member.
#' This is a two-column matrix and each line describes a merge of two
#' communities, the first line is the first merge and it creates community
#' \sQuote{\code{N}}, \code{N} is the number of initial communities in the
#' graph, the second line creates community \code{N+1}, etc.  }
#' \item{options}{Information about the underlying ARPACK computation, see
#' \code{\link{arpack}} for details.  }
#' @section Callback functions: The \code{callback} argument can be used to
#' supply a function that is called after each eigenvector calculation. The
#' following arguments are supplied to this function: \describe{
#'   \item{membership}{The actual membership vector, with zero-based indexing.}
#'   \item{community}{The community that the algorithm just tried to split,
#'     community numbering starts with zero here.}
#'   \item{value}{The eigenvalue belonging to the leading eigenvector the
#'     algorithm just found.}
#'   \item{vector}{The leading eigenvector the algorithm just found.}
#'   \item{multiplier}{An R function that can be used to multiple the actual
#'     modularity matrix with an arbitrary vector. Supply the vector as an
#'     argument to perform this multiplication. This function can be used
#'     with ARPACK.}
#'   \item{extra}{The \code{extra} argument that was passed to
#'     \code{cluster_leading_eigen}. }
#'   The callback function should return a scalar number. If this number
#'   is non-zero, then the clustering is terminated.
#' }
#' @author Gabor Csardi \email{csardi.gabor@@gmail.com}
#' @seealso \code{\link{modularity}}, \code{\link{cluster_walktrap}},
#' \code{\link{cluster_edge_betweenness}},
#' \code{\link{cluster_fast_greedy}}, \code{\link[stats]{as.dendrogram}}
#' @references MEJ Newman: Finding community structure using the eigenvectors
#' of matrices, Physical Review E 74 036104, 2006.
#' @export
#' @keywords graphs
#' @examples
#' 
#' g <- make_full_graph(5) %du% make_full_graph(5) %du% make_full_graph(5)
#' g <- add_edges(g, c(1,6, 1,11, 6, 11))
#' lec <- cluster_leading_eigen(g)
#' lec
#' 
#' cluster_leading_eigen(g, start=membership(lec))
#' 
cluster_leading_eigen <- function(graph, steps=-1, weights=NULL,
                                          start=NULL,
                                          options=arpack_defaults,
                                          callback=NULL, extra=NULL,
                                          env=parent.frame()){

  # Argument checks
  if (!is_igraph(graph)) { stop("Not a graph object") }
  steps <- as.integer(steps)
  if (is.null(weights) && "weight" %in% edge_attr_names(graph)) { 
    weights <- E(graph)$weight 
  } 
  if (!is.null(weights) && any(!is.na(weights))) { 
    weights <- as.numeric(weights) 
  } else { 
    weights <- NULL 
  }
  if (!is.null(start)) { start <- as.numeric(start)-1 }
  options.tmp <- arpack_defaults; options.tmp[ names(options) ] <- options ; options <- options.tmp
  
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  # Function call
  res <- .Call("R_igraph_community_leading_eigenvector", graph, steps,
               weights, options, start, callback, extra, env,
               environment(igraph.i.levc.arp),
               PACKAGE="igraph")
  if (igraph_opt("add.vertex.names") && is_named(graph)) {
    res$names <- V(graph)$name
  }
  res$algorithm <- "leading eigenvector"
  res$vcount <- vcount(graph)
  res$membership <- res$membership + 1
  res$merges <- res$merges + 1
  res$history <- res$history + 1
  class(res) <- "communities"
  res
}

#' Finding communities based on propagating labels
#' 
#' This is a fast, nearly linear time algorithm for detecting community
#' structure in networks. In works by labeling the vertices with unique labels
#' and then updating the labels by majority voting in the neighborhood of the
#' vertex.
#' 
#' This function implements the community detection method described in:
#' Raghavan, U.N. and Albert, R. and Kumara, S.: Near linear time algorithm to
#' detect community structures in large-scale networks. Phys Rev E 76, 036106.
#' (2007). This version extends the original method by the ability to take edge
#' weights into consideration and also by allowing some labels to be fixed.
#' 
#' From the abstract of the paper: \dQuote{In our algorithm every node is
#' initialized with a unique label and at every step each node adopts the label
#' that most of its neighbors currently have. In this iterative process densely
#' connected groups of nodes form a consensus on a unique label to form
#' communities.}
#'
#' @aliases label.propagation.community
#' @param graph The input graph, should be undirected to make sense.
#' @param weights An optional weight vector. It should contain a positive
#' weight for all the edges. The \sQuote{weight} edge attribute is used if
#' present. Supply \sQuote{\code{NA}} here if you want to ignore the
#' \sQuote{weight} edge attribute.
#' @param initial The initial state. If \code{NULL}, every vertex will have a
#' different label at the beginning. Otherwise it must be a vector with an
#' entry for each vertex. Non-negative values denote different labels, negative
#' entries denote vertices without labels.
#' @param fixed Logical vector denoting which labels are fixed. Of course this
#' makes sense only if you provided an initial state, otherwise this element
#' will be ignored. Also note that vertices without labels cannot be fixed.
#' @return \code{cluster_label_prop} returns a
#' \code{\link{communities}} object, please see the \code{\link{communities}}
#' manual page for details.
#' @author Tamas Nepusz \email{ntamas@@gmail.com} for the C implementation,
#' Gabor Csardi \email{csardi.gabor@@gmail.com} for this manual page.
#' @seealso \code{\link{communities}} for extracting the actual results.
#' 
#' \code{\link{cluster_fast_greedy}}, \code{\link{cluster_walktrap}} and
#' \code{\link{cluster_spinglass}} for other community detection methods.
#' @references Raghavan, U.N. and Albert, R. and Kumara, S.: Near linear time
#' algorithm to detect community structures in large-scale networks. \emph{Phys
#' Rev E} 76, 036106. (2007)
#' @export
#' @keywords graphs
#' @examples
#' 
#'   g <- sample_gnp(10, 5/10) %du% sample_gnp(9, 5/9)
#'   g <- add_edges(g, c(1, 12))
#'   cluster_label_prop(g)
#' 
cluster_label_prop <- function(graph, weights=NULL, initial=NULL,
                                        fixed=NULL) {
  # Argument checks
  if (!is_igraph(graph)) { stop("Not a graph object") }
  if (is.null(weights) && "weight" %in% edge_attr_names(graph)) { 
  weights <- E(graph)$weight 
  } 
  if (!is.null(weights) && any(!is.na(weights))) { 
  weights <- as.numeric(weights) 
  } else { 
  weights <- NULL 
  }
  if (!is.null(initial)) initial <- as.numeric(initial)
  if (!is.null(fixed)) fixed <- as.logical(fixed)

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  # Function call
  res <- .Call("R_igraph_community_label_propagation", graph, weights, initial, fixed,
        PACKAGE="igraph")
  if (igraph_opt("add.vertex.names") && is_named(graph)) {
    res$names <- V(graph)$name
  }
  res$vcount <- vcount(graph)
  res$algorithm <- "label propagation"
  res$membership <- res$membership + 1
  class(res) <- "communities"
  res
}



#' Finding community structure by multi-level optimization of modularity
#' 
#' This function implements the multi-level modularity optimization algorithm
#' for finding community structure, see references below. It is based on the
#' modularity measure and a hierarchial approach.
#' 
#' This function implements the multi-level modularity optimization algorithm
#' for finding community structure, see VD Blondel, J-L Guillaume, R Lambiotte
#' and E Lefebvre: Fast unfolding of community hierarchies in large networks,
#' \url{http://arxiv.org/abs/arXiv:0803.0476} for the details.
#' 
#' It is based on the modularity measure and a hierarchial approach.
#' Initially, each vertex is assigned to a community on its own. In every step,
#' vertices are re-assigned to communities in a local, greedy way: each vertex
#' is moved to the community with which it achieves the highest contribution to
#' modularity. When no vertices can be reassigned, each community is considered
#' a vertex on its own, and the process starts again with the merged
#' communities. The process stops when there is only a single vertex left or
#' when the modularity cannot be increased any more in a step.
#' 
#' This function was contributed by Tom Gregorovic.
#'
#' @aliases multilevel.community
#' @param graph The input graph.
#' @param weights Optional positive weight vector.  If the graph has a
#' \code{weight} edge attribute, then this is used by default. Supply \code{NA}
#' here if the graph has a \code{weight} edge attribute, but you want to ignore
#' it.
#' @return \code{cluster_louvain} returns a \code{\link{communities}}
#' object, please see the \code{\link{communities}} manual page for details.
#' @author Tom Gregorovic, Tamas Nepusz \email{ntamas@@gmail.com}
#' @seealso See \code{\link{communities}} for extracting the membership,
#' modularity scores, etc. from the results.
#' 
#' Other community detection algorithms: \code{\link{cluster_walktrap}},
#' \code{\link{cluster_spinglass}},
#' \code{\link{cluster_leading_eigen}},
#' \code{\link{cluster_edge_betweenness}},
#' \code{\link{cluster_fast_greedy}},
#' \code{\link{cluster_label_prop}}
#' @references Vincent D. Blondel, Jean-Loup Guillaume, Renaud Lambiotte,
#' Etienne Lefebvre: Fast unfolding of communities in large networks. J. Stat.
#' Mech. (2008) P10008
#' @export
#' @keywords graphs
#' @examples
#' 
#' # This is so simple that we will have only one level
#' g <- make_full_graph(5) %du% make_full_graph(5) %du% make_full_graph(5)
#' g <- add_edges(g, c(1,6, 1,11, 6, 11))
#' cluster_louvain(g)
#' 
cluster_louvain <- function(graph, weights=NULL) {
  # Argument checks
  if (!is_igraph(graph)) { stop("Not a graph object") }
  if (is.null(weights) && "weight" %in% edge_attr_names(graph)) { 
  weights <- E(graph)$weight 
  } 
  if (!is.null(weights) && any(!is.na(weights))) { 
  weights <- as.numeric(weights) 
  } else { 
  weights <- NULL 
  }

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  # Function call
  res <- .Call("R_igraph_community_multilevel", graph, weights,
        PACKAGE="igraph")
  if (igraph_opt("add.vertex.names") && is_named(graph)) {
    res$names <- V(graph)$name
  }
  res$vcount <- vcount(graph)
  res$algorithm <- "multi level"
  res$membership <- res$membership + 1
  res$memberships <- res$memberships + 1
  class(res) <- "communities"
  res
}



#' Optimal community structure
#' 
#' This function calculates the optimal community structure of a graph, by
#' maximizing the modularity measure over all possible partitions.
#' 
#' This function calculates the optimal community structure for a graph, in
#' terms of maximal modularity score.
#' 
#' The calculation is done by transforming the modularity maximization into an
#' integer programming problem, and then calling the GLPK library to solve
#' that. Please the reference below for details.
#' 
#' Note that modularity optimization is an NP-complete problem, and all known
#' algorithms for it have exponential time complexity. This means that you
#' probably don't want to run this function on larger graphs. Graphs with up to
#' fifty vertices should be fine, graphs with a couple of hundred vertices
#' might be possible.
#'
#' @aliases optimal.community
#' @param graph The input graph. Edge directions are ignored for directed
#' graphs.
#' @param weights Optional positive weight vector for optimizing weighted
#' modularity. If the graph has a \code{weight} edge attribute, then this is
#' used by default. Supply \code{NA} to ignore the weights of a weighted graph.
#' @return \code{cluster_optimal} returns a \code{\link{communities}} object,
#' please see the \code{\link{communities}} manual page for details.
#' @author Gabor Csardi \email{csardi.gabor@@gmail.com}
#' @seealso \code{\link{communities}} for the documentation of the result,
#' \code{\link{modularity}}. See also \code{\link{cluster_fast_greedy}} for a
#' fast greedy optimizer.
#' @references Ulrik Brandes, Daniel Delling, Marco Gaertler, Robert Gorke,
#' Martin Hoefer, Zoran Nikoloski, Dorothea Wagner: On Modularity Clustering,
#' \emph{IEEE Transactions on Knowledge and Data Engineering} 20(2):172-188,
#' 2008.
#' @export
#' @keywords graphs
#' @examples
#' 
#' ## Zachary's karate club
#' g <- make_graph("Zachary")
#' 
#' ## We put everything into a big 'try' block, in case 
#' ## igraph was compiled without GLPK support
#' 
#' try({
#'   ## The calculation only takes a couple of seconds
#'   oc <- cluster_optimal(g)
#' 
#'   ## Double check the result
#'   print(modularity(oc))
#'   print(modularity(g, membership(oc)))
#' 
#'   ## Compare to the greedy optimizer
#'   fc <- cluster_fast_greedy(g)
#'   print(modularity(fc))
#' }, silent=TRUE)
#' 
cluster_optimal <- function(graph, weights=NULL) {
  # Argument checks
  if (!is_igraph(graph)) { stop("Not a graph object") }
  if (is.null(weights) && "weight" %in% edge_attr_names(graph)) {
    weights <- E(graph)$weight
  }
  if (!is.null(weights) && any(!is.na(weights))) {
    weights <- as.numeric(weights)
  } else {
    weights <- NULL
  }

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  # Function call
  res <- .Call("R_igraph_community_optimal_modularity", graph, weights,
               PACKAGE="igraph")
  if (igraph_opt("add.vertex.names") && is_named(graph)) {
    res$names <- V(graph)$name
  }
  res$vcount <- vcount(graph)
  res$algorithm <- "optimal"
  res$membership <- res$membership + 1
  class(res) <- "communities"
  res
}



#' Infomap community finding
#' 
#' Find community structure that minimizes the expected description length of a
#' random walker trajectory
#' 
#' Please see the details of this method in the references given below.
#'
#' @aliases infomap.community
#' @param graph The input graph.
#' @param e.weights If not \code{NULL}, then a numeric vector of edge weights.
#' The length must match the number of edges in the graph.  By default the
#' \sQuote{\code{weight}} edge attribute is used as weights. If it is not
#' present, then all edges are considered to have the same weight.
#' @param v.weights If not \code{NULL}, then a numeric vector of vertex
#' weights. The length must match the number of vertices in the graph.  By
#' default the \sQuote{\code{weight}} vertex attribute is used as weights. If
#' it is not present, then all vertices are considered to have the same weight.
#' @param nb.trials The number of attempts to partition the network (can be any
#' integer value equal or larger than 1).
#' @param modularity Logical scalar, whether to calculate the modularity score
#' of the detected community structure.
#' @return \code{cluster_infomap} returns a \code{\link{communities}} object,
#' please see the \code{\link{communities}} manual page for details.
#' @author Martin Rosvall (\url{http://www.tp.umu.se/~rosvall/}) wrote the
#' original C++ code. This was ported to be more igraph-like by Emmanuel
#' Navarro (\url{http://www.irit.fr/~Emmanuel.Navarro/}).  The R interface and
#' some cosmetics was done by Gabor Csardi \email{csardi.gabor@@gmail.com}.
#' @seealso Other community finding methods and \code{\link{communities}}.
#' @references The original paper: M. Rosvall and C. T. Bergstrom, Maps of
#' information flow reveal community structure in complex networks, \emph{PNAS}
#' 105, 1118 (2008) \url{http://dx.doi.org/10.1073/pnas.0706851105},
#' \url{http://arxiv.org/abs/0707.0609}
#' 
#' A more detailed paper: M. Rosvall, D. Axelsson, and C. T. Bergstrom, The map
#' equation, \emph{Eur. Phys. J. Special Topics} 178, 13 (2009).
#' \url{http://dx.doi.org/10.1140/epjst/e2010-01179-1},
#' \url{http://arxiv.org/abs/0906.1405}.
#' @export
#' @keywords graphs
#' @examples
#' 
#' ## Zachary's karate club
#' g <- make_graph("Zachary")
#' 
#' imc <- cluster_infomap(g)
#' membership(imc)
#' communities(imc)
#' 
cluster_infomap <- function(graph, e.weights=NULL, v.weights=NULL,
                              nb.trials=10, modularity=TRUE) {
  
  # Argument checks
  if (!is_igraph(graph)) { stop("Not a graph object") }
  if (is.null(e.weights) && "weight" %in% edge_attr_names(graph)) { 
    e.weights <- E(graph)$weight 
  } 
  if (!is.null(e.weights) && any(!is.na(e.weights))) { 
    e.weights <- as.numeric(e.weights) 
  } else { 
    e.weights <- NULL 
  }
  if (is.null(v.weights) && "weight" %in% vertex_attr_names(graph)) { 
    v.weights <- V(graph)$weight 
  } 
  if (!is.null(v.weights) && any(!is.na(v.weights))) { 
    v.weights <- as.numeric(v.weights) 
  } else { 
    v.weights <- NULL 
  }
  nb.trials <- as.integer(nb.trials)
  
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  # Function call
  res <- .Call("R_igraph_community_infomap", graph, e.weights,
               v.weights, nb.trials,
               PACKAGE="igraph")

  if (igraph_opt("add.vertex.names") && is_named(graph)) {
    res$names <- V(graph)$name
  }
  res$vcount <- vcount(graph)
  res$algorithm <- "infomap"
  res$membership <- res$membership + 1
  if (modularity) {
    res$modularity <- modularity(graph, res$membership, weights=e.weights)
  }
  class(res) <- "communities"
  res
}

#' @rdname communities
#' @method plot communities
#' @export

plot.communities <- function(x, y,
                             colbar=rainbow(length(x)),
                             col=colbar[membership(x)],
                             mark.groups=communities(x),
                             edge.color=c("black", "red")[crossing(x,y)+1],
                             ...) {

  plot(y, vertex.color=col, mark.groups=mark.groups,
       edge.color=edge.color,
       ...)  
}



#' @rdname plot_dendrogram.communities
#' @aliases dendPlot
#' @export

plot_dendrogram <- function(x, mode=igraph_opt("dend.plot.type"), ...)
  UseMethod("plot_dendrogram")



#' Community structure dendrogram plots
#' 
#' Plot a hierarchical community structure as a dendrogram.
#' 
#' \code{plot_dendrogram} supports three different plotting functions, selected via
#' the \code{mode} argument. By default the plotting function is taken from the
#' \code{dend.plot.type} igraph option, and it has for possible values:
#' \itemize{ \item \code{auto} Choose automatically between the plotting
#' functions. As \code{plot.phylo} is the most sophisticated, that is choosen,
#' whenever the \code{ape} package is available. Otherwise \code{plot.hclust}
#' is used.  \item \code{phylo} Use \code{plot.phylo} from the \code{ape}
#' package.  \item \code{hclust} Use \code{plot.hclust} from the \code{stats}
#' package.  \item \code{dendrogram} Use \code{plot.dendrogram} from the
#' \code{stats} package.  }
#' 
#' The different plotting functions take different sets of arguments. When
#' using \code{plot.phylo} (\code{mode="phylo"}), we have the following syntax:
#' \preformatted{
#'     plot_dendrogram(x, mode="phylo",
#'             colbar = rainbow(11, start=0.7, end=0.1),
#'             edge.color = NULL, use.edge.length = FALSE, \dots)
#' } The extra arguments not documented above: \itemize{
#'   \item \code{colbar} Color bar for the edges.
#'   \item \code{edge.color} Edge colors. If \code{NULL}, then the
#'     \code{colbar} argument is used.
#'   \item \code{use.edge.length} Passed to \code{plot.phylo}.
#'   \item \code{dots} Attitional arguments to pass to \code{plot.phylo}.
#' }
#' 
#' The syntax for \code{plot.hclust} (\code{mode="hclust"}): \preformatted{
#'     plot_dendrogram(x, mode="hclust", rect = 0, colbar = rainbow(rect),
#'             hang = 0.01, ann = FALSE, main = "", sub = "", xlab = "",
#'             ylab = "", \dots)
#' } The extra arguments not documented above: \itemize{
#'   \item \code{rect} A numeric scalar, the number of groups to mark on
#'     the dendrogram. The dendrogram is cut into exactly \code{rect}
#'     groups and they are marked via the \code{rect.hclust} command. Set
#'     this to zero if you don't want to mark any groups.
#'   \item \code{colbar} The colors of the rectanges that mark the
#'     vertex groups via the \code{rect} argument.
#'   \item \code{hang} Where to put the leaf nodes, this corresponds to the
#'     \code{hang} argument of \code{plot.hclust}.
#'   \item \code{ann}  Whether to annotate the plot, the \code{ann}
#'     argument of \code{plot.hclust}.
#'   \item \code{main} The main title of the plot, the \code{main} argument
#'     of \code{plot.hclust}.
#'   \item \code{sub} The sub-title of the plot, the \code{sub} argument of
#'     \code{plot.hclust}.
#'   \item \code{xlab} The label on the horizontal axis, passed to
#'     \code{plot.hclust}.
#'   \item \code{ylab} The label on the vertical axis, passed to
#'     \code{plot.hclust}.
#'   \item \code{dots} Attitional arguments to pass to \code{plot.hclust}.
#' }
#' 
#' The syntax for \code{plot.dendrogram} (\code{mode="dendrogram"}):
#' \preformatted{
#'     plot_dendrogram(x, \dots)
#' } The extra arguments are simply passed to \code{as.dendrogram}.
#' 
#' @param x An object containing the community structure of a graph. See
#' \code{\link{communities}} for details.
#' @param mode Which dendrogram plotting function to use. See details below.
#' @param \dots Additional arguments to supply to the dendrogram plotting
#' function.
#' @param use.modularity Logical scalar, whether to use the modularity values
#' to define the height of the branches.
#' @return Returns whatever the return value was from the plotting function,
#' \code{plot.phylo}, \code{plot.dendrogram} or \code{plot.hclust}.
#' @author Gabor Csardi \email{csardi.gabor@@gmail.com}
#' @method plot_dendrogram communities
#' @export
#' @keywords graphs
#' @examples
#' 
#' karate <- make_graph("Zachary")
#' fc <- cluster_fast_greedy(karate)
#' plot_dendrogram(fc)
#' 
plot_dendrogram.communities <- function(x, 
                                 mode=igraph_opt("dend.plot.type"), ...,
                                 use.modularity=FALSE) {  
  mode <- igraph.match.arg(mode, c("auto", "phylo", "hclust", "dendrogram"))

  if (mode=="auto") {
    value <- tryCatch(suppressWarnings(library("ape", character.only=TRUE,
                                               logical.return=TRUE,
                                               warn.conflicts=FALSE,
                                               quietly=TRUE,
                                               pos="package:base")),
                      error=function(e) e)
    mode <- if (value) "phylo" else "hclust"
  }
  
  if (mode=="hclust") {
    dendPlotHclust(x, use.modularity=use.modularity, ...)
  } else if (mode=="dendrogram") {
    dendPlotDendrogram(x, use.modularity=use.modularity, ...)
  } else if (mode=="phylo") {
    dendPlotPhylo(x, use.modularity=use.modularity, ...)
  }
}

dendPlotHclust <- function(communities, rect=length(communities),
                           colbar=rainbow(rect), hang=-1, ann=FALSE,
                           main="", sub="", xlab="", ylab="", ...,
                           use.modularity=FALSE) {
  hc <- as.hclust(communities, hang=hang, use.modularity=use.modularity)
  ret <- plot(hc, hang=hang, ann=ann, main=main, sub=sub, xlab=xlab,
              ylab=ylab, ...)
  if (rect > 0) {
    rect.hclust(hc, k=rect, border=colbar)
  }
  invisible(ret)
}

dendPlotDendrogram <- function(communities, hang=-1, ...,
                               use.modularity=FALSE) {
  plot(as.dendrogram(communities, hang=hang, use.modularity=use.modularity),
       ...)
}

dendPlotPhylo <- function(communities, colbar=rainbow(length(communities)),
                          col=colbar[membership(communities)],
                          mark.groups=communities(communities),
                          use.modularity=FALSE, 
                          edge.color="#AAAAAAFF",
                          edge.lty=c(1,2), ...) {
  
  phy <- as_phylo(communities, use.modularity=use.modularity)

  getedges <- function(tip) {
    repeat {      
      ee <- which(! phy$edge[,1] %in% tip & phy$edge[,2] %in% tip)
      if (length(ee)<=1) { break }
      tip <- c(tip, unique(phy$edge[ee,1]))
    }
    ed <- which(phy$edge[,1] %in% tip & phy$edge[,2] %in% tip)
    eds <- phy$edge[ed, 1]
    good <- which(phy$edge[ed,1] %in% which(tabulate(eds) != 1))
    ed[good]
  }
  gredges <- lapply(mark.groups, getedges)

  if (length(mark.groups) > 0) {
    ecol <- rep(edge.color, nrow(phy$edge))
    for (gr in seq_along(gredges)) {
      ecol[gredges[[gr]]] <- colbar[gr]
    }
  } else {
    ecol <- edge.color
  }
  
  elty <- rep(edge.lty[2], nrow(phy$edge))
  elty[ unlist(gredges) ] <- edge.lty[1]
  
  plot(phy, edge.color=ecol, edge.lty=elty, tip.color=col, ...)
}

#' Compares community structures using various metrics
#' 
#' This function assesses the distance between two community structures.
#' 
#' 
#' @aliases compare.communities compare.membership compare
#' @param comm1 A \code{\link{communities}} object containing a community
#' structure; or a numeric vector, the membership vector of the first community
#' structure. The membership vector should contain the community id of each
#' vertex, the numbering of the communities starts with one.
#' @param comm2 A \code{\link{communities}} object containing a community
#' structure; or a numeric vector, the membership vector of the second
#' community structure, in the same format as for the previous argument.
#' @param method Character scalar, the comparison method to use. Possible
#' values: \sQuote{vi} is the variation of information (VI) metric of Meila
#' (2003), \sQuote{nmi} is the normalized mutual information measure proposed
#' by Danon et al. (2005), \sQuote{split.join} is the split-join distance of
#' can Dongen (2000), \sQuote{rand} is the Rand index of Rand (1971),
#' \sQuote{adjusted.rand} is the adjusted Rand index by Hubert and Arabie
#' (1985).
#' @return A real number.
#' @author Tamas Nepusz \email{ntamas@@gmail.com}
#' @seealso \code{\link{cluster_walktrap}},
#' \code{\link{cluster_edge_betweenness}},
#' \code{\link{cluster_fast_greedy}}, \code{\link{cluster_spinglass}} for
#' various community detection methods.
#' @references Meila M: Comparing clusterings by the variation of information.
#' In: Scholkopf B, Warmuth MK (eds.). \emph{Learning Theory and Kernel
#' Machines: 16th Annual Conference on Computational Learning Theory and 7th
#' Kernel Workshop}, COLT/Kernel 2003, Washington, DC, USA. Lecture Notes in
#' Computer Science, vol. 2777, Springer, 2003. ISBN: 978-3-540-40720-1.
#' 
#' Danon L, Diaz-Guilera A, Duch J, Arenas A: Comparing community structure
#' identification. \emph{J Stat Mech} P09008, 2005.
#' 
#' van Dongen S: Performance criteria for graph clustering and Markov cluster
#' experiments. Technical Report INS-R0012, National Research Institute for
#' Mathematics and Computer Science in the Netherlands, Amsterdam, May 2000.
#' 
#' Rand WM: Objective criteria for the evaluation of clustering methods.
#' \emph{J Am Stat Assoc} 66(336):846-850, 1971.
#' 
#' Hubert L and Arabie P: Comparing partitions. \emph{Journal of
#' Classification} 2:193-218, 1985.
#' @export
#' @keywords graphs
#' @examples
#' 
#' g <- make_graph("Zachary")
#' sg <- cluster_spinglass(g)
#' le <- cluster_leading_eigen(g)
#' compare(sg, le, method="rand")
#' compare(membership(sg), membership(le))
#' 
compare <- function(comm1, comm2, method=c("vi", "nmi",
                                    "split.join", "rand",
                                    "adjusted.rand"))
  UseMethod("compare")

#' @method compare communities
#' @export

compare.communities <- function(comm1, comm2,
                                method=c("vi", "nmi", "split.join", "rand",
                                  "adjusted.rand")) {

  i_compare(comm1, comm2, method)
}

#' @method compare membership
#' @export

compare.membership <- function(comm1, comm2,
                                method=c("vi", "nmi", "split.join", "rand",
                                  "adjusted.rand")) {

  i_compare(comm1, comm2, method)
}

#' @method compare default
#' @export

compare.default <- compare.membership

i_compare <- function (comm1, comm2, method=c("vi", "nmi", "split.join",
                                       "rand", "adjusted.rand")) {

  comm1 <- if (inherits(comm1, "communities")) {
    as.numeric(membership(comm1))
  } else {
    as.numeric(comm1)
  }
  comm2 <- if (inherits(comm2, "communities")) {
    as.numeric(membership(comm2))
  } else {
    as.numeric(comm2)
  }
  method <- switch(igraph.match.arg(method), vi = 0, nmi = 1, 
                   split.join = 2, rand = 3, adjusted.rand = 4)
  on.exit(.Call("R_igraph_finalizer", PACKAGE = "igraph"))
  res <- .Call("R_igraph_compare_communities", comm1, comm2, 
               method, PACKAGE = "igraph")
  res  
}


#' @export

split_join_distance <- function(comm1, comm2) {
  comm1 <- if (inherits(comm1, "communities")) {
    as.numeric(membership(comm1))
  } else {
    as.numeric(comm1)
  }
  comm2 <- if (inherits(comm2, "communities")) {
    as.numeric(membership(comm2))
  } else {
    as.numeric(comm2)
  }
  on.exit(.Call("R_igraph_finalizer", PACKAGE = "igraph"))
  res <- .Call("R_igraph_split_join_distance", comm1, comm2,
               PACKAGE = "igraph")
  unlist(res)
}

#' @export

groups <- function(x)
  UseMethod("groups")

#' @method groups default
#' @export

groups.default <- function(x) {
  vids <- names(x$membership)
  if (is.null(vids)) vids <- seq_along(x$membership)
  tapply(vids, x$membership, simplify=FALSE, function(x) x)
}

#' @method groups communities
#' @export

groups.communities <- function(x) {
  m <- membership(x)
  groups.default(list(membership = m))
}

#' @export

communities <- groups.communities

#' @method "[" communities
#' @export

`[.communities` <- function(x, i) {
  groups(x)[i]
}

#' @method "[[" communities
#' @export

`[[.communities` <- function(x, i) {
  groups(x)[[i]]
}
