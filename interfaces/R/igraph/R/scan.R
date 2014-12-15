## -----------------------------------------------------------------------
##
##   IGraph R package
##   Copyright (C) 2014  Gabor Csardi <csardi.gabor@gmail.com>
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
## -----------------------------------------------------------------------

#' Compute local scan statistics on graphs
#'
#' The scan statistic is a summary of the locality statistics that is
#' computed from the local neighborhood of each vertex. The
#' \code{local_scan} function computes the local statistics for each vertex
#' for a given neighborhood size and the statistic function.
#'
#' See the given reference below for the details on the local scan
#' statistics.
#'
#' \code{local_scan} calculates exact local scan statistics.
#'
#' If \code{graph.them} is \code{NULL}, then \code{local_scan} computes the
#' \sQuote{us} variant of the scan statistics.  Otherwise,
#' \code{graph.them} should be an igraph object and the \sQuote{them}
#' variant is computed using \code{graph.us} to extract the neighborhood
#' information, and applying \code{FUN} on these neighborhoods in
#' \code{graph.them}.
#'
#' @param graph.us,graph An igraph object, the graph for which the scan
#'   statistics will be computed
#' @param graph.them An igraph object or \code{NULL}, if not \code{NULL},
#'   then the \sQuote{them} statistics is computed, i.e. the neighborhoods
#'   calculated from \code{graph.us} are evaluated on \code{graph.them}.
#' @param k An integer scalar, the size of the local neighborhood for each
#'   vertex. Should be non-negative.
#' @param FUN Character, a function name, or a function object itself, for
#'   computing the local statistic in each neighborhood. If \code{NULL}(the
#'   default value), \code{ecount} is used for unweighted graphs (if
#'   \code{weighted=FALSE}) and a function that computes the sum of edge
#'   weights is used for weighted graphs (if \code{weighted=TRUE}). This
#'   argument is ignored if \code{k} is zero.
#' @param weighted Logical scalar, TRUE if the edge weights should be used
#'   for computation of the scan statistic. If TRUE, the graph should be
#'   weighted.  Note that this argument is ignored if \code{FUN} is not
#'   \code{NULL}, \code{"ecount"} and \code{"sumweights"}.
#' @param mode Character scalar, the kind of neighborhoods to use for the
#'   calculation. One of \sQuote{\code{out}}, \sQuote{\code{in}},
#'   \sQuote{\code{all}} or \sQuote{\code{total}}. This argument is ignored
#'   for undirected graphs.
#' @param neighborhoods A list of neighborhoods, one for each vertex, or
#'   \code{NULL}. If it is not \code{NULL}, then the function is evaluated on
#'   the induced subgraphs specified by these neighborhoods.
#'
#'   In theory this could be useful if the same \code{graph.us} graph is used
#'   for multiple \code{graph.them} arguments. Then the neighborhoods can be
#'   calculated on \code{graph.us} and used with multiple graphs. In
#'   practice, this is currently slower than simply using \code{graph.them}
#'   multiple times.
#' @param \dots Arguments passed to \code{FUN}, the function that computes
#'   the local statistics.
#' @return For \code{local_scan} typically a numeric vector containing the
#'   computed local statistics for each vertex. In general a list or vector
#'   of objects, as returned by \code{FUN}.
#'
#' @references Priebe, C. E., Conroy, J. M., Marchette, D. J., Park,
#'   Y. (2005).  Scan Statistics on Enron Graphs. \emph{Computational and
#'   Mathematical Organization Theory}.
#'
#' @family scan statistics
#' @export
#' @examples
#' pair <- sample_correlated_gnp_pair(n = 10^3, corr = 0.8, p = 0.1)
#' local_0_us <- local_scan(graph.us = pair$graph1, k = 0)
#' local_1_us <- local_scan(graph.us = pair$graph1, k = 1)
#'
#' local_0_them <- local_scan(graph.us = pair$graph1,
#'                            graph.them = pair$graph2, k = 0)
#' local_1_them <- local_scan(graph.us = pair$graph1,
#'                            graph.them = pair$graph2, k = 1)
#'
#' Neigh_1 <- neighborhood(pair$graph1, order = 1)
#' local_1_them_nhood <- local_scan(graph.us = pair$graph1,
#'                                  graph.them = pair$graph2,
#'                                  neighborhoods = Neigh_1)

local_scan <- function(graph.us, graph.them=NULL, k=1, FUN=NULL,
                       weighted=FALSE, mode=c("out", "in", "all"),
                       neighborhoods=NULL, ...) {

  ## Must be igraph object
  stopifnot(is.igraph(graph.us))

  ## Must be NULL or igraph object
  stopifnot(is.null(graph.them) || is.igraph(graph.them))

  ## If given, number of vertices must match
  stopifnot(is.null(graph.them) || vcount(graph.them) == vcount(graph.us))

  ## k must be non-negative integer
  stopifnot(length(k)==1, k >= 0, as.integer(k) == k)

  ## Must be NULL or a function
  stopifnot(is.null(FUN) || is.function(FUN) ||
            (is.character(FUN) && length(FUN) == 1))

  ## Logical scalar
  stopifnot(is.logical(weighted), length(weighted )== 1)

  ## If weighted, then the graph(s) must be weighted
  stopifnot(!weighted || (is.weighted(graph.us) && (is.null(graph.them) ||
                                                    is.weighted(graph.them))))

  ## Check if 'neighborhoods' makes sense
  if (!is.null(neighborhoods)) {
    stopifnot(is.list(neighborhoods))
    stopifnot(length(neighborhoods) == vcount(graph.us))
  }
  if (!is.null(neighborhoods) && k==0) {
    warning("`neighborhoods' ignored for k=0")
    neighborhoods <- NULL
  }

  ## Check mode argument
  mode <- igraph.match.arg(mode)
  cmode <- switch(mode, out = 1, `in` = 2, all = 3, total = 3)

  sumweights <- function(g) sum(E(g)$weight)

  if (is.null(FUN)) { FUN <- if (weighted) "sumweights" else "ecount" }

  res <- if (is.null(graph.them)) {

    if (!is.null(neighborhoods)) {
      if (is.character(FUN) && FUN %in% c("ecount", "sumweights")) {
        neighborhoods <- lapply(neighborhoods, function(x) {
          as.integer(x)-1L
        })
        on.exit(.Call("R_igraph_finalizer", PACKAGE = "igraph"))
        .Call("R_igraph_local_scan_neighborhood_ecount", graph.us,
              if (weighted) as.numeric(E(graph.us)$weight) else NULL,
              neighborhoods, PACKAGE="igraph")
      } else {
        sapply(lapply(neighborhoods, induced.subgraph, graph=graph.us),
               FUN, ...)
      }
    } else {
      ## scan-0
      if (k == 0) {
        on.exit(.Call("R_igraph_finalizer", PACKAGE = "igraph"))
        .Call("R_igraph_local_scan_0", graph.us,
              if (weighted) as.numeric(E(graph.us)$weight) else NULL, cmode,
              PACKAGE="igraph")

        ## scan-1, ecount
      } else if (k==1 && is.character(FUN) &&
                 FUN %in% c("ecount", "sumweights")) {
        on.exit(.Call("R_igraph_finalizer", PACKAGE = "igraph"))
        .Call("R_igraph_local_scan_1_ecount", graph.us,
              if (weighted) as.numeric(E(graph.us)$weight) else NULL, cmode,
              PACKAGE="igraph")

        ## scan-k, ecount
      } else if (is.character(FUN) && FUN %in% c("ecount", "sumweights")) {
        on.exit(.Call("R_igraph_finalizer", PACKAGE = "igraph"))
        .Call("R_igraph_local_scan_k_ecount", graph.us, as.integer(k),
              if (weighted) as.numeric(E(graph.us)$weight) else NULL, cmode,
              PACKAGE="igraph")

        ## General
      } else {
        sapply(graph.neighborhood(graph.us, order=k, V(graph.us), mode=mode),
               FUN, ...)
      }
    }
    
  } else {

    if (!is.null(neighborhoods)) {
      neighborhoods <- lapply(neighborhoods, as.vector)
      if (is.character(FUN) && FUN %in% c("ecount", "wumweights")) {
        neighborhoods <- lapply(neighborhoods, function(x) {
          as.integer(x)-1L
        })
        on.exit(.Call("R_igraph_finalizer", PACKAGE = "igraph"))
        .Call("R_igraph_local_scan_neighborhood_ecount", graph.them,
              if (weighted) as.numeric(E(graph.them)$weight) else NULL,
              neighborhoods, PACKAGE="igraph")
      } else {
        sapply(lapply(neighborhoods, induced.subgraph, graph=graph.them),
               FUN, ...)
      }
    } else {

      ## scan-0
      if (k == 0) {
        on.exit(.Call("R_igraph_finalizer", PACKAGE = "igraph"))
        .Call("R_igraph_local_scan_0_them", graph.us, graph.them,
              if (weighted) as.numeric(E(graph.them)$weight) else NULL,
              cmode, PACKAGE="igraph")

        ## scan-1, ecount
      } else if (k==1 && is.character(FUN) &&
                 FUN %in% c("ecount", "sumweights")) {
        on.exit(.Call("R_igraph_finalizer", PACKAGE = "igraph"))
        .Call("R_igraph_local_scan_1_ecount_them", graph.us, graph.them,
              if (weighted) as.numeric(E(graph.them)$weight) else NULL,
              cmode, PACKAGE="igraph")

        ## scan-k, ecount
      } else if (is.character(FUN) && FUN %in% c("ecount", "sumweights")) {
        on.exit(.Call("R_igraph_finalizer", PACKAGE = "igraph"))
        .Call("R_igraph_local_scan_k_ecount_them", graph.us, graph.them,
              as.integer(k),
              if (weighted) as.numeric(E(graph.them)$weight) else NULL,
              cmode, PACKAGE="igraph")

        ## general case
      } else {
        sapply(V(graph.us), function(x) {
          vei <- neighborhood(graph.us, order=k, nodes=x, mode=mode)[[1]]
          if (!is.function(FUN)) { FUN <- getFunction(FUN, where=environment()) }
          FUN(induced.subgraph(graph.them, vei), ...)
        })
      }
    }
  }

  res <- as.numeric(res)

  if (igraph_opt("add.vertex.names") && is_named(graph.us)) {
    names(res) <- V(graph.us)$name
  }

  res
}


#' Scan statistics on a time series of graphs
#'
#' Calculate scan statistics on a time series of graphs.
#' This is done by calculating the local scan statistics for
#' each graph and each vertex, and then normalizing across the
#' vertices and across the time steps.
#'
#' @param graphs A list of igraph graph objects. They must be all directed
#'   or all undirected and they must have the same number of vertices.
#' @param tau The number of previous time steps to consider for the
#'   time-dependent normalization for individual vertices.  In other words,
#'   the current locality statistics of each vertex will be compared to this
#'   many previous time steps of the same vertex to decide whether it is
#'   significantly larger.
#' @param ell The number of previous time steps to consider
#'   for the aggregated scan statistics. This is essentially a smoothing
#'   parameter.
#' @param locality Whether to calculate the \sWuote{us} or \sQuote{them}
#'   statistics.
#' @param ... Extra arguments are passed to \code{\link{local_scan}}.
#' @return A list with entries:
#'   \item{stat}
#'
#' @family scan statistics
#' @export
#' @examples
#' ## Generate a bunch of SBMs, with the last one being different
#' num_t <- 20
#' block_sizes <- c(10, 5, 5)
#' p_ij <- list(p = 0.1, h = 0.9, q = 0.9)
#'
#' P0 <- matrix(p_ij$p, 3, 3)
#' P0[2, 2] <- p_ij$h
#' PA <- P0
#' PA[3, 3] <- p_ij$q
#' num_v <- sum(block_sizes)
#'
#' tsg <- replicate(num_t - 1, P0, simplify = FALSE) %>%
#'   append(list(PA)) %>%
#'   lapply(sample_sbm, n = num_v, block.sizes = block_sizes, directed = TRUE)
#'
#' scan_stat(graphs = tsg, k = 1, tau = 4, ell = 2)
#' scan_stat(graphs = tsg, locality = "them", k = 1, tau = 4, ell = 2)

scan_stat <- function(graphs, tau = 1, ell = 0,
                      locality = c("us", "them"), ...) {

  ## List of igraph graphs, all have same directedness and
  ## weightedness
  stopifnot(is.list(graphs),
            length(graphs) > 0,
            all(sapply(graphs, is_igraph)),
            length(unique(sapply(graphs, is_directed))) == 1,
            length(unique(sapply(graphs, gorder))) == 1)

  ## tau must the a non-negative integer
  stopifnot(length(tau) == 1, tau >= 0, as.integer(tau) == tau)

  ## ell must the a non-negative integer
  stopifnot(length(ell) == 1, ell >= 0, as.integer(ell) == ell)

  locality <- igraph.match.arg(locality)

  ## number of time steps and number of vertices
  maxTime = length(graphs)
  nVertex = vcount(graphs[[1]])

  if (locality == 'us') {
    ## Underlying locality stat is us
    lstatPsi <- matrix(0, nrow = nVertex , ncol = maxTime)
    for (i in 1:maxTime) {
      ## locality statistics \Psi over all vertices at t=i
      lstatPsi[,i] <- local_scan(graphs[[i]], ...)
    }
    lstat <- lstatPsi

  } else if (locality == 'them') {
    ## Underlying locality stat is \Phi

    lstatPhi <- array(0, dim = c(nVertex, (tau + 1), maxTime))
    for (i in 1:maxTime) {
      if (i > tau) {
        ## graph to trace k-th order neighborhood
        g <- graphs[[i]]
        for (j in 0:tau) {
          ## locality statistics \Phi over all vertices with t=i and t'=i-tau+j
          lstatPhi[, (j + 1), i] <- local_scan(
            graph.us = graphs[[i]],
            graph.them= graphs[[i - tau + j]],
            ...
          )
        }
      }
    }
    lstat <- lstatPhi
  }

  ## vertex-dependent and temporal normalization
  scan_temp_norm(
    scan_vertex_norm(lstat, tau),
    tau,
    ell
  )
}

scan_vertex_norm <-function (input_stat, tau) {
  if (is.matrix(input_stat)) {
    n <- nrow(input_stat)
    nbins <- ncol(input_stat)
    nstat <- matrix(0, n, nbins)
    for (i in 1:nbins) {
      if (i > tau) {
        if (tau == 0) {
          nstat[,i] <- input_stat[, i]
        } else {
          muv <- apply(as.matrix(input_stat[, (i - tau):(i-1)]), 1, mean)
          sdv <- apply(as.matrix(input_stat[, (i - tau):(i-1)]), 1, sd)
          sdv[is.na(sdv)] <- 1
          nstat[, i] <- (input_stat[, i] - muv) / pmax(sdv, 1)
        }
      }
    }
  } else {
    dd <- dim(input_stat)
    n <- dd[1]
    nbins <- dd[3]
    nstat <- matrix(0, n, nbins)
    for (i in 1:nbins) {
      if (i > tau) {
        if (tau == 0) {
          nstat[, i] <- input_stat[, (tau + 1), i]
        } else {
          muv <- apply(as.matrix(input_stat[, (1 : tau), i]), 1, mean)
          sdv <- apply(as.matrix(input_stat[, (1 : tau), i]), 1, sd)
          sdv[is.na(sdv)] <- 1
          nstat[, i] <- (input_stat[, (tau + 1),i] - muv) / pmax(sdv, 1)
        }
      }
    }
  }
  return(nstat)
}

scan_temp_norm <- function (stat, tau, ell) {
  maxTime <- ncol(stat)
  Mtilde <- apply(stat, 2, max)
  argmaxV <- apply(stat, 2, which.max)

  if (ell == 0) {
    list(stat = Mtilde, arg_max_v = argmaxV)

  } else if(ell ==1 ) {
    list(stat = Mtilde - c(0, Mtilde[-maxTime]), arg_max_v = argmaxV)

  } else {
    muMtilde <- rep(0, maxTime)
    sdMtilde <- rep(1, maxTime)
    for (i in (ell + 1):maxTime) {
      muMtilde[i] <- mean(Mtilde[(i - ell):(i - 1)])
      sdMtilde[i] <- sd(Mtilde[(i - ell):(i - 1)])
    }
    sstat <- (Mtilde - muMtilde) / pmax(sdMtilde, 1)
    sstat[1:(tau + ell)] <- 0
    argmaxV[1:(tau + ell)] <- NA
    list(stat = sstat, arg_max_v = argmaxV)
  }
}
