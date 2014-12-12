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

#' @export

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
