
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

ba.game <- function(n, m=NULL, out.dist=NULL, out.seq=NULL,
                    out.pref=FALSE, directed=TRUE) {

  # Checks
  if (! is.null(out.seq) && (!is.null(m) || !is.null(out.dist))) {
    warning("if `out.seq' is given `m' and `out.dist' should be NULL")
    m <- out.dist <- NULL
  }
  if (is.null(out.seq) && !is.null(out.dist) && !is.null(m)) {
    warning("if `out.dist' is given `m' will be ignored")
    m <- NULL
  }
  if (!is.null(out.seq) && length(out.seq) != n) {
    stop("`out.seq' should be of length `n'")
  }
  if (!is.null(out.seq) && min(out.seq)<0) {
    stop("negative elements in `out.seq'");
  }
  if (!is.null(m) && m<0) {
    stop("`m' is negative")
  }
  if (!is.null(m) && m==0) {
    warning("`m' is zero, graph will be empty")
  }

  if (is.null(m) && is.null(out.dist) && is.null(out.seq)) {
    m <- 1
  }
  
  n <- as.numeric(n)
  if (!is.null(m)) { m <- as.numeric(m) }
  if (!is.null(out.dist)) { out.dist <- as.numeric(out.dist) }
  if (!is.null(out.seq)) { out.seq <- as.numeric(out.seq) }
  out.pref <- as.logical(out.pref)

  if (!is.null(out.dist)) {
    out.seq <- as.numeric(sample(0:(length(out.dist)-1), n,
                                 replace=TRUE, prob=out.dist))
  }

  if (is.null(out.seq)) {
    out.seq <- numeric()
  }

  .Call("R_igraph_barabasi_game", n, m, out.seq, out.pref, directed,
        PACKAGE="igraph")
}

barabasi.game <- ba.game

erdos.renyi.game <- function(n, p.or.m, type="gnp",
                             directed=FALSE, loops=FALSE, ...) {
  if (is.character(type)) {
    type <- switch(type, "gnp"=0, "gnm"=1)
  }

  .Call("R_igraph_erdos_renyi_game", as.numeric(n), as.numeric(type),
        as.numeric(p.or.m), as.logical(directed), as.logical(loops),
        PACKAGE="igraph")
}

random.graph.game <- erdos.renyi.game

degree.sequence.game <- function(out.deg, in.deg=numeric(0),
                                 method="simple", ...) {

  if (is.character(method)) {
    method <- switch(method, "simple"=0)
  }

  .Call("R_igraph_degree_sequence_game", as.numeric(out.deg),
        as.numeric(in.deg), as.numeric(method),
        PACKAGE="igraph")
}

growing.random.game <- function(n, m=1, directed=TRUE, citation=FALSE) {
  .Call("R_igraph_growing_random_game", as.numeric(n), as.numeric(m),
        as.logical(directed), as.logical(citation),
        PACKAGE="igraph")
}
  
aging.prefatt.game <- function(n, m=1, aging.type="exponential",
                               params=list(), ...) {

  if (! aging.type %in% c("exponential", "powerlaw")) {
    stop("Invalid aging type.")
  }
  
  probs <- rep(1, times=n)
  ind <- rep(0, times=n)
  born <- 1:n
  edges <- numeric( (n-1)*m )
  edgep <- 1

  if (aging.type=="exponential") {
    aging.exp <- params[["aging.exp"]]
    if (is.null(aging.exp)) { aging.exp <- 1 }
  } else if (aging.type=="powerlaw") {
    aging.exp <- params[["aging.exp"]]
    if (is.null(aging.exp)) { aging.exp <- 1 }
  }
  
  for (step in 2:n) {

    # choose neighbors
    newneis <- sample(1:(step-1), m, replace=TRUE,
                      prob=probs[1:(step-1)])
    
    # aging in probs
    if (aging.type=="exponential") {
      probs[1:(step-1)] <- probs[1:(step-1)] * exp(-aging.exp)
    } else if (aging.type=="powerlaw") {
      probs[1:(step-1)] <- probs[1:(step-1)] *
        ((step:2)/((step-1):1))^-aging.exp
    }
      
    # add the edges, recalculate probs
    for (nei in newneis) {
      if (aging.type=="exponential") {
        probs[nei] <- probs[nei] * ((ind[nei]+1)+1)/(ind[nei]+1)
      } else if (aging.type=="powerlaw") {
        probs[nei] <- probs[nei] * ((ind[nei]+1)+1)/(ind[nei]+1)
      }
      ind[nei] <- ind[nei]+1
      edges[edgep] <- step; edgep <- edgep + 1
      edges[edgep] <- nei ; edgep <- edgep + 1
    }
    
  }

  graph(edges-1, n=n, ...)
}

callaway.traits.game <- function(nodes, types, edge.per.step=1,
                                type.dist=rep(1, types),
                                pref.matrix=matrix(1, types, types),
                                directed=FALSE) {

  .Call("R_igraph_callaway_traits_game", as.double(nodes),
        as.double(types), as.double(edge.per.step),
        as.double(type.dist), matrix(as.double(pref.matrix), types, types),
        as.logical(directed))
}

establishment.game <- function(nodes, types, k=1, type.dist=rep(1, types),
                               pref.matrix=matrix(1, types, types),
                               directed=FALSE) {
  .Call("R_igraph_establishment_game", as.double(nodes),
        as.double(types), as.double(k), as.double(type.dist),
        matrix(as.double(pref.matrix), types, types),
        as.logical(directed))
}
        

                               
