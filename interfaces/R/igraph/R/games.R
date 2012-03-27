
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

ba.game <- function(n, power=1, m=NULL, out.dist=NULL, out.seq=NULL,
                    out.pref=FALSE, zero.appeal=1,
                    directed=TRUE, time.window=NULL) {

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
  if (!is.null(time.window) && time.window <= 0) {
    stop("time window size should be positive")
  }
  if (zero.appeal != 1 && power == 1) {
    warning("`zero.appeal' is set to 1 for traditional BA game")
  } else if (zero.appeal <= 0) {
    warning("`zero.appeal' is not positive")
  }
  if (!is.null(m) && m==0) {
    warning("`m' is zero, graph will be empty")
  }
  if (power < 0) {
    warning("`power' is negative")
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

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph0") )
  if (!is.null(time.window)) {
    .Call("R_igraph_recent_degree_game", n, power, time.window, m, out.seq,
          out.pref, as.numeric(zero.appeal), directed,
          PACKAGE="igraph0")
  } else if (power==1 && zero.appeal==1) {    
    .Call("R_igraph_barabasi_game", n, m, out.seq, out.pref, directed,
          PACKAGE="igraph0")
  } else {
    .Call("R_igraph_nonlinear_barabasi_game", n, power, m, out.seq,
          out.pref, as.numeric(zero.appeal), directed,
          PACKAGE="igraph0")
  }
}

barabasi.game <- ba.game

erdos.renyi.game <- function(n, p.or.m, type=c("gnp", "gnm"),
                             directed=FALSE, loops=FALSE, ...) {
  type <- igraph.match.arg(type)
  type <- switch(type, "gnp"=0, "gnm"=1)

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph0") )
  .Call("R_igraph_erdos_renyi_game", as.numeric(n), as.numeric(type),
        as.numeric(p.or.m), as.logical(directed), as.logical(loops),
        PACKAGE="igraph0")
}

random.graph.game <- erdos.renyi.game

degree.sequence.game <- function(out.deg, in.deg=numeric(0),
                                 method=c("simple", "vl"), ...) {

  method <- igraph.match.arg(method)
  method <- switch(method, "simple"=0, "vl"=1)

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph0") )
  .Call("R_igraph_degree_sequence_game", as.numeric(out.deg),
        as.numeric(in.deg), as.numeric(method),
        PACKAGE="igraph0")
}

growing.random.game <- function(n, m=1, directed=TRUE, citation=FALSE) {
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph0") )
  .Call("R_igraph_growing_random_game", as.numeric(n), as.numeric(m),
        as.logical(directed), as.logical(citation),
        PACKAGE="igraph0")
}

aging.prefatt.game <- function(n, pa.exp, aging.exp, m=NULL, aging.bin=300,
                               out.dist=NULL, out.seq=NULL,
                               out.pref=FALSE, directed=TRUE,
                               zero.deg.appeal=1, zero.age.appeal=0,
                               deg.coef=1, age.coef=1,
                               time.window=NULL) {
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
  if (!is.null(time.window) && time.window <= 0) {
    stop("time window size should be positive")
  }
  if (!is.null(m) && m==0) {
    warning("`m' is zero, graph will be empty")
  }
  if (pa.exp < 0) {
    warning("preferential attachment is negative")
  }
  if (aging.exp > 0) {
    warning("aging exponent is positive")
  }
  if (zero.deg.appeal <=0 ) {
    warning("initial attractiveness is not positive")
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

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph0") )
  if (is.null(time.window)) {
    .Call("R_igraph_barabasi_aging_game", as.numeric(n),
          as.numeric(pa.exp), as.numeric(aging.exp),
          as.numeric(aging.bin), m, out.seq,
          out.pref, as.numeric(zero.deg.appeal), as.numeric(zero.age.appeal),
          as.numeric(deg.coef), as.numeric(age.coef), directed, 
          PACKAGE="igraph0")
  } else {
    .Call("R_igraph_recent_degree_aging_game", as.numeric(n),
          as.numeric(pa.exp), as.numeric(aging.exp),
          as.numeric(aging.bin), m, out.seq, out.pref, as.numeric(zero.deg.appeal),
          directed,
          time.window,
          PACKAGE="igraph0")
  }
}

aging.barabasi.game <- aging.ba.game <- aging.prefatt.game

callaway.traits.game <- function(nodes, types, edge.per.step=1,
                                type.dist=rep(1, types),
                                pref.matrix=matrix(1, types, types),
                                directed=FALSE) {

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph0") )
  .Call("R_igraph_callaway_traits_game", as.double(nodes),
        as.double(types), as.double(edge.per.step),
        as.double(type.dist), matrix(as.double(pref.matrix), types, types),
        as.logical(directed),
        PACKAGE="igraph0")
}

establishment.game <- function(nodes, types, k=1, type.dist=rep(1, types),
                               pref.matrix=matrix(1, types, types),
                               directed=FALSE) {
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph0") )
  .Call("R_igraph_establishment_game", as.double(nodes),
        as.double(types), as.double(k), as.double(type.dist),
        matrix(as.double(pref.matrix), types, types),
        as.logical(directed),
        PACKAGE="igraph0")
}

grg.game <- function(nodes, radius, torus=FALSE, coords=FALSE) {
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph0") )
  res <- .Call("R_igraph_grg_game", as.double(nodes), as.double(radius),
               as.logical(torus), as.logical(coords),
               PACKAGE="igraph0")
  if (coords) {
    V(res[[1]])$x <- res[[2]]
    V(res[[1]])$y <- res[[3]]
  }
  res[[1]]
}

preference.game <- function(nodes, types, type.dist=rep(1, types),
                            pref.matrix=matrix(1, types, types),
                            directed=FALSE, loops=FALSE) {

  if (nrow(pref.matrix) != types || ncol(pref.matrix) != types) {
    stop("Invalid size for preference matrix")
  }
  
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph0") )
  .Call("R_igraph_preference_game", as.double(nodes), as.double(types),
        as.double(type.dist), matrix(as.double(pref.matrix), types, types),
        as.logical(directed), as.logical(loops),
        PACKAGE="igraph0")
}

asymmetric.preference.game <- function(nodes, types,
                                       type.dist.matrix=matrix(1, types,types),
                                       pref.matrix=matrix(1, types, types),
                                       loops=FALSE) {
  
  if (nrow(pref.matrix) != types || ncol(pref.matrix) != types) {
    stop("Invalid size for preference matrix")
  }
  if (nrow(type.dist.matrix) != types || ncol(type.dist.matrix) != types) {
    stop("Invalid size for type distribution matrix")
  }
  
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph0") )
  .Call("R_igraph_asymmetric_preference_game",
        as.double(nodes), as.double(types),
        matrix(as.double(type.dist.matrix), types, types),
        matrix(as.double(pref.matrix), types, types),
        as.logical(loops),
        PACKAGE="igraph0")
}

connect.neighborhood <- function(graph, order, mode=c("all", "out", "in", "total")) {
  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }
  mode <- igraph.match.arg(mode)
  mode <- switch(mode, "out"=1, "in"=2, "all"=3, "total"=3)

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph0") )
  .Call("R_igraph_connect_neighborhood", graph, as.numeric(order),
        as.numeric(mode),
        PACKAGE="igraph0")
}

rewire.edges <- function(graph, prob) {
  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph0") )
  .Call("R_igraph_rewire_edges", graph, as.numeric(prob),
        PACKAGE="igraph0")
}

watts.strogatz.game <- function(dim, size, nei, p) {
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph0") )
  .Call("R_igraph_watts_strogatz_game", as.numeric(dim), as.numeric(size),
        as.numeric(nei), as.numeric(p),
        PACKAGE="igraph0")
}
