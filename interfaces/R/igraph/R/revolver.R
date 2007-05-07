
#   IGraph R package
#   Copyright (C) 2007  Gabor Csardi <csardi@rmki.kfki.hu>
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

revolver.d <- function(graph, niter=5, sd=FALSE, norm=FALSE,
                      cites=FALSE, expected=FALSE, error=TRUE, debug=numeric(),
                      verbose=igraph.par("verbose")) {

  if (!is.igraph(graph)) {
    stop("Not a graph object!")
  }

  .Call("R_igraph_revolver_d", graph, as.numeric(niter), as.logical(sd),
        as.logical(norm), as.logical(cites), as.logical(expected),
        as.logical(error), as.numeric(debug), as.logical(verbose),
        PACKAGE="igraph")
}

revolver.ad <- function(graph, niter=5, agebins=max(vcount(graph)/7100, 10),
                       sd=FALSE, norm=FALSE, cites=FALSE, expected=FALSE, error=TRUE,
                       debug=matrix(nc=2, nr=0), verbose=igraph.par("verbose")) {
  if (!is.igraph(graph)) {
    stop("Not a graph object!")
  }

  .Call("R_igraph_revolver_ad", graph, as.numeric(niter), as.numeric(agebins),
        as.logical(sd), as.logical(norm), as.logical(cites),
        as.logical(expected), as.logical(error),
        structure(as.numeric(debug), dim=dim(debug)), as.logical(verbose),
        PACKAGE="igraph")
}

revolver.ade <- function(graph, cats, niter=5, agebins=max(vcount(graph)/7100, 10),
                        sd=FALSE, norm=FALSE, cites=FALSE, expected=FALSE,
                        error=TRUE, debug=matrix(nc=2, nr=0),
                        verbose=igraph.par("verbose")) {
  if (!is.igraph(graph)) {
    stop("Not a graph object!")
  }

  .Call("R_igraph_revolver_ade", graph, as.numeric(cats), as.numeric(niter),
        as.numeric(agebins), as.logical(sd), as.logical(norm), as.logical(cites),
        as.logical(expected), as.logical(error),
        structure(as.numeric(debug), dim=dim(debug)), as.logical(verbose),
        PACKAGE="igraph")
}

revolver.e <- function(graph, cats, niter=5,
                      sd=FALSE, norm=FALSE, cites=FALSE, expected=FALSE,
                      error=TRUE, debug=numeric(), verbose=igraph.par("verbose")) {
  if (!is.igraph(graph)) {
    stop("Not a graph object!")
  }

  .Call("R_igraph_revolver_e", graph, as.numeric(cats), as.numeric(niter),
        as.logical(sd), as.logical(norm), as.logical(cites), as.logical(expected),
        as.logical(error), as.numeric(debug), as.logical(verbose),
        PACKAGE="igraph")
}
                                                     
revolver.de <- function(graph, cats, niter=5,
                      sd=FALSE, norm=FALSE, cites=FALSE, expected=FALSE,
                      error=TRUE, debug=numeric(), verbose=igraph.par("verbose")) {
  if (!is.igraph(graph)) {
    stop("Not a graph object!")
  }

  .Call("R_igraph_revolver_de", graph, as.numeric(cats), as.numeric(niter),
        as.logical(sd), as.logical(norm), as.logical(cites), as.logical(expected),
        as.logical(error), as.numeric(debug), as.logical(verbose),
        PACKAGE="igraph")
}
                       
revolver.l <- function(graph, niter=5, agebins=max(vcount(graph)/7100, 10),
                      sd=FALSE, norm=FALSE, cites=FALSE, expected=FALSE,
                      error=TRUE, debug=numeric(), verbose=igraph.par("verbose")) {
  if (!is.igraph(graph)) {
    stop("Not a graph object!")
  }

  .Call("R_igraph_revolver_l", graph, as.numeric(niter), as.numeric(agebins),
        as.logical(sd), as.logical(norm), as.logical(cites), as.logical(expected),
        as.logical(error), as.numeric(debug), as.logical(verbose),
        PACKAGE="igraph")
}

revolver.dl <- function(graph, niter=5, agebins=max(vcount(graph)/7100, 10),
                       sd=FALSE, norm=FALSE, cites=FALSE, expected=FALSE,
                       error=TRUE, debug=numeric(), verbose=igraph.par("verbose")) {
  if (!is.igraph(graph)) {
    stop("Not a graph object!")
  }

  .Call("R_igraph_revolver_dl", graph, as.numeric(niter), as.numeric(agebins),
        as.logical(sd), as.logical(norm), as.logical(cites), as.logical(expected),
        as.logical(error), as.numeric(debug), as.logical(verbose),
        PACKAGE="igraph")
}

revolver.el <- function(graph, cats, niter=5, agebins=max(vcount(graph)/7100, 10),
                       sd=FALSE, norm=FALSE, cites=FALSE, expected=FALSE,
                       error=TRUE, debug=numeric(), verbose=igraph.par("verbose")) {
  if (!is.igraph(graph)) {
    stop("Not a graph object!")
  }

  .Call("R_igraph_revolver_el", graph, as.numeric(cats), as.numeric(niter),
        as.numeric(agebins),
        as.logical(sd), as.logical(norm), as.logical(cites), as.logical(expected),
        as.logical(error), as.numeric(debug), as.logical(verbose),
        PACKAGE="igraph")
}

revolver.r <- function(graph, window, niter=5, sd=FALSE, norm=FALSE,
                      cites=FALSE, expected=FALSE, error=TRUE, debug=numeric(),
                      verbose=igraph.par("verbose")) {

  if (!is.igraph(graph)) {
    stop("Not a graph object!")
  }

  .Call("R_igraph_revolver_r", graph, as.numeric(niter), as.numeric(window),
        as.logical(sd), as.logical(norm), as.logical(cites), as.logical(expected),
        as.logical(error), as.numeric(debug), as.logical(verbose),
        PACKAGE="igraph")
}

revolver.ar <- function(graph, window, niter=5, agebins=max(vcount(graph)/7100, 10),
                       sd=FALSE, norm=FALSE, cites=FALSE, expected=FALSE, error=TRUE,
                       debug=matrix(nc=2, nr=0), verbose=igraph.par("verbose")) {
  if (!is.igraph(graph)) {
    stop("Not a graph object!")
  }

  .Call("R_igraph_revolver_ar", graph, as.numeric(niter), as.numeric(agebins),
        as.numeric(window),
        as.logical(sd), as.logical(norm), as.logical(cites),
        as.logical(expected), as.logical(error),
        structure(as.numeric(debug), dim=dim(debug)), as.logical(verbose),
        PACKAGE="igraph")
}

revolver.di <- function(graph, cats, niter=5,
                      sd=FALSE, norm=FALSE, cites=FALSE, expected=FALSE,
                      error=TRUE, debug=numeric(), verbose=igraph.par("verbose")) {
  if (!is.igraph(graph)) {
    stop("Not a graph object!")
  }

  .Call("R_igraph_revolver_di", graph, as.numeric(cats), as.numeric(niter),
        as.logical(sd), as.logical(norm), as.logical(cites), as.logical(expected),
        as.logical(error), as.numeric(debug), as.logical(verbose),
        PACKAGE="igraph")
}

revolver.adi <- function(graph, cats, niter=5, agebins=max(vcount(graph)/7100, 10),
                        sd=FALSE, norm=FALSE, cites=FALSE, expected=FALSE,
                        error=TRUE, debug=matrix(nc=2, nr=0),
                        verbose=igraph.par("verbose")) {
  if (!is.igraph(graph)) {
    stop("Not a graph object!")
  }

  .Call("R_igraph_revolver_adi", graph, as.numeric(cats), as.numeric(niter),
        as.numeric(agebins), as.logical(sd), as.logical(norm), as.logical(cites),
        as.logical(expected), as.logical(error),
        structure(as.numeric(debug), dim=dim(debug)), as.logical(verbose),
        PACKAGE="igraph")
}

revolver.il <- function(graph, cats, niter=5, agebins=max(vcount(graph)/7100, 10),
                      sd=FALSE, norm=FALSE, cites=FALSE, expected=FALSE,
                      error=TRUE, debug=numeric(), verbose=igraph.par("verbose")) {
  if (!is.igraph(graph)) {
    stop("Not a graph object!")
  }

  .Call("R_igraph_revolver_il", graph, as.numeric(cats),
        as.numeric(niter), as.numeric(agebins),
        as.logical(sd), as.logical(norm), as.logical(cites), as.logical(expected),
        as.logical(error), as.numeric(debug), as.logical(verbose),
        PACKAGE="igraph")
}

revolver.ir <- function(graph, cats, window, niter=5,
                      sd=FALSE, norm=FALSE, cites=FALSE, expected=FALSE,
                      error=TRUE, debug=numeric(), verbose=igraph.par("verbose")) {
  if (!is.igraph(graph)) {
    stop("Not a graph object!")
  }

  .Call("R_igraph_revolver_ir", graph, as.numeric(cats), as.numeric(window),
        as.numeric(niter),
        as.logical(sd), as.logical(norm), as.logical(cites), as.logical(expected),
        as.logical(error), as.numeric(debug), as.logical(verbose),
        PACKAGE="igraph")
}

revolver.air <- function(graph, cats, window,
                        niter=5, agebins=max(vcount(graph)/7100, 10),
                        sd=FALSE, norm=FALSE, cites=FALSE, expected=FALSE,
                        error=TRUE, debug=matrix(nc=2, nr=0),
                        verbose=igraph.par("verbose")) {
  if (!is.igraph(graph)) {
    stop("Not a graph object!")
  }

  .Call("R_igraph_revolver_air", graph, as.numeric(cats), as.numeric(window),
        as.numeric(niter),
        as.numeric(agebins), as.logical(sd), as.logical(norm), as.logical(cites),
        as.logical(expected), as.logical(error),
        structure(as.numeric(debug), dim=dim(debug)), as.logical(verbose),
        PACKAGE="igraph")
}
