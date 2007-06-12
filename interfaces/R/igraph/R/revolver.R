
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

evolver.d <- function(nodes, kernel, outseq=NULL, outdist=NULL,
                      m=1, directed=TRUE, verbose=igraph.par("verbose")) {

  if (!is.null(outseq)) { outseq <- as.numeric(outseq) }
  if (!is.null(outdist)) { outdist <- as.numeric(outdist) }
  .Call("R_igraph_evolver_d", as.numeric(nodes), as.numeric(kernel),
        outseq, outdist, m, as.logical(directed), as.logical(verbose),
        PACKAGE="igraph")
}
  
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

revolver.error.d <- function(graph, kernel) {

  if (!is.igraph(graph)) {
    stop("Not a graph object!")
  }

  .Call("R_igraph_revolver_error2_d", graph, as.numeric(kernel),
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

revolver.error.ad <- function(graph, kernel) {

  if (!is.igraph(graph)) {
    stop("Not a graph object!")
  }

  kernel <- structure(as.numeric(kernel), dim=dim(kernel))
  
  .Call("R_igraph_revolver_error2_ad", graph, kernel,
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

revolver.error.ade <- function(graph, kernel, cats) {

  if (!is.igraph(graph)) {
    stop("Not a graph object!")
  }
  
  kernel <- structure(as.numeric(kernel), dim=dim(kernel))
  .Call("R_igraph_revolver_error2_ade", graph, kernel, as.numeric(cats),
        PACKAGE="igraph")
}

revolver.e <- function(graph, cats, niter=5, st=FALSE,
                      sd=FALSE, norm=FALSE, cites=FALSE, expected=FALSE,
                      error=TRUE, debug=numeric(), verbose=igraph.par("verbose")) {
  if (!is.igraph(graph)) {
    stop("Not a graph object!")
  }

  .Call("R_igraph_revolver_e", graph, as.numeric(cats), as.numeric(niter), as.logical(st),
        as.logical(sd), as.logical(norm), as.logical(cites), as.logical(expected),
        as.logical(error), as.numeric(debug), as.logical(verbose),
        PACKAGE="igraph")
}
                                                     
revolver.error.e <- function(graph, kernel, cats) {

  if (!is.igraph(graph)) {
    stop("Not a graph object!")
  }
  
  .Call("R_igraph_revolver_error2_e", graph, as.numeric(kernel), as.numeric(cats),
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
                       
revolver.error.de <- function(graph, kernel, cats) {

  if (!is.igraph(graph)) {
    stop("Not a graph object!")
  }

  kernel <- structure(as.numeric(kernel), dim=dim(kernel))
  .Call("R_igraph_revolver_error2_de", graph, kernel, as.numeric(cats),
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

revolver.error.l <- function(graph, kernel) {

  if (!is.igraph(graph)) {
    stop("Not a graph object!")
  }

  .Call("R_igraph_revolver_error2_l", graph, as.numeric(kernel),
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

revolver.error.dl <- function(graph, kernel) {

  if (!is.igraph(graph)) {
    stop("Not a graph object!")
  }

  kernel <- structure(as.numeric(kernel), dim=dim(kernel))
  .Call("R_igraph_revolver_error2_dl", graph, kernel,
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

revolver.error.el <- function(graph, kernel, cats) {

  if (!is.igraph(graph)) {
    stop("Not a graph object!")
  }

  kernel <- structure(as.numeric(kernel), dim=dim(kernel))
  .Call("R_igraph_revolver_error2_el", graph, kernel, as.numeric(cats),
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

revolver.error.r <- function(graph, kernel, window) {

  if (!is.igraph(graph)) {
    stop("Not a graph object!")
  }

  .Call("R_igraph_revolver_error2_r", graph, as.numeric(kernel), as.numeric(window),
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

revolver.error.ar <- function(graph, kernel, window) {

  if (!is.igraph(graph)) {
    stop("Not a graph object!")
  }

  kernel <- structure(as.numeric(kernel), dim=dim(kernel))
  .Call("R_igraph_revolver_error2_ar", graph, kernel, as.numeric(window),
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

revolver.error.di <- function(graph, kernel, cats) {

  if (!is.igraph(graph)) {
    stop("Not a graph object!")
  }

  kernel <- structure(as.numeric(kernel), dim=dim(kernel))
  .Call("R_igraph_revolver_error2_di", graph, kernel, as.numeric(cats),
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

revolver.error.adi <- function(graph, kernel, cats) {

  if (!is.igraph(graph)) {
    stop("Not a graph object!")
  }

  kernel <- structure(as.numeric(kernel), dim=dim(kernel))
  .Call("R_igraph_revolver_error2_adi", graph, kernel, as.numeric(cats),
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

revolver.error.il <- function(graph, kernel, cats) {

  if (!is.igraph(graph)) {
    stop("Not a graph object!")
  }

  kernel <- structure(as.numeric(kernel), dim=dim(kernel))
  .Call("R_igraph_revolver_error2_il", graph, kernel, as.numeric(cats),
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

revolver.error.ir <- function(graph, kernel, cats, window) {

  if (!is.igraph(graph)) {
    stop("Not a graph object!")
  }

  kernel <- structure(as.numeric(kernel), dim=dim(kernel))
  .Call("R_igraph_revolver_error2_ir", graph, kernel, as.numeric(cats),
        as.numeric(window),
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

revolver.error.air <- function(graph, kernel, cats, window) {

  if (!is.igraph(graph)) {
    stop("Not a graph object!")
  }

  kernel <- structure(as.numeric(kernel), dim=dim(kernel))
  .Call("R_igraph_revolver_error2_air", graph, kernel, as.numeric(cats),
        as.numeric(window),
        PACKAGE="igraph")
}

revolver.d.d <- function(graph, vtime=V(graph)$time, etime=E(graph)$time, niter=5,
                         sd=FALSE, norm=FALSE, cites=FALSE, expected=FALSE,
                         error=TRUE, debug=matrix(nc=2, nr=0),
                         verbose=igraph.par("verbose")) {
  if (!is.igraph(graph)) {
    stop("Not a graph object!")
  }

  if (is.null(vtime)) {
    stop("vtime missing")
  }
  if (is.null(etime)) {
    stop("etime missing")
  }

  .Call("R_igraph_revolver_d_d", graph, as.numeric(niter),
        as.numeric(vtime), as.numeric(etime), 
        as.logical(sd), as.logical(norm), as.logical(cites),
        as.logical(expected), as.logical(error),
        structure(as.numeric(debug), dim=dim(debug)), as.logical(verbose),
        PACKAGE="igraph")
}

revolver.p.p <- function(graph, events=get.graph.attribute(graph, "events"),
                         vtime=V(graph)$time, etime=E(graph)$time,
                         niter=5, sd=FALSE, norm=FALSE, cites=FALSE, expected=FALSE,
                         error=TRUE, debug=matrix(nc=2, nr=0),
                         verbose=igraph.par("verbose")) {
  if (!is.igraph(graph)) {
    stop("Not a graph object!")
  }
  if (is.null(events) || !is.list(events)) {
    stop("events missing or not a list")
  }
  if (is.null(vtime)) {
    ## TODO: calculate from events
    stop("vtime missing")
  }
  if (is.null(etime)) {
    ## TODO: calculate from events
    stop("etime missing")
  }
  
  authors <- unlist(events)
  eventsizes <- sapply(events, length)
  
  .Call("R_igraph_revolver_p_p", graph, as.numeric(niter),
        as.numeric(vtime), as.numeric(etime), as.numeric(authors),
        as.numeric(eventsizes), as.logical(sd), as.logical(norm), as.logical(cites),
        as.logical(expected), as.logical(error),
        structure(as.numeric(debug), dim=dim(debug)), as.logical(verbose),
        PACKAGE="igraph")
}
