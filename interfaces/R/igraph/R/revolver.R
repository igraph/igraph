
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

evolver.d <- function(graph, niter=5, sd=FALSE, norm=FALSE,
                      cites=TRUE, expected=TRUE, error=TRUE, debug=numeric()) {

  if (!is.igraph(graph)) {
    stop("Not a graph object!")
  }

  .Call("R_igraph_evolver_d", graph, as.numeric(niter), as.logical(sd),
        as.logical(norm), as.logical(cites), as.logical(expected),
        as.logical(error), as.numeric(debug),
        PACKAGE="igraph")
}

evolver.ad <- function(graph, niter=5, agebins=max(vcount(graph)/7100, 10),
                       sd=FALSE, norm=FALSE, cites=TRUE, expected=TRUE, error=TRUE,
                       debug=matrix(nc=2, nr=0)) {
  if (!is.igraph(graph)) {
    stop("Not a graph object!")
  }

  .Call("R_igraph_evolver_ad", graph, as.numeric(niter), as.numeric(agebins),
        as.logical(sd), as.logical(norm), as.logical(cites),
        as.logical(expected), as.logical(error),
        structure(as.numeric(debug), dim=dim(debug)),
        PACKAGE="igraph")
}

evolver.ade <- function(graph, cats, niter=5, agebins=max(vcount(graph)/7100, 10),
                        sd=FALSE, norm=FALSE, cites=FALSE, expected=FALSE,
                        error=TRUE, debug=matrix(nc=2, nr=0)) {
  if (!is.igraph(graph)) {
    stop("Not a graph object!")
  }

  .Call("R_igraph_evolver_ade", graph, as.numeric(cats), as.numeric(niter),
        as.numeric(agebins), as.logical(sd), as.logical(norm), as.logical(cites),
        as.logical(expected), as.logical(error),
        structure(as.numeric(debug), dim=dim(debug)),
        PACKAGE="igraph")
}

evolver.e <- function(graph, cats, niter=5,
                      sd=FALSE, norm=FALSE, cites=FALSE, expected=FALSE,
                      error=TRUE, debug=numeric()) {
  if (!is.igraph(graph)) {
    stop("Not a graph object!")
  }

  .Call("R_igraph_evolver_e", graph, as.numeric(cats), as.numeric(niter),
        as.logical(sd), as.logical(norm), as.logical(cites), as.logical(expected),
        as.logical(error), as.numeric(debug),
        PACKAGE="igraph")
}
                                                     
evolver.de <- function(graph, cats, niter=5,
                      sd=FALSE, norm=FALSE, cites=FALSE, expected=FALSE,
                      error=TRUE, debug=numeric()) {
  if (!is.igraph(graph)) {
    stop("Not a graph object!")
  }

  .Call("R_igraph_evolver_de", graph, as.numeric(cats), as.numeric(niter),
        as.logical(sd), as.logical(norm), as.logical(cites), as.logical(expected),
        as.logical(error), as.numeric(debug),
        PACKAGE="igraph")
}
                       
evolver.l <- function(graph, niter=5, agebins=max(vcount(graph)/7100, 10),
                      sd=FALSE, norm=FALSE, cites=FALSE, expected=FALSE,
                      error=TRUE, debug=numeric()) {
  if (!is.igraph(graph)) {
    stop("Not a graph object!")
  }

  .Call("R_igraph_evolver_l", graph, as.numeric(niter), as.numeric(agebins),
        as.logical(sd), as.logical(norm), as.logical(cites), as.logical(expected),
        as.logical(error), as.numeric(debug),
        PACKAGE="igraph")
}
