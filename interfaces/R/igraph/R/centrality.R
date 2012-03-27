
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

igraph.arpack.default <- list(bmat="I", n=0, which="XX", nev=1, tol=0.0,
                              ncv=3, ldv=0, ishift=1, maxiter=3000, nb=1,
                              mode=1, start=0, sigma=0.0, sigmai=0.0)

arpack <- function(func, extra=NULL, sym=FALSE, options=igraph.arpack.default,
                   env=parent.frame(), complex=!sym) {
  
  options.tmp <- igraph.arpack.default
  options.tmp[ names(options) ] <- options
  options <- options.tmp

  if (sym && complex) {
    complex <- FALSE
    warning("Symmetric matrix, setting `complex' to FALSE")
  }
  
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph0") )
  res <- .Call("R_igraph_arpack", func, extra, options, env, sym,
               PACKAGE="igraph0")

  if (complex) {
    rew <- arpack.unpack.complex(res$vectors, res$values, res$options$nev)
    res$vectors <- rew$vectors
    res$values <- rew$values

    res$values <- apply(res$values, 1, function(x) x[1]+x[2]*1i)
    dim(res$vectors) <- c(nrow(res$vectors)*2, ncol(res$vectors)/2)
    res$vectors <- apply(res$vectors, 2, function(x) {
      l <- length(x)/2
      x[1:l] + x[(l+1):length(x)]*1i
    })
  } else {
    if (is.matrix(res$values)) {
      if (!all(res$values[,2]==0)) {
        warning("Dropping imaginary parts of eigenvalues")
      }
      res$values <- res$values[,1]
    }
    res$vectors <- res$vectors[,1:length(res$values)]
  }
  
  res
}
