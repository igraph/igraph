
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

###################################################################
# Pit a power-law (khmm a Yule really) distribution,
# this is a common degree distribution in networks
###################################################################

power.law.fit <- function(x, xmin=NULL, start=2, ...) {

  if (length(x) == 0) {
    stop("zero length vector")
  }
  if (length(x) == 1) {
    stop("vector should be at least of length two")
  }  

  require(stats4)
  
  if (is.null(xmin)) { xmin <- min(x) }
  
  n <- length(x)
  x <- x[ x >= xmin]
  if (length(x) != n) {
    warning("too small values eliminated from vector")
    n <- length(x)
  }
  
#  mlogl <- function(alpha) {
#    if (xmin > 1) {
#      C <- 1/(1/(alpha-1)-sum(beta(1:(xmin-1), alpha)))
#    } else {
#      C <- alpha-1
#    }
#    -n*log(C)-sum(lbeta(x, alpha))
#  }

  mlogl <- function(alpha) {
     C <- 1/sum( (xmin:10000)^-alpha )
     -n*log(C)+alpha*sum(log(x))
  }

  alpha <- mle(mlogl, start=list(alpha=start), ...)

  alpha
}
