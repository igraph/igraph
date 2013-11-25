
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

running.mean <- function(v, binwidth) {

  v <- as.numeric(v)
  binwidth <- as.numeric(binwidth)
  if (length(v) < binwidth) {
    stop("Vector too short for this binwidth.")
  }
  
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  .Call("R_igraph_running_mean", v, binwidth,
       PACKAGE="igraph");
}

igraph.sample <- function(low, high, length) {
  if (length>high-low+1) {
    stop("length too big for this interval")
  }
  
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  .Call("R_igraph_random_sample", as.numeric(low), as.numeric(high),
        as.numeric(length),
        PACKAGE="igraph")
}

igraph.match.arg <- function(arg, choices, several.ok=FALSE) {
  if (missing(choices)) {
    formal.args <- formals(sys.function(sys.parent()))
    choices <- eval(formal.args[[deparse(substitute(arg))]])
  }

  arg <- tolower(arg)
  choices <- tolower(choices)

  match.arg(arg=arg, choices=choices, several.ok=several.ok)
}

igraph.i.spMatrix <- function(M) {
  if (M$type == "triplet") {
    Matrix::sparseMatrix(dims=M$dim, i=M$i+1L, j=M$p+1L, x=M$x)
  } else {
    new("dgCMatrix", Dim=M$dim, Dimnames=list(NULL, NULL),
        factors=list(), i=M$i, p=M$p, x=M$x)
  }
}

dimSelect <- function(sv) {

  n <- length(sv)
  if (n == 1) return(1)

  profile <- rep(0, n)
  for (i in 1:(n-1)) {
    sv.1 <- sv[1:i]
    sv.2 <- sv[(i+1):n]
    mean.1 <- mean(sv.1)
    mean.2 <- mean(sv.2)
    var.1 <- if (i == 1) 0 else var(sv.1)
    var.2 <- if (i == n-1) 0 else var(sv.2)
    sd <- sqrt(((i-1) * var.1 + (n-i-1) * var.2) / (n-2))

    profile[i] <- sum(dnorm(sv.1, mean.1, sd, log=TRUE)) +
      sum(dnorm(sv.2, mean.2, sd, log=TRUE))
  }
  mean <- mean(sv)
  sd <- sd(sv)
  profile[n] <- sum(dnorm(sv, mean, sd, log=TRUE))
  p <- which.max(profile)
  p
}
