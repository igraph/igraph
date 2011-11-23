
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
  require(Matrix)
  if (M$type == "triplet") {
    sparseMatrix(dims=M$dim, i=M$i+1L, j=M$p+1L, x=M$x)
  } else {
    new("dgCMatrix", Dim=M$dim, Dimnames=list(NULL, NULL),
        factors=list(), i=M$i, p=M$p, x=M$x)
  }
}
