
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

.igraph.pars <- list("print.vertex.attributes"=FALSE,
                    "print.edge.attributes"=FALSE,
                    "print.graph.attributes"=FALSE,
                    "verbose"=FALSE,
                    "vertex.attr.comb"=list(name="concat", "ignore"),
                    "edge.attr.comb"=list(weight="sum", name="concat",
                      "ignore"),
                     "sparsematrices"=TRUE,
                     "nexus.url"="http://nexus.igraph.org"
                    )

igraph.pars.set.verbose <- function(verbose) {
  if (is.logical(verbose)) {
    .Call("R_igraph_set_verbose", verbose, PACKAGE="igraph")
  } else if (is.character(verbose)) {
    if (!verbose %in% c("tk", "tkconsole")) {
      stop("Unknown 'verbose' value")
    }
    if (verbose %in% c("tk", "tkconsole")) {
      if (!capabilities()[["X11"]]) { stop("X11 not available")           }
      if (!require("tcltk"))        { stop("tcltk package not available") }
    }
    .Call("R_igraph_set_verbose", verbose, PACKAGE="igraph")
  } else {
    stop("'verbose' should be a logical or character scalar")
  }
  verbose
}

igraph.pars.callbacks <- list("verbose"=igraph.pars.set.verbose)

## This is based on 'sm.options' in the 'sm' package

igraph.options <- function(...) {
  if (nargs() == 0) return(.igraph.pars)
  current <- .igraph.pars
  temp <- list(...)
  if (length(temp) == 1 && is.null(names(temp))) {
    arg <- temp[[1]]
    switch(mode(arg),
           list = temp <- arg,
           character = return(.igraph.pars[arg]),
           stop("invalid argument: ", sQuote(arg)))
  }
  if (length(temp) == 0) return(current)
  n <- names(temp)
  if (is.null(n)) stop("options must be given by name")
  env <- asNamespace("igraph")
  cb <- intersect(names(igraph.pars.callbacks), n)
  for (cn in cb) {
    temp[[cn]] <- igraph.pars.callbacks[[cn]](temp[[cn]])
  }
  current <- .igraph.pars               # callback might have updated it
  current[n] <- temp
  assign(".igraph.pars", current, envir = env)
  invisible(current)
}

getIgraphOpt <- function(x, default=NULL) {
  if (missing(default)) 
    return(igraph.options(x)[[1L]])
  if (x %in% names(igraph.options())) 
    igraph.options(x)[[1L]]
  else default
}

## This is deprecated from 0.6

igraph.par <- function(parid, parvalue=NULL) {

  .Deprecated("igraph.options", package="igraph")
  
  if (is.null(parvalue)) {
    res <- .igraph.pars[[parid]]
    res
  } else {
    .igraph.pars[[parid]] <- parvalue
    invisible(parvalue)
  }
}
