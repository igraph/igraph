
#   IGraph R package
#   Copyright (C) 2003-2012  Gabor Csardi <csardi.gabor@gmail.com>
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

###################################################################
# Common functions for plot and tkplot
###################################################################

i.parse.plot.params <- function(graph, params) {
  
  ## store the arguments
  p <- list(vertex=list(), edge=list(), plot=list())
  for (n in names(params)) {
    if (substr(n, 1, 7)=="vertex.") {
      nn <- substring(n, 8)
      p[["vertex"]][[nn]] <- params[[n]]
    } else if (substr(n, 1, 5)=="edge.") {
      nn <- substring(n, 6)
      p[["edge"]][[nn]] <- params[[n]]
    } else {
      p[["plot"]][[n]] <- params[[n]]
    }
  }

  ## check that names are present

  mis <- ! names(p[["vertex"]]) %in% names(i.default.values$vertex) 
  if (any(mis)) {
    stop("Unknown vertex parameters: ",
         paste(sep=", ", collapse=", ", names(p[["vertex"]])[mis]))
  }
  mis <- ! names(p[["edge"]]) %in% names(i.default.values$edge) 
  if (any(mis)) {
    stop("Unknown edge parameters: ",
         paste(sep=", ", collapse=", ", names(p[["edge"]])[mis]))
  }
  mis <- ! names(p[["plot"]]) %in% names(i.default.values$plot) 
  if (any(mis)) {
    stop("Unknown plot parameters: ",
         paste(sep=", ", collapse=", ", names(p[["plot"]]) [ mis ]))
  }
  
  func <- function(type, name, range=NULL, dontcall=FALSE) {
    if (! type %in% names(p)) {
      stop("Invalid plot option type")
    }
    ret <- function() {
      v <- p[[type]][[name]]
      if (is.function(v) && !dontcall) {
        v <- v(graph)
      }
      if (is.null(range)) {
        return (v)        
      } else {
        if (length(v)==1) {
          return(rep(v, length(range)))
        } else {
          return (rep(v, length=max(range)+1)[[range+1]])
        }
      }
    }
    if (name %in% names(p[[type]])) {
      ## we already have the parameter
      return(ret())
    } else {
      ## we don't have the parameter, check attributes first
      if (type=="vertex" && name %in% list.vertex.attributes(graph)) {
        p[[type]][[name]] <- get.vertex.attribute(graph, name)
        return(ret())
      } else if (type=="edge" && name %in% list.edge.attributes(graph)) {
        p[[type]][[name]] <- get.edge.attribute(graph, name)
        return(ret())
      } else if (type=="plot" && name %in% list.graph.attributes(graph)) {
        p[[type]][[name]] <- get.graph.attribute(graph, name)
        return(ret())
      } else {
        ## no attributes either, check igraph parameters
        n <- paste(sep="", type, ".", name)
        v <- getIgraphOpt(n)
        if (!is.null(v)) {
          p[[type]][[name]] <- v
          return(ret())
        }
        ## no igraph parameter either, use default value
        p[[type]][[name]] <- i.default.values[[type]][[name]]
        return(ret())
      }
    }
    
  }

  return (func)
}

i.get.edge.labels <- function(graph, edge.labels=NULL) {

  if (is.null(edge.labels)) {
    edge.labels <- rep(NA, ecount(graph))
  }

  edge.labels
}

i.get.labels <- function(graph, labels=NULL) {

  if (is.null(labels)) {
    labels <- seq_len(vcount(graph))
  }
  labels
}

i.get.arrow.mode <- function(graph, arrow.mode=NULL) {

  if (is.character(arrow.mode) &&
      length(arrow.mode)==1 && substr(arrow.mode, 1, 2)=="a:") {
    arrow.mode <- get.vertex.attribute(graph, substring(arrow.mode,3))
  }

  if (is.character(arrow.mode)) {
    tmp <- numeric(length(arrow.mode))
    tmp[ arrow.mode %in% c("<", "<-") ] <- 1
    tmp[ arrow.mode %in% c(">", "->") ] <- 2
    tmp[ arrow.mode %in% c("<>", "<->") ] <- 3
    arrow.mode <- tmp
  }

  if (is.null(arrow.mode)) {
    if (is.directed(graph)) {
      arrow.mode <- 2
    } else {
      arrow.mode <- 0
    }
  }

  arrow.mode
}

igraph.check.shapes <- function(x) {
  xx <- unique(x)
  bad.shapes <- ! xx %in% names(.igraph.shapes)
  if (any(bad.shapes)) {
    bs <- paste(xx[bad.shapes], collapse=", ")
    stop("Bad vertex shape(s): ", bs, ".")
  }
  x
}

autocurve.edges <- function(graph, start=0.5) {
  cm <- count.multiple(graph)
  el <- apply(get.edgelist(graph, names=FALSE), 1, paste, collapse=":")
  ord <- order(el)
  res <- numeric(length(ord))

  p <- 1
  while (p <= length(res)) {
    m <- cm[ord[p]]
    idx <- p:(p+m-1)
    if (m==1) {
      r <- 0
    } else {
      r <- seq(-start, start, length=m)
    }
    res[ord[idx]] <- r
    p <- p + m
  }

  res
}

i.default.values <- list(vertex=list(color="SkyBlue2",
                           size=15,
                           size2=15,
                           label=i.get.labels,
                           label.degree=-pi/4,
                           label.color="darkblue",
                           label.dist=0,
                           label.family="serif",
                           label.font=1,
                           label.cex=1,
                           frame.color="black",
                           shape="circle",
                           pie=1,
                           pie.color=list(c("white", "lightblue", "mistyrose",
                             "lightcyan", "lavender", "cornsilk")),
                           pie.border=list(c("white", "lightblue","mistyrose",
                             "lightcyan", "lavender", "cornsilk")),
                           pie.angle=45,
                           pie.density=-1,
                           pie.lty=1),
                         edge=list(color="darkgrey",
                           label=i.get.edge.labels,
                           lty=1,
                           width=1,
                           loop.angle=0,
                           loop.angle2=0,
                           label.family="serif",
                           label.font=1,
                           label.cex=1,
                           label.color="darkblue",
                           arrow.size=1,
                           arrow.mode=i.get.arrow.mode,
                           curved=autocurve.edges,
                           arrow.width=1),
                         plot=list(layout=layout.auto,
                           margin=c(0,0,0,0),
                           rescale=TRUE,
                           asp=1,
                           frame=FALSE))
