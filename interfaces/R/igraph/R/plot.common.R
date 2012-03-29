
#   IGraph R package
#   Copyright (C) 2003, 2004, 2005  Gabor Csardi <csardi@rmki.kfki.hu>
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
        v <- igraph.par(n)
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

i.get.layout <- function(graph, layout, layout.par) {

  if (is.function(layout)) {
    layout <- layout(graph, layout.par)
  } else if (is.character(layout) && length(layout)==1 &&
             substr(layout, 1, 2)=="a:") {
    layout <- matrix(unlist(get.vertex.attribute(graph, substring(layout,3))),
                     nrow=vcount(graph), byrow=TRUE)[,1:2]
  }  

  layout
}

i.get.vertex.color <- function(graph, vertex.color) {

  if (length(vertex.color)==1 && substr(vertex.color, 1, 2)=="a:") {
    vertex.color <- unlist(get.vertex.attribute(graph,
                                                substring(vertex.color,3)))
  }

  if (is.numeric(vertex.color)) {
    vertex.color <- vertex.color %% length(palette())
    vertex.color[vertex.color==0] <- length(palette())
    vertex.color <- palette()[vertex.color]
  }

  vertex.color  
}

i.get.vertex.frame.color <- function(graph, vertex.frame.color) {

  if (length(vertex.frame.color)==1 &&
      substr(vertex.frame.color, 1, 2)=="a:") {
    vertex.frame.color <-
      unlist(get.vertex.attribute(graph,
                                  substring(vertex.frame.color,3)))
  }

  if (is.numeric(vertex.frame.color)) {
    vertex.frame.color <- vertex.frame.color %% length(palette())
    vertex.frame.color[vertex.frame.color==0] <- length(palette())
    vertex.frame.color <- palette()[vertex.frame.color]
  }

  vertex.frame.color  
}

i.get.vertex.size <- function(graph, vertex.size) {

  if (is.character(vertex.size) &&
      length(vertex.size)==1 && substr(vertex.size, 1, 2)=="a:") {
    vertex.size <- as.numeric(get.vertex.attribute
                                 (graph, substring(vertex.size,3)))
  }
  vertex.size
}

i.get.edge.color <- function(graph, edge.color) {

  if (length(edge.color)==1 && substr(edge.color, 1, 2)=="a:") {
    edge.color <- as.character(get.edge.attribute
                               (graph, substring(edge.color,3)))
  }

  if (is.numeric(edge.color)) {
    edge.color <- edge.color %% length(palette())
    edge.color[edge.color==0] <- length(palette())
    edge.color <- palette()[edge.color]
  }
  edge.color
}

i.get.edge.width <- function(graph, edge.width) {

  if (is.character(edge.width) &&
      length(edge.width)==1 && substr(edge.width, 1, 2)=="a:") {
    edge.width <- as.character(get.edge.attribute
                               (graph, substring(edge.width,3)))
  }
  edge.width
}

i.get.edge.labels <- function(graph, edge.labels=NULL) {

  if (is.null(edge.labels)) {
    edge.labels <- rep(NA, ecount(graph))
  }

  edge.labels
}

i.get.label.degree <- function(graph, label.degree) {

  if (is.character(label.degree) &&
      length(label.degree)==1 && substr(label.degree, 1, 2)=="a:") {
    label.degree <- as.numeric(get.vertex.attribute
                               (graph, substring(label.degree,3)))
  }
  label.degree
}

i.get.labels <- function(graph, labels=NULL) {

  if (is.null(labels)) {
    labels <- 0:(vcount(graph)-1)
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
                           shape="circle"),
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
                           curved=FALSE,
                           arrow.width=1),
                         plot=list(layout=layout.random,
                           margin=c(0,0,0,0),
                           rescale=TRUE,
                           asp=1,
                           frame=FALSE))
