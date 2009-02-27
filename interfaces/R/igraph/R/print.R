
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
# Convert graphs to human readable forms
###################################################################

print.igraph <- function(x,
                         graph.attributes=igraph.par("print.graph.attributes"),
                         vertex.attributes=igraph.par("print.vertex.attributes"),
                         edge.attributes=igraph.par("print.edge.attributes"),
                         names=TRUE,
                         quote.names=TRUE,
                         ...) {
  
  if (!is.igraph(x)) {
    stop("Not a graph object")
  }
  ec <- ecount(x)
  vc <- vcount(x)
  
  # From summary.graph
  cat("Vertices:", vc, "\n")
  cat("Edges:", ec, "\n")
  cat("Directed:", is.directed(x), "\n")

  # Graph attributes
  if (graph.attributes) {
    cat("Graph attributes:\n")
    list <- list.graph.attributes(x)
    lapply(list, function(n) {
      cat(sep="", "[[", n, "]]\n")
      print(get.graph.attribute(x, n))
    })
  }
  
  # Vertex attributes
  if (vertex.attributes) {
    list <- list.vertex.attributes(x)
    if (length(list) == 0) {
      cat("No vertex attributes\n")
    } else {
      cat("Vertex attributes:\n")
      mp <- getOption("max.print")
      options(max.print=1000000000)      # no built-in limit, we handle it by hand
      if (vc <= mp) {
        omitted.vertices <- 0
        ind <- as.numeric(V(x))
      } else {
        omitted.vertices <- vc-mp
        ind <- seq(length=mp)-1
      }
      if (vc==0 ||
          all(sapply(list, function(v) is.numeric(get.vertex.attribute(x, v)) |
                     is.character(get.vertex.attribute(x, v)) |
                     is.logical(get.vertex.attribute(x, v))))) {
        ## create a table
        tab <- data.frame(v=paste(sep="", "[", ind, "]"), row.names="v")
        for (i in list) {
          tab[i] <- get.vertex.attribute(x, i, ind)
        }
        print(tab)
      } else {
        for (i in ind) {
          cat(sep="", "[[", i, "]]\n")
          lapply(list, function(n) {
            cat(sep="", "[[", i, "]][[", n, "]]\n")
            print(get.vertex.attribute(x, n, i))})
        }
      }
      options(max.print=mp)
      if (omitted.vertices != 0) {
        cat(paste('[ reached getOption("max.print") -- omitted',
                  omitted.vertices, "vertices ]\n\n"))
      }      
    }
  }

  if (edge.attributes) {
    list <- list.edge.attributes(x)
  } else {
    list <- character()
  }
  
  arrow <- ifelse(is.directed(x), "->", "--")
  if (ec != 0) {
    if (!edge.attributes) {
      cat("Edges:\n")
    } else {
      cat("Edges and their attributes:\n")
    }
    if (names && ! "name" %in% list.vertex.attributes(x)) {
      names <- FALSE
    }
    if (names && "name" %in% list.vertex.attributes(x) &&
        !is.numeric(get.vertex.attribute(x, "name")) &&
        !is.character(get.vertex.attribute(x, "name")) &&
        !is.logical(get.vertex.attribute(x, "name"))) {
      warning("Can't print vertex names, complex `name' vertex attribute")
      names <- FALSE
    }
    el <- get.edgelist(x, names=names)
    if (names && quote.names) { el[] <- paste(sep="", "'", el, "'") }
    mp <- getOption("max.print")
    if (mp >= nrow(el)) {
      omitted.edges <- 0
    } else {
      omitted.edges <- nrow(el)-mp
      el <- el[1:mp,]
    }
    if (ec==0 || 
        all(sapply(list, function(v) is.numeric(get.edge.attribute(x, v)) |
                   is.character(get.edge.attribute(x,v)) |
                   is.logical(get.edge.attribute(x, v))))) {
      ## create a table
      tab <- data.frame(e=paste(sep="", "[", seq(length=nrow(el))-1, "]"), row.names="e")
      if (is.numeric(el)) { w <- nchar(max(el)) } else { w <- max(nchar(el)) }
      tab[" "] <- paste(format(el[,1], width=w), arrow, format(el[,2], width=w))
      for (i in list) {
        tab[i] <- get.edge.attribute(x, i)
      }
      print(tab)
    } else {
      i <- 0
      apply(el, 1, function(v) {
        cat(sep="", "[", i, "] ", v[1], " ", arrow, " ", v[2]);
        if (edge.attributes) {
          lapply(list, function(n) {
            cat(sep="", "\n[[", i, "]][[", n, "]]\n")
            print(get.edge.attribute(x, n, i))})
        }
        cat("\n")
        i <<- i+1
      })
    }
    if (omitted.edges != 0) {
      cat(paste('[ reached getOption("max.print") -- omitted', omitted.edges,
                'edges ]\n\n'))
    }    
  }
  
  invisible(x)
}

summary.igraph <- function(object, ...) {

  if (!is.igraph(object)) {
    stop("Not a graph object")
  }
  cat("Vertices:", vcount(object), "\n")
  cat("Edges:", ecount(object), "\n")
  cat("Directed:", is.directed(object), "\n")
  l <- list.graph.attributes(object)
  if (length(l)==0) {
    cat("No graph attributes.\n")
  } else {
    cat(sep="", "Graph attributes: ", paste(l, collapse=", "), ".\n")
  }
  l <- list.vertex.attributes(object)
  if (length(l)==0) {
    cat("No vertex attributes.\n")
  } else {
    cat(sep="", "Vertex attributes: ", paste(l, collapse=", "), ".\n")
  }
  l <- list.edge.attributes(object)
  if (length(l)==0) {
    cat("No edge attributes.\n")
  } else {
    cat(sep="", "Edge attributes: ", paste(l, collapse=", "), ".\n")
  }
  
  invisible(object)
}
