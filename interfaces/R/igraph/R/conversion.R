
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

get.adjacency <- function(graph, type=c("both", "upper", "lower"),
                          attr=NULL, names=TRUE,
                          binary=FALSE) {
  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }
  type <- igraph.match.arg(type)
  type <- switch(type, "upper"=0, "lower"=1, "both"=2)
  
  if (is.null(attr)) {    
    on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
    res <- .Call("R_igraph_get_adjacency", graph, as.numeric(type),
                 PACKAGE="igraph")
    if (binary) {
      res <- ifelse(res >= 1, 1, 0)
    }
  } else {
    attr <- as.character(attr)
    if (! attr %in% list.edge.attributes(graph)) {
      stop("no such edge attribute")
    }
    res <- matrix(0, nr=vcount(graph), nc=vcount(graph))
    if (is.directed(graph)) {
      for (i in seq(length=ecount(graph))-1) {
        e <- get.edge(graph, i)
        res[ e[1]+1, e[2]+1 ] <- get.edge.attribute(graph, attr, i)
      }
    } else {
      if (type==0) {
        ## upper
        for (i in seq(length=ecount(graph))-1) {
          e <- get.edge(graph, i)
          res[ min(e)+1, max(e)+1 ] <- get.edge.attribute(graph, attr, i)
        }        
      } else if (type==1) {
        ## lower
        for (i in seq(length=ecount(graph))-1) {
          e <- get.edge(graph, i)
          res[ max(e)+1, min(e)+1 ] <- get.edge.attribute(graph, attr, i)
        }        
      } else if (type==2) {
        ## both
        for (i in seq(length=ecount(graph))-1) {
          e <- get.edge(graph, i)
          res[ e[1]+1, e[2]+1 ] <- get.edge.attribute(graph, attr, i)
          if (e[1] != e[2]) {
            res[ e[2]+1, e[1]+1 ] <- get.edge.attribute(graph, attr, i)
          }
        }
      }
    }
  }

  if (names && "name" %in% list.vertex.attributes(graph)) {
    colnames(res) <- rownames(res) <- V(graph)$name
  }
  
  res
}

get.edgelist <- function(graph, names=TRUE) {
  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  res <- matrix(.Call("R_igraph_get_edgelist", graph, TRUE,
                      PACKAGE="igraph"), nc=2)
  if (names && "name" %in% list.vertex.attributes(graph)) {
    res <- matrix(V(graph)$name[ res+1 ], nc=2)
  }

  res
}

as.directed <- function(graph, mode=c("mutual", "arbitrary")) {
  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }

  mode <- igraph.match.arg(mode)
  mode <- switch(mode, "arbitrary"=0, "mutual"=1)
  
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  .Call("R_igraph_to_directed", graph, as.numeric(mode),
        PACKAGE="igraph")
}

as.undirected <- function(graph, mode=c("collapse", "each")) {
  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }

  mode <- igraph.match.arg(mode)
  mode <- switch(mode, "each"=0, "collapse"=1)
  
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  .Call("R_igraph_to_undirected", graph, as.numeric(mode),
        PACKAGE="igraph")  
}

get.adjlist <- function(graph, mode=c("all", "out", "in", "total")) {
  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }

  mode <- igraph.match.arg(mode)
  mode <- as.numeric(switch(mode, "out"=1, "in"=2, "all"=3, "total"=3))
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  .Call("R_igraph_get_adjlist", graph, mode,
        PACKAGE="igraph")
}

get.adjedgelist <- function(graph, mode=c("all", "out", "in", "total")) {
  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }

  mode <- igraph.match.arg(mode)
  mode <- as.numeric(switch(mode, "out"=1, "in"=2, "all"=3, "total"=3))
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  .Call("R_igraph_get_adjedgelist", graph, mode,
        PACKAGE="igraph")
}

igraph.from.graphNEL <- function(graphNEL, name=TRUE, weight=TRUE,
                                 unlist.attrs=TRUE) {

  require(graph)

  if (!is(graphNEL, "graphNEL")) {
    stop("Not a graphNEL graph")
  }
  
  al <- lapply(edgeL(graphNEL), "[[", "edges")
  al <- lapply(al, function(x) x-1)
  g <- graph.adjlist(al, directed= edgemode(graphNEL)=="directed")
  if (name) {
    V(g)$name <- nodes(graphNEL)
  }

  ## Graph attributes
  g.n <- names(graphNEL@graphData)
  g.n <- g.n [ g.n != "edgemode" ]
  for (n in g.n) {
    g <- set.graph.attribute(g, n, graphNEL@graphData[[n]])
  }
  
  ## Vertex attributes
  v.n <- names(nodeDataDefaults(graphNEL))
  for (n in v.n) {
    val <- unname(nodeData(graphNEL, attr=n))
    if (unlist.attrs && all(sapply(val, length)==1)) { val <- unlist(val) }
    g <- set.vertex.attribute(g, n, value=val)
  }

  ## Edge attributes
  e.n <- names(edgeDataDefaults(graphNEL))
  if (!weight) { e.n <- e.n [ e.n != "weight" ] }
  if (length(e.n) > 0) {
    el <- get.edgelist(g)
    el <- paste(sep="|", el[,1], el[,2])
    for (n in e.n) {
      val <- unname(edgeData(graphNEL, attr=n)[el])
      if (unlist.attrs && all(sapply(val, length)==1)) { val <- unlist(val) }
      g <- set.edge.attribute(g, n, value=val)
    }
  }
  
  g 
}

igraph.to.graphNEL <- function(graph) {

  if (!is.igraph(graph)) {
    stop("Not an igraph graph")
  }
  
  require(graph)

  if ("name" %in% list.vertex.attributes(graph) &&
      is.character(V(graph)$name)) {
    name <- V(graph)$name
  } else {
    name <- as.character(seq(vcount(graph))-1)    
  }

  edgemode <- if (is.directed(graph)) "directed" else "undirected"  

  if ("weight" %in% list.edge.attributes(graph) &&
      is.numeric(E(graph)$weight)) {
    al <- get.adjedgelist(graph, "out")
    for (i in seq(along=al)) {
      edges <- get.edges(graph, al[[i]])
      edges <- ifelse( edges[,2]==i-1, edges[,1], edges[,2])
      weights <- E(graph)$weight[al[[i]]+1]
      al[[i]] <- list(edges=edges+1, weights=weights)
    }
  } else {
    al <- get.adjlist(graph, "out")
    al <- lapply(al, function(x) list(edges=x+1))
  }  
  
  names(al) <- name
  res <- new("graphNEL", nodes=name, edgeL=al, edgemode=edgemode)

  ## Add graph attributes (other than 'directed')
  ## Are this "officially" supported at all?

  g.n <- list.graph.attributes(graph)
  if ("directed" %in% g.n) {
    warning("Cannot add graph attribute `directed'")
    g.n <- g.n[ g.n != "directed" ]
  }
  for (n in g.n) {
    res@graphData[[n]] <- get.graph.attribute(graph, n)
  }

  ## Add vertex attributes (other than 'name', that is already
  ## added as vertex names)
  
  v.n <- list.vertex.attributes(graph)
  v.n <- v.n[ v.n != "name" ]
  for (n in v.n) {
    nodeDataDefaults(res, attr=n) <- NA
    nodeData(res, attr=n) <- get.vertex.attribute(graph, n)
  }

  ## Add edge attributes (other than 'weight')
  
  e.n <- list.edge.attributes(graph)
  e.n <- e.n[ e.n != "weight" ]
  if (length(e.n) > 0) {
    el <- get.edgelist(graph)
    el <- paste(sep="|", el[,1], el[,2])
    for (n in e.n) {
      edgeDataDefaults(res, attr=n) <- NA
      res@edgeData@data[el] <- mapply(function(x,y) {
        xx <- c(x,y); names(xx)[length(xx)] <- n; xx },
                                      res@edgeData@data[el],
                                      get.edge.attribute(graph, n),
                                      SIMPLIFY=FALSE)
    }
  }
  
  res
}
