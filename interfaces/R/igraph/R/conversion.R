
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

get.adjacency.dense <- function(graph, type=c("both", "upper", "lower"),
                                attr=NULL, names=TRUE,
                                binary=FALSE) {
  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }
  
  type <- igraph.match.arg(type)
  type <- switch(type, "upper"=0, "lower"=1, "both"=2)
  
  if (is.null(attr)) {    
    on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph0") )
    res <- .Call("R_igraph_get_adjacency", graph, as.numeric(type),
                 PACKAGE="igraph0")
    if (binary) {
      res <- ifelse(res >= 1, 1, 0)
    }
  } else {
    attr <- as.character(attr)
    if (! attr %in% list.edge.attributes(graph)) {
      stop("no such edge attribute")
    }
    res <- matrix(0, nrow=vcount(graph), ncol=vcount(graph))
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

get.adjacency.sparse <- function(graph, type=c("both", "upper", "lower"),
                                 attr=NULL, names=TRUE,
                                 binary=FALSE) {
  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }

  require(Matrix)
  
  type <- igraph.match.arg(type)

  vc <- vcount(graph)
  
  el <- get.edgelist(graph, names=FALSE)
  if (!is.null(attr)) {
    attr <- as.character(attr)
    if (!attr %in% list.edge.attributes(graph)) {
      stop("no such edge attribute")
    }
    value <- get.edge.attribute(graph, name=attr)
  } else {
    value <- rep(1, nrow(el))
  }

  if (is.directed(graph)) {
    res <- spMatrix(vc, vc, i=el[,1]+1, j=el[,2]+1, x=value)
  } else {
    if (type=="upper") {
      ## upper
      res <- spMatrix(vc, vc, i=pmin(el[,1],el[,2])+1,
                      j=pmax(el[,1],el[,2])+1, x=value)
    } else if (type=="lower") {
      ## lower
      res <- spMatrix(vc, vc, i=pmax(el[,1],el[,2])+1,
                      j=pmin(el[,1],el[,2])+1, x=value)
    } else if (type=="both") {
      ## both
      i <- pn <- pmin(el[,1],el[,2])+1
      j <- px <- pmax(el[,1],el[,2])+1
      non.loop <- pn != px
      i <- c(pn, px[non.loop])
      j <- c(px, pn[non.loop])
      x <- c(value, value[non.loop])
      res <- spMatrix(vc, vc, i=i, j=j, x=x)
    }
  }

  if (names && "name" %in% list.vertex.attributes(graph)) {
    colnames(res) <- rownames(res) <- V(graph)$name
  }

  res
}

get.adjacency <- function(graph, type=c("both", "upper", "lower"),
                          attr=NULL, names=TRUE,
                          binary=FALSE, sparse=FALSE) {
  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }

  if (!sparse) {
    get.adjacency.dense(graph, type=type, attr=attr, names=names, binary=binary)
  } else {
    get.adjacency.sparse(graph, type=type, attr=attr, names=names, binary=binary)
  }  
}

get.edgelist <- function(graph, names=TRUE) {
  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph0") )
  res <- matrix(.Call("R_igraph_get_edgelist", graph, TRUE,
                      PACKAGE="igraph0"), ncol=2)
  if (names && "name" %in% list.vertex.attributes(graph)) {
    res <- matrix(V(graph)$name[ res+1 ], ncol=2)
  }

  res
}

as.directed <- function(graph, mode=c("mutual", "arbitrary")) {
  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }

  mode <- igraph.match.arg(mode)
  mode <- switch(mode, "arbitrary"=0, "mutual"=1)
  
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph0") )
  .Call("R_igraph_to_directed", graph, as.numeric(mode),
        PACKAGE="igraph0")
}

as.undirected <- function(graph, mode=c("collapse", "each")) {
  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }

  mode <- igraph.match.arg(mode)
  mode <- switch(mode, "each"=0, "collapse"=1)
  
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph0") )
  .Call("R_igraph_to_undirected", graph, as.numeric(mode),
        PACKAGE="igraph0")  
}

get.adjlist <- function(graph, mode=c("all", "out", "in", "total")) {
  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }

  mode <- igraph.match.arg(mode)
  mode <- as.numeric(switch(mode, "out"=1, "in"=2, "all"=3, "total"=3))
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph0") )
  .Call("R_igraph_get_adjlist", graph, mode,
        PACKAGE="igraph0")
}

get.adjedgelist <- function(graph, mode=c("all", "out", "in", "total")) {
  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }

  mode <- igraph.match.arg(mode)
  mode <- as.numeric(switch(mode, "out"=1, "in"=2, "all"=3, "total"=3))
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph0") )
  .Call("R_igraph_get_adjedgelist", graph, mode,
        PACKAGE="igraph0")
}

igraph.from.graphNEL <- function(graphNEL, name=TRUE, weight=TRUE,
                                 unlist.attrs=TRUE) {

  require(graph)

  if (!is(graphNEL, "graphNEL")) {
    stop("Not a graphNEL graph")
  }
  
  al <- lapply(edgeL(graphNEL), "[[", "edges")
  if (edgemode(graphNEL)=="undirected") {
    al <- mapply(seq_along(al), al, SIMPLIFY=FALSE, FUN=function(n, l) {
      c(l, rep(n, sum(l==n)))
    })
  }
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

get.incidence.dense <- function(graph, types, names, attr) {

  if (is.null(attr)) {
    on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph0") )
    ## Function call
    res <- .Call("R_igraph_get_incidence", graph, types,
                 PACKAGE="igraph0")

    if (names && "name" %in% list.vertex.attributes(graph)) {
      rownames(res$res) <- V(graph)$name[ res$row_ids+1 ]
      colnames(res$res) <- V(graph)$name[ res$col_ids+1 ]
    } else {
      rownames(res$res) <- res$row_ids
      colnames(res$res) <- res$col_ids
    }
    res$res
    
  } else {

    attr <- as.character(attr)
    if (!attr %in% list.edge.attributes(graph)) {
      stop("no such edge attribute")
    }

    vc <- vcount(graph)
    n1 <- sum(!types)
    n2 <- vc-n1    
    res <- matrix(0, n1, n2)

    recode <- numeric(vc)
    recode[!types] <- seq_len(n1)
    recode[types]  <- seq_len(n2)
    
    for (i in seq(length=ecount(graph))-1) {
      eo <- get.edge(graph, i)
      e <- recode[eo+1]
      if (!types[eo[1]+1]) {
        res[ e[1], e[2] ] <- get.edge.attribute(graph, attr, i)
      } else{
        res[ e[2], e[1] ] <- get.edge.attribute(graph, attr, i)
      }
    }

    if (names && "name" %in% list.vertex.attributes(graph)) {
      rownames(res) <- V(graph)$name[ which(!types) ]
      colnames(res) <- V(graph)$name[ which( types) ]
    } else {
      rownames(res) <- which(!types)-1
      colnames(res) <- which(types)
    }

    res
  }
}

get.incidence.sparse <- function(graph, types, names, attr) {

  vc <- vcount(graph)
  if (length(types) != vc) {
    stop("Invalid types vector")
  }
  
  require(Matrix)
  el <- get.edgelist(graph, names=FALSE)
  if (any(types[el[,1]+1] == types[el[,2]+1])) {
    stop("Invalid types vector, not a bipartite graph")
  }

  n1 <- sum(!types)
  n2 <- vc-n1

  recode <- numeric(vc)
  recode[!types] <- seq_len(n1)
  recode[types]  <- seq_len(n2) + n1

  el[,1] <- recode[el[,1]+1]
  el[,2] <- recode[el[,2]+1]

  change <- el[,1] > n1
  el[change,] <- el[change,2:1]
  el[,2] <- el[,2]-n1

  if (!is.null(attr)) {
    attr <- as.character(attr)
    if (!attr %in% list.edge.attributes(graph)) {
      stop("no such edge attribute")
    }
    value <- get.edge.attribute(graph, name=attr)
  } else { 
    value <- rep(1, nrow(el))
  }

  res <- spMatrix(n1, n2, i=el[,1], j=el[,2], x=value)

  if (names && "name" %in% list.vertex.attributes(graph)) {
    rownames(res) <- V(graph)$name[which(!types)]
    colnames(res) <- V(graph)$name[which(types)]
  } else {
    rownames(res) <- which(!types)-1
    colnames(res) <- which(types)-1
  }
  res
}

get.incidence <- function(graph, types=NULL, attr=NULL,
                          names=TRUE, sparse=FALSE) {
  # Argument checks
  if (!is.igraph(graph)) { stop("Not a graph object") }
  if (is.null(types) && "type" %in% list.vertex.attributes(graph)) { 
    types <- V(graph)$type 
  } 
  if (!is.null(types)) { 
    types <- as.logical(types) 
  } else { 
    stop("Not a bipartite graph, supply `types' argument") 
  }
  
  names <- as.logical(names)
  sparse <- as.logical(sparse)
  
  if (sparse) {
    get.incidence.sparse(graph, types=types, names=names, attr=attr)
  } else {
    get.incidence.dense(graph, types=types, names=names, attr=attr)
  }
}

