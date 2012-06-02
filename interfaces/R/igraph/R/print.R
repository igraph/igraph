
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

###################################################################
# Convert graphs to human readable forms
###################################################################

.print.header <- function(object) {

  if (!is.igraph(object)) {
    stop("Not a graph object")
  }

  title <- paste(sep="", "IGRAPH ",
                 c("U","D")[is.directed(object)+1],
                 c("-","N")[is.named(object)+1],
                 c("-","W")[is.weighted(object)+1],
                 c("-","B")[is.bipartite(object)+1], " ",
                 vcount(object), " ", ecount(object), " -- ")
  w <- getOption("width")
  if (nchar(title) < w && "name" %in% list.graph.attributes(object)) {
    title <- substring(paste(sep="", title, object$name), 1, w-1)
  }
  cat(title, "\n", sep="")

  aa <- function(names, code, fun) {
    if (length(names)==0) {
      ""
    } else {
      type <- sapply(names, function(x) mode(fun(object, x)))
      type <- sapply(type, switch, "numeric"="n", "character"="c", "x")
      paste(sep="", names, " (", code, "/", type, ")")
    }
  }
  ga <- aa(list.graph.attributes(object), "g", get.graph.attribute)
  va <- aa(list.vertex.attributes(object), "v", get.vertex.attribute)
  ea <- aa(list.edge.attributes(object), "e", get.edge.attribute)
  atxt <- c(ga, va, ea)
  atxt <- paste(atxt[atxt!=""], collapse=", ")
  if (atxt != "") {
    atxt <- strwrap(paste(sep="", "+ attr: ", atxt), exdent=2)
    cat(atxt, sep="\n")
  }
}

.print.graph.attributes <- function(x) {
  list <- list.graph.attributes(x)
  if (length(list)!=0) {
    cat("+ graph attributes:\n")
    lapply(list, function(n) {
      cat(sep="", "[[", n, "]]\n")
      print(get.graph.attribute(x, n))
    })
  }
}

.print.vertex.attributes <- function(x) {
  vc <- vcount(x)
  list <- list.vertex.attributes(x)
  if (length(list) != 0) {
    cat("+ vertex attributes:\n")
    mp <- getOption("max.print")
    options(max.print=1000000000)
    if (vc <= mp) {
      omitted.vertices <- 0
      ind <- as.numeric(V(x))
    } else {
      omitted.vertices <- vc-mp
      ind <- seq(length=mp)
    }
    if (vc==0 ||
        all(sapply(list, function(v)
                   is.numeric(get.vertex.attribute(x, v)) ||
                   is.character(get.vertex.attribute(x, v)) ||
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

.print.edges.edgelist <- function(x, names) {
  ec <- ecount(x)
  list <- list.edge.attributes(x)
  list <- list[list!="name"]
  arrow <- ifelse(is.directed(x), "->", "--")
  if (is.named(x)) {
    cat("+ edges (vertex names) and their attributes:\n")
  } else {
    cat("+ edges and their attributes:\n")
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
  
  mp <- getOption("max.print")
  if (mp >= ec) {
    omitted.edges <- 0
    el <- get.edgelist(x, names=names)
  } else {
    omitted.edges <- ec-mp
    el <- get.edges(x, seq_len(mp))
    if (names) { el[] <- V(x)$name[el] }
  }
  ename <- if ("name" %in% list.edge.attributes(x)) {
    paste(sep="", "'", E(x)$name, "'")
  } else {
    seq(length=nrow(el))
  }
  if (ec==0 || 
      all(sapply(list, function(v) is.numeric(get.edge.attribute(x, v)) |
                 is.character(get.edge.attribute(x,v)) |
                 is.logical(get.edge.attribute(x, v))))) {
    ## create a table
    tab <- data.frame(row.names=paste(sep="", "[", ename, "]"))
    if (is.numeric(el)) { w <- nchar(max(el)) } else { w <- max(nchar(el)) }
    tab["edge"] <- paste(sep="", format(el[,1], width=w),
                         arrow, format(el[,2], width=w))
    for (i in list) {
      tab[i] <- get.edge.attribute(x, i)
    }
    print(tab)
  } else {
    i <- 1
    apply(el, 1, function(v) {
      cat(sep="", "[", ename[i], "] ", v[1], " ", arrow, " ", v[2]);
      lapply(list, function(n) {
        cat(sep="", "\n[[", i, "]][[", n, "]]\n")
        print(get.edge.attribute(x, n, i))})
      cat("\n")
      i <<- i+1
    })
  }
  if (omitted.edges != 0) {
    cat(paste('[ reached getOption("max.print") -- omitted', omitted.edges,
              'edges ]\n\n'))
  }    
} 

.print.edges.compressed <- function(x, names) {
  ## TODO: getOption("max.print")
  if (is.named(x)) cat("+ edges (vertex names):\n") else cat("+ edges:\n")
  el <- get.edgelist(x, names=names)
  arrow <- c("--", "->")[is.directed(x)+1]
  edges <- paste(sep="", format(el[,1]), arrow, format(el[,2]))
  print(edges, quote=FALSE)
}

.print.edges.adjlist <- function(x) {
  ## TODO: getOption("max.print")
  cat("+ edges:\n")
  vc <- vcount(x)
  arrow <- c(" -- ", " -> ")[is.directed(x)+1]
  al <- get.adjlist(x, mode="out")
  w <- nchar(max(which(degree(x, mode="in") != 0)))
  mpl <- trunc((getOption("width")-nchar(arrow)-nchar(vc)) / (w+1))
  if (any(sapply(al, length) > mpl)) {
    ## Wrapping needed
    mw <- nchar(vcount(x))
    sm <- paste(collapse="", rep(" ", mw+4))
    alstr <- lapply(seq_along(al), function(x) {
      len <- length(al[[x]])
      fac <- rep(1:(len/mpl+1), each=mpl, length=len)
      nei <- tapply(format(al[[x]], width=mw), fac, paste, collapse=" ")
      mark <- paste(sep="", format(x, width=mw), arrow)
      mark <- c(mark, rep(sm, max(0, length(nei)-1)))
      paste(sep="", mark, nei)
    })
    cat(unlist(alstr), sep="\n")
  } else { 
    alstr <- sapply(al, function(x) {
      paste(format(x, width=w), collapse=" ")
    })
    mark <- paste(sep="", format(seq_len(vc)), arrow)
    alstr <- paste(sep="", mark, alstr)
    maxw <- max(nchar(alstr))
    sep <- "   "
    ncol <- trunc((getOption("width")-1+nchar(sep)) / (maxw+nchar(sep)))
    if (ncol > 1) {
      alstr <- format(alstr, width=maxw, justify="left")
      fac <- rep(1:(vc/ncol+1), each=ncol, length=vc)
      alstr <- tapply(alstr, fac, paste, collapse=sep)
    }
    cat(alstr, sep="\n")
  }
}

.print.edges.adjlist.named <- function(x) {
  ## TODO getOption("max.print")
  cat("+ edges (vertex names):\n")

  arrow <- c(" -- ", " -> ")[is.directed(x)+1]
  vn <- V(x)$name

  al <- get.adjlist(x, mode="out")
  alstr <- sapply(al, function(x) { paste(collapse=", ", vn[x]) })
  alstr <- paste(sep="", format(vn), arrow, alstr)
  alstr <- strwrap(alstr, exdent=max(nchar(vn))+nchar(arrow))
  cat(alstr, sep="\n")
}

str.igraph <- function(object, ...) {
  print.igraph(object, full=TRUE, ...)
}

print.igraph <- function(x, full=getIgraphOpt("print.full"),
                graph.attributes=getIgraphOpt("print.graph.attributes"),
                vertex.attributes=getIgraphOpt("print.vertex.attributes"),
                edge.attributes=getIgraphOpt("print.edge.attributes"),
                names=TRUE, ...) {
  
  if (!is.igraph(x)) {
    stop("Not a graph object")
  }

  .print.header(x)
  if (full) { 
    if (graph.attributes)  .print.graph.attributes(x)
    if (vertex.attributes) .print.vertex.attributes(x)
    if (ecount(x)==0) {
      ## Do nothing
    } else if (edge.attributes && length(list.edge.attributes(x)) !=0 ) {
      .print.edges.edgelist(x, names)
    } else if (median(degree(x, mode="out")) < 3) {
      .print.edges.compressed(x, names)
    } else if (is.named(x)) {
      .print.edges.adjlist.named(x)
    } else {
      .print.edges.adjlist(x)
    }
  }
  
  invisible(x)
}

summary.igraph <- function(object, ...) {

  if (!is.igraph(object)) {
    stop("Not a graph object")
  }

  title <- paste(sep="", "IGRAPH ",
                 c("U","D")[is.directed(object)+1],
                 c("-","N")[is.named(object)+1],
                 c("-","W")[is.weighted(object)+1],
                 c("-","B")[is.bipartite(object)+1], " ",
                 vcount(object), " ", ecount(object), " -- ")
  w <- getOption("width")
  if (nchar(title) < w && "name" %in% list.graph.attributes(object)) {
    title <- substring(paste(sep="", title, object$name), 1, w-1)
  }
  cat(title, "\n", sep="")

  aa <- function(names, code, fun) {
    if (length(names)==0) {
      ""
    } else {
      type <- sapply(names, function(x) mode(fun(object, x)))
      type <- sapply(type, switch, "numeric"="n", "character"="c", "x")
      paste(sep="", names, " (", code, "/", type, ")")
    }
  }
  ga <- aa(list.graph.attributes(object), "g", get.graph.attribute)
  va <- aa(list.vertex.attributes(object), "v", get.vertex.attribute)
  ea <- aa(list.edge.attributes(object), "e", get.edge.attribute)
  atxt <- c(ga, va, ea)
  atxt <- paste(atxt[atxt!=""], collapse=", ")
  if (atxt != "") {
    atxt <- strwrap(paste(sep="", "attr: ", atxt), exdent=2)
    cat(atxt, sep="\n")
  }

  invisible(object)
}

"
####################################################################
## Various designs for printing graphs

## Summary

IGRAPH UNW- 5 5 -- A ring
Attr: name (g/c), name (v/c), weight (e/n)

IGRAPH D-W- 100 200 -- Gnm random graph

## Printing, edge list

IGRAPH-UNW--V5-E5----------------------------------------- A ring -
+ attributes: name (g), name (v), weight (e).
+ edges:
     edge  weight              
[1]' a--b       1
[2]' b--c       2
[3]' c--d      -1
[4]' d--e     0.5
[5]' a--e       1

## Compressed edge list

IGRAPH UNW- 5 10 -- A ring 
+ attributes: name (g/c), name (v/n), weight (e/n)
+ edges:
[1]' 1--2 2--3 3--4 4--5 1--5 2--5 5--1
[8]' 1--4 4--2 1--3

## This is good if vertices are named

IGRAPH UNW- 10 18 -- Krackhardt kite
+ attributes: name (g/c), name (v/c), weight (e/n)
+ edges:
Andre    -- [1] Beverly, Carol, Diane, Fernando
Beverly  -- [1] Andre, Diane, Ed, Garth
Carol    -- [1] Andre, Diane, Fernando
Diane    -- [1] Andre, Beverly, Carol, Diane, Ed
         -- [6] Garth
Ed       -- [1] Beverly, Diane, Garth
Fernando -- [1] Andre, Carol, Diane, Garth
Garth    -- [1] Beverly, Diane, Ed, Fernando
Heather  -- [1] Fernando, Garth
Ike      -- [1] Heather, Jane
Jane     -- [1] Ike

IGRAPH UNW- 10 18 -- Krackhardt kite
+ attributes: name (g/c), name (v/c), weight (e/n)
+ edges:
Andre    -- Beverly, Carol, Diane, Fernando
Beverly  -- Andre, Diane, Ed, Garth
Carol    -- Andre, Diane, Fernando
Diane    -- Andre, Beverly, Carol, Diane, Ed, Garth
Ed       -- Beverly, Diane, Garth
Fernando -- Andre, Carol, Diane, Garth
Garth    -- Beverly, Diane, Ed, Fernando
Heather  -- Fernando, Garth
Ike      -- Heather, Jane
Jane     -- Ike

## This is the good one if vertices are not named

IGRAPH U--- 100 200 -- Gnm random graph 
+ edges:
[  1] 28 46 89 90                 [  2] 47 69 72 89
[  3] 29                          [  4] 17 20
[  5] 11 40 42 51 78 89           [  6] 27 32 70 87 93
[  7] 18 27 87                    [  8] 18 24 82
[  9] 18 20 85 94                 [ 10] 24 70 77 91
[ 11]  5 12 34 61 62              [ 12] 11 41 44 61 65 80
...

## Alternative designs, summary

IGRAPH-UNW--V5-E5,---------------------------------------- A ring -
+ attributes: name (g/c), name (v/c), weight (e/n)

IGRAPH. |V|=5, |E|=5, undirected, named, weighted.
Attributes: name (g/c), name (v/c), weight (e/n)

IGRAPH: 'A ring'
Graph attributes: |V|=5, |E|=5, undirected, name.
Vertex attributes: name.
Edge attributes: weight.

## Alternative designs, printing

IGRAPH-UNW--V5-E5----------------------------------------- A ring -
'- attributes: name (g), name (v), weight (e).
'         edge  weight              
[1] 'a' -- 'b'       1
[2] 'b' -- 'c'       2
[3] 'c' -- 'd'      -1
[4] 'd' -- 'e'     0.5
[5] 'a' -- 'e'       1

IGRAPH-UNW--V-5-E-10-------------------------------------- A ring -
|- attributes: name (g), name (v), weight (e).
|- edges:
[1] 'a'--'b'  'b'--'c'  'c'--'d'  'd'--'e'  'a'--'e'  'b'-'e'
[7] 'e'--'a'  'a'--'d'  'd'--'b'  'a'--'c'


IGRAPH-UNW--V-5-E-10-------------------------------------- A ring -
+ attributes: name (g), name (v), weight (e).
+ vertices:
|     name
| [1]    a
| [2]    b
| [3]    c
| [4]    d
| [5]    e
+ edges:
[1] 'a'--'b'  'b'--'c'  'c'--'d'  'd'--'e'  'a'--'e'  'b'-'e'
[7] 'e'--'a'  'a'--'d'  'd'--'b'  'a'--'c'




IGRAPH-UNW--V-5-E-10-------------------------------------- A ring -
+ graph attributes: name
+ vertex attributes: name
+ edge attributes: weight
+ vertices:
|   name
|1]    a
|2]    b
|3]    c
|4]    d
|5]    e
+ edges:
|1] a--b  b--c  c--d  d--e  a--e  b-e
|7] e--a  a--d  d--b  a--c



IGRAPH-UNW--V-5-E-10-------------------------------------- A ring -
+ graph attributes:  name (c)
+ vertex attributes: name (c)
+ edge attributes:   weight (n)
+ edges:
[1] a--b  b--c  c--d  d--e  a--e  b-e
[7] e--a  a--d  d--b  a--c


IGRAPH-UNW--V-5-E-10-------------------------------------- A ring -
+ attributes: name (g/c), name (v/c), weight (e/n)
+ edges:
[ 1] a--b b--c c--d d--e a--e b--e e--a a--d d--b
[10] a--c

IGRAPH-DNW--V-5-E-10-------------------------------------- A ring -
+ attributes: name (g/c), name (v/n), weight (e/n)
+ edges:
[1]' 1->2 2->3 3->4 4->5 1->5 2->5 5->1
[8]' 1->4 4->2 1->3


IGRAPH-UNW--V-5-E-20-------------------------------------- A ring -
+ attributes: name (g/c), name (v/c), weight (e/n)
+ edges:
[ 1] a-b b-c c-d d-e a-e b-e e-a a-d d-b a-c
[11] a-b b-c c-d d-e a-e b-e e-a a-d d-b a-c


IGRAPH-UNW--V-8-E-10-------------------------------------- A ring -
+ attributes: name (g/c), name (v/c), weight (e/n)
+ edges:
[a] b c e f h
[b] a c e
[c] a b d
[d] a b c h
[e] a b d
[f] a
[g]
[h] a d

IGRAPH-UNW--V-10-E-18------------------------------------- A ring -
+ attributes: name (g/c), name (v/c), weight (e/n)
+ edges:
[a] a--{b,c,e,f,h}  b--{a,c,e}  c--{a,b,d}  d--{a,b,c,h}
[e] e--{a,b,d}      f--{a}      g--{}       h--{a,d}


IGRAPH-UNW--V10-E18------------------------------Krackhardt kite--
+ attributes: name (g/c), name (v/c), weight (e/n)
+ edges:
[   Andre][1] Beverly  Carol    Diane    Fernando
[ Beverly][1] Andre    Diane    Ed       Garth
[   Carol][1] Andre    Diane    Fernando
[   Diane][1] Andre    Beverly  Carol    Diane    Ed
[   Diane][6] Garth
[      Ed][1] Beverly  Diane    Garth
[Fernando][1] Andre    Carol    Diane    Garth
[   Garth][1] Beverly  Diane    Ed       Fernando
[ Heather][1] Fernando Garth
[     Ike][1] Heather  Jane
[    Jane][1] Ike

IGRAPH-UNW--V10-E18-------------------------------Krackhardt kite--
+ attributes: name (g/c), name (v/c), weight (e/n)
+ edges:
[   Andre][1] Beverly/1  Carol/3    Diane/3    Fernando/1
[ Beverly][1] Andre/1    Diane/1    Ed/2       Garth/2
[   Carol][1] Andre/2    Diane/2    Fernando/1
[   Diane][1] Andre/5    Beverly/1  Carol/0.4  Diane/2
[   Diane][5] Ed/1.5     Garth/2.5
[      Ed][1] Beverly/-1 Diane/1.5  Garth/2
[Fernando][1] Andre/1    Carol/2    Diane/1    Garth/1
[   Garth][1] Beverly/2  Diane/3    Ed/1       Fernando/-1
[ Heather][1] Fernando/3 Garth/1
[     Ike][1] Heather/1  Jane/-1
[    Jane][1] Ike/-2


IGRAPH-UNW--V10-E18-------------------------------Krackhardt kite--
+ attributes: name (g/c), name (v/c), weight (e/n)
+ edges:
[   Andre][1] Beverly (1)  Carol (3)    Diane (3)    Fernando (1)
[ Beverly][1] Andre (1)    Diane (1)    Ed (2)       Garth (2)
[   Carol][1] Andre (2)    Diane (2)    Fernando (1)
[   Diane][1] Andre (5)    Beverly (1)  Carol (0.5)  Diane (2)
[   Diane][5] Ed (1.5)     Garth (2.5)
[      Ed][1] Beverly (-1) Diane (1.5)  Garth (2)
[Fernando][1] Andre (1)    Carol (2)    Diane (1)    Garth (1)
[   Garth][1] Beverly (2)  Diane (3)    Ed (1)       Fernando (-1)
[ Heather][1] Fernando (3) Garth (1)
[     Ike][1] Heather (1)  Jane (-1)
[    Jane][1] Ike (-2)

IGRAPH UNW- V10 E18 -- Krackhardt kite
+ attr: name (g/c), name (v/c), weight (e/n)
+ edges:
[   Andre][1] Beverly (1)  Carol (3)    Diane (3)    Fernando (1)
[ Beverly][1] Andre (1)    Diane (1)    Ed (2)       Garth (2)
[   Carol][1] Andre (2)    Diane (2)    Fernando (1)
[   Diane][1] Andre (5)    Beverly (1)  Carol (0.5)  Diane (2)
[   Diane][5] Ed (1.5)     Garth (2.5)
[      Ed][1] Beverly (-1) Diane (1.5)  Garth (2)
[Fernando][1] Andre (1)    Carol (2)    Diane (1)    Garth (1)
[   Garth][1] Beverly (2)  Diane (3)    Ed (1)       Fernando (-1)
[ Heather][1] Fernando (3) Garth (1)
[     Ike][1] Heather (1)  Jane (-1)
[    Jane][1] Ike (-2)



IGRAPH-U----V100-E200----------------------------Gnm random graph--
+ edges:
[  1] 28 46 89 90
[  2] 47 69 72 89
[  3] 29
[  4] 17 20
[  5] 11 40 42 51 78 89
[  6] 27 32 70 87 93
[  7] 18 27 87
[  8] 18 24 82
[  9] 18 20 85 94
[ 10] 24 70 77 91
[ 11]  5 12 34 61 62
[ 12] 11 41 44 61 65 80
...

IGRAPH-U----100-200------------------------------Gnm random graph--
+ edges:
[  1] 28 46 89 90                 [  2] 47 69 72 89
[  3] 29                          [  4] 17 20
[  5] 11 40 42 51 78 89           [  6] 27 32 70 87 93
[  7] 18 27 87                    [  8] 18 24 82
[  9] 18 20 85 94                 [ 10] 24 70 77 91
[ 11]  5 12 34 61 62              [ 12] 11 41 44 61 65 80
...



"
