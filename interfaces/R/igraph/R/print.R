
## ----------------------------------------------------------------------
##
##   IGraph R package
##   Copyright (C) 2005-2014  Gabor Csardi <csardi.gabor@gmail.com>
##   334 Harvard street, Cambridge, MA 02139 USA
##
##   This program is free software; you can redistribute it and/or modify
##   it under the terms of the GNU General Public License as published by
##   the Free Software Foundation; either version 2 of the License, or
##   (at your option) any later version.
##
##   This program is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY; without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##   GNU General Public License for more details.
##
##   You should have received a copy of the GNU General Public License
##   along with this program; if not, write to the Free Software
##   Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA
##   02110-1301 USA
##
## ----------------------------------------------------------------------

###################################################################
# Convert graphs to human readable forms
###################################################################

.get.attr.codes <- function(object) {
  ga <- va <- ea <- ""
  gal <- graph_attr_names(object)
  if (length(gal) != 0) {
    ga <- paste(sep="", gal, " (g/",
                .Call("R_igraph_get_attr_mode", object, 2L, PACKAGE="igraph"),
                ")")
  }
  val <- vertex_attr_names(object)
  if (length(val) != 0) {
    va <- paste(sep="", val, " (v/",
                .Call("R_igraph_get_attr_mode", object, 3L, PACKAGE="igraph"),
                ")")
  }
  eal <- edge_attr_names(object)
  if (length(eal) != 0) {
    ea <- paste(sep="", edge_attr_names(object), " (e/",
                .Call("R_igraph_get_attr_mode", object, 4L, PACKAGE="igraph"),
                ")")
  }
  c(ga, va, ea)
}

.print.header <- function(object) {

  if (!is_igraph(object)) {
    stop("Not a graph object")
  }

  title <- paste(sep="", "IGRAPH ",
                 c("U","D")[is_directed(object)+1],
                 c("-","N")[is_named(object)+1],
                 c("-","W")[is_weighted(object)+1],
                 c("-","B")[is_bipartite(object)+1], " ",
                 vcount(object), " ", ecount(object), " -- ")
  w <- getOption("width")
  if (nchar(title) < w && "name" %in% graph_attr_names(object)) {
    title <- substring(paste(sep="", title,
                             as.character(object$name)[1]), 1, w-1)
  }
  cat(title, "\n", sep="")

  atxt <- .get.attr.codes(object)
  atxt <- paste(atxt[atxt!=""], collapse=", ")
  if (atxt != "") {
    atxt <- strwrap(paste(sep="", "+ attr: ", atxt), prefix = "| ",
                    initial = "")
    cat(atxt, sep="\n")
  }
  1 + if (length(atxt) == 1 && atxt == "") 0 else length(atxt)
}

#' @importFrom printr indent_print

.print.graph.attributes <- function(x, full, max.lines) {
  list <- graph_attr_names(x)
  if (length(list)!=0) {
    cat("+ graph attributes:\n")
    out <- capture.output({
      lapply(list, function(n) {
        cat(sep="", "+ ", n, ":\n")
        indent_print(graph_attr(x, n), .indent = "  ")
      })
      invisible(NULL)
    })
    indent_print(out, sep = "\n", .indent = "| ", .printer = cat)
    length(out) + 1
  } else {
    0
  }
}

## IGRAPH U--- 10 10 -- Ring graph
## + attr: name (g/c), mutual (g/l), circular (g/l)
## + graph attributes:
## | + name:
## |   [1] "Ring graph"
## | + mutual:
## |   [1] FALSE
## | + circular=
## |   [1] TRUE
## | + layout =
## |            [,1]          [,2]
## |    [1,]  0.000000  0.000000e+00
## |    [2,]  1.000000  0.000000e+00
## |    [3,]  0.809017  5.877853e-01
## |    [4,]  0.309017  9.510565e-01
## |    [5,] -0.309017  9.510565e-01
## |    [6,] -0.809017  5.877853e-01
## |    [7,] -1.000000  1.224647e-16
## |    [8,] -0.809017 -5.877853e-01
## |    [9,] -0.309017 -9.510565e-01
## |   [10,]  0.309017 -9.510565e-01
## |   [11,]  0.809017 -5.877853e-01
## + edges:
##  [1] 1-- 2 2-- 3 3-- 4 4-- 5 5-- 6 6-- 7 7-- 8 8-- 9 9--10 1--10

#' @importFrom printr indent_print

.print.vertex.attributes <- function(x, full, max.lines) {
  pf <- function(x) .print.vertex.attributes.old(x, full, max.lines)
  if (length(vertex_attr_names(x))) cat("+ vertex attributes:\n")
  indent_print(x, .indent = "| ", .printer = pf)
}

.print.vertex.attributes.old <- function(x, full, max.lines) {
  vc <- vcount(x)
  list <- vertex_attr_names(x)
  if (length(list) != 0) {
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
                   is.numeric(vertex_attr(x, v)) ||
                   is.character(vertex_attr(x, v)) ||
                   is.logical(vertex_attr(x, v))))) {
      ## create a table
      tab <- data.frame(v=paste(sep="", "[", ind, "]"), row.names="v")
      for (i in list) {
        tab[i] <- vertex_attr(x, i, ind)
      }
      print(tab)
    } else {
      for (i in ind) {
        cat(sep="", "[[", i, "]]\n")
        lapply(list, function(n) {
          cat(sep="", "[[", i, "]][[", n, "]]\n")
          print(vertex_attr(x, n, i))})
      }
    }
    options(max.print=mp)
    if (omitted.vertices != 0) {
      cat(paste('[ reached getOption("max.print") -- omitted',
                omitted.vertices, "vertices ]\n\n"))
    }      
  }
}

.print.edges.edgelist <- function(x, edges = E(x), names) {
  ec <- length(edges)
  list <- edge_attr_names(x)
  list <- list[list!="name"]
  arrow <- ifelse(is_directed(x), "->", "--")
  if (is_named(x)) {
    cat("+ edges (vertex names) and their attributes:\n")
  } else {
    cat("+ edges and their attributes:\n")
  }
  if (names && ! "name" %in% vertex_attr_names(x)) {
    names <- FALSE
  }
  if (names && "name" %in% vertex_attr_names(x) &&
      !is.numeric(vertex_attr(x, "name")) &&
      !is.character(vertex_attr(x, "name")) &&
      !is.logical(vertex_attr(x, "name"))) {
    warning("Can't print vertex names, complex `name' vertex attribute")
    names <- FALSE
  }
  
  mp <- getOption("max.print")
  if (mp >= ec) {
    omitted.edges <- 0
    el <- ends(x, edges, names=names)
  } else {
    omitted.edges <- ec-mp
    el <- ends(x, ends[seq_len(mp)])
    if (names) { el[] <- V(x)$name[el] }
  }
  ename <- if ("name" %in% edge_attr_names(x)) {
    paste(sep="", "'", E(x)$name, "'")
  } else {
    seq(length=nrow(el))
  }
  if (ec==0 || 
      all(sapply(list, function(v) is.numeric(edge_attr(x, v)) |
                 is.character(edge_attr(x,v)) |
                 is.logical(edge_attr(x, v))))) {
    ## create a table
    tab <- data.frame(row.names=paste(sep="", "[", ename, "]"))
    if (is.numeric(el)) { w <- nchar(max(el)) } else { w <- max(nchar(el)) }
    tab["edge"] <- paste(sep="", format(el[,1], width=w),
                         arrow, format(el[,2], width=w))
    for (i in list) {
      tab[i] <- edge_attr(x, i)
    }
    print(tab)
  } else {
    i <- 1
    apply(el, 1, function(v) {
      cat(sep="", "[", ename[i], "] ", v[1], " ", arrow, " ", v[2]);
      lapply(list, function(n) {
        cat(sep="", "\n[[", i, "]][[", n, "]]\n")
        print(edge_attr(x, n, i))})
      cat("\n")
      i <<- i+1
    })
  }
  if (omitted.edges != 0) {
    cat(paste('[ reached getOption("max.print") -- omitted', omitted.edges,
              'edges ]\n\n'))
  }    
}

#' @importFrom printr head_print printer_callback

.print.edges.compressed <- function(x, edges = E(x), names, num = FALSE,
                                      max.lines = igraph_opt("auto.print.lines")) {

  title <- "+" %+%
    (if (num) " " %+% chr(length(edges)) %+% "/" %+%
       (if (is.null(x)) "?" else chr(ecount(x))) else "") %+%
    " edges" %+%
    (if (is.null(x)) ", unknown graph" else "") %+%
    (if (is.null(attr(edges, "vnames"))) "" else " (vertex names)") %+%
    ":\n"
  cat(title)

  if (is.null(max.lines)) {
    .print.edges.compressed.all(x, edges, names)
  } else {
    .print.edges.compressed.limit(x, edges, names, max.lines)
  }
}

.print.edges.compressed.all <- function(x, edges, names) {

  arrow <- c("--", "->")[is_directed(x)+1]

  if (!is.null(x)) {
    el <- ends(x, edges, names=names)
    pr <- paste(sep="", format(el[,1]), arrow, format(el[,2]))
    print(pr, quote=FALSE)
  } else {
    if (!is.null(attr(edges, "vnames"))) {
      print(as.vector(attr(edges, "vnames")), quote = FALSE)
    } else if (!is.null(names(edges))) {
      print(names(edges), quote = FALSE)
    } else {
      print(as.vector(edges))
    }
  }

}

.print.edges.compressed.limit <- function(x, edges, names, max.lines) {

  if (!is.null(x)) {

    arrow <- c("--", "->")[is_directed(x)+1]

    can_max <- NA
    el <- NA

    fun <- function(q, no) {
      if (q == "length") {
        length(edges)
      } else if (q == "min_width") {
        5
      } else if (q == "width") {
        el <<- ends(x, edges[seq_len(no)], names = names)
        cummax(nchar(el[,1])) + nchar(arrow) + cummax(nchar(el[,2])) + 1
      } else if (q == "print") {
        el <<- el[seq_len(no), , drop = FALSE]
        out <- paste(sep="", format(el[,1]), arrow, format(el[,2]))
        capture.output(print(out, quote = FALSE))
      } else if (q == "max") {
        can_max <<- no
      } else if (q == "done") {
        if (no["tried_items"] < length(edges) ||
            no["printed_lines"] < no["tried_lines"]) {
          cat("+ ... omitted several edges\n")
        }
      }
    }

    fun <- printer_callback(fun)
    head_print(fun, max_lines = max.lines)
  } else {
    if (!is.null(attr(edges, "vnames"))) {
      head_print(as.vector(attr(edges, "vnames")), quote = FALSE)
    } else if (!is.null(names(edges))) {
      head_print(names(edges), quote = FALSE)
    } else {
      head_print(as.vector(edges))
    }
  }
}

.print.edges.adjlist <- function(x) {
  ## TODO: getOption("max.print")
  cat("+ edges:\n")
  vc <- vcount(x)
  arrow <- c(" -- ", " -> ")[is_directed(x)+1]
  al <- as_adj_list(x, mode="out")
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

.print.edges.adjlist.named <- function(x, edges = E(x)) {
  ## TODO getOption("max.print")
  cat("+ edges (vertex names):\n")

  arrow <- c(" -- ", " -> ")[is_directed(x)+1]
  vn <- V(x)$name

  al <- as_adj_list(x, mode="out")
  alstr <- sapply(al, function(x) { paste(collapse=", ", vn[x]) })
  alstr <- paste(sep="", format(vn), arrow, alstr)
  alstr <- strwrap(alstr, exdent=max(nchar(vn))+nchar(arrow))
  cat(alstr, sep="\n")
}

#' @method str igraph
#' @export

str.igraph <- function(object, ...) {
  print.igraph(object, full=TRUE, ...)
}



#' Print graphs to the terminal
#' 
#' These functions attempt to print a graph to the terminal in a human readable
#' form.
#' 
#' \code{summary.igraph} prints the number of vertices, edges and whether the
#' graph is directed.
#' 
#' \code{str.igraph} prints the same information, and also lists the edges, and
#' optionally graph, vertex and/or edge attributes.
#' 
#' \code{print.igraph} behaves either as \code{summary.igraph} or
#' \code{str.igraph} depending on the \code{full} argument. See also the
#' \sQuote{print.full} igraph option and \code{\link{igraph_opt}}.
#' 
#' The graph summary printed by \code{summary.igraph} (and \code{print.igraph}
#' and \code{str.igraph}) consists one or more lines. The first line contains
#' the basic properties of the graph, and the rest contains its attributes.
#' Here is an example, a small star graph with weighed directed edges and named
#' vertices: \preformatted{    IGRAPH DNW- 10 9 -- In-star
#'     + attr: name (g/c), mode (g/c), center (g/n), name (v/c),
#'       weight (e/n) }
#' The first line always
#' starts with \code{IGRAPH}, showing you that the object is an igraph graph.
#' Then a four letter long code string is printed. The first letter
#' distinguishes between directed (\sQuote{\code{D}}) and undirected
#' (\sQuote{\code{U}}) graphs. The second letter is \sQuote{\code{N}} for named
#' graphs, i.e. graphs with the \code{name} vertex attribute set. The third
#' letter is \sQuote{\code{W}} for weighted graphs, i.e. graphs with the
#' \code{weight} edge attribute set. The fourth letter is \sQuote{\code{B}} for
#' bipartite graphs, i.e. for graphs with the \code{type} vertex attribute set.
#' 
#' Then, after two dashes, the name of the graph is printed, if it has one,
#' i.e. if the \code{name} graph attribute is set.
#' 
#' From the second line, the attributes of the graph are listed, separated by a
#' comma. After the attribute names, the kind of the attribute -- graph
#' (\sQuote{\code{g}}), vertex (\sQuote{\code{v}}) or edge (\sQuote{\code{e}})
#' -- is denoted, and the type of the attribute as well, character
#' (\sQuote{\code{c}}), numeric (\sQuote{\code{n}}), logical
#' (\sQuote{\code{l}}), or other (\sQuote{\code{x}}).
#' 
#' As of igraph 0.4 \code{str.igraph} and \code{print.igraph} use the
#' \code{max.print} option, see \code{\link[base]{options}} for details.
#' 
#' @aliases print.igraph str.igraph summary.igraph
#' @param x The graph to print.
#' @param full Logical scalar, whether to print the graph structure itself as
#' well.
#' @param graph.attributes Logical constant, whether to print graph attributes.
#' @param vertex.attributes Logical constant, whether to print vertex
#' attributes.
#' @param edge.attributes Logical constant, whether to print edge attributes.
#' @param names Logical constant, whether to print symbolic vertex names (ie.
#' the \code{name} vertex attribute) or vertex ids.
#' @param object The graph of which the summary will be printed.
#' @param \dots Additional agruments.
#' @return All these functions return the graph invisibly.
#' @author Gabor Csardi \email{csardi.gabor@@gmail.com}
#' @method print igraph
#' @export
#' @keywords graphs
#' @examples
#' 
#' g <- make_ring(10)
#' g
#' summary(g)
#' 
print.igraph <- function(x, full=igraph_opt("print.full"),
                graph.attributes=igraph_opt("print.graph.attributes"),
                vertex.attributes=igraph_opt("print.vertex.attributes"),
                edge.attributes=igraph_opt("print.edge.attributes"),
                names=TRUE, max.lines = igraph_opt("auto.print.lines"), ...) {
  
  if (!is_igraph(x)) {
    stop("Not a graph object")
  }

  head_lines <- .print.header(x)
  if (is.logical(full) && full) {
    if (graph.attributes) {
      head_lines <- head_lines + .print.graph.attributes(x, full, max.lines)
    }
    if (vertex.attributes) {
      head_lines <- head_lines + .print.vertex.attributes(x, full, max.lines)
    }
    if (ecount(x)==0) {
      ## Do nothing
    } else if (edge.attributes && length(edge_attr_names(x)) != 0 ) {
      .print.edges.edgelist(x, names = names)
    } else if (median(degree(x, mode="out")) < 3) {
      .print.edges.compressed(x, names = names)
    } else if (is_named(x)) {
      .print.edges.adjlist.named(x)
    } else {
      .print.edges.adjlist(x)
    }
  } else if (full == "auto") {
    .print.edges.compressed(x, names = names, max.lines =
                              max.lines - head_lines)
  }
  
  invisible(x)
}

#' @rdname print.igraph
#' @method summary igraph
#' @export

summary.igraph <- function(object, ...) {
  .print.header(object)
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
