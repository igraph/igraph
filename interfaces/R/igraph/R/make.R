
## ----------------------------------------------------------------
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
## -----------------------------------------------------------------

#' Make a new graph
#'
#' This is is generic function for creating graphs.
#'
#' @details
#' TODO
#'
#'
#' @param ... Parameters, see details below.
#'
#' @importFrom lazyeval lazy_eval
#' @export
#' @examples
#' r <- make_(ring(10))
#' l <- make_(lattice(c(3, 3, 3)))

make_ <- function(...) {

  me <- attr(sys.function(), "name") %||% "construct"
  args <- list(...)
  cidx <- vapply(args, inherits, TRUE, what = "igraph_constructor_spec")
  if (sum(cidx) == 0) {
    stop("Don't know how to ", me, ", nothing given")
  }
  if (sum(cidx) > 1) {
    stop("Don't know how to ", me, ", multiple constructors given")
  }
  cons <- args[ cidx][[1]]
  args <- args[!cidx]

  args2 <- if (cons$lazy) lapply(cons$args, "[[", "expr") else lazy_eval(cons$args)

  do_call(cons$fun, args2, args)

}

#' Sample from a random graph model
#'
#' Generic function for sampling from network models.
#'
#' @details
#' TODO
#'
#' @param ... Parameters, see details below.
#'
#' @export
#' @examples
#' pref_matrix <- cbind(c(0.8, 0.1), c(0.1, 0.7))
#' blocky <- sample_(sbm(n = 20, pref.matrix = pref_matrix,
#'   block.sizes = c(10, 10)))
#'
#' blocky2 <- pref_matrix %>%
#'   sample_sbm(n = 20, block.sizes = c(10, 10))
#'
#' ## Arguments are passed on from sample_ to sample_sbm
#' blocky3 <- pref_matrix %>%
#'   sample_(sbm(), n = 20, block.sizes = c(10, 10))

sample_ <- make_

#' Convert object to a graph
#'
#' This is a generic function to convert R objects to igraph graphs.
#'
#' @details
#' TODO
#'
#' @param ... Parameters, see details below.
#'
#' @importFrom lazyeval lazy_dots
#' @export
#' @examples
#' ## These are equivalent
#' graph_(cbind(1:5,2:6), from_edgelist(directed = FALSE))
#' graph_(cbind(1:5,2:6), from_edgelist(), directed = FALSE)

graph_ <- make_

attr(make_, "name") <- "make_"
attr(sample_, "name") <- "sample_"
attr(graph_, "name") <- "graph_"

constructor_spec <- function(fun, ..., .lazy = FALSE) {
  structure(
    list(
      fun = fun,
      args = lazy_dots(...),
      lazy = .lazy
    ),
    class = "igraph_constructor_spec"
  )
}

## -----------------------------------------------------------------

#' Create an igraph graph from a list of edges, or a notable graph
#'
#' @section Notable graphs:
#'
#' \code{make_graph} can create some notable graphs. The name of the
#' graph (case insensitive), a character scalar must be suppliced as
#' the \code{edges} argument, and other arguments are ignored. (A warning
#' is given is they are specified.)
#'
#' \code{make_graph} knows the following graphs: \describe{
#'   \item{Bull}{The bull graph, 5 vertices, 5 edges, resembles to the head
#'     of a bull if drawn properly.}
#'   \item{Chvatal}{This is the smallest triangle-free graph that is
#'     both 4-chromatic and 4-regular. According to the Grunbaum conjecture there
#'     exists an m-regular, m-chromatic graph with n vertices for every m>1 and
#'     n>2. The Chvatal graph is an example for m=4 and n=12. It has 24 edges.}
#'   \item{Coxeter}{A non-Hamiltonian cubic symmetric graph with 28 vertices and
#'     42 edges.}
#'   \item{Cubical}{The Platonic graph of the cube. A convex regular
#'     polyhedron with 8 vertices and 12 edges.}
#'   \item{Diamond}{A graph with 4 vertices and 5 edges, resembles to a
#'     schematic diamond if drawn properly.}
#'   \item{Dodecahedral, Dodecahedron}{Another Platonic solid with 20 vertices
#'     and 30 edges.}
#'   \item{Folkman}{The semisymmetric graph with minimum number of
#'     vertices, 20 and 40 edges. A semisymmetric graph is regular, edge transitive
#'     and not vertex transitive.}
#'   \item{Franklin}{This is a graph whose embedding
#'     to the Klein bottle can be colored with six colors, it is a counterexample
#'     to the neccessity of the Heawood conjecture on a Klein bottle. It has 12
#'     vertices and 18 edges.}
#'   \item{Frucht}{The Frucht Graph is the smallest
#'     cubical graph whose automorphism group consists only of the identity
#'     element. It has 12 vertices and 18 edges.}
#'   \item{Grotzsch}{The Grötzsch
#'     graph is a triangle-free graph with 11 vertices, 20 edges, and chromatic
#'     number 4. It is named after German mathematician Herbert Grötzsch, and its
#'     existence demonstrates that the assumption of planarity is necessary in
#'     Grötzsch's theorem that every triangle-free planar graph is 3-colorable.}
#'   \item{Heawood}{The Heawood graph is an undirected graph with 14 vertices and
#'     21 edges. The graph is cubic, and all cycles in the graph have six or more
#'     edges. Every smaller cubic graph has shorter cycles, so this graph is the
#'     6-cage, the smallest cubic graph of girth 6.}
#'   \item{Herschel}{The Herschel
#'     graph is the smallest nonhamiltonian polyhedral graph. It is the unique such
#'     graph on 11 nodes, and has 18 edges.}
#'   \item{House}{The house graph is a
#'     5-vertex, 6-edge graph, the schematic draw of a house if drawn properly,
#'     basicly a triangle of the top of a square.}
#'   \item{HouseX}{The same as the
#'     house graph with an X in the square. 5 vertices and 8 edges.}
#'   \item{Icosahedral, Icosahedron}{A Platonic solid with 12 vertices and 30
#'     edges.}
#'   \item{Krackhardt kite}{A social network with 10 vertices and 18
#'     edges.  Krackhardt, D. Assessing the Political Landscape: Structure,
#'     Cognition, and Power in Organizations.  Admin. Sci. Quart. 35, 342-369,
#'     1990.}
#'   \item{Levi}{The graph is a 4-arc transitive cubic graph, it has 30
#'     vertices and 45 edges.}
#'   \item{McGee}{The McGee graph is the unique 3-regular
#'     7-cage graph, it has 24 vertices and 36 edges.}
#'   \item{Meredith}{The Meredith
#'     graph is a quartic graph on 70 nodes and 140 edges that is a counterexample
#'     to the conjecture that every 4-regular 4-connected graph is Hamiltonian.}
#'   \item{Noperfectmatching}{A connected graph with 16 vertices and 27 edges
#'     containing no perfect matching. A matching in a graph is a set of pairwise
#'     non-adjacent edges; that is, no two edges share a common vertex. A perfect
#'     matching is a matching which covers all vertices of the graph.}
#'   \item{Nonline}{A graph whose connected components are the 9 graphs whose
#'     presence as a vertex-induced subgraph in a graph makes a nonline graph. It
#'     has 50 vertices and 72 edges.}
#'   \item{Octahedral, Octahedron}{Platonic solid
#'     with 6 vertices and 12 edges.}
#'   \item{Petersen}{A 3-regular graph with 10
#'     vertices and 15 edges. It is the smallest hypohamiltonian graph, ie. it is
#'     non-hamiltonian but removing any single vertex from it makes it
#'     Hamiltonian.}
#'   \item{Robertson}{The unique (4,5)-cage graph, ie. a 4-regular
#'    graph of girth 5. It has 19 vertices and 38 edges.}
#'   \item{Smallestcyclicgroup}{A smallest nontrivial graph whose automorphism
#'     group is cyclic. It has 9 vertices and 15 edges.}
#'   \item{Tetrahedral,
#'     Tetrahedron}{Platonic solid with 4 vertices and 6 edges.}
#'   \item{Thomassen}{The smallest hypotraceable graph, on 34 vertices and 52
#'     edges. A hypotracable graph does not contain a Hamiltonian path but after
#'     removing any single vertex from it the remainder always contains a
#'     Hamiltonian path. A graph containing a Hamiltonian path is called tracable.}
#'   \item{Tutte}{Tait's Hamiltonian graph conjecture states that every
#'     3-connected 3-regular planar graph is Hamiltonian.  This graph is a
#'     counterexample. It has 46 vertices and 69 edges.}
#'   \item{Uniquely3colorable}{Returns a 12-vertex, triangle-free graph with
#'     chromatic number 3 that is uniquely 3-colorable.}
#'   \item{Walther}{An identity
#'     graph with 25 vertices and 31 edges. An identity graph has a single graph
#'     automorphism, the trivial one.}
#'   \item{Zachary}{Social network of friendships
#'     between 34 members of a karate club at a US university in the 1970s. See W.
#'     W. Zachary, An information flow model for conflict and fission in small
#'     groups, Journal of Anthropological Research 33, 452-473 (1977).  } }
#'
#' @encoding UTF-8
#' @aliases graph.famous graph
#' @param edges A vector defining the edges, the first edge points
#'   from the first element to the second, the second edge from the third
#'   to the fourth, etc. For a numeric vector, these are interpreted
#'   as internal vertex ids. For character vectors, they are interpreted
#'   as vertex names.
#'
#'   Alternatively, this can be a character scalar, the name of a
#'   notable graph. See Notable graphs below. The name is case
#'   insensitive.
#' @param n The number of vertices in the graph. This argument is
#'   ignored (with a warning) if \code{edges} are symbolic vertex names. It
#'   is also ignored if there is a bigger vertex id in \code{edges}. This
#'   means that for this function it is safe to supply zero here if the
#'   vertex with the largest id is not an isolate.
#' @param isolates Character vector, names of isolate vertices,
#'   for symbolic edge lists. It is ignored for numeric edge lists.
#' @param directed Whether to create a directed graph.
#' @return An igraph graph.
#'
#' @family determimistic constructors
#' @export
#' @examples
#' make_graph(c(1, 2, 2, 3, 3, 4, 5, 6), directed = FALSE)
#' make_graph(c("A", "B", "B", "C", "C", "D"), directed = FALSE)
#'
#' solids <- list(make_graph("Tetrahedron"),
#'                make_graph("Cubical"),
#'                make_graph("Octahedron"),
#'                make_graph("Dodecahedron"),
#'                make_graph("Icosahedron"))

make_graph <- function(edges, n=max(edges), isolates = NULL, directed=TRUE) {

  if (is.character(edges) && length(edges) == 1) {
    if (!missing(n)) warning("'n' is ignored for the '", edges, "' graph")
    if (!missing(isolates)) {
      warning("'isolates' is ignored for the '", edges, "' graph")
    }
    if (!missing(directed)) {
      warning("'directed' is ignored for the '", edges, "' graph")
    }
    make_famous_graph(edges)

  } else if (is.numeric(edges)) {
    if (!is.null(isolates)) {
      warning("'isolates' ignored for numeric edge list")
    }
    on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
    .Call("R_igraph_create", as.numeric(edges)-1, as.numeric(n),
          as.logical(directed),
          PACKAGE="igraph")

  } else if (is.character(edges)) {
    if (!missing(n)) {
      warning("'n' is ignored for edge list with vertex names")
    }
    el <- matrix(edges, ncol = 2, byrow = TRUE)
    res <- graph_from_edgelist(el, directed = directed)
    if (!is.null(isolates)) {
      isolates <- as.character(isolates)
      res <- res + vertices(isolates)
    }
    res

  } else {
    stop("'edges' must be numeric of character")
  }
}

make_famous_graph <- function(name) {

  name <- gsub("\\s", "_", name)

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  res <- .Call("R_igraph_famous", as.character(name),
               PACKAGE="igraph")
  if (igraph_opt("add.params")) {
    res$name <- capitalize(name)
  }
  res
}

#' @rdname make_graph
#' @export

make_directed_graph <- function(edges, n = max(edges)) {
  if (missing(n)) {
    make_graph(edges, directed = TRUE)
  } else {
    make_graph(edges, n = n, directed = TRUE)
  }
}

#' @rdname make_graph
#' @export

make_undirected_graph <- function(edges, n = max(edges)) {
  if (missing(n)) {
    make_graph(edges, directed = FALSE)
  } else {
    make_graph(edges, n = n, directed = FALSE)
  }
}

#' @rdname make_graph
#' @param ... Passed to \code{make_directed_graph} or
#'   \code{make_undirected_graph}.
#' @export

directed_graph <- function(...) constructor_spec(make_directed_graph, ...)

#' @rdname make_graph
#' @export

undirected_graph <- function(...) constructor_spec(make_undirected_graph, ...)

## -----------------------------------------------------------------

#' A graph with no edges
#'
#' @aliases graph.empty
#' @concept Empty graph.
#' @param n Number of vertices.
#' @param directed Whether to create a directed graph.
#' @return An igraph graph.
#'
#' @family determimistic constructors
#' @export
#' @examples
#' make_empty_graph(n = 10)
#' make_empty_graph(n = 5, directed = FALSE)

make_empty_graph <- function(n=0, directed=TRUE) {
  # Argument checks
  n <- as.integer(n)
  directed <- as.logical(directed)

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  # Function call
  res <- .Call("R_igraph_empty", n, directed,
        PACKAGE="igraph")

  res
}

#' @rdname make_empty_graph
#' @param ... Passed to \code{make_graph_empty}.
#' @export

empty_graph <- function(...) constructor_spec(make_empty_graph, ...)

## -----------------------------------------------------------------


#' Creating (small) graphs via a simple interface
#'
#' This function is useful if you want to create a small (named) graph
#' quickly, it works for both directed and undirected graphs.
#'
#' @details
#' \code{graph_from_literal} is very handy for creating small graphs quickly.
#' You need to supply one or more R expressions giving the structure of
#' the graph. The expressions consist of vertex names and edge
#' operators. An edge operator is a sequence of \sQuote{\code{-}} and
#' \sQuote{\code{+}} characters, the former is for the edges and the
#' latter is used for arrow heads. The edges can be arbitrarily long,
#' ie. you may use as many \sQuote{\code{-}} characters to \dQuote{draw}
#' them as you like.
#'
#' If all edge operators consist of only \sQuote{\code{-}} characters
#' then the graph will be undirected, whereas a single \sQuote{\code{+}}
#' character implies a directed graph.
#'
#' Let us see some simple examples. Without arguments the function
#' creates an empty graph:
#' \preformatted{  graph_from_literal()
#' }
#'
#' A simple undirected graph with two vertices called \sQuote{A} and
#' \sQuote{B} and one edge only:
#' \preformatted{  graph_from_literal(A-B)
#' }
#'
#' Remember that the length of the edges does not matter, so we could
#' have written the following, this creates the same graph:
#' \preformatted{  graph_from_literal( A-----B )
#' }
#'
#' If you have many disconnected components in the graph, separate them
#' with commas. You can also give isolate vertices.
#' \preformatted{  graph_from_literal( A--B, C--D, E--F, G--H, I, J, K )
#' }
#'
#' The \sQuote{\code{:}} operator can be used to define vertex sets. If
#' an edge operator connects two vertex sets then every vertex from the
#' first set will be connected to every vertex in the second set. The
#' following form creates a full graph, including loop edges:
#' \preformatted{  graph_from_literal( A:B:C:D -- A:B:C:D )
#' }
#'
#' In directed graphs, edges will be created only if the edge operator
#' includes a arrow head (\sQuote{+}) \emph{at the end} of the edge:
#' \preformatted{  graph_from_literal( A -+ B -+ C )
#'   graph_from_literal( A +- B -+ C )
#'   graph_from_literal( A +- B -- C )
#' }
#' Thus in the third example no edge is created between vertices \code{B}
#' and \code{C}.
#'
#' Mutual edges can be also created with a simple edge operator:
#' \preformatted{  graph_from_literal( A +-+ B +---+ C ++ D + E)
#' }
#' Note again that the length of the edge operators is arbitrary,
#' \sQuote{\code{+}}, \sQuote{\code{++}} and \sQuote{\code{+-----+}} have
#' exactly the same meaning.
#'
#' If the vertex names include spaces or other special characters then
#' you need to quote them:
#' \preformatted{  graph_from_literal( "this is" +- "a silly" -+ "graph here" )
#' }
#' You can include any character in the vertex names this way, even
#' \sQuote{+} and \sQuote{-} characters.
#'
#' See more examples below.
#'
#' @aliases graph.formula
#' @param ... For \code{graph_from_literal} the formulae giving the
#'   structure of the graph, see details below. For \code{from_literal}
#'   all arguments are passed to \code{graph_from_literal}.
#' @param simplify Logical scalar, whether to call \code{\link{simplify}}
#'   on the created graph. By default the graph is simplified, loop and
#'   multiple edges are removed.
#' @return An igraph graph
#'
#' @family determimistic constructors
#' @export
#' @examples
#' # A simple undirected graph
#' g <- graph_from_literal( Alice-Bob-Cecil-Alice, Daniel-Cecil-Eugene,
#'                      Cecil-Gordon )
#' g
#'
#' # Another undirected graph, ":" notation
#' g2 <- graph_from_literal( Alice-Bob:Cecil:Daniel, Cecil:Daniel-Eugene:Gordon )
#' g2
#'
#' # A directed graph
#' g3 <- graph_from_literal( Alice +-+ Bob --+ Cecil +-- Daniel,
#'                      Eugene --+ Gordon:Helen )
#' g3
#'
#' # A graph with isolate vertices
#' g4 <- graph_from_literal( Alice -- Bob -- Daniel, Cecil:Gordon, Helen )
#' g4
#' V(g4)$name
#'
#' # "Arrows" can be arbitrarily long
#' g5 <- graph_from_literal( Alice +---------+ Bob )
#' g5
#'
#' # Special vertex names
#' g6 <- graph_from_literal( "+" -- "-", "*" -- "/", "%%" -- "%/%" )
#' g6
#'

graph_from_literal <- function(..., simplify=TRUE) {
  mf <- as.list(match.call())[-1]

  ## In case 'simplify' is given
  if ('simplify' %in% names(mf)) {
    w <- which(names(mf)=='simplify')
    if (length(w) > 1) { stop("'simplify' specified multiple times") }
    mf <- mf[-w]
  }

  ## Operators first
  f <- function(x) {
    if (is.call(x)) {
      return (list(as.character(x[[1]]), lapply(x[-1], f)))
    } else {
      return (NULL)
    }
  }
  ops <- unlist(lapply(mf, f))
  if (all(ops %in% c("-", ":"))) {
    directed <- FALSE
  } else if (all(ops %in% c("-", "+", ":"))) {
    directed <- TRUE
  } else {
    stop("Invalid operator in formula")
  }

  f <- function(x) {
    if (is.call(x)) {
      if (length(x)==3) {
        return( list(f(x[[2]]), op=as.character(x[[1]]), f(x[[3]])) )
      } else {
        return( list(op=as.character(x[[1]]), f(x[[2]])) )
      }
    } else {
      return( c(sym=as.character(x)) )
    }
  }

  ret <- lapply(mf, function(x) unlist(f(x)))

  v <- unique(unlist(lapply(ret, function(x) { x[ names(x)=="sym" ] })))

  ## Merge symbols for ":"
  ret <- lapply(ret, function(x) {
    res <- list()
    for (i in seq(along=x)) {
      if (x[i]==":" && names(x)[i]=="op") {
        ## SKIP
      } else if (i>1 && x[i-1]==":" && names(x)[i-1]=="op") {
        res[[length(res)]] <- c(res[[length(res)]], unname(x[i]))
      } else {
        res <- c(res, x[i])
      }
    }
    res
  })

  ## Ok, create the edges
  edges <- numeric()
  for (i in seq(along=ret)) {
    prev.sym <- character()
    lhead <- rhead <- character()
    for (j in seq(along=ret[[i]])) {
      act <- ret[[i]][[j]]
      if (names(ret[[i]])[j]=="op") {
        if (length(lhead)==0) {
          lhead <- rhead <- act
        } else {
          rhead <- act
        }
      } else if (names(ret[[i]])[j]=="sym") {
        for (ps in prev.sym) {
          for (ps2 in act) {
            if (lhead=="+") {
              edges <- c(edges, unname(c(ps2, ps)))
            }
            if (!directed || rhead=="+") {
              edges <- c(edges, unname(c(ps, ps2)))
            }
          }
        }
        lhead <- rhead <- character()
        prev.sym <- act
      }
    }
  }

  ids <- seq(along=v)
  names(ids) <- v
  res <- graph( unname(ids[edges]), n=length(v), directed=directed)
  if (simplify) res <- simplify(res)
  res <- set_vertex_attr(res, "name", value=v)
  res
}

#' @rdname graph_from_literal
#' @export

from_literal <- function(...)
  constructor_spec(graph_from_literal, ..., .lazy = TRUE)

## -----------------------------------------------------------------

#' Create a star graph, a tree with n vertices and n - 1 leaves
#'
#' \code{star} creates a star graph, in this every single vertex is
#' connected to the center vertex and nobody else.
#'
#' @aliases graph.star
#' @concept Star graph
#' @param n Number of vertices.
#' @param mode It defines the direction of the
#'   edges, \code{in}: the edges point \emph{to} the center, \code{out}:
#'   the edges point \emph{from} the center, \code{mutual}: a directed
#'   star is created with mutual edges, \code{undirected}: the edges
#'   are undirected.
#' @param center ID of the center vertex.
#' @return An igraph graph.
#'
#' @family determimistic constructors
#' @export
#' @examples
#' make_star(10, mode = "out")
#' make_star(5, mode = "undirected")

make_star <- function(n, mode=c("in", "out", "mutual", "undirected"),
                   center=1 ) {

  mode <- igraph.match.arg(mode)
  mode1 <- switch(mode, "out"=0, "in"=1, "undirected"=2, "mutual"=3)

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  res <- .Call("R_igraph_star", as.numeric(n), as.numeric(mode1),
               as.numeric(center)-1,
               PACKAGE="igraph")
  if (igraph_opt("add.params")) {
    res$name <- switch(mode, "in"="In-star", "out"="Out-star", "Star")
    res$mode <- mode
    res$center <- center
  }
  res
}

#' @rdname make_star
#' @param ... Passed to \code{make_star}.
#' @export

star <- function(...) constructor_spec(make_star, ...)

## -----------------------------------------------------------------

#' Create a full graph
#'
#' @aliases graph.full
#' @concept Full graph
#' @param n Number of vertices.
#' @param directed Whether to create a directed graph.
#' @param loops Whether to add self-loops to the graph.
#' @return An igraph graph
#'
#' @family determimistic constructors
#' @export
#' @examples
#' make_full_graph(5)
#' str(make_full_graph(4, directed = TRUE))

make_full_graph <- function(n, directed=FALSE, loops=FALSE) {
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  res <- .Call("R_igraph_full", as.numeric(n), as.logical(directed),
               as.logical(loops),
               PACKAGE="igraph")
  if (igraph_opt("add.params")) {
    res$name <- "Full graph"
    res$loops <- loops
  }
  res
}

#' @rdname make_full_graph
#' @param ... Passed to \code{make_full_graph}.
#' @export

full_graph <- function(...) constructor_spec(make_full_graph, ...)

## -----------------------------------------------------------------

#' Create a lattice graph
#'
#' \code{make_lattice} is a flexible function, it can create lattices of
#' arbitrary dimensions, periodic or unperiodic ones. It has two
#' forms. In the first form you only supply \code{dimvector}, but not
#' \code{length} and \code{dim}. In the second form you omit
#' \code{dimvector} and supply \code{length} and \code{dim}.
#'
#' @aliases graph.lattice
#' @concept Lattice
#' @param dimvector A vector giving the size of the lattice in each
#'   dimension.
#' @param length Integer constant, for regular lattices, the size of the
#'   lattice in each dimension.
#' @param dim Integer constant, the dimension of the lattice.
#' @param nei The distance within which (inclusive) the neighbors on the
#'   lattice will be connected. This parameter is not used right now.
#' @param directed Whether to create a directed lattice.
#' @param mutual Logical, if \code{TRUE} directed lattices will be
#'   mutually connected.
#' @param circular Logical, if \code{TRUE} the lattice or ring will be
#'   circular.
#' @return An igraph graph.
#'
#' @family determimistic constructors
#' @export
#' @examples
#' make_lattice(c(5, 5, 5))
#' make_lattice(length = 5, dim = 3)

make_lattice <- function(dimvector = NULL, length = NULL, dim = NULL,
                          nei = 1, directed = FALSE, mutual = FALSE,
                          circular=FALSE) {

  if (is.null(dimvector)) {
    dimvector <- rep(length, dim)
  }

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  res <- .Call("R_igraph_lattice", as.numeric(dimvector), as.numeric(nei),
               as.logical(directed), as.logical(mutual),
               as.logical(circular),
               PACKAGE="igraph")
  if (igraph_opt("add.params")) {
    res$name <- "Lattice graph"
    res$dimvector <- dimvector
    res$nei <- nei
    res$mutual <- mutual
    res$circular <- circular
  }
  res
}

#' @rdname make_lattice
#' @param ... Passed to \code{make_lattice}.
#' @export

lattice <- function(...) constructor_spec(make_lattice, ...)

## -----------------------------------------------------------------

#' Create a ring graph
#'
#' A ring is a one-dimensional lattice and this function is a special case
#' of \code{\link{make_lattice}}.
#'
#' @aliases make_ring graph.ring
#' @param n Number of vertices.
#' @param directed Whether the graph is directed.
#' @param mutual Whether directed edges are mutual. It is ignored in
#'   undirected graphs.
#' @param circular Whether to create a circular ring. A non-circular
#'   ring is essentially a \dQuote{line}: a tree where every non-leaf
#'   vertex has one child.
#' @return An igraph graph.
#'
#' @family determimistic constructors
#' @export
#' @examples
#' str(make_ring(10))
#' str(make_ring(10, directed = TRUE, mutual = TRUE))

make_ring <- function(n, directed=FALSE, mutual=FALSE, circular=TRUE) {
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  res <- .Call("R_igraph_ring", as.numeric(n), as.logical(directed),
               as.logical(mutual), as.logical(circular),
               PACKAGE="igraph")
  if (igraph_opt("add.params")) {
    res$name <- "Ring graph"
    res$mutual <- mutual
    res$circular <- circular
  }
  res
}

#' @rdname make_ring
#' @param ... Passed to \code{make_ring}.
#' @export

ring <- function(...) constructor_spec(make_ring, ...)

## -----------------------------------------------------------------

#' Create tree graphs
#'
#' Create a regular tree graph.
#'
#' @aliases graph.tree
#' @concept Trees.
#' @param n Number of vertices.
#' @param children Integer scalar, the number of children of a vertex
#'   (except for leafs)
#' @param mode Defines the direction of the
#'   edges. \code{out} indicates that the edges point from the parent to
#'   the children, \code{in} indicates that they point from the children
#'   to their parents, while \code{undirected} creates an undirected
#'   graph.
#' @return An igraph graph
#'
#' @family determimistic constructors
#' @export
#' @examples
#' make_tree(10, 2)
#' make_tree(10, 3, mode = "undirected")

make_tree <- function(n, children=2, mode=c("out", "in", "undirected")) {

  mode <- igraph.match.arg(mode)
  mode1 <- switch(mode, "out"=0, "in"=1, "undirected"=2);

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  res <- .Call("R_igraph_tree", as.numeric(n), as.numeric(children),
               as.numeric(mode1),
               PACKAGE="igraph")
  if (igraph_opt("add.params")) {
    res$name <- "Tree"
    res$children <- children
    res$mode <- mode
  }
  res
}

#' @rdname make_tree
#' @param ... Passed to \code{make_tree}.
#' @export

tree <- function(...) constructor_spec(make_tree, ...)

## -----------------------------------------------------------------

#' Create a graph from the Graph Atlas
#'
#' \code{graph_from_atlas} creates graphs from the book
#' \sQuote{An Atlas of Graphs} by
#' Roland C. Read and Robin J. Wilson. The atlas contains all undirected
#' graphs with up to seven vertices, numbered from 0 up to 1252. The
#' graphs are listed:
#' \enumerate{
#'    \item in increasing order of number of nodes;
#'    \item for a fixed number of nodes, in increasing order of the number
#'      of edges;
#'    \item for fixed numbers of nodes and edges, in increasing order of
#'      the degree sequence, for example 111223 < 112222;
#'    \item for fixed degree sequence, in increasing number of
#'      automorphisms.
#' }
#'
#' @aliases graph.atlas
#' @concept Graph Atlas.
#' @param n The id of the graph to create.
#' @return An igraph graph.
#'
#' @family determimistic constructors
#' @export
#' @examples
#' ## Some randomly picked graphs from the atlas
#' graph_from_atlas(sample(0:1252, 1))
#' graph_from_atlas(sample(0:1252, 1))


graph_from_atlas <- function(n) {

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  res <- .Call("R_igraph_atlas", as.numeric(n),
               PACKAGE="igraph")
  if (igraph_opt("add.params")) {
    res$name <- sprintf("Graph from the Atlas #%i", n)
    res$n <- n
  }
  res
}

#' @rdname graph_from_atlas
#'  @param ... Passed to \code{graph_from_atlas}.
#' @export

atlas <- function(...) constructor_spec(graph_from_atlas, ...)

## -----------------------------------------------------------------

#' Create an extended chordal ring graph
#'
#' \code{make_chordal_ring} creates an extended chordal ring.
#' An extended chordal ring is regular graph, each node has the same
#' degree. It can be obtained from a simple ring by adding some extra
#' edges specified by a matrix. Let p denote the number of columns in
#' the \sQuote{\code{W}} matrix. The extra edges of vertex \code{i}
#' are added according to column \code{i mod p} in
#' \sQuote{\code{W}}. The number of extra edges is the number
#' of rows in \sQuote{\code{W}}: for each row \code{j} an edge
#' \code{i->i+w[ij]} is added if \code{i+w[ij]} is less than the number
#' of total nodes. See also Kotsis, G: Interconnection Topologies for
#' Parallel Processing Systems, PARS Mitteilungen 11, 1-6, 1993.
#'
#' @aliases graph.extended.chordal.ring
#' @param n The number of vertices.
#' @param w A matrix which specifies the extended chordal ring. See
#'   details below.
#' @return An igraph graph.
#'
#' @family determimistic constructors
#' @export
#' @examples
#' chord <- make_chordal_ring(15,
#'     matrix(c(3, 12, 4, 7, 8, 11), nr = 2))

make_chordal_ring <- function(n, w) {

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  res <- .Call("R_igraph_extended_chordal_ring", as.numeric(n),
               as.matrix(w),
               PACKAGE="igraph")
  if (igraph_opt("add.params")) {
    res$name <- "Extended chordal ring"
    res$w <- w
  }
  res
}

#' @rdname make_chordal_ring
#' @param ... Passed to \code{make_chordal_ring}.
#' @export

chordal_ring <- function(...) constructor_spec(make_chordal_ring, ...)

## -----------------------------------------------------------------

#' Line graph of a graph
#'
#' This function calculates the line graph of another graph.
#'
#' The line graph \code{L(G)} of a \code{G} undirected graph is defined as
#' follows. \code{L(G)} has one vertex for each edge in \code{G} and two
#' vertices in \code{L(G)} are connected by an edge if their corresponding
#' edges share an end point.
#'
#' The line graph \code{L(G)} of a \code{G} directed graph is slightly
#' different, \code{L(G)} has one vertex for each edge in \code{G} and two
#' vertices in \code{L(G)} are connected by a directed edge if the target of
#' the first vertex's corresponding edge is the same as the source of the
#' second vertex's corresponding edge.
#'
#' @aliases line.graph
#' @param graph The input graph, it can be directed or undirected.
#' @return A new graph object.
#' @author Gabor Csardi \email{csardi.gabor@@gmail.com}, the first version of
#' the C code was written by Vincent Matossian.
#' @keywords graphs
#' @examples
#'
#' # generate the first De-Bruijn graphs
#' g <- make_full_graph(2, directed=TRUE, loops=TRUE)
#' make_line_graph(g)
#' make_line_graph(make_line_graph(g))
#' make_line_graph(make_line_graph(make_line_graph(g)))
#'
make_line_graph <- function(graph) {

  if (!is_igraph(graph)) {
    stop("Not a graph object")
  }

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  res <- .Call("R_igraph_line_graph", graph,
               PACKAGE="igraph")
  if (igraph_opt("add.params")) {
    res$name <- "Line graph"
  }
  res
}

#' @rdname make_line_graph
#' @param ... Passed to \code{make_line_graph}.
#' @export

line_graph <- function(...) constructor_spec(make_line_graph, ...)

## -----------------------------------------------------------------

#' De Bruijn graphs
#'
#' De Bruijn graphs are labeled graphs representing the overlap of strings.
#'
#' A de Bruijn graph represents relationships between strings. An alphabet of
#' \code{m} letters are used and strings of length \code{n} are considered.  A
#' vertex corresponds to every possible string and there is a directed edge
#' from vertex \code{v} to vertex \code{w} if the string of \code{v} can be
#' transformed into the string of \code{w} by removing its first letter and
#' appending a letter to it.
#'
#' Please note that the graph will have \code{m} to the power \code{n} vertices
#' and even more edges, so probably you don't want to supply too big numbers
#' for \code{m} and \code{n}.
#'
#' De Bruijn graphs have some interesting properties, please see another
#' source, eg. Wikipedia for details.
#'
#' @aliases graph.de.bruijn
#' @param m Integer scalar, the size of the alphabet. See details below.
#' @param n Integer scalar, the length of the labels. See details below.
#' @return A graph object.
#' @author Gabor Csardi <csardi.gabor@@gmail.com>
#' @seealso \code{\link{make_kautz_graph}}, \code{\link{make_line_graph}}
#' @keywords graphs
#' @export
#' @examples
#'
#' # de Bruijn graphs can be created recursively by line graphs as well
#' g <- make_de_bruijn_graph(2,1)
#' make_de_bruijn_graph(2,2)
#' make_line_graph(g)

make_de_bruijn_graph <- function(m, n) {

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  res <- .Call("R_igraph_de_bruijn", as.numeric(m), as.numeric(n),
               PACKAGE="igraph")
  if (igraph_opt("add.params")) {
    res$name <- sprintf("De-Bruijn graph %i-%i", m, n)
    res$m <- m
    res$n <- n
  }
  res
}

#' @rdname make_de_bruijn_graph
#' @param ... Passed to \code{make_de_bruijn_graph}.
#' @export

de_bruijn_graph <- function(...) constructor_spec(make_de_bruijn_graph, ...)

## -----------------------------------------------------------------

#' Kautz graphs
#'
#' Kautz graphs are labeled graphs representing the overlap of strings.
#'
#' A Kautz graph is a labeled graph, vertices are labeled by strings of length
#' \code{n+1} above an alphabet with \code{m+1} letters, with the restriction
#' that every two consecutive letters in the string must be different. There is
#' a directed edge from a vertex \code{v} to another vertex \code{w} if it is
#' possible to transform the string of \code{v} into the string of \code{w} by
#' removing the first letter and appending a letter to it.
#'
#' Kautz graphs have some interesting properties, see eg. Wikipedia for
#' details.
#'
#' @aliases graph.kautz
#' @param m Integer scalar, the size of the alphabet. See details below.
#' @param n Integer scalar, the length of the labels. See details below.
#' @return A graph object.
#' @author Gabor Csardi <csardi.gabor@@gmail.com>, the first version in R was
#' written by Vincent Matossian.
#' @seealso \code{\link{make_de_bruijn_graph}}, \code{\link{make_line_graph}}
#' @keywords graphs
#' @export
#' @examples
#'
#' make_line_graph(make_kautz_graph(2,1))
#' make_kautz_graph(2,2)
#'
make_kautz_graph <- function(m, n) {

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  res <- .Call("R_igraph_kautz", as.numeric(m), as.numeric(n),
               PACKAGE="igraph")
  if (igraph_opt("add.params")) {
    res$name <- sprintf("Kautz graph %i-%i", m, n)
    res$m <- m
    res$n <- n
  }
  res
}

#' @rdname make_kautz_graph
#' @param ... Passed to \code{make_kautz_graph}.
#' @export

kautz_graph <- function(...) constructor_spec(make_kautz_graph, ...)

## -----------------------------------------------------------------

#' Create a full bipartite graph
#'
#' Bipartite graphs are also called two-mode by some. This function creates a
#' bipartite graph in which every possible edge is present.
#'
#' Bipartite graphs have a \sQuote{\code{type}} vertex attribute in igraph,
#' this is boolean and \code{FALSE} for the vertices of the first kind and
#' \code{TRUE} for vertices of the second kind.
#'
#' @aliases graph.full.bipartite
#' @param n1 The number of vertices of the first kind.
#' @param n2 The number of vertices of the second kind.
#' @param directed Logical scalar, whether the graphs is directed.
#' @param mode Scalar giving the kind of edges to create for directed graphs.
#' If this is \sQuote{\code{out}} then all vertices of the first kind are
#' connected to the others; \sQuote{\code{in}} specifies the opposite
#' direction; \sQuote{\code{all}} creates mutual edges. This argument is
#' ignored for undirected graphs.x
#' @return An igraph graph, with the \sQuote{\code{type}} vertex attribute set.
#' @author Gabor Csardi \email{csardi.gabor@@gmail.com}
#' @seealso \code{\link{make_full_graph}} for creating one-mode full graphs
#' @keywords graphs
#' @examples
#'
#' g <- make_full_bipartite_graph(2, 3)
#' g2 <- make_full_bipartite_graph(2, 3, dir=TRUE)
#' g3 <- make_full_bipartite_graph(2, 3, dir=TRUE, mode="in")
#' g4 <- make_full_bipartite_graph(2, 3, dir=TRUE, mode="all")
#'
make_full_bipartite_graph <- function(n1, n2, directed=FALSE,
                       mode=c("all", "out", "in")) {

  n1 <- as.integer(n1)
  n2 <- as.integer(n2)
  directed <- as.logical(directed)
  mode1 <- switch(igraph.match.arg(mode), "out"=1, "in"=2, "all"=3, "total"=3)

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  res <- .Call("R_igraph_full_bipartite", n1, n2, as.logical(directed), mode1,
               PACKAGE="igraph")
  if (igraph_opt("add.params")) {
    res$graph$name <- "Full bipartite graph"
    res$n1 <- n1
    res$n2 <- n2
    res$mode <- mode
  }
  set_vertex_attr(res$graph, "type", value=res$types)
}

#' @rdname make_full_bipartite_graph
#' @param ... Passed to \code{make_full_bipartite_graph}.
#' @export

full_bipartite_graph <- function(...) constructor_spec(make_full_bipartite_graph, ...)

## -----------------------------------------------------------------

#' Create a bipartite graph
#'
#' A bipartite graph has two kinds of vertices and connections are only allowed
#' between different kinds.
#'
#' Bipartite graphs have a \code{type} vertex attribute in igraph, this is
#' boolean and \code{FALSE} for the vertices of the first kind and \code{TRUE}
#' for vertices of the second kind.
#'
#' \code{make_bipartite_graph} basically does three things. First it checks tha
#' \code{edges} vector against the vertex \code{types}. Then it creates a graph
#' using the \code{edges} vector and finally it adds the \code{types} vector as
#' a vertex attribute called \code{type}.
#'
#' \code{is_bipartite} checks whether the graph is bipartite or not. It just
#' checks whether the graph has a vertex attribute called \code{type}.
#'
#' @aliases make_bipartite_graph graph.bipartite is.bipartite is_bipartite
#' @param types A vector giving the vertex types. It will be coerced into
#' boolean. The length of the vector gives the number of vertices in the graph.
#' @param edges A vector giving the edges of the graph, the same way as for the
#' regular \code{\link{graph}} function. It is checked that the edges indeed
#' connect vertices of different kind, accoding to the supplied \code{types}
#' vector.
#' @param directed Whether to create a directed graph, boolean constant. Note
#' that by default undirected graphs are created, as this is more common for
#' bipartite graphs.
#' @param graph The input graph.
#' @return \code{make_bipartite_graph} returns a bipartite igraph graph. In other
#' words, an igraph graph that has a vertex attribute named \code{type}.
#'
#' \code{is_bipartite} returns a logical scalar.
#' @author Gabor Csardi \email{csardi.gabor@@gmail.com}
#' @seealso \code{\link{graph}} to create one-mode networks
#' @keywords graphs
#' @examples
#'
#' g <- make_bipartite_graph( rep(0:1,length=10), c(1:10))
#' print(g, v=TRUE)
#'
make_bipartite_graph <- function(types, edges, directed=FALSE) {

  types <- as.logical(types)
  edges <- as.numeric(edges)-1
  directed <- as.logical(directed)

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  res <- .Call("R_igraph_create_bipartite", types, edges, directed,
               PACKAGE="igraph")
  set_vertex_attr(res, "type", value=types)
}

#' @rdname make_bipartite_graph
#' @param ... Passed to \code{make_bipartite_graph}.
#' @export

bipartite_graph <- function(...) constructor_spec(make_bipartite_graph, ...)

## -----------------------------------------------------------------

#' Create a complete (full) citation graph
#'
#' \code{make_full_citation_graph} creates a full citation graph. This is a
#' directed graph, where every \code{i->j} edge is present if and only if
#' \eqn{j<i}. If \code{directed=FALSE} then the graph is just a full graph.
#'
#' @aliases graph.full.citation
#' @param n The number of vertices.
#' @param directed Whether to create a directed graph.
#' @return An igraph graph.
#'
#' @family determimistic constructors
#' @export
#' @examples
#' str(make_full_citation_graph(10))

make_full_citation_graph <- function(n, directed=TRUE) {
  # Argument checks
  n <- as.integer(n)
  directed <- as.logical(directed)

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  # Function call
  res <- .Call("R_igraph_full_citation", n, directed,
        PACKAGE="igraph")

  res <- set_graph_attr(res, 'name', 'Full citation graph')
  res
}

#' @rdname make_full_citation_graph
#' @param ... Passed to \code{make_full_citation_graph}.
#' @export

full_citation_graph <- function(...) constructor_spec(make_full_citation_graph, ...)

## -----------------------------------------------------------------
