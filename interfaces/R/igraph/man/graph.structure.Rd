\name{graph.structure}
\alias{add.edges}
\alias{add.vertices}
\alias{delete.edges}
\alias{delete.vertices}
\alias{[.igraph}
\alias{[[.igraph}
\alias{[<-.igraph}
\alias{+.igraph}
\alias{-.igraph}
\alias{edge}
\alias{edges}
\alias{vertex}
\alias{vertices}
\alias{path}
\title{Method for structural manipulation of graphs}
\description{These are the methods for simple manipulation of graphs:
  adding and deleting edges and vertices.}
\usage{
\method{[}{igraph}(x, i, j, \dots, from, to,
                   sparse=getIgraphOpt("sparsematrices"),
                   edges=FALSE, drop=TRUE,
                   attr=if (is.weighted(x)) "weight" else NULL)
\method{[[}{igraph}(x, i, j, \dots, directed=TRUE, edges=FALSE, exact=TRUE)
\method{[}{igraph}(x, i, j, \dots, from, to,
                   attr=if (is.weighted(x)) "weight" else NULL) <- value 

\method{+}{igraph}(e1, e2)
\method{-}{igraph}(e1, e2)
vertex(\dots)
vertices(\dots)
edge(\dots)
edges(\dots)
path(\dots)

add.edges(graph, edges, \dots, attr=list())
add.vertices(graph, nv, \dots, attr=list())
delete.edges(graph, edges)
delete.vertices(graph, v)
}
\arguments{
  \item{x,graph,e1}{The graph to work on.}
  \item{i,j}{Vertex ids or names or logical vectors. See details below.}
  \item{\dots}{These are currently ignored for the indexing operators.
    For \code{vertex}, \code{vertices}, \code{edge}, \code{edges} and
    \code{path} see details below.
    For \code{add.edges} and \code{add.vertices} these additional
    parameters will be added as edge/vertex attributes. Note that these
    arguments have to be named.}
  \item{from}{A numeric or character vector giving vertex ids or names.
    Together with the \code{to} argument, it can be used to query/set
    a sequence of edges. See details below.

    This argument cannot be present together with any of the
    \code{i} and \code{j} arguments and if it is present, then the
    \code{to} argument must be present as well.}
  \item{to}{A numeric or character vector giving vertex ids or names.
    Together with the \code{from} argument, it can be used to query/set
    a sequence of edges. See details below.

    This argument cannot be present together with any of the
    \code{i} and \code{j} arguments and if it is present, then the
    \code{from} argument must be present as well.}
  \item{sparse}{Logical scalar, whether to use sparse matrix.}
  \item{directed}{Logical scalar, whether to consider edge directions
    in directed graphs. It is ignored for undirected graphs.}
  \item{edges}{Logical scalar, whether to return edge ids.}
  \item{drop,exact}{These arguments are ignored.}
  \item{value}{A logical or numeric scalar or \code{NULL}. If
     \code{FALSE}, \code{NULL} or zero, then the specified edges will be
     deleted. If \code{TRUE} or a non-zero numeric value, then the
     specified edges will be added. (Only if they don't yet exist.)}
  \item{e2}{See details below.}
  \item{attr}{For the indexing operators: if not \code{NULL}, then it
     should be the name of an edge attribute. This attribute is queried,
     or updated to the given value. For \code{add.edges} and
     \code{add.vertices}: additional edge/vertex attributes to add. This
     will be concatenated to the other supplied attributes.}
  \item{nv}{Numeric constant, the number of vertices to add.}
  \item{v}{Vector sequence, the vertices to remove.}
}
\details{
  There are, by and large, three ways to manipulate the structure
  of a graph in igraph. The first way is using the \sQuote{\code{[}} and
  \sQuote{\code{[[}} indexing operators on the graph object, very
  much like the graph was an adjacency matrix (\code{[}) or an adjacency
  list (\code{[}). The single bracket indexes the (possibly weighted)
  adjacency matrix of the graph. The double bracket operator is
  similar, but queries the adjacencly list of the graph. The details
  on how to use the indexing operators are 
  discussed below.

  The addition (\sQuote{\code{+}}) and division (\sQuote{\code{-}})
  operators can also be used to add and remove vertices and edges. This
  form is sometimes more readable, and is usually the best if the user 
  also wants to add attributes, together with the new vertices/edges.
  Please see the details below.

  In addition, the four functions, \code{add.vertices}, \code{add.edges},
  \code{delete.vertices} and \code{delete.edges} can also be used
  to manipulate the structure.
}
\section{The indexing operators}{
  The one-bracket (\sQuote{\code{[}}) and two-brackets
  (\sQuote{\code{[[}}) indexing operators allow relatively
  straightforward query and update operations on graphs. The one bracket
  operator works on the (imaginary) adjacency matrix of the graph.
  Here is what you can do with it:

  \enumerate{
    \item Check whether there is an edge between two vertices (\eqn{v}
      and \eqn{w}) in the graph: \preformatted{  graph[v, w]}
      A numeric scalar is returned, one if the edge exists, zero
      otherwise.
    \item Extract the (sparse) adjacency matrix of the graph, or part of
      it: \preformatted{  graph[]
  graph[1:3,5:6]
  graph[c(1,3,5),]}
      The first variants returns the full adjacency matrix, the other
      two return part of it.
    \item The \code{from} and \code{to} arguments can be used to check
      the existence of many edges. In this case, both \code{from} and
      \code{to} must be present and they must have the same length. They
      must contain vertex ids or names. A numeric vector is returned, of
      the same length as \code{from} and \code{to}, it contains ones
      for existing edges edges and zeros for non-existing ones.
      Example: \preformatted{  graph[from=1:3, to=c(2,3,5)]}.
    \item For weighted graphs, the \code{[} operator returns the edge
      weights. For non-esistent edges zero weights are returned. Other
      edge attributes can be queried as well, by giving the \code{attr}
      argument. 
    \item Querying edge ids instead of the existance of edges or edge
      attributes. E.g. \preformatted{  graph[1, 2, edges=TRUE]}
      returns the id of the edge between vertices 1 and 2, or zero if
      there is no such edge.
    \item Adding one or more edges to a graph. For this the element(s) of
      the imaginary adjacency matrix must be set to a non-zero numeric
      value (or \code{TRUE}): \preformatted{  graph[1, 2] <- 1
  graph[1:3,1] <- 1
  graph[from=1:3, to=c(2,3,5)] <- TRUE}
      This does not affect edges that are already present in the graph,
      i.e. no multiple edges are created.
    \item Adding weighted edges to a graph. The \code{attr} argument
      contains the name of the edge attribute to set, so it does not
      have to be \sQuote{weight}: \preformatted{  graph[1, 2, attr="weight"]<- 5
  graph[from=1:3, to=c(2,3,5)] <- c(1,-1,4)}
      If an edge is already present in the network, then only its
      weigths or other attribute are updated. If the graph is already
      weighted, then the \code{attr="weight"} setting is implicit, and
      one does not need to give it explicitly.
    \item Deleting edges. The replacement syntax allow the deletion of
      edges, by specifying \code{FALSE} or \code{NULL} as the
      replacement value: \preformatted{  graph[v, w] <- FALSE}
      removes the edge from vertex \eqn{v} to vertex \eqn{w}.
      As this can be used to delete edges between two sets of vertices,
      either pairwise: \preformatted{  graph[from=v, to=w] <- FALSE}
      or not: \preformatted{  graph[v, w] <- FALSE }
      if \eqn{v} and \eqn{w} are vectors of edge ids or names.
  }

  The double bracket operator indexes the (imaginary) adjacency list
  of the graph. This can used for the following operations:
  \enumerate{
    \item Querying the adjacent vertices for one or more
      vertices: \preformatted{  graph[[1:3,]]
  graph[[,1:3]]}
      The first form gives the successors, the second the predessors
      or the 1:3 vertices. (For undirected graphs they are equivalent.)
    \item Querying the incident edges for one or more vertices,
      if the \code{edges} argument is set to
      \code{TRUE}: \preformatted{  graph[[1:3, , edges=TRUE]]
  graph[[, 1:3, edges=TRUE]]}
    \item Querying the edge ids between two sets or vertices,
      if both indices are used. E.g. \preformatted{  graph[[v, w, edges=TRUE]]}
      gives the edge ids of all the edges that exist from vertices
      \eqn{v} to vertices \eqn{w}.
  }

  Both the \sQuote{\code{[}} and \sQuote{\code{[[}} operators allow
  logical indices and negative indices as well, with the usual R
  semantics. E.g. \preformatted{  graph[degree(graph)==0, 1] <- 1}
  adds an edge from every isolate vertex to vertex one,
  and \preformatted{  G <- graph.empty(10)
  G[-1,1] <- TRUE}
  creates a star graph.

  Of course, the indexing operators support vertex names,
  so instead of a numeric vertex id a vertex can also be given to
  \sQuote{\code{[}} and \sQuote{\code{[[}}.
}
\section{The plus operator for adding vertices and edges}{
  The plus operator can be used to add vertices or edges to graph.
  The actual operation that is performed depends on the type of the
  right hand side argument.
  \itemize{
  \item If it is another igraph graph object, then the disjoint union of
    the two graphs is calculated, see \code{\link{graph.disjoint.union}}.
  \item If it is a numeric scalar, then the specified number of vertices
    are added to the graph.
  \item If it is a character scalar or vector, then it is interpreted as
    the names of the vertices to add to the graph.
  \item If it is an object created with the \code{vertex} or
    \code{vertices} function, then new vertices are added to the
    graph. This form is appropriate when one wants to add some vertex
    attributes as well. The operands of the \code{vertices} function
    specifies the number of vertices to add and their attributes as
    well.

    The unnamed arguments of \code{vertices} are concatenated and
    used as the \sQuote{\code{name}} vertex attribute (i.e. vertex
    names), the named arguments will be added as additional vertex
    attributes. Examples: \preformatted{  g <- g + vertex(shape="circle", color="red")
  g <- g + vertex("foo", color="blue")
  g <- g + vertex("bar", "foobar")
  g <- g + vertices("bar2", "foobar2", color=1:2, shape="rectangle")}
    See more examples below.

    \code{vertex} is just an alias to \code{vertices}, and it is
    provided for readability. The user should use it if a single vertex
    is added to the graph.

  \item If it is an object created with the \code{edge} or \code{edges}
    function, then new edges will be added to the graph. The new edges
    and possibly their attributes can be specified as the arguments of
    the \code{edges} function.

    The unnamed arguments of \code{edges} are concatenated and used
    as vertex ids of the end points of the new edges. The named
    arguments will be added as edge attributes.

    Examples: \preformatted{  g <- graph.empty() + vertices(letters[1:10]) +
  vertices("foo", "bar", "bar2", "foobar2")
  g <- g + edge("a", "b")
  g <- g + edges("foo", "bar", "bar2", "foobar2")
  g <- g + edges(c("bar", "foo", "foobar2", "bar2"), color="red", weight=1:2)}
    See more examples below.

    \code{edge} is just an alias to \code{edges} and it is provided
    for readability. The user should use it if a single edge is added to
    the graph.

  \item If it is an object created with the \code{path} function, then
    new edges that form a path are added. The edges and possibly their
    attributes are specified as the arguments to the \code{path}
    function. The non-named arguments are concatenated and interpreted
    as the vertex ids along the path. The remaining arguments are added
    as edge attributes.

    Examples: \preformatted{  g <- graph.empty() + vertices(letters[1:10])
  g <- g + path("a", "b", "c", "d")
  g <- g + path("e", "f", "g", weight=1:2, color="red")
  g <- g + path(c("f", "c", "j", "d"), width=1:3, color="green")}
  }

  It is important to note that, although the plus operator is
  commutative, i.e. is possible to write \preformatted{  graph <- "foo" + graph.empty()}
  it is not associative, e.g. \preformatted{  graph <- "foo" + "bar" + graph.empty()}
  results a syntax error, unless parentheses are used: \preformatted{  graph <- "foo" + ( "bar" + graph.empty() )}
  For clarity, we suggest to always put the graph object on the left
  hand side of the operator: \preformatted{  graph <- graph.empty() + "foo" + "bar"}  
}
\section{The minus operator for deleting vertices and edges}{
  The minus operator (\sQuote{\code{-}}) can be used to remove vertices
  or edges from the graph. The operation performed is selected based on
  the type of the right hand side argument:
  \itemize{
  \item If it is an igraph graph object, then the difference of the
    two graphs is calculated, see \code{\link{graph.difference}}.
  \item If it is a numeric or character vector, then it is interpreted
    as a vector of vertex ids and the specified vertices will be
    deleted from the graph. Example: \preformatted{  g <- graph.ring(10)
  V(g)$name <- letters[1:10]
  g <- g - c("a", "b")}
  \item If \code{e2} is a vertex sequence (e.g. created by the
    \code{\link{V}} function), then these vertices will be deleted from
    the graph.
  \item If it is an edge sequence (e.g. created by the \code{\link{E}}
    function), then these edges will be deleted from the graph.
  \item If it is an object created with the \code{vertex} (or the
    \code{vertices}) function, then all arguments of \code{vertices} are
    concatenated and the result is interpreted as a vector of vertex
    ids. These vertices will be removed from the graph.
  \item If it is an object created with the \code{edge} (or the
    \code{edges}) function, then all arguments of \code{edges} are
    concatenated and then interpreted as edges to be removed from the
    graph.
    Example: \preformatted{  g <- graph.ring(10)
  V(g)$name <- letters[1:10]
  E(g)$name <- LETTERS[1:10]
  g <- g - edge("e|f")
  g <- g - edge("H")}
  \item If it is an object created with the \code{path} function,
    then all \code{path} arguments are concatenated and then interpreted
    as a path along which edges will be removed from the graph.
    Example: \preformatted{  g <- graph.ring(10)
  V(g)$name <- letters[1:10]
  g <- g - path("a", "b", "c", "d")}
  }
}
\section{More functions to manipulate graph structure}{
  \code{add.edges} adds the specified edges to the graph. The ids of the
  vertices are preserved. The additionally supplied named arguments will
  be added as edge attributes for the new edges. If an attribute was not
  present in the original graph, its value for the original edges will
  be \code{NA}.

  \code{add.vertices} adds the specified number of isolate vertices to
  the graph. The ids of the old vertices are preserved. The additionally
  supplied named arguments will be added as vertex attributes for the
  new vertices. If an attribute was not present in the original graph,
  its value is set to \code{NA} for the original vertices.

  \code{delete.edges} removes the specified edges from the graph. If a
  specified edge is not present, the function gives an error message,
  and the original graph remains unchanged.
  The ids of the vertices are preserved.

  \code{delete.vertices} removes the specified vertices from the graph
  together with their adjacent edges. The ids of the vertices are
  \emph{not} preserved.
}
\value{For the indexing operators see the description above. The other
  functions return a new graph.}
% \references{}
\author{Gabor Csardi \email{csardi.gabor@gmail.com}}
% \seealso{}
\examples{
# 10 vertices named a,b,c,... and no edges
g <- graph.empty() + vertices(letters[1:10])

# Add edges to make it a ring
g <- g + path(letters[1:10], letters[1], color="grey")

# Add some extra random edges
g <- g + edges(sample(V(g), 10, replace=TRUE), color="red")
g$layout <- layout.circle
if (interactive()) {
  plot(g)
}

# The old-style operations
g <- graph.ring(10)
add.edges(g, c(2,6,3,7) )
delete.edges(g, E(g, P=c(1,10, 2,3)) )
delete.vertices(g, c(2,7,8) )
}
\keyword{graphs}
