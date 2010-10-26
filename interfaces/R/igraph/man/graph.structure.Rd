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
\alias{vertex}
\alias{path}
\title{Method for structural manipulation of graphs}
\description{These are the methods for simple manipulation of graphs:
  adding and deleting edges and vertices.}
\usage{
\method{[}{igraph}(x, i, j, \dots, sparse=getIgraphOpt("sparsematrices"),
                   edges=FALSE, drop=TRUE,
                   attr=if (is.weighted(x)) "weight" else NULL)
\method{[[}{igraph}(x, i, j, \dots, directed=TRUE, edges=FALSE, exact=TRUE)
\method{[}{igraph}(x, i, j, \dots,
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
  much like the graph was an adjacency matrix. The single bracket
  indexes the (possibly weighted) adjacency matrix of the graph. The
  double bracket operator is similar, but queries the adjacencly list of
  the graph. The details on how to use the indexing operators are
  discussed below.

  The addition (\sQuote{\code{+}}) and division (\sQuote{\code{-}})
  operators can also be used to add and remove vertices and edges. This
  form is sometimes more readable, and is usually the best if the user 
  also wants to add attributes, together with the new vertices/edges.
  Please see the details below.

  In addition, the four functions \code{add.vertices}, \code{add.edges},
  \code{delete.vertices} and \code{delete.edges} can also be used
  to manipulate the structure.
}
\section{The indexing operators}{
  The one-bracket (\sQuote{\code{[}}) and two-brackets
  (\sQuote{\code{[[}}) indexing operators allow relatively
  straightforward query and update operations on graphs. Here is
  a little summary of what they can do.

  TODO

  Of course, the indexing operators support vertex names,
  so instead of a numeric vertex id a vertex can also be given to
  \sQuote{\code{[}} and \sQuote{\code{[[}}. The same is true for
  edge names as well.
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
\author{Gabor Csardi \email{csardi@rmki.kfki.hu}}
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
