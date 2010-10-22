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
\method{[}{igraph}(x, i, j, \dots, directed=TRUE, weighted=TRUE,
            multi=FALSE, edges=NULL, simplify=TRUE, unique=TRUE,
            drop=TRUE)
\method{[[}{igraph}(x, i, j, \dots, directed=TRUE, multi=FALSE,
            simplify=TRUE, unique=TRUE, exact=TRUE)
\method{[}{igraph}(x, i, j, \dots, add=FALSE) <- value 

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
  \item{i,j}{Vertex ids or names. See details below.}
  \item{\dots}{These are currently ignored for the indexing operators.
    For \code{vertex}, \code{vertices}, \code{edge}, \code{edges} and
    \code{path} see details below.
    For \code{add.edges} and \code{add.vertices} these additional
    parameters will be added as edge/vertex attributes. Note that these
    arguments have to be named.} 
  \item{directed}{Logical scalar, whether to consider edge directions
    in directed graphs. It is ignored for undirected graphs.}
  \item{weighted}{Logical scalar, whether to use the
    \sQuote{\code{weight}} edge attribute as edge weights.}
  \item{multi}{Logical scalar, whether to handle multiple edges. If it
     is true and the same pairs of vertices are given multiple times in
     \sQuote{\code{i}} and \sQuote{\code{j}}, then multiple edges are
     properly reported (if they exist). Otherwise an arbitrary edge is
     reported between the pairs.}
  \item{edges}{Edge sequence. For the indexing operator it gives the
    edges whose incident vertices are returned. For \code{delete.edges}
    these are the edges to remove.}
  \item{simplify}{Logical scalar. Whether to put together
    vertices/edges into a single vector, when querying adjacenct
    vertices or incident edges.}
  \item{unique}{Logical scalar. It is only used if
    \code{simplify=TRUE}. Whether to make the result vector unique.}
  \item{drop,exact}{These arguments are ignored.}
  \item{add}{Logical scalar, whether to add a new edge to the graph,
     even if an edge already exists between the specified vertices.}
  \item{value}{A logical or numeric vector or scalar. If \code{FALSE},
     then the specified edges will be deleted. If \code{TRUE}, then the
     specified edges will be added. (By default only if they don't yet
     exist, see also the \code{add} argument.) If numeric, then the
     given values will be added as edge weights, together with the
     edges if they don't yet exist.}
  \item{e2}{See details below.}
  \item{attr}{Additional edge/vertex attributes to add. This will be
    concatenated to the other supplied attributes.}
  \item{nv}{Numeric constant, the number of vertices to add.}
  \item{v}{Vector sequence, the vertices to remove.}
}
\details{
  There are, by and large, three ways to manipulate the structure
  of a graph in igraph. The first way is using the \sQuote{\code{[}} and
  \sQuote{\code{[[}} indexing operators on the graph object, very
  much like the graph was an adjacency matrix. The single bracket
  indexes the (possibly weighted) adjacency matrix of the graph. The
  double bracket operator is similar, but indexes a matrix that contains
  the ids of the edges, for connected pairs of vertices. (And \code{NA}
  otherwise.) The details on how to use the indexing operators are
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

  \enumerate{
    \item Check whether there is an edge between two vertices (\eqn{v}
      and \eqn{w}) in the graph: \preformatted{  graph[v, w]}
      A logical scalar is returned.

      The \sQuote{\code{directed}} option can be used to ignore edge
      direction in directed graphs, i.e.
      \preformatted{  graph[v, w, directed=FALSE]}
      returns \code{TRUE}
      if an edge from \eqn{v} to \eqn{w}, or an edge from \eqn{w} to
      \eqn{v} exists.
    \item The same syntax can be used to check multiple vertex pairs,
      in other words, \eqn{v} and \eqn{w} may be vectors, and in this
      case a logical vector is returned.
    \item Weighted graphs behave somewhat differently. Instead of a
      logical vector, a numeric vector of edge weights is returned.
      For non-existant edges the vector contains zero weights.
    \item Query adjacent vertices of one or more vertices. To do
      this, you need to give only one of the \sQuote{\code{i}} and
      \sQuote{\code{j}} arguments: \preformatted{  graph[v, ]
  graph[, w]}
      Here the first line list all vertices reachable from
      \eqn{v} via a directed edge, and the second line lists vertices
      from which \eqn{w} is reachable by a directed edge.
    \item The same syntax works for multiple vertices as well. By
      default all adjacency vertices are put into a single vector,
      but the \sQuote{\code{simplify}} and \sQuote{\code{unique}}
      arguments change this behavior. If \sQuote{\code{unique}} is
      \code{FALSE}, then a non-unique vector is returned. If
      \sQuote{\code{simplify}} is \code{FALSE}, then a list of vectors
      is returned, one vector for each vertex.
    \item Incident vertices for an edge (or more edges). For this you
      need to give the \sQuote{\code{edges}}
      argument: \preformatted{  graph[edges=e]}
      Whether a single vector is returned, or a list of vectors, is
      controlled by the \sQuote{\code{simplify}} and
      \sQuote{\code{unique}} arguments, just like in the previous
      case.
    \item The double-brackets indexing operator can be used to
      query edge ids. E.g. \preformatted{  graph[[v, w]]}
      returns the id of the edge from \eqn{v} to \eqn{w}.
    \item To query the ids of incident edges, give only one of the
      arguments: \preformatted{  graph[[v, ]]
  graph[[, w]]}
      The first form queries outgoing, the second incoming edges in
      directed graphs. As before, the \sQuote{\code{simplify}} and
      \sQuote{\code{unique}} arguments define whether to return a
      list of edge id vectors, separately for each vertex, or a
      (possibly unique) vector of all of them.
    \item Adding edges. The single bracket (\sQuote{\code{[}}) operator
      can be used for adding edges to the graph: \preformatted{  graph[v, w] <-TRUE}
      adds a new (possibly directed) edge from vertex \eqn{v} to vertex
      \eqn{w}. By default the new edge is only added, if it does not yet
      exist an edge between the two vertices. If the \sQuote{\code{add}}
      argument is set to \code{TRUE}, then a new edge is always added.
    \item Adding several edges to the graph. The syntax in the previous
      point can be used to add several edge to the graph, \code{v} and
      \code{w} can be vectors of the same length, and then their
      vertices will be pairwise connected.
    \item Adding weighted edges. For weighted graphs, weighted edges
      can be added to the graph, if the replacement value is not
      logical, but a numeric vector: \preformatted{  graph[v, w] <- 1/2}
      This adds a new edge from vertex \eqn{v} to vertex \eqn{w}, and
      sets its weight to \eqn{1/2}. As previously, \code{v} and \code{w}
      can be vectors of vertex ids, to add many edges at the same time.
      In this case the replacement value can be a vector of the same
      length as \code{v} and \code{w} to set different edge weights for
      the new edges; or it can be a single number, to set the same
      weight for all new edges.
    \item Updating edge weights. The same syntax can be used to update
      the edge weights of already existing edges. Even mixing the
      addition of edges and the updating of weights is possibly, if some
      of the specified edges already exist, but some don't.
    \item Deleting edges. The replacement syntax allow the deletion of
      edges, by specifying \code{FALSE} or \code{NULL} as the
      replacement value: \preformatted{  graph[v, w] <- FALSE}
      removes the edge from vertex \eqn{v} to vertex \eqn{w}.
      As befre, the \code{v} and \code{w} arguments can be vectors of
      vertex ids, to remove two or more edges at the same time.
  }

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
