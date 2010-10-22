\name{graph.structure}
\alias{add.edges}
\alias{add.vertices}
\alias{delete.edges}
\alias{delete.vertices}
\alias{[.igraph}
\alias{[[.igraph}
\alias{[<-.igraph}
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

add.edges(graph, edges, \dots, attr=list())
add.vertices(graph, nv, \dots, attr=list())
delete.edges(graph, edges)
delete.vertices(graph, v)
}
\arguments{
  \item{x,graph}{The graph to work on.}
  \item{i,j}{Vertex ids or names. See details below.}
  \item{\dots}{These are currently ignored for the indexing operators.
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
  \item{attr}{Additional edge/vertex attributes to add. This will be
    concatenated to the other supplied attributes.}
  \item{nv}{Numeric constant, the number of vertices to add.}
  \item{v}{Vector sequence, the vertices to remove.}
}
\details{
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
      For non-existant edges the vector contains \code{NA} values.
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
g <- graph.ring(10)
add.edges(g, c(2,6,3,7) )
delete.edges(g, E(g, P=c(1,10, 2,3)) )
delete.vertices(g, c(2,7,8) )
}
\keyword{graphs}
