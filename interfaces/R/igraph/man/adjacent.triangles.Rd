\name{triangles}
\alias{count_triangles}
\alias{adjacent.triangles}
\alias{triangles}
\title{Find triangles in graphs}
\description{Count how many triangles a vertex is part of, in a graph,
  or just list the triangles of a graph.}
\usage{
triangles (graph)
count_triangles (graph, vids = V(graph))
}
\arguments{
  \item{graph}{The input graph. It might be directed, but edge
    directions are ignored.}
  \item{vids}{The vertices to query, all of them by default. This might
    be a vector of numeric ids, or a character vector of symbolic vertex
    names for named graphs.}
}
\details{
  \code{triangles} lists all triangles of a graph. For efficiency, all
  triangles are returned in a single vector. The first three vertices
  belong to the first triangle, etc.
  
  \code{count_triangles} counts how many triangles a vertex is part
  of.
}
\value{
  For \code{triangles} a numeric vector of vertex ids, the first three
  vertices belong to the first triangle found, etc.
  
  For \code{triangles} a numeric vector, the number of triangles for all
  vertices queried.
}
% \references{}
\author{ Gabor Csardi \email{csardi.gabor@gmail.com} }
\seealso{\code{\link{transitivity}}}
\examples{
## A small graph
kite <- make_graph("Krackhardt_Kite")
plot(kite)
matrix(triangles(kite), nrow=3)

## Adjacenct triangles
atri <- count_triangles(kite)
plot(kite, vertex.label=atri)

## Always true
sum(count_triangles(kite)) == length(triangles(kite))

## Should match, local transitivity is the
## number of adjacent triangles divided by the number
## of adjacency triples
transitivity(kite, type="local")
count_triangles(kite) / (degree(kite) * (degree(kite)-1)/2)
}
\keyword{graphs}
