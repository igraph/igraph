\name{adjacent.triangles}
\alias{adjacent.triangles}
\title{Count adjacenct triangles}
\description{Count how many triangles a vertex is part of, in a graph.}
\usage{
adjacent.triangles (graph, vids = V(graph))
}
\arguments{
  \item{graph}{The input graph. It might be directed, but edge
    directions are ignored.}
  \item{vids}{The vertices to query, all of them by default. This might
    be a vector of numeric ids, or a character vector of symbolic vertex
    names for named graphs.}
}
\details{
  Count how many triangles a vertex is part of.
}
\value{
  A numeric vector, the number of triangles for all vertices queried.
}
% \references{}
\author{ Gabor Csardi \email{csardi.gabor@gmail.com} }
\seealso{\code{\link{transitivity}}}
\examples{
## A small graph
kite <- graph.famous("Krackhardt_Kite")
atri <- adjacent.triangles(kite)
plot(kite, vertex.label=atri)

## Should match, local transitivity is the
## number of adjacent triangles divided by the number
## of adjacency triples
transitivity(kite, type="local")
adjacent.triangles(kite) / (degree(kite) * (degree(kite)-1)/2)
}
\keyword{graphs}
