\name{is.named}
\alias{is.named}
\concept{Named graphs}
\title{Named graphs}
\description{An igraph graph is named, if there is a symbolic name
  associated with its vertices.}
\usage{
is.named(graph)
}
\arguments{
  \item{graph}{The input graph.}
}
\details{
  In igraph vertices can always be identified and specified via their
  numeric vertex ids. This is, however, not always convenient, and in
  many cases there exist symbolic ids that correspond to the
  vertices. To allow this more flexible identification of vertices, one
  can assign a vertex attribute called \sQuote{name} to an igraph
  graph. After doing this, the symbolic vertex names can be used in all
  igraph functions, instead of the numeric ids.

  Note that the uniqueness of vertex names are currently not enforced in
  igraph, you have to check that for yourself, when assigning the vertex
  names.
}
\value{A logical scalar.}
\author{Gabor Csardi \email{csardi.gabor@gmail.com}}
\examples{
g <- graph.ring(10)
is.named(g)
V(g)$name <- letters[1:10]
is.named(g)
neighbors(g, "a")
}
\keyword{graphs}
