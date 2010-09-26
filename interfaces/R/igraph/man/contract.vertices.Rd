\name{contract.vertices}
\alias{contract.vertices}
\title{Contract several vertices into a single one}
\description{
  This function creates a new graph, by merging several 
  vertices into one. The vertices in the new graph correspond
  to sets of vertices in the input graph.
}
\usage{
contract.vertices(graph, mapping,
                  vertex.attr.comb=igraph.par("vertex.attr.comb"))
}
\arguments{
  \item{graph}{The input graph, it can be directed or undirected.}
  \item{mapping}{A numeric vector that specifies the mapping. Its
    elements correspond to the vertices, and for each element the id in
    the new graph is given.}
  \item{vertex.attr.comb}{Specifies how to combine the vertex attributes
    in the new graph. Please see \code{\link{attribute.combination}} for
    details.
  }
}
\details{
  The attributes of the graph are kept. Graph and edge attributes are
  unchanged, vertex attributes are combined, according to the
  \code{vertex.attr.comb} parameter.
}
\value{
  A new graph object.
}
% \references{}
\author{Gabor Csardi \email{csardi.gabor@gmail.com}}
% \seealso{}
\examples{
g <- graph.ring(10)
g$name <- "Ring"
V(g)$name <- letters[1:vcount(g)]
E(g)$weight <- runif(ecount(g))

g2 <- contract.vertices(g, rep(1:5, each=2),
                        vertex.attr.comb=toString)

## graph and edge attributes are kept, vertex attributes are
## combined using the 'toString' function.
print(g2, g=TRUE, v=TRUE, e=TRUE)

}
\keyword{graphs}

  
  