\name{graph.difference}
\alias{graph.difference}
\alias{\%m\%}
\concept{Graph operators}
\title{Difference of graphs}
\description{The difference of two graphs are created.}
\usage{
graph.difference(big, small, byname = "auto")
}
\arguments{
  \item{big}{The left hand side argument of the minus operator. A
    directed or undirected graph.}
  \item{small}{The right hand side argument of the minus operator. A
    directed ot undirected graph.}
  \item{byname}{A logical scalar, or the character scalar
    \code{auto}. Whether to perform the operation based on symbolic
    vertex names. If it is \code{auto}, that means \code{TRUE} if both
    graphs are named and \code{FALSE} otherwise. A warning is generated
    if \code{auto} and one graph, but not both graphs are named.}
}
\details{
  \code{graph.difference} creates the difference of two graphs. Only
  edges present in the first graph but not in the second will be be
  included in the new graph. The corresponding operator is \%m\%.

  If the \code{byname} argument is \code{TRUE} (or \code{auto} and the
  graphs are all named), then the operation is performed based on
  symbolic vertex names. Otherwise numeric vertex ids are used.
  
  \code{graph.difference} keeps all attributes
  (graph, vertex and edge) of the first graph.

  Note that \code{big} and \code{small} must both be directed or both be
  undirected, otherwise an error message is given.
}
\value{
  A new graph object.
}
% \references{}
\author{Gabor Csardi \email{csardi.gabor@gmail.com}}
% \seealso{}
\examples{
## Create a wheel graph
wheel <- graph.union(graph.ring(10),
                     graph.star(11, center=11, mode="undirected"))
V(wheel)$name <- letters[seq_len(vcount(wheel))]

## Subtract a star graph from it
sstar <- graph.star(6, center=6, mode="undirected")
V(sstar)$name <- letters[c(1,3,5,7,9,11)]
G <- wheel \%m\% sstar
str(G)
plot(G, layout=layout.auto(wheel))
}
\keyword{graphs}
