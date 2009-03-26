\name{cliques}
\alias{cliques}
\alias{largest.cliques}
\alias{maximal.cliques}
\alias{clique.number}
\concept{Clique}
\concept{Maximal clique}
\concept{Largest clique}
\title{The functions find cliques, ie. complete subgraphs in a graph}
\description{These functions find all, the largest or all the maximal
  cliques in an undirected graph. The size of the largest clique can
  also be calculated.}
\usage{
cliques(graph, min=NULL, max=NULL)
largest.cliques(graph)
maximal.cliques(graph)
clique.number(graph)
}
\arguments{
  \item{graph}{The input graph, directed graphs will be considered as
    undirected ones, multiple edges and loops are ignored.}
  \item{min}{Numeric constant, lower limit on the size of the cliques to
    find. \code{NULL} means no limit, ie. it is the same as 0.}
  \item{max}{Numeric constant, upper limit on the size of the cliques to
    find. \code{NULL} means no limit.}
}
\details{
  \code{cliques} find all complete subgraphs in the input graph, obeying
  the size limitations given in the \code{min} and \code{max} arguments.

  \code{largest.cliques} finds all largest cliques in the input graph. A
  clique is largest if there is no other clique including more vertices.

  \code{maximal.cliques} finds all maximal cliques in the input graph.
  A clique in maximal if it cannot be extended to a larger clique. The
  largest cliques are always maximal, but a maximal clique is not
  neccessarily the largest.

  \code{clique.number} calculates the size of the largest clique(s).

  The current implementation of these functions searches
  for maximal independent vertex sets (see
  \code{\link{independent.vertex.sets}}) in the complementer graph.
}
\value{
  \code{cliques}, \code{largest.cliques} and \code{clique.number} return
  a list containing numeric vectors of vertex ids. Each list element is
  a clique.

  \code{clique.number} returns an integer constant.
}
%\references{TODO}
\author{Tamas Nepusz \email{ntamas@rmki.kfki.hu} and Gabor Csardi
  \email{csardi@rmki.kfki.hu} for the R interface and the manual page. }
\seealso{\code{\link{independent.vertex.sets}}}
\examples{
# this usually contains cliques of size six
g <- erdos.renyi.game(100, 0.3)
clique.number(g)
cliques(g, min=6)
largest.cliques(g)

# To have a bit less maximal cliques, about 100-200 usually
g <- erdos.renyi.game(100, 0.03)
maximal.cliques(g)

}
\keyword{graphs}
