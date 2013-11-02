\name{cliques}
\alias{cliques}
\alias{largest.cliques}
\alias{maximal.cliques}
\alias{maximal.cliques.count}
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
maximal.cliques(graph, min=NULL, max=NULL, subset=NULL, file=NULL)
maximal.cliques.count(graph, min=NULL, max=NULL, subset=NULL)
clique.number(graph)
}
\arguments{
  \item{graph}{The input graph, directed graphs will be considered as
    undirected ones, multiple edges and loops are ignored.}
  \item{min}{Numeric constant, lower limit on the size of the cliques to
    find. \code{NULL} means no limit, ie. it is the same as 0.}
  \item{max}{Numeric constant, upper limit on the size of the cliques to
    find. \code{NULL} means no limit.}
  \item{subset}{If not \code{NULL}, then it must be a vector of vertex
    ids, numeric or symbolic if the graph is named. The algorithm is run
    from these vertices only, so only a subset of all maximal cliques is
    returned. See the Eppstein paper for details. This argument makes it
    possible to easily parallelize the finding of maximal cliques.}
  \item{file}{If not \code{NULL}, then it must be a file name, i.e. a
    character scalar. The output of the algorithm is written to this
    file. (If it exists, then it will be overwritten.) Each clique will
    be a separate line in the file, given with the numeric ids of its
    vertices, separated by whitespace.}
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

  \code{maximal.cliques.count} counts the maximal cliques.

  \code{clique.number} calculates the size of the largest clique(s).

  The current implementation of these functions searches
  for maximal independent vertex sets (see
  \code{\link{independent.vertex.sets}}) in the complementer graph.
}
\value{
  \code{cliques}, \code{largest.cliques} and \code{clique.number} return
  a list containing numeric vectors of vertex ids. Each list element is
  a clique.

  \code{maximal.cliques} returns \code{NULL}, invisibly, if its
  \code{file} argument is not \code{NULL}. The output is written to the
  specified file in this case.
  
  \code{clique.number} and \code{maximal.cliques.count} return an
  integer scalar.
}
\references{
  For maximal cliques the following algorithm is implemented:
  David Eppstein, Maarten Loffler, Darren Strash:
  Listing All Maximal Cliques in Sparse Graphs in Near-optimal Time.
  \url{http://arxiv.org/abs/1006.5440}
}
\author{Tamas Nepusz \email{ntamas@gmail.com} and Gabor Csardi
  \email{csardi.gabor@gmail.com}}
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
