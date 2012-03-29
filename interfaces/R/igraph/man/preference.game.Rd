\name{preference.game}
\alias{preference.game}
\alias{asymmetric.preference.game}
\title{Trait-based random generation}
\description{Generation of random graphs based on different vertex types.}
\usage{
preference.game(nodes, types, type.dist=rep(1, types),
    pref.matrix=matrix(1, types, types), directed=FALSE, loops=FALSE) 
asymmetric.preference.game(nodes, types,
    type.dist.matrix=matrix(1,types,types),
    pref.matrix = matrix(1, types, types), loops=FALSE) 
}
\arguments{
  \item{nodes}{The number of vertices in the graphs.}
  \item{types}{The number of different vertex types.}
  \item{type.dist}{The distribution of the vertex types, a numeric
    vector of length \sQuote{types} containing non-negative numbers. The
    vector will be normed to obtain probabilities.}
  \item{type.dist.matrix}{The joint distribution of the in- and
    out-vertex types.}
  \item{pref.matrix}{A square matrix giving the preferences of the
    vertex types. The matrix has \sQuote{types} rows and columns.}
  \item{directed}{Logical constant, whether to create a directed graph.}
  \item{loops}{Logical constant, whether self-loops are allowed in the
    graph.}
}
\details{
  Both models generate random graphs with given vertex types. For
  \code{preference.game} the probability that two vertices will be
  connected depends on their type and is given by the
  \sQuote{pref.matrix} argument. This matrix should be symmetric to make
  sense but this is not checked. The distribution of the different
  vertes types is given by the \sQuote{type.dist} vector.

  For \code{asymmetric.preference.game} each vertex has an in-type and
  an out-type and a directed graph is created. The probability that a
  directed edge is realized from a vertex with a given out-type to a
  vertex with a given in-type is given in the \sQuote{pref.matrix}
  argument, which can be asymmetric. The joint distribution for the in-
  and out-types is given in the \sQuote{type.dist.matrix} argument.
}
\value{An igraph graph.}
%\references{}
\author{Tamas Nepusz \email{ntamas@rmki.kfki.hu} and Gabor Csardi
  \email{csardi@rmki.kfki.hu} for the R interface} 
\examples{
pf <- matrix( c(1, 0, 0, 1), nrow=2)
g <- preference.game(20, 2, pref.matrix=pf)
\dontrun{tkplot(g, layout=layout.fruchterman.reingold)}

pf <- matrix( c(0, 1, 0, 0), nrow=2)
g <- asymmetric.preference.game(20, 2, pref.matrix=pf)
\dontrun{tkplot(g, layout=layout.circle)}
}
\seealso{\code{\link{establishment.game}}.
  \code{\link{callaway.traits.game}}
}
\keyword{graphs}
