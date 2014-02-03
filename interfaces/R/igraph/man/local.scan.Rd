\name{local.scan}
\concept{Scan Statistics on Graphs}
\title{Compute local scan statistics on static and time-varying graphs}
\description{The scan statistic is a summary of the locality statistics that is computed from the local neighborhood of each vertex  }
\usage{
local.scan( graph.us)
}
\arguments{
  \item{graph.us}{}
  \item{graph.them}{}
  \item{k}{ The size of the local neighborhood for each vertex}
  \item{FUN}{}
  \item{weighted}{}
  \item{mode}{}
  \item{neighborhoods}{}
  \item{\dots}
}
\details{
  In G(n,p) graphs, the graph has \sQuote{n} vertices and for each edge
  the probability that it is present in the graph is \sQuote{p}.

  In G(n,m) graphs, the graph has \sQuote{n} vertices and \sQuote{m}
  edges, and the \sQuote{m} edges are chosen uniformly randomly from the
  set of all possible edges. This set includes loop edges as well if the
  \code{loops} parameter is TRUE.
  
  \code{random.graph.game} is an alias to this function.
}
\value{A graph object.}
\references{ Erdos, P. and Renyi, A., On random graphs, \emph{Publicationes
  Mathematicae} 6, 290--297 (1959).
}
\author{Gabor Csardi \email{csardi.gabor@gmail.com}}
\seealso{\code{\link{barabasi.game}}}
\examples{
g <- erdos.renyi.game(1000, 1/1000)
degree.distribution(g)
}
\keyword{graphs}
