\name{sgm}
\concept{Seeded Graph Matching}
\title{Match Graphs given a seeding of vertex correspondences}
\description{The approximate graph matching problem is to find a bijection between the vertices of two graphs, such that the number of edge disagreements between the corresponding vertex pairs are minimized. Given a part of the bijection that consist of known correspondences (the seeds), \code{sgm} completes the bijection by computing the permutation matrix that permutes the rows and columns of the adjacency matrix  }
\usage{
sgm( \dots)
}
\arguments{
  \item{n}{The number of vertices in the graph.}
  
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
