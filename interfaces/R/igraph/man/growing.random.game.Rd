\name{growing.random.game}
\alias{growing.random.game}
\concept{Random graph model}
\title{Growing random graph generation}
\description{This function creates a random graph by simulating its
  stochastic evolution.}
\usage{
growing.random.game(n, m=1, directed=TRUE, citation=FALSE)
}
\arguments{
  \item{n}{Numeric constant, number of vertices in the graph.}
  \item{m}{Numeric constant, number of edges added in each time step.}
  \item{directed}{Logical, whether to create a directed graph.}
  \item{citation}{Logical. If \code{TRUE} a citation graph is created,
    ie. in each time step the added edges are originating from the new
    vertex. }
}
\details{
  This is discrete time step model, in each time step a new vertex is
  added to the graph and \code{m} new edges are created. If
  \code{citation} is \code{FALSE} these edges are connecting two
  uniformly randomly chosen vertices, otherwise the edges are connecting 
  new vertex to uniformly randomly chosen old vertices.
}
\value{
  A new graph object.
}
% \references{}
\author{Gabor Csardi \email{csardi@rmki.kfki.hu}}
\seealso{\code{\link{barabasi.game}}, \code{\link{erdos.renyi.game}}}
\examples{
g <- growing.random.game(500, citation=FALSE)
g2 <- growing.random.game(500, citation=TRUE)
}
\keyword{graphs}
