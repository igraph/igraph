\name{growing.random.game}
\alias{growing.random.game}
\title{Growing random graph generation}
\description{This function creates a random graph by simulating its
  stochastic evolution.}
\usage{
growing.random.game(n, m=1, citation=FALSE, \dots)
}
\arguments{
  \item{n}{Numeric constant, number of vertices in the graph.}
  \item{m}{Numeric constant, number of edges added in each time step.}
  \item{citation}{Logical. If \code{TRUE} a citation graph is created,
    ie. in each time step the added edges are originating from the new
    vertex. }
  \item{\dots}{Additional arguments, passed to the graph constructor.}
}
\details{
  This is discrete time step model, in each time step a new vertex is
  added to the graph and \code{m} new edges are created. If
  \code{citation} is \code{FALSE} these edges are connecting two
  uniformly randomly chosen vertices, otherwise the edges are connecting 
  new vertex to uniformly randomly chosen old vertices. For directed
  graphs, the edges are always originating from the new vertex.
}
\value{
  A new graph object.
}
% \references{}
\author{Gabor Csardi \email{csardi@rmki.kfki.hu}}
\seealso{\code{\link{ba.game}}, \code{\link{erdos.renyi.game}}}
\examples{
g <- growing.random.game(500, citation=FALSE)
g2 <- growing.random.game(500, citation=TRUE)
}
\keyword{graphs}
