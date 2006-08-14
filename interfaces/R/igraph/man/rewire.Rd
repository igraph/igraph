\name{rewire}
\alias{rewire}
\title{Graph rewiring}
\description{Randomly rewires a graph while preserving the degree
  distribution.}
\usage{
rewire(graph, mode = "simple", niter = 100) 
}
\arguments{
  \item{graph}{The graph to be rewired.}
  \item{mode}{The rewiring algorithm to be used. It can be one of the
    following: \code{simple}: simple rewiring algorithm which
    chooses two arbitrary edges in each step (namely (a,b) and (c,d))
    and substitutes them with (a,d) and (c,b) if they don't yet exist.}
  \item{niter}{Number of rewiring trials to perform.}
}
\details{
This function generates a new graph based on the original one by
randomly rewiring edges while preserving the original graph's degree
distribution.
}
\value{A new graph object.}
% \references{}
\author{Tamas Nepusz \email{ntamas@rmki.kfki.hu} and
  Gabor Csardi \email{csardi@rmki.kfki.hu}}
\seealso{\code{\link{degree.sequence.game}}}
\examples{
g <- graph.ring(20)
g2 <- rewire(g, niter=3)
}
\keyword{graphs}
