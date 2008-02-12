\name{degree}
\alias{degree}
\concept{Vertex degree}
\concept{Degree distribution}
\alias{degree.distribution}
\title{Degree and degree distribution of the vertices}
\description{The degree of a vertex is its most basic structural
  property, the number of its adjacent edges.}
\usage{
degree(graph, v=V(graph), mode = c("all", "out", "in", "total"), loops = TRUE) 
degree.distribution(graph, cumulative = FALSE, \dots) 
}
\arguments{
  \item{graph}{The graph to analyze.}
  \item{v}{The ids of vertices of which the degree will be calculated.}
  \item{mode}{Character string, \dQuote{out} for out-degree, \dQuote{in} for
    in-degree or \dQuote{total} for the sum of the two. For undirected
    graphs this argument is ignored.}
  \item{loops}{Logical; whether the loop edges are also counted.}
  \item{cumulative}{Logical; whether the cumulative degree distribution
    is to be calculated.}
  \item{\dots}{Additional arguments to pass to \code{degree},
    eg. \code{mode} is useful but also \code{v} and \code{loops} make
    sense.}
}
% \details{}
\value{
  For \code{degree} a numeric vector of the same length as argument
  \code{v}.

  For \code{degree.distribution} a numeric vector of the same length as
  the maximum degree plus one. The first element is the relative
  frequency zero degree vertices, the second vertices with degree one, 
  etc. 
}
% \references{}
\author{Gabor Csardi \email{csardi@rmki.kfki.hu}}
% \seealso{}
\examples{
g <- graph.ring(10)
degree(g)
g2 <- erdos.renyi.game(1000, 10/1000)
degree.distribution(g2)
}
\keyword{graphs}
