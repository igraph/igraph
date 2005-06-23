\name{degree.sequence.game}
\alias{degree.sequence.game}
\title{Generate random graphs with a given degree sequence}
\description{It is often to create a graph with given vertex
  degrees. This is exactly what \code{degree.sequence.game} does.}
\usage{
degree.sequence.game(out.deg, in.deg = NULL, method = "simple", \dots) 
}
\arguments{
  \item{out.deg}{Numeric vector, the sequence of degrees (for undirected
  graphs) or out-degrees (for directed graphs). For undirected graphs
  its sum should be even. For directed graphs its sum should be the same
  as the sum of \code{in.deg}.}
  \item{in.deg}{For directed graph, the in-degree sequence.}
  \item{method}{Character, the method for generating the graph. Right
    now only \dquote{simple} is implemented.}
}
\details{The dquote{simple} method simply connects the out-stubs of the
  edges (undirected graphs) or the out-stubs and in-stubs (directed
  graphs) together. This way loops edges and also multiple edges will be
  generated.}
\value{The new graph object.}
%\references{}
\author{Gabor Csardi \email{csardi@rmki.kfki.hu}}
\seealso{\code{\link{erdos.renyi.game}}, \code{\link{barabasi.game}}}
\examples{
g <- degree.sequence.game(rep(2,100))
degree(g)
g2 <- degree.sequence.game(1:10, 10:1)
degree(g2, mode="out")
degree(g2, mode="in")
}
\keyword{graphs}
