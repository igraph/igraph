\name{as.igraph}
\alias{as.igraph}
\alias{as.igraph.igraphHRG}
\concept{Hierarchical random graphs}
\title{Conversion to igraph}
\description{These fucntions convert various objects to igraph graphs.}
\usage{
\method{as.igraph}{igraphHRG}(x, \dots)
}
\arguments{
  \item{x}{The object to convert.}
  \item{\dots}{Additional arguments. None currently.}
}
\details{
  You can use \code{as.igraph} to convert various objects to igraph
  graphs. Right now the following objects are supported:
  \itemize{
    \item code{igraphHRG} These objects are created by the
    \code{\link{hrg.fit}} and \code{\link{hrg.consensus}} functions.
  }
}
\value{
  All these functions return an igraph graph.
}
% \references{}
\author{Gabor Csardi \email{csardi.gabor@gmail.com}.}
% \seealso{}
\examples{
g <- graph.full(5) + graph.full(5)
hrg <- hrg.fit(g)
as.igraph(hrg)
}
\keyword{graphs}
