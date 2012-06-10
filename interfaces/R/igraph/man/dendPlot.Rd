\name{dendPlot}
\alias{dendPlot}
\concept{Dendrograms}
\title{Plot dendrograms}
\description{
  This is generic function that can plot various objects as dendrograms.
}
% \usage
% \arguments{}
\details{
  Currently the function is defined for \code{communities} (see
  \code{\link{dendPlot.communities}}) and \code{igraphHRG}
  (see \code{\link{dendPlot.igraphHRG}}) objects.
}
% \value{}
%\references{}
\author{Gabor Csardi \email{csardi.gabor@gmail.com}}
% \seealso{}
\examples{
karate <- graph.famous("Zachary")
fc <- fastgreedy.community(karate)
dendPlot(fc)
}
\keyword{graphs}
