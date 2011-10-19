\name{plot.igraph}
\alias{plot.igraph}
\concept{Visualization}
\title{Plotting of graphs}
\description{\code{plot.graph} is able to plot graphs to any R
  device. It is the non-interactive companion of the \code{tkplot}
  function.}
\usage{
\method{plot}{igraph}(x, axes=FALSE, xlab="", ylab="", add=FALSE,
      xlim=c(-1,1), ylim=c(-1,1), main="", sub="", \dots)
}
\arguments{
  \item{x}{The graph to plot.}
  \item{axes}{Logical, whether to plot axes, defaults to FALSE.}
  \item{xlab}{The label of the horizontal axis. Defaults to the empty string.}
  \item{ylab}{The label of the vertical axis. Defaults to the empty
    string.}
  \item{add}{Logical scalar, whether to add the plot to the current
    device, or delete the device's current contents first.}
  \item{xlim}{The limits for the horizontal axis, it is unlikely that
    you want to modify this.}
  \item{ylim}{The limits for the vertical axis, it is unlikely that you
    want to modify this.}
  \item{main}{Main title.}
  \item{sub}{Subtitle.}
  \item{\dots}{Additional arguments, passed to \code{\link{plot}}.}
}
\details{
  One convenient way to plot graphs is to plot with \code{\link{tkplot}}
  first, handtune the placement of the vertices, query the coordinates
  by the \code{\link{tkplot.getcoords}} function and use them with
  \code{plot} to plot the graph to any R device.}
\value{Returns \code{NULL}, invisibly.}
% \references{}
\author{Gabor Csardi \email{csardi@rmki.kfki.hu}}
\seealso{\code{\link{layout}} for different layouts,
  \code{\link{igraph.plotting}} for the detailed description of the
  plotting parameters and \code{\link{tkplot}} and \code{\link{rglplot}}
  for other graph plotting functions. }
\examples{
g <- graph.ring(10)
\dontrun{plot(g, layout=layout.kamada.kawai, vertex.color="green")}
}
\keyword{graphs}
