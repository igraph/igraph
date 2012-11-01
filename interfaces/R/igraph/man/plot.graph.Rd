\name{plot.igraph}
\alias{plot.igraph}
\concept{Visualization}
\title{Plotting of graphs}
\description{\code{plot.graph} is able to plot graphs to any R
  device. It is the non-interactive companion of the \code{tkplot}
  function.}
\usage{
\S3method{plot}{igraph}(x, axes=FALSE, xlab="", ylab="", add=FALSE,
      xlim=c(-1,1), ylim=c(-1,1), main="", sub="",
      mark.groups=list(), mark.shape=1/2, 
      mark.col=rainbow(length(mark.groups), alpha=0.3),
      mark.border=rainbow(length(mark.groups), alpha=1)
      mark.expand=15, 
      \dots)
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
  \item{mark.groups}{A list of vertex id vectors. It is interpreted as a
    set of vertex groups. Each vertex group is highlighted, by plotting
    a colored smoothed polygon around and \dQuote{under} it. See the
    arguments below to control the look of the polygons.}
  \item{mark.shape}{A numeric scalar or vector. Controls the smoothness
    of the vertex group marking polygons. This is basically the \sQuote{shape}
    parameter of the \code{\link[graphics]{xspline}} function, its
    possible values are between -1 and 1. If it is a vector, then a
    different value is used for the different vertex groups.}
  \item{mark.col}{A scalar or vector giving the colors of marking the
    polygons, in any format accepted by
    \code{\link[graphics]{xspline}}; e.g. numeric color ids, symbolic
    color names, or colors in RGB.}
  \item{mark.border}{A scalar or vector giving the colors of the borders
    of the vertex group marking polygons. If it is \code{NA}, then no
    border is drawn.}
  \item{mark.expand}{A numeric scalar or vector, the size of the border
    around the marked vertex groups. It is in the same units as the vertex
    sizes. If a vector is given, then different values are used for the
    different vertex groups.} 
  \item{\dots}{Additional plotting parameters. See
    \link{igraph.plotting} for the complete list.}
}
\details{
  One convenient way to plot graphs is to plot with \code{\link{tkplot}}
  first, handtune the placement of the vertices, query the coordinates
  by the \code{\link{tkplot.getcoords}} function and use them with
  \code{plot} to plot the graph to any R device.}
\value{Returns \code{NULL}, invisibly.}
% \references{}
\author{Gabor Csardi \email{csardi.gabor@gmail.com}}
\seealso{\code{\link{layout}} for different layouts,
  \code{\link{igraph.plotting}} for the detailed description of the
  plotting parameters and \code{\link{tkplot}} and \code{\link{rglplot}}
  for other graph plotting functions. }
\examples{
g <- graph.ring(10)
\dontrun{plot(g, layout=layout.kamada.kawai, vertex.color="green")}
}
\keyword{graphs}
