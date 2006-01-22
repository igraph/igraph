\name{plot.igraph}
\alias{plot.igraph}
\title{Plotting of graphs}
\description{\code{plot.graph} is able to plot graphs to any R
  device. It is the non-interactive companion of the \code{tkplot}
  function.}
\usage{
plot.igraph(x, layout = layout.random, layout.par = list(), 
    labels = NULL, label.color = "darkblue", label.font = NULL, 
    label.degree = -pi/4, label.dist = 0, vertex.color = "SkyBlue2", 
    vertex.size = 15, edge.color = "darkgrey", edge.width = 1,
    edge.labels = NA, axes = FALSE, xlab = "", ylab = "", \dots)
}
\arguments{
  \item{x}{The graph to plot.}
  \item{layout}{The coordinates of the vertices, or the method to
    calculate them. See the same argument of \code{\link{tkplot}}.}
  \item{layout.par}{List, the parameters of the layout function. See the
    same argument of \code{\link{tkplot}}.}
  \item{labels}{The labels of the vertices. See the same argument of
    \code{\link{tkplot}}.}
  \item{label.color}{The color(s) of the labels. See the same argument
    of \code{\link{tkplot}}.}
  \item{label.font}{The font to be used for the labels. Right now all
    labels use the same font.}
  \item{label.degree}{Specifies the placement of the labels of the
    vertices. See the same argument of \code{\link{tkplot}}.}
  \item{label.dist}{Specifies the placement of the vertex labels. See
    the same argument of \code{\link{tkplot}}.}
  \item{vertex.color}{The color(s) of the vertices. See the same
    argument of \code{\link{tkplot}}.}
  \item{vertex.size}{The size(s) of the vertices. See the same argument
    of \code{\link{tkplot}}.}
  \item{edge.color}{The color(s) of the edges. See the same argument of
    \code{\link{tkplot}}. }
  \item{edge.width}{The width(s) of the edges. See the same argument of
    \code{\link{tkplot}}.}
  \item{edge.labels}{Labels for the edges, this will be converted to a
    character vector. The default is \code{NA} which suppresses edge
    labels. The vector might contain \code{NA}'s as well to suppress
    labels for some edges.}  
  \item{axes}{Logical, whether to plot axes, defaults to FALSE.}
  \item{xlab}{The label of the horizontal axis. Defaults to the empty string.}
  \item{ylab}{The label of the vertical axis. Defaults to the empty
    string.}
  \item{\dots}{Additional arguments, passed to \code{\link{plot}}.}
}
\details{
  One convenient way to plot graphs is to plot with \code{\link{tkplot}}
  first, handtune the placement of the vertices, query the coordinates
  by the \code{\link{tkplot.getcoords}} function and use them with
  \code{plot} to plot the graph to any R device.
  
  The additional arguments specified are passed to the \code{plot}
  function.} 
\value{Returns \code{NULL}, invisibly.}
% \references{}
\author{Gabor Csardi \email{csardi@rmki.kfki.hu}}
\seealso{\code{\link{layout}} for different layouts,
  \code{\link{tkplot}} for the better description of the parameters and
  an interactive version.}
\examples{
g <- graph.ring(10)
\dontrun{plot(g, layout=layout.kamada.kawai, vertex.color="green")}
}
\keyword{graphs}
