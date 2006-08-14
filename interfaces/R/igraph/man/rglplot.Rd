\name{rglplot}
\alias{rglplot}
\alias{rglplot.igraph}
\title{3D plotting of graphs with OpenGL}
\description{Using the \code{rgl} package, \code{rglplot} plots a graph
  in 3D. The plot can be zoomed, rotated, shifted, etc. but the
  coordinates of the vertices is fixed.}
\usage{
rglplot (x, layout = layout.random, layout.par = list(), labels = NULL, 
    label.color = "darkblue", label.font = NULL, label.degree = -pi/4, 
    label.dist = 0, vertex.color = "SkyBlue2", vertex.size = 15, 
    edge.color = "darkgrey", edge.width = 1, edge.labels = NA, 
    \dots)
}
\arguments{
  \item{x}{The graph to plot.}
  \item{layout}{A function or matrix giving the coordinates of the
    vertices. See the same argument of \code{\link{tkplot}}.}
  \item{layout.par}{
    A named list specifying the optional parameters of the layout
    function, only used if the \code{layout} argument is a
    function. The parameters themselves are layout function dependent,
    they are discussed in the \code{\link{layout}} manual page. 
  }
  \item{labels}{Gives the labels to use for the vertices. It can be a
    numeric or character vector. If \code{NA} the plot will have no
    labels, if \code{NULL} numeric ids starting from one will be
    used. If a character vector of length one with prefix
    \dQuote{\code{a:}} is given then the corresponding vertex attribute
    (without the \dQuote{\code{a:}} prefix) is used for the labels.}
  \item{label.color}{The color of the labels. If it is a character
    scalar with prefix \dQuote{\code{a:}} then the corresponding vertex
    attribute will be used. Otherwise if it is a vector of length one
    then all labels will be the same color. It can also be a character
    or numeric vector of the same length as the number of
    vertices. Colors can be given as strings or integers, in the latter
    case the actual color palette will be used.}
  \item{label.font}{The font for the labels, currently not used.}
  \item{label.degree}{Governs the placement of the vertex labels. 0
    is right, pi/2 is down, pi is left. This argument can be a constant
    (for all labels) or a numeric vector. If it is a constant string
    with prefix \dQuote{\code{a:}} then the corresponding vertex
    argument will be used.} 
  \item{label.dist}{Numeric constant, gives the distance between the
    labels and the vertices. 1 is some \dQuote{ideal} distance,
    based on the vertex and label sizes. Use 0 for drawing the labels on
    the vertices.}
  \item{vertex.color}{The color of the vertices, either a character
    constant with prefix \dQuote{\code{a:}} for using an attribute or
    without it: a constant color for all vertices, or a numeric or
    character or numeric vector giving the different colors for the
    vertices.} 
  \item{vertex.size}{The size of the vertices. Either a numeric
    constant, or a numeric vector, or a string with prefix
    \quote{\code{a:}} giving the vertex attribute to use.}
  \item{edge.color}{The color of the edges. Either a character vector
    with prefix \dQuote{\code{a:}} naming the edge attribute to use or
    one without it, or a numeric or character vector giving the colors
    for each edge.
  }
  \item{edge.width}{The width of the edges. Either a numeric constant or
  a numeric vector, or a string with prefix \dQuote{\code{a:}} giving
  the edge attribute to use.}
  \item{edge.labels}{Labels for the edges, this will be converted to a
    character vector. The default is \code{NA} which suppresses edge
    labels. The vector might contain \code{NA}'s as well to suppress
    labels for some edges.}
  \item{\dots}{Additional arguments, none right now.}
}
\details{
  Note that \code{rglplot} is considered to be highly experimental. It
  is not very useful either.
}
\value{\code{NULL}, invisibly.}
% \references{}
\author{ Gabor Csardi \email{csardi@rmki.kfki.hu}}
\seealso{ \code{\link{plot.igraph}} for the 2D version,
  \code{\link{tkplot}} for interactive graph drawing in 2D.}
\examples{
\dontrun{
g <- graph.lattice( c(5,5,5) )
coords <- layout.fruchterman.reingold(g, dim=3)
rglplot(g, layout=coords)
}
}
\keyword{graphs}
