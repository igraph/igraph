\name{tkplot}
\alias{tkplot}
\alias{tkplot.close}
\alias{tkplot.off}
\alias{tkplot.fit.to.screen}
\alias{tkplot.reshape}
\alias{tkplot.export.postscript}
\alias{tkplot.getcoords}
\alias{tkplot.center}
\alias{tkplot.rotate}
\title{Interactive plotting of graphs}
\description{\code{tkplot} and its companion functions serve as an
  interactive graph drawing facility. Not all parameters of the plot can
  be changed interactively right now though, eg. the colors of vertices,
  edges, and also others have to be pre-defined.}
\usage{
tkplot (graph, layout = layout.random, layout.par = list(), 
    labels = NULL, label.color = "darkblue", label.font = NULL, 
    label.degree = -pi/4, label.dist = 0, vertex.color = "SkyBlue2", 
    vertex.size = 15, edge.color = "darkgrey", edge.width = 1,
    edge.labels=NA, \dots)

tkplot.close(tkp.id, window.close = TRUE) 
tkplot.off()
tkplot.fit.to.screen(tkp.id, width = NULL, height = NULL)
tkplot.reshape(tkp.id, newlayout, \dots)
tkplot.export.postscript(tkp.id)
tkplot.getcoords(tkp.id, norm = FALSE)
tkplot.center(tkp.id)
tkplot.rotate(tkp.id, degree = NULL, rad = NULL)
}
\arguments{
  \item{graph}{The \code{graph} to plot.}
  \item{layout}{Either a function, or a character string or a numeric
    matrix. It specifies how the vertices will be placed on the plot.

    If it is a function, this function will be called with the
    \code{graph} and \code{layout.param} as parameters to determine the
    actual coordinates.

    If it is a character sting starting with \dQuote{\code{a:}} then the
    remaining part of the string specifies the name of the vertex
    attribute to use as coordinates.

    If it is a numeric matrix, then the matrix has to have one line for
    each vertex, specifying its coordinates. The matrix should have at
    least two columns, for the \code{x} and \code{y} coordinates, but it
    can also have more (eg. for a future \code{z} coordinate).
  }
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
  \item{label.font}{The font for the labels, a tkfont object, by default
    16pt Helvetica bold is used.}
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
  \item{tkp.id}{The id of the tkplot window to close.}
  \item{window.close}{Leave this on the default value.}
  \item{width}{The width of the rectangle for generating new
    coordinates.}
  \item{height}{The height of the rectangle for generating new
    coordinates.}
  \item{newlayout}{The new layout, see the \code{layout} parameter of
    tkplot.}
  \item{norm}{Logical, should we norm the coordinates.}
  \item{degree}{The degree to rotate the plot.}
  \item{rad}{The degree to rotate the plot, in radian.}
  \item{\dots}{Not used right now.}
}
\details{
  \code{tkplot} is an interactive graph drawing facility. It is not very
  well developed at this stage, but it should be still useful.

  It's hangling should be quite straightforward most of the time, here
  are some remarks and hints.

  There are different popup menus, activated by the right mouse button,
  for vertices and edges. Both operate on the current selection if the
  vertex/edge under the cursor is part of the selection and operate on
  the vertex/edge under the cursor if it is not.

  One selection can be active at a time, either a vertex or an edge
  selection. A vertex/edge can be added to a selection by holding the
  \code{control} key while clicking on it with the left mouse
  button. Doing this again deselect the vertex/edge.

  Selections can be made also from the \code{Select} menu. The `Select
  some vertices' dialog allows to give an expression for the vertices to
  be selected: this can be a list of numeric R expessions separated by
  commas, like `\code{1,2:10,12,14,15}' for example. Similarly in the
  `Select some edges' dialog two such lists can be given and all edges
  connecting a vertex in the first list to one in the second list will
  be selected.

  In the color dialog a color name like 'orange' or RGB notation can
  also be used.
}
\value{
  \code{tkplot} returns an integer, the id of the plot, this can be used
  to manipulate it from the command line.

  \code{tkplot.close} does not return a value (?).

  \code{tkplot.off} returns TRUE.

  \code{tkplot.fit.to.screen} returns TRUE.

  \code{tkplot.reshape} does not return a value (?).

  \code{tkplot.export.postscript} does not return a value (?).

  \code{tkplot.getcoords} returns a matrix with the coordinates.

  \code{tkplot.center} returns TRUE.

  \code{tkplot.rotate} returns TRUE.
}
% \references{}
\author{Gabor Csardi \email{csardi@rmki.kfki.hu}}
\seealso{\code{\link{plot.igraph}}, \code{\link{layout}}}
\examples{
g <- graph.ring(10)
\dontrun{tkplot(g)}
}
\keyword{graphs}
