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
\concept{Visualization}
\title{Interactive plotting of graphs}
\description{\code{tkplot} and its companion functions serve as an
  interactive graph drawing facility. Not all parameters of the plot can
  be changed interactively right now though, eg. the colors of vertices,
  edges, and also others have to be pre-defined.}
\usage{
tkplot(graph, canvas.width=450, canvas.height=450, \dots)

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
  \item{canvas.width,canvas.height}{The size of the tkplot
    drawing area.}
  \item{tkp.id}{The id of the tkplot window to close/reshape/etc.}
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

  It's handling should be quite straightforward most of the time, here
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

  The \code{tkplot} command creates a new Tk window with the graphical
  representation of \code{graph}. The command returns an integer number,
  the tkplot id. The other commands utilize this id to be able to
  query or manipulate the plot.

  \code{tkplot.close} closes the Tk plot with id \code{tkp.id}.

  \code{tkplot.off} closes all Tk plots.

  \code{tkplot.fit.to.screen} fits the plot to the given rectange
  (\code{width} and \code{height}), if some of these are \code{NULL}
  the actual phisical width od height of the plot window is used.

  \code{tkplot.reshape} applies a new layout to the plot, its optional
  parameters will be collected to a list analogous to \code{layout.par}.

  \code{tkplot.export.postscript} creates a dialog window for saving the
  plot in postscript format.

  \code{tkplot.getcoords} returns the coordinates of the vertices in a
  matrix. Each row corresponds to one vertex.

  \code{tkplot.center} shifts the figure to the center of its plot window.

  \code{tkplot.rotate} rotates the figure, its parameter can be given
  either in degrees or in radians.
}
\value{
  \code{tkplot} returns an integer, the id of the plot, this can be used
  to manipulate it from the command line.

  \code{tkplot.getcoords} returns a matrix with the coordinates.

  \code{tkplot.close}, \code{tkplot.off}, \code{tkplot.fit.to.screen},
  \code{tkplot.reshape}, \code{tkplot.export.postscript},
  \code{tkplot.center} and \code{tkplot.rotate} return \code{NULL}
  invisibly. 
}
% \references{}
\author{Gabor Csardi \email{csardi@rmki.kfki.hu}}
\seealso{\code{\link{plot.igraph}}, \code{\link{layout}}}
\examples{
g <- graph.ring(10)
\dontrun{tkplot(g)}
}
\keyword{graphs}
