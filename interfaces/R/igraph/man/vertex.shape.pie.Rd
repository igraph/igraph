\name{Pie charts as vertices}
\alias{vertex.shape.pie}
\concept{Vertex shapes}
\title{Using pie charts as vertices in graph plots}
\description{
  More complex vertex images can be used to express addtional
  information about vertices. E.g. pie charts can be used as vertices,
  to denote vertex classes, fuzzy classification of vertices, etc.
}
\details{
  The vertex shape \sQuote{pie} makes igraph draw a pie chart for every
  vertex. There are some extra graphical vertex parameters that specify
  how the pie charts will look like:
  \describe{
    \item{pie}{Numeric vector, gives the sizes of the pie slices.}
    \item{pie.color}{A list of color vectors to use for the pies. If it
      is a list of a single vector, then this is used for all pies. It
      the color vector is shorter than the number of areas in a pie,
      then it is recycled.}
    \item{pie.border}{The color of the border line of the pie charts, in
      the same format as \code{pie.color}.}
    \item{pie.angle}{The slope of shading lines, given as an angle in
      degrees (counter-clockwise).}
    \item{pie.density}{The density of the shading lines, in lines per
      inch. Non-positive values inhibit the drawing of shading lines.}
    \item{pie.lty}{The line type of the border of the slices.}
  }
}
\author{ Gabor Csardi \email{csardi.gabor@gmail.com} }
\seealso{ \code{\link{igraph.plotting}}, \code{\link{plot.igraph}} }
\examples{
g <- make_ring(10)
values <- lapply(1:10, function(x) sample(1:10,3))
if (interactive()) {
  plot(g, vertex.shape="pie", vertex.pie=values,
       vertex.pie.color=list(heat.colors(5)),
       vertex.size=seq(10,30,length=10), vertex.label=NA)
}
}
\keyword{graphs}
