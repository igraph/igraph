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
\description{}
\usage{
tkplot (graph, layout = layout.random, layout.par = list(), 
    labels = NULL, label.color = "darkblue", label.font = NULL, 
    label.degree = -pi/4, label.dist = 0, vertex.color = "SkyBlue2", 
    vertex.size = 15, edge.color = "darkgrey", edge.width = 1, \dots)

tkplot.close(tkp.id, window.close = TRUE) 
tkplot.off()
tkplot.fit.to.screen(tkp.id, width = NULL, height = NULL)
tkplot.reshape(tkp.id, newlayout, \dots)
tkplot.export.postcript(tkp.id)
tkplot.getcoords(tkp.id, norm = FALSE)
tkplot.center(tkp.id)
tkplot.rotate(tkp.id, degree = NULL, rad = NULL)
}
\arguments{
  \item{graph}{}
  \item{layout}{}
  \item{layout.par}{}
  \item{labels}{}
  \item{label.color}{}
  \item{label.font}{}
  \item{label.degree}{}
  \item{label.dist}{}
  \item{vertex.color}{}
  \item{vertex.size}{}
  \item{edge.color}{}
  \item{edge.width}{}
  \item{tkp.id}{}
  \item{window.close}{}
  \item{width}{}
  \item{height}{}
  \item{newlayout}{}
  \item{norm}{}
  \item{degree}{}
  \item{rad}{}
}
\details{}
\value{}
\references{}
\author{Gabor Csardi \email{csardi@rmki.kfki.hu}}
\seealso{}
\examples{}
\keyword{graphs}
