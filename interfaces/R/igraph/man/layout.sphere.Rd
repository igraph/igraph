\name{layout.sphere}
\alias{layout.sphere}
\concept{Graph layout}
\title{Graph layout with vertices on the surface of a sphere}
\description{
  Place vertices on a sphere, approximately uniformly,
  in the order of their vertex ids.
}
\usage{
layout.sphere (graph)
}
\arguments{
  \item{graph}{The input graph.}
}
\details{
  \code{layout.sphere} places the vertices (approximately) uniformly on
  the surface of a sphere, this is thus a 3d layout. It is not clear
  however what \dQuote{uniformly on a sphere} means.
  
  If you want to order the vertices differently, then permute them using
  the \code{\link{permute.vertices}} function.
}
\value{
  A numeric matrix with three columns, and one row for each vertex.
}
% \references{}
\author{Gabor Csardi \email{csardi.gabor@gmail.com}}
% \seealso{}
% \examples{}
\keyword{graphs}
