\name{layout.circle}
\alias{layout.circle}
\concept{Graph layout}
\title{Graph layout with vertices on a circle.}
\description{
  Place vertices on a circle, in the order of their vertex ids.
}
\usage{
layout.circle (graph)
}
\arguments{
  \item{graph}{The input graph.}
}
\details{
  If you want to order the vertices differently, then permute them using
  the \code{\link{permute.vertices}} function.
}
\value{
  A numeric matrix with two columns, and one row for each vertex.
}
% \references{}
\author{Gabor Csardi \email{csardi.gabor@gmail.com}}
% \seealso{}
% \examples{}
\keyword{graphs}
