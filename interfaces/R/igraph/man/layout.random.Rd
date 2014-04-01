\name{layout.random}
\alias{layout.random}
\concept{Graph layout}
\title{Randomly place vertices on a plane or in 3d space}
\description{
  This function uniformly randomly places the vertices of the graph
  in two or three dimensions.}
\usage{
layout.random (graph, dim=2)
}
\arguments{
  \item{graph}{The input graph.}
  \item{dim}{Integer scalar, the dimension of the space to use. It must
    be 2 or 3.}
}
\details{
  Randomly places vertices on a [-1,1] square (in 2d) or in a cube (in
  3d). It is probably a useless layout, but it can use as a starting
  point for other layout generators.}
\value{
  A numeric matrix with two or three columns.
}
% \references{}
\author{Gabor Csardi \email{csardi.gabor@gmail.com}}
% \seealso{}
\keyword{graphs}
