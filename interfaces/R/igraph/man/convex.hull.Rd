\name{convex.hull}
\alias{convex.hull}
\concept{Convex hull}
\title{Convex hull of a set of vertices}
\description{
  Calculate the convex hull of a set of points, i.e. the covering
  polygon that has the smallest area.
}
\usage{
convex.hull(data)
}
\arguments{
  \item{data}{The data points, a numeric matrix with two columns.}
}
%\details{}
\value{
  A named list with components:
  \item{resverts}{The indices of the input vertices that constritute the
    convex hull.}
  \item{rescoords}{The coordinates of the corners of the convex hull.}
}
\references{
 Thomas H. Cormen, Charles E. Leiserson, Ronald L. Rivest, and Clifford
 Stein. Introduction to Algorithms, Second Edition. MIT Press and
 McGraw-Hill, 2001. ISBN 0262032937. Pages 949-955 of section 33.3:
 Finding the convex hull.
}
\author{Tamas Nepusz \email{ntamas@gmail.com}}
% \seealso
\examples{
M <- cbind( runif(100), runif(100) )
convex.hull(M)
}
\keyword{graphs}
