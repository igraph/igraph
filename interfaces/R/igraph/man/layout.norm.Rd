\name{layout.norm}
\alias{layout.norm}
\concept{Graph layout}
\title{Normalize coordinates for plotting graphs}
\description{Rescale coordinates linearly to be within given bounds.}
\usage{
layout.norm(layout, xmin = -1, xmax = 1, ymin = -1, ymax = 1,
      zmin = -1, zmax = 1)
}
\arguments{
  \item{layout}{A matrix with two or three columns, the layout to
    normalize.}
  \item{xmin,xmax}{The limits for the first coordinate, if one of them
    or both are \code{NULL} then no normalization is performed along
    this direction.}
  \item{ymin,ymax}{The limits for the second coordinate, if one of them
    or both are \code{NULL} then no normalization is performed along
    this direction.}
  \item{zmin,zmax}{The limits for the third coordinate, if one of them
    or both are \code{NULL} then no normalization is performed along
    this direction.}  
}
\details{
  \code{layout.norm} normalizes a layout, it linearly transforms each
  coordinate separately to fit into the given limits.
}
\value{A numeric matrix with at the same dimension as \code{layout}.}
%\references{}
\author{Gabor Csardi \email{csardi.gabor@gmail.com}}
% \seealso{}
% \examples{}
\keyword{graphs}
