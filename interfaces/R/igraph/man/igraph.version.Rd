\name{igraph.version}
\alias{igraph.version}
\title{Query igraph's version string}
\description{Queries igraph's original version string. See details
  below.}
\usage{
igraph.version()
}
\details{
  The igraph version string is the same as the version of the R package
  for all realeased igraph versions. For development versions and
  nightly builds, they might differ however.

  The reason for this is, that R package version numbers are not
  flexible enough to cover in-between releases versions, e.g. alpha and
  beta versions, release candidates, etc.
}
\value{
  A character scalar, the igraph version string.
}
% \references{}
\author{ Gabor Csardi \email{csardi.gabor@gmail.com} }
% \seealso{}
\examples{
## Compare to the package version
packageDescription("igraph")$Version
igraph.version()

}
\keyword{graphs}
