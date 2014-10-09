\name{groups}
\alias{groups}
\alias{groups.default}
\alias{groups.communities}
\title{Groups of a vertex partitioning}
\description{Create a list of vertex groups from some graph clustering
  or community structure.}
\usage{
groups(x)
\method{groups}{default}(x)
\method{groups}{communities}(x)
}
\arguments{
  \item{x}{Some object that represents a grouping of the vertices. See
    details below.}
}
\details{
  Currently two methods are defined for this function. The default
  method works on the output of \code{\link{components}}. (In fact it
  works on any object that is a list with an entry called
  \code{membership}.)

  The second method works on \code{\link{communities}} objects.
}
\value{
  A named list of numeric or character vectors. The names are just numbers that
  refer to the groups. The vectors themselves are numeric or symbolic
  vertex ids.
}
\author{
  Gabor Csardi \email{csardi.gabor@gmail.com}
}
\seealso{
  \code{\link{components}} and the various community finding functions.
}
\examples{
g <- make_graph("Zachary")
fgc <- cluster_fast_greedy(g)
groups(fgc)

g2 <- ring(10) + full_graph(5)
groups(components(g2))
}
\keyword{graphs}
