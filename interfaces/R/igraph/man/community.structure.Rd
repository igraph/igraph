\name{communities}
\alias{communities}
\alias{community.to.membership}
\concept{Community structure}
\title{Common functions supporting community detection algorithms}
\description{
  \code{community.to.membership} takes a merge matrix, a typical result
  of community structure detection algorithms and creates a membership
  vector by performing a given number of merges in the merge matrix.
}
\usage{
community.to.membership(graph, merges, steps, membership=TRUE, csize=TRUE)
}
\arguments{
  \item{graph}{The graph to which the merge matrix belongs.}
  \item{merges}{The merge matrix, see
    e.g. \code{\link{walktrap.community}} for the exact format.}
  \item{steps}{The number of steps, ie. merges to be performed.}
  \item{membership}{Logical scalar, whether to include the membership
    vector in the result.}
  \item{csize}{Logical scalar, whether to include the sizes of the
    communities in the result.}
}
%\details{}
\value{
  A named list with two members:
  \item{membership}{The membership vector.}
  \item{csize}{A numeric vector giving the sizes of the communities.}
}
%\references{}
\author{Gabor Csardi \email{csardi@rmki.kfki.hu}}
\seealso{
  \code{\link{walktrap.community}},
  \code{\link{edge.betweenness.community}},
  \code{\link{fastgreedy.community}},
  \code{\link{spinglass.community}} for various community detection
  methods.
}
\examples{
g <- graph.full(5) \%du\% graph.full(5) \%du\% graph.full(5)
g <- add.edges(g, c(0,5, 0,10, 5, 10))
wtc <- walktrap.community(g)
community.to.membership(g, wtc$merges, steps=12)
}
\keyword{graphs}
