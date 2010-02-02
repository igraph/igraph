\name{compare.communities}
\alias{compare.communities}
\concept{Community structure}
\title{Compares community structures using various metrics}
\description{
  This function assesses the distance between two community structures.
}
\usage{
compare.communities (comm1, comm2, method = c("vi"))
}
\arguments{
  \item{comm1}{Numeric vector. The membership vector of the first
    community structure. It should contain the community id of each
    vertex, the numbering of the communities starts with zero.} 
  \item{comm2}{Numeric vector. The membership vector of the second
    community structure, in the same format as for the previous
    argument.} 
  \item{method}{Character scalar, the comparison method to use, right
    now only the variation of information (VI) metric of Meila (2003) is
    implemented.}
}
% \details{}
\value{A real number.}
\author{Tamas Nepusz \email{ntamas@gmail.com}}
\seealso{
  \code{\link{walktrap.community}},
  \code{\link{edge.betweenness.community}},
  \code{\link{fastgreedy.community}},
  \code{\link{spinglass.community}} for various community detection
  methods.
}
\examples{
g <- graph.famous("Zachary")
sg <- spinglass.community(g)
le <- leading.eigenvector.community(g)
compare.communities(sg$membership, le$membership)
}
\keyword{graphs}
