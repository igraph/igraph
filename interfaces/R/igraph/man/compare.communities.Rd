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
  \item{method}{Character scalar, the comparison method to use. Possible
    values: \sQuote{vi} is the variation of information (VI) metric of
    Meila (2003), \sQuote{nmi} is the normalized mutual information
    measure proposed by Danon et al. (2005).}
}
% \details{}
\value{A real number.}
\references{
  Meila M: Comparing clusterings by the variation of information.
  In: Schölkopf B, Warmuth MK (eds.). Learning Theory and Kernel Machines:
  16th Annual Conference on Computational Learning Theory and 7th Kernel
  Workshop, COLT/Kernel 2003, Washington, DC, USA. Lecture Notes in Computer
  Science, vol. 2777, Springer, 2003. ISBN: 978-3-540-40720-1.

  Danon L, Diaz-Guilera A, Duch J, Arenas A: Comparing community structure
  identification. J Stat Mech P09008,   
}
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
