\name{compare.communities}
\alias{compare.communities}
\concept{Community structure}
\title{Compares community structures using various metrics}
\description{
  This function assesses the distance between two community structures.
}
\usage{
compare.communities (comm1, comm2, method = c("vi", "nmi", "split.join",
                                      "rand", "adjusted.rand"))
}
\arguments{
  \item{comm1}{A \code{\link{communities}} object containing a community
    structure; or a numeric vector, the membership vector of the first
    community structure. The membership vector should contain the
    community id of each vertex, the numbering of the communities starts
    with one.}
  \item{comm2}{A \code{\link{communities}} object containing a community
    structure; or a numeric vector, the membership vector of the second
    community structure, in the same format as for the previous
    argument.} 
  \item{method}{Character scalar, the comparison method to use. Possible
    values: \sQuote{vi} is the variation of information (VI) metric of
    Meila (2003), \sQuote{nmi} is the normalized mutual information
    measure proposed by Danon et al. (2005), \sQuote{split.join} is the
    split-join distance of can Dongen (2000), \sQuote{rand} is the Rand
    index of Rand (1971), \sQuote{adjusted.rand} is the adjusted Rand
    index by Hubert and Arabie (1985).}
}
% \details{}
\value{A real number.}
\references{
  Meila M: Comparing clusterings by the variation of information.
  In: Scholkopf B, Warmuth MK (eds.). \emph{Learning Theory and Kernel
    Machines: 16th Annual Conference on Computational Learning Theory
    and 7th Kernel Workshop}, COLT/Kernel 2003, Washington, DC,
  USA. Lecture Notes in Computer Science, vol. 2777, Springer,
  2003. ISBN: 978-3-540-40720-1.

  Danon L, Diaz-Guilera A, Duch J, Arenas A: Comparing community structure
  identification. \emph{J Stat Mech} P09008, 2005.

  van Dongen S: Performance criteria for graph clustering and Markov cluster
  experiments. Technical Report INS-R0012, National Research Institute for
  Mathematics and Computer Science in the Netherlands, Amsterdam, May 2000.
 
  Rand WM: Objective criteria for the evaluation of clustering methods.
  \emph{J Am Stat Assoc} 66(336):846-850, 1971.

  Hubert L and Arabie P: Comparing partitions. \emph{Journal of
    Classification} 2:193-218, 1985.
  
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
compare(sg, le, method="rand")
compare(membership(sg), membership(le))
}
\keyword{graphs}
