\name{minimum.size.separators}
\alias{minimum.size.separators}
\concept{Minimum size vertex separator}
\concept{Vertex separator}
\title{Minimum size vertex separators}
\description{Find all vertex sets of minimal size whose removal
  separates the graph into more components}
\usage{
minimum.size.separators (graph)
}
\arguments{
  \item{graph}{The input graph. It may be directed, but edge directions
    are ignored.}
}
\details{
  This function implements the Kanevsky algorithm for finding all
  minimal-size vertex separators in an undirected graph. See the
  reference below for the details.

  In the special case of a fully connected input graph with \eqn{n}
  vertices, all subsets of size \eqn{n-1} are listed as the result.
}
\value{
  A list of numeric vectors. Each numeric vector is a vertex separator.
}
\references{
  Arkady Kanevsky: Finding all minimum-size separating vertex sets in a
  graph. \emph{Networks} 23 533--541, 1993.

  JS Provan and DR Shier: A Paradigm for listing
  (s,t)-cuts in graphs, \emph{Algorithmica} 15, 351--372, 1996.

  J. Moody and D. R. White. Structural cohesion and embeddedness: A
  hierarchical concept of social groups. \emph{American Sociological
    Review}, 68 103--127, Feb 2003.  
}
\author{Gabor Csardi \email{csardi.gabor@gmail.com}}
\seealso{
  \code{\link{is.separator}}
}
\examples{
# The graph from the Moody-White paper
mw <- graph.formula(1-2:3:4:5:6, 2-3:4:5:7, 3-4:6:7, 4-5:6:7,
                    5-6:7:21, 6-7, 7-8:11:14:19, 8-9:11:14, 9-10,
                    10-12:13, 11-12:14, 12-16, 13-16, 14-15, 15-16,
                    17-18:19:20, 18-20:21, 19-20:22:23, 20-21,
                    21-22:23, 22-23)

# Cohesive subgraphs
mw1 <- induced.subgraph(mw, as.character(c(1:7, 17:23)))
mw2 <- induced.subgraph(mw, as.character(7:16))
mw3 <- induced.subgraph(mw, as.character(17:23))
mw4 <- induced.subgraph(mw, as.character(c(7,8,11,14)))
mw5 <- induced.subgraph(mw, as.character(1:7))

minimum.size.separators(mw)
minimum.size.separators(mw1)
minimum.size.separators(mw2)
minimum.size.separators(mw3)
minimum.size.separators(mw4)
minimum.size.separators(mw5)

# Another example, the science camp network
camp <- graph.formula(Harry:Steve:Don:Bert - Harry:Steve:Don:Bert,
                      Pam:Brazey:Carol:Pat - Pam:Brazey:Carol:Pat,
                      Holly   - Carol:Pat:Pam:Jennie:Bill,
                      Bill    - Pauline:Michael:Lee:Holly,
                      Pauline - Bill:Jennie:Ann,
                      Jennie  - Holly:Michael:Lee:Ann:Pauline,
                      Michael - Bill:Jennie:Ann:Lee:John,
                      Ann     - Michael:Jennie:Pauline,
                      Lee     - Michael:Bill:Jennie,
                      Gery    - Pat:Steve:Russ:John,
                      Russ    - Steve:Bert:Gery:John,
                      John    - Gery:Russ:Michael)
lapply(minimum.size.separators(camp), function(x) V(camp)[x])
}
\keyword{graphs}
