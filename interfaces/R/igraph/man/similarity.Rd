\name{similarity}
\alias{similarity.jaccard}
\alias{similarity.dice}
\alias{similarity.invlogweighted}
\concept{Vertex similarity}
\title{Similarity measures of two vertices}
\description{
  These functions calculates similarity scores for vertices based on
  their connection patterns.
}
\usage{
similarity.jaccard(graph, vids = V(graph), mode = c("all", "out", "in",
      "total"), loops = FALSE) 
similarity.dice(graph, vids = V(graph), mode = c("all", "out", "in",
      "total"), loops = FALSE)
similarity.invlogweighted(graph, vids = V(graph),
       mode = c("all", "out", "in", "total"))
}
\arguments{
  \item{graph}{The input graph.}
  \item{vids}{The vertex ids for which the similarity is calculated.}
  \item{mode}{The type of neighboring vertices to use for the
    calculation, possible values: \sQuote{\code{out}},
    \sQuote{\code{in}}, \sQuote{\code{all}}.}
  \item{loops}{Whether to include vertices themselves in the neighbor
    sets.}
}
\details{
 The Jaccard similarity coefficient of two vertices is the number of common
 neighbors divided by the number of vertices that are neighbors of at
 least one of the two vertices being
 considered. \code{similarity.jaccard} calculates
 the pairwise Jaccard similarities for some (or all) of the vertices.

 The Dice similarity coefficient of two vertices is twice the number of common
 neighbors divided by the sum of the degrees of the
 vertices. \code{similarity.dice} calculates the pairwise Dice
 similarities for some (or all) of the vertices.

 The inverse log-weighted similarity of two vertices is the number of
 their common neighbors, weighted by the inverse logarithm of their degrees.
 It is based on the assumption that two vertices should be considered
 more similar if they share a low-degree common neighbor, since high-degree
 common neighbors are more likely to appear even by pure chance.
 Isolated vertices will have zero similarity to any other vertex.
 Self-similarities are not calculated.
 See the following paper for more details: Lada A. Adamic and Eytan Adar:
 Friends and neighbors on the Web. Social Networks, 25(3):211-230, 2003.
}
\value{
  A \code{length(vids)} by \code{length(vids)} numeric matrix containing
  the similarity scores.
}
\references{Lada A. Adamic and Eytan Adar:
 Friends and neighbors on the Web. \emph{Social Networks}, 25(3):211-230, 2003.}
\author{Tamas Nepusz \email{ntamas@rmki.kfki.hu} and Gabor Csardi
  \email{csardi@rmki.kfki.hu} for the manual page.}
\seealso{\code{\link{cocitation}} and \code{\link{bibcoupling}}}
\examples{
g <- graph.ring(5)
similarity.dice(g)
similarity.jaccard(g)
}
\keyword{graphs}
