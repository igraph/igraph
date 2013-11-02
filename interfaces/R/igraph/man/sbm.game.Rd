\name{sbm.game}
\alias{sbm.game}
\concept{Stochastic block model}
\concept{Random graph model}
\title{Sample stochastic block model}
\description{Sampling from the stochastic block model of networks}
\usage{
sbm.game (n, pref.matrix, block.sizes, directed = FALSE, loops = FALSE) 
}
\arguments{
  \item{n}{Number of vertices in the graph.}
  \item{pref.matrix}{The matrix giving the Bernoulli rates.
    This is a \eqn{K\times K}{KxK} matrix, where \eqn{K} is the number
    of groups. The probability of creating an edge between vertices from
    groups \eqn{i} and \eqn{j} is given by element \eqn{(i,j)}. For
    undirected graphs, this matrix must be symmetric.}
  \item{block.sizes}{Numeric vector giving the number of vertices in
    each group. The sum of the vector must match the number of vertices.}
  \item{directed}{Logical scalar, whether to generate a directed
    graph.}
  \item{loops}{Logical scalar, whether self-loops are allowed in the
    graph.}
}
\details{
  This function samples graphs from a stochastic block
  model by (doing the equivalent of) Bernoulli
  trials for each potential edge with the probabilities
  given by the Bernoulli rate matrix, \code{pref.matrix}.
}
\value{An igraph graph.}
\references{
  Faust, K., & Wasserman, S. (1992a). Blockmodels:
  Interpretation and evaluation. \emph{Social Networks}, 14, 5--61.
}
\author{Gabor Csardi \email{csardi.gabor@gmail.com}}
\seealso{\code{\link{random.graph.game}}}
\examples{
## Two groups with not only few connection between groups
pm <- cbind( c(.1, .001), c(.001, .05) )
g <- sbm.game(1000, pref.matrix=pm, block.sizes=c(300,700))
g
}
\keyword{graphs}
