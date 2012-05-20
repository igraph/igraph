\name{assortativity}
\alias{assortativity}
\alias{assortativity.degree}
\alias{assortativity.nominal}
\concept{Assortativity coefficient}
\title{Assortativity coefficient}
\description{
  The assortativity coefficient is positive is similar vertices (based
  on some external property) tend to connect to each, and negative
  otherwise.
}
\usage{
assortativity (graph, types1, types2 = NULL, directed = TRUE) 
assortativity.nominal (graph, types, directed = TRUE) 
assortativity.degree (graph, directed = TRUE) 
}
\arguments{
  \item{graph}{The input graph, it can be directed or undirected.}
  \item{types}{Vector giving the vertex types. They as assumed to be
    integer numbers, starting with one. Non-integer values are
    converted to integers with \code{\link{as.integer}}.}
  \item{types1}{The vertex values, these can be arbitrary numeric
    values.}
  \item{types2}{A second value vector to be using for the incoming
    edges when calculating assortativity for a directed graph. 
    Supply \code{NULL} here if you want to use the same values
    for outgoing and incoming edges. This argument is ignored
    (with a warning) if it is not \code{NULL} and undirected
    assortativity coefficient is being calculated.}
  \item{directed}{Logical scalar, whether to consider edge directions
    for directed graphs. This argument is ignored for undirected
    graphs. Supply \code{TRUE} here to do the natural thing, i.e. use
    directed version of the measure for directed graphs and the
    undirected version for undirected graphs.}
}
\details{
  The assortativity coefficient measures the level of homophyly of the
  graph, based on some vertex labeling or values assigned to
  vertices. If the coefficient is high, that means that connected
  vertices tend to have the same labels or similar assigned values.
  
  M.E.J. Newman defined two kinds of assortativity coefficients, the
  first one is for categorical labels of
  vertices. \code{assortativity.nominal} calculates this measure. It is
  defines as

  \deqn{r=\frac{\sum_i e_{ii}-\sum_i a_i b_i}{1-\sum_i a_i b_i}}{
    r=(sum(e(i,i), i) - sum(a(i)b(i), i)) / (1 - sum(a(i)b(i), i))}

  where \eqn{e_{ij}}{e(i,j)} is the fraction of edges connecting vertices
  of type \eqn{i} and \eqn{j},
  \eqn{a_i=\sum_j e_{ij}}{a(i)=sum(e(i,j), j)} and
  \eqn{b_j=\sum_i e_{ij}}{b(j)=sum(e(i,j), i)}.
  
  The second assortativity variant is based on values assigned to the
  vertices. \code{assortativity} calculates this measure. It is defined
  as

  \deqn{r=\frac1{\sigma_q^2}\sum_{jk} jk(e_{jk}-q_j q_k)}{
    sum(jk(e(j,k)-q(j)q(k)), j, k) / sigma(q)^2}

  for undirected graphs (\eqn{q_i=\sum_j e_{ij}}{q(i)=sum(e(i,j), j)})
  and as 

  \deqn{r=\frac1{\sigma_o\sigma_i}\sum_{jk}jk(e_{jk}-q_j^o q_k^i)}{
    sum(jk(e(j,k)-qout(j)qin(k)), j, k) / sigma(qin) / sigma(qout)
  }

  for directed ones. Here
  \eqn{q_i^o=\sum_j e_{ij}}{qout(i)=sum(e(i,j), j)},
  \eqn{q_i^i=\sum_j e_{ji}}{qin(i)=sum(e(j,i), j)},
  moreover, \eqn{\sigma_q}{sigma(q)},
  \eqn{sigma_o}{sigma(qout)} and
  \eqn{sigma_i}{sigma(qin)} are the standard deviations of \eqn{q},
  \eqn{q^o}{qout} and \eqn{q^i}{qin}, respectively.

  The reason of the difference is that in directed
  networks the relationship is not symmetric, so it is possible to
  assign different values to the outgoing and the incoming end of the
  edges.

  \code{assortativity.degree} uses vertex degree (minus one) as vertex
  values and calls \code{assortativity}.
}
\value{
  A single real number.
}
\references{
 M. E. J. Newman: Mixing patterns in networks, \emph{Phys. Rev. E} 67,
 026126 (2003) \url{http://arxiv.org/abs/cond-mat/0209450}

 M. E. J. Newman: Assortative mixing in networks,
 \emph{Phys. Rev. Lett.} 89, 208701 (2002)
 \url{http://arxiv.org/abs/cond-mat/0205405/}
}
\author{
  Gabor Csardi \email{csardi.gabor@gmail.com}
}
% \seealso{}
\examples{
# random network, close to zero
assortativity.degree(erdos.renyi.game(10000,3/10000))

# BA model, tends to be dissortative
assortativity.degree(ba.game(10000, m=4))
}
\keyword{graphs}
