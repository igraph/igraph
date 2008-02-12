\name{aging.prefatt.game}
\alias{aging.prefatt.game}
\alias{aging.barabasi.game}
\alias{aging.ba.game}
\concept{Preferential attachment}
\concept{Aging of vertices}
\concept{Random graph model}
\title{Generate an evolving random graph with preferential attachment
  and aging}
\description{This function creates a random graph by simulating its
  evolution. Each time a new vertex is added it creates a number of links
  to old vertices and the probability that an old vertex is cited depends on
  its in-degree (preferential attachment) and age.}
\usage{
aging.prefatt.game (n, pa.exp, aging.exp, m = NULL, aging.bin = 300,
     out.dist = NULL, out.seq = NULL, out.pref = FALSE,
     directed = TRUE, zero.deg.appeal = 1, zero.age.appeal = 0,
     deg.coef = 1, age.coef = 1, time.window = NULL)
}
\arguments{
  \item{n}{The number of vertices in the graph.}
  \item{pa.exp}{The preferantial attachment exponent, see the details
    below.}
  \item{aging.exp}{The exponent of the aging, usually a non-positive
    number, see details below.}
  \item{m}{The number of edges each new vertex creates (except the very
    first vertex). This argument is used only if both the
    \code{out.dist} and \code{out.seq} arguments are NULL.}
  \item{aging.bin}{The number of bins to use for measuring the age of
    vertices, see details below.}
  \item{out.dist}{The discrete distribution to generate the number of
    edges to add in each time step if \code{out.seq} is NULL. See
    details below.}
  \item{out.seq}{The number of edges to add in each time step, a vector
    containing as many elements as the number of vertices. See details
    below.}
  \item{out.pref}{Logical constant, whether to include edges not
    initiated by the vertex as a basis of preferential attachment. See
    details below.}
  \item{directed}{Logical constant, whether to generate a directed
    graph. See details below.}
  \item{zero.deg.appeal}{The degree-dependent part of the
    \sQuote{attractiveness} of the vertices with no adjacent edges. See
    also details below.} 
  \item{zero.age.appeal}{The age-dependent part of the
    \sQuote{attrativeness} of the vertices with age zero. It is usually
    zero, see details below.} 
  \item{deg.coef}{The coefficient of the degree-dependent
    \sQuote{attractiveness}. See details below.}
  \item{age.coef}{The coefficient of the age-dependent part of the
    \sQuote{attractiveness}. See details below.}
  \item{time.window}{Integer constant, if NULL only adjacent added in
    the last \code{time.windows} time steps are counted as a basis of
    the preferential attachment. See also details below.}
}
\details{
  This is a discrete time step model of a growing graph. We start with a
  network containing a single vertex (and no edges) in the first time
  step. Then in each time step (starting with the second) a new vertex
  is added and it initiates a number of edges to the old vertices in the
  network. The probability that an old vertex is connected to is
  proportional to 
  \deqn{P[i] \sim (c\cdot k_i^\alpha+a)(d\cdot l_i^\beta+b)\cdot }{%
        P[i] ~ (c k[i]^alpha + a) (d l[i]^beta + a)}

  Here \eqn{k_i}{k[i]} is the in-degree of vertex \eqn{i} in the current
  time step and \eqn{l_i}{l[i]} is the age of vertex \eqn{i}. The age is
  simply defined as the number of time steps passed since the vertex is
  added, with the extension that vertex age is divided to be in
  \code{aging.bin} bins.

 \eqn{c}, \eqn{\alpha}{alpha}, \eqn{a}, \eqn{d}, \eqn{\beta}{beta} and
 \eqn{b} are parameters and they can be set via the following
 arguments: \code{pa.exp} (\eqn{\alpha}{alpha}, mandatory argument),
 \code{aging.exp} (\eqn{\beta}{beta}, mandatory argument),
 \code{zero.deg.appeal} (\eqn{a}, optional, the default value is 1),
 \code{zero.age.appeal} (\eqn{b}, optional, the default is 0),
 \code{deg.coef} (\eqn{c}, optional, the default is 1),
 and \code{age.coef} (\eqn{d}, optional, the default is 1).

 The number of edges initiated in each time step is governed by the
 \code{m}, \code{out.seq} and \code{out.pref} parameters. If
 \code{out.seq} is given then it is interpreted as a vector giving the
 number of edges to be added in each time step. It should be of length
 \code{n} (the number of vertices), and its first element will be
 ignored. If \code{out.seq} is not given (or NULL) and \code{out.dist}
 is given then it will be used as a discrete probability distribution to
 generate the number of edges. Its first element gives the probability
 that zero edges are added at a time step, the second element is the
 probability that one edge is added, etc. (\code{out.seq} should contain
 non-negative numbers, but if they don't sum up to 1, they will be
 normalized to sum up to 1. This behavior is similar to the \code{prob}
 argument of the \code{sample} command.)

 By default a directed graph is generated, but it \code{directed} is set
 to \code{FALSE} then an undirected is created. Even if an undirected
 graph is generaed \eqn{k_i}{k[i]} denotes only the adjacent edges not
 initiated by the vertex itself except if \code{out.pref} is set to
 \code{TRUE}.

 If the \code{time.window} argument is given (and not NULL) then
 \eqn{k_i}{k[i]} means only the adjacent edges added in the previous
 \code{time.window} time steps.

 This function might generate graphs with multiple edges.
}
\value{A new graph.}
% \references{}
\author{Gabor Csardi \email{csardi@rmki.kfki.hu}}
\seealso{\code{\link{barabasi.game}}, \code{\link{erdos.renyi.game}}}
\examples{
# The maximum degree for graph with different aging exponents
g1 <- aging.prefatt.game(10000, pa.exp=1, aging.exp=0, aging.bin=1000)
g2 <- aging.prefatt.game(10000, pa.exp=1, aging.exp=-1,   aging.bin=1000)
g3 <- aging.prefatt.game(10000, pa.exp=1, aging.exp=-3,   aging.bin=1000)
max(degree(g1))
max(degree(g2))
max(degree(g3))
}
\keyword{graphs}
