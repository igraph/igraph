\name{barabasi.game}
\alias{barabasi.game}
\alias{ba.game}
\concept{Preferential attachment model}
\concept{Random graph model}
\title{Generate scale-free graphs according to the Barabasi-Albert
  model}
\description{The BA-model is a very simple stochastic algorithm for
  building a graph.}
\usage{
barabasi.game(n, power = 1, m = NULL, out.dist = NULL, out.seq = NULL, 
    out.pref = FALSE, zero.appeal = 1, directed = TRUE, time.window = NULL) 
}
\arguments{
  \item{n}{Number of vertices.}
  \item{power}{The power of the preferential attachment, the default is
    one, ie. linear preferential attachment.}
  \item{m}{Numeric constant, the number of edges to add in each time
    step This argument is only used if both \code{out.dist} and
    \code{out.seq} are omitted or NULL.}
  \item{out.dist}{Numeric vector, the distribution of the number of
    edges to add in each time step. This argument is only used if the
    \code{out.seq} argument is omitted or NULL.}
  \item{out.seq}{Numeric vector giving the number of edges to add in
    each time step. Its first element is ignored as no edges are added
    in the first time step.}
  \item{out.pref}{Logical, if true the total degree is used for
    calculating the citation probability, otherwise the in-degree is
    used. }
  \item{zero.appeal}{The \sQuote{attractiveness} of the vertices with no
    adjacent edges. See details below.}
  \item{directed}{Whether to create a directed graph.}
  \item{time.window}{Integer constant, if given only edges added in the
    previous \code{time.window} time steps are counted as the basis of
    preferential attachment.}
}
\details{
  This is a simple stochastic algorithm to generate a graph. It is a
  discrete time step model and in each time step a single vertex is
  added.

  We start with a single vertex and no edges in the first time
  step. Then we add one vertex in each time step and the new vertex
  initiates some edges to old vertices. The probability that an old vertex
  is chosen is given by
  \deqn{P[i] \sim k_i^\alpha+a}{P[i] ~ k[i]^alpha + a}
  where \eqn{k_i}{k[i]} is the in-degree of vertex \eqn{i} in the
  current time step (more precisely the number of adjacent edges of
  \eqn{i} which were not initiated by \eqn{i} itself) and
  \eqn{\alpha}{alpha} and \eqn{a} are parameters given by the
  \code{power} and \code{zero.appeal} arguments.

  The number of edges initiated in a time step is given by the \code{m},
  \code{out.dist} and \code{out.seq} arguments. If \code{out.seq} is
  given and not NULL then it gives the number of edges to add in a
  vector, the first element is ignored, the second is the number of
  edges to add in the second time step and so on. If \code{out.seq} is
  not given or null and \code{out.dist} is given and not NULL then it is
  used as a discrete distribution to generate the number of edges in
  each time step. Its first element is the probability that no edges
  will be added, the second is the probability that one edge is added,
  etc. (\code{out.dist} does not need to sum up to one, it normalized
  automatically.) \code{out.dist} should contain non-negative numbers
  and at east one element should be positive.

  If both \code{out.seq} and \code{out.dist} are omitted or NULL then
  \code{m} will be used, it should be a positive integer constant and
  \code{m} edges will be added in each time step.

  \code{barabasi.game} generates a directed graph by default, set
  \code{directed} to \code{FALSE} to generate an undirected graph. Note
  that even if an undirected graph is generated \eqn{k_i}{k[i]} denotes
  the number of adjacent edges not initiated by the vertex itself and
  not the total (in- + out-) degree of the vertex, unless the
  \code{out.pref} argument is set to \code{TRUE}.

  If the \code{time.window} argument is not NULL then \eqn{k_i}{k[i]} is
  the number of adjacent edges of \eqn{i} added in the previous
  \code{time.window} time steps.

  Note that \code{barabasi.game} might generate graphs with multiple
  edges.
}
\value{A graph object.}
\references{  Barabasi, A.-L. and Albert R. 1999. Emergence of scaling
  in random networks \emph{Science}, 286 509--512.}
\author{Gabor Csardi \email{csardi@rmki.kfki.hu}}
\seealso{\code{\link{random.graph.game}}}
\examples{
g <- barabasi.game(10000)
degree.distribution(g)
}
\keyword{graphs}
