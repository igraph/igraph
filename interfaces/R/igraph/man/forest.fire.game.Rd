\name{forest.fire.game}
\alias{forest.fire.game}
\concept{Random graph model}
\concept{Forest fire model}
\title{Forest Fire Network Model}
\description{This is a growing network model, which resembles of how the
  forest fire spreads by igniting trees close by.}
\usage{
forest.fire.game (nodes, fw.prob, bw.factor = 1, ambs = 1, directed = TRUE, 
    verbose = igraph.par("verbose")) 
}
\arguments{
  \item{nodes}{The number of vertices in the graph.}
  \item{fw.prob}{The forward burning probability, see details below.}
  \item{bw.factor}{The backward burning ratio. The backward burning
    probability is calculated as \code{bw.factor*fw.prob}.}
  \item{ambs}{The number of ambassador vertices.}
  \item{directed}{Logical scalar, whether to create a directed graph.}
  \item{verbose}{Logical scalar, whether to \dQuote{draw} a progress bar.}
}
\details{
  The forest fire model intends to reproduce the following network
  characteristics, observed in real networks:
  \itemize{
    \item Heavy-tailed in-degree distribution. 
    \item Heavy-tailed out-degree distribution.
    \item Communities.
    \item Densification power-law. The network is densifying in time,
      according to a power-law rule.
    \item Shrinking diameter. The diameter of the network decreases in
      time.
  }

  The network is generated in the following way. One vertex is added at
  a time. This vertex connects to (cites) \code{ambs} vertices already
  present in the network, chosen uniformly random. Now, for each cited
  vertex \eqn{v} we do the following procedure:
  \enumerate{
    \item We generate two random number, \eqn{x} and \eqn{y}, that are
      geometrically distributed with means \eqn{p/(1-p)} and
      \eqn{rp(1-rp)}. (\eqn{p} is \code{fw.prob}, \eqn{r} is
      \code{bw.factor}.) The new vertex cites \eqn{x} outgoing neighbors
      and \eqn{y} incoming neighbors of \eqn{v}, from those which are
      not yet cited by the new vertex. If there are less than \eqn{x} or
      \eqn{y} such vertices available then we cite all of them.
    \item The same procedure is applied to all the newly cited
      vertices.
  }
}
\note{
  The version of the model in the published paper is incorrect
  in the sense that it cannot generate the kind of graphs the authors
  claim. A corrected version is available from
  \url{http://www.cs.cmu.edu/~jure/pubs/powergrowth-tkdd.pdf}, our
  implementation is based on this.
}
\value{A simple graph, possibly directed if the \code{directed} argument
  is \code{TRUE}.}
\references{
  Jure Leskovec, Jon Kleinberg and Christos Faloutsos. Graphs over time:
  densification laws, shrinking diameters and possible explanations.
  \emph{KDD '05: Proceeding of the eleventh ACM SIGKDD international
    conference on Knowledge discovery in data mining}, 177--187, 2005.
}
\author{
  Gabor Csardi \email{csardi@rmki.kfki.hu}
}
\seealso{
  \code{\link{barabasi.game}} for the basic preferential attachment
  model.
}
\examples{
g <- forest.fire.game(10000, fw.prob=0.37, bw.factor=0.32/0.37)
dd1 <- degree.distribution(g, mode="in")
dd2 <- degree.distribution(g, mode="out")
if (interactive()) {
  plot(seq(along=dd1)-1, dd1, log="xy")
  points(seq(along=dd2)-1, dd2, col=2, pch=2)
}
}
\keyword{graphs}
