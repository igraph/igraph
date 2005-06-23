\name{ba.game}
\alias{ba.game}
\alias{barabasi.game}
\title{Generate scale-free graphs according to the Barabasi-Albert
  model}
\description{The BA-model is a very simple stochastic algorithm for
  building a graph.}
\usage{
ba.game(n, m = NULL, out.dist = NULL, out.seq = NULL, out.pref = FALSE, 
        \dots) 
}
\arguments{
  \item{n}{Number of vertices.}
  \item{m}{Numeric constant, the contant out-degree of the vertices.}
  \item{out.dist}{Numeric vector, the distiribution of the out-degree.}
  \item{out.seq}{Numeric vector, the out-degree sequence of the vertices}
  \item{out.pref}{Logical, if true the total degree is used for
    calculating the citation probability, otherwise the in-degree is
    used. }
}
\details{
  The out-degree of the vertices is determined by the \code{m},
  \code{out.dist} and \code{out.seq} parameters, only at most one of
  these should be not NULL. If all of these are NULL, then constant
  out-degree one is used.
  
  If \code{m} is not NULL, it specifies the constant out-degree.

  Otherwise if \code{out.dist} is not NULL, it is used to generate an
  out-degree sequence, the maximum out-degree is
  \code{length(out.dist)-1} and the probability of having out-degree
  \eqn{k}{\code{k}} is \code{out.dist[k+1]/sum(out.dist)}.

  Otherwise if \code{out.seq} is not NULL it is used directly as the
  out-degrees of the vertices. Its length should be exactly \code{n}.
}
\value{A graph object.}
\references{  Barabási, A.-L. and Albert R. 1999. ``Emergence of scaling
  in random networks'' Science, 286 509--512.}
\author{Gabor Csardi \email{csardi@rmki.kfki.hu}}
\seealso{\code{\link{random.graph.game}}}
\examples{
g <- ba.game(1000)
degree.distribution(g)
}
\keyword{graphs}
