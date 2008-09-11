\name{degree.sequence.game}
\alias{degree.sequence.game}
\concept{Degree sequence}
\concept{Configuration model}
\title{Generate random graphs with a given degree sequence}
\description{It is often useful to create a graph with given vertex
  degrees. This is exactly what \code{degree.sequence.game} does.}
\usage{
degree.sequence.game(out.deg, in.deg = numeric(0),
     method = c("simple", "vl"), \dots) 
}
\arguments{
  \item{out.deg}{Numeric vector, the sequence of degrees (for undirected
  graphs) or out-degrees (for directed graphs). For undirected graphs
  its sum should be even. For directed graphs its sum should be the same
  as the sum of \code{in.deg}.}
  \item{in.deg}{For directed graph, the in-degree sequence.}
  \item{method}{Character, the method for generating the graph. Right
    now the \dQuote{simple} and \dQuote{vl} methods are implemented.}
  \item{\dots}{Additional arguments, these are used as graph
    attributes.}
}
\details{The \dQuote{simple} method connects the out-stubs of the
  edges (undirected graphs) or the out-stubs and in-stubs
  (directed graphs) together. This way loop edges and also multiple
  edges may be generated.

  This method is not adequate if one needs to
  generate simple graphs with a given degree sequence. The multiple and
  loop edges can be deleted, but then the degree sequence is distorted
  and there is nothing to ensure that the graphs are sampled uniformly.

  THe \dQuote{vl} method is a more sophisticated generator. The
  algorithm and the implementation was done by Fabien Viger and Matthieu
  Latapy. This generator always generates undirected, connected simple
  graphs, it is an error to pass the \code{in.deg} argument to it.
  The algorithm relies on first creating an initial (possibly
  unconnected) simple undirected graph with the given degree sequence
  (if this is possible at all). Then some rewiring is done to make the
  graph connected. Finally a Monte-Carlo algorithm is used to randomize
  the graph. The \dQuote{vl} samples from the undirected, connected
  simple graphs unformly. See
  \url{http://www-rp.lip6.fr/~latapy/FV/generation.html} for details.
}
\value{The new graph object.}
%\references{}
\author{Gabor Csardi \email{csardi@rmki.kfki.hu}}
\seealso{\code{\link{erdos.renyi.game}}, \code{\link{barabasi.game}},
  \code{\link{simplify}} to get rid of the multiple and/or loops edges.}
\examples{
## The simple generator
g <- degree.sequence.game(rep(2,100))
degree(g)
is.simple(g)   # sometimes TRUE, but can be FALSE
g2 <- degree.sequence.game(1:10, 10:1)
degree(g2, mode="out")
degree(g2, mode="in")

## The vl generator
g3 <- degree.sequence.game(rep(2,100), method="vl")
degree(g3)
is.simple(g3)  # always TRUE

## Exponential degree distribution
## Note, that we correct the degree sequence if its sum is odd
degs <- sample(1:100, 100, replace=TRUE, prob=exp(-0.5*(1:100)))
if (sum(degs) \%\% 2 != 0) { degs[1] <- degs[1] + 1 }
g4 <- degree.sequence.game(degs, method="vl")
all(degree(g4) == degs)

## Power-law degree distribution
## Note, that we correct the degree sequence if its sum is odd
degs <- sample(1:100, 100, replace=TRUE, prob=(1:100)^-2)
if (sum(degs) \%\% 2 != 0) { degs[1] <- degs[1] + 1 }
g5 <- degree.sequence.game(degs, method="vl")
all(degree(g5) == degs)
}
\keyword{graphs}
