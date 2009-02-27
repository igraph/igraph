\name{dyad.census}
\alias{dyad.census}
\concept{Dyad census}
\title{Dyad census of a graph}
\description{Classify dyads in a directed graphs. The relationship
  between each pair of vertices is measured. It can be in three states:
  mutual, asymmetric or non-existent.}
\usage{
dyad.census(graph)
}
\arguments{
  \item{graph}{The input graph. A warning is given if it is not directed.}
}
%\details{}
\value{A named numeric vector with three elements:
  \item{mut}{The number of pairs with mutual connections.}
  \item{asym}{The number of pairs with non-mutual connections.}
  \item{null}{The number of pairs with no connection between them.}
}
\references{
  Holland, P.W. and Leinhardt, S. A Method for Detecting
  Structure in Sociometric Data. \emph{American Journal of Sociology}, 76,
  492--513. 1970.

  Wasserman, S., and Faust, K. \emph{Social Network Analysis: Methods
  and Applications.} Cambridge: Cambridge University Press. 1994.
}
\author{
  Gabor Csardi \email{csardi@rmki.kfki.hu}
}
\seealso{
  \code{\link{triad.census}} for the same classification, but with triples.
}
\examples{
g <- ba.game(100)
dyad.census(g)
}
\keyword{graphs}
