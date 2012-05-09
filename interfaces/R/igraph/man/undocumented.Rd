\name{igraph.undocumented}
\alias{cited.type.game}
\alias{citing.cited.type.game}
\alias{lastcit.game}
\title{Undocumented and unsupportted igraph functions}
\description{These functions are still in the alpha stage or
  their arguments are expected to change, so they're not documented yet.
  They are also not very useful for the general audience.
}
\usage{
lastcit.game(n, edges=1, agebins=n/7100, pref=(1:(agebins+1))^-3,
     directed=TRUE)
cited.type.game(n, edges=1, types=rep(0, n),
     pref=rep(1, length(types)),
     directed=TRUE, attr=TRUE)
citing.cited.type.game(n, edges=1, types=rep(0, n),
    pref=matrix(1, nrow=length(types), ncol=length(types)),
    directed=TRUE, attr=TRUE)
}
\arguments{
  \item{n}{Number of vertices.}
  \item{edges}{Number of edges per step.}
  \item{agebins}{Number of aging bins.}
  \item{pref}{Vector (\code{lastcit.game} and \code{cited.type.game} or
    matrix (\code{citing.cited.type.game}) giving the (unnormalized)
    citation probabilities for the different vertex types.} 
  \item{directed}{Logical scalar, whether to generate directed
    networks.}
  \item{types}{Vector of length \sQuote{\code{n}}, the types of the
    vertices. Types are numbered from zero.}
  \item{attr}{Logical scalar, whether to add the vertex types to the
    generated graph as a vertex attribute called \sQuote{\code{type}}. }
}
%\details{}
\value{A new graph.}
%\references{}
\author{Gabor Csardi \email{csardi@rmki.kfki.hu}}
%\seealso{}
%\examples{}
\keyword{graphs}
