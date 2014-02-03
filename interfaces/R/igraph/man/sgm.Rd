\name{sgm}
\alias{sgm}
\concept{Seeded Graph Matching}
\title{Match Graphs given a seeding of vertex correspondences}
\description{The approximate graph matching problem is to find a
  bijection between the vertices of two graphs, such that the number of
  edge disagreements between the corresponding vertex pairs are
  minimized.}
\usage{
sgm(A, B, m, start, iteration)
}
\arguments{
  \item{A}{TODO}
  \item{B}{TODO}
  \item{m}{TODO}
  \item{start}{TODO}
  \item{iteration}{TODO}
}
\details{
  The approximate graph matching problem is to find a
  bijection between the vertices of two graphs, such that the number of
  edge disagreements between the corresponding vertex pairs are
  minimized.

  Given a part of the bijection that consist of known
  correspondences (the seeds), \code{sgm} completes the bijection by
  computing the permutation matrix that permutes the rows and columns of
  the adjacency matrix.
}
\value{TODO}
\references{TODO}
\author{Vince Lyzinski \url{http://www.ams.jhu.edu/~lyzinski/}}
\seealso{TODO}
\examples{TODO}
\keyword{graphs}
