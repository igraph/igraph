\name{igraph.sample}
\alias{igraph.sample}
\title{Sampling a random integer sequence}
\description{This function provides a very efficient way to pull an
  integer random sample sequence from an integer interval.}
\usage{
igraph.sample(low, high, length)
}
\arguments{
  \item{low}{The lower limit of the interval (inclusive).}
  \item{high}{The higher limit of the interval (inclusive).}
  \item{length}{The length of the sample.}
}
\details{The algorithm runs in \code{O(length)} expected time, even if
  \code{high-low} is big. It is much faster (but of course less general)
  than the builtin \code{sample} function of R.}
\value{
  An increasing numeric vector containing integers, the sample.
}
\references{Jeffrey Scott Vitter: An Efficient Algorithm for Sequential
  Random Sampling, \emph{ACM Transactions on Mathematical Software}, 13/1,
  58--67.}
\author{Gabor Csardi \email{csardi@rmki.kfki.hu}}
\examples{
rs <- igraph.sample(1, 100000000, 10)
rs
}
\keyword{datagen}
    