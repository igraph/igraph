\name{running.mean}
\alias{running.mean}
\title{Running mean of a time series}
\description{\code{running.mean} calculates the running mean in a vector
  with the given bin width.}
\usage{
running.mean(v, binwidth)
}
\arguments{
  \item{v}{The numeric vector.}
  \item{binwidth}{Numeric constant, the size of the bin, should be
    meaningful, ie. smaller than the length of \code{v}. }
}
\details{The running mean of \code{v} is a \code{w}
  vector of length \code{length(v)-binwidth+1}. The first element of
  \code{w} id the average of the first \code{binwidth} elements of
  \code{v}, the second element of \code{w} is the average of elements
  \code{2:(binwidth+1)}, etc.}
\value{
  A numeric vector of length \code{length(v)-binwidth+1}
}
% \references{}
\author{Gabor Csardi \email{csardi@rmki.kfki.hu}}
% \seealso{}
\examples{
running.mean(1:100, 10)
}
\keyword{manip}
