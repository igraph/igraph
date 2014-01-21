\name{srand}
\alias{srand}
\title{Set random seed of the C library's RNG}
\description{Set the random seed of the C library's RNG, for a new
  sequence of pseudo-random numbers.}
\usage{
srand(seed)
}
\arguments{
  \item{seed}{Numeric scalar, the new random seed. It must be
    non-negative and will be converted to an integer.}
}
\details{
  Note that this function has nothing to do with R's random number
  generator, see \code{set.seed} for that.

  Some package (e.g. ngspatial) use internal C code and generate random
  numbers using the standard C library's built-in random number
  generator instead of using R's RNGs. The \code{srand} function is
  provided to set the random seed for these packages. It simply calls
  the standard C function \code{srand}, with the supplied integer seed
  value.

  Note that the standard C library's RNGs are typically of very bad
  quality, and also slower than R's RNGs. It is not worth using
  them, really, other than taking over some legacy C code that already
  uses them, and that would be difficult to rewrite to use R's RNGs.
}
\value{
  \code{NULL}, invisibly.
}
% \references{}
\author{Gabor Csardi \email{csardi.gabor@gmail.com}}
% \seealso{}
% \examples{}
