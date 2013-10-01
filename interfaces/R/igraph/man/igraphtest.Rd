\name{igraphtest}
\alias{igraphtest}
\title{Run package tests}
\description{Runs all package tests.}
\usage{
igraphtest()
}
\details{
	The \code{testthat} package is needed to run all tests. The location
	tests themselves can be extracted from the package via
	\code{system.file("tests", package="igraph")}.

	This function simply calls the \code{test_dir} function from the
	\code{testthat} package on the test directory.
}
\value{Whatever is returned by \code{test_dir} from the \code{testthat}
	package.}
% \references{}
\author{ Gabor Csardi \email{csardi.gabor@gmail.com} }
% \seealso{}
% \examples{}
\keyword{graphs}
