\name{igraph console}
\alias{igraph.console}
\concept{The igraph console}
\title{The igraph console}
\description{The igraph console is a GUI windows that shows what the
  currently running igraph function is doing.}
\usage{
igraph.console()
}
%\arguments{}
\details{
  The console can be started by calling the \code{igraph.console}
  function. Then it stays open, until the user closes it.

  Another way to start it to set the \code{verbose} igraph option to
  \dQuote{tkconsole} via \code{igraph.options}. Then the console
  (re)opens each time an igraph function supporting it starts; to close
  it, set the \code{verbose} option to another value.

  The console is written in Tcl/Tk and required the \code{tcltk}
  package.
}
\value{\code{NULL}, invisibly.}
% \references
\author{Gabor Csardi \email{csardi.gabor@gmail.com}}
\seealso{\code{\link{igraph.options}} and the \code{verbose} option.}
%\examples{}
\keyword{graphs}
