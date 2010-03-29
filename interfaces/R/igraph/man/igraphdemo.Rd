\name{igraphdemo}
\alias{igraphdemo}
\title{Run igraph demos, step by step}
\description{Run one of the accompanying igraph demos, somewhat
  interactively, using a Tk window. }
\usage{
igraphdemo(which)
}
\arguments{
  \item{which}{If not given, then the names of the available demos are
    listed. Otherwise it should be either a filename or the name of an
    igraph demo.}
}
\details{
  This function provides a somewhat nicer interface to igraph demos that
  come with the package, than the standard \code{\link{demo}}
  function. Igraph demos are divided into chunks and \code{igraphdemo}
  runs them chunk by chunk, with the possibility of inspecting the
  workspace between two chunks.

  The \code{tcltk} package is needed for \code{igraphdemo}.
}
\value{Returns \code{NULL}, invisibly.}
%\references{}
\author{Gabor Csardi \email{Gabor.Csardi@unil.ch}}
\seealso{\code{\link{demo}}}
\examples{
igraphdemo()
if (interactive()) {
  igraphdemo("centrality")
}
}
\keyword{graphs}
