\name{igraph options}
\alias{igraph.options}
\alias{getIgraphOpt}
\alias{igraph.par}
\title{Parameters for the igraph package}
\description{igraph has some parameters which (usually) affect the
  behavior of many functions. These can be set for the whole session via
  \code{igraph.options}.
}
\usage{
igraph.options(\dots)
getIgraphOpt(x, default = NULL)
igraph.par(parid, parvalue = NULL)
}
\arguments{
  \item{\dots}{A list may be given as the only argument, or any number
    of arguments may be in the \code{name=value} form, or no argument at
    all may be given. See the Value and Details sections for explanation.}
  \item{x}{A character string holding an option name.}
  \item{default}{If the specified option is not set in the options list,
    this value is returned. This facilitates retrieving an option and
    checking whether it is set and setting it separately if not.}
  \item{parid}{The name of the parameter. See the currently used
    parameters below.}
  \item{parvalue}{The new value of the parameter. If \code{NULL} then
    the current value of the parameter is listed.}
}
\details{
  From igraph version 0.6, \code{igraph.par} is deprecated. Please use
  the more flexible \code{igraph.options} and \code{getIgraphOpt}
  functions instead.

  The parameter values set via a call to the \code{igraph.options}
  function will remain in effect for the rest of the session, affecting
  the subsequent behaviour of the other functions of the \code{igraph}
  package for which the given parameters are relevant.

  This offers the possibility of customizing the functioning of the
  \code{igraph} package, for instance by insertions of appropriate calls
  to \code{igraph.options} in a load hook for package \pkg{igraph}.
  
  The currently used parameters in alphabetical order:
  \describe{
    \item{edge.attr.comb:}{Specifies what to do with the edge
      attributes if the graph is modified. The default value is
      \code{list(weight="sum", name="concat", "ignore")}. See
      \code{\link{attribute.combination}} for details on this.}
    \item{print.edge.attributes:}{Logical constant, whether to print edge
      attributes when printing graphs. Defaults to \code{FALSE}.}
    \item{print.graph.attributes:}{Logical constant, whether to print
      graph attributes when printing graphs. Defaults to \code{FALSE}.}
    \item{print.vertex.attriubutes:}{Logical constant, whether to print
      vertex attributes when printing graphs. Defaults to \code{FALSE}.}
    \item{verbose:}{Logical constant, whether igraph functions should
      talk more than minimal. Eg. if \code{TRUE} thne some functions
      will use progress bars while computing. Defaults to \code{FALSE}.}
    \item{vertex.attr.comb:}{Specifies what to do with the vertex
      attributes if the graph is modified. The default value is
      \code{list(name="concat", "ignore")} See
      \code{\link{attribute.combination}} for details on this.}
  }
}
\value{
  \code{igraph.options} returns a list with  the updated values of the
  parameters. If the argument list is not empty, the returned list is
  invisible.

  For \code{getIgraphOpt}, the current value set for option \code{x}, or
  \code{NULL} if the option is unset.
  
  If \code{parvalue} is \code{NULL} then \code{igraph.par} returns the
  current value of the parameter. Otherwise the new value of the
  parameter is returned invisibly.
}
% \references{}
\author{ Gabor Csardi \email{csardi@rmki.kfki.hu}}
\seealso{\code{igraph.options} is similar to \code{\link{options}} and
  \code{getIgraphOpt} is similar to \code{\link{getOption}}.}
\examples{
oldval <- getIgraphOpt("verbose")
igraph.options(verbose=TRUE)
layout.kamada.kawai(graph.ring(10))
igraph.options(verbose=oldval)
}
\keyword{graphs}
