\name{igraph options}
\alias{igraph.options}
\alias{igraph_options}
\alias{getIgraphOpt}
\alias{igraph_opt}
\title{Parameters for the igraph package}
\description{igraph has some parameters which (usually) affect the
  behavior of many functions. These can be set for the whole session via
  \code{igraph_options}.
}
\usage{
igraph_options(\dots)
igraph_opt(x, default = NULL)
}
\arguments{
  \item{\dots}{A list may be given as the only argument, or any number
    of arguments may be in the \code{name=value} form, or no argument at
    all may be given. See the Value and Details sections for explanation.}
  \item{x}{A character string holding an option name.}
  \item{default}{If the specified option is not set in the options list,
    this value is returned. This facilitates retrieving an option and
    checking whether it is set and setting it separately if not.}
}
\details{
  The parameter values set via a call to the \code{igraph_options}
  function will remain in effect for the rest of the session, affecting
  the subsequent behaviour of the other functions of the \code{igraph}
  package for which the given parameters are relevant.

  This offers the possibility of customizing the functioning of the
  \code{igraph} package, for instance by insertions of appropriate calls
  to \code{igraph_options} in a load hook for package \pkg{igraph}.
  
  The currently used parameters in alphabetical order:
  \describe{
    \item{add.params}{Logical scalar, whether to add model
      parameter to the graphs that are created by the various
      graph constructors. By default it is \code{TRUE}.}
    \item{add.vertex.names}{Logical scalar, whether to add
      vertex names to node level indices, like degree, betweenness
      scores, etc. By default it is \code{TRUE}.}
    \item{annotate.plot}{Logical scalar, whether to annotate igraph
      plots with the graph's name (\code{name} graph attribute, if
      present) as \code{main}, and with the number of vertices and edges
      as \code{xlab}. Defaults to \code{FALSE}.}
    \item{dend.plot.type}{The plotting function to use when plotting
      community structure dendrograms via
      \code{\link{plot_dendrogram}}}. Possible values are \sQuote{auto} (the
      default), \sQuote{phylo}, \sQuote{hclust} and
      \sQuote{dendrogram}. See \code{\link{plot_dendrogram}} for details.
    \item{edge.attr.comb}{Specifies what to do with the edge
      attributes if the graph is modified. The default value is
      \code{list(weight="sum", name="concat", "ignore")}. See
      \code{\link{attribute.combination}} for details on this.}
    \item{nexus.url}{The base URL of the default Nexus server. See
      \code{\link{nexus}} for details.}
    \item{print.edge.attributes}{Logical constant, whether to print edge
      attributes when printing graphs. Defaults to \code{FALSE}.}
    \item{print.full}{Logical scalar, whether \code{\link{print.igraph}}
      should show the graph structure as well, or only a summary of the
      graph.}
    \item{print.graph.attributes}{Logical constant, whether to print
      graph attributes when printing graphs. Defaults to \code{FALSE}.}
    \item{print.vertex.attributes}{Logical constant, whether to print
      vertex attributes when printing graphs. Defaults to \code{FALSE}.}
    \item{sparsematrices}{Whether to use the \code{Matrix} package for
      (sparse) matrices. It is recommended, if the user works with
      larger graphs.}
    \item{verbose}{Logical constant, whether igraph functions should
      talk more than minimal. Eg. if \code{TRUE} thne some functions
      will use progress bars while computing. Defaults to \code{FALSE}.}
    \item{vertex.attr.comb}{Specifies what to do with the vertex
      attributes if the graph is modified. The default value is
      \code{list(name="concat", "ignore")} See
      \code{\link{attribute.combination}} for details on this.}
  }
}
\value{
  \code{igraph_options} returns a list with  the updated values of the
  parameters. If the argument list is not empty, the returned list is
  invisible.

  For \code{igraph_opt}, the current value set for option \code{x}, or
  \code{NULL} if the option is unset.
}
% \references{}
\author{ Gabor Csardi \email{csardi.gabor@gmail.com}}
\seealso{\code{igraph_options} is similar to \code{\link{options}} and
  \code{igraph_opt} is similar to \code{\link{getOption}}.}
\examples{
oldval <- igraph_opt("verbose")
igraph_options(verbose=TRUE)
layout_with_kk(make_ring(10))
igraph_options(verbose=oldval)
}
\keyword{graphs}
