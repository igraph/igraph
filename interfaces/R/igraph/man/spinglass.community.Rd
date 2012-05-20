\name{spinglass.community}
\alias{spinglass.community}
\concept{Statistical mechanics}
\concept{Spin-glass}
\concept{Community structure}
\title{Finding communities in graphs based on statistical meachanics}
\description{This function tries to find communities in graphs via
  a spin-glass model and simulated annealing.}
\usage{
spinglass.community(graph, weights=NULL, vertex=NULL, spins=25,
                    parupdate=FALSE, start.temp=1, stop.temp=0.01,
                    cool.fact=0.99, update.rule=c("config", "random",
                    "simple"), gamma=1, implementation=c("orig", "neg"),
                    gamma.minus=1)
}
\arguments{
  \item{graph}{The input graph, can be directed but the direction of the
    edges is neglected.}
  \item{weights}{The weights of the edges. Either a numeric vector or
    \code{NULL}. If it is null and the input graph has a \sQuote{weight}
    edge attribute then that will be used. If \code{NULL} and no such
    attribute is present then the edges will have equal weights. Set
    this to \code{NA} if the graph was a \sQuote{weight} edge attribute,
    but you don't want to use it for community detection.}
  \item{vertex}{This parameter can be used to calculate the community of
      a given vertex without calculating all communities. Note that if
      this argument is present then some other arguments are ignored.}
  \item{spins}{Integer constant, the number of spins to use. This is the
    upper limit for the number of communities. It is not a problem to
    supply a (reasonably) big number here, in which case some
    spin states will be unpopulated. }
  \item{parupdate}{Logical constant, whether to update the spins of the
    vertices in parallel (synchronously) or not. This argument is
    ignored if the second form of the function is used (ie. the
    \sQuote{\code{vertex}} argument is present). }
  \item{start.temp}{Real constant, the start temperature.
    This argument is ignored if the second form of the
    function is used (ie. the \sQuote{\code{vertex}} argument is
    present). } 
  \item{stop.temp}{Real constant, the stop temperature. The simulation
    terminates if the temperature lowers below this level.
    This argument is ignored if the second form of the
    function is used (ie. the \sQuote{\code{vertex}} argument is
    present). } 
  \item{cool.fact}{Cooling factor for the simulated annealing.
    This argument is ignored if the second form of the
    function is used (ie. the \sQuote{\code{vertex}} argument is
    present). }     
  \item{update.rule}{Character constant giving the \sQuote{null-model}
    of the simulation. Possible values: \dQuote{simple} and
    \dQuote{config}. \dQuote{simple} uses a random graph with the same
    number of edges as the baseline probability and \dQuote{config} uses
    a random graph with the same vertex degrees as the input graph.}
  \item{gamma}{Real constant, the gamma argument of the algorithm. This
    specifies the balance between the importance of present and
    non-present edges in a community. Roughly, a comunity is a set of
    vertices having many edges inside the community and few edges
    outside the community. The default 1.0 value makes existing and
    non-existing links equally important. Smaller values make the
    existing links, greater values the missing links more important.}
  \item{implementation}{Character scalar. Currently igraph contains two
    implementations for the Spin-glass community finding algorithm. The
    faster original implementation is the default. The other
    implementation, that takes into account negative weights, can be
    chosen by supplying \sQuote{neg} here.}
  \item{gamma.minus}{Real constant, the gamma.minus parameter of
    the algorithm. This specifies the balance between the importance of
    present and non-present negative weighted edges in a
    community. Smaller values of gamma.minus, leads to communities with
    lesser negative intra-connectivity. If this argument is set to zero,
    the algorithm reduces to a graph coloring algorithm, using the
    number of spins as the number of colors. This argument is ignored if
    the \sQuote{orig} implementation is chosen.}
}
\details{
  This function tries to find communities in a graph. A community is a
  set of nodes with many edges inside the community and  few edges
  between outside it (i.e. between the community itself and the rest of
  the graph.)

  This idea is reversed for edges having a negative weight, 
  ie. few negative edges inside a community and many negative edges
  between communities. Note that only the \sQuote{neg} implementation
  supports negative edge weights.

  The \code{spinglass.cummunity} function can solve two problems related
  to community detection. If the \code{vertex} argument is not given (or
  it is \code{NULL}), then the regular community detection problem is
  solved (approximately), i.e. partitioning the vertices into
  communities, by optimizing the an energy function.

  If the \code{vertex} argument is given and it is not \code{NULL}, then
  it must be a vertex id, and the same energy function is used to find
  the community of the the given vertex. See also the examples below.
}
\value{
  If the \code{vertex} argument is not given, ie. the first form is used
  then a \code{\link{spinglass.community}} returns a
  \code{\link{communities}} object. 

  If the \code{vertex} argument is present, ie. the second form is used
  then a named list is returned with the following components:
  \item{community}{Numeric vector giving the ids of the vertices in
    the same community as \code{vertex}.}
  \item{cohesion}{The cohesion score of the result, see references.}
  \item{adhesion}{The adhesion score of the result, see references.}
  \item{inner.links}{The number of edges within the community of
    \code{vertex}.}
  \item{outer.links}{The number of edges between the community of
    \code{vertex} and the rest of the graph. }
}
\references{
  J. Reichardt and S. Bornholdt: Statistical Mechanics of Community
  Detection, \emph{Phys. Rev. E}, 74, 016110 (2006),
  \url{http://arxiv.org/abs/cond-mat/0603718}

  M. E. J. Newman and M. Girvan: Finding and evaluating community
  structure in networks, \emph{Phys. Rev. E} 69, 026113 (2004)

  V.A. Traag and Jeroen Bruggeman: Community detection in networks 
  with positive and negative links, \url{http://arxiv.org/abs/0811.2329}
  (2008).
  
}
\author{Jorg Reichardt \email{lastname@physik.uni-wuerzburg.de} for the
  original code and Gabor Csardi \email{csardi.gabor@gmail.com} for the
  igraph glue code.
  
  Changes to the original function for including the
  possibility of negative ties were implemented 
  by Vincent Traag \email{vtraag@f-m.fm}.}
\seealso{\code{\link{communities}}, \code{\link{clusters}}}
\examples{
  g <- erdos.renyi.game(10, 5/10) \%du\% erdos.renyi.game(9, 5/9)
  g <- add.edges(g, c(1, 12))
  g <- induced.subgraph(g, subcomponent(g, 1))
  spinglass.community(g, spins=2)
  spinglass.community(g, vertex=1)
}
\keyword{graphs}

    