\name{spinglass}
\alias{spinglass}
\concept{Statistical mechanics}
\concept{Spin-glass}
\concept{Community structure (roles)}
\title{Finding communities in graphs based on statistical meachanics, including negative weights.}
\description{This function tries to find communities in graphs via
  a spin-glass model and simulated annealing.}
\synopsis{
spinglass(graph, weights=NULL, spins=25,
                    parupdate=FALSE, start.temp=1, stop.temp=0.01,
                    cool.fact=0.99, update.rule=c("config", "random",
                    "simple"), gamma=1, lambda=1)
}
\usage{
spinglass(graph, weights=NULL, spins=25,
                    parupdate=FALSE, start.temp=1, stop.temp=0.01,
                    cool.fact=0.99, update.rule=c("config", "random",
                    "simple"), gamma=1, lambda=1)
}
\arguments{
  \item{graph}{The input graph, can be directed.}
  \item{weights}{The weights of the edges, which can be negative. Either a numeric vector or
    \code{NULL}. If it is null and the input graph has a \sQuote{weight}
    edge attribute then that will be used. If \code{NULL} and no such
    attribute is present then the edges will have equal weights.}
  \item{spins}{Integer constant, the number of spins to use. This is the
    upper limit for the number of communities. It is not a problem to
    supply a (reasonably) big number here, in which case some
    spin states will be unpopulated. }
  \item{parupdate}{Logical constant, whether to update the spins of the
    vertices in parallel (synchronously) or not. This argument is currently
	not used }
  \item{start.temp}{Real constant, the start temperature. } 
  \item{stop.temp}{Real constant, the stop temperature. } 
  \item{cool.fact}{Cooling factor for the simulated annealing. }     
  \item{update.rule}{Character constant giving the \sQuote{null-model}
    of the simulation. Possible values: \dQuote{simple} and
    \dQuote{config}. \dQuote{simple} uses a random graph with the same
    number of edges as the baseline probability and \dQuote{config} uses
    a random graph with the same vertex degrees as the input graph. This
	is currently not used.}
  \item{gamma}{Real constant, the gamma argument of the algorithm. This
    specifies the balance between the importance of present and
    non-present edges in a community. Larger values of gamma, leads to
	communities with an higher intraconnectivity.}
  \item{lambda}{Real constant, the ;ambda argument of the algorithm. This
    specifies the balance between the importance of present and
    non-present negative weighted edges in a community. Smaller values of 
	lambda, leads to communities with lesser negative intraconnectivity.
	If this argument is set to zero, the algorithm reduces to a graph
	coloring algorithm, using the number of spins as the number of colors.}	
}
\details{
  This function tries to find communities in a graph. A community is a
  set of nodes with many edges inside the community and  few edges
  between outside it (ie. between the community itself and the rest of
  the graph. This idea is reversed for edges having a negative weight, 
  ie. few negative edges inside a community and many negative edges
  between communities.)
}
\value{
  The function returns a named list with the following slots:
  \itemize{
    \item{membership}{Integer vector giving the communities found. The
      communities have ids starting from zero and for each graph vertex
      ids community id is given in this vector.}
    \item{csize}{The sizes of the communities in the order of their ids.}
    \item{modularity}{The modularity score of the result, as defined by
      Newman and Girvan, see references.}  	  
    \item{temperature}{The temperature of the system when the algorithm
      terminated.}
  }
}
\references{
  J. Reichardt and S. Bornholdt: Statistical Mechanics of Community
  Detection, \emph{Phys. Rev. E}, 74, 016110 (2006),
  \url{http://arxiv.org/abs/cond-mat/0603718}

  M. E. J. Newman and M. Girvan: Finding and evaluating community
  structure in networks, \emph{Phys. Rev. E} 69, 026113 (2004) 
  
  V.A. Traag and Jeroen Bruggeman: Community detection in networks 
  with positive and negative links, \url{http://arxiv.org/abs/0811.2329}  (2008).
}
\author{Jorg Reichardt \email{lastname@physik.uni-wuerzburg.de} for the
  original code and Gabor Csardi \email{csardi@rmki.kfki.hu} for the
  igraph glue code. Changes to the original function for including the
  possibility of negative ties (and directed networks) were implemented 
  by Vincent Traag \email{vtraag@f-m.fm}.}
\seealso{\code{\link{clusters}}}
\examples{
  g <- erdos.renyi.game(10, 5/10) \%du\% erdos.renyi.game(9, 5/9)
  g <- add.edges(g, c(0, 11))
  g <- subgraph(g, subcomponent(g, 0))
  spinglass(g, spins=2)
}
\keyword{graphs}

    