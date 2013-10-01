\name{Hierarchical random graphs}
\alias{HRG}
\alias{hrg.consensus}
\alias{hrg.game}
\alias{hrg.create}
\alias{hrg.predict}
\alias{hrg.dendrogram}
\alias{hrg.fit}
\alias{print.igraphHRG}
\alias{print.igraphHRGConsensus}
\concept{Hierarchical random graphs}
\title{Hierarchical random graphs}
\description{Fitting and sampling hierarchical random graph models.}
\usage{
hrg.fit (graph, hrg = NULL, start = FALSE, steps = 0)
hrg.consensus (graph, hrg = NULL, start = FALSE, num.samples = 10000)

hrg.create (graph, prob)
hrg.dendrogram (hrg)

hrg.game (hrg)

hrg.predict (graph, hrg = NULL, start = FALSE, num.samples = 10000,
             num.bins = 25)

\method{print}{igraphHRG}(x, type=c("auto", "tree", "plain"),
           level = 3, \dots)
\method{print}{igraphHRGConsensus}(x, \dots)
}
\arguments{
  \item{graph}{The graph to fit the model to. Edge directions are
    ignored in directed graphs.}
  \item{hrg}{A hierarchical random graph model, in the form of an
    \code{igraphHRG} object. \code{hrg.fit} allows this to be
    \code{NULL}, in which case a random starting point is used for the
    fitting. \code{hrg.consensus} and \code{hrg.predict} allow this to
    be \code{NULL} as well, then a HRG is fitted to the graph first,
    from a random starting point.}
  \item{start}{Logical, whether to start the fitting/sampling from the
    supplied \code{igraphHRG} object, or from a random starting point.}
  \item{steps}{The number of MCMC steps to make. If this is zero, then
    the MCMC procedure is performed until convergence.}
  \item{num.samples}{Number of samples to use for consensus generation
    or missing edge prediction.}
  \item{prob}{A vector of probabilities, one for each vertex, in the
    order of vertex ids.}
  \item{num.bins}{Number of bins for the edge probabilities. Give a
    higher number for a more accurate prediction.}
  \item{x}{\code{igraphHRG} or \code{igraphHRGConsensus} object to
    print.}
  \item{type}{How to print the dendrogram, see details below.}
  \item{level}{The number of top levels to print from the dendrogram.}
  \item{\dots}{Additional arguments, not used currently.}
}
\details{
 A hierarchical random graph is an ensemble of undirected
 graphs with \eqn{n} vertices. It is defined via a binary tree with
 \eqn{n} leaf and \eqn{n-1} internal vertices, where the
 internal vertices are labeled with probabilities.
 The probability that two vertices are connected in the random graph
 is given by the probability label at their closest common
 ancestor.

 Please see references below for more about hierarchical random graphs. 
 
 igraph contains functions for fitting HRG models to a given network
 (\code{hrg.fit}, for generating networks from a given HRG
 ensemble (\code{hrg.game}), converting an igraph graph to a
 HRG and back (\code{hrg.create}, \code{hrg.dendrogram}),
 for calculating a consensus tree from a set of sampled HRGs
 (\code{hrg.consensus}) and for predicting missing edges in a
 network based on its HRG models (\code{hrg.predict}). 
 
 The igraph HRG implementation is heavily based on the code
 published by Aaron Clauset, at his website,
 \url{http://tuvalu.santafe.edu/~aaronc/hierarchy/}.

 \code{hrg.fit} fits a HRG to a given graph. It takes the specified
 \code{steps} number of MCMC steps to perform the fitting, or a
 convergence criteria if the specified number of steps is
 zero. \code{hrg.fit} can start from a given HRG, if this is given in
 the \code{hrg} argument and the \code{start} argument is \code{TRUE}.

 \code{hrg.consensus} creates a consensus tree from several fitted
 hierarchical random graph models, using phylogeny methods. If the
 \code{hrg} argument is given and \code{start} is set to \code{TRUE},
 then it starts sampling from the given HRG. Otherwise it optimizes the
 HRG log-likelihood first, and then samples starting from the optimum.

 \code{hrg.create} creates a HRG from an igraph graph. The igraph graph
 must be a directed binary tree, with \eqn{n-1} internal and \eqn{n}
 leaf vertices. The \code{prob} argument contains the HRG probability
 labels for each vertex; these are ignored for leaf vertices.

 \code{hrg.dendrogram} creates the corresponsing igraph tree of a
 hierarchical random graph model.

 \code{hrg.game} samples a graph from a given hierarchical random graph
 model.

 \code{hrg.predict} uses a hierarchical random graph model to predict
 missing edges from a network. This is done by sampling hierarchical
 models around the optimum model, proportionally to their
 likelihood. The MCMC sampling is stated from \code{hrg}, if it is given
 and the \code{start} argument is set to \code{TRUE}. Otherwise a HRG is
 fitted to the graph first.
}
\section{Printing HRGs to the screen}{
  \code{igraphHRG} objects can be printed to the screen in two forms: as
  a tree or as a list, depending on the \code{type} argument of the
  print function. By default the \code{auto} type is used, which selects
  \code{tree} for small graphs and \code{simple} (=list) for bigger
  ones. The \code{tree} format looks like
  this: \preformatted{Hierarchical random graph, at level 3:
g1        p=   0  
'- g15    p=0.33  1 
   '- g13 p=0.88  6  3  9  4  2  10 7  5  8 
'- g8     p= 0.5  
   '- g16 p= 0.2  20 14 17 19 11 15 16 13
   '- g5  p=   0  12 18  }
  This is a graph with 20 vertices, and the
  top three levels of the fitted hierarchical random graph are
  printed. The root node of the HRG is always vertex group #1
  (\sQuote{\code{g1}} in the the printout). Vertex pairs in the left
  subtree of \code{g1} connect to vertices in the right subtree with
  probability zero, according to the fitted model. \code{g1} has two
  subgroups, \code{g15} and \code{g8}. \code{g15} has a subgroup of a
  single vertex (vertex 1), and another larger subgroup that contains
  vertices 6, 3, etc. on lower levels, etc.

  The \code{plain} printing is simpler and faster to produce, but less
  visual: \preformatted{Hierarchical random graph:
g1  p=0.0 -> g12 g10   g2  p=1.0 -> 7 10      g3  p=1.0 -> g18 14    
g4  p=1.0 -> g17 15    g5  p=0.4 -> g15 17    g6  p=0.0 -> 1 4       
g7  p=1.0 -> 11 16     g8  p=0.1 -> g9 3      g9  p=0.3 -> g11 g16   
g10 p=0.2 -> g4 g5     g11 p=1.0 -> g6 5      g12 p=0.8 -> g8 8      
g13 p=0.0 -> g14 9     g14 p=1.0 -> 2 6       g15 p=0.2 -> g19 18    
g16 p=1.0 -> g13 g2    g17 p=0.5 -> g7 13     g18 p=1.0 -> 12 19     
g19 p=0.7 -> g3 20}
  It lists the two subgroups of each internal node, in
  as many columns as the screen width allows.

  Consensus dendrograms (\code{igraphHRGConsensus} objects) are printed
  simply by listing the children of each internal node of the
  dendrogram: \preformatted{HRG consensus tree:
g1 -> 11 12 13 14 15 16 17 18 19 20
g2 -> 1  2  3  4  5  6  7  8  9  10   
g3 -> g1 g2}
  The root of the dendrogram is \code{g3} (because it has no incoming
  edges), and it has two subgroups, \code{g1} and \code{g2}.
}
\value{
  \code{hrg.fit} returns an \code{igraphHRG} object. This is a list with
  the following members:
  \item{left}{Vector that contains the left children of the internal
    tree vertices. The first vertex is always the root vertex, so
    the first element of the vector is the left child of the root
    vertex. Internal vertices are denoted with negative numbers,
    starting from -1 and going down, i.e. the root vertex is
    -1. Leaf vertices are denoted by non-negative number, starting
    from zero and up.}
  \item{right}{Vector that contains the right children of the
    vertices, with the same encoding as the \code{left} vector.}
  \item{prob}{The connection probabilities attached to the internal
    vertices, the first number belongs to the root vertex
    (i.e. internal vertex -1), the second to internal vertex -2,
    etc.}
  \item{edges}{The number of edges in the subtree below the given
    internal vertex.}
  \item{vertices}{The number of vertices in the subtree below the
    given internal vertex, including itself.}

  \code{hrg.consensus} returns a list of two objects. The first is an
  \code{igraphHRGConsensus} object, the second is an \code{igraphHRG}
  object. The \code{igraphHRGConsensus} object has the following
  members:
  \item{parents}{For each vertex, the id of its parent vertex is stored,
    or zero, if the vertex is the root vertex in the tree. The first n
    vertex ids (from 0) refer to the original vertices of the graph,
    the other ids refer to vertex groups.}
  \item{weights}{Numeric vector, counts the number of times a given
    tree split occured in the generated network samples, for each
    internal vertices. The order is the same as in the \code{parents}
    vector.}

  \code{hrg.create} returns an \code{igraphHRG} object.

  \code{hrg.dendrogram} returns an igraph graph.

  \code{hrg.game} returns an igraph graph.
}
\references{
  A. Clauset, C. Moore, and M.E.J. Newman. Hierarchical structure and the
  prediction of missing links in networks. \emph{Nature} 453, 98--101
  (2008); 

  A. Clauset, C. Moore, and M.E.J. Newman. Structural Inference of
  Hierarchies in Networks. In E. M. Airoldi  et al. (Eds.): ICML 2006 Ws,
  \emph{Lecture Notes in Computer Science} 4503, 1--13. Springer-Verlag,
  Berlin Heidelberg (2007). 
}
\author{
  Gabor Csardi \email{csardi.gabor@gmail.com}, based on code from Aaron
  Clauset, thanks Aaron!
}
% \seealso{}
\examples{

## We are not running these examples any more, because they
## take a long time (~15 seconds) to run and this is against the CRAN
## repository policy. Copy and paste them by hand to your R prompt if
## you want to run them.

\dontrun{
## A graph with two dense groups
g <- erdos.renyi.game(10, p=1/2) + erdos.renyi.game(10, p=1/2)
hrg <- hrg.fit(g)
hrg

## The consensus tree for it
hrg.consensus(g, hrg=hrg, start=TRUE)

## Prediction of missing edges
g2 <- graph.full(4) + (graph.full(4) - path(1,2))
hrg.predict(g2)
}
}
\keyword{graphs}
