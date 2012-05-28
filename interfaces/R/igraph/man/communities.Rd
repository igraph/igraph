\name{communities}
\alias{communities}
\alias{membership}
\alias{algorithm}
\alias{crossing}
\alias{cutat}
\alias{merges}
\alias{sizes}
\alias{is.hierarchical}
\alias{print.communities}
\alias{plot.communities}
\alias{length.communities}
\alias{modularity.communities}
\alias{as.dendrogram.communities}
\alias{as.hclust.communities}
\alias{asPhylo}
\alias{asPhylo.communities}
\alias{showtrace}
\alias{code.length}
\concept{Community structure}
\title{Functions to deal with the result of network community detection}
\description{igraph community detection functions return their results
  as an object from the \code{communities} class. This manual page
  describes the operations of this class.}
\usage{
\method{print}{communities}(x, \dots)

\method{length}{communities}(x)
sizes(communities)
membership(communities)
\method{modularity}{communities}(x, \dots)
algorithm(communities)
crossing(communities, graph)

is.hierarchical(communities, full = FALSE)
merges(communities)
cutat(communities, no, steps)
\method{as.dendrogram}{communities}(object, hang=-1,
   use.modularity=FALSE, \dots)
\method{as.hclust}{communities}(object, use.modularity=FALSE, \dots)
\method{asPhylo}{communitues}(x)
showtrace(communities)

code.length(communities)

\method{plot}{communities}(x, y,
   colbar=rainbow(length(x)),
   col=colbar[membership(x)],
   mark.groups=communities(x),
   edge.color=c("black", "red")[crossing(x,y)+1],
   \dots)
}
\arguments{
  \item{communities,x,object}{A \code{communities} object, the result of
    an igraph community detection function.}
  \item{graph}{An igraph graph object, corresponding to
    \code{communities}.}
  \item{full}{Logical scalar, if \code{TRUE}, then
    \code{is.hierarchical} only returns \code{TRUE} for fully
   hierarchical algorithms. The \sQuote{leading eigenvector} algorithm
   is hierarchical, it gives a hierarchy of groups, but not a full
   dendrogram with all vertices, so it is not fully hierarchical.}
  \item{y}{An igraph graph object, corresponding to the communities in
    \code{x}.}
  \item{no}{Integer scalar, the desired number of communities. If too
    low or two high, then an error message is given. Exactly one of
    \code{no} and \code{steps} must be supplied.}
  \item{steps}{The number of merge operations to perform to produce the
    communities. Exactly one of \code{no} and \code{steps} must be
    supplied.}
  \item{colbar}{A vector of colors, in any format that is accepted by
    the regular R plotting methods. E.g. it may be an integer vector, a
    character vector of color names, a character vector of RGB
    colors. This vector gives the color bar for the vertices. The length
    of the vector should be the same as the number of communities.}
  \item{col}{A vector of colors, in any format that is accepted by the
    regular R plotting methods. This vector gives the colors of the
    vertices explicitly.}
  \item{mark.groups}{A list of numeric vectors. The communities can be
    highlighted using colored polygons. The groups for which the
    polygons are drawn are given here. The default is to use the groups
    given by the communities. Supply \code{NULL} here if you do not want
    to highlight any groups.}
  \item{edge.color}{The colors of the edges. By default the edges within
    communities are colored green and other edges are red.}
  \item{hang}{Numeric scalar indicating how the height of leaves should
    be computed from the heights of their parents; see
   \code{\link{plot.hclust}}.}
 \item{use.modularity}{Logical scalar, whether to use the modularity
   values to define the height of the branches.}
 \item{\dots}{Additional arguments. \code{plot.communities} passes
   these to \code{\link{plot.igraph}}. The other functions silently
   ignore them.}
}
\details{
  Community structure detection algorithms try to find dense subgraphs
  in directed or undirected graphs, by optimizing some criteria, and
  usually using heuristics.

  igraph implements a number of commmunity detection methods (see them
  below), all of which return an object of the class
  \code{communities}. Because the community structure detection
  algorithms are different, \code{communities} objects do not always
  have the same structure. Nevertheless, they have some common
  operations, these are documented here.

  The \code{print} generic function is defined for \code{communities},
  it prints a short summary.

  The \code{length} generic function call be called on
  \code{communities} and returns the number of communities.

  The \code{sizes} function returns the community sizes, in the order of
  their ids.
  
  \code{membership} gives the division of the vertices, into
  communities. It returns a numeric vector, one value for each vertex,
  the id of its community. Community ids start from one. Note that some
  algorithms calculate the complete (or incomplete) hierarchical
  structure of the communities, and not just a single
  partitioning. For these algorithms typically the membership for the
  highest modularity value is returned, but see also the manual pages of
  the individual algorithms.
  
  \code{modularity} gives the modularity score of the partitioning. (See
  \code{\link{modularity.igraph}} for details. For algorithms that do
  not result a single partitioning, the highest modularity value is
  returned.

  \code{algorithm} gives the name of the algorithm that was used to
  calculate the community structure.

  \code{crossing} returns a logical vector, with one value for each
  edge, ordered according to the edge ids. The value is \code{TRUE} iff
  the edge connects two different communities, according to the (best)
  membership vector, as returned by \code{membership()}.

  \code{is.hierarchical} checks whether a hierarchical algorithm was
  used to find the community structure. Some functions only make sense
  for hierarchical methods (e.g. \code{merges}, \code{cutat} and
  \code{as.dendrogram}).
  
  \code{merges} returns the merge matrix for hierarchical methods. An
  error message is given, if a non-hierarchical method was used to find
  the community structure. You can check this by calling
  \code{is.hierarchical} on the \code{communities} object.

  \code{cutat} cuts the merge tree of a hierarchical community finding
  method, at the desired place and returns a membership vector. The
  desired place can be expressed as the desired number of communities or
  as the number of merge steps to make. The function gives an error
  message, if called with a non-hierarchical method.

  \code{as.dendrogram} converts a hierarchical community structure to a
  \code{dendrogram} object. It only works for hierarchical methods, and
  gives an error message to others. See \code{\link[stats]{dendrogram}}
  for details.

  \code{as.hclust} is similar to \code{as.dendrogram}, but converts a
  hierarchical community structure to a \code{hclust} object.

  \code{asPhylo} converts a hierarchical community structure to
  a \code{phylo} object, you will need the \code{ape} package for this.
  
  \code{showtrace} works (currently) only for communities found by the
  leading eigenvector method
  (\code{\link{leading.eigenvector.community}}), and returns a character
  vector that gives the steps performed by the algorithm while finding
  the communities.

  \code{code.length} is defined for the InfoMAP method
  (\code{\link{infomap.community}} and returns the code length of the
  partition.
  
  It is possibly to call the \code{plot} function on \code{communities}
  objects. This will plot the graph (and uses \code{\link{plot.igraph}}
  internally), with the communities shown. By default it colores the
  vertices according to their communities, and also marks the vertex
  groups corresponding to the communities. It passes additional
  arguments to \code{\link{plot.igraph}}, please see that and also
  \code{\link{igraph.plotting}} on how to change the plot.
}
\value{
  \code{print} returns the \code{communities} object itself, invisibly.

  \code{length} returns an integer scalar.

  \code{sizes} returns a numeric vector.
  
  \code{membership} returns a numeric vector, one number for each vertex
  in the graph that was the input of the community detection.
  
  \code{modularity} returns a numeric scalar.

  \code{algorithm} returns a character scalar.

  \code{crossing} returns a logical vector.

  \code{is.hierarchical} returns a logical scalar.

  \code{merges} returns a two-column numeric matrix.

  \code{cutat} returns a numeric vector, the membership vector of the
  vertices.

  \code{as.dendrogram} returns a \code{\link[stats]{dendrogram}}
  object.

  \code{showtrace} returns a character vector.

  \code{code.length} returns a numeric scalar for communities found with
  the InfoMAP method and \code{NULL} for other methods.
  
  \code{plot} for \code{communities} objects returns \code{NULL},
  invisibly.
}
% \references{}
\author{Gabor Csardi \email{csardi.gabor@gmail.com}}
\seealso{
  \code{\link{compare.communities}} for comparing two community
  structures on the same graph.
  
  The different methods for finding communities, they all return a
  \code{communities} object:
  \code{\link{edge.betweenness.community}},
  \code{\link{fastgreedy.community}},
  \code{\link{label.propagation.community}},
  \code{\link{leading.eigenvector.community}},
  \code{\link{multilevel.community}},
  \code{\link{optimal.community}},
  \code{\link{spinglass.community}},
  \code{\link{walktrap.community}}.
}
\examples{
karate <- graph.famous("Zachary")
wc <- walktrap.community(karate)
modularity(wc)
membership(wc)
plot(wc, karate)
}
\keyword{graphs}

