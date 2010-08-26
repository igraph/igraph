\name{dominator.tree}
\alias{dominator.tree}
\concept{Domintor tree}
\title{Dominator tree}
\description{Dominator tree of a directed graph.
}
\usage{
dominator.tree (graph, root, mode = c("out", "in"))
}
\arguments{
  \item{graph}{A directed graph. If it is not a flowgraph, and it
    contains some vertices not reachable from the root vertex,
    then these vertices will be collected and returned as part of the
    result.}
  \item{root}{The id of the root (or source) vertex, this will be the
    root of the tree.}
  \item{mode}{Constant, must be \sQuote{\code{in}} or
    \sQuote{\code{out}}. IF it is \sQuote{\code{in}}, then all
    directions are considered as opposite to the original one in the
    input graph.}
}
\details{
 A flowgraph is a directed graph with a distinguished start (or
 root) vertex \eqn{r}, such that for any vertex \eqn{v}, there is a path
 from \eqn{r} to \eqn{v}. A vertex \eqn{v} dominates another vertex
 \eqn{w} (not equal to \eqn{v}), if every path from \eqn{r} to \eqn{w}
 contains \eqn{v}. Vertex \eqn{v} is the immediate dominator or \eqn{w},
 \eqn{v=\textrm{idom}(w)}{v=idom(w)}, if \eqn{v} dominates \eqn{w} and
 every other dominator of \eqn{w} dominates \eqn{v}. The edges
 \eqn{{(\textrm{idom}(w), w)| w \ne r}}{{(idom(w),w)| w is not r}}
 form a directed tree, rooted at \eqn{r}, called the dominator tree of
 the graph. Vertex \eqn{v} dominates vertex \eqn{w} if and only if
 \eqn{v} is an ancestor of \eqn{w} in the dominator tree.

 This function implements the Lengauer-Tarjan algorithm
 to construct the dominator tree of a directed graph. For details see
 the reference below.
}
\value{
  A list with components:
  \item{dom}{
    A numeric vector giving the immediate dominators for each
    vertex. For vertices that are unreachable from the root, it contains
    \code{NaN}. For the root vertex itself it contains minus one.
  }
  \item{domtree}{
    A graph object, the dominator tree. Its vertex ids are the as the
    vertex ids of the input graph. Isolate vertices are the ones that
    are unreachable from the root.
  }
  \item{leftout}{
    A numeric vector containing the vertex ids that are unreachable from
    the root.
  }
}
\references{
  Thomas Lengauer, Robert Endre Tarjan: A fast algorithm
  for finding dominators in a flowgraph, \emph{ACM Transactions on
  Programming Languages and Systems (TOPLAS)} I/1, 121--141, 1979.
}
\author{Gabor Csardi \email{csardi.gabor@gmail.com}}
% \seealso{}
\examples{
## The example from the paper
g <- graph.formula(R-+A:B:C, A-+D, B-+A:D:E, C-+F:G, D-+L,
                   E-+H, F-+I, G-+I:J, H-+E:K, I-+K, J-+I,
                   K-+I:R, L-+H)
dtree <- dominator.tree(g, root="R")
layout <- layout.reingold.tilford(dtree$domtree, root="R")
layout[,2] <- -layout[,2]
if (interactive()) {
  plot(dtree$domtree, layout=layout, vertex.label=V(dtree$domtree)$name)
}
}
\keyword{graphs}
