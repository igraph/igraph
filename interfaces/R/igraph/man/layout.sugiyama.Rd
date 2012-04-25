\name{layout.sugiyama}
\alias{layout.sugiyama}
\concept{Graph layout}
\title{The Sugiyama graph layout generator}
\description{
  Sugiyama layout algorithm for layered directed acyclic graphs. The
  algorithm minimized edge crossings.
}
\usage{
layout.sugiyama (graph, layers = NULL, hgap = 1, vgap = 1, maxiter = 100, 
    weights = NULL) 
}
\arguments{
  \item{graph}{The input graph.}
  \item{layers}{A numeric vector or \code{NULL}. If not \code{NULL},
    then it should specify the layer index of the vertices. If
    \code{NULL}, then igraph calculates the layers automatically.}
  \item{hgap}{Real scalar, the minimum horizontal gap between vertices
    in the same layer.}
  \item{vgap}{Real scalar, the distance between layers.}
  \item{maxiter}{Integer scalar, the maximum number of iterations in the
    crossing minimization stage. 100 is a reasonable default; if you
    feel that you have too many edge crossings, increase this.}
  \item{weights}{Optional edge weight vector. If \code{NULL}, then the
    'weight' edge attribute is used, if there is one. Supply \code{NA}
    here and igraph ignores the edge weights.}
}
\details{
  This layout algorithm is designed for directed acyclic graphs where each
  vertex is assigned to a layer. Layers are indexed from zero, and vertices
  of the same layer will be placed on the same horizontal line. The X
  coordinates of vertices within each layer are decided by the heuristic
  proposed by Sugiyama et al. to minimize edge crossings.

  You can also try to lay out undirected graphs, graphs containing cycles, or
  graphs without an a priori layered assignment with this algorithm. igraph
  will try to eliminate cycles and assign vertices to layers, but there is no
  guarantee on the quality of the layout in such cases.

  The Sugiyama layout may introduce \dQuote{bends} on the edges in order
  to obtain a visually more pleasing layout. This is achieved by adding
  dummy nodes to edges spanning more than one layer. The resulting
  layout assigns coordinates not only to the nodes of the original graph
  but also to the dummy nodes. The layout algorithm will also return the
  extended graph with the dummy nodes.

  For more details, see the reference below.
}
\value{
  A list with the components:
  \item{res}{The layout, a two-column matrix, including coordinates for
    the dummy vertices. Use the top rows as coordinates for the original
    vertices.}
  \item{extd_graph}{The original graph, extended with dummy vertices.}
  \item{extd_to_orig_eids}{TODO}
}
\references{
  K. Sugiyama, S. Tagawa and M. Toda, "Methods for Visual Understanding
  of Hierarchical Systems". IEEE Transactions on Systems, Man and
  Cybernetics 11(2):109-125, 1981.
}
\author{
  Tamas Nepusz \email{ntamas@gmail.com}
}
\examples{
# TODO
}
\keyword{graphs}
