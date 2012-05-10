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
    then it should specify the layer index of the vertices. Layers are
    numbered from one. If \code{NULL}, then igraph calculates the layers
    automatically.}
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
  \item{extd_to_orig_eids}{Mapping from the edges of the extended graph
    to the edges of the original graph.}
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
## Data taken from http://tehnick-8.narod.ru/dc_clients/
DC <- graph.formula("DC++" -+
                "LinuxDC++":"BCDC++":"EiskaltDC++":"StrongDC++":"DiCe!++",
                "LinuxDC++" -+ "FreeDC++", "BCDC++" -+ "StrongDC++",
                "FreeDC++" -+ "BMDC++":"EiskaltDC++",
                "StrongDC++" -+ "AirDC++":"zK++":"ApexDC++":"TkDC++",
                "StrongDC++" -+ "StrongDC++ SQLite":"RSX++",
                "ApexDC++" -+ "FlylinkDC++ ver <= 4xx",
                "ApexDC++" -+ "ApexDC++ Speed-Mod":"DiCe!++",
                "StrongDC++ SQLite" -+ "FlylinkDC++ ver >= 5xx",
                "ApexDC++ Speed-Mod" -+ "FlylinkDC++ ver <= 4xx",
                "ApexDC++ Speed-Mod" -+ "GreylinkDC++",
                "FlylinkDC++ ver <= 4xx" -+ "FlylinkDC++ ver >= 5xx",
                "FlylinkDC++ ver <= 4xx" -+ AvaLink,
                "GreylinkDC++" -+ AvaLink:"RayLinkDC++":"SparkDC++":PeLink)

## Layers, as on the plot
layers <- list(c("DC++"),
               c("LinuxDC++", "BCDC++"),
               c("FreeDC++", "StrongDC++"),
               c("BMDC++", "EiskaltDC++", "AirDC++", "zK++", "ApexDC++",
                 "TkDC++", "RSX++"),
               c("StrongDC++ SQLite", "ApexDC++ Speed-Mod", "DiCe!++"),
               c("FlylinkDC++ ver <= 4xx", "GreylinkDC++"),
               c("FlylinkDC++ ver >= 5xx", "AvaLink", "RayLinkDC++",
                 "SparkDC++", "PeLink"))

## Check that we have all nodes
all(sort(unlist(layers)) == sort(V(DC)$name))

## Create a similar layout using the predefined layers
lay1 <- layout.sugiyama(DC, layers=apply(sapply(layers,
                        function(x) V(DC)$name \%in\% x), 1, which))

## Flip vertically
lay1$res[,2] <- max(lay1$res[,2]) - lay1$res[,2] + 1

## Simple plot, not vey nice
labels <- lapply(V(DC)$name, function(x)
          paste(strwrap(x, 12), collapse="\n"))
par(mar=rep(.1, 4))
plot(DC, layout=lay1$res[1:vcount(DC),], vertex.label=labels,
     edge.arrow.size=.5, vertex.color="white",
     vertex.shape="rectangle", vertex.label.cex=0.5,
     vertex.size=20, vertex.size2=10)

## Sugiyama plot
origvert <- c(rep(TRUE, vcount(DC)), rep(FALSE, nrow(lay1$res)-vcount(DC)))
realedge <- get.edgelist(lay1$extd_graph)[,2] <= vcount(DC)
plot(lay1$extd_graph, layout=lay1$res, vertex.label.cex=0.5,
     edge.arrow.size=.5, vertex.color="white",
     vertex.size=ifelse(origvert, 20, 0),
     vertex.size2=ifelse(origvert, 10, 0),
     vertex.shape=ifelse(origvert, "rectangle", "none"),
     vertex.label=ifelse(origvert, labels, ""),
     edge.arrow.mode=ifelse(realedge, 2, 0))

## The same with automatic layer calculation
lay2 <- layout.sugiyama(DC)
lay2$res[,2] <- max(lay2$res[,2]) - lay2$res[,2] + 1
origvert <- c(rep(TRUE, vcount(DC)), rep(FALSE, nrow(lay2$res)-vcount(DC)))
realedge <- get.edgelist(lay2$extd_graph)[,2] <= vcount(DC)
plot(lay2$extd_graph, layout=lay2$res, vertex.label.cex=0.5,
     edge.arrow.size=.5, vertex.color="white",
     vertex.size=ifelse(origvert, 20, 0),
     vertex.size2=ifelse(origvert, 10, 0),
     vertex.shape=ifelse(origvert, "rectangle", "none"),
     vertex.label=ifelse(origvert, labels, ""),
     edge.arrow.mode=ifelse(realedge, 2, 0))

## Use edge types
E(DC)$lty <- 1
E(DC)["BCDC++" \%->\% "StrongDC++"]$lty <- 2
E(DC)["FreeDC++" \%->\% "EiskaltDC++"]$lty <- 2
E(DC)["ApexDC++" \%->\% "FlylinkDC++ ver <= 4xx"]$lty <- 2
E(DC)["ApexDC++" \%->\% "DiCe!++"]$lty <- 2
E(DC)["StrongDC++ SQLite" \%->\% "FlylinkDC++ ver >= 5xx"]$lty <- 2
E(DC)["GreylinkDC++" \%->\% "AvaLink"]$lty <- 2

plot(lay2$extd_graph, layout=lay2$res, vertex.label.cex=0.5,
     edge.arrow.size=.5, vertex.color="white",
     vertex.size=ifelse(origvert, 20, 0),
     vertex.size2=ifelse(origvert, 10, 0),
     vertex.shape=ifelse(origvert, "rectangle", "none"),
     vertex.label=ifelse(origvert, labels, ""),
     edge.arrow.mode=ifelse(realedge, 2, 0),
     edge.lty=E(DC)$lty[lay2$extd_to_orig_eids])

}
\keyword{graphs}
