\name{dendPlot}
\alias{dendPlot}
\alias{dendPlot.communities}
\concept{Dendrograms}
\concept{Community structure}
\title{Plot fully hierarchical community structures as dendrograms}
\description{
  Plot a dendrogram that corresponds to a fully hierarchical community
  structure.
}
\usage{
\method{dendPlot}{communities}(communities,
        mode = getIgraphOpt("dend.plot.type"), \dots,
        use.modularity = FALSE)
}
\arguments{
  \item{communities}{A \code{communities} object, the community
    structure, as returned by some igraph community detection function.}
  \item{mode}{Which dendrogram plotting function to use. See details
    below.}
  \item{\dots}{Additional arguments to supply to the dendrogram plotting
    function.}
  \item{use.modularity}{Logical scalar, whether to use the modularity
    scores calculated by the community detection algorithm, to set the
    heights of the dendrogram branches. Right now this must be
    \code{FALSE} for undirected graphs, for technical readons.}
}
\details{
  \code{dendPlot} supports three different plotting functions, selected
  via the \code{mode} argument. By default the plotting function is
  taken from the \code{dend.plot.type} igraph option, and it has for
  possible values: \itemize{
  \item \code{auto} Choose automatically between the plotting
    functions. As \code{plot.phylo} is the most sophisticated, that is
    choosen, whenever the \code{ape} package is available. Otherwise
    \code{plot.hclust} is used.
  \item \code{phylo} Use \code{plot.phylo} from the \code{ape}
    package.
  \item \code{hclust} Use \code{plot.hclust} from the \code{stats}
    package.
  \item \code{dendrogram} Use \code{plot.dendrogram} from the
    \code{stats} package.
  }

  The different plotting functions take different sets of
  arguments. When using \code{plot.phylo} (\code{mode="phylo"}), we have
  the following syntax: \preformatted{
    dendPlot(communities, mode="phylo",
         colbar = rainbow(length(communities)),
         col = colbar[membership(communities)],
         mark.groups = communities(communities),
         edge.color = "#AAAAAAFF", edge.lty = c(1,2), \dots,
         use.modularity = FALSE)
  } The extra arguments not documented above: \itemize{
  \item \code{colbar} This is a vector of colors, in any format
    supported by R. These colors will be used to mark the communities on
    the dendrogram, if \code{mark.groups} are given.
  \item \code{col} The colors of the vertex labels, by default the
    community membership is color-coded here, using the colors from the
    \code{colbar} argument.
  \item \code{mark.groups} A list of groups to mark on the tree. Groups
    must correspond to disjunct subtrees of the dendrogram, otherwise
    you might get strange results. This argument is interpreted as a
    list of vertex numeric ids.
  \item \code{edge.color} The default branch color. This gives the color
    of the dendrogram branches that are not marked according to the
    \code{mark.groups} argument.
  \item \code{edge.lty} The line types of the branches, a vector of two
    line types, the first value is for branches marked via
    \code{mark.groups}, the second value is for unmarked branches. See
    the documentation of \code{lty} at the manual page of \code{par} for
    possible values. By default unmarked branches use dashed, marked
    branches use solid lines.
  \item \code{\dots} Additional arguments to pass to \code{plot.phylo}.
  }

  The syntax for \code{plot.hclust} (\code{mode="hclust"}): \preformatted{
    dendPlot(communities, mode = "hclust", rect = length(communities),
         colbar = rainbow(rect), hang = -1, ann = FALSE, main = "",
	 sub = "", xlab = "", ylab = "", \dots, use.modularity = FALSE)
  } The extra arguments not documented above: \itemize{
  \item \code{rect} A numeric constant, the number of groups to mark on
    the dendrogram. The dendrogram is cut into exactly \code{rect}
    groups and they are marked via the \code{rect.hclust} command. Set
    this to zero if you don't want to mark any groups.
  \item \code{colbar} The colors of the rectanges that mark the
    communities via the \code{rect} argument.
  \item \code{hang} Where to put the leaf nodes, this corresponds to the
    \code{hang} argument of \code{plot.hclust}.
  \item \code{ann} Whether to annotate the plot, the \code{ann} argument
    of \code{plot.hclust}.
  \item \code{main} The main title of the plot, the \code{main} argument
    of \code{plot.hclust}.
  \item \code{sub} The sub-title of the plot, the \code{sub} argument
    of \code{plot.hclust}.
  \item \code{xlab} The label on the horizontal axis, passed to
    \code{plot.hclust}. 
  \item \code{ylab} The label on the vertical axis, passed to
    \code{plot.hclust}. 
  \item \code{dots} Attitional arguments to pass to \code{plot.hclust}.
  }

  The syntax for \code{plot.dendrogram} (\code{mode="dendrogram"}): \preformatted{
    dendPlot(communities, mode="dendrogram", hang = -1, \dots,
    use.modularity = FALSE)
  } The extra argument not documented above: \itemize{
  \item \code{hang} Where to place the leaf nodes on the plot, this is
    passed to \code{as.dendrogram}.
  }
}
\value{
  Returns whatever the return value was from the plotting function,
  \code{plot.phylo}, \code{plot.dendrogram} or \code{plot.hclust}.
}
%\references{}
\author{Gabor Csardi \email{csardi.gabor@gmail.com}}
\seealso{ \code{\link{communities}} for other operations on
  \code{communities} obejcts.
}
\examples{
karate <- graph.famous("Zachary")
fc <- fastgreedy.community(karate)
dendPlot(fc)
}
\keyword{graphs}
