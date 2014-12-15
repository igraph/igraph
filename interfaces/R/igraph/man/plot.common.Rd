\name{Drawing graphs}
\alias{igraph.plotting}
\concept{Visualization}
\title{Drawing graphs}
\description{The common bits of the three plotting functions
  \code{plot.igraph}, \code{tkplot} and \code{rglplot} are discussed in
  this manual page}
\details{
  There are currently three different functions in the igraph package
  which can draw graph in various ways:

  \code{plot.igraph} does simple non-interactive 2D plotting to R devices.
  Actually it is an implementation of the \code{\link[graphics]{plot}}
  generic function, 
  so you can write \code{plot(graph)} instead of
  \code{plot.igraph(graph)}. As it used the standard R devices it
  supports every output format for which R has an output device. The
  list is quite impressing: PostScript, PDF files, XFig files, SVG
  files, JPG, PNG and of course you can plot to the screen as well using
  the default devices, or the good-looking anti-aliased Cairo device.
  See \code{\link{plot.igraph}} for some more information.

  \code{\link{tkplot}} does interactive 2D plotting using the \code{tcltk}
  package. It can only handle graphs of moderate size, a thousend
  vertices is probably already too many. Some parameters of the plotted
  graph can be changed interactively after issuing the \code{tkplot}
  command: the position, color and size of the vertices and the color
  and width of the edges. See \code{\link{tkplot}} for details.
    
  \code{\link{rglplot}} is an experimental function to draw graphs in 3D
  using OpenGL. See \code{\link{rglplot}} for some more information.

  Please also check the examples below.
}
\section{How to specify graphical parameters}{
  There are three ways to give values to the parameters described below,
  in section 'Parameters'. We give these three ways here in the order of
  their precedence. 

  The first method is to supply named arguments to the plotting commands:
  \code{\link{plot.igraph}}, \code{\link{tkplot}} or
  \code{\link{rglplot}}. Parameters for vertices start with prefix
  \sQuote{\code{vertex.}}, parameters for edges have prefix
  \sQuote{\code{edge.}}, and global parameters have no prefix. Eg. the
  color of the vertices can be given via argument \code{vertex.color},
  whereas \code{edge.color} sets the color of the edges. \code{layout}
  gives the layout of the graphs.

  The second way is to assign vertex, edge and graph attributes to the
  graph. These attributes have no prefix, ie. the color of the vertices
  is taken from the \code{color} vertex attribute and the color of the
  edges from the \code{color} edge attribute. The layout of the graph is
  given by the \code{layout} graph attribute. (Always assuming that the
  corresponding command argument is not present.) Setting vertex and
  edge attributes are handy if you want to assign a given \sQuote{look}
  to a graph, attributes are saved with the graph is you save it with
  \code{\link[base]{save}} or in GraphML format with
  \code{\link{write_graph}}, so the graph will have the same look after
  loading it again.

  If a parameter is not given in the command line, and the corresponding
  vertex/edge/graph attribute is also missing then the general igraph
  parameters handled by \code{\link{igraph_options}} are also
  checked. Vertex parameters have prefix \sQuote{\code{vertex.}}, edge
  parameters are prefixed with \sQuote{\code{edge.}}, general parameters
  like \code{layout} are prefixed with \sQuote{\code{plot}}.
  These parameters are useful if you want
  all or most of your graphs to have the same look, vertex size, vertex
  color, etc. Then you don't need to set these at every plotting, and
  you also don't need to assign vertex/edge attributes to every graph.

  If the value of a parameter is not specified by any of the three ways
  described here, its default valued is used, as given in the source
  code.

  Different parameters can have different type, eg. vertex colors can be
  given as a character vector with color names, or as an integer vector
  with the color numbers from the current palette. Different types are
  valid for different parameters, this is discussed in detail in the
  next section. It is however always true that the parameter can always
  be a function object in which it will be called with the graph as its
  single argument to get the \dQuote{proper} value of the parameter.
  (If the function returns another function object that will \emph{not}
  be called again\dots)
}
\section{The list of parameters}{
  Vertex parameters first, note that the \sQuote{\code{vertex.}} prefix
  needs to be added if they are used as an argument or when setting via
  \code{\link{igraph_options}}. The value of the parameter may be scalar
  valid for every vertex or a vector with a separate value for each
  vertex. (Shorter vectors are recycled.)
  \describe{
    \item{size}{The size of the vertex, a numeric scalar or vector, in
      the latter case each vertex sizes may differ. This vertex sizes
      are scaled in order have about the same size of vertices for a
      given value for all three plotting commands. It does not need to
      be an integer number.

      The default value is 15. This is big enough to place short labels
      on vertices.}
    \item{size2}{The \dQuote{other} size of the vertex, for some vertex
      shapes. For the various rectangle shapes this gives the height of
      the vertices, whereas \code{size} gives the width. It is ignored
      by shapes for which the size can be specified with a single
      number.

      The default is 15.
    }
    \item{color}{The fill color of the vertex. If it is numeric then
      the current palette is used, see
      \code{\link[grDevices]{palette}}. If it is a character vector then
      it may either contain integer values, named colors or RGB
      specified colors with three or four bytes. All strings starting
      with \sQuote{\code{#}} are assumed to be RGB color
      specifications. It is possible to mix named color and RGB
      colors. Note that \code{\link{tkplot}} ignores the fourth byte
      (alpha channel) in the RGB color specification.

      For \code{plot.igraph} and integer values, the default igraph
      palette is used (see the \sQuote{palette} parameter below. Note
      that this is different from the R palette.

      If you don't want (some) vertices to have any color, supply
      \code{NA} as the color name.
      
      The default value is \dQuote{\code{SkyBlue2}}.
    }
    \item{frame.color}{The color of the frame of the vertices, the same
      formats are allowed as for the fill color.

      If you don't want vertices to have a frame, supply \code{NA} as
      the color name.

      By default it is \dQuote{black}.
    }
    \item{shape}{The shape of the vertex, currently
      \dQuote{\code{circle}}, \dQuote{\code{square}},
      \dQuote{\code{csquare}}, \dQuote{\code{rectangle}},
      \dQuote{\code{crectangle}}, \dQuote{\code{vrectangle}},
      \dQuote{\code{pie}} (see \link{vertex.shape.pie}),
      \sQuote{\code{sphere}}, and \dQuote{\code{none}} are supported,
      and only by the \code{\link{plot.igraph}}
      command. \dQuote{\code{none}} does not draw the vertices at all,
      although vertex label are plotted (if given). See
      \code{\link{shapes}} for details about vertex
      shapes and \code{\link{vertex.shape.pie}} for using pie charts as
      vertices.

      The \dQuote{\code{sphere}} vertex shape plots vertices as 3D
      ray-traced spheres, in the given color and size. This produces a
      raster image and it is only supported with some graphics
      devices. On some devices raster transparency is not supported and
      the spheres do not have a transparent background. See
      \code{\link{dev.capabilities}} and the \sQuote{\code{rasterImage}}
      capability to check that your device is supported.
      
      By default vertices are drawn as circles.
    }
    \item{label}{The vertex labels. They will be converted to
      character. Specify \code{NA} to omit vertex labels.

      The default vertex labels are the vertex ids.
    }
    \item{label.family}{The font family to be used for vertex labels.
      As different plotting commands can used different fonts, they
      interpret this parameter different ways. The basic notation is,
      however, understood by both \code{\link{plot.igraph}} and
      \code{\link{tkplot}}. \code{\link{rglplot}} does not support fonts
      at all right now, it ignores this parameter completely.

      For \code{\link{plot.igraph}} this parameter is simply passed to
      \code{\link[graphics]{text}} as argument \code{family}.

      For \code{\link{tkplot}} some
      conversion is performed. If this parameter is the name of an
      exixting Tk font, then that font is used and the \code{label.font}
      and \code{label.cex} parameters are ignored complerely. If it is
      one of the base families (serif, sans, mono) then Times,
      Helvetica or Courier fonts are used, there are guaranteed to exist
      on all systems. For the \sQuote{symbol} base family we used the
      symbol font is available, otherwise the first font which has
      \sQuote{symbol} in its name. If the parameter is not a name of the
      base families and it is also not a named Tk font then we pass it
      to \code{\link[tcltk]{tkfont.create}} and hope the user knows what
      she is doing. The \code{label.font} and \code{label.cex}
      parameters are also passed to \code{\link[tcltk]{tkfont.create}}
      in this case.

      The default value is \sQuote{serif}.
    }
    \item{label.font}{The font within the font family to use for the
      vertex labels. It is interpreted the same way as the the
      \code{font} graphical parameter: 1 is plain text, 2 is bold face,
      3 is italic, 4 is bold and italic and 5 specifies the symbol
      font.

      For \code{\link{plot.igraph}} this parameter is simply passed to
      \code{\link[graphics]{text}}.

      For \code{\link{tkplot}}, if the \code{label.family} parameter is
      not the name of a Tk font then this parameter is used to set
      whether the newly created font should be italic and/or
      boldface. Otherwise it is ignored.

      For \code{\link{rglplot}} it is ignored.

      The default value is 1.
    }
    \item{label.cex}{The font size for vertex labels. It is interpreted
      as a multiplication factor of some device-dependent base font
      size.

      For \code{\link{plot.igraph}} it is simply passed to
      \code{\link[graphics]{text}} as argument \code{cex}.

      For \code{\link{tkplot}} it is multiplied by 12 and then used as
      the \code{size} argument for \code{\link[tcltk]{tkfont.create}}.
      The base font is thus 12 for tkplot.

      For \code{\link{rglplot}} it is ignored.

      The default value is 1.
    }
    \item{label.dist}{
      The distance of the label from the center of the vertex. If it is 0
      then the label is centered on the vertex. If it is 1 then the
      label is displayed beside the vertex.

      The default value is 0.
    }
    \item{label.degree}{
      It defines the position of the vertex labels, relative to the
      center of the vertices. It is interpreted as an angle in radian,
      zero means \sQuote{to the right}, and \sQuote{\code{pi}} means to
      the left, up is \code{-pi/2} and down is \code{pi/2}.

      The default value is \code{-pi/4}.
    }
    \item{label.color}{The color of the labels, see the \code{color}
      vertex parameter discussed earlier for the possible values.

      The default value is \code{black}.
    }    
  }

  Edge parameters require to add the \sQuote{\code{edge.}} prefix when
  used as arguments or set by \code{\link{igraph_options}}. The edge
  parameters:
  \describe{
    \item{color}{The color of the edges, see the \code{color} vertex
      parameter for the possible values.

      By default this parameter is \code{darkgrey}.
    }
    \item{width}{The width of the edges.

      The default value is 1.
    }
    \item{arrow.size}{The size of the arrows. Currently this is a
      constant, so it is the same for every edge. If a vector is
      submitted then only the first element is used, ie. if this is
      taken from an edge attribute then only the attribute of the first
      edge is used for all arrows. This will likely change in the future.

      The default value is 1.
    }
    \item{arrow.width}{The width of the arrows. Currently this is a
      constant, so it is the same for every edge. If a vector is
      submitted then only the first element is used, ie. if this is
      taken from an edge attribute then only the attribute of the first
      edge is used for all arrows. This will likely change in the future.

      This argument is currently only used by \code{\link{plot.igraph}}.
      
      The default value is 1, which gives the same width as before
      this option appeared in igraph.
    }
    \item{lty}{The line type for the edges. Almost the same format is
      accepted as for the standard graphics \code{\link[graphics]{par}},
      0 and \dQuote{blank} mean no edges, 1 and \dQuote{solid} are for
      solid lines, the other possible values are: 2 (\dQuote{dashed}),
      3 (\dQuote{dotted}), 4 (\dQuote{dotdash}), 5 (\dQuote{longdash}),
      6 (\dQuote{twodash}).

      \code{\link{tkplot}} also accepts standard Tk line type strings,
      it does not however support \dQuote{blank} lines, instead of type
      \sQuote{0} type \sQuote{1}, ie. solid lines will be drawn.
      
      This argument is ignored for \code{\link{rglplot}}.

      The default value is type 1, a solid line.
    }
    \item{label}{The edge labels. They will be converted to
      character. Specify \code{NA} to omit edge labels.

      Edge labels are omitted by default.}
    \item{label.family}{Font family of the edge labels. See the vertex
      parameter with the same name for the details.}
    \item{label.font}{The font for the edge labels. See the
      corresponding vertex parameter discussed earlier for details.}
    \item{label.cex}{The font size for the edge labels, see the
      corresponding vertex parameter for details.}
    \item{label.color}{The color of the edge labels, see the
      \code{color} vertex parameters on how to specify colors. }
    \item{label.x}{The horizontal coordinates of the edge labels might
      be given here, explicitly. The \code{NA} elements will be
      replaced by automatically calculated coordinates. If \code{NULL},
      then all edge horizontal coordinates are calculated
      automatically. This parameter is only supported by
      \code{plot.igraph}.}
    \item{label.y}{The same as \code{label.x}, but for vertical
      coordinates.}
    \item{curved}{Specifies whether to draw curved edges, or not. This
      can be a logical or a numeric vector or scalar. 

      First the vector is replicated to have the same length as the
      number of edges in the graph. Then it is interpreted for each edge
      separately. A numeric value specifies the curvature of the edge;
      zero curvature means straight edges, negative values means the
      edge bends clockwise, positive values the opposite. \code{TRUE}
      means curvature 0.5, \code{FALSE} means curvature zero.

      By default the vector specifying the curvatire is calculated via a
      call to the \code{\link{curve_multiple}} function. This function makes
      sure that multiple edges are curved and are all visible. This
      parameter is ignored for loop edges. 
      
      The default value is \code{FALSE}.
      
      This parameter is currently ignored by \code{\link{rglplot}}.}
    \item{arrow.mode}{This parameter can be used to specify for which
      edges should arrows be drawn. If this parameter is given by the
      user (in either of the three ways) then it specifies which edges
      will have forward, backward arrows, or both, or no arrows at all.
      As usual, this parameter can be a vector or a scalar value.
      It can be an integer or character type. If it is integer then 0
      means no arrows, 1 means backward arrows, 2 is for forward arrows
      and 3 for both. If it is a character vector then \dQuote{<} and
      \dQuote{<-} specify backward, \dQuote{>} and \dQuote{->} forward
      arrows and \dQuote{<>} and \dQuote{<->} stands for both arrows.
      All other values mean no arrows, perhaps you should use \dQuote{-}
      or \dQuote{--} to specify no arrows.

      Hint: this parameter can be used as a \sQuote{cheap} solution for
      drawing \dQuote{mixed} graphs: graphs in which some edges are
      directed some are not. If you want do this, then please create a
      \emph{directed} graph, because as of version 0.4 the vertex pairs
      in the edge lists can be swapped in undirected graphs.

      By default, no arrows will be drawn for undirected graphs, and for 
      directed graphs, an arrow will be drawn for each edge, according
      to its direction. This is not very surprising, it is the expected
      behavior.
    }
    \item{loop.angle}{Gives the angle in radian for plotting loop
      edges. See the \code{label.dist} vertex parameter to see how this
      is interpreted.

      The default value is 0.
    }    
    \item{loop.angle2}{Gives the second angle in radian for plotting
      loop edges. This is only used in 3D, \code{loop.angle} is enough
      in 2D.

      The default value is 0.
    }
  }

  Other parameters:
  \describe{
    \item{layout}{
      Either a function or a numeric
      matrix. It specifies how the vertices will be placed on the plot.
      
      If it is a numeric matrix, then the matrix has to have one line for
      each vertex, specifying its coordinates. The matrix should have at
      least two columns, for the \code{x} and \code{y} coordinates, and
      it can also have third column, this will be the \code{z}
      coordinate for 3D plots and it is ignored for 2D plots.

      If a two column matrix is given for the 3D plotting function
      \code{\link{rglplot}} then the third column is assumed to be 1 for
      each vertex.

      If \code{layout} is a function, this function will be called with the
      \code{graph} as the single parameter to determine the
      actual coordinates. The function should return a matrix with two
      or three columns. For the 2D plots the third column is ignored.

      The default value is \code{layout_nicely}, a smart function that
      chooses a layouter based on the graph.}
    \item{margin}{The amount of empty space below, over, at the left and
      right of the plot, it is a numeric vector of length four. Usually
      values between 0 and 0.5 are meaningful, but negative values are
      also possible, that will make the plot zoom in to a part of the
      graph. If it is shorter than four then it is recycled.

      \code{\link{rglplot}} does not support this parameter, as it can
      zoom in and out the graph in a more flexible way.
      
      Its default value is 0.
    }
    \item{palette}{The color palette to use for vertex color.
      The default is \code{\link{categorical_pal}}, which is a
      color-blind friendly categorical palette. See its manual page
      for details and other palettes. This parameters is only supported
      by \code{plot}, and not by \code{tkplot} and \code{rglplot}.
    }
    \item{rescale}{Logical constant, whether to rescale the coordinates
      to the [-1,1]x[-1,1](x[-1,1]) interval. This parameter is not
      implemented for \code{tkplot}.

      Defaults to \code{TRUE}, the layout will be rescaled.
    }
    \item{asp}{A numeric constant, it gives the \code{asp} parameter
      for \code{\link{plot}}, the aspect ratio. Supply 0 here if you
      don't want to give an aspect ratio. It is ignored by \code{tkplot}
      and \code{rglplot}.

      Defaults to 1.
    }
    \item{frame}{Boolean, whether to plot a frame around the graph. It
      is ignored by \code{tkplot} and \code{rglplot}.

      Defaults to \code{FALSE}.
    }
    \item{main}{Overall title for the main plot. The default is empty if
      the \code{annotate.plot} igraph option is \code{FALSE}, and the
      graph's \code{name} attribute otherwise. See the same argument of
      the base \code{plot} function. Only supported by \code{plot}.}
    \item{sub}{Subtitle of the main plot, the default is empty. Only
      supported by \code{plot}.}
    \item{xlab}{Title for the x axis, the default is empty if the
      \code{annotate.plot} igraph option is \code{FALSE}, and the number
      of vertices and edges, if it is \code{TRUE}. Only supported by
      \code{plot}.}
    \item{ylab}{Title for the y axis, the default is empty. Only
      supported by \code{plot}.}
  }
}
\author{Gabor Csardi \email{csardi.gabor@gmail.com}}
\seealso{ \code{\link{plot.igraph}}, \code{\link{tkplot}},
  \code{\link{rglplot}}, \code{\link{igraph_options}}}
\examples{
\dontrun{

# plotting a simple ring graph, all default parameters, except the layout
g <- make_ring(10)
g$layout <- layout_in_circle
plot(g)
tkplot(g)
rglplot(g)

# plotting a random graph, set the parameters in the command arguments
g <- barabasi.game(100)
plot(g, layout=layout_with_fr, vertex.size=4,
     vertex.label.dist=0.5, vertex.color="red", edge.arrow.size=0.5)

# plot a random graph, different color for each component
g <- sample_gnp(100, 1/100)
comps <- components(g)$membership
colbar <- rainbow(max(comps)+1)
V(g)$color <- colbar[comps+1]
plot(g, layout=layout_with_fr, vertex.size=5, vertex.label=NA)

# plot communities in a graph
g <- make_full_graph(5) \%du\% make_full_graph(5) \%du\% make_full_graph(5)
g <- add_edges(g, c(1,6, 1,11, 6,11))
com <- cluster_spinglass(g, spins=5)
V(g)$color <- com$membership+1
g <- set_graph_attr(g, "layout", layout_with_kk(g))
plot(g, vertex.label.dist=1.5)

# draw a bunch of trees, fix layout
igraph_options(plot.layout=layout_as_tree)
plot(make_tree(20, 2))
plot(make_tree(50, 3), vertex.size=3, vertex.label=NA)
tkplot(make_tree(50, 2, mode="undirected"), vertex.size=10,
vertex.color="green")
}
}
\keyword{graphs}
