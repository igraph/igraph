\name{Combining attributes}
\alias{attribute.combination}
\concept{Vertex/edge/graph attributes}
\title{How igraph functions handle attributes when the graph changes}
\description{Many times, when the structure of a graph is modified,
  vertices/edges map of the original graph map to vertices/edges in the
  newly created (modified) graph. For example \code{\link{simplify}}
  maps multiple edges to single edges. igraph provides a flexible
  mechanism to specify what to do with the vertex/edge attributes in
  these cases.}
\details{
  The functions that support the combination of attributes have one or
  two extra arguments called \code{vertex.attr.comb} and/or
  \code{edge.attr.comb} that specify how to perform the mapping of the
  attributes. E.g. \code{\link{contract.vertices}} contracts many
  vertices into a single one, the attributes of the vertices can be
  combined and stores as the vertex attributes of the new graph.

  The specification of the combination of (vertex or edge) attributes
  can be given as
  \enumerate{
    \item{a character scalar,}
    \item{a function object or}
    \item{a list of character scalars and/or function objects.}
  }

  If it is a character scalar, then it refers to one of the predefined
  combinations, see their list below.

  If it is a function, then the given function is expected to perform
  the combination. It will be called once for each new vertex/edge in
  the graph, with a single argument: the attribute values of the
  vertices that map to that single vertex.

  The third option, a list can be used to specify different combination
  methods for different attributes. A named entry of the list
  corresponds to the attribute with the same name. An unnamed entry
  (i.e. if the name is the empty string) of the list specifies the
  default combination method. I.e.
  \preformatted{list(weight="sum", "ignore")}
  specifies that the weight of the new edge should be sum of the weights
  of the corresponding edges in the old graph; and that the rest of the
  attributes should be ignored (=dropped).
}
\section{Predefined combination functions}{
  The following combination behaviors are predefined:
  \describe{
    \item{\sQuote{ignore}}{The attribute is ignored and dropped.}
    \item{\sQuote{sum}}{The sum of the attributes is calculated. This
      does not work for character attributes and works for complex
      attributes only if they have a \code{sum} generic
      defined. (E.g. it works for sparse matrices from the \code{Matrix}
      package, because they have a \code{sum} method.)}
    \item{\sQuote{prod}}{The product of the attributes is
      calculated. This does not work for character attributes and 
      works for complex attributes only if they have a \code{prod}
      function defined.}
    \item{\sQuote{min}}{The minimum of the attributes is calculated and
      returned. For character and complex attributes the standard R
      \code{min} function is used.}
    \item{\sQuote{max}}{The maximum of the attributes is calculated and
      returned. For character and complex attributes the standard R
      \code{max} function is used.}
    \item{\sQuote{random}}{Chooses one of the supplied attribute
      values, uniformly randomly. For character and complex attributes
      this is implemented by calling \code{sample}.}
    \item{\sQuote{first}}{Always chooses the first attribute value. It
      is implemented by calling the \code{head} function.}
    \item{\sQuote{last}}{Always chooses the last attribute value. It is
      implemented by calling the \code{tail} function.}
    \item{\sQuote{mean}}{The mean of the attributes is calculated and
      returned. For character and complex attributes this simply calls
      the \code{mean} function.}
    \item{\sQuote{median}}{The median of the attributes is
      selected. Calls the R \code{median} function for all attribute
      types.}
    \item{\sQuote{concat}}{Concatenate the attributes, using the
      \code{c} function. This results almost always a complex
      attribute.}
  }
}
\section{Specifying combination methods for all graphs}{
  The are two standard igraph parameters that define the default
  behavior when combining vertices and edges: \code{vertex.attr.comb}
  specifies how to combine vertices by default, \code{edge.attr.comb}
  does the same for edges.

  E.g. if you want to drop all vertex attributes when combining
  vertices, you can specify
  \preformatted{igraph.options(vertex.attr.comb="ignore")}

  As another example, if -- when combining edges -- you want to keep the
  mean weight of the edges, concatenate their names into a single
  character scalar, and drop everything else, then use
  \preformatted{igraph.options(edge.attr.comb=list(weight="mean",
    name=toString, "ignore")}
}
\section{Simple and complex attributes}{
  An attribute is simple if (for all vertices/edges) it can be specified
  as an atomic vector. Character and numeric attributes are always
  simple. E.g. a vertex attribute that is a numeric vector of arbitrary
  length for each vertex, is a complex attribute.

  Combination of attributes might turn a complex attribute into
  a single one, and the opposite is possible, too. E.g. when
  contatenating attribute values to form the new attribute value, the
  result will be typically a complex attribute.

  See also examples below.
}
\author{Gabor Csardi \email{csardi.gabor@gmail.com}}
\seealso{\code{\link[igraph]{attributes}} on how to use
  graph/vertex/edges attributes in general. \code{\link{igraph.options}} on
  igraph parameters.}
\examples{
g <- graph( c(1,2, 1,2, 1,2, 2,3, 3,4) )
E(g)$weight <- 1:5

## print attribute values with the graph
igraph.options(print.graph.attributes=TRUE)
igraph.options(print.vertex.attributes=TRUE)
igraph.options(print.edge.attributes=TRUE)

## new attribute is the sum of the old ones
simplify(g, edge.attr.comb="sum")

## collect attributes into a string
simplify(g, edge.attr.comb=toString)

## concatenate them into a vector, this creates a complex
## attribute
simplify(g, edge.attr.comb="concat")

E(g)$name <- letters[seq_len(ecount(g))]

## both attributes are collected into strings
simplify(g, edge.attr.comb=toString)

## harmonic average of weights, names are dropped
simplify(g, edge.attr.comb=list(weight=function(x) length(x)/sum(1/x),
                                name="ignore"))

}
\keyword{graphs}
