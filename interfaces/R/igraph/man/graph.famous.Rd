\encoding{UTF-8}
\name{graph.famous}
\alias{graph.famous}
\title{Creating named graphs}
\description{There are some famous, named graphs, sometimes
  counterexamples to some conjecture or unique graphs with given
  features. These can be created with this function}
\usage{
graph.famous(name)
}
\arguments{
  \item{name}{Character constant giving the name of the graph. It is
    case insensitive.}
}
\details{
  \code{graph.famous} knows the following graphs:
  \itemize{
    \item{Bull}{The bull graph, 5 vertices, 5 edges, resembles to the
      head of a bull if drawn properly.}
    \item{Chvatal}{This is the smallest triangle-free graph that is
      both 4-chromatic and 4-regular. According to the Grunbaum
      conjecture there exists an m-regular, m-chromatic graph
      with n vertices for every m>1 and n>2. The Chvatal graph
      is an example for m=4 and n=12. It has 24 edges.}
    \item{Coxeter}{A non-Hamiltonian cubic symmetric graph with 28
      vertices and 42 edges.}
    \item{Cubical}{The Platonic graph of the cube. A convex regular
      polyhedron with 8 vertices and 12 edges.}
    \item{Diamond}{A graph with 4 vertices and 5 edges, resembles to a
      schematic diamond if drawn properly.}
    \item{Dodecahedral, Dodecahedron}{Another Platonic solid
      with 20 vertices and 30 edges.}
    \item{Folkman}{The semisymmetric graph with minimum number of
      vertices, 20 and 40 edges. A semisymmetric graph is
      regular, edge transitive and not vertex transitive.}
    \item{Franklin}{This is a graph whose embedding to the Klein
      bottle can be colored with six colors, it is a
      counterexample to the neccessity of the Heawood
      conjecture on a Klein bottle. It has 12 vertices and 18
      edges.}
    \item{Frucht}{The Frucht Graph is the smallest cubical graph
      whose automorphism group consists only of the identity
      element. It has 12 vertices and 18 edges.}
    \item{Grotzsch}{The Grötzsch graph is a triangle-free graph with
      11 vertices, 20 edges, and chromatic number 4. It is named after
      German mathematician Herbert Grötzsch, and its existence
      demonstrates that the assumption of planarity is necessary in
      Grötzsch's theorem that every triangle-free planar
      graph is 3-colorable.}
    \item{Heawood}{The Heawood graph is an undirected graph with 14
      vertices and 21 edges. The graph is cubic, and all cycles in the
      graph have six or more edges. Every smaller cubic graph has shorter
      cycles, so this graph is the 6-cage, the smallest cubic graph of
      girth 6.}
    \item{Herschel}{The Herschel graph is the smallest
      nonhamiltonian polyhedral graph. It is the
      unique such graph on 11 nodes, and has 18 edges.}
    \item{House}{The house graph is a 5-vertex, 6-edge graph, the
      schematic draw of a house if drawn properly, basicly a 
      triangle of the top of a square.}
    \item{HouseX}{The same as the house graph with an X in the square. 5
      vertices and 8 edges.}
    \item{Icosahedral, Icosahedron}{A Platonic solid with 12
      vertices and 30 edges.}
    \item{Krackhardt\_Kite}{A social network with 10 vertices and 18 edges.
      Krackhardt, D. Assessing the Political Landscape:
      Structure, Cognition, and Power in Organizations.
      Admin. Sci. Quart. 35, 342-369, 1990.}
    \item{Levi}{The graph is a 4-arc transitive cubic graph, it has
      30 vertices and 45 edges.}
    \item{McGee}{The McGee graph is the unique 3-regular 7-cage
      graph, it has 24 vertices and 36 edges.}
    \item{Meredith}{The Meredith graph is a quartic graph on 70
      nodes and 140 edges that is a counterexample to the conjecture that
      every 4-regular 4-connected graph is Hamiltonian.}
    \item{Noperfectmatching}{A connected graph with 16 vertices and
      27 edges containing no perfect matching. A matching in a graph
      is a set of pairwise non-adjacent edges; that is, no two edges
      share a common vertex. A perfect matching is a matching
      which covers all vertices of the graph.}
    \item{Nonline}{A graph whose connected components are the 9
      graphs whose presence as a vertex-induced subgraph in a
      graph makes a nonline graph. It has 50 vertices and 72 edges.}
    \item{Octahedral, Octahedron}{Platonic solid with 6
      vertices and 12 edges.} 
    \item{Petersen}{A 3-regular graph with 10 vertices and 15 edges. It is
      the smallest hypohamiltonian graph, ie. it is
      non-hamiltonian but removing any single vertex from it makes it
      Hamiltonian.}
    \item{Robertson}{The unique (4,5)-cage graph, ie. a 4-regular
      graph of girth 5. It has 19 vertices and 38 edges.}
    \item{Smallestcyclicgroup}{A smallest nontrivial graph
      whose automorphism group is cyclic. It has 9 vertices and
      15 edges.}
    \item{Tetrahedral, Tetrahedron}{Platonic solid with 4
      vertices and 6 edges.}
    \item{Thomassen}{The smallest hypotraceable graph,
      on 34 vertices and 52 edges. A hypotracable graph does
      not contain a Hamiltonian path but after removing any
      single vertex from it the remainder always contains a
      Hamiltonian path. A graph containing a Hamiltonian path
      is called tracable.}
    \item{Tutte}{Tait's Hamiltonian graph conjecture states that
      every 3-connected 3-regular planar graph is Hamiltonian.
      This graph is a counterexample. It has 46 vertices and 69
      edges.}
    \item{Uniquely3colorable}{Returns a 12-vertex, triangle-free
      graph with chromatic number 3 that is uniquely
      3-colorable.}
    \item{Walther}{An identity graph with 25 vertices and 31
      edges. An identity graph has a single graph automorphism,
      the trivial one.}
  }
}
\value{A graph object.}
\author{Gabor Csardi \email{csardi@rmki.kfki.hu}}
\seealso{ \code{\link{graph}} can create arbitrary graphs, see also the other
  functions on the its manual page for creating special graphs.
}
\examples{
solids <- list(graph.famous("Tetrahedron"),
               graph.famous("Cubical"),
               graph.famous("Octahedron"),
               graph.famous("Dodecahedron"),
               graph.famous("Icosahedron"))
}
\keyword{graphs}
