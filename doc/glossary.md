<!--
  The glossary is generated from these Markdown sources using
    pandoc glossary.md --to docbook > glossary.xml
  and manually updated to fit into the documentation system.
-->

# Glossary

This glossary defines common terms used throughout the igraph documentation.

 - **attribute**: A piece of data associated with a vertex, an edge, or the graph itself. The igraph C library currently supports numeric, string and Boolean attribute values, and provides a means for implementing attribute handlers that support custom types.
 - **adjacent**: Two vertices are called **adjacent** if there is an edge connecting them. This term describes a vertex-to-vertex relation.
 - **adjacency list**: A data structure that associates a list of neighbours (i.e. adjacent vertices) to each vertex.
 - **adjacency matrix**: A representation of a graph as a square matrix. `A_ij` gives the number of edge endpoints connecting from the `i`th vertex to the `j`th vertex. Conventionally, the diagonal of the adjacency matrix of an undirected graph contains _twice_ the number of self-loops. All igraph functions follow this convention unless noted otherwise.
 - **biadjacency matrix**: Analogous to the adjacency matrix, but used for bipartite graphs. Element `B_ij` gives the number of edges from the `i`th vertex of the first group to the `j`th vertex of the second group.
 - **bipartite graph**: A graph whose vertices can be partitioned into two groups in such a way that connections are present only between members of different groups.
 - **complete graph**: Also called **full graph** within the context of igraph, a graph in which all pairs of vertices are connected to each other.
 - **connected graph**: A connected graph consists of a single component, in which any vertex is reachable from any other. In igraph, the null graph is not considered connected, as it has not one, but zero components.
 - **edge**: A **connection** between two vertices, also called a **link**. In igraph, edges are referred to by integer indices called **edge IDs**.
 - **finalizer stack**: A global stack used internally by igraph to keep track of currently allocated objects and their destructors, so that they can be automatically destroyed in case of an error.
 - **game**: Within igraph, this term is used for stochastic graph generators, i.e. functions that sample from random graph models.
 - **graph** or **network**: A set of vertices with connections between them. In igraph, graphs may carry associated data in the form of vertex, edge or graph attributes.
 - **incident**: An edge is called **incident** to the vertices that are its endpoints. This term describes a vertex-to-edge relation.
 - **incidence list**: A data structure that associates a list of incident edges to each vertex.
 - **incidence matrix**: A matrix describing the incidence relation between vertices (rows) and edges (columns).
 - **membership vector**: Membership vectors are a means of encoding a partitioning of items, usually vertices, into several groups. The `i`th element of the vector gives an integer identifier of the group the `i`th vertex belongs to. Membership vectors are typically used to describe a vertex clustering obtained through community detection, or by identifying the connected components of a graph.
 - **multi-edges** or **parallel edges**: More than one edge connecting the same two vertices. In a directed graph, `a -> b, a -> b` are considered parallel edges, but `a -> b, a <- b` are not.
 - **null graph**: A graph with no vertices (and no edges).
 - **self-loop**, **self-edge**, or simply **loop**: An edge that connects a vertex to itself.
 - **simple graph**: A graph that does not have self-loops or multi-edges.
 - **singleton graph**: A graph having a single vertex. This term usually refers to a single vertex with no edges, but note that self-loops may in principle be present.
 - **vertex**: Graphs consist of vertices, also called **nodes**, that are connected to each other. In igraph, vertices are referred to by integer indices called **vertex IDs**.
