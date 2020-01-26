# igraph C library changelog

## [0.8] - 2020-01-28

### Added

#### Trees

 - `igraph_to_prufer()` and `igraph_from_prufer()` for converting trees to/from Prüfer sequences
 - `igraph_tree_game()` to sample uniformly from the set of trees
 - `igraph_is_tree()` to check if a graph is a tree
 - `igraph_random_spanning_tree()` to pick a spanning tree of a graph uniformly at random
 - `igraph_random_edge_walk()` returns the indices of edges traversed by a random walk; useful for multigraphs

#### Community detection

 - `igraph_community_fluid_communities()`, community detection based on interacting fluids
 - `igraph_community_leiden()`, Leiden community detection method

#### Cliques

 - `igraph_maximal_cliques_hist()` to count maximal cliques of each size
 - `igraph_maximal_cliques_callback()` to call a function for each maximal clique
 - `igraph_clique_size_hist()` to count cliques of each size
 - `igraph_cliques_callback()` to call a function for each clique
 - `igraph_weighted_cliques()` to find weighted cliques in graphs with integer vertex weights
 - `igraph_weighted_clique_number()` to find the weighted clique number
 - `igraph_largest_weighted_cliques()` to find the largest weighted cliques

#### Other

 - `igraph_simplify_and_colorize()` to encode edge and self-loop multiplicities into edge and vertex colors
 - `igraph_bridges()` to find edges whose removal would disconnect a graph
 - `igraph_realize_degree_sequence()` to create a single graph with a given degree sequence (Havel-Hakimi algorithm)
 - `igraph_vertex_coloring_greedy()`
 - `igraph_rewire_directed_edges()` randomly rewires only the starting points or only the endpoints of directed edges
 - `igraph_malloc()`, to be paired with the existing `igraph_free()`

### Changed 

 - `igraph_degree_sequence_game()`: new method added for uniform sampling: `IGRAPH_DEGSEQ_SIMPLE_NO_MULTIPLE_UNIFORM`
 - `igraph_modularity_matrix()`: removed `membership` argument (PR #1194)
 - `igraph_avg_nearest_neighbor_degree()`: added `mode` and `neighbor_degree_mode` arguments (PR #1214).
 - `igraph_get_all_simple_paths()`: added `cutoff` argument (PR #1232).
 - `igraph_unfold_tree()`: no longer preserves edge ordering of original graph
 - `igraph_decompose()`: support strongly connected components

### Other

 - The [Bliss library](http://www.tcs.hut.fi/Software/bliss/) was updated to version 0.73
 - igraph now uses the high-performance [Cliquer library](https://users.aalto.fi/~pat/cliquer.html) to find (non-maximal) cliques 
