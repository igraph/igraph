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

 #### Games
 - `igraph_hsbm_game` for a hierarchical stochastic block model
 - `igraph_hsbm_list_game` for a more general hierarchical stochastic block model
 - `igraph_correlated_game` generates pairs of correlated random graphs by perturbing existing adjacency matrix
 - `igraph_correlated_pair_game` generates pairs of correlated random graphs
 - `igraph_tree_game` generates a random tree with the given number of nodes
 - `igraph_dot_product_game` generates a random dot product graph
 - `igraph_sample_sphere_surface` samples points uniformly from the surface of a sphere
 - `igraph_sample_sphere_volume` samples points uniformly from the volume of a sphere
 - `igraph_sample_dirichlet` samples points from a Dirichlet distribution

#### Other

 - `igraph_simplify_and_colorize()` to encode edge and self-loop multiplicities into edge and vertex colors
 - `igraph_bridges()` to find edges whose removal would disconnect a graph
 - `igraph_realize_degree_sequence()` to create a single graph with a given degree sequence (Havel-Hakimi algorithm)
 - `igraph_vertex_coloring_greedy()` computes a vertex coloring using a greedy algorithm
 - `igraph_rewire_directed_edges()` randomly rewires only the starting points or only the endpoints of directed edges
- `igraph_adjacency_spectral_embedding` and `igraph_laplacian_spectral_embedding` provide graph embedddings and `igraph_dim_select` provides dimensionality selection for singular values using profile likelihood.
- Various `igraph_local_scan_*` functions provide local counts and statistics of neighborhoods.
- `igraph_solve_lsap` solves the linear sum assignment problem
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
 - Provide proper support for Windows, using `__declspec(dllexport)` and `__declspec(dllimport)` for `DLL`s and static usage by using `#define IGRAPH_STATIC 1`.
 - Provided integer versions of `dqueue` and `stack` data types.