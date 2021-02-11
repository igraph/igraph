# igraph C library changelog

## [Unreleased]

### Added

 - Eulerian paths/cycles (PR #1346):
   * `igraph_is_eulerian()` finds out whether an Eulerian path/cycle exists.
   * `igraph_eulerian_path()` returns an Eulerian path.
   * `igraph_eulerian_cycle()` returns an Eulerian cycle.
 - Efficiency (PR #1344):
   * `igraph_global_efficiency()` computes the global efficiency of a network.
   * `igraph_local_efficiency()` computes the local efficiency around each vertex.
   * `igraph_average_local_efficiency()` computes the mean local efficiency.
 - Degree sequences (PR #1445):
   * `igraph_is_graphical()` checks if a degree sequence has a realization as a simple or multigraph, with or without self-loops.
   * `igraph_is_bigraphical()` checks if two degree sequences have a realization as a bipartite graph.
   * `igraph_realize_degree_sequence()` now supports constructing non-simple graphs as well.
 - There is a new fatal error handling mechanism (PR #1548):
   * `igraph_set_fatal_handler()` sets the fatal error handler. It is the only function in this functionality group that is relevant to end users.
   * The macro `IGRAPH_FATAL()` and the functions `igraph_fatal()` and `igraph_fatalf()` raise a fatal error. These are for internal use.
   * `IGRAPH_ASSERT()` is a replacement for the `assert()` macro. It is for internal use.
   * `igraph_fatal_handler_abort()` is the default fatal error handler.
 - The new `IGRAPH_WARNINGF`, `IGRAPH_ERRORF` and `IGRAPH_FATALF` macros provide warning/error reporting with `printf`-like syntax. (PR #1627, thanks to Daniel Noom!)
 - `igraph_average_path_length_dijkstra()` computes the mean shortest path length in weighted graphs (PR #1344).
 - `igraph_is_same_graph()` cheks that two labelled graphs are the same (PR #1604).
 - Harmonic centrality (PR #1583):
   * `igraph_harmonic_centrality()` computes the harmonic centrality of vertices.
   * `igraph_harmonic_centrality_cutoff()` computes the range-limited harmonic centrality.
 - Range-limited centralities, currently equivalent to the old functions with names ending in `_estimate` (PR #1583):
   * `igraph_closeness_cutoff()`.
   * `igraph_betweenness_cutoff()`.
   * `igraph_edge_betweenness_cutoff()`.
 - `igraph_vector_is_any_nan()` checks if any elements of an `igraph_vector_t` is NaN.
 - `igraph_inclist_size()` returns the number of vertices in an incidence list.
 - `igraph_lazy_adjlist_size()` returns the number of vertices in a lazy adjacency list.
 - `igraph_lazy_inclist_size()` returns the number of vertices in a lazy incidence list.
 - `igraph_bfs_simple()` now provides a simpler interface to the breadth-first search functionality.

### Changed

 - igraph now uses a CMake-based build sysyem.
 - GMP support can no longer be disabled. When GMP is not present on the system, igraph will use an embedded copy of Mini-GMP (PR #1549).
 - Bliss has been updated to version 0.75. Bliss functions are now interruptible. Thanks to Tommi Junttila for making this possible!
 - Adjacency and incidence lists:
   * `igraph_adjlist_init()` and `igraph_lazy_adjlist_init()` now require the caller to specify what to do with loop and multiple edges.
   * `igraph_inclist_init()` and `igraph_lazy_inclist_init()` now require the caller to specify what to do with loop edges.
   * Adjacency and incidence lists now use `igraph_vector_int_t` consistently.
 - Community detection:
   * `igraph_community_multilevel()`: added resolution parameter.
   * `igraph_community_fluid_communities()`: graphs with no vertices or with one vertex only are now supported; they return a trivial partition.
 - Modularity:
   * `igraph_modularity()` and `igraph_modularity_matrix()`: added resolution parameter.
   * `igraph_modularity()` and `igraph_modularity_matrix()` now support the directed version of modularity.
   * `igraph_modularity()` returns NaN for graphs with no edges to indicate that the modularity is not well-defined for such graphs.
 - Centralities:
   * `cutoff=0` is no longer interpreted as infinity (i.e. no cutoff) in `betweenness`, `edge_betweenness` and `closeness`. If no cutoff is desired, use a negative value such as `cutoff=-1`.
   * The `nobigint` argument has been removed from `igraph_betweenness()`, `igraph_betweenness_estimate()` and `igraph_centralization_betweenness()`, as it is not longer needed. The current implementation is more accurate than the old one using big integers.
   * `igraph_closeness()` now considers only reachable vertices during the calculation (i.e. the closeness is calculated per-component in the undirected case) (PR #1630).
   * `igraph_closeness()` gained two additional output parameters, `reachable_count` and `all_reachable`, returning the number of reached vertices from each vertex, as well as whether all vertices were reachable. This allows for computing various generalizations of closeness for disconnected graphs (PR #1630).
   * `igraph_pagerank()`, `igraph_personalized_pagerank()` and `igraph_personalized_pagerank_vs()` no longer support the `IGRAPH_PAGERANK_ALGO_POWER` method. Their `options` argument now has type `igraph_arpack_options_t *` instead of `void *`.
 - Shortest paths (PR #1344):
   * `igraph_average_path_length()` now returns the number of disconnected vertex pairs in the new `unconn_pairs` output argument.
   * `igraph_diameter()` now return the result as an `igraph_real_t` instead of an `igraph_integer_t`.
   * `igraph_average_path_length()`  and `igraph_diameter()` now return `IGRAPH_INFINITY` when `unconn=FALSE` and the graph is not connected. Previously they returned the number of vertices.
 - `igraph_subisomorphic_lad()` now supports graphs with self-loops.
 - `igraph_is_chordal()` and `igraph_maximum_cardinality_search()` now support non-simple graphs and directed graphs.
 - `igraph_realize_degree_sequence()` has an additional argument controlling whether multi-edges or self-loops are allowed.   
 - `igraph_is_connected()` now returns false for the null graph; see https://github.com/igraph/igraph/issues/1538 for the reasoning behind this decision.
 - `igraph_lapack_ddot()` is renamed to `igraph_blas_ddot()`.
 - `igraph_to_directed()`: added RANDOM and ACYCLIC modes (PR #1511).
 - `igraph_topological_sorting()` now issues an error if the input graph is not acyclic. Previously it issued a warning.
 - `igraph_vector_(which_)(min|max|minmax)()` now handles NaN elements.
 - `igraph_i_set_attribute_table()` is renamed to `igraph_set_attribute_table()`.
 - `igraph_i_sparsemat_view()` is renamed to `igraph_sparsemat_view()`.

### Deprecated

 - `igraph_is_degree_sequence()` and `igraph_is_graphical_degree_sequence()` are deprecated in favour of the newly added `igraph_is_graphical()`.
 - `igraph_closeness_estimate()` is deprecated in favour of the newly added `igraph_closeness_cutoff()`.
 - `igraph_betweenness_estimate()` and `igraph_edge_betweenness_estimate()` are deprecated in favour of the newly added `igraph_betweenness_cutoff()` and `igraph_edge_betweenness_cutoff()`.
 - `igraph_adjlist_remove_duplicate()` and `igraph_inclist_remove_duplicate()` are now deprecated in favour of the new constructor arguments in `igraph_adjlist_init()` and `igraph_inclist_init()`.

### Removed

 - The following functions, all deprecated in igraph 0.6, have been removed (PR #1562):
   * `igraph_adjedgelist_init()`, `igraph_adjedgelist_destroy()`, `igraph_adjedgelist_get()`, `igraph_adjedgelist_print()`, `igraph_adjedgelist_remove_duplicate()`.
   * `igraph_lazy_adjedgelist_init()`, `igraph_lazy_adjedgelist_destroy()`, `igraph_lazy_adjedgelist_get()`, `igraph_lazy_adjedgelist_get_real()`.
   * `igraph_adjacent()`.
   * `igraph_es_adj()`.
   * `igraph_subgraph()`.
 - `igraph_pagerank_old()`, deprecated in 0.7, has been removed.

### Fixed

 - Betweenness calculations are no longer at risk from integer overflow.
 - The actual cutoff distance used in closeness calculation was one smaller than the `cutoff` parameter. This is corrected (PR #1630).
 - `igraph_layout_gem()` was not interruptible; now it is.
 - `igraph_barabasi_aging_game()` now checks its parameters more carefully.
 - `igraph_callaway_traits_game()` now checks its parameters.
 - `igraph_lastcit_game()` checks its parameters more carefully, and no longer crashes with zero vertices (PR #1625).
 - `igraph_residual_graph()` now returns the correct _residual_ capacities; previously it wrongly returned the original capacities (PR #1598).
 - `igraph_psumtree_update()` now checks for negative values and NaN.
 - `igraph_communities_spinglass()`: fixed several memory leaks in the `IGRAPH_SPINCOMM_IMP_NEG` implementation.
 - `igraph_incident()` now returns edges in the same order as `igraph_neighbors()`.
 - `igraph_modularity_matrix()` returned incorrect results for weighted graphs. This is now fixed. (PR #1649, thanks to Daniel Noom!)
 - PageRank (PR #1640):
   * `igraph_(personalized_)pagerank(_vs)()` now check their parameters more carefully.
   * `igraph_personalized_pagerank()` no longer modifies its `reset` parameter.
   * `igraph_(personalized_)pagerank(_vs)`: the `IGRAPH_PAGERANK_ALGO_ARPACK` method now handles self-loops correctly.
   * `igraph_personalized_pagerank(_vs)()`: the result retuned for edgeless or all-zero-weight graphs with the `IGRAPH_PAGERANK_ALGO_ARPACK` ignored the personalization vector. This is now corrected.
   * `igraph_personalized_pagerank(_vs)()` with a non-uniform personalization vector, a disconnected graph and the `IGRAPH_PAGERANK_ALGO_PRPACK` method would return results that were inconsistent with `IGRAPH_PAGERANK_ALGO_ARPACK`. This happened because PRPACK always used a uniform reset distribution when the random walk got stuck in a sink vertex. Now it uses the user-specified reset distribution for this case as well.
 - Fixed crashes in several functions when passing a weighted graph with zero edges (due to `vector_min` being called on the zero-length weight vector).
 - Weighted betweenness, closeness, PageRank and shortest path calculations, as well as random walk functions now check if any weights are NaN.
 - Compatibility with the PGI compiler.

### Other

 - Documentation improvements.
 - Improved error and warning messages.
 - More robust error handling.
 - General code cleanup to reduce the number of compiler warnings.
 - igraph's source files have been re-organized for better maintainability.
 - igraph can now be built with an external CXSparse library.
 - The references to igraph source files in error and warning messages are now always relative to igraph's base directory.
 - When igraph is built as a shared library, only public symbols are exported even on Linux and macOS.

### Acknowledgments

 - Thanks to Daniel Noom for significantly expanding igraph's test coverage and exposing several issues in the process!

## [0.8.5] - 2020-12-07

### Changed

 - `igraph_write_graph_pajek()`: the function now always uses the platform-native line endings (CRLF on Windows, LF on Unix and macOS). Earlier versions tried to enforce Windows line endings, but this was error-prone, and since all recent versions of Pajek support both line endings, enforcing Windows line endings is not necessary any more.

### Fixed

 - Fixed several compilation issues with MINGW32/64 (PR #1554)
 - `igraph_layout_davidson_harel()` was not interruptible; now it is.
 - Added a missing memory cleanup call in `igraph_i_cattribute_combine_vertices()`.
 - Fixed a few memory leaks in test cases.

## [0.8.4] - 2020-11-24

### Fixed

 - `igraph_i_cattribute_combine_vertices()`: fixed invalid cleanup code that eventually filled up the "finally" stack when combining vertices with attributes extensively.
 - `igraph_hrg_sample()`: fixed incorrect function prototype
 - `igraph_is_posinf()` and `igraph_is_neginf()`: fixed incorrect result on platforms where the sign of the result of `isinf()` is not indicative of the sign of the input.
 - Fixed building with vendored LAPACK and external BLAS
 - Fixed building with XCode 12.2 on macOS

### Other

 - Documentation improvements
 - General code cleanup to reduce the number of compiler warnings

## [0.8.3] - 2020-10-02

### Added

 - `igraph_vector_binsearch_slice()` performs binary search on a sorted slice of a vector.

### Changed

 - `igraph_eigenvector_centrality()` assumes the adjacency matrix of undirected graphs to have twice the number of self-loops for each vertex on the diagonal. This makes the results consistent between an undirected graph and its directed equivalent when each edge is replaced by a mutual edge pair.

### Fixed

 - `igraph_isomorphic()` now verifies that the input graphs have no multi-edges (PR #1464).
 - `igraph_difference()` was creating superfluous self loops (#597).
 - `igraph_count_multiple()` was giving incorrect results for self-loops in directed graph (PR #1399).
 - `igraph_betweenness_estimate()`: fixed incorrect results with finite cutoff (PR #1392).
 - `igraph_count_multiple()` was giving incorrect results for self-loops in directed graph (PR #1399).
 - `igraph_eigen_matrix_symmetric()`: fixed incorrect matrix multiplication (PR #1379).
 - Corrected several issues that could arise during an error condition (PRs #1405, #1406, #1438).
 - `igraph_realize_degree_sequence()` did not correctly detect some non-graphical inputs.
 - `igraph_is_graphical_degree_sequence()`: fixed incorrect results in undirected case (PR #1441).
 - `igraph_community_leiden()`: fixed incorrect result when self-loops are present (PR #1476).
 - `igraph_eigenvector_centrality()`: fixed incorrect value for isolated vertices in weighted graphs.
 - `igraph_eigenvector_centrality()`: corrected the handling of self-loops.
 - `igraph_layout_reingold_tilford()`: fixed an issue where branches of the tree would sometimes overlap.

### Other

 - `igraph_degree_sequence_game()`: improved performance with `IGRAPH_DEGSEQ_SIMPLE_NO_MULTIPLE_UNIFORM` method.
 - Improved the robustness of the test suite.
 - Documentation improvements.
 - Improved error and warning messages.
 - Improved compatibility with recent versions of Microsoft Visual C.

## [0.8.2] - 2020-04-28

### Changed

 - Improved argument checking: `igraph_all_st_mincuts()` and `igraph_sir()`
 - Improved interruptibility: `igraph_sir()`

### Fixed

 - `igraph_community_leiden()`: fixed crash when interrupting
 - The tests are now more robust. Some incorrect test failures were fixed when
   running on i386 architecture, or when using different versions of external
   dependencies.

### Other

 - Improved error messages from `igraph_sir()`.
 - Improved compatibility with more recent versions of Microsoft Visual C.

## [0.8.1] - 2020-03-13

### Changed

 - Improved interruptability: `igraph_degree_sequence_game()`
 - Improved argument checking: `igraph_forest_fire_game()`
 - Updated the plfit library to version 0.8.1

### Fixed

 - `igraph_community_edge_betweenness()`: fix for graphs with no edges (PR #1312)
 - `igraph_bridges()` now handles multigraphs correctly (PR #1335)
 - `igraph_avg_nearest_neighbor_degree()`: fix for memory leak in weighted case (PR #1339)
 - `igraph_community_leiden()`: fix crash bug (PR #1357)

### Other

 - Included `ACKOWLEDGEMENTS.md`
 - Documentation improvements

## [0.8.0] - 2020-01-29

### Added

 * Trees

   - `igraph_to_prufer()` and `igraph_from_prufer()` convert labelled trees to/from Prüfer sequences
   - `igraph_tree_game()` samples uniformly from the set of labelled trees
   - `igraph_is_tree()` checks if a graph is a tree
   - `igraph_random_spanning_tree()` picks a spanning tree of a graph uniformly at random
   - `igraph_random_edge_walk()` returns the indices of edges traversed by a random walk; useful for multigraphs

 * Community detection

   - `igraph_community_fluid_communities()` detects communities based on interacting fluids
   - `igraph_community_leiden()` detects communities with the Leiden method

 * Cliques

   - `igraph_maximal_cliques_hist()` counts maximal cliques of each size
   - `igraph_maximal_cliques_callback()` calls a function for each maximal clique
   - `igraph_clique_size_hist()` counts cliques of each size
   - `igraph_cliques_callback()` calls a function for each clique
   - `igraph_weighted_cliques()` finds weighted cliques in graphs with integer vertex weights
   - `igraph_weighted_clique_number()` computes the weighted clique number
   - `igraph_largest_weighted_cliques()` finds the largest weighted cliques

 * Graph generators

   - `igraph_hsbm_game()` for a hierarchical stochastic block model
   - `igraph_hsbm_list_game()` for a more general hierarchical stochastic block model
   - `igraph_correlated_game()` generates pairs of correlated random graphs by perturbing existing adjacency matrix
   - `igraph_correlated_pair_game()` generates pairs of correlated random graphs
   - `igraph_tree_game()` samples uniformly from the set of labelled trees
   - `igraph_dot_product_game()` generates a random dot product graph
   - `igraph_realize_degree_sequence()` creates a single graph with a given degree sequence (Havel-Hakimi algorithm)

 * Graph embeddings

   - `igraph_adjacency_spectral_embedding()` and `igraph_laplacian_spectral_embedding()` provide graph embedddings
   - `igraph_dim_select()` provides dimensionality selection for singular values using profile likelihood

 * Isomorphism

   - `igraph_automorphism_group()` computes the generators of the automorphism group of a simple graph
   - `igraph_simplify_and_colorize()` encodes edge and self-loop multiplicities into edge and vertex colors; use in conjunction with VF2 to test isomorphism of non-simple graphs

 * Other

   - `igraph_bridges()` finds edges whose removal would disconnect a graph
   - `igraph_vertex_coloring_greedy()` computes a vertex coloring using a greedy algorithm
   - `igraph_rewire_directed_edges()` randomly rewires only the starting points or only the endpoints of directed edges
   - Various `igraph_local_scan_*` functions provide local counts and statistics of neighborhoods
   - `igraph_sample_sphere_surface()` samples points uniformly from the surface of a sphere
   - `igraph_sample_sphere_volume()` samples points uniformly from the volume of a sphere
   - `igraph_sample_dirichlet()` samples points from a Dirichlet distribution
   - `igraph_malloc()`, to be paired with the existing `igraph_free()`

### Changed

 - `igraph_degree_sequence_game()`: new method added for uniform sampling: `IGRAPH_DEGSEQ_SIMPLE_NO_MULTIPLE_UNIFORM`
 - `igraph_modularity_matrix()`: removed `membership` argument (PR #1194)
 - `igraph_avg_nearest_neighbor_degree()`: added `mode` and `neighbor_degree_mode` arguments (PR #1214).
 - `igraph_get_all_simple_paths()`: added `cutoff` argument (PR #1232).
 - `igraph_unfold_tree()`: no longer preserves edge ordering of original graph
 - `igraph_decompose()`: support strongly connected components
 - `igraph_isomorphic_bliss()`, `igraph_canonical_permutation()`, `igraph_automorphisms()`: added additional arguments to support vertex colored graphs (PR #873)
 - `igraph_extended_chordal_ring`: added argument to support direction (PR #1096), and fixed issue #1093.

### Other

 - The [Bliss isomorphism library](http://www.tcs.hut.fi/Software/bliss/) was updated to version 0.73. This version adds support for vertex colored and directed graphs.
 - igraph now uses the high-performance [Cliquer library](https://users.aalto.fi/~pat/cliquer.html) to find (non-maximal) cliques
 - Provide proper support for Windows, using `__declspec(dllexport)` and `__declspec(dllimport)` for `DLL`s and static usage by using `#define IGRAPH_STATIC 1`.
 - Provided integer versions of `dqueue` and `stack` data types.
