# igraph C library changelog

## [0.9.4] - 2021-05-31

### Changed

 - Unweighted transitivity (i.e. clustering coefficient) calculations now ignore multi-edges and edge directions instead of rejecting multigraphs and directed graphs.
 - `igraph_transitivity_barrat()` now returns an error code if the input graph has multiple edges (which is not handled correctly by the implementation yet).

### Fixed

 - `igraph_local_scan_k_ecount()` now handles loops correctly.
 - `igraph_transitivity_avglocal_undirected()` is no longer slower than `igraph_transitivity_local_undirected()`.
 - Worked around an invalid warning issued by Clang 9.0 when compiling with OpenMP.

### Other

 - Documentation improvements.

## [0.9.3] - 2021-05-05

### Added

 - OpenMP is now enabled and used by certain functions (notably PageRank calculation) when the compiler supports it. Set `IGRAPH_OPENMP_SUPPORT=OFF` at configuration time to disable this.

### Fixed

 - `igraph_get_incidence()` no longer reads and writes out of bounds when given a non-bipartite graph, but gives a warning and ignores edges within a part.
 - `igraph_dyad_census()` no longer reports an overflow on singleton graphs, and handles loops and multigraphs correctly. Undirected graphs are handled consistently and will no longer give a warning.
 - `igraph_vector_lex_cmp()` and `igraph_vector_colex_cmp()` dereferenced their arguments only once instead of twice, and therefore did not work with `igraph_vector_ptr_sort()`.
 - `igraph_maximal_cliques_subset()` and `igraph_transitivity_barrat()` corrupted the error handling stack ("finally stack") under some circumstances.
 - CMake package files did not respect `CMAKE_INSTALL_LIBDIR`. This only affected Linux distributions which install into `lib64` or other locations instead of `lib`.
 - The parser sources could not be generated when igraph was in a location that contained spaces in its path.
 - igraph no longer links to the math library (`libm`) when this is not necessary.
 - `_CRT_SECURE_NO_WARNINGS` is now defined during compilation to enable compatibility with UWP.
 - Fixed a compilation issue on MSYS / MinGW when link-time optimization was enabled and the `MSYS Makefiles` CMake generator was used. Some source files in igraph were renamed as a consequence, but these should not affect users of the library.

### Deprecated

 - `igraph_rng_min()` is now deprecated; assume a constant zero as its return value if you used this function in your own code.

### Other

 - Updated the vendored CXSparse library to version 3.2.0

## [0.9.2] - 2021-04-14

### Added

 - CMake package files are now installed with igraph. This allows `find_package(igraph)` to find igraph and detect the appropriate compilation options for projects that link to it.

### Fixed

 - igraph can now be used as a CMake subproject in other CMake-based projects.
 - The documentaton can now be built from the release tarball.
 - Configuration will no longer fail when the release tarball is extracted into a subdirectory of an unrelated git repository.
 - The generated pkg-config file was incorrect when `CMAKE_INSTALL_<dir>` variables were absolute paths.
 - On Unix-like systems, the library name is now `libigraph.so.0.0.0`, as it used to be for igraph 0.8 and earlier.
 - Fixed a return type mismatch in parser sources, and fixed some warnings with recent versions of gcc.
 - Fixed a bug in `igraph_get_shortest_paths_dijkstra()` and `igraph_get_shortest_paths_bellman_ford()` that returned incorrect results for unreachable vertices.

### Other

 - Improved installation instructions and tutorial.

## [0.9.1] - 2021-03-23

### Added

 - `igraph_vector_lex_cmp()` and `igraph_vector_colex_cmp()` for lexicographic
   and colexicographic comparison of vectors. These functions may also be used
   for sorting.

### Changed

 - `igraph_community_multilevel()` is now randomized (PR #1696, thanks to Daniel Noom).

### Fixed

 - CMake settings that controlled the library installation directory name, such as `CMAKE_INSTALL_LIBDIR`, were not respected.
 - Under some conditions, the generated pkg-config file contained an incorrect include directory path.
 - The following functions were not exported from the shared library: `igraph_subcomponent()`, `igraph_stack_ptr_free_all()`, `igraph_stack_ptr_destroy_all()`, `igraph_status_handler_stderr()`, `igraph_progress_handler_stderr()`.
 - Built-in random number generators (`igraph_rngtype_mt19937`, `igraph_rngtype_rand`, `igraph_rngtype_glibc2`) were not exported from the shared library.
 - `igraph_layout_graphopt()` no longer rounds the `spring_length` parameter to an integer.
 - `igraph_get_all_shortest_paths_dijkstra()` no longer modifies the `res` vector's item destructor.
 - `igraph_get_shortest_path_bellman_ford()` did not work correctly when calculating paths to all vertices.
 - `igraph_arpack_rnsolve()` checks its parameters more carefully.
 - `igraph_community_to_membership()` does not crash anymore when `csize` is requested but `membership` is not.
 - `igraph_citing_cited_type_game()`: fixed memory leaks (PR #1700, thanks to Daniel Noom).
 - `igraph_transitivity_undirected()`, `igraph_transitivity_avglocal_undirected()` and `igraph_transitivity_barrat()` no longer trigger an assertion failure when used with the null graph (PRs #1709, #1710).
 - `igraph_(personalized_)pagerank()` would return incorrect results for weighted multigraphs with fewer than 128 vertices when using `IGRAPH_PAGERANK_ALGO_PRPACK`.
 - `igraph_diversity()` now checks its input more carefully, and throws an error when the input graph has multi-edges or is directed.
 - `igraph_shortest_paths_johnson()` would return incorrect results when the `to` argument differed from `from` (thanks to Daniel Noom).
 - `igraph_is_graphical()` would fail to set the result variable for certain special degree sequences in the undirected simple graph case.
 - Non-maximal clique finding functions would sometimes return incomplete results when finding more than 2147483647 (i.e. 2^31 - 1) cliques.
 - GLPK internal errors no longer crash igraph.
 - Fixed some potential memory leaks that could happen on error conditions or when certain functions were interrupted.
 - When testing a DLL build on Windows, the `PATH` was sometimes not set correctly, causing the tests to fail (PR #1692).
 - When compiling from the git repository (as opposed to the release tarball), the build would fail with recent versions of `bison` and `flex`.

### Other

 - Documentation improvements.
 - Much faster documentation builds.
 - Allow using a pre-generated `arith.h` header for f2c when cross-compiling; see the Installation section of the documentation.
 - The `IGRAPH_ENABLE_LTO` build option now supports the `AUTO` value, which uses LTO only if the compiler supports it. Warning: CMake may not always be able to detect that LTO is not fully supported. Therefore, the default setting is `OFF`.
 - The following functions are now interruptible: `igraph_grg_game()`, `igraph_sbm_game()`, `igraph_barabasi_game()`, `igraph_barabasi_aging_game()`.
 - Functions that use GLPK, such as `igraph_feedback_arc_set()` and `igraph_community_optimal_modularity()` are now interruptible.
 - Add support for older versions of Clang that do not recognize the `-Wno-varargs` flag.

### Acknowledgments

 - Big thanks to Daniel Noom for continuing to expand the test suite and discovering and fixing several bugs in the process!

## [0.9.0] - 2021-02-16

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
 - `igraph_get_shortest_paths_bellman_ford()` computes the shortest paths (including the vertex and edge IDs along the paths) using the Bellman-Ford algorithm (PR #1642, thanks to Guy Rozenberg). This makes it possible to calculate the shortest paths on graphs with negative edge weights, which was not possible before with Dijkstra's algorithm.
 - `igraph_get_shortest_path_bellman_ford()` is a wrapper for `igraph_get_shortest_paths_bellman_ford()` for the single path case.
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
 - Trait-based random graph generators:
   * `igraph_callaway_traits_game()` and `igraph_establishment_game()` now have an optional output argument to retrieve the generated vertex types.
   * `igraph_callaway_traits_game()` and `igraph_establishment_game()` now allow omitting the type distribution vector, in which case they assume a uniform distribution.
   * `igraph_asymmetric_preference_game()` now accept a different number of in-types and out-types.
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
 - `igraph_vector_bool` and `igraph_matrix_bool` functions that relied on inequality-comparing `igraph_bool_t` values are removed.

### Fixed

 - Betweenness calculations are no longer at risk from integer overflow.
 - The actual cutoff distance used in closeness calculation was one smaller than the `cutoff` parameter. This is corrected (PR #1630).
 - `igraph_layout_gem()` was not interruptible; now it is.
 - `igraph_barabasi_aging_game()` now checks its parameters more carefully.
 - `igraph_callaway_traits_game()` and `igraph_establishment_game()` now check their parameters.
 - `igraph_lastcit_game()` checks its parameters more carefully, and no longer crashes with zero vertices (PR #1625).
 - `igraph_cited_type_game()` incorrectly rounded the attractivity vector entries to integers.
 - `igraph_residual_graph()` now returns the correct _residual_ capacities; previously it wrongly returned the original capacities (PR #1598).
 - `igraph_psumtree_update()` now checks for negative values and NaN.
 - `igraph_communities_spinglass()`: fixed several memory leaks in the `IGRAPH_SPINCOMM_IMP_NEG` implementation.
 - `igraph_incident()` now returns edges in the same order as `igraph_neighbors()`.
 - `igraph_modularity_matrix()` returned incorrect results for weighted graphs. This is now fixed. (PR #1649, thanks to Daniel Noom!)
 - `igraph_lapack_dgetrf()` would crash when passing `NULL` for its `ipiv` argument (thanks for the fix to Daniel Noom).
 - Some `igraph_matrix` functions would fail to report errors on out-of-memory conditions.
 - `igraph_maxdegree()` now returns 0 for the null graph or empty vector set. Previously, it did not handle this case.
 - `igraph_vector_bool_all_e()` now considers all nonzero (i.e. "true") values to be the same.
 - PageRank (PR #1640):
   * `igraph_(personalized_)pagerank(_vs)()` now check their parameters more carefully.
   * `igraph_personalized_pagerank()` no longer modifies its `reset` parameter.
   * `igraph_(personalized_)pagerank(_vs)`: the `IGRAPH_PAGERANK_ALGO_ARPACK` method now handles self-loops correctly.
   * `igraph_personalized_pagerank(_vs)()`: the result retuned for edgeless or all-zero-weight graphs with the `IGRAPH_PAGERANK_ALGO_ARPACK` ignored the personalization vector. This is now corrected.
   * `igraph_personalized_pagerank(_vs)()` with a non-uniform personalization vector, a disconnected graph and the `IGRAPH_PAGERANK_ALGO_PRPACK` method would return results that were inconsistent with `IGRAPH_PAGERANK_ALGO_ARPACK`. This happened because PRPACK always used a uniform reset distribution when the random walk got stuck in a sink vertex. Now it uses the user-specified reset distribution for this case as well.
 - Fixed crashes in several functions when passing a weighted graph with zero edges (due to `vector_min` being called on the zero-length weight vector).
 - Fixed problems in several functions when passing in a graph with zero vertices.
 - Weighted betweenness, closeness, PageRank, shortest path calculations and random walk functions now check if any weights are NaN.
 - Many functions now reject input arguments containing NaN values.
 - Compatibility with the PGI compiler.

### Other

 - Documentation improvements.
 - Improved error and warning messages.
 - More robust error handling.
 - General code cleanup to reduce the number of compiler warnings.
 - igraph's source files have been re-organized for better maintainability.
 - Debugging aid: When igraph is build with AddressSanitizer, the default error handler prints a stack trace before exiting.
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
