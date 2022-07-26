# igraph C library changelog

## [Unreleased]

### Release notes

This release focuses on infrastructural improvements, stability, and making the igraph interface more consistent, more predictable and easier to use. It contains many API-breaking changes and function renamings, in preparation for a future 1.0 release, at which point the API will become stable. Changes in this direction are likely to continue through a 0.11 release. It is recommended that you migrate your code from 0.9 to 0.10 soon, to make the eventual transition to 1.0 easier.

Some of the highlights are:

 - A consistent use of `igraph_integer_t` for all indices and most integer quantities, both in the API and internally. This type is 64-bit by default on all 64-bit systems, bringing support for very large graphs with more than 2 billion vertices. Previously, vertex and edge indices were often represented as `igraph_real_t`. The move to an `igraph_integer_t` also implies a change from `igraph_vector_t` to `igraph_vector_int_t` in many functions.

 - The random number generation framework has been overhauled. Sampling from the full range of `igraph_integer_t` is now possible. Similarly, the sampling of random reals has been improved to utilize almost the full range of the mantissa of an `igraph_real_t`.

 - There is a new fully memory-managed container type for lists of vectors (`igraph_vector_list_t`), replacing most prevous uses of the non-managed `igraph_vector_ptr_t`.

 - File format readers are much more robust and more tolerant of invalid input.

 - Many improvements to robustness and reliability, made possible by internal refactorings.

### Breaking changes

 - igraph now requires CMake 3.18 or later.

 - In order to facilitate the usage of graphs with more than 2 billion vertices
   and edges, we have made the size of the `igraph_integer_t` data type to be
   32 bits on 32-bit platforms and 64 bits on 64-bit platforms by default. You
   also have the option to compile a 32-bit igraph variant on a 64-bit platform
   by changing the `IGRAPH_INTEGER_SIZE` build variable in CMake to 32.

 - The random number generator interface, `igraph_rng_type_t`, has been overhauled.
   Check the declaration of the type for details.

 - Since `igraph_integer_t` aims to be the largest integer size that is feasible
   on a particular platform, there is no need for generic data types based on
   `long int` any more. The `long` variants of generic data types (e.g.,
   `igraph_vector_long_t`) are therefore removed; you should use the corresponding
   `int` variant instead, whose elements are of type `igraph_integer_t`.

 - Several igraph functions that used to take a `long int` or return a
   `long int` now takes or returns an `igraph_integer_t` instead to make the
   APIs more consistent. Similarly, igraph functions that used `igraph_vector_t`
   for arguments that take or return _integral_ vectors (e.g., vertex or
   edge indices) now take `igraph_vector_int_t` instead. Graph-related functions
   where the API was changed due to this reason are listed below, one by one.

 - Similarly, igraph functions that used to accept the `long` variant of a
   generic igraph data type (e.g., `igraph_vector_long_t`) now take the `int`
   variant of the same data type.

 - Error handlers should no longer perform a `longjmp()`. Doing so will introduce
   memory leaks, as resource cleanup is now done in multiple stages, through
   multiple calls to the error handler. Thus, the error handler should either
   abort execution immediately (as the default handler does), or report the error,
   call `IGRAPH_FINALLY_FREE()`, and return normally.

 - `igraph_add_edges()` now uses an `igraph_vector_int_t` for its `edges`
   parameter.

 - `igraph_all_minimal_st_separators()` now returns the separators in
   an `igraph_vector_int_list_t` containing `igraph_vector_int_t` vectors.

 - `igraph_all_st_cuts()` and `igraph_all_st_mincuts()` now return the cuts in
   an `igraph_vector_int_list_t` containing `igraph_vector_int_t` vectors.

 - `igraph_arpack_unpack_complex()` now uses `igraph_integer_t` for its `nev`
   argument instead of `long int`.

 - `igraph_articulation_points()` now uses an `igraph_vector_int_t` to return
   the list of articulation points, not an `igraph_vector_t`.

 - `igraph_assortativity_nominal()` now accepts vertex types in an `igraph_vector_int_t`
   instead of an `igraph_vector_t`.

 - `igraph_asymmetric_preferennce_game()` now uses an `igraph_vector_int_t` to
   return the types of the nodes in the generated graph.

 - `igraph_automorphism_group()` now returns the generators in an `igraph_vector_int_list_t`
   instead of a pointer vector containing `igraph_vector_t` objects.

 - `igraph_barabasi_game()`, `igraph_barabasi_aging_game()`,
   `igraph_recent_degree_game()` and `igraph_recent_degree_aging_game()` now use
   an `igraph_vector_int_t` for the out-degree sequence of the nodes being
   generated instead of an `igraph_vector_t`.

 - `igraph_bfs()` now takes an `igraph_vector_int_t` for its `roots`,
   `restricted`, `order`, `father`, `pred`, `succ` and `dist` arguments instead
   of an `igraph_vector_t`.

 - `igraph_bfs_simple()` now takes `igraph_vector_int_t` for its `vids`,
   `layers` and `parents` arguments instead of an `igraph_vector_t`.

 - `igraph_bfs_simple()` now returns -1 in `parents` for the root node of the
   traversal, and -2 for unreachable vertices. This is now consistent with other
   functions that return a parent vector.

 - `igraph_biconnected_components()` now uses an `igraph_vector_int_t` to return
   the list of articulation points, not an `igraph_vector_t`. Also, the container
   used for the edges and vertices of the components is now an `igraph_vector_int_list_t`
   instead of a pointer vector containing `igraph_vector_t` objects.

 - `igraph_bipartite_projection()` now uses `igraph_vector_int_t` to return
   `multiplicity1` and `multiplicity2`, not `igraph_vector_t`.

 - `igraph_bridges()` now uses an `igraph_vector_int_t` to return the list of
   bridges, not an `igraph_vector_t`.

 - `igraph_callaway_traits_game()` returns the node types in an `igraph_vector_int_t`
   instead of an `igraph_vector_t`.

 - `igraph_cattribute_list()` now uses `igraph_vector_int_t`
   to return `gtypes`, `vtypes` and `etypes`.

 - `igraph_canonical_permutation()` now uses an `igraph_vector_int_t` for its
   labeling parameter.

 - `igraph_cited_type_game()` now uses an `igraph_vector_int_t` for its
   types parameter.

 - `igraph_citing_cited_type_game()` now uses an `igraph_vector_int_t` for its
   types parameter.

 - `igraph_clique_handler_t` now uses an `igraph_vector_int_t` for its
   `clique` parameter, and must return an `igraph_error_t`. Use `IGRAPH_STOP`
   as the return code to terminate the search prematurely. The vector that the
   handler receives is owned by the clique search routine. If you want to hold
   on to the vector for a longer period of time, you need to make a copy of it
   in the handler. Cliques passed to the callback are marked as `const` as a
   reminder to this change.

 - The `res` parameter of `igraph_cliques()` is now an `igraph_vector_int_list_t`.

 - Callbacks used by `igraph_cliques_callback()` need to be updated to account
   for the fact that the callback does not own the clique passed to it any more;
   the callback needs to make a copy if it wants to hold on to the clique for a
   longer period of time. If the callback does not need to store the clique, it
   does not need to do anything any more, and it must not destroy or free the
   clique.

 - `igraph_closeness()` and `igraph_closeness_cutoff()` now use an
   `igraph_vector_int_t` to return `reachable_count`, not an `igraph_vector_t`.

 - `igraph_cohesive_blocks()` now uses an `igraph_vector_int_t` to return the
   mapping from block indices to parent block indices, and the `cohesion`; also,
   it uses an `igraph_vector_int_list_t` to return the blocks themselves instead of
   a pointer vector of `igraph_vector_t`.

 - `igraph_community_walktrap()`, `igraph_community_edge_betweenness()`,
   `igraph_community_eb_get_merges()`, `igraph_community_fastgreedy()`,
   `igraph_community_to_membership()`, `igraph_le_community_to_membership()`,
   `igraph_community_leading_eigenvector()` now use an `igraph_vector_int_t`
   for their `merges` parameter.

 - `igraph_community_walktrap()` now uses `igraph_integer_t` for its `steps` argument.

 - `igraph_community_infomap()` now uses `igraph_integer_t` for its `nb_trials` argument.

 - `igraph_community_edge_betweenness()` now uses an `igraph_vector_int_t` to
   return the edge IDs in the order of their removal as well as the list of edge
   IDs whose removal broke a single component into two.

 - `igraph_community_fluid_communities()` does not provide the modularity in a
   separate output argument any more; use `igraph_modularity()` to retrieve the
   modularity if you need it.

 - `igraph_community_label_propagation()` now uses an `igraph_vector_int_t` for its
   `initial` parameter. It also takes a `mode` argument that specifies how
   labels should be propagated along edges (forward, backward or ignoring edge
   directions).

 - `igraph_community_label_propagation()` does not provide the modularity in a
   separate output argument any more; use `igraph_modularity()` to retrieve the
   modularity if you need it.

 - `igraph_coreness()` now uses an `igraph_vector_int_t` to return the coreness
   values.

 - `igraph_convex_hull()` now uses an `igraph_vector_int_t` to return the
   indices of the input vertices that were chosen to be in the convex hull.

 - `igraph_correlated_game()` and `igraph_correlated_pair_game()` now take an
   `igraph_vector_int_t` as the permutation vector, not an `igraph_vector_t`.

 - `igraph_create()` now uses an `igraph_vector_int_t` for its
   `edges` parameter.

 - `igraph_create_bipartite()` now uses an `igraph_vector_int_t` for its
   `edges` parameter.

 - `igraph_compose()` now returns the edge maps in an `igraph_vector_int_t`
   instead of an `igraph_vector_t`.

 - `igraph_count_multiple()` now returns the multiplicities in an
   `igraph_vector_int_t` instead of an `igraph_vector_t`.

 - `igraph_decompose()` now uses an `igraph_integer_t` for its `maxcompno` and
   `minelements` arguments instead of a `long int`.

 - `igraph_degree()` now uses an `igraph_vector_int_t` to return the degrees.
   If you need the degrees in a vector containing floating-point numbers
   instead (e.g., because you want to pass them on to some other function that
   takes an `igraph_vector_t`), use `igraph_strength()` instead with a null
   weight vector.

 - `igraph_degree_sequence_game()` now takes degree sequences represented as
   `igraph_vector_int_t` instead of `igraph_vector_t`.

 - `igraph_degseq_t`, used by `igraph_degree_sequence_game()`, uses new names
   for its constants. The old names are deprecated, but retained for compatibility.
   See `igraph_constants.h` to see which new name corresponds to which old one.

 - `igraph_delete_vertices_idx()` now uses `igraph_vector_int_t` vectors to
   return the mapping and the inverse mapping of old vertex IDs to new ones.

 - `igraph_deterministic_optimal_imitation()` now expects the list of strategies
   in an `igraph_vector_int_t` instead of an `igraph_int_t`.

 - `igraph_dfs()` now takes an `igraph_vector_int_t` for its `order`, `order_out`,
   `father` and `dist` arguments instead of an `igraph_vector_t`. Furthermore,
   these vectors will contain -2 for vertices that have not been visited; in
   earlier versions, they used to contain NaN instead. Note that -1 is still
   used in the `father` vector to indicate the root of a DFS tree.

 - `igraph_diameter()` and `igraph_diameter_dijkstra()` now use `igraph_vector_int_t`
   vectors to return the list of vertex and edge IDs in the diameter.

 - `igraph_dnorm()` was removed. This is not really a breaking change as the
   function was never documented, but it was exposed from one of the headers.

 - `igraph_dominator_tree()` now takes an `igraph_vector_int_t` for its
   `dom` and `leftout` arguments instead of an `igraph_vector_t`.

 - `igraph_edges()` now takes an `igraph_vector_int_t` for its
   `edges` argument instead of an `igraph_vector_t`.

 - `igraph_es_multipairs()` was removed; you can use the newly added
   `igraph_es_all_between()` instead.

 - `igraph_establishment_game()` now takes an `igraph_vector_int_t` for its
   `node_type_vec` argument instead of an `igraph_vector_t`.

 - `igraph_eulerian_path()` and `igraph_eulerian_cycle()` now use
   `igraph_vector_int_t` to return the list of edge and vertex IDs participating
   in an Eulerian path or cycle instead of an `igraph_vector_t`.

 - `igraph_feedback_arc_set()` now uses an `igraph_vector_int_t` to return the
   IDs of the edges in the feedback arc set instead of an `igraph_vector_t`.

 - `igraph_get_adjacency()` no longer has the `eids` argument, which would produce
   an adjacency matrix where non-zero values were 1-based (not 0-based) edge IDs.
   If you need a matrix with edge IDs, create it manually.

 - `igraph_get_adjacency_sparse()` now returns the sparse adjacency matrix in
   an `igraph_sparsemat_t` structure, and it assumes that the input matrix is
   _initialized_ for sake of consistency with other igraph functions.

 - `igraph_get_adjacency()` and `igraph_get_adjacency_sparse()` now has a
   `loops` argument that lets the user specify how loop edges should be handled.

 - `igraph_get_edgelist()` now uses an `igraph_vector_int_t` for its
   `res` parameter.

 - `igraph_get_eids()` now uses `igraph_vector_int_t` to return lists of edge IDs
   and to receive lists of vertex IDs.

 - The `path` argument of `igraph_get_eids()` was removed. You can replicate the
   old behaviour by constructing the list of vertex IDs explicitly from the
   path by duplicating each vertex in the path except the first and last ones.
   A helper function called `igraph_expand_path_to_pairs()` is provided to ease
   the transition.

 - `igraph_get_eids_multi()` was removed as its design was fundamentally broken;
   there was no way to retrieve the IDs of all edges between a specific pair of
   vertices without knowing in advance how many such edges there are in the graph.

 - `igraph_get_incidence()` now returns the vertex IDs corresponding to the
   rows and columns of the incidence matrix as `igraph_vector_int_t`.

 - `igraph_get_shortest_path()`, `igraph_get_shortest_path_bellman_ford()` and
   `igraph_get_shortest_path_dijkstra()` now use `igraph_vector_int_t` vectors
   to return the list of vertex and edge IDs in the shortest path.

 - `igraph_get_shortest_paths()`, `igraph_get_shortest_paths_dijkstra()` and
   `igraph_get_shortest_paths_bellman_ford()` now use an `igraph_vector_int_t`
   to return the predecessors and inbound edges instead of an
   `igraph_vector_long_t`.

 - The functions `igraph_get_all_shortest_paths()`,
   `igraph_get_all_shortest_paths_dijkstra()`, `igraph_get_shortest_paths()`,
   `igraph_get_shortest_paths_bellman_ford()` and
   `igraph_get_shortest_paths_dijkstra()` now return paths in an
   `igraph_vector_int_list_t` instead of a pointer vector containing
   `igraph_vector_t` objects.

 - The `maps` parameters in `igraph_get_isomorphisms_vf2()` and
   `igraph_get_subisomorphisms_vf2()` are now of type `igraph_vector_int_list_t`.

 - `igraph_get_stochastic()` now has an additional `weights` argument for edge
   weights.

 - `igraph_get_stochastic_sparse()` now returns the sparse adjacency matrix in
   an `igraph_sparsemat_t` structure, and it assumes that the input matrix is
   _initialized_ for sake of consistency with other igraph functions. It also
   received an additional `weights` argument for edge weights.

 - `igraph_girth()` now uses an `igraph_vector_int_t` for its
   `circle` parameter.

 - The `igraph_vector_ptr_t` cliques parameter in `igraph_graphlets()`
   now contains `igraph_vector_int_t`, not `igraph_vector_t`.

 - The `igraph_vector_ptr_t` cliques parameter in `igraph_graphlets_candidate_basis()`
   now contains `igraph_vector_int_t`, not `igraph_vector_t`.

 - The `igraph_vector_ptr_t` cliques parameter in `igraph_graphlets_project()`
   now contains `igraph_vector_int_t`, not `igraph_vector_t`.

 - `igraph_hrg_init()` and `igraph_hrg_resize()` now takes an `igraph_integer_t`
   as their size arguments instead of an `int`.

 - `igraph_hrg_consensus()` now returns the parent vector in an `igraph_vector_int_t`
   instead of an `igraph_vector_t`.

 - `igraph_hrg_predict()` now uses an `igraph_vector_int_t` for its
   `edges` parameter.

 - `igraph_hrg_sample()` now always samples a single graph only. Use
   `igraph_hrg_sample_many()` if you need more than one sample, and call
   `igraph_hrg_fit()` beforehand if you do not have a HRG model but only a
   single input graph.

 - `igraph_hrg_size()` now returns an `igraph_integer_t` instead of an `int`.

 - `igraph_incident()` now uses an `igraph_vector_int_t` for its
   `eids` parameter.

 - The `res` parameter in `igraph_independent_vertex_sets()` is now an
   `igraph_vector_int_list_t`.

 - `igraph_induced_subgraph_map()` now uses `igraph_vector_int_t` vectors to
   return the mapping and the inverse mapping of old vertex IDs to new ones.

 - `igraph_intersection()` now uses an `igraph_vector_int_t` for its
   `edge_map1` and `edge_map2` parameters.

 - The `edgemaps` parameter of `igraph_intersection_many()` is now an
   `igraph_vector_int_list_t` instead of a pointer vector.

 - `igraph_is_chordal()` now uses an `igraph_vector_int_t` for its
   `alpha`, `alpham1` and `fill_in` parameters.

 - `igraph_is_graphical()` and `igraph_is_bigraphical()` now take degree
   sequences represented as `igraph_vector_int_t` instead of `igraph_vector_t`.

 - `igraph_is_matching()`, `igraph_is_maximal_matching()`,
   `igraph_maximum_bipartite_matching` and `igraph_maximum_matching()` now
   use an `igraph_vector_int_t` to return the matching instead of an
   `igraph_vector_long_t`.

 - The `vids` parameter for `igraph_isoclass_subgraph()` is now an
   `igraph_vector_int_t` instead of `igraph_vector_t`.

 - `igraph_isomorphic_vf2()`, `igraph_get_isomorphisms_vf2_callback()` (which
   used to be called `igraph_isomorphic_function_vf2()`) and
   `igraph_isohandler_t` now all use `igraph_vector_int_t` for their `map12` and
   `map21` parameters.

 - The `igraph_vector_ptr_t` cliques parameter in `igraph_largest_cliques()`
   now contains `igraph_vector_int_t`, not `igraph_vector_t`.

 - The `igraph_vector_ptr_t` res parameter in `igraph_largest_independent_vertex_sets()`
   now contains `igraph_vector_int_t`, not `igraph_vector_t`.

 - The `igraph_vector_ptr_t` res parameter in `igraph_largest_weighted_cliques()`
   now contains `igraph_vector_int_t`, not `igraph_vector_t`.

 - The dimension vector parameter for `igraph_square_lattice()` (used to be
   `igraph_lattice()`) is now an `igraph_vector_int_t` instead of `igraph_vector_t`.

 - The maxiter parameter of `igraph_layout_bipartite()` is now an `igraph_integer_t`
   instead of `long int`.

 - The fixed parameter of `igraph_layout_drl()` and `igraph_layout_drl_3d()`
   was removed as it has never been implemented properly.

 - The width parameter of `igraph_layout_grid()` is now an `igraph_integer_t`
   instead of `long int`.

 - The width and height parameters of `igraph_layout_grid_3d()` are now
   `igraph_integer_t` instead of `long int`.

 - The dimension parameter of `igraph_layout_mds()` is now an `igraph_integer_t`
   instead of `long int`.

 - The `roots` and `rootlevel` parameters of `igraph_layout_reingold_tilford()`
   are now `igraph_vector_int_t` instead of `igraph_vector_t`.

 - The `roots` and `rootlevel` parameters of `igraph_layout_reingold_tilford_circular()`
   are now `igraph_vector_int_t` instead of `igraph_vector_t`.

 - The order parameter of `igraph_layout_star()` is now an `igraph_vector_int_t`
   instead of an `igraph_vector_t`.

 - The maxiter parameter of `igraph_layout_sugiyama()` is now an `igraph_integer_t`
   instead of `long int`. Also, the function now uses an `igraph_vector_int_t`
   for its `extd_to_orig_eids` parameter.

 - The shifts parameter of `igraph_lcf_vector()` is now an `igraph_vector_int_t`
   instead of an `igraph_vector_t`.

 - `igraph_maxflow()` now uses an `igraph_vector_int_t` for its
   `cut`, `partition` and `partition2` parameters.

 - The `igraph_maxflow_stats_t` struct now contains `igraph_integer_t` values instead of `int` ones.

 - The `res` parameters in `igraph_maximal_cliques()` and `igraph_maximal_cliques_subset()`
   are now of type `igraph_vector_int_list_t`.

 - Callbacks used by `igraph_maximal_cliques_callback()` need to be updated to account
   for the fact that the callback does not own the clique passed to it any more;
   the callback needs to make a copy if it wants to hold on to the clique for a
   longer period of time. If the callback does not need to store the clique, it
   does not need to do anything any more, and it must not destroy or free the
   clique.

 - The `res` parameter in `igraph_maximal_independent_vertex_sets()` is now
   an `igraph_vector_int_list_t`.

 - `igraph_maximum_cardinality_search()` now uses an `igraph_vector_int_t` for
   its `alpha` and `alpham1` arguments.

 - `igraph_mincut()` now uses an `igraph_vector_int_t` for its
   `cut`, `partition` and `partition2` parameters.

 - `igraph_moran_process()` now expects the list of strategies in an
   `igraph_vector_int_t` instead of an `igraph_int_t`.

 - Motif callbacks of type `igraph_motifs_handler_t` now take an `igraph_vector_int_t`
   with the vertex IDs instead of an `igraph_vector_t`, and use `igraph_integer_t`
   for the isoclass parameter.

 - Motif functions now use `igraph_integer_t` instead of `int` for their `size`
   parameter.

 - `igraph_neighborhood_size()` now uses an `igraph_vector_int_t` for its
   `res` parameter.

 - The `res` parameter of `igraph_neighborhood()` is now an `igraph_vector_int_list_t`.

 - `igraph_neighbors()` now uses an `igraph_vector_int_t` for its
   `neis` parameter.

 - `igraph_permute_vertices()` now takes an `igraph_vector_int_t` as the
   permutation vector.

 - `igraph_power_law_fit()` does not calculate the p-value automatically any
   more because the previous estimation method did not match the results from
   the original paper of Clauset, Shalizi and Newman (2009) and the
   implementation of the method outlined in the paper runs slower than the
   previous naive estimate. A separate function named
   `igraph_plfit_result_calculate_p_value()` is now provided for calculating
   the p-value. The automatic selection of the `x_min` cutoff also uses a
   different method than earlier versions. As a consequence, results might be
   slightly different if you used tests where the `x_min` cutoff was selected
   automatically. The new behaviour is now consistent with the defaults of the
   underlying `plfit` library.

 - `igraph_preference_game()` now uses an `igraph_vector_int_t` to return the
   types of the nodes in the generated graph.

 - `igraph_random_walk()` now uses an `igraph_vector_int_t` for its
   results. Also, the function now takes both vertices and edges as
   parameters. It can return IDs of vertices and/or edges on the walk.
   The function now takes weights as a parameter to support weighted graphs.

 - `igraph_random_edge_walk()` now uses an `igraph_vector_int_t` for its
   `edgewalk` parameter.

 - `igraph_read_graph_dimacs_flow()` now uses an `igraph_vector_int_t` for its
   label parameter.

 - `igraph_realize_degree_sequence()` now uses an `igraph_vector_int_t` for its
   `outdeg` and `indeg` parameters.

 - `igraph_reindex_membership()` now uses an `igraph_vector_int_t` for its
   `new_to_old` parameter.

 - `igraph_roulette_wheel_imitation()` now expects the list of strategies
   in an `igraph_vector_int_t` instead of an `igraph_int_t`.

 - `igraph_similarity_dice_pairs()` now uses an `igraph_vector_int_t` for its
   `pairs` parameter.

 - `igraph_similarity_jaccard_pairs()` now uses an `igraph_vector_int_t` for its
   `pairs` parameter.

 - `igraph_simple_interconnected_islands_game()` does not generate multi-edges
   between islands any more.

 - `igraph_sort_vertex_ids_by_degree()` and `igraph_topological_sorting()` now
   use an `igraph_vector_int_t` to return the vertex IDs instead of an
   `igraph_vector_t`.

 - `igraph_spanner()` now uses an `igraph_vector_int_t` to return the vector
   of edge IDs in the spanner instead of an `igraph_vector_t`.

 - `igraph_spanning_tree()`, `igraph_minimum_spanning_tree()` and
   `igraph_random_spanning_tree()` now all use an `igraph_vector_int_t` to
   return the vector of edge IDs in the spanning tree instead of an
   `igraph_vector_t`.

 - `igraph_spmatrix_t` and related functions were removed as they mostly
   duplicated functionality that was already present in `igraph_sparsemat_t`.
   Functions that used `igraph_spmatrix_t` in the library now use
   `igraph_sparsemat_t`.

 - `igraph_stochastic_imitation()` now expects the list of strategies
   in an `igraph_vector_int_t` instead of an `igraph_int_t`.

 - `igraph_st_mincut()` now uses an `igraph_vector_int_t` for its
   `cut`, `partition` and `partition2` parameters.

 - `igraph_strvector_get()` now returns strings in the return value, not in an
   output argument.

 - `igraph_subcomponent()` now uses an `igraph_integer_t` for the seed vertex
   instead of an `igraph_real_t`. It also uses an `igraph_vector_int_t` to
   return the list of vertices in the same component as the seed vertex instead
   of an `igraph_vector_t`.

 - `igraph_subisomorphic_vf2()`, `igraph_get_subisomorphisms_vf2_callback()`
   (which used to be called `igraph_subisomorphic_function_vf2()`) and
   `igraph_isomorphic_bliss()` now all use `igraph_vector_int_t` for their `map12`
   and `map21` parameters.

 - The `maps` parameters in `igraph_subisomorphic_lad()`,
   `igraph_get_isomorphisms_vf2()` and `igraph_get_subisomorphisms_vf2()` are
   now of type `igraph_vector_int_list_t`.

 - `igraph_subisomorphic_lad()` now uses an `igraph_vector_int_t` for its `map`
   parameter. Also, its `domains` parameter is now an `igraph_vector_int_list_t`
   instead of a pointer vector containing `igraph_vector_t` objects.

 - `igraph_unfold_tree()` now uses an `igraph_vector_int_t` for its `vertex_index`
   and `roots` parameters.

 - `igraph_union()` now uses an `igraph_vector_int_t` for its
   `edge_map1` and `edge_map2` parameters.

 - The `edgemaps` parameter of `igraph_union_many()` is now an
   `igraph_vector_int_list_t` instead of a pointer vector.

 - `igraph_vector_init_copy()` was refactored to take _another_ vector that the
   newly initialized vector should copy. The old array-based initialization
   function is now called `igraph_vector_init_array()`.

 - `igraph_vector_ptr_init_copy()` was renamed to `igraph_vector_ptr_init_array()`
   for sake of consistency.

 - `igraph_vs_vector()`, `igraph_vss_vector()` and `igraph_vs_vector_copy()` now
   all take an `igraph_vector_int_t` as the vector of vertex IDs, not an
   `igraph_vector_t`. Similarly, `igraph_vs_as_vector()` now returns the vector
   of matched vertex IDs in an `igraph_vector_int_t`, not an `igraph_vector_t`.

 - The `res` parameter of `igraph_weighted_cliques()` is now an
   `igraph_vector_int_list_t`.

 - `igraph_write_graph_dimacs_flow()` now uses `igraph_integer_t` for the source and
   target vertex index instead of a `long int`.

 - `igraph_vector_*()`, `igraph_matrix_*()`, `igraph_stack_*()`, `igraph_array_*()`
   and several other generic igraph data types now use `igraph_integer_t` for
   indexing, _not_ `long int`. Please refer to the headers for the exact details;
   the list of affected functions is too large to include here.

 - `igraph_vector_order()` was removed; use `igraph_vector_int_pair_order()` instead.
   (The original function worked for vectors containing integers only).

 - `igraph_vector_qsort_ind()` and its variants now take an `igraph_order_t` enum
   instead of a boolean to denote whether the order should be ascending or
   descending.

 - The `add_edges()` function in the attribute handler now takes an
   `igraph_vector_int_t` for its `edges` parameter instead of an
   `igraph_vector_t`. The `add_vertices()` function now takes an
   `igraph_integer_t` for the vertex count instead of a `long int`. The
   `combine_vertices()` and `combine_edges()` functions now take an
   `igraph_vector_ptr_t` containing vectors of type `igraph_vector_int_t` in
   their `merges` parameters. The `get_info()` function now uses
   `igraph_vector_int_t` to return the types of the graph, vertex and edge
   attribute types. The `permute_vertices()` and `permute_edges()` functions in
   the attribute handler tables now take an `igraph_vector_int_t` instead of an
   `igraph_vector_t` for the index vectors. These are relevant only to maintainers
   of higher level interfaces to igraph; they should update their attribute
   handlers accordingly.

 - igraph functions that interface with external libraries such as BLAS or LAPACK
   may now fail if the underlying BLAS or LAPACK implementation cannot handle
   the size of input vectors or matrices (BLAS and LAPACK are usually limited to
   vectors whose size fits in an `int`). `igraph_blas_dgemv()` and
   `igraph_blas_dgemv_array()` thus now return an `igraph_error_t`, which may be
   set to `IGRAPH_EOVERFLOW` if the input vectors or matrices are too large.

 - The `values` argument of `igraph_sparsemat_transpose()` was removed; now the
   function always copies the values over to the transposed matrix.

 - `igraph_sparsemat_cholsol()`, `igraph_sparsemat_lusol()`, `igraph_sparsemat_symbqr`
   and `igraph_sparsemat_symblu` now take an `igraph_integer_t` as their `order` parameter.

 - `igraph_sparsemat_count_nonzero()` and `igraph_sparsemat_count_nonzerotol()` now
   return an `igraph_integer_t`.

 - Functions that used an `igraph_vector_t` to represent cluster size
   and cluster membership now use an `igraph_vector_int_t` instead. These are:
   - `igraph_connected_components()` (used to be `igraph_clusters()` in 0.9 and before)
   - `igraph_community_eb_get_merges()`
   - `igraph_community_edge_betweenness()`
   - `igraph_community_fastgreedy()`
   - `igraph_community_fluid_communities()`
   - `igraph_community_infomap()`
   - `igraph_community_label_propagation()`
   - `igraph_community_leading_eigenvector()`
   - `igraph_community_leiden()`
   - `igraph_community_multilevel()`
   - `igraph_community_optimal_modularity()`
   - `igraph_community_spinglass()`
   - `igraph_community_spinglass_single()`
   - `igraph_community_to_membership()`
   - `igraph_community_walktrap()`
   - `igraph_compare_communities()`
   - `igraph_le_community_to_membership()`
   - `igraph_modularity()`
   - `igraph_reindex_membership()`
   - `igraph_split_join_distance()`
   `igraph_community_multilevel()` additionaly uses a `igraph_matrix_int_t`
   instead of `igraph_matrix_t()` for its memberships parameter.

 - `IGRAPH_TOTAL` was removed from the `igraph_neimode_t` enum; use the equivalent
   `IGRAPH_ALL` instead.

 - `igraph_vector_resize_min()` and `igraph_matrix_resize_min()` no longer return an
   error code (return type is now `void`). The vector or matrix is always left in
   a consistent state by these functions, with all data intact, even if releasing
   unused storage is not successful.

 - `igraph_atlas()` now uses `igraph_integer_t` for its `number` argument.

 - `igraph_read_graph_graphml()` now uses `igraph_integer_t` for its `index` argument.

 - `igraph_read_graph_pajek()` now creates a Boolean `type` attribute for bipartite graphs.
   Previously it created a numeric attribute.

 - `igraph_adjacency()` no longer accepts a negative number of edges in its
   adjacency matrix. When negative entries are found, an error is generated.

 - `igraph_adjacency()` gained an additional `loops` argument that lets you
   specify whether the diagonal entries should be ignored or should be interpreted
   as raw edge counts or _twice_ the number of edges (which is common in linear
   algebra contexts).

 - `igraph_weighted_adjacency()` now returns the weights in a separate vector
   instead of storing it in a vertex attribute. The reason is twofold: first,
   the previous solution worked only with the C attribute handler (not the ones
   from the higher-level interfaces), and second, it wasn't consistent with
   other igraph functions that use weights provided as separate arguments.

 - The `loops` argument of `igraph_weighted_adjacency()` was converted to an
   `igraph_loops_t` for sake of consistency with `igraph_adjacency()` and
   `igraph_get_adjacency()`.

 - `igraph_st_vertex_connectivity()` now ignores edges between source and target for `IGRAPH_VCONN_NEI_IGNORE`

 - The `igraph_community_eb_get_merges()` bridges parameter now starts the indices into the
   edge removal vector at 0, not 1.

 - The `igraph_community_eb_get_merges()` now reports an error when not all edges in the graph are
   removed, instead of a nonsensical result.

 - `igraph_incidence()` does not accept negative incidences anymore.

 - `igraph_write_graph_gml()` takes an additional bitfield parameter controlling some aspects of writing
   the GML file.

 - The vector of parents in `igraph_get_shortest_paths()`, `igraph_get_shortest_paths_bellman_ford()`
   and `igraph_get_shortest_paths_dijkstra()` now use -1 to represent the starting vertex, and -2
   for unreachable vertices.

 - `igraph_rng_seed()` now requires an `igraph_uint_t` as its seed arguments. RNG implementations are free to use only the lower bits of the seed if they do not support 64-bit seeds.

 - `igraph_rngtype_rand` (i.e. the RNG that is based on BSD `rand()`) was removed due to poor statistical properties that sometimes resulted in weird artifacts like all-even "random" numbers when igraph's usage patterns happened to line up with the shortcomings of the `rand()` generator in a certain way.

 - The default random number generator has been changed from Mersenne Twister to PCG32.

 - `igraph_vector_minmax()` and `igraph_vector_which_minmax()` no longer return an error code. The return type is now `void`. These functions never fail.

 - `igraph_matrix_minmax()`, `igraph_matrix_which_minmax()`, `igraph_matrix_which_min()` and `igraph_matrix_which_max()` no longer return an error code. The return type is now `void`. These functions never fail.

 - `igraph_is_mutual()` has an additional parameter which controls whether directed self-loops are considered mutual.

### Added

 - A new integer type, `igraph_uint_t` has been added. This is the unsigned pair of `igraph_integer_t` and they are always consistent in size.
 - A new container type, `igraph_vector_list_t` has been added, replacing most uses of `igraph_vector_ptr_t` in the API. It contains `igraph_vector_t` objects, and it is fully memory managed (i.e. its contents do not need to be allocated and destroyed manually). There are specializations for all vector types, such as for `igraph_vector_int_list_t`.
 - `igraph_adjlist_init_from_inclist()` to create an adjacency list from an already existing incidence list by resolving edge IDs to their corresponding endpoints. This function is useful for algorithms when both an adjacency and an incidence list is needed and they should be in the same order.
 - `igraph_vector_*_permute()` functions to permute a vector based on an index vector.
 - `igraph_vector_*_remove_fast()` functions to remove an item from a vector by swapping it with the last element and then popping it off. It allows one to remove an item from a vector in constant time if the order of items does not matter.
 - `igraph_vector_ptr_sort_ind()` to obtain an index vector that would sort a vector of pointers based on some comparison function.
 - `igraph_hub_and_authority_scores()` calculates the hub and authority scores of a graph as a matching pair.
 - `igraph_generalized_petersen()` to create generalized Petersen graphs (#1844, thanks to @alexsyou).
 - `igraph_circulant()` to create circulant graphs (#1856, thanks to @Gomango999).
 - `igraph_symmetric_tree()` to create a tree with the specified number of branches at each level (#1859, thanks to @YuliYudith and @DoruntinaM).
 - `igraph_regular_tree()` creates a regular tree where all internal vertices have the same total degree.
 - `igraph_is_forest()` to check whether a graph is a forest (#1888, thanks to @rohitt28).
 - `igraph_is_acyclic()` to check whether a graph is acyclic (#1945, thanks to @borsgeorgica).
 - `igraph_es_all_between()` to create an edge selector that selects all edges between a pair of vertices.
 - `igraph_blas_dgemm()` to multiply two matrices.
 - `igraph_wheel()` to create a wheel graph (#1938, thanks to @kwofach).
 - `igraph_stack_capacity()` to query the capacity of a stack.
 - `igraph_almost_equals()` and `igraph_cmp_epsilon()` to compare floating point numbers with a relative tolerance.
 - `igraph_complex_almost_equals()` to compare complex numbers with a relative tolerance.
 - `igraph_vector_all_almost_e()`, `igraph_vector_complex_all_almost_e()`, `igraph_matrix_all_almost_e()`, `igraph_matrix_complex_all_almost_e()` for elementwise comparisons of floating point vector and matrices with a relative tolerance.
 - `igraph_vector_range()` to fill an existing vector with a range of increasing numbers.
 - `igraph_roots_for_tree_layout()` computes a set of roots suitable for a nice tree layout.
 - `igraph_fundamental_cycles()` computes a fundamental cycle basis (experimental).
 - `igraph_minimum_cycle_basis()` computes an unweighted minimum cycle basis (experimental).
 - `igraph_strvector_merge()` moves all strings from one string vectors to the end of another without re-allocating them.
 - `igraph_get_k_shortest_paths()` finds the k shortest paths between a source and a target vertex (#1763, thanks to @GroteGnoom)
 - `igraph_get_widest_path()`, `igraph_get_widest_paths()`, `igraph_widest_path_widths_dijkstra()` and `igraph_widest_path_widths_floyd_warshall()` to find widest paths (#1893, thanks to @Gomango999).
 - `igraph_get_laplacian()` and `igraph_get_laplacian_sparse()` return the Laplacian matrix of the graph as a dense or sparse matrix, with various kinds of normalizations. They replace the now-deprecated `igraph_laplacian()`. This makes the API consistent with `igraph_get_adjacency()` and `igraph_get_adjacency_sparse()`.
 - `igraph_enter_safelocale()` and `igraph_exit_safelocale()` for temporarily setting the locale to C. Foreign format readers and writers require a locale which uses a decimal point instead of decimal comma.
 - `igraph_vertex_path_from_edge_path()` converts a sequence of edge IDs representing a path to an equivalent sequence of vertex IDs that represent the vertices the path travelled through.
 - `igraph_graph_count()` gives the number of unlabelled graphs on a given number of vertices. It is meant to find the maximum isoclass value.
 - `igraph_rngtype_pcg32` and `igraph_rngtype_pcg64` implement 32-bit and 64-bit variants of the PCG random number generator.
 - `igraph_rng_get_pois()` generates random variates from the Poisson distribution.
 - `igraph_sparse_adjacency()` and `igraph_sparse_weighted_adjacency()` constructs graphs from (weighted) sparse matrices.
 - `igraph_full_multipartite()` generates full multipartite graphs (a generalization of bipartite graphs to multiple groups).
 - `igraph_turan()` generates Turán graphs.
 - `igraph_has_mutual()` checks if a directed graph has any mutual edges.
 - `igraph_vs_range()`, `igraph_vss_range()`, `igraph_es_range()` and `igraph_ess_range()` creates vertex and edge sequences from C-style intervals (closed from the left, open from the right).

### Removed

 - The type `igraph_stack_ptr_t`, `igraph_stack_ptr_free_all()` and `igraph_stack_ptr_destroy_all()` are removed. Use `igraph_vector_ptr_t` and associated functions instead.

### Changed

 - `igraph_version()` no longer returns an error code.
 - `igraph_write_graph_ncol()` now preserves the edge ordering of the graph when writing an NCOL file.
 - The Pajek parser is now less strict and accepts more files.
 - `igraph_degree_sequence_game()` now supports an additional method, `IGRAPH_DEGSEQ_EDGE_SWITCHING_SIMPLE`,
    and edge-switching MCMC sampler.
 - `igraph_ring()` no longer simplifies its result when generating a one- or two-vertex graph. The one-cycle has a self-loop and the undirected two-cycle has parallel edges.
 - igraph functions that take an ARPACK options object now also accept `NULL` in place of an options object, and they will fall back to using a default object provided by `igraph_arpack_options_get_default()`.
 - `igraph_vector_view()` now allows `data` to be `NULL` in the special case when `length == 0`.
 - `igraph_write_graph_gml()` uses the `creator` parameter in a different way: the supplied
   string is now written into the Creator line as-is instead of being appended to a default
   value.
 - `igraph_read_graph_gml()` now supports graph attributes (in addition to vertex and edge attributes).
 - `igraph_read_graph_gml()` now uses NaN as the default numerical attribute values instead of 0.
   `igraph_write_graph_gml()` skips writing NaN values. These two changes ensure consistent round-tripping.
 - `igraph_write_graph_gml()` and `igraph_read_graph_gml()` now have limited support for entity encoding.
 - Foreign format readers now present more informative error messages.
 - `igraph_get_adjacency()` and `igraph_get_adjacency_sparse()` now counts loop edges _twice_ in undirected graphs when using `IGRAPH_GET_ADJACENCY_BOTH`. This is to ensure consistency with `IGRAPH_GET_ADJACENCY_UPPER` and `IGRAPH_GET_ADJACENCY_LOWER` such that the sum of the upper and the lower triangle matrix is equal to the full adjacency matrix even in the presence of loop edges.
 - It is now possible to overide the uniform integer and the Poisson samplers in the random number generator interface.

### Fixed

 - When an error occurs during parsing DL, GML, NCOL, LGL or Pajek files, line numbers are now reported correctly.
 - The GraphML parser does not print to stderr any more in case of encoding
   errors and other error conditions originating from the underlying `libxml2`
   library.
 - The GraphML parser would omit some edges and vertices when reading files with custom attribute types, such as those produced by yEd. This is now corrected.
 - The GML parser no longer mixes up Inf and NaN and -Inf now works.
 - The GML parser now supports nodes with no id field.
 - The GML parser now performs more stringent checks on the input file, such as verifying that `id`, `source`, `target` and `directed` fields are not duplicated.
 - `igraph_write_graph_gml()` no longer produces corrupt output when some string attribute values contain `"` characters.
 - Graphs no longer lose all their attributes after calling `igraph_contract_vertices()`.
 - `igraph_matrix_complex_create()` and `igraph_matrix_complex_create_polar()` now set their sizes correctly.
 - The core data structures (vector, etc.) have overflow checks now.
 - Deterministic graph generators have overflow checks now.
 - `igraph_sparsemat_getelements_sorted()` did not sort the elements for triplet matrices correctly; this is fixed now.

### Deprecated

 - `igraph_complex_eq_tol()` is now deprecated in favour of `igraph_complex_almost_equals()`.

 - `igraph_clusters()` has been renamed to `igraph_connected_components()`; the
   old name is deprecated and will be removed in 0.11.

 - `igraph_get_sparsemat()` is deprecated in favour of `igraph_get_adjacency_sparse()`,
   and will be removed in 0.11. Note that `igraph_get_adjacency_sparse()` takes an
   _initialized_ sparse matrix as input, unlike `igraph_get_sparsemat()` which takes
   an uninitialized one.

 - `igraph_get_stochastic_sparsemat()` is deprecated in favour of `igraph_get_stochastic_sparse()`,
   and will be removed in 0.11. Note that `igraph_get_stochastic_sparse()` takes an
   _initialized_ sparse matrix as input, unlike `igraph_get_stochastic_sparsemat()` which
   takes an uninitialized one.

 - `igraph_isomorphic_34()` has been deprecated in favour of `igraph_isomorphic()`.
   Note that `igraph_isomorphic()` calls an optimized version for directed graphs
   of size 3 and 4, and undirected graphs with 3-6 vertices, so there is no need
   for a separate function.

 - `igraph_lattice()` has been renamed to `igraph_square_lattice()` to indicate
   that this function generates square lattices only. The old name is deprecated
   and will either be removed in 0.11 or will be changed to become a generic
   lattice generator that also supports other types of lattices.

 - `igraph_laplacian()` is now deprecated; use `igraph_get_laplacian()` or
   `igraph_get_laplacian_sparse()` depending on whether you need a dense or a
   sparse matrix.

 - `igraph_matrix_all_e_tol()` is now deprecated in favour of `igraph_matrix_all_almost_e()`.

 - `igraph_matrix_copy()` is now deprecated; use `igraph_matrix_init_copy()`
   instead. The new name emphasizes that the function _initializes_ the first
   argument instead of expecting an already-initialized target matrix. The old
   name will be removed in 0.11.

 - `igraph_matrix_e()` and `igraph_matrix_e_ptr()` have been renamed to
   `igraph_matrix_get()` and `igraph_matrix_get_ptr()`. The old names are
   deprecated and will be removed in 0.11.

- `igraph_random_edge_walk()` has been deprecated by `igraph_random_walk()`
   to support edges and/or vertices for the random walk in a single function.
   It will be removed in 0.11.

 - `igraph_read_graph_dimacs()` has been renamed to `igraph_read_graph_dimacs_flow()`;
   the old name is deprecated and might be re-used as a generic DIMACS reader
   in the future. Also, the function now uses `igraph_integer_t` as the source
   and target vertex IDs instead of a `long int`.

 - `igraph_shortest_paths()` and related functions were renamed to `igraph_distances()`;
   the old name was unfortunate because these functions calculated _path lengths_
   only and not the paths themselves. The old names are deprecated and will be
   removed in 0.11.

 - `igraph_strvector_add()` has been renamed to `igraph_strvector_push_back()`
   for sake of consistency with other vector-like data structures; the old name
   is deprecated and will be removed in 0.11.

 - `igraph_strvector_copy()` has been renamed to `igraph_strvector_init_copy()`
   for sake of consistency with other vector-like data structures; the old name
   is deprecated and will be removed in 0.11.

 - `igraph_strvector_get()` now returns a `const char*` and not a `char*` to
   indicate that you are not supposed to modify the string in the vector
   directly. If you do want to modify it and you are aware of the implications
   (i.e. the new string must not be longer than the original one), you can
   cast away the constness of the return value before modifying it.

 - `igraph_strvector_set2()` has been renamed to `igraph_strvector_set_len()`;
   the old name is deprecated and will be removed in 0.11.

 - `igraph_tree()` has been renamed to `igraph_kary_tree()`; the old name is
   deprecated and will be removed in 0.11.

 - `igraph_vector_e()` and `igraph_vector_e_ptr()` have been renamed to
   `igraph_vector_get()` and `igraph_vector_get_ptr()`. The old names are
   deprecated and will be removed in 0.11.

 - `igraph_vector_e_tol()` is now deprecated in favour of `igraph_vector_all_almost_e()`.

 - `igraph_vector_copy()` is now deprecated; use `igraph_vector_init_copy()`
   instead. The new name emphasizes that the function _initializes_ the first
   argument instead of expecting an already-initialized target vector. The old
   name will be removed in 0.11.

 - `igraph_vector_init_seq()` is now deprecated in favour of
   `igraph_vector_init_range()`, which uses C-style intervals (closed from the
   left and open from the right).

 - `igraph_write_graph_dimacs()` has been renamed to `igraph_write_graph_dimacs_flow()`;
   the old name is deprecated and might be re-used as a generic DIMACS writer
   in the future. Also, the function now uses `igraph_integer_t` as the source
   and target vertex IDs instead of a `long int`.

 - The macros `igraph_Calloc`, `igraph_Realloc` and `igraph_Free` have been
   deprecated in favour of `IGRAPH_CALLOC`, `IGRAPH_REALLOC` and `IGRAPH_FREE`
   to simplify the API. The deprecated variants will be removed in 0.11.
 - `igraph_vs_seq()`, `igraph_vss_seq()`, `igraph_es_seq()` and `igraph_ess_seq()`
   are now deprecated in favour of `igraph_vs_range()`, `igraph_vss_range()`,
   `igraph_es_range()` and `igraph_ess_range()` because these use C-style
   intervals (closed from the left, open from the right).

### Other

 - Documentation improvements

## [Unreleased 0.9]

### Added

 - `igraph_reverse_edges()` reverses the specified edges in the graph while preserving all attributes.

### Changes

 - The `IGRAPH_ARPACK_PROD` error code is no longer used. Instead, the specific error encountered while doing matrix multiplication is reported.

### Fixed

 - Fixed incorrect results from `igraph_local_scan_1_ecount()` when the graph was directed but the mode was `IGRAPH_ALL` and some nodes had loop edges. See issue #2092.
 - In some rare edge cases, `igraph_pagerank()` with the ARPACK method and `igraph_hub_score()` / `igraph_authority_score()` could return incorrect results. The problem could be detected by checking that the returned eigenvalue is not negative. See issue #2090.
 - `igraph_permute_vertices()` now checks for out-of-range indices and duplicates in the permutation vector.
 - `igraph_create()` now checks for non-finite vertex indices in the edges vector.
 - `igraph_eigenvector_centrality()` would return incorrect scores when some weights were negative.
 - `igraph_es_seq()` and `igraph_ess_seq()` did not include the `to` vertex in the sequence.
 - `igraph_eit_create()` and `igraph_vit_create()` now check that all edge/vertex indices are in range when creating iterators from sequence-type selectors.
 - `igraph_grg_game()` now validates its arguments.
 - `igraph_layout_drl()` and its 3D version now validate their inputs.
 - `igraph_layout_kamada_kawai()`, `igraph_layout_fruchterman_reingold()`, `igraph_layout_drl()`, as well as their 3D versions now check for non-positive weights.
 - `igraph_asymmetric_preference_game()` interpreted its `type_dist_matrix` argument incorrectly.

### Other

 - Documentation improvements.

## [0.9.9] - 2022-06-04

### Changed

 - `igraph_community_walktrap()` now uses double precision floating point operations internally instead of single precision.
 - In `igraph_community_leiden()`, the `nb_clusters` output parameter is now optional (i.e. it can be `NULL`).
 - `igraph_read_graph_graphml()` no longer attempts to temporarily set the C locale, and will therefore not work correctly if the current locale uses a decimal comma.

### Fixed

 - `igraph_community_walktrap()` would return an invalid `modularity` vector when the `merges` matrix was not requested.
 - `igraph_community_walktrap()` would return a `modularity` vector that was too long for disconnected graphs. This would cause a failure in some weighted graphs when the `membership` vector was requested.
 - `igraph_community_walktrap()` now checks the weight vector: only non-negative weights are accepted, and all vertices must have non-zero strength.
 - `igraph_community_walktrap()` now returns a modularity score of NaN for graphs with no edges.
 - `igraph_community_fast_greedy()` now returns a modularity score of NaN for graphs with no edges.
 - `igraph_community_edge_betweenness()` now returns a modularity vector with a single NaN entry for graph with no edges. Previously it returned a zero-length vector.
 - `igraph_community_leading_eigenvector()` does not ignore non-ARPACK-related errors from `igraph_arpack_rssolve()` any more.
 - `igraph_preference_game()` now works correctly when `fixed_size` is true and
   `type_dist` is not given; earlier versions had a bug where more than half of
   the vertices mistakenly ended up in group 0.
 - Fixed a memory leak in `igraph_hrg_fit()` when using `start=1`.
 - `igraph_write_graph_dot()` now outputs NaN values unchanged.
 - `igraph_write_graph_dot()` no longer produces invalid DOT files when empty string attributes are present.
 - `igraph_layout_fruchterman_reingold()` and `igraph_layout_kamada_kawai()`, as well as their 3D versions, did not respect vertex coordinate bounds (`xmin`, `xmax`, etc.) when minimum values were large or maximum values were small. This is now fixed.
 - The initial coordinates of the Kamada-Kawai layout (`igraph_layout_kamada_kawai()` and `igraph_layout_kamada_kawai_3d()`) are chosen to be more in line with the original publication, improving the stability of the result. See isse #963. This changes the output of the function for the same graph, compared with previous versions. To obtain the same layout, initialize coordinates with `igraph_layout_circle()` (in 2D) or `igraph_layout_sphere()` (in 3D).
 - Improved numerical stability in Kamada-Kawai layout.
 - Corrected a problem in the calculation of displacements in `igraph_layout_fruchterman_reingold()` and its 3D version. This fixes using the "grid" variant of the algorithm on disconnected graphs.
 - `igraph_sumtree_search()` would consider search intervals open on the left and closed on the right, contrary to the documentation. This is now corrected to closed on the left and open on the right. In some cases this lead to a zero-weight element being returned for a zero search value. See issue #2080.

### Other

 - Greatly improved error reporting from foregin format parsers.
 - Documentation improvements.

## [0.9.8] - 2022-04-08

### Fixed

 - Assertion failure in `igraph_bfs()` when an empty `roots` or `restricted` vector was provided.
 - `igraph_diversity()` now returns 0 for degree-1 vertices. Previously it incorrectly returned NaN or +-Inf depending on roundoff errors.
 - `igraph_community_walktrap()` does not crash any more when provided with
   `modularity=NULL` and `membership=NULL`.

### Other

 - Documentation improvements.

## [0.9.7] - 2022-03-16

### Changed

 - `igraph_get_all_shortest_paths_dijsktra()` now uses tolerances when comparing path
   lengths, and is thus robust to numerical roundoff errors.
 - `igraph_vector_*_swap` and `igraph_matrix_swap` now take O(1) instead of O(n) and accept all sizes.

### Fixed

 - NCOL and LGL format writers no longer accept "name" and "weight" attributes
   of invalid types.
 - The LGL writer could not access numerical weight attributes, potentially leading
   to crashes.
 - External PLFIT libraries and their headers are now detected at their standard
   installation location.
 - `igraph_vector_init()` no longer accepts negative vector sizes.
 - `igraph_assortativity_nominal()` crashed on the null graph.
 - Label propagation now ensures that all labels are dominant.
 - Fixed incorrect partition results for walktrap algorithm (issue #1927)
 - Negative values returned by `igraph_rng_get_integer()` and `RNG_INTEGER()` were incorrect,
   one larger than they should have been.
 - `igraph_community_walktrap()` now checks its `steps` input argument.
 - The first modularity value reported by `igraph_community_walktrap()` was
   incorrect (it was always zero). This is now fixed.
 - `igraph_correlated_game()` would return incorrect results, or exhaust the memory,
    for most input graphs that were not generated with `igraph_erdos_renyi_game_gnp()`.

### Other

 - The C attribute handler now verifies attribute types when retrieving attributes.
 - Documentation improvements.

## [0.9.6] - 2022-01-05

 - Isomorphism class functions (`igraph_isoclass()`, `igraph_isoclass_subgraph()`,
   `igraph_isoclass_create`) and motif finder functions (`igraph_motifs_randesu()`,
   `igraph_motifs_randesu_estimate()`, `igraph_motifs_randesu_callback()`) now
   support undirected (sub)graphs of sizes 5 and 6. Previsouly only sizes 3 and 4
   were supported.

### Fixed

 - igraph would not build with MinGW when using the vendored GLPK and enabling TLS.
 - Removed some uses of `abort()` from vendored libraries, which could unexpectedly
   shut down the host language of igraph's high-level interfaces.
 - `igraph_community_label_propagation()` no longer leaves any vertices unlabeled
   when they were not reachable from any labeled ones, i.e. the returned membership
   vector is guaranteed not to contain negative values (#1853).
 - The Kamada-Kawai layout is now interruptible.
 - The Fruchterman-Reingold layout is now interruptible.
 - Fixed a bug in `igraph_cmp_epsilon()` that resulted in incorrect results for
   edge betweenness calculations in certain rare cases with x87 floating point
   math when LTO was also enabled (#1894).
 - Weighted clique related functions now fall back to the unweighted variants
   when a null vertex weight vector is given to them.
 - `igraph_erdos_renyi_game_(gnm|gnp)` would not produce self-loops for the singleton
   graph.
 - Fixed a bug in `igraph_local_efficiency()` that sometimes erroneously
   reported zero as the local efficiency of a vertex in directed graphs.
 - `igraph_vector_update()` (and its type-specific variants) did not check for
   memory allocation failure.
 - Fixed a potential crash in the GraphML reader that would be triggered by some
   invalid GraphML files.

### Other

 - `igraph_is_tree()` has improved performance and memory usage.
 - `igraph_is_connected()` has improved performance when checking weak connectedness.
 - Improved error handling in `igraph_maximal_cliques()` and related functions.
 - The build system now checks that GLPK is of a compatible version (4.57 or later).
 - The vendored `plfit` package was updated to 0.9.3.
 - You can now build igraph with an external `plfit` instead of the vendored one.
 - Documentation improvements.

## [0.9.5] - 2021-11-11

### Fixed

 - `igraph_reindex_membership()` does not allow negative membership indices any more.

 - `igraph_rewire_directed_edges()` now generates multigraphs when edge directions
   are ignored, to make it consistent with the directed case.

 - Fixed a bug in `igraph_gomory_hu_tree()` that returned only the equivalent flow
   tree instead of the cut tree (#1810).

 - Fixed a bug in the `IGRAPH_TO_UNDIRECTED_COLLAPSE` mode of
   `igraph_to_undirected()` that provided an incorrect merge vector to the
   attribute handler, leading to problems when edge attributes were merged
   using an attribute combination (#1814).

 - Fixed the behaviour of the `IGRAPH_ENABLE_LTO` option when it was set to
   `AUTO`; earlier versions had a bug where `AUTO` simply checked whether LTO
   is supported but then did not use LTO even if it was supported.

 - When using igraph from a CMake project, it is now checked that the project has
   the C++ language enabled. This is necessary for linking to igraph with CMake.

### Other

 - Improved the root selection method for disconnected graphs in the
   Reingold-Tilford layout (#1836). The new root selection method provides
   niceer results if the graph is not a tree, although it is still recommended
   to use the Sugiyama layout instead, unless the input graph is _almost_ a
   tree, in which case Reingold-Tilfold may still be preferred.

 - `igraph_decompose()` is now much faster for large graphs containing many
   isolates or small components (#960).

 - `igraph_largest_cliques()` and `igraph_clique_number()` were re-written to
   use `igraph_maximal_cliques_callback()` so they are much faster now (#804).

 - The vendored GLPK has been upgraded to GLPK 5.0.

 - Documentation improvements.

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

 - `igraph_trussness()` calculates the trussness of each edge in the graph (PR #1034, thanks to Alex Perrone and Fabio Zanini)
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

[Unreleased]: https://github.com/igraph/igraph/compare/0.9.9..HEAD
[0.9.9]: https://github.com/igraph/igraph/compare/0.9.8...0.9.9
[0.9.8]: https://github.com/igraph/igraph/compare/0.9.7...0.9.8
[0.9.7]: https://github.com/igraph/igraph/compare/0.9.6...0.9.7
[0.9.6]: https://github.com/igraph/igraph/compare/0.9.5...0.9.6
[0.9.5]: https://github.com/igraph/igraph/compare/0.9.4...0.9.5
[0.9.4]: https://github.com/igraph/igraph/compare/0.9.3...0.9.4
[0.9.3]: https://github.com/igraph/igraph/compare/0.9.2...0.9.3
[0.9.2]: https://github.com/igraph/igraph/compare/0.9.1...0.9.2
[0.9.1]: https://github.com/igraph/igraph/compare/0.9.0...0.9.1
[0.9.0]: https://github.com/igraph/igraph/compare/0.8.5...0.9.0
[0.8.5]: https://github.com/igraph/igraph/compare/0.8.4...0.8.5
[0.8.4]: https://github.com/igraph/igraph/compare/0.8.3...0.8.4
[0.8.3]: https://github.com/igraph/igraph/compare/0.8.2...0.8.3
[0.8.2]: https://github.com/igraph/igraph/compare/0.8.1...0.8.2
[0.8.1]: https://github.com/igraph/igraph/compare/0.8.0...0.8.1
[0.8.0]: https://github.com/igraph/igraph/releases/tag/0.8.0
