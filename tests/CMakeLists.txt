include(test_helpers)
include(benchmark_helpers)

# Add a compatibility alias to the "test" target so it can also be invoked as
# "make check" - for people who have it in their muscle memories from autotools
add_custom_target(build_tests)
add_custom_target(
  check
  COMMAND ${CMAKE_CTEST_COMMAND} --progress --output-on-failure -C $<CONFIG>
  COMMENT "Executing unit tests..."
  USES_TERMINAL
)
add_dependencies(check build_tests)

# Add a custom target for benchmarks and another one for building them
# TODO(ntamas)
add_custom_target(build_benchmarks)
add_custom_target(
  benchmark
  COMMAND true
  COMMENT "Running benchmarks..."
)
add_dependencies(benchmark build_benchmarks)

# Some newer gcc version have --enable-new-dtags on by default. This then leads
# to using RUNPATH instead of RPATH. Since RUNPATH is only considered after
# LD_LIBRARY_PATH, if another version of igraph is installed somewhere it will
# be linked to that library.
include(CheckCCompilerFlag)
set(ORIG_CMAKE_REQUIRED_FLAGS ${CMAKE_REQUIRED_FLAGS})
set(CMAKE_REQUIRED_FLAGS "-Wl,--enable-new-dtags")
check_c_compiler_flag("" HAVE_ENABLE_NEW_DTAGS)
set(CMAKE_REQUIRED_FLAGS ${ORIG_CMAKE_REQUIRED_FLAGS})
if (HAVE_ENABLE_NEW_DTAGS AND BUILD_SHARED_LIBS)
  message(STATUS "Disabling new dtags for testing to use RPATH to ensure the correct library is found.")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,--disable-new-dtags")
endif()

# tutorial examples and other snippets from the documentation
add_examples(
  FOLDER examples/simple NAMES
  igraph_free
)
add_examples(
  FOLDER examples/tutorial NAMES
  tutorial1
  tutorial2
  tutorial3
)

# version.at
add_examples(
  FOLDER examples/simple NAMES
  igraph_version
)

# types.at
add_examples(
  FOLDER examples/simple NAMES
  dqueue
  igraph_sparsemat
  igraph_sparsemat3
  igraph_sparsemat4
  igraph_sparsemat6
  igraph_sparsemat7
  igraph_sparsemat8
  igraph_strvector
  igraph_vector_ptr_sort
)

add_legacy_tests(
  FOLDER tests/unit NAMES
  heap
  igraph_array
  igraph_complex
  igraph_psumtree
  igraph_sparsemat5
  igraph_sparsemat9
  igraph_sparsemat_droptol
  igraph_sparsemat_fkeep
  igraph_sparsemat_getelements_sorted
  igraph_sparsemat_is_symmetric
  igraph_sparsemat_iterator_idx
  igraph_sparsemat_minmax
  igraph_sparsemat_nonzero_storage
  igraph_sparsemat_view
  igraph_sparsemat_which_minmax
  igraph_spmatrix_add_col_values
  igraph_vector_floor
  igraph_vector_lex_cmp
  matrix
  matrix2
  matrix3
  spmatrix
  spmatrix_clear
  stack
  strvector_set2_remove_print
  vector
  vector2
  vector3
  vector_ptr
)

if ((NOT BUILD_SHARED_LIBS) OR (NOT BLAS_IS_VENDORED AND NOT ARPACK_IS_VENDORED))
  add_legacy_tests(
    FOLDER tests/unit NAMES
    igraph_sparsemat2 # Uses ARPACK and BLAS functions which are not publicly available when building with internal ARPACK/BLAS
  )
endif()

add_legacy_tests(
  FOLDER tests/unit NAMES
  2wheap
  cutheap
  d_indheap
  hashtable
  marked_queue
  set
  trie
)

# basic.at
add_examples(
  FOLDER examples/simple NAMES
  igraph_add_edges
  igraph_add_vertices
  igraph_copy
  igraph_degree
  igraph_delete_edges
  igraph_delete_vertices
  igraph_empty
  igraph_get_eid
  igraph_get_eids
  igraph_is_directed
  igraph_neighbors
)

add_legacy_tests(
  FOLDER tests/unit NAMES
  igraph_is_same_graph
  igraph_i_incident
  igraph_i_neighbors
)

# iterators.at
add_examples(
  FOLDER examples/simple NAMES
  igraph_es_pairs
#  igraph_es_fromto
  igraph_vs_nonadj
  igraph_vs_seq
  igraph_vs_vector
)

add_legacy_tests(
  FOLDER tests/unit NAMES
  edge_selectors
  igraph_es_path
  vertex_selectors
)

# structure_generators.at
add_examples(
  FOLDER examples/simple NAMES
  igraph_adjacency
  igraph_atlas
  igraph_barabasi_game
  igraph_barabasi_game2
  igraph_create
  igraph_degree_sequence_game
  igraph_erdos_renyi_game
  igraph_full
  igraph_grg_game
  igraph_lcf
  igraph_ring
  igraph_small
  igraph_star
  igraph_tree
  igraph_weighted_adjacency
)

add_legacy_tests(
  FOLDER tests/unit NAMES
  erdos_renyi_game
  full
  igraph_barabasi_aging_game
  igraph_bipartite_game
  igraph_callaway_traits_game
  igraph_cited_type_game
  igraph_citing_cited_type_game
  igraph_correlated_game
  igraph_degree_sequence_game
  igraph_establishment_game
  igraph_extended_chordal_ring
  igraph_forest_fire_game
  igraph_from_prufer
  igraph_full_citation
  igraph_grg_game
  igraph_growing_random_game
  igraph_k_regular_game
  igraph_lastcit_game
  igraph_lattice
  igraph_linegraph
  igraph_kautz
  igraph_preference_game
  igraph_realize_degree_sequence
  igraph_recent_degree_aging_game
  igraph_recent_degree_game
  igraph_sbm_game
  igraph_simple_interconnected_islands_game
  igraph_static_power_law_game
  tree
  tree_game
  ring
  watts_strogatz_game
)

# structural_properties.at
add_examples(
  FOLDER examples/simple NAMES
  assortativity
  bellman_ford
  dijkstra
  igraph_average_path_length
  igraph_cocitation
  igraph_diameter
  igraph_eccentricity
  igraph_feedback_arc_set
  igraph_feedback_arc_set_ip
  igraph_get_all_shortest_paths_dijkstra
  igraph_get_shortest_paths
  igraph_get_shortest_paths_dijkstra
  igraph_girth
  igraph_has_multiple
  igraph_is_loop
  igraph_is_multiple
  igraph_knn
  igraph_minimum_spanning_tree
  igraph_pagerank
  igraph_radius
  igraph_reciprocity
  igraph_similarity
  igraph_simplify
  igraph_topological_sorting
  igraph_transitivity
)

add_legacy_tests(
  FOLDER tests/unit NAMES
  igraph_rewire # Uses internal igraph_i_rewire
)

add_legacy_tests(
  FOLDER tests/unit NAMES
  all_shortest_paths
  efficiency
  global_transitivity
  igraph_adjacent_triangles
  igraph_are_connected
  igraph_average_path_length
  igraph_average_path_length_dijkstra
  igraph_betweenness
  igraph_closeness
  igraph_convergence_degree
  igraph_count_multiple
  igraph_density
  igraph_diversity
  igraph_eccentricity
  igraph_edge_betweenness
  igraph_get_all_simple_paths
  igraph_get_shortest_paths2
  igraph_get_shortest_paths_bellman_ford
  igraph_is_bipartite
  igraph_is_chordal
  igraph_is_mutual
  igraph_is_tree
  igraph_list_triangles
  igraph_local_scan_k_ecount
  igraph_local_transitivity
  igraph_neighborhood
  igraph_neighborhood_graphs
  igraph_neighborhood_size
  igraph_pagerank
  igraph_shortest_paths_johnson
  igraph_transitive_closure_dag
  igraph_transitivity_avglocal_undirected
  igraph_transitivity_barrat
  harmonic_centrality
  random_spanning_tree
  single_target_shortest_path
  topological_sorting
)

add_legacy_tests(
  FOLDER tests/regression NAMES
  bug_1760
)

# components.at
add_examples(
  FOLDER examples/simple NAMES
  igraph_biconnected_components
  igraph_decompose
)

add_legacy_tests(
  FOLDER tests/unit NAMES
  igraph_bridges
  igraph_decompose_strong
  igraph_subcomponent
)

# layout.at
add_examples(
  FOLDER examples/simple NAMES
  igraph_layout_reingold_tilford
)

add_legacy_tests(
  FOLDER tests/unit NAMES
  igraph_layout_drl
  igraph_layout_bipartite
  igraph_layout_fruchterman_reingold
  igraph_layout_graphopt
  igraph_layout_grid
  igraph_layout_lgl
  igraph_layout_mds
  igraph_layout_merge2
  igraph_layout_merge3
  igraph_layout_random_3d
  igraph_layout_reingold_tilford_circular
  igraph_layout_reingold_tilford_extended
  igraph_layout_star
  igraph_layout_sugiyama
)
add_legacy_tests(
  FOLDER tests/regression NAMES
  igraph_layout_kamada_kawai_3d_bug_1462
  igraph_layout_reingold_tilford_bug_879
)

add_legacy_tests(
  FOLDER tests/unit NAMES
  igraph_i_layout_sphere
  igraph_layout_davidson_harel # Uses igraph_i_layout_segments_intersect and igraph_i_layout_point_segment_dist2
  igraph_layout_merge # Uses igraph_i_layout_merge functions
)

# visitors.at
add_examples(
  FOLDER examples/simple NAMES
  igraph_bfs
  igraph_bfs_callback
  igraph_bfs_simple
)
add_legacy_tests(
  FOLDER tests/unit NAMES
  bfs
  bfs_simple
  igraph_random_walk
)

# topology.at
add_examples(
  FOLDER examples/simple NAMES
  igraph_isomorphic_vf2
  igraph_subisomorphic_lad
)

add_legacy_tests(
  FOLDER tests/unit NAMES
  simplify_and_colorize
  bliss_automorphisms
  igraph_get_isomorphisms_vf2
  igraph_get_subisomorphisms_vf2
  igraph_subisomorphic
  igraph_isomorphic_bliss
  isomorphism_test
  isoclasses
  VF2-compat
)

# coloring.at
add_examples(
  FOLDER examples/simple NAMES
  igraph_coloring
)

# motifs.at
add_examples(
  FOLDER examples/simple NAMES
  igraph_motifs_randesu
)

add_legacy_tests(
  FOLDER tests/unit NAMES
  igraph_dyad_census
  igraph_motifs_randesu_estimate
  igraph_motifs_randesu_no
  triad_census
)

# foreign.at
add_examples(
  FOLDER examples/simple NAMES
  dot
  foreign
  gml
  graphml
  igraph_read_graph_dl
  igraph_read_graph_graphdb
  igraph_read_graph_lgl
  igraph_write_graph_lgl
  igraph_write_graph_pajek
)

add_legacy_tests(
  FOLDER tests/unit NAMES
  igraph_write_graph_leda
  igraph_write_graph_dimacs
  lineendings
  pajek
  pajek2
  pajek_bipartite
  pajek_bipartite2
  pajek_signed
)

# other.at
add_examples(
  FOLDER examples/simple NAMES
  igraph_convex_hull
  igraph_power_law_fit
)

# operators.at
add_examples(
  FOLDER examples/simple NAMES
  igraph_complementer
  igraph_compose
  igraph_difference
  igraph_disjoint_union
  igraph_intersection
  igraph_union
)

add_legacy_tests(
  FOLDER tests/unit NAMES
  igraph_induced_subgraph
  igraph_induced_subgraph_map
  igraph_intersection2
  igraph_rewire_directed_edges
)

# conversion.at
add_examples(
  FOLDER examples/simple NAMES
  adjlist
  igraph_laplacian
  igraph_to_undirected
)

add_legacy_tests(
  FOLDER tests/unit NAMES
  adjlist
  igraph_get_adjacency_sparse
  igraph_to_directed
  igraph_to_prufer
  inclist
)

# flow.at
add_examples(
  FOLDER examples/simple NAMES
  dominator_tree
  even_tarjan
  flow
  flow2
  igraph_all_st_mincuts
  igraph_mincut
)

add_legacy_tests(
  FOLDER tests/unit NAMES
  igraph_st_mincut_value
  igraph_vertex_disjoint_paths
  igraph_adhesion
  igraph_cohesion
  igraph_residual_graph
  igraph_edge_disjoint_paths
  igraph_st_edge_connectivity
  igraph_st_mincut
)

add_legacy_tests(
  FOLDER tests/unit NAMES
  igraph_all_st_cuts # Uses igraph_marked_queue, which is internal.
)

add_legacy_tests(
  FOLDER tests/unit NAMES
  igraph_gomory_hu_tree
)

# community.at
add_examples(
  FOLDER examples/simple NAMES
  igraph_community_edge_betweenness
  igraph_community_fastgreedy
  igraph_community_fluid_communities
  igraph_community_label_propagation
  igraph_community_leading_eigenvector
  igraph_community_leiden
  igraph_community_multilevel
  igraph_community_optimal_modularity
  walktrap
)

add_legacy_tests(
  FOLDER tests/unit NAMES
  community_leiden
  community_label_propagation
  igraph_community_infomap
  igraph_community_leading_eigenvector2
  igraph_compare_communities
  igraph_le_community_to_membership
  igraph_modularity
  igraph_modularity_matrix
  igraph_split_join_distance
  levc-stress
  spinglass
)
add_legacy_tests(
  FOLDER tests/regression NAMES
  bug-1149658
)

# use a higher test timeout for the Infomap algorithm
set_tests_properties("test::igraph_community_infomap" PROPERTIES TIMEOUT 150)

# cliques.at
add_examples(
  FOLDER examples/simple NAMES
  igraph_cliques
  igraph_independent_sets
  igraph_maximal_cliques
)

add_legacy_tests(
  FOLDER tests/unit NAMES
  igraph_clique_size_hist
  igraph_maximal_cliques2
  igraph_maximal_cliques3
  igraph_maximal_cliques4
  igraph_maximal_cliques_file
  igraph_weighted_cliques
  maximal_cliques_callback
  maximal_cliques_hist
)

# eigen.at
add_legacy_tests(
  FOLDER tests/unit NAMES
  igraph_eigen_matrix
  igraph_eigen_matrix2
  igraph_eigen_matrix3
  igraph_eigen_matrix4
  igraph_eigen_matrix_symmetric
  igraph_eigen_matrix_symmetric_arpack
)

# attributes.at
add_examples(
  FOLDER examples/simple NAMES
  cattributes
  cattributes2
  cattributes3
  cattributes4
)

add_legacy_tests(
  FOLDER tests/unit NAMES
  igraph_attribute_combination_remove
  cattributes5
)
add_legacy_tests(
  FOLDER tests/regression NAMES
  cattr_bool_bug
  cattr_bool_bug2
)

# arpack.at
add_examples(
  FOLDER examples/simple NAMES
  blas
  eigenvector_centrality
  igraph_lapack_dgeev
  igraph_lapack_dgeevx
  igraph_lapack_dgesv
  igraph_lapack_dsyevr
)

add_legacy_tests(
  FOLDER tests/unit NAMES
  dgemv
  igraph_arpack_rnsolve
  igraph_arpack_unpack_complex
  igraph_lapack_dgehrd
  igraph_lapack_dgetrf
  igraph_lapack_dgetrs
)

# bipartite.at
add_examples(
  FOLDER examples/simple NAMES
  igraph_bipartite_create
  igraph_bipartite_projection
)

add_legacy_tests(
  FOLDER tests/unit NAMES
  igraph_get_incidence
)

# centralization.at
add_examples(
  FOLDER examples/simple NAMES
  centralization
)

# eulerian.at
add_legacy_tests(
  FOLDER tests/unit NAMES
  igraph_is_eulerian
  igraph_eulerian_cycle
  igraph_eulerian_path
)

# separators.at
add_examples(
  FOLDER examples/simple NAMES
  cohesive_blocks
  igraph_is_minimal_separator
  igraph_is_separator
  igraph_minimal_separators
  igraph_minimum_size_separators
)
add_legacy_tests(
  FOLDER tests/regression NAMES
  bug-1033045
)

# hrg.at
add_legacy_tests(
  FOLDER tests/unit NAMES
  igraph_hrg
  igraph_hrg2
  igraph_hrg3
)

# microscopic.at
add_examples(
  FOLDER examples/simple NAMES
  igraph_deterministic_optimal_imitation
  igraph_roulette_wheel_imitation
  igraph_stochastic_imitation
)

add_legacy_tests(
  FOLDER tests/unit NAMES
  igraph_moran_process
)

# mt.at -- only if we have pthreads
if(CMAKE_USE_PTHREADS_INIT)
  add_legacy_tests(
    FOLDER tests/unit NAMES tls1
    LIBRARIES Threads::Threads
  )

  # tls2 should be added only if we use vendored ARPACK because a non-vendored
  # ARPACK is not guaranteed to be thread-safe
  if(ARPACK_IS_VENDORED AND NOT BUILD_SHARED_LIBS)
    add_legacy_tests(
      FOLDER tests/unit NAMES tls2
      LIBRARIES Threads::Threads
    )
  endif()
endif()

# scg.at
add_examples(
  FOLDER examples/simple NAMES
  igraph_scg_grouping
  igraph_scg_grouping2
  igraph_scg_grouping3
  igraph_scg_grouping4
  igraph_scg_semiprojectors
  igraph_scg_semiprojectors2
  igraph_scg_semiprojectors3
  scg
)

add_legacy_tests(
  FOLDER tests/unit NAMES
  scg2
  scg3
)

# random.at
add_examples(
  FOLDER examples/simple NAMES
  igraph_fisher_yates_shuffle
  igraph_random_sample
  random_seed
)

add_legacy_tests(
  FOLDER tests/unit NAMES
  igraph_rng_get_exp
  igraph_rng_get_integer
  mt
  rng_reproducibility
  rng_init_destroy_max_min_name_set_default
)

# qsort.at
add_legacy_tests(
  FOLDER tests/unit NAMES
  igraph_qsort
  igraph_qsort_r
)

# matching.at
add_examples(
  FOLDER examples/simple NAMES
  igraph_maximum_bipartite_matching
)

# embedding.at
add_legacy_tests(
  FOLDER tests/unit NAMES
  igraph_adjacency_spectral_embedding
)

# graphicality

add_legacy_tests(
  FOLDER tests/unit NAMES
  igraph_is_graphical
  igraph_is_bigraphical
)

# handlers

add_legacy_tests(
  FOLDER tests/unit NAMES
  fatal_handler
  igraph_progress_handler_stderr
  igraph_set_progress_handler
)

# error output

add_legacy_tests(
  FOLDER tests/unit NAMES
  error_macros
)

# GLPK

add_legacy_tests(
  FOLDER tests/unit NAMES
  glpk_error
)

# regression and fuzzing tests

add_legacy_tests(
  FOLDER tests/regression NAMES
  igraph_read_graph_gml_invalid_inputs
)

# non-graph

add_legacy_tests(
  FOLDER tests/unit NAMES
  igraph_running_mean
  igraph_solve_lsap
)

# memory allocation

add_legacy_tests(
  FOLDER tests/unit NAMES
  zero_allocs
)

# simulation

add_legacy_tests(
  FOLDER tests/unit NAMES
  igraph_sir
)


# benchmarks
add_benchmarks(
  NAMES
  igraph_betweenness
  igraph_betweenness_weighted
  igraph_cliques
  igraph_closeness_weighted
  igraph_coloring
  igraph_maximal_cliques
  igraph_pagerank
  igraph_pagerank_weighted
  igraph_power_law_fit
  igraph_random_walk
  igraph_transitivity
)
