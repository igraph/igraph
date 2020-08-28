include(test_helpers)

# Add a compatibility alias to the "test" target so it can also be invoked as
# "make check" - for people who have it in their muscle memories from autotools
add_custom_target(build_tests)
add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND} -C ${CMAKE_BUILD_TYPE})
add_dependencies(check build_tests)

# version.at
add_legacy_tests(
  FOLDER simple NAMES
  igraph_version
)

# types.at
add_legacy_tests(
  FOLDER simple NAMES
  2wheap
  biguint
  d_indheap
  dqueue
  heap
  igraph_array
  igraph_complex
  igraph_hashtable
  igraph_i_cutheap
  igraph_marked_queue
  igraph_psumtree
  igraph_set
  igraph_sparsemat
  igraph_sparsemat2
  igraph_sparsemat3
  igraph_sparsemat4
  # igraph_sparsemat5
  igraph_sparsemat6
  igraph_sparsemat7
  igraph_sparsemat8
  igraph_sparsemat9
  igraph_sparsemat_is_symmetric
  igraph_sparsemat_minmax
  igraph_sparsemat_which_minmax
  igraph_strvector
  igraph_trie
  indheap
  matrix
  matrix2
  matrix3
  spmatrix
  stack
  vector
  vector2
  vector3
  vector_ptr
)

# basic.at
add_legacy_tests(
  FOLDER simple NAMES
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

# iterators.at
add_legacy_tests(
  FOLDER simple NAMES
  igraph_es_pairs
  igraph_es_path
  igraph_vs_nonadj
  igraph_vs_seq
  igraph_vs_vector
)

# structure_generators.at
add_legacy_tests(
  FOLDER simple NAMES
  igraph_adjacency
  igraph_atlas
  igraph_barabasi_game
  igraph_barabasi_game2
  igraph_correlated_game
  igraph_create
  igraph_degree_sequence_game
  igraph_erdos_renyi_game
  igraph_full
  igraph_grg_game
  igraph_growing_random_game
  igraph_k_regular_game
  igraph_lattice
  igraph_lcf
  igraph_preference_game
  igraph_ring
  igraph_small
  igraph_star
  igraph_tree
  igraph_weighted_adjacency
  watts_strogatz_game
)

add_legacy_tests(
  FOLDER tests NAMES
  tree
  igraph_from_prufer
  igraph_realize_degree_sequence
)

# structural_properties.at
add_legacy_tests(
  FOLDER simple NAMES
  assortativity
  bellman_ford
  biguint_betweenness
  dijkstra
  igraph_are_connected
  igraph_average_path_length
  igraph_betweenness
  igraph_cocitation
  igraph_convergence_degree
  igraph_density
  igraph_diameter
  igraph_eccentricity
  igraph_edge_betweenness
  igraph_feedback_arc_set
  igraph_feedback_arc_set_ip
  igraph_get_all_shortest_paths_dijkstra
  igraph_get_all_simple_paths
  igraph_get_shortest_paths
  igraph_get_shortest_paths2
  igraph_get_shortest_paths_dijkstra
  igraph_girth
  igraph_has_multiple
  igraph_is_loop
  igraph_is_multiple
  igraph_is_tree
  igraph_knn
  igraph_local_transitivity
  igraph_minimum_spanning_tree
  # igraph_pagerank
  igraph_radius
  igraph_reciprocity
  igraph_rewire
  igraph_similarity
  igraph_simplify
  igraph_transitive_closure_dag
  igraph_transitivity
  single_target_shortest_path
)

add_legacy_tests(
  FOLDER tests NAMES
  igraph_closeness
  igraph_topological_sorting
)

# components.at
add_legacy_tests(
  FOLDER simple NAMES
  igraph_biconnected_components
  igraph_bridges
  igraph_decompose
)

add_legacy_tests(
  FOLDER tests NAMES
  igraph_decompose_strong
)

# layout.at
add_legacy_tests(
  FOLDER simple NAMES
  igraph_i_layout_sphere
  igraph_layout_davidson_harel
  igraph_layout_grid
  igraph_layout_lgl
  igraph_layout_mds
  igraph_layout_merge
  igraph_layout_merge2
  igraph_layout_merge3
  igraph_layout_reingold_tilford
  igraph_layout_sugiyama
)

add_legacy_tests(
  FOLDER tests NAMES
  igraph_layout_reingold_tilford_extended
)

# visitors.at
add_legacy_tests(
  FOLDER simple NAMES
  igraph_bfs
  igraph_bfs2
  igraph_random_walk
)

# topology.at
add_legacy_tests(
  FOLDER simple NAMES
  VF2-compat
  igraph_isomorphic_bliss
  igraph_isomorphic_vf2
  igraph_subisomorphic_lad
  isomorphism_test
  topology
)

add_legacy_tests(
  FOLDER tests NAMES
  simplify_and_colorize
)

# coloring.at
add_legacy_tests(
  FOLDER simple NAMES
  igraph_coloring
)

# motifs.at
add_legacy_tests(
  FOLDER simple NAMES
  igraph_motifs_randesu
  triad_census
)

# foreign.at
add_legacy_tests(
  FOLDER simple NAMES
  dot
  foreign
  gml
  graphml
  igraph_read_graph_dl
  igraph_read_graph_graphdb
  igraph_read_graph_lgl
  igraph_write_graph_leda
  igraph_write_graph_lgl
  igraph_write_graph_pajek
  lineendings
  pajek
  pajek2
  pajek_bipartite
  pajek_bipartite2
  pajek_signed
)

# other.at
add_legacy_tests(
  FOLDER simple NAMES
  igraph_convex_hull
  igraph_power_law_fit
)

# operators.at
add_legacy_tests(
  FOLDER simple NAMES
  igraph_complementer
  igraph_compose
  igraph_difference
  igraph_disjoint_union
  igraph_intersection
  igraph_intersection2
  igraph_union
)

# conversion.at
add_legacy_tests(
  FOLDER simple NAMES
  adjlist
  igraph_laplacian
  igraph_to_prufer
  igraph_to_undirected
)

# flow.at
add_legacy_tests(
  FOLDER simple NAMES
  dominator_tree
  even_tarjan
  flow
  flow2
  igraph_all_st_cuts
  igraph_all_st_mincuts
  igraph_gomory_hu_tree
  igraph_mincut
)

# community.at
add_legacy_tests(
  FOLDER simple NAMES
  bug-1149658
  igraph_community_edge_betweenness
  igraph_community_fastgreedy
  igraph_community_infomap
  igraph_community_leading_eigenvector
  igraph_community_leading_eigenvector2
  igraph_community_multilevel
  igraph_community_optimal_modularity
  levc-stress
  spinglass
  walktrap
)

add_legacy_tests(
  FOLDER tests NAMES
  igraph_community_fluid_communities
  igraph_community_label_propagation
  igraph_community_leiden
  igraph_modularity
)

# cliques.at
add_legacy_tests(
  FOLDER simple NAMES
  igraph_cliques
  igraph_independent_sets
  igraph_maximal_cliques
  igraph_maximal_cliques2
  igraph_maximal_cliques3
  igraph_maximal_cliques4
  igraph_weighted_cliques
)

add_legacy_tests(
  FOLDER tests NAMES
  maximal_cliques_callback
  maximal_cliques_hist
)

# eigen.at
add_legacy_tests(
  FOLDER simple NAMES
  igraph_eigen_matrix
  igraph_eigen_matrix2
  igraph_eigen_matrix3
  igraph_eigen_matrix4
  igraph_eigen_matrix_symmetric
  igraph_eigen_matrix_symmetric_arpack
)

# attributes.at
add_legacy_tests(
  FOLDER simple NAMES
  cattr_bool_bug
  cattributes
  cattributes2
  cattributes3
  cattributes4
  cattributes5
)

add_legacy_tests(
  FOLDER tests NAMES
  cattr_bool_bug2
)

# arpack.at
add_legacy_tests(
  FOLDER simple NAMES
  blas
  eigenvector_centrality
  igraph_arpack_rnsolve
  igraph_lapack_dgeev
  igraph_lapack_dgeevx
  igraph_lapack_dgehrd
  igraph_lapack_dgesv
  igraph_lapack_dsyevr
)

# bipartite.at
add_legacy_tests(
  FOLDER simple NAMES
  igraph_bipartite_create
  igraph_bipartite_projection
)

# centralization.at
add_legacy_tests(
  FOLDER simple NAMES
  centralization
)

# separators.at
add_legacy_tests(
  FOLDER simple NAMES
  bug-1033045
  cohesive_blocks
  igraph_is_minimal_separator
  igraph_is_separator
  igraph_minimal_separators
  igraph_minimum_size_separators
)

# hrg.at
add_legacy_tests(
  FOLDER simple NAMES
  igraph_hrg
  igraph_hrg2
  igraph_hrg3
)

# microscopic.at
add_legacy_tests(
  FOLDER simple NAMES
  igraph_deterministic_optimal_imitation
  igraph_moran_process
  igraph_roulette_wheel_imitation
  igraph_stochastic_imitation
)

# mt.at
if (NOT(MSVC))
  add_legacy_tests(
    FOLDER simple NAMES
    tls1
    tls2
  )
endif()

# scg.at
add_legacy_tests(
  FOLDER simple NAMES
  igraph_scg_grouping
  igraph_scg_grouping2
  igraph_scg_grouping3
  igraph_scg_grouping4
  igraph_scg_semiprojectors
  igraph_scg_semiprojectors2
  igraph_scg_semiprojectors3
  scg
  scg2
  scg3
)

# random.at
add_legacy_tests(
  FOLDER simple NAMES
  igraph_fisher_yates_shuffle
  igraph_random_sample
  igraph_rng_get_exp
  mt
  random_seed
)

add_legacy_tests(
  FOLDER tests NAMES
  rng_reproducibility
)

# qsort.at
add_legacy_tests(
  FOLDER simple NAMES
  igraph_qsort
  igraph_qsort_r
)

# matching.at
add_legacy_tests(
  FOLDER simple NAMES
  igraph_maximum_bipartite_matching
)

# embedding.at
add_legacy_tests(
  FOLDER simple NAMES
  igraph_adjacency_spectral_embedding
)
