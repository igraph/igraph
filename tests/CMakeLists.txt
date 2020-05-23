include(test_helpers)

# version.at
add_simple_tests(
  igraph_version
)

# types.at
add_simple_tests(
  vector
  vector2
  vector3
  matrix
  matrix2
  matrix3
  dqueue
  vector_ptr
  stack
  heap
  indheap
  d_indheap
  igraph_strvector
  igraph_trie
  igraph_psumtree
  igraph_array
  igraph_hashtable
  igraph_i_cutheap
  igraph_set
  2wheap
  igraph_sparsemat
  igraph_sparsemat2
  igraph_sparsemat3
  igraph_sparsemat4
  igraph_sparsemat5
  igraph_sparsemat6
  igraph_sparsemat7
  igraph_sparsemat8
  igraph_sparsemat9
  igraph_sparsemat_is_symmetric
  igraph_sparsemat_minmax
  igraph_sparsemat_which_minmax
  spmatrix
  biguint
  igraph_marked_queue
  igraph_complex
)

# basic.at
add_simple_tests(
  igraph_empty
  igraph_copy
  igraph_add_edges
  igraph_add_vertices
  igraph_delete_edges
  igraph_delete_vertices
  igraph_neighbors
  igraph_is_directed
  igraph_degree
  igraph_get_eid
  igraph_get_eids
)

# iterators.at
add_simple_tests(
  igraph_vs_vector
  igraph_vs_nonadj
  igraph_vs_seq
  # igraph_es_adj
  # igraph_es_fromto
  igraph_es_pairs
  igraph_es_path
)

# structure_generators.at
add_simple_tests(
  igraph_create
  igraph_barabasi_game
  igraph_barabasi_game2
  igraph_erdos_renyi_game
  igraph_degree_sequence_game
  igraph_k_regular_game
  igraph_growing_random_game
  igraph_preference_game
  igraph_adjacency
  igraph_weighted_adjacency
  igraph_star
  igraph_lattice
  igraph_ring
  igraph_tree
  igraph_from_prufer
  igraph_full
  igraph_atlas
  igraph_small
  igraph_grg_game
  igraph_lcf
  watts_strogatz_game
  igraph_correlated_game
  igraph_realize_degree_sequence
)

# structural_properties.at
add_simple_tests(
  igraph_are_connected
  igraph_density
  igraph_diameter
  igraph_average_path_length
  igraph_pagerank
  igraph_rewire
  igraph_get_shortest_paths
  igraph_get_shortest_paths2
  dijkstra
  bellman_ford
  igraph_get_shortest_paths_dijkstra
  igraph_get_all_shortest_paths_dijkstra
  igraph_get_all_simple_paths
  single_target_shortest_path
  igraph_betweenness
  biguint_betweenness
  igraph_edge_betweenness
  igraph_transitivity
  igraph_local_transitivity
  igraph_reciprocity
  igraph_minimum_spanning_tree
  igraph_cocitation
  igraph_similarity
  igraph_simplify
  igraph_topological_sorting
  igraph_feedback_arc_set
  igraph_feedback_arc_set_ip
  igraph_is_loop
  igraph_is_multiple
  igraph_has_multiple
  igraph_is_tree
  igraph_girth
  igraph_convergence_degree
  # assortativity
  igraph_knn
  igraph_transitive_closure_dag
  igraph_eccentricity
  igraph_radius
)

# components.at
add_simple_tests(
  igraph_decompose
  igraph_biconnected_components
  igraph_bridges
)

# layout.at
add_simple_tests(
  igraph_layout_grid
  igraph_layout_lgl
  igraph_layout_reingold_tilford
  igraph_layout_sugiyama
  igraph_layout_mds
  igraph_i_layout_sphere
  igraph_layout_merge
  igraph_layout_merge2
  igraph_layout_merge3
  igraph_layout_davidson_harel
)

# visitors.at
add_simple_tests(
  igraph_bfs
  igraph_bfs2
  igraph_random_walk
)

# topology.at
add_simple_tests(
  topology
  igraph_isomorphic_vf2
  igraph_isomorphic_bliss
  VF2-compat
  igraph_subisomorphic_lad
  isomorphism_test
  igraph_is_degree_sequence
)

# coloring.at
add_simple_tests(
  igraph_coloring
)

# motifs.at
add_simple_tests(
  igraph_motifs_randesu
  triad_census
)

# foreign.at
add_simple_tests(
  foreign
  graphml
  igraph_write_graph_pajek
  pajek
  pajek2
  pajek_bipartite2
  pajek_signed
  pajek_bipartite
  igraph_read_graph_lgl
  igraph_write_graph_lgl
  igraph_read_graph_graphdb
  gml
  dot
  lineendings
  igraph_read_graph_dl
  igraph_write_graph_leda
)

# other.at
add_simple_tests(
  igraph_convex_hull
  igraph_power_law_fit
)

# operators.at
add_simple_tests(
  igraph_disjoint_union
  igraph_union
  igraph_intersection
  igraph_intersection2
  igraph_difference
  igraph_complementer
  igraph_compose
)

# conversion.at
add_simple_tests(
  igraph_to_undirected
  adjlist
  igraph_laplacian
  igraph_to_prufer
)

# flow.at
add_simple_tests(
  flow
  flow2
  igraph_mincut
  even_tarjan
  dominator_tree
  igraph_all_st_cuts
  igraph_all_st_mincuts
  igraph_gomory_hu_tree
)

# community.at
add_simple_tests(
  spinglass
  walktrap
  igraph_community_edge_betweenness
  igraph_community_fastgreedy
  igraph_community_leading_eigenvector
  igraph_community_leading_eigenvector2
  # levc-stress
  igraph_community_multilevel
  bug-1149658
  igraph_community_optimal_modularity
  igraph_community_infomap
)

# cliques.at
add_simple_tests(
  igraph_cliques
  igraph_maximal_cliques
  igraph_maximal_cliques2
  igraph_maximal_cliques3
  igraph_maximal_cliques4
  igraph_weighted_cliques
  igraph_independent_sets
)

# eigen.at
add_simple_tests(
  igraph_eigen_matrix_symmetric
  igraph_eigen_matrix_symmetric_arpack
  igraph_eigen_matrix
  igraph_eigen_matrix2
  igraph_eigen_matrix4
  igraph_eigen_matrix3
)

# attributes.at
add_simple_tests(
  cattributes
  cattributes2
  cattributes3
  cattributes4
  cattributes5
  cattr_bool_bug
)

# arpack.at
add_simple_tests(
  blas
  igraph_lapack_dsyevr
  igraph_lapack_dgeev
  igraph_lapack_dgeevx
  igraph_lapack_dgesv
  igraph_lapack_dgehrd
  eigenvector_centrality
  igraph_arpack_rnsolve
)

# bipartite.at
add_simple_tests(
  igraph_bipartite_create
  igraph_bipartite_projection
)

# centralization.at
add_simple_tests(
  centralization
)

# separators.at
add_simple_tests(
  igraph_is_separator
  igraph_is_minimal_separator
  igraph_minimal_separators
  bug-1033045
  igraph_minimum_size_separators
  cohesive_blocks
)

# hrg.at
add_simple_tests(
  igraph_hrg
  igraph_hrg2
  igraph_hrg3
)

# microscopic.at
add_simple_tests(
  igraph_deterministic_optimal_imitation
  igraph_stochastic_imitation
  igraph_roulette_wheel_imitation
  igraph_moran_process
)

# mt.at
# add_simple_tests(
  # tls1
  # tls2
# )

# scg.at
add_simple_tests(
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
add_simple_tests(
  random_seed
  mt
  igraph_rng_get_exp
  igraph_random_sample
  igraph_fisher_yates_shuffle
)

# qsort.at
add_simple_tests(
  igraph_qsort
  igraph_qsort_r
)

# matching.at
add_simple_tests(
  igraph_maximum_bipartite_matching
)

# embedding.at
add_simple_tests(
  igraph_adjacency_spectral_embedding
)
