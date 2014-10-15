
#   IGraph R package
#   Copyright (C) 2014  Gabor Csardi <csardi.gabor@gmail.com>
#   334 Harvard street, Cambridge, MA 02139 USA
#   
#   This program is free software; you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation; either version 2 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#   
#   You should have received a copy of the GNU General Public License
#   along with this program; if not, write to the Free Software
#   Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA
#   02110-1301 USA
#
###################################################################

#' @include attributes.R auto.R basic.R bipartite.R centrality.R
#' @include cliques.R cocitation.R cohesive.blocks.R community.R
#' @include components.R console.R conversion.R decomposition.R demo.R
#' @include epi.R fit.R flow.R foreign.R games.R glet.R hrg.R indexing.R
#' @include interface.R iterators.R layout.R minimum.spanning.tree.R
#' @include motifs.R nexus.R operators.R other.R package.R par.R plot.R
#' @include plot.common.R plot.shapes.R pp.R print.R scg.R socnet.R
#' @include sparsedf.R structural.properties.R structure.generators.R
#' @include structure.info.R test.R tkplot.R topology.R layout_drl.R
NULL

## For the future, right now, we do not warn or even message

deprecated <- function(old, new) {
  assign(old, new, envir = asNamespace(packageName()))
}

#' @export add.edges
deprecated("add.edges", add_edges)
#' @export add.vertex.shape
deprecated("add.vertex.shape", add_shape)
#' @export add.vertices
deprecated("add.vertices", add_vertices)
#' @export adjacent.triangles
deprecated("adjacent.triangles", count_triangles)
#' @export articulation.points
deprecated("articulation.points", articulation_points)
#' @export aging.prefatt.game
deprecated("aging.prefatt.game", sample_pa_age)
#' @export aging.ba.game
deprecated("aging.ba.game", sample_pa_age)
#' @export aging.barabasi.game
deprecated("aging.barabasi.game", sample_pa_age)
#' @export alpha.centrality
deprecated("alpha.centrality", alpha_centrality)
#' @export are.connected
deprecated("are.connected", is_connected_to)
#' @export asPhylo
deprecated("asPhylo", as_phylo)
#' @export asPhylo.communities
deprecated("asPhylo.communities", as_phylo.communities)
#' @export asPhylo.igraphHRG
deprecated("asPhylo.igraphHRG", as_phylo.igraphHRG)
#' @export assortativity.degree
deprecated("assortativity.degree", assortativity_degree)
#' @export assortativity.nominal
deprecated("assortativity.nominal", assortativity_nominal)
#' @export asymmetric.preference.game
deprecated("asymmetric.preference.game", sample_asym_pref)
#' @export authority.score
deprecated("authority.score", authority_score)
#' @export autocurve.edges
deprecated("autocurve.edges", curve_multiple)
#' @export average.path.length
deprecated("average.path.length", mean_distance)

#' @export ba.game
deprecated("ba.game", sample_pa)
#' @export barabasi.game
deprecated("barabasi.game", sample_pa)
#' @export betweenness.estimate
deprecated("betweenness.estimate", estimate_betweenness)
#' @export biconnected.components
deprecated("biconnected.components", biconnected_components)
#' @export bipartite.mapping
deprecated("bipartite.mapping", bipartite_mapping)
#' @export bipartite.projection
deprecated("bipartite.projection", bipartite_projection)
#' @export bipartite.projection.size
deprecated("bipartite.projection.size", bipartite_projection_size)
#' @export bipartite.random.game
deprecated("bipartite.random.game", sample_bipartite)
#' @export blockGraphs
deprecated("blockGraphs", graphs_from_cohesive_blocks)
#' @export bonpow
deprecated("bonpow", power_centrality)

#' @export callaway.traits.game
deprecated("callaway.traits.game", sample_traits_callaway)
#' @export canonical.permutation
deprecated("canonical.permutation", canonical_permutation)
#' @export centralization.betweenness
deprecated("centralization.betweenness", centr_betw)
#' @export centralization.betweenness.tmax
deprecated("centralization.betweenness.tmax", centr_betw_tmax)
#' @export centralization.closeness
deprecated("centralization.closeness", centr_clo)
#' @export centralization.closeness.tmax
deprecated("centralization.closeness.tmax", centr_clo_tmax)
#' @export centralization.degree
deprecated("centralization.degree", centr_degree)
#' @export centralization.degree.tmax
deprecated("centralization.degree.tmax", centr_degree_tmax)
#' @export centralization.evcent
deprecated("centralization.evcent", centr_eigen)
#' @export centralization.evcent.tmax
deprecated("centralization.evcent.tmax", centr_eigen_tmax)
#' @export centralize.scores
deprecated("centralize.scores", centralize)
#' @export cited.type.game
deprecated("cited.type.game", sample_cit_types)
#' @export citing.cited.type.game
deprecated("citing.cited.type.game", sample_cit_cit_types)
#' @export clique.number
deprecated("clique.number", clique_num)
#' @export closeness.estimate
deprecated("closeness.estimate", estimate_closeness)
#' @export cluster.distribution
deprecated("cluster.distribution", component_distribution)
#' @export clusters
deprecated("clusters", components)
#' @export code.length
deprecated("code.length", code_len)
#' @export cohesive.blocks
deprecated("cohesive.blocks", cohesive_blocks)
#' @export compare.communities
deprecated("compare.communities", compare)
#' @export connect.neighborhood
deprecated("connect.neighborhood", connect)
#' @export contract.vertices
deprecated("contract.vertices", contract)
#' @export convex.hull
deprecated("convex.hull", convex_hull)
#' @export count.multiple
deprecated("count.multiple", count_multiple)
#' @export cutat
deprecated("cutat", cut_at)

#' @export decompose.graph
deprecated("decompose.graph", decompose)
#' @export degree.distribution
deprecated("degree.distribution", degree_distribution)
#' @export degree.sequence.game
deprecated("degree.sequence.game", sample_degseq)
#' @export delete.edges
deprecated("delete.edges", delete_edges)
#' @export delete.vertices
deprecated("delete.vertices", delete_vertices)
#' @export dendPlot
deprecated("dendPlot", plot_dendrogram)
#' @export dendPlot.communities
deprecated("dendPlot.communities", plot_dendrogram.communities)
#' @export dendPlot.igraphHRG
deprecated("dendPlot.igraphHRG", plot_dendrogram.igraphHRG)
#' @export dominator.tree
deprecated("dominator.tree", dominator_tree)
#' @export dyad.census
deprecated("dyad.census", dyad_census)

#' @export ecount
deprecated("ecount", gsize)
#' @export edge.betweenness
deprecated("edge.betweenness", edge_betweenness)
#' @export edge.betweenness.community
deprecated("edge.betweenness.community", cluster_edge_betweenness)
#' @export edge.betweenness.estimate
deprecated("edge.betweenness.estimate", estimate_edge_betweenness)
#' @export edge.connectivity
deprecated("edge.connectivity", edge_connectivity)
#' @export edge.disjoint.paths
deprecated("edge.disjoint.paths", edge_disjoint_paths)
#' @export establishment.game
deprecated("establishment.game", sample_traits)
#' @export evcent
deprecated("evcent", eigen_centrality)

#' @export farthest.nodes
deprecated("farthest.nodes", farthest_vertices)
#' @export fastgreedy.community
deprecated("fastgreedy.community", cluster_fast_greedy)
#' @export forest.fire.game
deprecated("forest.fire.game", sample_forestfire)

#' @export get.adjedgelist
deprecated("get.adjedgelist", as_adj_edge_list)
#' @export get.adjlist
deprecated("get.adjlist", as_adj_list)
#' @export get.adjacency
deprecated("get.adjacency", as_adj)
#' @export get.data.frame
deprecated("get.data.frame", as_data_frame)
#' @export get.edge.attribute
deprecated("get.edge.attribute", edge_attr)
#' @export get.edgelist
deprecated("get.edgelist", as_edgelist)
#' @export get.graph.attribute
deprecated("get.graph.attribute", graph_attr)
#' @export get.incidence
deprecated("get.incidence", as_incidence_matrix)
#' @export get.stochastic
deprecated("get.stochastic", stochastic_matrix)
#' @export get.vertex.attribute
deprecated("get.vertex.attribute", vertex_attr)
#' @export graph.adhesion
deprecated("graph.adhesion", adhesion)
#' @export graph.adjacency
deprecated("graph.adjacency", graph_from_adjacency_matrix)
#' @export graph.adjlist
deprecated("graph.adjlist", graph_from_adj_list)
#' @export graph.atlas
deprecated("graph.atlas", graph_atlas)
#' @export graph.automorphisms
deprecated("graph.automorphisms", automorphisms)
#' @export graph.bfs
deprecated("graph.bfs", bfs)
#' @export graph.bipartite
deprecated("graph.bipartite", bipartite_graph)
#' @export graph.cohesion
deprecated("graph.cohesion", cohesion)
#' @export graph.complementer
deprecated("graph.complementer", complementer)
#' @export graph.compose
deprecated("graph.compose", compose)
#' @export graph.coreness
deprecated("graph.coreness", coreness)
#' @export graph.data.frame
deprecated("graph.data.frame", graph_from_data_frame)
#' @export graph.de.bruijn
deprecated("graph.de.bruijn", de_bruijn_graph)
#' @export graph.density
deprecated("graph.density", density)
#' @export graph.disjoint.union
deprecated("graph.disjoint.union", disjoint_union)
#' @export graph.dfs
deprecated("graph.dfs", dfs)
#' @export graph.difference
deprecated("graph.difference", difference)
#' @export graph.diversity
deprecated("graph.diversity", diversity)
#' @export graph.edgelist
deprecated("graph.edgelist", graph_from_edgelist)
#' @export graph.eigen
deprecated("graph.eigen", spectrum)
#' @export graph.empty
deprecated("graph.empty", empty_graph)
#' @export graph.extended.chordal.ring
deprecated("graph.extended.chordal.ring", chordal_ring)
#' @export graph.formula
deprecated("graph.formula", graph_from_formula)
#' @export graph.full
deprecated("graph.full", full_graph)
#' @export graph.full.bipartite
deprecated("graph.full.bipartite", full_bipartite_graph)
#' @export graph.full.citation
deprecated("graph.full.citation", full_citation_graph)
#' @export graph.graphdb
deprecated("graph.graphdb", graph_from_graphdb)
#' @export graph.incidence
deprecated("graph.incidence", graph_from_incidence_matrix)
#' @export graph.isoclass
deprecated("graph.isoclass", iso_class)
#' @export graph.isocreate
deprecated("graph.isocreate", graph_from_iso_class)
#' @export graph.kautz
deprecated("graph.kautz", kautz_graph)
#' @export graph.knn
deprecated("graph.knn", knn)
#' @export graph.laplacian
deprecated("graph.laplacian", laplacian_matrix)
#' @export graph.lattice
deprecated("graph.lattice", lattice)
#' @export graph.lcf
deprecated("graph.lcf", graph_from_lcf)
#' @export graph.maxflow
deprecated("graph.maxflow", max_flow)
#' @export graph.mincut
deprecated("graph.mincut", min_cut)
#' @export graph.motifs
deprecated("graph.motifs", motifs)
#' @export graph.motifs.est
deprecated("graph.motifs.est", sample_motifs)
#' @export graph.motifs.no
deprecated("graph.motifs.no", count_motifs)
#' @export graph.neighborhood
deprecated("graph.neighborhood", ego_graph)
#' @export graph.star
deprecated("graph.star", star)
#' @export graph.strength
deprecated("graph.strength", strength)
#' @export graph.tree
deprecated("graph.tree", tree)
#' @export graph.union
deprecated("graph.union", union)
#' @export graph.ring
deprecated("graph.ring", ring)
#' @export graphlets.candidate.basis
deprecated("graphlets.candidate.basis", graphlet_basis)
#' @export graphlets.project
deprecated("graphlets.project", graphlet_proj)
#' @export growing.random.game
deprecated("growing.random.game", sample_growing)
#' @export grg.game
deprecated("grg.game", sample_grg)

#' @export has.multiple
deprecated("has.multiple", any_multiple)
#' @export hrg.consensus
deprecated("hrg.consensus", consensus_tree)
#' @export hrg.create
deprecated("hrg.create", hrg)
#' @export hrg.dendrogram
deprecated("hrg.dendrogram", hrg_tree)
#' @export hrg.game
deprecated("hrg.game", sample_hrg)
#' @export hrg.fit
deprecated("hrg.fit", fit_hrg)
#' @export hrg.predict
deprecated("hrg.predict", predict_edges)
#' @export hub.score
deprecated("hub.score", hub_score)

#' @export igraph.arpack.default
deprecated("igraph.arpack.default", arpack_defaults)
#' @export igraph.console
deprecated("igraph.console", console)
#' @export igraph.eigen.default
deprecated("igraph.eigen.default", eigen_defaults)
#' @export igraph.sample
deprecated("igraph.sample", sample_seq)
#' @export igraph.version
deprecated("igraph.version", igraph_version)
#' @export igraphdemo
deprecated("igraphdemo", igraph_demo)
#' @export igraphtest
deprecated("igraphtest", igraph_test)
#' @export independence.number
deprecated("independence.number", ivs_size)
#' @export independent.vertex.sets
deprecated("independent.vertex.sets", ivs)
#' @export infomap.community
deprecated("infomap.community", cluster_infomap)
#' @export induced.subgraph
deprecated("induced.subgraph", induced_subgraph)
#' @export interconnected.islands.game
deprecated("interconnected.islands.game", sample_islands)
#' @export is.bipartite
deprecated("is.bipartite", is_bipartite)
#' @export is.chordal
deprecated("is.chordal", is_chordal)
#' @export is.connected
deprecated("is.connected", is_connected)
#' @export is.dag
deprecated("is.dag", is_dag)
#' @export is.degree.sequence
deprecated("is.degree.sequence", is_degseq)
#' @export is.directed
deprecated("is.directed", is_directed)
#' @export is.graphical.degree.sequence
deprecated("is.graphical.degree.sequence", is_graphical)
#' @export is.hierarchical
deprecated("is.hierarchical", is_hierarchical)
#' @export is.igraph
deprecated("is.igraph", is_igraph)
#' @export is.loop
deprecated("is.loop", which_loop)
#' @export is.matching
deprecated("is.matching", is_matching)
#' @export is.maximal.matching
deprecated("is.maximal.matching", is_max_matching)
#' @export is.minimal.separator
deprecated("is.minimal.separator", is_min_separator)
#' @export is.multiple
deprecated("is.multiple", which_multiple)
#' @export is.mutual
deprecated("is.mutual", which_mutual)
#' @export is.named
deprecated("is.named", is_named)
#' @export is.separator
deprecated("is.separator", is_separator)
#' @export is.simple
deprecated("is.simple", is_simple)
#' @export is.weighted
deprecated("is.weighted", is_weighted)

#' @export k.regular.game
deprecated("k.regular.game", sample_k_regular)

#' @export label.propagation.community
deprecated("label.propagation.community", cluster_label_prop)
#' @export largest.cliques
deprecated("largest.cliques", largest_cliques)
#' @export largest.independent.vertex.sets
deprecated("largest.independent.vertex.sets", largest_ivs)
#' @export lastcit.game
deprecated("lastcit.game", sample_last_cit)
#' @export layout.auto
deprecated("layout.auto", layout_nicely)
#' @export layout.bipartite
deprecated("layout.bipartite", layout_as_bipartite)
#' @export layout.circle
deprecated("layout.circle", layout_in_circle)
#' @export layout.davidson.harel
deprecated("layout.davidson.harel", layout_with_dh)
#' @export layout.drl
deprecated("layout.drl", layout_with_drl)
#' @export layout.fruchterman.reingold
deprecated("layout.fruchterman.reingold", layout_with_fr)
#' @export layout.gem
deprecated("layout.gem", layout_with_gem)
#' @export layout.graphopt
deprecated("layout.graphopt", layout_with_graphopt)
#' @export layout.grid
deprecated("layout.grid", layout_on_grid)
#' @export layout.kamada.kawai
deprecated("layout.kamada.kawai", layout_with_kk)
#' @export layout.lgl
deprecated("layout.lgl", layout_with_lgl)
#' @export layout.mds
deprecated("layout.mds", layout_with_mds)
#' @export layout.merge
deprecated("layout.merge", merge_coords)
#' @export layout.norm
deprecated("layout.norm", norm_coords)
#' @export layout.random
deprecated("layout.random", layout_randomly)
#' @export layout.reingold.tilford
deprecated("layout.reingold.tilford", layout_as_tree)
#' @export layout.sphere
deprecated("layout.sphere", layout_on_sphere)
#' @export layout.star
deprecated("layout.star", layout_as_star)
#' @export layout.sugiyama
deprecated("layout.sugiyama", layout_with_sugiyama)
#' @export leading.eigenvector.community
deprecated("leading.eigenvector.community", cluster_leading_eigen)
#' @export line.graph
deprecated("line.graph", line_graph)
#' @export list.edge.attributes
deprecated("list.edge.attributes", edge_attr_names)
#' @export list.graph.attributes
deprecated("list.graph.attributes", graph_attr_names)
#' @export list.vertex.attributes
deprecated("list.vertex.attributes", vertex_attr_names)

#' @export maxcohesion
deprecated("maxcohesion", max_cohesion)
#' @export maximal.cliques
deprecated("maximal.cliques", max_cliques)
#' @export maximal.cliques.count
deprecated("maximal.cliques.count", count_max_cliques)
#' @export maximal.independent.vertex.sets
deprecated("maximal.independent.vertex.sets", maximal_ivs)
#' @export minimal.st.separators
deprecated("minimal.st.separators", min_st_separators)
#' @export maximum.bipartite.matching
deprecated("maximum.bipartite.matching", max_bipartite_match)
#' @export maximum.cardinality.search
deprecated("maximum.cardinality.search", max_cardinality)
#' @export minimum.size.separators
deprecated("minimum.size.separators", min_separators)
#' @export minimum.spanning.tree
deprecated("minimum.spanning.tree", mst)
#' @export mod.matrix
deprecated("mod.matrix", modularity_matrix)
#' @export multilevel.community
deprecated("multilevel.community", cluster_louvain)

#' @export neighborhood
deprecated("neighborhood", ego)
#' @export neighborhood.size
deprecated("neighborhood.size", ego_size)
#' @export nexus.get
deprecated("nexus.get", nexus_get)
#' @export nexus.info
deprecated("nexus.info", nexus_info)
#' @export nexus.list
deprecated("nexus.list", nexus_list)
#' @export nexus.search
deprecated("nexus.search", nexus_search)
#' @export no.clusters
deprecated("no.clusters", count_components)

#' @export optimal.community
deprecated("optimal.community", cluster_optimal)

#' @export page.rank
deprecated("page.rank", page_rank)
#' @export page.rank.old
deprecated("page.rank.old", page_rank_old)
#' @export path.length.hist
deprecated("path.length.hist", distance_table)
#' @export permute.vertices
deprecated("permute.vertices", permute)
#' @export piecewise.layout
deprecated("piecewise.layout", layout_components)
#' @export plotHierarchy
deprecated("plotHierarchy", plot_hierarchy)
#' @export power.law.fit
deprecated("power.law.fit", fit_power_law)
#' @export preference.game
deprecated("preference.game", sample_pref)

#' @export read.graph
deprecated("read.graph", read_graph)
#' @export remove.edge.attribute
deprecated("remove.edge.attribute", delete_edge_attr)
#' @export remove.graph.attribute
deprecated("remove.graph.attribute", delete_graph_attr)
#' @export remove.vertex.attribute
deprecated("remove.vertex.attribute", delete_vertex_attr)
#' @export running.mean
deprecated("running.mean", running_mean)

#' @export sbm.game
deprecated("sbm.game", sample_sbm)
#' @export scgGrouping
deprecated("scgGrouping", scg_group)
#' @export scgNormEps
deprecated("scgNormEps", scg_eps)
#' @export scgSemiProjectors
deprecated("scgSemiProjectors", scg_semi_proj)
#' @export set.edge.attribute
deprecated("set.edge.attribute", set_edge_attr)
#' @export set.graph.attribute
deprecated("set.graph.attribute", set_graph_attr)
#' @export set.vertex.attribute
deprecated("set.vertex.attribute", set_vertex_attr)
#' @export shortest.paths
deprecated("shortest.paths", distances)
#' @export showtrace
deprecated("showtrace", show_trace)
#' @export spinglass.community
deprecated("spinglass.community", cluster_spinglass)
#' @export stCuts
deprecated("stCuts", st_cuts)
#' @export stMincuts
deprecated("stMincuts", st_min_cuts)
#' @export static.fitness.game
deprecated("static.fitness.game", sample_fitness)
#' @export static.power.law.game
deprecated("static.power.law.game", sample_fitness_pl)
#' @export subgraph.centrality
deprecated("subgraph.centrality", subgraph_centrality)

#' @export tkplot.canvas
deprecated("tkplot.canvas", tk_canvas)
#' @export tkplot.center
deprecated("tkplot.center", tk_center)
#' @export tkplot.close
deprecated("tkplot.close", tk_close)
#' @export tkplot.export.postscript
deprecated("tkplot.export.postscript", tk_postscript)
#' @export tkplot.fit.to.screen
deprecated("tkplot.fit.to.screen", tk_fit)
#' @export tkplot.getcoords
deprecated("tkplot.getcoords", tk_coords)
#' @export tkplot.off
deprecated("tkplot.off", tk_off)
#' @export tkplot.reshape
deprecated("tkplot.reshape", tk_reshape)
#' @export tkplot.rotate
deprecated("tkplot.rotate", tk_rotate)
#' @export tkplot.setcoords
deprecated("tkplot.setcoords", tk_set_coords)

#' @export topological.sort
deprecated("topological.sort", topo_sort)
#' @export triad.census
deprecated("triad.census", triad_census)

#' @export unfold.tree
deprecated("unfold.tree", unfold_tree)

#' @export vcount
deprecated("vcount", gorder)
#' @export vertex.connectivity
deprecated("vertex.connectivity", vertex_connectivity)
#' @export vertex.disjoint.paths
deprecated("vertex.disjoint.paths", vertex_disjoint_paths)

#' @export walktrap.community
deprecated("walktrap.community", cluster_walktrap)
#' @export watts.strogatz.game
deprecated("watts.strogatz.game", sample_smallworld)
#' @export write.graph
deprecated("write.graph", write_graph)
#' @export graph.famous
deprecated("graph.famous", make_graph)
#' @export igraph.from.graphNEL
deprecated("igraph.from.graphNEL", graph_from_graphnel)
#' @export igraph.to.graphNEL
deprecated("igraph.to.graphNEL", as_graphnel)
#' @export getIgraphOpt
deprecated("getIgraphOpt", igraph_opt)
#' @export igraph.options
deprecated("igraph.options", igraph_options)
#' @export graph.intersection
deprecated("graph.intersection", intersection)
#' @export exportPajek
deprecated("exportPajek", export_pajek)
#' @export get.diameter
deprecated("get.diameter", get_diameter)
#' @export get.all.shortest.paths
deprecated("get.all.shortest.paths", all_shortest_paths)
#' @export get.shortest.paths
deprecated("get.shortest.paths", shortest_paths)
