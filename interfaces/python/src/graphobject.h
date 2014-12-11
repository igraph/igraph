/* -*- mode: C -*-  */
/* 
   IGraph library.
   Copyright (C) 2006-2012  Tamas Nepusz <ntamas@gmail.com>
   
   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
   
   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
   
   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA 
   02110-1301 USA

*/

#ifndef PYTHON_GRAPHOBJECT_H
#define PYTHON_GRAPHOBJECT_H

#include <Python.h>
#include <igraph.h>
#include "structmember.h"
#include "common.h"

extern PyTypeObject igraphmodule_GraphType;

/**
 * \ingroup python_interface
 * \brief A structure containing all the fields required to access an igraph from Python
 */
typedef struct 
{
  PyObject_HEAD
  // The graph object
  igraph_t g;
  // Python object to be called upon destruction
  PyObject* destructor;
  // Python object representing the sequence of vertices
  PyObject* vseq;
  // Python object representing the sequence of edges
  PyObject* eseq;
  // Python object of the weak reference list
  PyObject* weakreflist;
} igraphmodule_GraphObject;

void igraphmodule_Graph_init_internal(igraphmodule_GraphObject *self);
PyObject* igraphmodule_Graph_new(PyTypeObject *type, PyObject *args, PyObject *kwds);
int igraphmodule_Graph_clear(igraphmodule_GraphObject *self);
int igraphmodule_Graph_traverse(igraphmodule_GraphObject *self, visitproc visit, void *arg);
void igraphmodule_Graph_dealloc(igraphmodule_GraphObject* self);
int igraphmodule_Graph_init(igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_from_igraph_t(igraph_t *graph);
PyObject* igraphmodule_Graph_str(igraphmodule_GraphObject *self);

PyObject* igraphmodule_Graph_vcount(igraphmodule_GraphObject *self);
PyObject* igraphmodule_Graph_ecount(igraphmodule_GraphObject *self);
PyObject* igraphmodule_Graph_is_dag(igraphmodule_GraphObject *self);
PyObject* igraphmodule_Graph_is_directed(igraphmodule_GraphObject *self);
PyObject* igraphmodule_Graph_is_simple(igraphmodule_GraphObject *self);
PyObject* igraphmodule_Graph_add_vertices(igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_delete_vertices(igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_add_edges(igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_delete_edges(igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_degree(igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_is_loop(igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_count_multiple(igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_neighbors(igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_successors(igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_predecessors(igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_get_eid(igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds);

PyObject* igraphmodule_Graph_Adjacency(PyTypeObject *type, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_Asymmetric_Preference(PyTypeObject *type, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_Atlas(PyTypeObject *type, PyObject *args);
PyObject* igraphmodule_Graph_Barabasi(PyTypeObject *type, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_Degree_Sequence(PyTypeObject *type, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_Establishment(PyTypeObject *type, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_Erdos_Renyi(PyTypeObject *type, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_Famous(PyTypeObject *type, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_Forest_Fire(PyTypeObject *type, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_Full_Citation(PyTypeObject *type, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_Full(PyTypeObject *type, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_GRG(PyTypeObject *type, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_Growing_Random(PyTypeObject *type, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_Isoclass(PyTypeObject *type, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_Lattice(PyTypeObject *type, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_LCF(PyTypeObject *type, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_Preference(PyTypeObject *type, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_Recent_Degree(PyTypeObject *type, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_Ring(PyTypeObject *type, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_SBM(PyTypeObject *type, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_Star(PyTypeObject *type, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_Tree(PyTypeObject *type, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_Watts_Strogatz(PyTypeObject *type, PyObject *args, PyObject *kwds);

PyObject* igraphmodule_Graph_is_connected(igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_are_connected(igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_adjacency_spectral_embedding(igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_articulation_points(igraphmodule_GraphObject *self);
PyObject* igraphmodule_Graph_average_path_length(igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_betweenness(igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_bibcoupling(igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_closeness(igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_clusters(igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_cocitation(igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_constraint(igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_copy(igraphmodule_GraphObject *self);
PyObject* igraphmodule_Graph_decompose(igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_density(igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_diameter(igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_edge_betweenness(igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_eigen_adjacency(igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_get_shortest_paths(igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_get_all_shortest_paths(igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_maxdegree(igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_pagerank(igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_path_length_hist(igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_reciprocity(igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_rewire(igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_shortest_paths(igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_spanning_tree(igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_simplify(igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_subcomponent(igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_subgraph(igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_transitivity_undirected(igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_transitivity_local_undirected(igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds);

PyObject* igraphmodule_Graph_scan1(igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds);

PyObject* igraphmodule_Graph_layout_circle(igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_layout_sphere(igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_layout_random(igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_layout_random_3d(igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_layout_kamada_kawai(igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_layout_kamada_kawai_3d(igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_layout_drl(igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_layout_fruchterman_reingold(igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_layout_fruchterman_reingold_3d(igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_layout_grid_fruchterman_reingold(igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_layout_lgl(igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_layout_reingold_tilford(igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds);

PyObject* igraphmodule_Graph_get_adjacency(igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_get_edgelist(igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_to_undirected(igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_to_directed(igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds);

PyObject* igraphmodule_Graph_laplacian(igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds);

PyObject* igraphmodule_Graph_Read_DIMACS(PyTypeObject *type, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_Read_Edgelist(PyTypeObject *type, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_Read_GML(PyTypeObject *type, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_Read_Ncol(PyTypeObject *type, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_Read_Lgl(PyTypeObject *type, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_Read_Pajek(PyTypeObject *type, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_Read_GraphML(PyTypeObject *type, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_write_dimacs(igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_write_dot(igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_write_edgelist(igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_write_ncol(igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_write_lgl(igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_write_gml(igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_write_graphml(igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds);

PyObject* igraphmodule_Graph_isoclass(igraphmodule_GraphObject* self, PyObject* args, PyObject* kwds);
PyObject* igraphmodule_Graph_isomorphic(igraphmodule_GraphObject* self, PyObject* args, PyObject* kwds);
PyObject* igraphmodule_Graph_count_isomorphisms(igraphmodule_GraphObject* self, PyObject* args, PyObject* kwds);
PyObject* igraphmodule_Graph_get_isomorphisms(igraphmodule_GraphObject* self, PyObject* args, PyObject* kwds);
PyObject* igraphmodule_Graph_subisomorphic(igraphmodule_GraphObject* self, PyObject* args, PyObject* kwds);
PyObject* igraphmodule_Graph_count_subisomorphisms(igraphmodule_GraphObject* self, PyObject* args, PyObject* kwds);
PyObject* igraphmodule_Graph_get_subisomorphisms(igraphmodule_GraphObject* self, PyObject* args, PyObject* kwds);

Py_ssize_t igraphmodule_Graph_attribute_count(igraphmodule_GraphObject* self);
PyObject* igraphmodule_Graph_get_attribute(igraphmodule_GraphObject* self, PyObject* s);
int igraphmodule_Graph_set_attribute(igraphmodule_GraphObject* self, PyObject* k, PyObject* v);
PyObject* igraphmodule_Graph_attributes(igraphmodule_GraphObject* self);
PyObject* igraphmodule_Graph_vertex_attributes(igraphmodule_GraphObject* self);
PyObject* igraphmodule_Graph_edge_attributes(igraphmodule_GraphObject* self);

PyObject* igraphmodule_Graph_get_vertices(igraphmodule_GraphObject* self, void* closure);
PyObject* igraphmodule_Graph_get_edges(igraphmodule_GraphObject* self, void* closure);

PyObject* igraphmodule_Graph_complementer(igraphmodule_GraphObject* self, PyObject* args);
PyObject* igraphmodule_Graph_complementer_op(igraphmodule_GraphObject* self);
PyObject* igraphmodule_Graph_compose(igraphmodule_GraphObject* self, PyObject* other);
PyObject* igraphmodule_Graph_difference(igraphmodule_GraphObject* self, PyObject* other);
PyObject* igraphmodule_Graph_disjoint_union(igraphmodule_GraphObject* self, PyObject* other);
PyObject* igraphmodule_Graph_intersection(igraphmodule_GraphObject* self, PyObject* other);
PyObject* igraphmodule_Graph_union(igraphmodule_GraphObject* self, PyObject* other);

PyObject* igraphmodule_Graph_bfs(igraphmodule_GraphObject* self, PyObject* args, PyObject* kwds);
PyObject* igraphmodule_Graph_bfsiter(igraphmodule_GraphObject* self, PyObject* args, PyObject* kwds);

PyObject* igraphmodule_Graph_maxflow(igraphmodule_GraphObject* self, PyObject* args, PyObject* kwds);
PyObject* igraphmodule_Graph_maxflow_value(igraphmodule_GraphObject* self, PyObject* args, PyObject* kwds);
PyObject* igraphmodule_Graph_mincut(igraphmodule_GraphObject* self, PyObject* args, PyObject* kwds);
PyObject* igraphmodule_Graph_mincut_value(igraphmodule_GraphObject* self, PyObject* args, PyObject* kwds);

PyObject* igraphmodule_Graph_cliques(igraphmodule_GraphObject* self, PyObject* args, PyObject* kwds);
PyObject* igraphmodule_Graph_maximal_cliques(igraphmodule_GraphObject* self, PyObject* args, PyObject* kwds);
PyObject* igraphmodule_Graph_largest_cliques(igraphmodule_GraphObject* self);
PyObject* igraphmodule_Graph_clique_number(igraphmodule_GraphObject* self);
PyObject* igraphmodule_Graph_independent_sets(igraphmodule_GraphObject* self, PyObject* args, PyObject* kwds);
PyObject* igraphmodule_Graph_maximal_independent_sets(igraphmodule_GraphObject* self);
PyObject* igraphmodule_Graph_largest_independent_sets(igraphmodule_GraphObject* self);
PyObject* igraphmodule_Graph_independence_number(igraphmodule_GraphObject* self);

PyObject* igraphmodule_Graph_community_edge_betweenness(igraphmodule_GraphObject* self, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_community_fastgreedy(igraphmodule_GraphObject* self, PyObject *args, PyObject *kwds);
PyObject *igraphmodule_Graph_community_infomap(igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_community_label_propagation(igraphmodule_GraphObject* self, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_community_leading_eigenvector(igraphmodule_GraphObject* self, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_community_leading_eigenvector_naive(igraphmodule_GraphObject* self, PyObject *args, PyObject *kwds);
PyObject *igraphmodule_Graph_community_multilevel(igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds);
PyObject *igraphmodule_Graph_community_optimal_modularity(igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_community_spinglass(igraphmodule_GraphObject* self, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_community_walktrap(igraphmodule_GraphObject* self, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_modularity(igraphmodule_GraphObject* self, PyObject *args, PyObject *kwds);

PyObject *igraphmodule_Graph_is_bipartite(igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds);

PyObject* igraphmodule_Graph___graph_as_cobject__(igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph___register_destructor__(igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds);

#endif
