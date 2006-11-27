/* -*- mode: C -*-  */
/* 
   IGraph library.
   Copyright (C) 2006  Gabor Csardi <csardi@rmki.kfki.hu>
   MTA RMKI, Konkoly-Thege Miklos st. 29-33, Budapest 1121, Hungary
   
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
#include "igraph.h"
#include "types.h"
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
PyObject* igraphmodule_Graph_str(igraphmodule_GraphObject *self);

PyObject* igraphmodule_Graph_vcount(igraphmodule_GraphObject *self);
PyObject* igraphmodule_Graph_ecount(igraphmodule_GraphObject *self);
PyObject* igraphmodule_Graph_is_directed(igraphmodule_GraphObject *self);
PyObject* igraphmodule_Graph_add_vertices(igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_delete_vertices(igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_add_edges(igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_delete_edges(igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_degree(igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_neighbors(igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_successors(igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_predecessors(igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_get_eid(igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds);

PyObject* igraphmodule_Graph_Adjacency(PyTypeObject *type, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_Atlas(PyTypeObject *type, PyObject *args);
PyObject* igraphmodule_Graph_Barabasi(PyTypeObject *type, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_Establishment(PyTypeObject *type, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_Erdos_Renyi(PyTypeObject *type, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_Full(PyTypeObject *type, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_Growing_Random(PyTypeObject *type, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_Lattice(PyTypeObject *type, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_Star(PyTypeObject *type, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_Ring(PyTypeObject *type, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_Tree(PyTypeObject *type, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_Degree_Sequence(PyTypeObject *type, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_Isoclass(PyTypeObject *type, PyObject *args, PyObject *kwds);

PyObject* igraphmodule_Graph_is_connected(igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_are_connected(igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_average_path_length(igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_betweenness(igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_bibcoupling(igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_closeness(igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_clusters(igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_cocitation(igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_copy(igraphmodule_GraphObject *self);
PyObject* igraphmodule_Graph_decompose(igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_density(igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_diameter(igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_edge_betweenness(igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_get_shortest_paths(igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_get_all_shortest_paths(igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_maxdegree(igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_pagerank(igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_reciprocity(igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_rewire(igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_shortest_paths(igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_spanning_tree(igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_simplify(igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_subcomponent(igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_subgraph(igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_transitivity_undirected(igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_transitivity_local_undirected(igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds);

PyObject* igraphmodule_Graph_layout_circle(igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_layout_sphere(igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_layout_random(igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_layout_random_3d(igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_layout_kamada_kawai(igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_layout_kamada_kawai_3d(igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_layout_fruchterman_reingold(igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_layout_fruchterman_reingold_3d(igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_layout_grid_fruchterman_reingold(igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_layout_lgl(igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_layout_reingold_tilford(igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds);

PyObject* igraphmodule_Graph_get_adjacency(igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_get_edgelist(igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_to_undirected(igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_to_directed(igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds);

PyObject* igraphmodule_Graph_Read_Edgelist(PyTypeObject *type, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_Read_Ncol(PyTypeObject *type, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_Read_Lgl(PyTypeObject *type, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_Read_Pajek(PyTypeObject *type, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_Read_GraphML(PyTypeObject *type, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_write_edgelist(igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_write_ncol(igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_write_lgl(igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_write_graphml(igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds);

PyObject* igraphmodule_Graph_isoclass(igraphmodule_GraphObject* self, PyObject* args, PyObject* kwds);
PyObject* igraphmodule_Graph_isomorphic(igraphmodule_GraphObject* self, PyObject* args, PyObject* kwds);

int igraphmodule_Graph_attribute_count(igraphmodule_GraphObject* self);
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

PyObject* igraphmodule_Graph_maxflow_value(igraphmodule_GraphObject* self, PyObject* args, PyObject* kwds);
PyObject* igraphmodule_Graph_mincut_value(igraphmodule_GraphObject* self, PyObject* args, PyObject* kwds);

PyObject* igraphmodule_Graph___graph_as_cobject__(igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph___register_destructor__(igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds);

/** \ingroup python_interface
 * \brief Member list of the \c igraph.Graph object type
 */
#define OFF(x) offsetof(igraphmodule_GraphObject, x)

static PyGetSetDef igraphmodule_Graph_getseters[] = {
  {"vs", (getter)igraphmodule_Graph_get_vertices, NULL,
      "The sequence of vertices in the graph.", NULL
  },
  {"es", (getter)igraphmodule_Graph_get_edges, NULL,
      "The sequence of edges in the graph.", NULL
  },
  {NULL}
};
/*
static PyMemberDef igraphmodule_Graph_members[] = {
  {"vertices", T_OBJECT, OFF(vseq), RO,
      "Sequence of vertices in the graph"
  },
  {"vs", T_OBJECT, OFF(vseq), RO,
      "Sequence of vertices in the graph. Alias for 'vertices'."
  },
  {"nodes", T_OBJECT, OFF(vseq), RO,
      "Sequence of vertices in the graph. Alias for 'nodes'."
  },
  {NULL}
};*/

/** \ingroup python_interface
 * \brief Method list of the \c igraph.Graph object type
 */
static PyMethodDef igraphmodule_Graph_methods[] = 
{
  ////////////////////////////
  // BASIC IGRAPH INTERFACE //
  ////////////////////////////
  
  // interface to igraph_vcount
  {"vcount", (PyCFunction)igraphmodule_Graph_vcount,
      METH_NOARGS,
      "vcount()\n\n"
      "Counts the number of vertices.\n"
      "@return: the number of vertices in the graph.\n"
      "@rtype: integer"
  },
   
  // interface to igraph_ecount
  {"ecount", (PyCFunction)igraphmodule_Graph_ecount,
      METH_NOARGS,
      "ecount()\n\n"
      "Counts the number of edges.\n"
      "@return: the number of edges in the graph.\n"
      "@rtype: integer"
  },
   
  // interface to igraph_is_directed
  {"is_directed", (PyCFunction)igraphmodule_Graph_is_directed,
      METH_NOARGS,
      "is_directed()\n\n"
      "Checks whether the graph is directed."
      "@return: C{True} if it is directed, C{False} otherwise.\n"
      "@rtype: boolean"
  },
   
  // interface to igraph_add_vertices
  {"add_vertices", (PyCFunction)igraphmodule_Graph_add_vertices,
      METH_VARARGS,
      "add_vertices(n)\n\n"
      "Adds vertices to the graph.\n\n"
      "@param n: the number of vertices to be added\n"
      "@return: the same graph object\n"
  },
   
  // interface to igraph_delete_vertices
  {"delete_vertices", (PyCFunction)igraphmodule_Graph_delete_vertices,
      METH_VARARGS,
      "delete_vertices(vs)\n\n"
      "Deletes vertices and all its edges from the graph.\n\n"
      "@param vs: a single vertex ID or the list of vertex IDs\n"
      "  to be deleted.\n"
      "@return: the same graph object\n"
  },
   
  // interface to igraph_add_edges
  {"add_edges", (PyCFunction)igraphmodule_Graph_add_edges,
      METH_VARARGS,
      "add_edges(es)\n\n"
      "Adds edges to the graph.\n\n"
      "@param es: the list of edges to be added. Every edge is\n"
      "  represented with a tuple, containing the vertex IDs of the\n"
      "  two endpoints. Vertices are enumerated from zero. It is\n"
      "  allowed to provide a single pair instead of a list consisting\n"
      "  of only one pair.\n"
      "@return: the same graph object\n"
  },
   
  // interface to igraph_delete_edges
  {"delete_edges", (PyCFunction)igraphmodule_Graph_delete_edges,
      METH_VARARGS,
      "delete_edges(es)\n\n"
      "Removes edges from the graph.\n\n"
      "@param es: the list of edges to be removed. Every edge is\n"
      "  represented with a tuple, containing the vertex IDs of the\n"
      "  two endpoints. Vertices are enumerated from zero. It is\n"
      "  allowed to provide a single pair instead of a list consisting\n"
      "  of only one pair. Nonexistent edges will be silently ignored.\n"
      "  All vertices will be kept, even if they lose all their edges.\n"
      "@return: the same graph object\n"
  },
   
  // interface to igraph_degree
  {"degree", (PyCFunction)igraphmodule_Graph_degree,
      METH_VARARGS | METH_KEYWORDS,
      "degree(vertices, type=ALL, loops=False)\n\n"
      "Returns some vertex degrees from the graph.\n\n"
      "This method accepts a single vertex ID or a list of vertex IDs as a\n"
      "parameter, and returns the degree of the given vertices (in the\n"
      "form of a single integer or a list, depending on the input\n"
      "parameter).\n"
      "\n"
      "@param vertices: a single vertex ID or a list of vertex IDs\n"
      "@param type: the type of degree to be returned (L{OUT} for\n"
      "  out-degrees, L{IN} IN for in-degrees or L{ALL} for the sum of\n"
      "  them).\n"
      "@param loops: whether self-loops should be counted.\n"
  },
   
  // interfaces to igraph_neighbors
  {"neighbors", (PyCFunction)igraphmodule_Graph_neighbors,
      METH_VARARGS | METH_KEYWORDS,
      "neighbors(vertex, type=ALL)\n\n"
      "Returns adjacent vertices to a given vertex.\n\n"
      "@param vertex: a vertex ID\n"
      "@param type: whether to return only predecessors (L{OUT}),\n"
      "  successors (L{OUT}) or both (L{ALL}). Ignored for undirected\n"
      "  graphs."
  },
   
  {"successors", (PyCFunction)igraphmodule_Graph_successors,
      METH_VARARGS | METH_KEYWORDS,
      "successors(vertex)\n\n"
      "Returns the successors of a given vertex.\n\n"
      "Equivalent to calling the L{Graph.neighbors} method with type=L{OUT}."
  },
   
  {"predecessors", (PyCFunction)igraphmodule_Graph_predecessors,
      METH_VARARGS | METH_KEYWORDS,
      "predecessors(vertex)\n\n"
      "Returns the predecessors of a given vertex.\n\n"
      "Equivalent to calling the L{Graph.neighbors} method with type=L{IN}."
  },

  /* interface to igraph_get_eid */
  {"get_eid", (PyCFunction)igraphmodule_Graph_get_eid,
   METH_VARARGS | METH_KEYWORDS,
   "get_eid(v1, v2)\n\n"
   "Returns the edge ID of an arbitrary edge between vertices v1 and v2\n\n"
   "@param v1: the first vertex ID\n"
   "@param v2: the second vertex ID\n"
   "@return: the edge ID of an arbitrary edge between vertices v1 and v2\n"
  },

  //////////////////////
  // GRAPH GENERATORS //
  //////////////////////

  // interface to igraph_adjacency
  {"Adjacency", (PyCFunction)igraphmodule_Graph_Adjacency,
      METH_CLASS | METH_VARARGS | METH_KEYWORDS,
      "Adjacency(matrix, mode=ADJ_DIRECTED)\n\n"
      "Generates a graph from its adjacency matrix.\n\n"
      "@param matrix: the adjacency matrix\n"
      "@param mode: the mode to be used. Possible values are:\n"
      "\n"
      "  - C{ADJ_DIRECTED} - the graph will be directed and a matrix\n"
      "    element gives the number of edges between two vertex.\n"
      "  - C{ADJ_UNDIRECTED} - alias to C{ADJ_MAX} for convenience.\n"
      "  - C{ADJ_MAX}   - undirected graph will be created and the number of\n"
      "    edges between vertex M{i} and M{j} is M{max(A(i,j), A(j,i))}\n"
      "  - C{ADJ_MIN}   - like C{ADJ_MAX}, but with M{min(A(i,j), A(j,i))}\n"
      "  - C{ADJ_PLUS}  - like C{ADJ_MAX}, but with M{A(i,j) + A(j,i)}\n"
      "  - C{ADJ_UPPER} - undirected graph with the upper right triangle of\n"
      "    the matrix (including the diagonal)\n"
      "  - C{ADJ_LOWER} - undirected graph with the lower left triangle of\n"
      "    the matrix (including the diagonal)\n"
      " Optional, defaults to ADJ_DIRECTED.\n"
  },
  
  // interface to igraph_atlas
  {"Atlas", (PyCFunction)igraphmodule_Graph_Atlas,
      METH_CLASS | METH_KEYWORDS,
      "Atlas(idx)\n\n"
      "Generates a graph from the Graph Atlas.\n\n"
      "@param idx: The index of the graph to be generated.\n"
      "  Indices start from zero, graphs are listed:\n\n"
      "    1. in increasing order of number of nodes;\n"
      "    2. for a fixed number of nodes, in increasing order of the\n"
      "       number of edges;\n"
      "    3. for fixed numbers of nodes and edges, in increasing order\n"
      "       of the degree sequence, for example 111223 < 112222;\n"
      "    4. for fixed degree sequence, in increasing number of automorphisms.\n"
  },
	
  // interface to igraph_barabasi_game
  {"Barabasi", (PyCFunction)igraphmodule_Graph_Barabasi,
      METH_VARARGS | METH_CLASS | METH_KEYWORDS,
      "Barabasi(n, m, outpref=False, directed=False, power=1)\n\n"
      "Generates a graph based on the Barabasi-Albert model.\n\n"
      "@param n: the number of vertices\n"
      "@param m: either the number of outgoing edges generated for\n"
      "  each vertex or a list containing the number of outgoing\n"
      "  edges for each vertex explicitly.\n"
      "@param outpref: C{True} if the out-degree of a given vertex\n"
      "  should also increase its citation probability (as well as\n"
      "  its in-degree), but it defaults to C{False}.\n"
      "@param directed: C{True} if the generated graph should be\n"
      "  directed (default: C{False}).\n"
      "@param power: the power constant of the nonlinear model.\n"
      "  It can be omitted, and in this case the usual linear model\n"
      "  will be used.\n"
  },

  // interface to igraph_establishment_game
  {"Establishment", (PyCFunction)igraphmodule_Graph_Establishment,
      METH_VARARGS | METH_CLASS | METH_KEYWORDS,
      "Establishment(n, k, type_dist, pref_matrix, directed=False)\n\n"
      "Generates a graph based on a simple growing model with vertex types.\n\n"
      "A single vertex is added at each time step. This new vertex tries to\n"
      "connect to k vertices in the graph. The probability that such a\n"
      "connection is realized depends on the types of the vertices involved.\n"
      "\n"
      "@param n: the number of vertices in the graph\n"
      "@param k: the number of connections tried in each step\n"
      "@param type_dist: list giving the distribution of vertex types\n"
      "@param pref_matrix: matrix (list of lists) giving the connection\n"
      "  probabilities for different vertex types\n"
      "@param directed: whether to generate a directed graph.\n"
  },
  
  // interface to igraph_erdos_renyi_game
  {"Erdos_Renyi", (PyCFunction)igraphmodule_Graph_Erdos_Renyi,
      METH_VARARGS | METH_CLASS | METH_KEYWORDS,
      "Erdos_Renyi(n, p, m, directed=False, loops=False)\n\n"
      "Generates a graph based on the Erdos-Renyi model.\n\n"
      "@param n: the number of vertices.\n"
      "@param p: the probability of edges. If given, C{m} must be missing.\n"
      "@param m: the number of edges. If given, C{p} must be missing.\n"
      "@param directed: whether to generate a directed graph.\n"
      "@param loops: whether self-loops are allowed.\n"
  },
  
  // interface to igraph_full_game
  {"Full", (PyCFunction)igraphmodule_Graph_Full,
      METH_VARARGS | METH_CLASS | METH_KEYWORDS,
      "Full(n, directed=False, loops=False)\n\n"
      "Generates a full graph (directed or undirected, with or without loops).\n\n"
      "@param n: the number of vertices.\n"
      "@param directed: whether to generate a directed graph.\n"
      "@param loops: whether self-loops are allowed.\n"
  },
  
  // interface to igraph_growing_random_game
  {"Growing_Random", (PyCFunction)igraphmodule_Graph_Growing_Random,
      METH_VARARGS | METH_CLASS | METH_KEYWORDS,
      "Growing_Random(n, m, directed=False, citation=False)\n\n"
      "Generates a growing random graph.\n\n"
      "@param n: The number of vertices in the graph\n"
      "@param m: The number of edges to add in each step (after adding a new vertex)\n"
      "@param directed: whether the graph should be directed.\n"
      "@param citation: whether the new edges should originate from the most\n"
      "   recently added vertex.\n"
  },
  
  // interface to igraph_star
  {"Star", (PyCFunction)igraphmodule_Graph_Star,
      METH_VARARGS | METH_CLASS | METH_KEYWORDS,
      "Star(n, mode=STAR_UNDIRECTED, center=0)\n\n"
      "Generates a star graph.\n\n"
      "@param n: the number of vertices in the graph\n"
      "@param mode: Gives the type of the star graph to create. Should be\n"
      "  one of the constants C{STAR_OUT}, C{STAR_IN} and C{STAR_UNDIRECTED}.\n"
      "@param center: Vertex ID for the central vertex in the star.\n"
  },
  
  // interface to igraph_lattice
  {"Lattice", (PyCFunction)igraphmodule_Graph_Lattice,
      METH_VARARGS | METH_CLASS | METH_KEYWORDS,
      "Lattice(dim, nei=1, directed=False, mutual=True, circular=True)\n\n"
      "Generates a regular lattice.\n\n"
      "@param dim: list with the dimensions of the lattice\n"
      "@param nei: value giving the distance (number of steps) within which\n"
      "   two vertices will be connected. Not implemented yet.\n"
      "@param directed: whether to create a directed graph.\n"
      "@param mutual: whether to create all connections as mutual\n"
      "    in case of a directed graph.\n"
      "@param circular: whether the generated lattice is periodic.\n"
  },
  
  // interface to igraph_ring
  {"Ring", (PyCFunction)igraphmodule_Graph_Ring,
      METH_VARARGS | METH_CLASS | METH_KEYWORDS,
      "Ring(n, directed=False, mutual=False, circular=True)\n\n"
      "Generates a ring graph.\n\n"
      "@param n: the number of vertices in the ring\n"
      "@param directed: whether to create a directed ring.\n"
      "@param mutual: whether to create mutual edges in a directed ring.\n"
      "@param circular: whether to create a closed ring.\n"
  },
  
  // interface to igraph_tree
  {"Tree", (PyCFunction)igraphmodule_Graph_Tree,
      METH_VARARGS | METH_CLASS | METH_KEYWORDS,
      "Tree(n, children, type=TREE_UNDIRECTED)\n\n"
      "Generates a tree in which almost all vertices have the same number of children.\n\n"
      "@param n: the number of vertices in the graph\n"
      "@param children: the number of children of a vertex in the graph\n"
      "@param type: determines whether the tree should be directed, and if\n"
      "  this is the case, also its orientation. Must be one of\n"
      "  C{TREE_IN}, C{TREE_OUT} and C{TREE_UNDIRECTED}.\n"
  },
  
  // interface to igraph_degree_sequence_game
  {"Degree_Sequence", (PyCFunction)igraphmodule_Graph_Degree_Sequence,
      METH_VARARGS | METH_CLASS | METH_KEYWORDS,
      "Degree_Sequence(out, in=None)\n\n"
      "Generates a graph with a given degree sequence.\n\n"
      "@param out: the out-degree sequence for a directed graph. If the\n"
      "  in-degree sequence is omitted, the generated graph\n"
      "  will be undirected, so this will be the in-degree\n"
      "  sequence as well\n"
      "@param in: the in-degree sequence for a directed graph.\n"
      "   If omitted, the generated graph will be undirected.\n"
  },
  
  // interface to igraph_isoclass_create
  {"Isoclass", (PyCFunction)igraphmodule_Graph_Isoclass,
      METH_VARARGS | METH_CLASS | METH_KEYWORDS,
      "Isoclass(n, class, directed=False)\n\n"
      "Generates a graph with a given isomorphy class.\n\n"
      "@param n: the number of vertices in the graph (3 or 4)\n"
      "@param class: the isomorphy class\n"
      "@param directed: whether the graph should be directed.\n"
  },
  
  /////////////////////////////////////
  // STRUCTURAL PROPERTIES OF GRAPHS //
  /////////////////////////////////////
  
  // interface to igraph_are_connected
  {"are_connected", (PyCFunction)igraphmodule_Graph_are_connected,
      METH_VARARGS | METH_KEYWORDS,
      "are_connected(v1, v2)\n\n"
      "Decides whether two given vertices are directly connected.\n\n"
      "@param v1: the first vertex\n"
      "@param v2: the second vertex\n"
      "@return: C{True} if there exists an edge from v1 to v2, C{False}\n"
      "  otherwise.\n"
  },
  
  // interface to igraph_average_path_length
  {"average_path_length", (PyCFunction)igraphmodule_Graph_average_path_length,
      METH_VARARGS | METH_KEYWORDS,
      "average_path_length(directed=True, unconn=True)\n\n"
      "Calculates the average path length in a graph.\n\n"
      "@param directed: whether to consider directed paths in case of a\n"
      "  directed graph. Ignored for undirected graphs.\n"
      "@param unconn: what to do when the graph is unconnected. If C{True},\n"
      "  the average of the geodesic lengths in the components is\n"
      "  calculated. Otherwise for all unconnected vertex pairs,\n"
      "  a path length equal to the number of vertices is used.\n"
      "@return: the average path length in the graph\n"
  },
  
  // interface to igraph_betweenness
  {"betweenness", (PyCFunction)igraphmodule_Graph_betweenness,
      METH_VARARGS | METH_KEYWORDS,
      "betweenness(vertices=None, directed=True)\n\n"
      "Calculates the betweenness of nodes in a graph.\n\n"
      "Keyword arguments:\n"
      "@param vertices: the vertices for which the betweennesses must be returned.\n"
      "  If C{None}, assumes all of the vertices in the graph.\n"
      "@param directed: whether to consider directed paths.\n"
      "@return: the betweenness of the given nodes in a list\n"
  },
  
  // interface to igraph_bibcoupling
  {"bibcoupling", (PyCFunction)igraphmodule_Graph_bibcoupling,
      METH_VARARGS | METH_KEYWORDS,
      "bibcoupling(vertices)\n\n"
      "Calculates bibliographic coupling values for given vertices\n"
      "in a graph.\n\n"
      "@param vertices: the vertices to be analysed.\n"
      "@return: bibliographic coupling values for all given\n"
      "  vertices in a matrix.\n"
  },
  
  // interface to igraph_closeness
  {"closeness", (PyCFunction)igraphmodule_Graph_closeness,
      METH_VARARGS | METH_KEYWORDS,
      "closeness(vertices=None, mode=ALL)\n\n"
      "Calculates the closeness centralities of given nodes in a graph.\n\n"
      "The closeness centerality of a vertex measures how easily other\n"
      "vertices can be reached from it (or the other way: how easily it\n"
      "can be reached from the other vertices). It is defined as the\n"
      "number of the number of vertices minus one divided by the sum of\n"
      "the lengths of all geodesics from/to the given vertex.\n\n"
      "If the graph is not connected, and there is no path between two\n"
      "vertices, the number of vertices is used instead the length of\n"
      "the geodesic. This is always longer than the longest possible\n"
      "geodesic.\n\n"
      "@param vertices: the vertices for which the closenesses must\n"
      "  be returned. If C{None}, uses all of the vertices in the graph.\n"
      "@param mode: must be one of C{IN}, C{OUT} and C{ALL}. C{IN} means\n"
      "  that the length of the incoming paths, C{OUT} means that the\n"
      "  length of the outgoing paths must be calculated. C{ALL} means\n"
      "  that both of them must be calculated.\n"
      "@return: the calculated closenesses in a list\n"
  },   	    
  
  // interface to igraph_clusters
  {"clusters", (PyCFunction)igraphmodule_Graph_clusters,
      METH_VARARGS | METH_KEYWORDS,
      "clusters(mode=STRONG)\n\n"
      "Calculates the (strong or weak) clusters for a given graph.\n\n"
      "@param mode: must be either C{STRONG} or C{WEAK}, depending on\n"
      "  the clusters being sought. Optional, defaults to C{STRONG}.\n"
      "@return: the component index for every node in the graph.\n"
  },
  {"components", (PyCFunction)igraphmodule_Graph_clusters,
      METH_VARARGS | METH_KEYWORDS,
      "components(mode=STRONG)\n\n"
      "Alias for L{Graph.clusters}.\n\n"
      "See the documentation of L{Graph.clusters} for details."
  },
  {"copy", (PyCFunction)igraphmodule_Graph_copy,
      METH_NOARGS,
      "copy()\n\n"
      "Creates an exact deep copy of the graph."
  },
  {"decompose", (PyCFunction)igraphmodule_Graph_decompose,
      METH_VARARGS | METH_KEYWORDS,
      "decompose(mode=STRONG, maxcompno=None, minelements=1)\n\n"
      "Decomposes the graph into subgraphs.\n\n"
      "@param mode: must be either STRONG or WEAK, depending on the\n"
      "  clusters being sought.\n"
      "@param maxcompno: maximum number of components to return.\n"
      "  C{None} means all possible components.\n"
      "@param minelements: minimum number of vertices in a component.\n"
      "  By setting this to 2, isolated vertices are not returned\n"
      "  as separate components.\n"
      "@return: a list of the subgraphs. Every returned subgraph is a\n"
      "  copy of the original.\n"
  },
  
  // interface to igraph_cocitation
  {"cocitation", (PyCFunction)igraphmodule_Graph_cocitation,
      METH_VARARGS | METH_KEYWORDS,
      "cocitation(vertices)\n\n"
      "Calculates cocitation scores for given vertices in a graph.\n\n"
      "@param vertices: the vertices to be analysed.\n"
      "@return: cocitation scores for all given vertices in a matrix."
  },
  
  /* interface to igraph_density */
  {"density", (PyCFunction)igraphmodule_Graph_density,
      METH_VARARGS | METH_KEYWORDS,
       "density(loops=False)\n\n"
       "Calculates the density of the graph.\n\n"
       "@param loops: whether to take loops into consideration. If C{True},\n"
       "  the algorithm assumes that there might be some loops in the graph\n"
       "  and calculates the density accordingly. If C{False}, the algorithm\n"
       "  assumes that there can't be any loops.\n"
      "@return: the reciprocity of the graph."
  },
  
  // interface to igraph_diameter
  {"diameter", (PyCFunction)igraphmodule_Graph_diameter,
      METH_VARARGS | METH_KEYWORDS,
      "diameter(directed=True, unconn=True)\n\n"
      "Calculates the diameter of the graph.\n\n"
      "@param directed: whether to consider directed paths.\n"
      "@param unconn: if C{True} and the graph is undirected, the\n"
      "  longest geodesic within a component will be returned. If\n"
      "  C{False} and the graph is undirected, the result is the\n"
      "  number of vertices."
  },
  
  // interface to igraph_edge_betweenness
  {"edge_betweenness", (PyCFunction)igraphmodule_Graph_edge_betweenness,
      METH_VARARGS | METH_KEYWORDS,
      "edge_betweenness(directed=True)\n\n"
      "Calculates the edge betweennesses in a graph.\n\n"
      "@param directed: whether to consider directed paths.\n"
      "@return: a list with the edge betweennesses of all specified edges.\n"
  },
  
  // interface to igraph_get_shortest_paths
  {"get_shortest_paths", (PyCFunction)igraphmodule_Graph_get_shortest_paths,
      METH_VARARGS | METH_KEYWORDS,
      "get_shortest_paths(v, mode=OUT)\n\n"
      "Calculates the shortest paths from/to a given node in a graph.\n\n"
      "@param v: the source/destination for the calculated paths\n"
      "@param mode: the directionality of the paths. C{IN} means to\n"
      "  calculate incoming paths, C{OUT} means to calculate outgoing\n"
      "  paths, C{ALL} means to calculate both ones.\n"
      "@return: at most one shortest path for every node in the graph in a\n"
      "list. For unconnected graphs, some of the list elements will be\n"
      "empty lists. Note that in case of mode=C{IN}, the nodes in a path\n"
      "are returned in reversed order!"
  },
  
  // interface to igraph_get_all_shortest_paths
  {"get_all_shortest_paths", (PyCFunction)igraphmodule_Graph_get_all_shortest_paths,
      METH_VARARGS | METH_KEYWORDS,
      "get_all_shortest_paths(v, mode=OUT)\n\n"
      "Calculates all of the shortest paths from/to a given node in a graph.\n\n"
      "@param v: the source/destination for the calculated paths\n"
      "@param mode: the directionality of the paths. C{IN} means to calculate\n"
      "  incoming paths, C{OUT} means to calculate outgoing paths,\n"
      "  C{ALL} means to calculate both ones.\n"
      "@return: all of the shortest path from the given node to every other\n"
      "reachable node in the graph in a list. Note that in case of mode=C{IN},\n"
      "the nodes in a path are returned in reversed order!"
  },
  
  // interface to igraph_is_connected
  {"is_connected", (PyCFunction)igraphmodule_Graph_is_connected,
      METH_VARARGS | METH_KEYWORDS,
      "is_connected(mode=STRONG)\n\n"
      "Decides whether a graph is connected.\n\n"
      "@param mode: whether we should calculate strong or weak connectivity.\n"
      "@return: C{True} if the graph is connected, C{False} otherwise.\n"
  },
  
  /* interface to igraph_maxdegree */
  {"maxdegree", (PyCFunction)igraphmodule_Graph_maxdegree,
      METH_VARARGS | METH_KEYWORDS,
      "maxdegree(vertices=None, type=ALL, loops=False)\n\n"
      "Returns the maximum degree of a vertex set in the graph.\n\n"
      "This method accepts a single vertex ID or a list of vertex IDs as a\n"
      "parameter, and returns the degree of the given vertices (in the\n"
      "form of a single integer or a list, depending on the input\n"
      "parameter).\n"
      "\n"
      "@param vertices: a single vertex ID or a list of vertex IDs or\n"
      "  C{None} meaning all the vertices in the graph.\n"
      "@param type: the type of degree to be returned (L{OUT} for\n"
      "  out-degrees, L{IN} IN for in-degrees or L{ALL} for the sum of\n"
      "  them).\n"
      "@param loops: whether self-loops should be counted.\n"
  },
   
  // interface to igraph_pagerank
  {"pagerank", (PyCFunction)igraphmodule_Graph_pagerank,
      METH_VARARGS | METH_KEYWORDS,
      "pagerank(vertices=None, directed=True, niter=1000, eps=0.001, damping=0.85)\n\n"
      "Calculates the Google PageRank values of a graph.\n\n"
      "@param vertices: the indices of the vertices being queried.\n"
      "  C{None} means all of the vertices.\n"
      "@param directed: whether to consider directed paths.\n"
      "@param niter: the maximum number of iterations to be performed.\n"
      "@param eps: the iteration stops if all of the PageRank values change\n"
      "  less than M{eps} between two iterations.\n"
      "@param damping: the damping factor.\n"
      "  M{1-damping} is the PageRank value for nodes with no\n"
      "  incoming links.\n"
      "@return: a list with the Google PageRank values of the specified\n"
      "  vertices."
  },

  // interface to igraph_reciprocity
  {"reciprocity", (PyCFunction)igraphmodule_Graph_reciprocity,
      METH_VARARGS | METH_KEYWORDS,
      "reciprocity()\n\n"
      "@return: the reciprocity of the graph."
  },
  
  // interface to igraph_rewire
  {"rewire", (PyCFunction)igraphmodule_Graph_rewire,
      METH_VARARGS | METH_KEYWORDS,
      "rewire(n=1000, mode=REWIRING_SIMPLE)\n\n"
      "Randomly rewires the graph while preserving the degree distribution.\n\n"
      "Please note that the rewiring is done \"in-place\", so the original\n"
      "graph will be modified. If you want to preserve the original graph,\n"
      "use the L{copy} method before.\n\n"
      "@param n: the number of rewiring trials.\n"
      "@param mode: the rewiring algorithm to use. As for now, only\n"
      "  C{REWIRING_SIMPLE} is supported.\n"
      "@return: the modified graph.\n"
  },
  
  // interface to igraph_shortest_paths
  {"shortest_paths", (PyCFunction)igraphmodule_Graph_shortest_paths,
      METH_VARARGS | METH_KEYWORDS,
      "shortest_paths(vertices, mode=OUT)\n\n"
      "Calculates shortest path lengths for given nodes in a graph.\n\n"
      "@param vertices: a list containing the vertex IDs which should be\n"
      "  included in the result.\n"
      "@param mode: the type of shortest paths to be used for the\n"
      "  calculation in directed graphs. C{OUT} means only outgoing,\n"
      "  C{IN} means only incoming paths. C{ALL} means to consider\n"
      "  the directed graph as an undirected one.\n"
      "@return: the shortest path lengths for given nodes in a matrix\n"
  },
  
  // interface to igraph_simplify
  {"simplify", (PyCFunction)igraphmodule_Graph_simplify,
      METH_VARARGS | METH_KEYWORDS,
      "simplify(multiple=True, loops=True)\n\n"
      "Simplifies a graph by removing self-loops and/or multiple edges.\n\n"
      "@param multiple: whether to remove multiple edges.\n"
      "@param loops: whether to remove loops.\n"
  },
  
  // interface to igraph_minimum_spanning_tree_unweighted and
  // igraph_minimum_spanning_tree_prim
  {"spanning_tree", (PyCFunction)igraphmodule_Graph_spanning_tree,
      METH_VARARGS | METH_KEYWORDS,
      "spanning_tree(weights=None)\n\n"
      "Calculates a minimum spanning tree for a graph (weighted or unweighted)\n\n"
      "@param weights: a vector containing weights for every edge in\n"
      "  the graph. C{None} means that the graph is unweighted.\n"
      "@return: the spanning tree as an igraph.Graph object."
  },
  
  // interface to igraph_subcomponent
  {"subcomponent", (PyCFunction)igraphmodule_Graph_subcomponent,
      METH_VARARGS | METH_KEYWORDS,
      "subcomponent(v, mode=ALL)\n\n"
      "Determines the indices of vertices which are in the same component as a given vertex.\n\n"
      "@param v: the index of the vertex used as the source/destination\n"
      "@param mode: if equals to C{IN}, returns the vertex IDs from\n"
      "  where the given vertex can be reached. If equals to C{OUT},\n"
      "  returns the vertex IDs which are reachable from the given\n"
      "  vertex. If equals to C{ALL}, returns all vertices within the\n"
      "  same component as the given vertex, ignoring edge directions.\n"
      "  Note that this is not equal to calculating the union of the \n"
      "  results of C{IN} and C{OUT}.\n"
      "@return: the indices of vertices which are in the same component as a given vertex.\n"
  },
  
  // interface to igraph_subgraph
  {"subgraph", (PyCFunction)igraphmodule_Graph_subgraph,
      METH_VARARGS | METH_KEYWORDS,
      "subgraph(vertices)\n\n"
      "Returns a subgraph based on the given vertices.\n\n"
      "@param vertices: a list containing the vertex IDs which\n"
      "  should be included in the result.\n"
      "@return: a copy of the subgraph\n"
  },
  
  // interface to igraph_transitivity_undirected
  {"transitivity_undirected", (PyCFunction)igraphmodule_Graph_transitivity_undirected,
      METH_VARARGS | METH_KEYWORDS,
      "transitivity_undirected()\n\n"
      "Calculates the transitivity (clustering coefficient) of the graph.\n\n"
      "@return: the transitivity\n"
  },
  
  // interface to igraph_transitivity_local_undirected
  {"transitivity_local_undirected", (PyCFunction)igraphmodule_Graph_transitivity_local_undirected,
      METH_VARARGS | METH_KEYWORDS,
   "transitivity_local_undirected(vertices=None)\n\n"
   "Calculates the local transitivity of given vertices in the graph.\n\n"
   "@param vertices: a list containing the vertex IDs which should be\n"
   "  included in the result. C{None} means all of the vertices.\n"
   "@return: the transitivities for the given vertices in a list\n"
  },
  
  //////////////////////
  // LAYOUT FUNCTIONS //
  //////////////////////
  
  // interface to igraph_layout_circle
  {"layout_circle", (PyCFunction)igraphmodule_Graph_layout_circle,
      METH_VARARGS | METH_KEYWORDS,
      "layout_circle()\n\n"
      "Places the vertices of the graph uniformly on a circle.\n\n"
      "@return: the calculated coordinate pairs in a list."
  },
  
  // interface to igraph_layout_sphere
  {"layout_sphere", (PyCFunction)igraphmodule_Graph_layout_sphere,
      METH_VARARGS | METH_KEYWORDS,
      "layout_sphere()\n\n"
      "Places the vertices of the graph uniformly on a sphere.\n\n"
      "@return: the calculated coordinate triplets in a list."
  },
  
  // interface to igraph_layout_kamada_kawai
  {"layout_kamada_kawai", (PyCFunction)igraphmodule_Graph_layout_kamada_kawai,
      METH_VARARGS | METH_KEYWORDS,
      "layout_kamada_kawai(maxiter=1000, sigma=None, initemp=10, coolexp=0.99, kkconst=None)\n\n"
      "Places the vertices on a plane according to the Kamada-Kawai algorithm.\n\n"
      "This is a force directed layout, see Kamada, T. and Kawai, S.:\n"
      "An Algorithm for Drawing General Undirected Graphs.\n"
      "Information Processing Letters, 31/1, 7--15, 1989.\n\n"
      "@param maxiter: the number of iterations to perform.\n"
      "@param sigma: the standard base deviation of the position\n"
      "  change proposals. C{None} means the number of vertices * 0.25\n"
      "@param initemp: initial temperature of the simulated annealing.\n"
      "@param coolexp: cooling exponent of the simulated annealing.\n"
      "@param kkconst: the Kamada-Kawai vertex attraction constant.\n"
      "  C{None} means the square of the number of vertices.\n"
      "@return: the calculated coordinate pairs in a list."
  },
  
  // interface to igraph_layout_kamada_kawai_3d
  {"layout_kamada_kawai_3d", (PyCFunction)igraphmodule_Graph_layout_kamada_kawai_3d,
      METH_VARARGS | METH_KEYWORDS,
      "layout_kamada_kawai_3d(maxiter=1000, sigma=None, initemp=10, coolexp=0.99, kkconst=None)\n\n"
      "Places the vertices in the 3D space according to the Kamada-Kawai algorithm.\n\n"
      "This is a force directed layout, see Kamada, T. and Kawai, S.:\n"
      "An Algorithm for Drawing General Undirected Graphs.\n"
      "Information Processing Letters, 31/1, 7--15, 1989.\n\n"
      "@param maxiter: the number of iterations to perform.\n"
      "@param sigma: the standard base deviation of the position\n"
      "  change proposals. C{None} means the number of vertices * 0.25\n"
      "@param initemp: initial temperature of the simulated annealing.\n"
      "@param coolexp: cooling exponent of the simulated annealing.\n"
      "@param kkconst: the Kamada-Kawai vertex attraction constant.\n"
      "  C{None} means the square of the number of vertices.\n"
      "@return: the calculated coordinate triplets in a list."
  },
  
  // interface to igraph_layout_fruchterman_reingold
  {"layout_fruchterman_reingold", (PyCFunction)igraphmodule_Graph_layout_fruchterman_reingold,
      METH_VARARGS | METH_KEYWORDS,
      "layout_fruchterman_reingold(maxiter=500, maxdelta=None, area=None, coolexp=0.99, repulserad=maxiter*maxdelta)\n\n"
      "Places the vertices on a 2D plane according to the Fruchterman-Reingold algorithm.\n\n"
      "This is a force directed layout, see Fruchterman, T. M. J. and Reingold, E. M.:\n"
      "Graph Drawing by Force-directed Placement.\n"
      "Software -- Practice and Experience, 21/11, 1129--1164, 1991\n\n"
      "@param maxiter: the number of iterations to perform.\n"
      "@param maxdelta: the maximum distance to move a vertex in\n"
      "  an iteration. C{None} means the number of vertices.\n"
      "@param area: the area of the square on which the vertices\n"
      "  will be placed. C{None} means the square of M{maxdelta}.\n"
      "@param coolexp: the cooling exponent of the simulated annealing.\n"
      "@param repulserad: determines the radius at which vertex-vertex\n"
      "  repulsion cancels out attraction of adjacent vertices.\n"
      "  C{None} means M{maxiter*maxdelta}.\n"
      "@return: the calculated coordinate pairs in a list."
  },
  
  // interface to igraph_layout_fruchterman_reingold_3d
  {"layout_fruchterman_reingold_3d", (PyCFunction)igraphmodule_Graph_layout_fruchterman_reingold_3d,
      METH_VARARGS | METH_KEYWORDS,
      "layout_fruchterman_reingold_3d(maxiter=500, maxdelta=None, area=None, coolexp=0.99, repulserad=maxiter*maxdelta)\n\n"
      "Places the vertices in the 3D space according to the Fruchterman-Reingold grid algorithm.\n\n"
      "This is a force directed layout, see Fruchterman, T. M. J. and Reingold, E. M.:\n"
      "Graph Drawing by Force-directed Placement.\n"
      "Software -- Practice and Experience, 21/11, 1129--1164, 1991\n\n"
      "@param maxiter: the number of iterations to perform.\n"
      "@param maxdelta: the maximum distance to move a vertex in\n"
      "  an iteration. C{None} means the number of vertices.\n"
      "@param area: the area of the square on which the vertices\n"
      "  will be placed. C{None} means the square of M{maxdelta}.\n"
      "@param coolexp: the cooling exponent of the simulated annealing.\n"
      "@param repulserad: determines the radius at which vertex-vertex\n"
      "  repulsion cancels out attraction of adjacent vertices.\n"
      "  C{None} means M{maxiter*maxdelta}.\n"
      "@return: the calculated coordinate triplets in a list."
  },
  
  // interface to igraph_layout_grid_fruchterman_reingold
  {"layout_grid_fruchterman_reingold", (PyCFunction)igraphmodule_Graph_layout_grid_fruchterman_reingold,
      METH_VARARGS | METH_KEYWORDS,
      "layout_grid_fruchterman_reingold(maxiter=500, maxdelta=None, area=None, coolexp=0.99, repulserad=maxiter*maxdelta, cellsize=1.0)\n\n"
      "Places the vertices on a 2D plane according to the Fruchterman-Reingold grid algorithm.\n\n"
      "This is a modified version of a force directed layout, see\n"
      "Fruchterman, T. M. J. and Reingold, E. M.:\n"
      "Graph Drawing by Force-directed Placement.\n"
      "Software -- Practice and Experience, 21/11, 1129--1164, 1991.\n"
      "The algorithm partitions the 2D space to a grid and vertex\n"
      "repulsion is then calculated only for vertices nearby.\n\n"
      "@param maxiter: the number of iterations to perform.\n"
      "@param maxdelta: the maximum distance to move a vertex in\n"
      "  an iteration. C{None} means the number of vertices.\n"
      "@param area: the area of the square on which the vertices\n"
      "  will be placed. C{None} means the square of M{maxdelta}.\n"
      "@param coolexp: the cooling exponent of the simulated annealing.\n"
      "@param repulserad: determines the radius at which vertex-vertex\n"
      "  repulsion cancels out attraction of adjacent vertices.\n"
      "  C{None} means M{maxiter*maxdelta}.\n"
      "@param cellsize: the size of the grid cells.\n"
      "@return: the calculated coordinate pairs in a list."
  },
  
  // interface to igraph_layout_lgl
  {"layout_lgl", (PyCFunction)igraphmodule_Graph_layout_lgl,
      METH_VARARGS | METH_KEYWORDS,
      "layout_lgl(maxiter=500, maxdelta=None, area=None, coolexp=0.99, repulserad=maxiter*maxdelta, cellsize=1.0, proot=None)\n\n"
      "Places the vertices on a 2D plane according to the Large Graph Layout.\n\n"
      "@param maxiter: the number of iterations to perform.\n"
      "@param maxdelta: the maximum distance to move a vertex in\n"
      "  an iteration. C{None} means the number of vertices.\n"
      "@param area: the area of the square on which the vertices\n"
      "  will be placed. C{None} means the square of M{maxdelta}.\n"
      "@param coolexp: the cooling exponent of the simulated annealing.\n"
      "@param repulserad: determines the radius at which vertex-vertex\n"
      "  repulsion cancels out attraction of adjacent vertices.\n"
      "  C{None} means M{maxiter*maxdelta}.\n"
      "@param cellsize: the size of the grid cells.\n"
      "@param proot: the root vertex, this is placed first, its neighbors\n"
      "  in the first iteration, second neighbors in the second,\n"
      "  etc. C{None} means a random vertex.\n"
      "@return: the calculated coordinate pairs in a list."
  },
  
  // interface to igraph_layout_reingold_tilford
  {"layout_reingold_tilford", (PyCFunction)igraphmodule_Graph_layout_reingold_tilford,
      METH_VARARGS | METH_KEYWORDS,
      "layout_reingold_tilford(root)\n"
      "Places the vertices on a 2D plane according to the Reingold-Tilford\n"
      "layout algorithm. See the following reference for details:\n"
      "EM Reingold, JS Tilford: Tidier Drawings of Trees.\n"
      "IEEE Transactions on Software Engineering 7:22, 223-228, 1981.\n\n"
      "@param root: the root of the tree.\n"
      "@return: the calculated coordinate pairs in a list."
  },
  
  // interface to igraph_layout_random
  {"layout_random", (PyCFunction)igraphmodule_Graph_layout_random,
      METH_VARARGS | METH_KEYWORDS,
      "layout_random()\n"
      "Places the vertices of the graph randomly in a 2D space.\n\n"
      "@return: the \"calculated\" coordinate pairs in a list."
  },
   
  // interface to igraph_layout_random_3d
  {"layout_random_3d", (PyCFunction)igraphmodule_Graph_layout_random_3d,
      METH_VARARGS | METH_KEYWORDS,
      "layout_random_3d()\n"
      "Places the vertices of the graph randomly in a 3D space.\n\n"
      "@return: the \"calculated\" coordinate triplets in a list."
  },

  ////////////////////////////
  // VISITOR-LIKE FUNCTIONS //
  ////////////////////////////
  {"bfs", (PyCFunction)igraphmodule_Graph_bfs,
      METH_VARARGS | METH_KEYWORDS,
      "bfs(vid, mode=OUT)\n\n"
      "Conducts a breadth first search (BFS) on the graph.\n\n"
      "@param vid: the root vertex ID\n"
      "@param mode: either C{IN} or C{OUT} or C{ALL}, ignored\n"
      "  for undirected graphs.\n"
      "@return: a tuple with the following items:\n"
      "   - The vertex IDs visited (in order)\n"
      "   - The start indices of the layers in the vertex list\n"
      "   - The parent of every vertex in the BFS\n"
  },
  {"bfsiter", (PyCFunction)igraphmodule_Graph_bfsiter,
      METH_VARARGS | METH_KEYWORDS,
      "bfsiter(vid, mode=OUT, advanced=False)\n\n"
      "Constructs a breadth first search (BFS) iterator of the graph.\n\n"
      "@param vid: the root vertex ID\n"
      "@param mode: either C{IN} or C{OUT} or C{ALL}.\n"
      "@param advanced: if C{False}, the iterator returns the next\n"
      "  vertex in BFS order in every step. If C{True}, the iterator\n"
      "  returns the distance of the vertex from the root and the\n"
      "  parent of the vertex in the BFS tree as well.\n"
      "@return: the BFS iterator as an L{igraph.BFSIter} object.\n"
  },
  
  /////////////////
  // CONVERSIONS //
  /////////////////
  
  // interface to igraph_get_adjacency
  {"get_adjacency", (PyCFunction)igraphmodule_Graph_get_adjacency,
      METH_VARARGS | METH_KEYWORDS,
      "get_adjacency(type=GET_ADJACENCY_BOTH)\n\n"
      "Returns the adjacency matrix of a graph.\n\n"
      "@param type: either C{GET_ADJACENCY_LOWER} (uses the\n"
      "  lower triangle of the matrix) or C{GET_ADJACENCY_UPPER}\n"
      "  (uses the upper triangle) or C{GET_ADJACENCY_BOTH}\n"
      "  (uses both parts). Ignored for directed graphs.\n"
      "@return: the adjacency matrix.\n"
  },
  
  // interface to igraph_get_edgelist
  {"get_edgelist", (PyCFunction)igraphmodule_Graph_get_edgelist,
      METH_NOARGS,
      "get_edgelist()\n\n"
      "Returns the edge list of a graph."
  },
  
  // interface to igraph_to_directed
  {"to_directed", (PyCFunction)igraphmodule_Graph_to_directed,
   METH_VARARGS | METH_KEYWORDS,
   "to_directed(mutual=True)\n\n"
   "Converts an undirected graph to directed.\n\n"
   "@param mutual: C{True} if mutual directed edges should be\n"
   "  created for every undirected edge. If C{False}, a directed\n"
   "  edge with arbitrary direction is created.\n"
  },

  // interface to igraph_to_undirected
  {"to_undirected", (PyCFunction)igraphmodule_Graph_to_undirected,
   METH_VARARGS | METH_KEYWORDS,
   "to_undirected(collapse=True)\n\n"
   "Converts a directed graph to undirected.\n\n"
   "@param collapse: C{True} if only a single edge should be\n"
   "  created from multiple directed edges going between the\n"
   "  same vertex pair. If C{False}, the edge count is kept constant.\n"
  },

  ///////////////////////////////
  // LOADING AND SAVING GRAPHS //
  ///////////////////////////////
  
  // interface to igraph_read_graph_edgelist
  {"Read_Edgelist", (PyCFunction)igraphmodule_Graph_Read_Edgelist,
      METH_VARARGS | METH_KEYWORDS | METH_CLASS,
      "Read_Edgelist(f, directed=True)\n\n"
      "Reads an edge list from a file and creates a graph based on it.\n\n"
      "Please note that the vertex indices are zero-based.\n\n"
      "@param f: the name of the file\n"
      "@param directed: whether the generated graph should be directed.\n"
  },
  // interface to igraph_read_graph_ncol
  {"Read_Ncol", (PyCFunction)igraphmodule_Graph_Read_Ncol,
      METH_VARARGS | METH_KEYWORDS | METH_CLASS,
      "Read_Ncol(f, names=True, weights=True)\n\n"
      "Reads an .ncol file used by LGL.\n\n"
      "It is also useful for creating graphs from \"named\" (and\n"
      "optionally weighted) edge lists.\n\n"
      "This format is used by the Large Graph Layout program. See the\n"
      "U{documentation of LGL <http://bioinformatics.icmb.utexas.edu/bgl/>}\n"
      "regarding the exact format description.\n\n"
      "LGL originally cannot deal with graphs containing multiple or loop\n"
      "edges, but this condition is not checked here, as igraph is happy\n"
      "with these.\n\n"
      "@param f: the name of the file\n"
      "@param names: If C{True}, the vertex names are added as a\n"
      "  vertex attribute called 'name'.\n"
      "@param weights: If True, the edge weights are added as an\n"
      "  edge attribute called 'weight'.\n"
  },
  // interface to igraph_read_graph_lgl
  {"Read_Lgl", (PyCFunction)igraphmodule_Graph_Read_Lgl,
      METH_VARARGS | METH_KEYWORDS | METH_CLASS,
      "Read_Lgl(f, names=True, weights=True)\n\n"
      "Reads an .lgl file used by LGL.\n\n"
      "It is also useful for creating graphs from \"named\" (and\n"
      "optionally weighted) edge lists.\n\n"
      "This format is used by the Large Graph Layout program. See the\n"
      "U{documentation of LGL <http://bioinformatics.icmb.utexas.edu/bgl/>}\n"
      "regarding the exact format description.\n\n"
      "LGL originally cannot deal with graphs containing multiple or loop\n"
      "edges, but this condition is not checked here, as igraph is happy\n"
      "with these.\n\n"
      "@param f: the name of the file\n"
      "@param names: If C{True}, the vertex names are added as a\n"
      "  vertex attribute called 'name'.\n"
      "@param weights: If True, the edge weights are added as an\n"
      "  edge attribute called 'weight'.\n"
  },
  // interface to igraph_read_graph_pajek
  {"Read_Pajek", (PyCFunction)igraphmodule_Graph_Read_Pajek,
      METH_VARARGS | METH_KEYWORDS | METH_CLASS,
      "Read_Pajek(f)\n\n"
      "Reads a Pajek format file and creates a graph based on it.\n\n"
      "@param f: the name of the file\n"
  },
  // interface to igraph_read_graph_graphml
  {"Read_GraphML", (PyCFunction)igraphmodule_Graph_Read_GraphML,
      METH_VARARGS | METH_KEYWORDS | METH_CLASS,
      "Read_GraphML(f, directed=True, index=0)\n\n"
      "Reads a GraphML format file and creates a graph based on it.\n\n"
      "@param f: the name of the file\n"
      "@param index: if the GraphML file contains multiple graphs,\n"
      "  specifies the one that should be loaded. Graph indices\n"
      "  start from zero, so if you want to load the first graph,\n"
      "  specify 0 here.\n"
  },
  // interface to igraph_write_graph_edgelist
  {"write_edgelist", (PyCFunction)igraphmodule_Graph_write_edgelist,
      METH_VARARGS | METH_KEYWORDS,
      "write_edgelist(f)\n\n"
      "Writes the edge list of a graph to a file.\n\n"
      "Directed edges are written in (from, to) order.\n\n"
      "@param f: the name of the file to be written\n"
  },
  // interface to igraph_write_graph_ncol
  {"write_ncol", (PyCFunction)igraphmodule_Graph_write_ncol,
      METH_VARARGS | METH_KEYWORDS,
      "write_ncol(f, names=\"name\", weights=\"weights\")\n\n"
      "Writes the edge list of a graph to a file in .ncol format.\n\n"
      "Note that multiple edges and/or loops break the LGL software,\n"
      "but igraph does not check for this condition. Unless you know\n"
      "that the graph does not have multiple edges and/or loops, it\n"
      "is wise to call L{simplify()} before saving.\n\n"
      "@param f: the name of the file to be written\n"
      "@param names: the name of the vertex attribute containing the name\n"
      "  of the vertices. If you don't want to store vertex names,\n"
      "  supply C{None} here.\n"
      "@param weights: the name of the edge attribute containing the weight\n"
      "  of the vertices. If you don't want to store weights,\n"
      "  supply C{None} here.\n"
  },
  // interface to igraph_write_graph_lgl
  {"write_lgl", (PyCFunction)igraphmodule_Graph_write_lgl,
      METH_VARARGS | METH_KEYWORDS,
      "write_lgl(f, names=\"name\", weights=\"weights\", isolates=True)\n\n"
      "Writes the edge list of a graph to a file in .lgl format.\n\n"
      "Note that multiple edges and/or loops break the LGL software,\n"
      "but igraph does not check for this condition. Unless you know\n"
      "that the graph does not have multiple edges and/or loops, it\n"
      "is wise to call L{simplify()} before saving.\n\n"
      "@param f: the name of the file to be written\n"
      "@param names: the name of the vertex attribute containing the name\n"
      "  of the vertices. If you don't want to store vertex names,\n"
      "  supply C{None} here.\n"
      "@param weights: the name of the edge attribute containing the weight\n"
      "  of the vertices. If you don't want to store weights,\n"
      "  supply C{None} here.\n"
      "@param isolates: whether to include isolated vertices in the output.\n"
  },
  // interface to igraph_write_graph_edgelist
  {"write_graphml", (PyCFunction)igraphmodule_Graph_write_graphml,
      METH_VARARGS | METH_KEYWORDS,
      "write_graphml(f)\n\n"
      "Writes the graph to a GraphML file.\n\n"
      "@param f: the name of the file to be written\n"
  },

  /////////////////
  // ISOMORPHISM //
  /////////////////
  {"isoclass", (PyCFunction)igraphmodule_Graph_isoclass,
      METH_VARARGS | METH_KEYWORDS,
      "isoclass(vertices)\n\n"
      "Returns the isomorphy class of the graph or its subgraph.\n\n"
      "Isomorphy class calculations are implemented only for graphs with\n"
      "3 or 4 nodes.\n\n"
      "@param vertices: a list of vertices if we want to calculate the\n"
      "  isomorphy class for only a subset of vertices. C{None} means to\n"
      "  use the full graph.\n"
      "@return: the isomorphy class of the (sub)graph\n\n"
  },
  {"isomorphic", (PyCFunction)igraphmodule_Graph_isomorphic,
      METH_VARARGS | METH_KEYWORDS,
      "isomorphic(other)\n\n"
      "Checks whether the graph is isomorphic with another graph.\n\n"
      "Works only for graphs with 3 or 4 vertices.\n\n"
      "@param other: the other graph with which we want to compare the graph.\n"
      "@return: C{True} if the graphs are isomorphic, C{False} if not.\n"
  },
  
  ////////////////////////
  // ATTRIBUTE HANDLING //
  ////////////////////////
  {"attributes", (PyCFunction)igraphmodule_Graph_attributes,
      METH_NOARGS,
      "attributes()\n\n"
      "@return: the attribute name list of the graph\n"
  },
  {"vertex_attributes", (PyCFunction)igraphmodule_Graph_vertex_attributes,
      METH_NOARGS,
      "vertex_attributes()\n\n"
      "@return: the attribute name list of the graph's vertices\n"
  },
  {"edge_attributes", (PyCFunction)igraphmodule_Graph_edge_attributes,
      METH_NOARGS,
      "edge_attributes()\n\n"
      "@return: the attribute name list of the graph's edges\n"
  },

  ///////////////
  // OPERATORS //
  ///////////////
  {"complementer", (PyCFunction)igraphmodule_Graph_complementer,
      METH_VARARGS,
      "complementer(loops=False)\n\n"
      "Returns the complementer of the graph\n\n"
      "@param loops: whether to include loop edges in the complementer.\n"
      "@return: the complementer of the graph\n"
  },
  {"compose", (PyCFunction)igraphmodule_Graph_compose,
      METH_O, "compose(other)\n\nReturns the composition of two graphs."
  },
  {"difference", (PyCFunction)igraphmodule_Graph_difference,
      METH_O, "difference(other)\n\nSubtracts the given graph from the original"
  },
  {"disjoint_union", (PyCFunction)igraphmodule_Graph_disjoint_union,
      METH_O,
      "disjoint_union(graphs)\n\n"
      "Creates the disjoint union of two (or more) graphs.\n\n"
      "@param graphs: the list of graphs to be united with the current one.\n"
  },
  {"intersection", (PyCFunction)igraphmodule_Graph_intersection,
      METH_O,
      "intersection(graphs)\n\n"
      "Creates the intersection of two (or more) graphs.\n\n"
      "@param graphs: the list of graphs to be intersected with\n"
      "  the current one.\n"
  },
  {"union", (PyCFunction)igraphmodule_Graph_union,
      METH_O,
      "union(graphs)\n\n"
      "Creates the union of two (or more) graphs.\n\n"
      "@param graphs: the list of graphs to be intersected with\n"
      "  the current one.\n"
  },

  /**************************/
  /* FLOW RELATED FUNCTIONS */  
  /**************************/
  {"maxflow_value", (PyCFunction)igraphmodule_Graph_maxflow_value,
   METH_VARARGS | METH_KEYWORDS,
   "maxflow_value(source, target, capacity=None)\n\n"
   "Returns the maximum flow between the source and target vertices.\n\n"
   "@param source: the source vertex ID\n"
   "@param target: the target vertex ID\n"
   "@param capacity: the capacity of the edges. It must be a list or a valid\n"
   "  attribute name or C{None}. In the latter case, every edge will have the\n"
   "  same capacity.\n"
   "@return: the value of the maximum flow between the given vertices\n"
  },

  {"mincut_value", (PyCFunction)igraphmodule_Graph_mincut_value,
   METH_VARARGS | METH_KEYWORDS,
   "mincut_value(source=-1, target=-1, capacity=None)\n\n"
   "Returns the minimum cut between the source and target vertices.\n\n"
   "@param source: the source vertex ID. If negative, the calculation is\n"
   "  done for every vertex except the target and the minimum is returned.\n"
   "@param target: the target vertex ID. If negative, the calculation is\n"
   "  done for every vertex except the source and the minimum is returned.\n"
   "@param capacity: the capacity of the edges. It must be a list or a valid\n"
   "  attribute name or C{None}. In the latter case, every edge will have the\n"
   "  same capacity.\n"
   "@return: the value of the minimum cut between the given vertices\n"
  },

  ////////////////////////////////////
  // INTERNAL/DEVELOPMENT FUNCTIONS //
  ////////////////////////////////////
  {"__graph_as_cobject", (PyCFunction)igraphmodule_Graph___graph_as_cobject__,
      METH_VARARGS | METH_KEYWORDS,
      "__graph_as_cobject()\n\n"
      "Returns the igraph graph encapsulated by the Python object as\n"
      "a PyCObject\n\n."
      "A PyObject is barely a regular C pointer. This function\n"
      "should not be used directly by igraph users, it is useful only\n"
      "in the case when the underlying igraph object must be passed to\n"
      "another C code through Python.\n\n"
  },
  {"__register_destructor", (PyCFunction)igraphmodule_Graph___register_destructor__,
      METH_VARARGS | METH_KEYWORDS,
      "__register_destructor(destructor)\n\n"
      "Registers a destructor to be called when the object is freed by\n"
      "Python. This function should not be used directly by igraph users."
  },
  {NULL}

}
;

#endif
