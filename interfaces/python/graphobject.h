#ifndef PYTHON_GRAPHOBJECT_H
#define PYTHON_GRAPHOBJECT_H

#include <Python.h>
#include "igraph.h"
#include "types.h"
#include "structmember.h"
#include "common.h"

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

PyObject* igraphmodule_Graph_Atlas(PyTypeObject *type, PyObject *args);
PyObject* igraphmodule_Graph_Barabasi(PyTypeObject *type, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_Erdos_Renyi(PyTypeObject *type, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_Full(PyTypeObject *type, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_Growing_Random(PyTypeObject *type, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_Star(PyTypeObject *type, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_Ring(PyTypeObject *type, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_Tree(PyTypeObject *type, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_Degree_Sequence(PyTypeObject *type, PyObject *args, PyObject *kwds);

PyObject* igraphmodule_Graph_diameter(igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_is_connected(igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_are_connected(igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_average_path_length(igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_betweenness(igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_bibcoupling(igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_closeness(igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_clusters(igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_cocitation(igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_decompose(igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_edge_betweenness(igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_get_shortest_paths(igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_get_all_shortest_paths(igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_pagerank(igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_rewire(igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_shortest_paths(igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_spanning_tree(igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_simplify(igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_subcomponent(igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_subgraph(igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds);

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

PyObject* igraphmodule_Graph_get_adjacency(igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_get_edgelist(igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds);

PyObject* igraphmodule_Graph_Read_Edgelist(PyTypeObject *type, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_Read_Ncol(PyTypeObject *type, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_Read_Lgl(PyTypeObject *type, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_Read_Pajek(PyTypeObject *type, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_Read_GraphML(PyTypeObject *type, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_write_edgelist(igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_write_ncol(igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_write_lgl(igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph_write_graphml(igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds);

int igraphmodule_Graph_attribute_count(igraphmodule_GraphObject* self);
PyObject* igraphmodule_Graph_get_attribute(igraphmodule_GraphObject* self, PyObject* s);
int igraphmodule_Graph_set_attribute(igraphmodule_GraphObject* self, PyObject* k, PyObject* v);
PyObject* igraphmodule_Graph_attributes(igraphmodule_GraphObject* self, PyObject* args, PyObject* kwds);
PyObject* igraphmodule_Graph_vertex_attributes(igraphmodule_GraphObject* self, PyObject* args, PyObject* kwds);
PyObject* igraphmodule_Graph_edge_attributes(igraphmodule_GraphObject* self, PyObject* args, PyObject* kwds);

PyObject* igraphmodule_Graph_get_vertices(igraphmodule_GraphObject* self, void* closure);
PyObject* igraphmodule_Graph_get_edges(igraphmodule_GraphObject* self, void* closure);

PyObject* igraphmodule_Graph___graph_as_cobject__(igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds);
PyObject* igraphmodule_Graph___register_destructor__(igraphmodule_GraphObject *self, PyObject *args, PyObject *kwds);

PyTypeObject igraphmodule_GraphType;

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
      "Returns the number of vertices in the graph"
  },
  // interface to igraph_ecount
  {"ecount", (PyCFunction)igraphmodule_Graph_ecount,
      METH_NOARGS,
      "Returns the number of edges in the graph"
  },
  // interface to igraph_is_directed
  {"is_directed", (PyCFunction)igraphmodule_Graph_is_directed,
      METH_NOARGS,
      "Checks whether the graph is directed"
  },
  // interface to igraph_add_vertices
  {"add_vertices", (PyCFunction)igraphmodule_Graph_add_vertices,
      METH_VARARGS,
      "Adds vertices to the graph. The only parameter is the number of "
      "vertices to be added"
  },
  // interface to igraph_delete_vertices
  {"delete_vertices", (PyCFunction)igraphmodule_Graph_delete_vertices,
      METH_VARARGS,
      "Deletes vertices and all its edges from the graph. The only "
      "parameter is a list of the vertices to be added. It is allowed "
      "to provide a single integer instead of a list consisting of only "
      "one integer."
  },
  // interface to igraph_add_edges
  {"add_edges", (PyCFunction)igraphmodule_Graph_add_edges,
      METH_VARARGS,
      "Adds edges to the graph. The only parameter is a list of "
      "edges to be added. Every edge is represented with a tuple, "
      "containing the vertex IDs of the two endpoints. Vertices are "
      "enumerated from zero. It is allowed to provide a single pair "
      "instead of a list consisting of only one pair."
  },
  // interface to igraph_delete_edges
  {"delete_edges", (PyCFunction)igraphmodule_Graph_delete_edges,
      METH_VARARGS,
      "Removes edges from the graph. The only parameter is a list of "
      "edges to be removed. Every edge is represented with a tuple, "
      "containing the vertex IDs of the two endpoints. Vertices are "
      "enumerated from zero. It is allowed to provide a single pair "
      "instead of a list consisting of only one pair. Nonexistent "
      "edges will be silently ignored. All vertices will be kept, even "
      "if they lose all their edges."
  },
  // interface to igraph_degree
  {"degree", (PyCFunction)igraphmodule_Graph_degree,
      METH_VARARGS | METH_KEYWORDS,
      "Returns some vertex degrees from the graph.\n"
      "This method accepts a single vertex ID or a list of vertex IDs as a "
      "parameter, and returns the degree of the given vertices (in the form of "
      "a single integer or a list, depending on the input parameter). A "
      "second and a third argument may be passed as well, the second one "
      "meaning the type of degree to be returned (OUT for out-degrees, "
      "IN for in-degrees or ALL for the sum of them) and the third one "
      "meaning whether self-loops should be counted. The default for them is "
      "ALL and False. The type of degree is ignored for undirected graphs."
  },
  // interfaces to igraph_neighbors
  {"neighbors", (PyCFunction)igraphmodule_Graph_neighbors,
      METH_VARARGS | METH_KEYWORDS,
      "Returns adjacent vertices to a given vertex.\n"
      "This method accepts a single vertex ID as an argument, "
      "and returns the adjacent vertices of that vertex. An optional "
      "second argument allows the user to limit the result to only "
      "predecessors (IN), only successors (OUT) or both of them (ALL). "
      "The default behaviour is the latter one. "
      "This argument is ignored for undirected graphs."
  },
  {"successors", (PyCFunction)igraphmodule_Graph_successors,
      METH_VARARGS | METH_KEYWORDS,
      "Returns the successors of a given vertex.\n"
      "Equivalent to calling the neighbors method with type=OUT."
  },
  {"predecessors", (PyCFunction)igraphmodule_Graph_predecessors,
      METH_VARARGS | METH_KEYWORDS,
      "Returns the predecessors of a given vertex.\n"
      "Equivalent to calling the neighbors method with type=IN."
  },
  
  //////////////////////
  // GRAPH GENERATORS //
  //////////////////////
  
  // interface to igraph_atlas
  {"Atlas", (PyCFunction)igraphmodule_Graph_Atlas,
      METH_CLASS | METH_KEYWORDS,
      "Generates a graph from the Graph Atlas.\n"
      "The only argument denotes the index of the graph to be generated.\n"
      "Indices start from zero, graphs are listed:\n"
      "1. in increasing order of number of nodes;\n"
      "2. for a fixed number of nodes, in increasing order of the\n"
      "   number of edges;\n"
      "3. for fixed numbers of nodes and edges, in increasing order\n"
      "   of the degree sequence, for example 111223 < 112222;\n"
      "4. for fixed degree sequence, in increasing number of automorphisms.\n"
  },
	
  // interface to igraph_barabasi_game
  {"Barabasi", (PyCFunction)igraphmodule_Graph_Barabasi,
      METH_VARARGS | METH_CLASS | METH_KEYWORDS,
      "Generates a graph based on the Barabási-Albert model.\n"
      "The first two arguments are mandatory: the first one is the "
      "number of vertices, the second one is either the number of "
      "outgoing edges generated for each vertex or a list containing the "
      "number of outgoing edges for each vertex explicitly. The third "
      "argument is True if the out-degree of a given vertex should also "
      "increase its citation probability (as well as its in-degree), but "
      "it defaults to False. The fourth argument is True if the generated "
      "graph should be directed (default: False).\n\n"
      "Keywords for the arguments: n, m, outpref, directed"
  },
  
  // interface to igraph_erdos_renyi_game
  {"Erdos_Renyi", (PyCFunction)igraphmodule_Graph_Erdos_Renyi,
      METH_VARARGS | METH_CLASS | METH_KEYWORDS,
      "Generates a graph based on the Erdõs-Rényi model.\n"
      "There are a total of five possible arguments, two of them are "
      "mutually exclusive. The first argument (keyword: n) is the number "
      "of vertices. The second and the third (keywords: p and m) define "
      "the density of the graph: if p is missing, there will be m edges; "
      "if m is missing, every edge will be present with a probability of "
      "p. These two arguments influence the same graph property (the "
      "number of edges) in two different ways, so either p or m must be "
      "present (but not both of them). The remaining two arguments are "
      "optional. The fourth argument (keyword: directed) is True if the "
      "generated graph should be directed (default: False), the fifth "
      "(keyword: loops) is True if self-loops are allowed (default: False)."
  },
  
  // interface to igraph_full_game
  {"Full", (PyCFunction)igraphmodule_Graph_Full,
      METH_VARARGS | METH_CLASS | METH_KEYWORDS,
      "Generates a full graph (directed or undirected, with or without loops).\n"
      "The only mandatory argument (keyword: n) is the number "
      "of vertices. The remaining two arguments are optional. "
      "The second argument (keyword: directed) is True if the "
      "generated graph should be directed (default: False), the third "
      "(keyword: loops) is True if self-loops are allowed (default: False)."
  },
  
  // interface to igraph_growing_random_game
  {"Growing_Random", (PyCFunction)igraphmodule_Graph_Growing_Random,
      METH_VARARGS | METH_CLASS | METH_KEYWORDS,
      "Generates a growing random graph.\n\n"
      "Keyword arguments:\n"
      "n -- The number of vertices in the graph\n"
      "m -- The number of edges to add in each step (after adding a new vertex)\n"
      "directed -- whether the graph should be directed.\n"
      "            Optional, defaults to False.\n"
      "citation -- whether the new edges should originate from the most\n"
      "            recently added vertex.\n"
      "            Optional, defaults to False."
  },
  
  // interface to igraph_star
  {"Star", (PyCFunction)igraphmodule_Graph_Star,
      METH_VARARGS | METH_CLASS | METH_KEYWORDS,
      "Generates a star graph.\n\n"
      "Keyword arguments:\n"
      "n -- The number of vertices in the graph\n"
      "mode -- Gives the type of the star graph to create. Should be\n"
      "        one of the constants STAR_OUT, STAR_IN and STAR_UNDIRECTED.\n"
      "        Optional, defaults to STAR_UNDIRECTED.\n"
      "center -- Vertex ID for the central vertex in the star.\n"
      "          Optional, defaults to zero.\n"
  },
  
  // interface to igraph_lattice
  {"Lattice", (PyCFunction)igraphmodule_unimplemented,
      METH_VARARGS | METH_CLASS | METH_KEYWORDS,
      "Generates a lattice. This function is yet unimplemented.\n\n"
      "Throws a NotImplementedError."
  },
  
  // interface to igraph_ring
  {"Ring", (PyCFunction)igraphmodule_Graph_Ring,
      METH_VARARGS | METH_CLASS | METH_KEYWORDS,
      "Generates a ring graph.\n\n"
      "Keyword arguments:\n"
      "n -- The number of vertices in the ring\n"
      "directed -- whether to create a directed ring.\n"
      "            Optional, defaults to False.\n"
      "mutual -- whether to create mutual edges in a directed ring.\n"
      "          Optional, defaults to False.\n"
      "          Ignored for undirected graphs.\n"
      "circular -- whether to create a closed ring.\n"
      "            Optional, defaults to True."
  },
  
  // interface to igraph_tree
  {"Tree", (PyCFunction)igraphmodule_Graph_Tree,
      METH_VARARGS | METH_CLASS | METH_KEYWORDS,
      "Generates a tree in which almost all vertices have the same number of children.\n\n"
      "Keyword arguments:\n"
      "n -- The number of vertices in the graph\n"
      "children -- The number of children of a vertex in the graph\n"
      "type -- determines whether the tree should be directed, and if\n"
      "        this is the case, also its orientation. Must be one of\n"
      "        TREE_IN, TREE_OUT and TREE_UNDIRECTED.\n"
      "        Optional, defaults to TREE_UNDIRECTED\n"
  },
  
  // interface to igraph_adjacency
  {"Adjacency", (PyCFunction)igraphmodule_unimplemented,
      METH_VARARGS | METH_CLASS | METH_KEYWORDS,
      "Generates a graph from an adjacency matrix. This function is yet unimplemented.\n\n"
      "Throws a NotImplementedError."
  },
  
  // interface to igraph_degree_sequence_game
  {"Degree_Sequence", (PyCFunction)igraphmodule_Graph_Degree_Sequence,
      METH_VARARGS | METH_CLASS | METH_KEYWORDS,
      "Generates a graph with a given degree sequence.\n\n"
      "Keyword arguments:\n"
      "out -- the out-degree sequence for a directed graph. If the\n"
      "       in-degree sequence is omitted, the generated graph\n"
      "       will be undirected, so this will be the in-degree\n"
      "       sequence as well\n"
      "in  -- the in-degree sequence for a directed graph. Optional,\n"
      "       if omitted, the generated graph will be undirected.\n"
  },
  
  /////////////////////////////////////
  // STRUCTURAL PROPERTIES OF GRAPHS //
  /////////////////////////////////////
  
  // interface to igraph_are_connected
  {"are_connected", (PyCFunction)igraphmodule_Graph_are_connected,
      METH_VARARGS | METH_KEYWORDS,
      "Decides whether two given vertices are directly connected.\n\n"
      "Keyword arguments:\n"
      "v1 -- the first vertex\n"
      "v2 -- the second vertex\n"
      "Returns true if there exists an edge from v1 to v2."
  },
  
  // interface to igraph_average_path_length
  {"average_path_length", (PyCFunction)igraphmodule_Graph_average_path_length,
      METH_VARARGS | METH_KEYWORDS,
      "Calculates the average path length in a graph.\n\n"
      "Keyword arguments:\n"
      "directed -- whether to consider directed paths.\n"
      "            Ignored for undirected graphs. Optional, defaults to True.\n"
      "unconn -- what to do when the graph is unconnected. If True, the\n"
      "          average of the geodesic lengths in the components is\n"
      "          calculated. Otherwise for all unconnected vertex pairs,\n"
      "          a path length equal to the number of vertices is used.\n"
  },
  
  // interface to igraph_betweenness
  {"betweenness", (PyCFunction)igraphmodule_Graph_betweenness,
      METH_VARARGS | METH_KEYWORDS,
      "Calculates the betweenness of nodes in a graph.\n\n"
      "Keyword arguments:\n"
      "vertices -- the vertices for which the betweennesses must be returned.\n"
      "            Optional, defaults to all of the vertices in the graph.\n"
      "directed -- whether to consider directed paths.\n"
      "            Ignored for undirected graphs. Optional, defaults to True.\n"
  },
  
  // interface to igraph_bibcoupling
  {"bibcoupling", (PyCFunction)igraphmodule_Graph_bibcoupling,
      METH_VARARGS | METH_KEYWORDS,
      "Calculates bibliographic coupling values for given vertices in a graph.\n\n"
      "Keyword arguments:\n"
      "vertices -- the vertices to be analysed.\n"
      "Returns bibliographic coupling values for all given vertices in a matrix."
  },
  
  // interface to igraph_closeness
  {"closeness", (PyCFunction)igraphmodule_Graph_closeness,
      METH_VARARGS | METH_KEYWORDS,
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
      "Keyword arguments:\n"
      "vertices -- the vertices for which the betweennesses must be returned.\n"
      "            Optional, defaults to all of the vertices in the graph.\n"
      "mode -- must be one of IN, OUT and ALL. IN means that the length of\n"
      "        incoming paths, OUT means that the length of the outgoing\n"
      "        paths must be calculated. ALL means that both of them must\n"
      "        be calculated. Optional, defaults to ALL.\n"
  },   	    
  
  // interface to igraph_clusters
  {"clusters", (PyCFunction)igraphmodule_Graph_clusters,
      METH_VARARGS | METH_KEYWORDS,
      "Calculates the (strong or weak) clusters for a given graph.\n\n"
      "Keyword arguments:\n"
      "mode -- must be either STRONG or WEAK, depending on the clusters\n"
      "        being sought. Optional, defaults to STRONG.\n"
      "Returns the component index for every node in the graph."
  },
  {"components", (PyCFunction)igraphmodule_Graph_clusters,
      METH_VARARGS | METH_KEYWORDS,
      "Alias for 'clusters'.\n\n"
      "See the documentation of 'clusters' for details."
  },
  {"decompose", (PyCFunction)igraphmodule_Graph_decompose,
      METH_VARARGS | METH_KEYWORDS,
      "Decomposes the graph into subgraphs.\n\n"
      "Keyword arguments:\n"
      "mode -- must be either STRONG or WEAK, depending on the clusters\n"
      "        being sought. Optional, defaults to STRONG.\n"
      "maxcompno -- maximum number of components to return. Optional,\n"
      "             defaults to all possible components.\n"
      "minelements -- minimum number of vertices in a component. Optional,\n"
      "               defaults to 1. By setting this to 2, isolated vertices\n"
      "               are not returned as separate components.\n"
      "A list of the subgraphs is returned. Every returned subgraph is a\n"
      "copy of the original.\n"
  },
  
  // interface to igraph_cocitation
  {"cocitation", (PyCFunction)igraphmodule_Graph_cocitation,
      METH_VARARGS | METH_KEYWORDS,
      "Calculates cocitation scores for given vertices in a graph.\n\n"
      "Keyword arguments:\n"
      "vertices -- the vertices to be analysed.\n"
      "Returns cocitation scores for all given vertices in a matrix."
  },
  
  // interface to igraph_diameter
  {"diameter", (PyCFunction)igraphmodule_Graph_diameter,
      METH_VARARGS | METH_KEYWORDS,
      "Calculates the diameter of the graph.\n\n"
      "Keyword arguments:\n"
      "directed -- whether to consider directed paths.\n"
      "            Ignored for undirected graphs. Optional, defaults to True.\n"
      "unconn -- if True and the graph is undirected, the longest geodesic "
      "          within a component will be returned. If False and the "
      "          graph is undirected, the result is the number of vertices."
  },
  
  // interface to igraph_edge_betweenness
  {"edge_betweenness", (PyCFunction)igraphmodule_Graph_edge_betweenness,
      METH_VARARGS | METH_KEYWORDS,
      "Calculates the edge betweennesses in a graph.\n\n"
      "Keyword arguments:\n"
      "directed -- whether to consider directed paths.\n"
      "            Ignored for undirected graphs. Optional, defaults to True.\n"
      "Returns a list with the edge betweennesses of all the edges.\n"
  },
  
  // interface to igraph_get_shortest_paths
  {"get_shortest_paths", (PyCFunction)igraphmodule_Graph_get_shortest_paths,
      METH_VARARGS | METH_KEYWORDS,
      "Calculates the shortest paths from/to a given node in a graph.\n\n"
      "Keyword arguments:\n"
      "v    -- the source/destination for the calculated paths\n"
      "mode -- the directionality of the paths. IN means to calculate\n"
      "        incoming paths, OUT means to calculate outgoing paths,\n"
      "        ALL means to calculate both ones. Ignored for undirected\n"
      "        graphs. Optional, defaults to ALL\n"
      "Returns at most one shortest path for every node in the graph in a\n"
      "list. For unconnected graphs, some of the list elements will be\n"
      "an empty list. Note that in case of mode=IN, the nodes in a path\n"
      "are returned in reversed order!"
  },
  
  // interface to igraph_get_all_shortest_paths
  {"get_all_shortest_paths", (PyCFunction)igraphmodule_Graph_get_all_shortest_paths,
      METH_VARARGS | METH_KEYWORDS,
      "Calculates all of the shortest paths from/to a given node in a graph.\n\n"
      "Keyword arguments:\n"
      "v    -- the source/destination for the calculated paths\n"
      "mode -- the directionality of the paths. IN means to calculate\n"
      "        incoming paths, OUT means to calculate outgoing paths,\n"
      "        ALL means to calculate both ones. Ignored for undirected\n"
      "        graphs. Optional, defaults to ALL\n"
      "Returns all of the shortest path from the given node to every other\n"
      "reachable node in the graph in a list. Note that in case of mode=IN,\n"
      "the nodes in a path are returned in reversed order!"
  },
  
  // interface to igraph_is_connected
  {"is_connected", (PyCFunction)igraphmodule_Graph_is_connected,
      METH_VARARGS | METH_KEYWORDS,
      "Decides whether a graph is connected.\n\n"
      "Keyword arguments:\n"
      "mode -- whether we should calculate strong or weak connectivity.\n"
      "        Ignored for undirected graphs. Optional, defaults to\n"
      "        STRONG."
  },
  
  // interface to igraph_pagerank
  {"pagerank", (PyCFunction)igraphmodule_Graph_pagerank,
      METH_VARARGS | METH_KEYWORDS,
      "Calculates the Google PageRank values of a graph.\n\n"
      "Keyword arguments:\n"
      "vertices -- the indices of the vertices being queried. Optional,"
      "            defaults to all of the vertices.\n"
      "directed -- whether to consider directed paths.\n"
      "            Ignored for undirected graphs. Optional, defaults to True.\n"
      "niter    -- the maximum number of iterations to be performed.\n"
      "            Optional, defaults to 1000.\n"
      "eps      -- the iteration stops if all of the PageRank values change\n"
      "            less than eps between two iterations. Optional, defaults\n"
      "            to 0.001\n"
      "damping  -- the damping factor. Optional, defaults to 0.85.\n"
      "            1-damping is the PageRank value for nodes with no\n"
      "            incoming links.\n"
      "Returns a list with the Google PageRank values of the specified\n"
      "vertices."
  },

  // interface to igraph_rewire
  {"rewire", (PyCFunction)igraphmodule_Graph_rewire,
      METH_VARARGS | METH_KEYWORDS,
      "Randomly rewires the graph while preserving the degree distribution.\n\n"
      "Please note that the rewiring is done \"in-place\", so the original\n"
      "graph will be modified. If you want to preserve the original graph,\n"
      "use the copy method before.\n\n"
      "Keyword arguments:\n"
      "n        -- the number of rewiring trials. Optional,\n"
      "            defaults to 1000.\n"
      "mode     -- the rewiring algorithm to use. Optional,\n"
      "            defaults to REWIRING_SIMPLE.\n"
      "Returns the modified graph.\n"
  },
  
  // interface to igraph_shortest_paths
  {"shortest_paths", (PyCFunction)igraphmodule_Graph_shortest_paths,
      METH_VARARGS | METH_KEYWORDS,
      "Calculates shortest path lengths for given nodes in a graph.\n\n"
      "Keyword arguments:\n"
      "vertices -- a list containing the vertex IDs which should be included in the result.\n"
      "mode -- the type of shortest paths to be used for the calculation in directed graphs.\n"
      "        OUT -- outgoing paths, IN -- incoming paths, ALL -- the directed graph\n"
      "        is considered as an undirected one."
  },
  
  // interface to igraph_simplify
  {"simplify", (PyCFunction)igraphmodule_Graph_simplify,
      METH_VARARGS | METH_KEYWORDS,
      "Simplifies a graph by removing self-loops and/or multiple edges.\n"
      "Keywords arguments:\n"
      "multiple -- whether to remove multiple edges. Optional, defaults\n"
      "            to True.\n"
      "loops -- whether to remove loops. Optional, defaults to True.\n"
  },
  
  // interface to igraph_minimum_spanning_tree_unweighted and
  // igraph_minimum_spanning_tree_prim
  {"spanning_tree", (PyCFunction)igraphmodule_Graph_spanning_tree,
      METH_VARARGS | METH_KEYWORDS,
      "Calculates a minimum spanning tree for a graph (weighted or unweighted)\n\n"
      "Keyword arguments:\n"
      "weights -- a vector containing weights for every edge in the graph.\n"
      "           If omitted, every edge is assumed to have an equal weight.\n"
      "Returns the spanning tree as an igraph.Graph object."
  },
  
  // interface to igraph_subcomponent
  {"subcomponent", (PyCFunction)igraphmodule_Graph_subcomponent,
      METH_VARARGS | METH_KEYWORDS,
      "Returns the indices of vertices which are in the same component as a given vertex.\n\n"
      "Keyword arguments:\n"
      "v -- the index of the vertex used as the source/destination\n"
      "mode -- if equals to IN, returns the vertex IDs from where the\n"
      "        given vertex can be reached. If equals to OUT, returns the\n"
      "        vertex IDs which are reachable from the given vertex. If\n"
      "        equals to ALL, returns all vertices within the same component\n"
      "        as the given vertex, ignoring edge directions. Note that this\n"
      "        not equals to calculating the union of the results of IN and OUT.\n"
  },
  
  // interface to igraph_subgraph
  {"subgraph", (PyCFunction)igraphmodule_Graph_subgraph,
      METH_VARARGS | METH_KEYWORDS,
      "Returns a subgraph based on the given vertices.\n\n"
      "Keyword arguments:\n"
      "vertices -- a list containing the vertex IDs which should be included in the result.\n"
  },
  
  //////////////////////
  // LAYOUT FUNCTIONS //
  //////////////////////
  
  // interface to igraph_layout_circle
  {"layout_circle", (PyCFunction)igraphmodule_Graph_layout_circle,
      METH_VARARGS | METH_KEYWORDS,
      "Places the vertices of the graph uniformly on a circle.\n\n"
      "Returns the calculated coordinate pairs in a vector."
  },
  
  // interface to igraph_layout_sphere
  {"layout_sphere", (PyCFunction)igraphmodule_Graph_layout_sphere,
      METH_VARARGS | METH_KEYWORDS,
      "Places the vertices of the graph uniformly on a sphere.\n\n"
      "Returns the calculated coordinate pairs in a vector."
  },
  
  // interface to igraph_layout_kamada_kawai
  {"layout_kamada_kawai", (PyCFunction)igraphmodule_Graph_layout_kamada_kawai,
      METH_VARARGS | METH_KEYWORDS,
      "Places the vertices on a plane according to the Kamada-Kawai algorithm.\n\n"
      "This is a force directed layout, see Kamada, T. and Kawai, S.:\n"
      "An Algorithm for Drawing General Undirected Graphs.\n"
      "Information Processing Letters, 31/1, 7--15, 1989.\n\n"
      "Keyword arguments:\n"
      "maxiter -- the number of iterations to perform. Optional, defaults to 1000.\n"
      "sigma   -- the standard base deviation of the position change proposals.\n"
      "           Optional, defaults to the number of vertices * 0.25\n"
      "initemp -- initial temperature of the simulated annealing.\n"
      "           Optional, defaults to 10.\n"
      "coolexp -- cooling exponent of the simulated annealing.\n"
      "           Optional, defaults to 0.99\n"
      "kkconst -- the Kamada-Kawai vertex attraction constant. Optional,\n"
      "           defaults to the square of the number of vertices.\n"
  },
  
  // interface to igraph_layout_kamada_kawai_3d
  {"layout_kamada_kawai_3d", (PyCFunction)igraphmodule_Graph_layout_kamada_kawai_3d,
      METH_VARARGS | METH_KEYWORDS,
      "Places the vertices in the 3D space according to the Kamada-Kawai algorithm.\n\n"
      "For argument list, see Graph.layout_kamada_kawai"
  },
  
  // interface to igraph_layout_fruchterman_reingold
  {"layout_fruchterman_reingold", (PyCFunction)igraphmodule_Graph_layout_fruchterman_reingold,
      METH_VARARGS | METH_KEYWORDS,
      "Places the vertices on a 2D plane according to the Fruchterman-Reingold algorithm.\n\n"
      "This is a force directed layout, see Fruchterman, T. M. J. and Reingold, E. M.:\n"
      "Graph Drawing by Force-directed Placement.\n"
      "Software -- Practice and Experience, 21/11, 1129--1164, 1991\n\n"
      "Keyword arguments:\n"
      "maxiter    -- the number of iterations to perform. Optional, defaults to 500.\n"
      "maxdelta   -- the maximum distance to move a vertex in an iteration\n"
      "              Optional, defaults to the number of vertices\n"
      "area       -- the area of the square on which the vertices\n"
      "              will be placed. Optional, defaults to the square of\n"
      "              maxdelta.\n"
      "coolexp    -- the cooling exponent of the simulated annealing.\n"
      "              Optional, defaults to 0.99\n"
      "repulserad -- Determines the radius at which vertex-vertex repulsion\n"
      "              cancels out attraction of adjacent vertices. Optional,\n"
      "              defaults to maxiter * maxdelta\n"
  },
  
  // interface to igraph_layout_fruchterman_reingold_3d
  {"layout_fruchterman_reingold_3d", (PyCFunction)igraphmodule_Graph_layout_fruchterman_reingold_3d,
      METH_VARARGS | METH_KEYWORDS,
      "Places the vertices in the 3D space according to the Fruchterman-Reingold grid algorithm.\n\n"
      "For argument list, see Graph.layout_fruchterman_reingold"
  },
  
  // interface to igraph_layout_grid_fruchterman_reingold
  {"layout_grid_fruchterman_reingold", (PyCFunction)igraphmodule_Graph_layout_grid_fruchterman_reingold,
      METH_VARARGS | METH_KEYWORDS,
      "Places the vertices on a 2D plane according to the Fruchterman-Reingold grid algorithm.\n\n"
      "This is a modified version of a force directed layout, see\n"
      "Fruchterman, T. M. J. and Reingold, E. M.:\n"
      "Graph Drawing by Force-directed Placement.\n"
      "Software -- Practice and Experience, 21/11, 1129--1164, 1991\n"
      "The algorithm partitions the 2D space to a grid and vertex\n"
      "repulsion is then calculated only for vertices nearby.\n\n"
      "Keyword arguments:\n"
      "maxiter    -- the number of iterations to perform. Optional, defaults to 500.\n"
      "maxdelta   -- the maximum distance to move a vertex in an iteration\n"
      "              Optional, defaults to the number of vertices\n"
      "area       -- the area of the square on which the vertices\n"
      "              will be placed. Optional, defaults to the square of\n"
      "              maxdelta.\n"
      "coolexp    -- the cooling exponent of the simulated annealing.\n"
      "              Optional, defaults to 0.99\n"
      "repulserad -- Determines the radius at which vertex-vertex repulsion\n"
      "              cancels out attraction of adjacent vertices. Optional,\n"
      "              defaults to maxiter * maxdelta\n"
      "cellsize   -- the size of the grid cells.\n"
  },
  
  // interface to igraph_layout_lgl
  {"layout_lgl", (PyCFunction)igraphmodule_Graph_layout_lgl,
      METH_VARARGS | METH_KEYWORDS,
      "Places the vertices on a 2D plane according to the Large Graph Layout.\n\n"
      "Keyword arguments:\n"
      "maxiter    -- the number of iterations to perform. Optional, defaults to 500.\n"
      "maxdelta   -- the maximum distance to move a vertex in an iteration\n"
      "              Optional, defaults to the number of vertices\n"
      "area       -- the area of the square on which the vertices\n"
      "              will be placed. Optional, defaults to the square of\n"
      "              maxdelta.\n"
      "coolexp    -- the cooling exponent of the simulated annealing.\n"
      "              Optional, defaults to 0.99\n"
      "repulserad -- Determines the radius at which vertex-vertex repulsion\n"
      "              cancels out attraction of adjacent vertices. Optional,\n"
      "              defaults to maxiter * maxdelta\n"
      "cellsize   -- The size of the grid cells, one side of the square.\n"
      "              Optional.\n"
      "proot      -- The root vertex, this is placed first, its neighbors\n"
      "              in the first iteration, second neighbors in the second,\n"
      "              etc. Optional, defaults to a random vertex.\n"
  },
  
  // interface to igraph_layout_random
  {"layout_random", (PyCFunction)igraphmodule_Graph_layout_random,
      METH_VARARGS | METH_KEYWORDS,
      "Places the vertices of the graph randomly in a 2D space.\n\n"
      "Returns the \"calculated\" coordinate pairs in a vector."
  },
   
  // interface to igraph_layout_random_3d
  {"layout_random_3d", (PyCFunction)igraphmodule_Graph_layout_random_3d,
      METH_VARARGS | METH_KEYWORDS,
      "Places the vertices of the graph randomly in a 3D space.\n\n"
      "Returns the \"calculated\" coordinate triplets in a vector."
  },
   
  //////////////////////////////////////////////////////
  // CONVERT A GRAPH TO EDGE LIST OR ADJACENCY MATRIX //
  //////////////////////////////////////////////////////
  
  // interface to igraph_get_edgelist
  {"get_adjacency", (PyCFunction)igraphmodule_Graph_get_adjacency,
      METH_VARARGS | METH_KEYWORDS,
      "Returns the edge list of a graph.\n\n"
      "Keyword arguments:\n"
      "type -- either GET_ADJACENCY_LOWER (uses the lower triangle\n"
      "        of the matrix) or GET_ADJACENCY_UPPER (uses the upper\n"
      "        triangle) or GET_ADJACENCY_BOTH (uses both parts).\n"
      "        Optional, defaults to GET_ADJACENCY_BOTH, ignored for\n"
      "        directed graphs."
  },
  
  // interface to igraph_get_edgelist
  {"get_edgelist", (PyCFunction)igraphmodule_Graph_get_edgelist,
      METH_NOARGS,
      "Returns the edge list of a graph."
  },
  
  ///////////////////////////////
  // LOADING AND SAVING GRAPHS //
  ///////////////////////////////
  
  // interface to igraph_read_graph_edgelist
  {"Read_Edgelist", (PyCFunction)igraphmodule_Graph_Read_Edgelist,
      METH_VARARGS | METH_KEYWORDS | METH_CLASS,
      "Reads an edge list from a file and creates a graph based on it.\n"
      "Please note that the vertex indices are zero-based.\n\n"
      "Keyword arguments:\n"
      "f        -- the name of the file\n"
      "directed -- whether the generated graph should be directed.\n"
      "            Optional, defaults to True.\n\n"
  },
  // interface to igraph_read_graph_ncol
  {"Read_Ncol", (PyCFunction)igraphmodule_Graph_Read_Ncol,
      METH_VARARGS | METH_KEYWORDS | METH_CLASS,
      "Reads an .ncol file used by LGL, also useful for creating graphs\n"
      "from \"named\" (and optionally weighted) edge lists.\n\n"
      "This format is used by the Large Graph Layout program. See the\n"
      "documentation of LGL regarding the exact format description:\n"
      "http://bioinformatics.icmb.utexas.edu/bgl\n\n"
      "LGL originally cannot deal with graphs containing multiple or loop\n"
      "edges, but this condition is not checked here, as igraph is happy\n"
      "with these.\n\n"
      "Keyword arguments:\n"
      "f       -- the name of the file\n"
      "names   -- logical value. If True, the vertex names are added as a\n"
      "           vertex attribute called 'name'. Optional, defaults to\n"
      "           True.\n"
      "weights -- logical value. If True, the edge weights are added as an\n"
      "           edge attribute called 'weight'. Optional, defaults to\n"
      "           True.\n"
  },
  // interface to igraph_read_graph_lgl
  {"Read_Lgl", (PyCFunction)igraphmodule_Graph_Read_Lgl,
      METH_VARARGS | METH_KEYWORDS | METH_CLASS,
      "Reads an .lgl file used by LGL, also useful for creating graphs\n"
      "from \"named\" (and optionally weighted) edge lists.\n\n"
      "This format is used by the Large Graph Layout program. See the\n"
      "documentation of LGL regarding the exact format description:\n"
      "http://bioinformatics.icmb.utexas.edu/bgl\n\n"
      "LGL originally cannot deal with graphs containing multiple or loop\n"
      "edges, but this condition is not checked here, as igraph is happy\n"
      "with these.\n\n"
      "Keyword arguments:\n"
      "f       -- the name of the file\n"
      "names   -- logical value. If True, the vertex names are added as a\n"
      "           vertex attribute called 'name'. Optional, defaults to\n"
      "           True.\n"
      "weights -- logical value. If True, the edge weights are added as an\n"
      "           edge attribute called 'weight'. Optional, defaults to\n"
      "           True.\n"
  },
  // interface to igraph_read_graph_pajek
  {"Read_Pajek", (PyCFunction)igraphmodule_Graph_Read_Pajek,
      METH_VARARGS | METH_KEYWORDS | METH_CLASS,
      "Reads a Pajek format file and creates a graph based on it.\n"
      "Keyword arguments:\n"
      "f -- the name of the file\n"
  },
  // interface to igraph_read_graph_graphml
  {"Read_GraphML", (PyCFunction)igraphmodule_Graph_Read_GraphML,
      METH_VARARGS | METH_KEYWORDS | METH_CLASS,
      "Reads a GraphML format file and creates a graph based on it.\n"
      "Keyword arguments:\n"
      "f        -- the name of the file\n"
      "directed -- whether the graph should be directed. Please note that\n"
      "            if you ask for a directed graph, but the GraphML file\n"
      "            contains an undirected graph (with some optional directed\n"
      "            edges), then the undirected edges will be added in both\n"
      "            directions and the resulting Graph object will be a\n"
      "            directed graph. However, if you ask for an undirected\n"
      "            graph and the GraphML file contains a directed one,\n"
      "            all of the edges will be undirected (so your request\n"
      "            takes precedence over what the file states). Optional,\n"
      "            defaults to true\n"
      "index    -- if the GraphML file contains multiple graphs, specifies\n"
      "            the one which should be loaded. Graph indices start from\n"
      "            zero, so if you want to load the first graph, specify 0\n"
      "            here. Optional, defaults to zero.\n"
  },
  // interface to igraph_write_graph_edgelist
  {"write_edgelist", (PyCFunction)igraphmodule_Graph_write_edgelist,
      METH_VARARGS | METH_KEYWORDS,
      "Writes the edge list of a graph to a file. Directed edges are\n"
      "written in (from, to) order.\n\n"
      "Keyword arguments:\n"
      "f -- the name of the file to be written\n"
  },
  // interface to igraph_write_graph_ncol
  {"write_ncol", (PyCFunction)igraphmodule_Graph_write_ncol,
      METH_VARARGS | METH_KEYWORDS,
      "Writes the edge list of a graph to a file in .ncol format.\n"
      "Note that multiple edges and/or loops break the LGL software,\n"
      "but igraph does not check for this condition. Unless you know\n"
      "that the graph does not have multiple edges and/or loops, it\n"
      "is wise to call simplify() before saving.\n\n"
      "Keyword arguments:\n"
      "f       -- the name of the file to be written\n"
      "names   -- the name of the vertex attribute containing the name\n"
      "           of the vertices. Optional, defaults to 'name'. If you\n"
      "           don't want to store vertex names, supply None here.\n"
      "weights -- the name of the edge attribute containing the weight\n"
      "           of the vertices. Optional, defaults to 'weight'. If you\n"
      "           don't want to store weights, supply None here.\n"
  },
  // interface to igraph_write_graph_lgl
  {"write_lgl", (PyCFunction)igraphmodule_Graph_write_lgl,
      METH_VARARGS | METH_KEYWORDS,
      "Writes the edge list of a graph to a file in .lgl format.\n"
      "Note that multiple edges and/or loops break the LGL software,\n"
      "but igraph does not check for this condition. Unless you know\n"
      "that the graph does not have multiple edges and/or loops, it\n"
      "is wise to call simplify() before saving.\n\n"
      "Keyword arguments:\n"
      "f        -- the name of the file to be written\n"
      "names    -- the name of the vertex attribute containing the name\n"
      "            of the vertices. Optional, defaults to 'name'. If you\n"
      "            don't want to store vertex names, supply None here.\n"
      "weights  -- the name of the edge attribute containing the weight\n"
      "            of the vertices. Optional, defaults to 'weight'. If you\n"
      "            don't want to store weights, supply None here.\n"
      "isolates -- whether to include isolated vertices in the output.\n"
      "            Optional, defaults to True.\n"
  },
  // interface to igraph_write_graph_edgelist
  {"write_graphml", (PyCFunction)igraphmodule_Graph_write_graphml,
      METH_VARARGS | METH_KEYWORDS,
      "Writes the graph to a GraphML file.\n\n"
      "Keyword arguments:\n"
      "f -- the name of the file to be written\n"
  },

  ////////////////////////
  // ATTRIBUTE HANDLING //
  ////////////////////////
  {"attributes", (PyCFunction)igraphmodule_Graph_attributes,
      METH_NOARGS,
      "Returns the attribute list of the graph\n"
  },
  {"vertex_attributes", (PyCFunction)igraphmodule_Graph_vertex_attributes,
      METH_NOARGS,
      "Returns the attribute list of the graph's vertices\n"
  },
  {"edge_attributes", (PyCFunction)igraphmodule_Graph_edge_attributes,
      METH_NOARGS,
      "Returns the attribute list of the graph's edges\n"
  },
  
  ////////////////////////////////////
  // INTERNAL/DEVELOPMENT FUNCTIONS //
  ////////////////////////////////////
  {"__graph_as_cobject__", (PyCFunction)igraphmodule_Graph___graph_as_cobject__,
      METH_VARARGS | METH_KEYWORDS,
      "Returns the igraph graph encapsulated by the Python object as\n"
      "a PyCObject (which is barely a regular C pointer). This function\n"
      "should not be used directly by igraph users, it is useful only\n"
      "in the case when the underlying igraph object must be passed to\n"
      "another C code through Python.\n\n"
      /*"Keyword arguments:\n"
      "ref -- increases the reference count of the graph when True.\n"
      "       Optional, defaults to False.\n"*/
  },
  {"__register_destructor__", (PyCFunction)igraphmodule_Graph___register_destructor__,
      METH_VARARGS | METH_KEYWORDS,
      "Registers a destructor to be called when the object is freed by "
      "Python. This function should not be used directly by igraph users."
  },
  {NULL}

}
;

#endif
