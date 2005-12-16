/* -*- mode: C -*-  */
/* 
   IGraph library.
   Copyright (C) 2003, 2004, 2005  Gabor Csardi <csardi@rmki.kfki.hu>
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
   Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

*/

#ifndef RESTGAME_H
#define RESTGAME_H

#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
# define __BEGIN_DECLS extern "C" {
# define __END_DECLS }
#else
# define __BEGIN_DECLS /* empty */
# define __END_DECLS /* empty */
#endif

__BEGIN_DECLS

#ifndef _GNU_SOURCE
# define _GNU_SOURCE
#endif

#include "types.h"
#include "attributes.h"
#include "error.h"

#include <stdio.h> 		/* FILE */

/**
 * \defgroup errorhandling Error handling
 * \brief How igraph handles errors.
 * 
 * This very first section describes the \a igraph error handling.
 * The error handling is based on the code of the GNU GSL. 
 * The following base principles are used:
 * -# All functions which may run into an error return an
 *   <code>int</code>, the code of the error.
 * -# A function which notices an error should call the error handler
 *   and (if the error handler returns) return to its caller with the
 *   error code. Functions performing their jobs successfully return
 *   with IGRAPH_SUCCESS.
 * -# The error handler is called by using the igraph_error() function.
 *   There is also an IGRAPH_ERROR macro for convenience, it is
 *   recommended to use this. This macro calles the error handler and
 *   returns with the error code.
 * -# Functions returning anything else than IGRAPH_SUCCESS always
 *   leave their parameters <i>>semantically</i> unchanged.
 * -# The \a igraph error codes and their explanation:
 *   - <b>IGRAPH_SUCCESS</b>: the function performed its task
 *     successfully. It can be assumed that this error code is the
 *     same as integer zero.
 *   - <b>IGRAPH_FAILURE</b>: some unknown error happened. 
 *   - <b>IGRAPH_ENOMEM</b>: a memory allocation operation failed.
 *   - <b>IGRAPH_PARSEERROR</b>: a syntactic error was found in a file
 *     while parsing.
 *   - <b>IGRAPH_EINVAL</b>: invalid value of a parameter.
 *   - <b>IGRAPH_EXISTS</b>: the specified attribute already exists so
 *     it cannot be added.
 *   - <b>IGRAPH_EINVEVECTOR</b>: invalid edge vector, it either
 *     contains a too high or negative value.
 *   - <b>IGRAPH_EINVVID</b>: invalid vertex id. It is either too big
 *     or negative.
 *   - <b>IGRAPH_NONSQUARE</b>: non-square matrix was passed when a
 *     square matrix was expected.
 *   - <b>IGRAPH_EINVMODE</b>: invalid <code>mode</code> parameter.
 *   - <b>IGRAPH_EFILE</b>: a file operation (like opening a file)
 *     failed. 
 * -# It is possible to define an error handler function which will be 
 *   called after an error is noticed. This is possible with the
 *   igraph_set_error_handler() function, see this for details. The
 *   default error handler prints an error message and aborts the
 *   program. 
 */

/** \defgroup types Basic data types.
 * \brief Important data types for vectors, matrices and graphs.
 */ 

/** \defgroup interface The basic igraph interface 
 * \brief This is the very minimal API in \a igraph. All the other
 * functions use this minimal set for creating and manipulating the
 * graphs. 
 * 
 * This is a very important principle since it makes possible to
 * implement other data representations by implementing only this
 * minimal set.
 */

/** \defgroup generators Graph generators 
 * \brief Graph generators create graphs.
 * 
 * Almost all functions which create graph objects are documented
 * here. The exceptions are <code>subgraph</code> and alike, these 
 * create graphs based on another graph.
 */

/** \defgroup structural Structural properties of graphs 
 * \brief These functions usually calculate some structural property
 * of a graph, like its diameter, the degree of the nodes, etc.
 */

/**
 * \defgroup attributes Graph, vertex and edge attributes
 * \brief Attributes are numberic or textual labels assigned to 
 * graphs, vertices or edges.
 * 
 * Attributes are associated values of graph, vertices or edges. A
 * graph attribute can contain the date of its creations, vertex
 * attributes describe the color to be used for the vertices when the
 * graph is plotted, edge attributes can simply be the weights of the
 * edges in a weighted graph.
 * 
 * The name space for graph, vertex and edge attributes are different,
 * thus the ``color'' vertex attribute has nothing to do with the
 * ``color'' edge attribute.
 *
 * In order to use an attribute, first it has to be added by the
 * igraph_add_graph_attribute(), igraph_add_vertex_attribute() or
 * igraph_add_edge_attribute() function. The value(s) of the added
 * attribute is/are undefined.
 * 
 * The value of the attribute can be set after it was added, and also
 * it can be requested, see the documentation of the specific
 * functions.
 * 
 * The attribute can be removed id it is not needed any more.
 * 
 * The attributes are handled properly if vertices and/or edges are
 * added or removed to or from a graph, the attributes are also copied
 * if the graph is copied by igraph_copy().
 * 
 * There are two types of attributes, numeric
 * (<b>IGRAPH_ATTRIBUTE_NUM</b>) and  string
 * (<b>IGRAPH_ATTRIBUTE_STR</b>). These types cannot be mixed for an
 * attribute, ie. if the "id" vertex attribute is a string attribute
 * then it is a string attribute for all vertices.
 *
 * Attribute handling does not change the time complexity of any
 * functions in the igraph library.
 */

/**
 * \defgroup layout Layout functions for graph drawing
 * \brief These functions place the vertices of a graph in 2D and 3D,
 * trying to please the eye as much as possible.
 */

/**
 * \defgroup loadsave Reading and writing graphs from and to files
 * \brief These functions can write a graph to a file, or read a graph
 * from a file.
 * 
 * Note that as \a igraph uses the traditional C streams it is
 * possible to read/write files from/to memory, at least on GNU
 * operating systems supporting "non-standard" streams.
 */

/**
 * \defgroup conversion Convert a graph to an edge list or adjacency matrix
 * \brief 
 */

/**
 * \defgroup nongraph Non-graph related functions
 * \brief These functions are not directly graph related but they 
 * are either used by the graph related routines or are just useful 
 * themselves.
 */

/** \defgroup internal Internal functions and data types
 * \brief Everything documented here should be of interest of
 * \a igraph developers only.
 */


/** 
 * \mainpage The <em>igraph</em> library manual
 *
 * \section toc Table of Contents
 * - \ref intro
 * - \ref features
 * - \ref installation
 * - \ref license
 * - \ref tutorial
 *
 * \section intro Introduction
 * This is another library for creating and manipulating graphs. Its
 * aim is performance above all. There is simply no
 * other graph library out there which can be handle graphs of the
 * size the author is confronted with efficiently. Whenever possible \a igraph
 * tries to be also user friendly and portable. 
 *
 * \a igraph started as an additional package to the GNU R
 * statistical environment, and still some functions which are hard to
 * implement in C are available only in R (like interactive
 * graphics). Most functions are however now written in C and they can
 * be compiled without R as a separate library. Still the author
 * highly recommends R, as it is a sophisticated, well desinged,
 * modern system.
 *
 * The advantage of the C version is that it is portable, can
 * be used with other mathematical software (as soon as somebody
 * writes the interfaces to it), and can be used for large scale
 * scientific computation in the Condor environment we are
 * using here at KFKI.
 *
 * \section features Features
 * - \a igraph contains functions for generating regular and
 *   random graphs according to known algorithms and models in the
 *   network theory literature.
 * - \a igraph provides routines for manipulating graphs,
 *   adding and removing edges and vertices. 
 * - a set of structural property calculation functions like
 *   degree, betweenness, etc. are also included. 
 * - force based layout generators are included for smaller graphs,
 *   another method is expected to be added for large graphs soon.
 * - a set of conversion functions are also included and will be
 *   extended shortly.
 * - \a igraph iterators provide a simple and efficient way
 *   of walking through graphs.
 * - the documentation will be finished shortly, for every function
 *   its time complexity is stated.
 *
 * \section installation Installation
 * First download the latest version of \a igraph and
 * uncompress it to a temporary directory:
 * \verbatim wget http://geza.kzoo.edu/~csardi/igraph_latest.tar.gz
 tar xzf igraph_latest.tar.gz 
 cd igraph-latest 
\endverbatim
 * \subsection installation-r Installation as an R package
 * This is very simple, the \a igraph directory in the
 * tarball contains a complete R package. Eg. in Linux you can install
 * it like
 * \verbatim R CMD INSTALL -l ~/.R/library igraph \endverbatim
 * 
 * Alternatively you can download the R package version of \a igraph only.
 *
 * Also, consult your R documentation.
 * 
 * \subsection installation-lib Installation as a library
 * Simply typing
 * \verbatim ./configure
 make
 make install \endverbatim 
 * (the latter as root) should work on most systems. You can try
 * \verbatim ./configure --help \endverbatim
 * to see installations options, and reading the <code>INSTALL</code> 
 * file. 
 *
 * \section license License
 *
 * Copyright (C) 2003, 2004, 2005  Gabor Csardi <csardi@rmki.kfki.hu>
 * MTA RMKI, Konkoly-Thege Miklos st. 29-33, Budapest 1121, Hungary
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 *
 * \section tutorial A short tutorial
 * This is under development. Until it is done, you can just read
 * through the documentation of each modules in the Modules section.
 *
 */

/**
 * \ingroup internal
 * \struct igraph_t
 * \brief The internal data structure for storing graphs.
 *
 * It is simple and efficient. It has the following members:
 * - <b>n</b> The number of vertices, reduntant.
 * - <b>directed</b> Whether the graph is directed.
 * - <b>from</b> The first column of the edge list.
 * - <b>to</b> The second column of the edge list.
 * - <b>oi</b> The index of the edge list by the first column. Thus
 *   the first edge according to this order goes from
 *   \c from[oi[0]] to \c to[oi[0]]. The length of
 *   this vector is the same as the number of edges in the graph.
 * - <b>ii</b> The index of the edge list by the second column. 
 *   The length of this vector is the same as the number of edges.
 * - <b>os</b> Contains pointers to the edgelist (\c from
 *   and \c to for every vertex. The first edge \em from
 *   vertex \c v is edge no. \c from[oi[os[v]]] if 
 *   \c os[v]<os[v+1]. If \c os[v]==os[v+1] then
 *   there are no edges \em from node \c v. Its length is
 *   the number of vertices plus one, the last element is always the 
 *   same as the number of edges and is contained only to ease the
 *   queries.
 * - <b>is</b> This is basically the same as <b>os</b>, but this time
 *   for the incoming edges.
 * 
 * For undirected graph, the same edge list is stored, ie. an
 * undirected edge is stored only once, and for checking whether there
 * is an undirected edge from \c v1 to \c v2 one
 * should search for both \c from=v1, \c to=v2 and 
 * \c from=v2, \c to=v1.
 *
 * The storage requirements for a graph with \c |V| vertices
 * and \c |E| edges is \c O(|E|+|V|).
 */
typedef struct igraph_s {
  integer_t n;
  bool_t directed;
  vector_t from;
  vector_t to;
  vector_t oi;
  vector_t ii;
  vector_t os;
  vector_t is;
  igraph_attribute_list_t gal;
  igraph_attribute_list_t val;
  igraph_attribute_list_t eal;
} igraph_t;

/* -------------------------------------------------- */
/* Iterators                                          */
/* -------------------------------------------------- */

/* Vertex iterators */
#define IGRAPH_ITERATOR_VS_ALL          0
#define IGRAPH_ITERATOR_VS_NONE         1
#define IGRAPH_ITERATOR_VS_ADJ          2
#define IGRAPH_ITERATOR_VS_VECTOR       3
#define IGRAPH_ITERATOR_VS_1            4
#define IGRAPH_ITERATOR_VS_SEQ          5
#define IGRAPH_ITERATOR_VS_RW           6
#define IGRAPH_ITERATOR_VS_RW1          7

#define IGRAPH_ITERATOR_ES_ALL          0
#define IGRAPH_ITERATOR_ES_FROMORDER    1
#define IGRAPH_ITERATOR_ES_ADJ          2
#define IGRAPH_ITERATOR_ES_NONE         3
#define IGRAPH_ITERATOR_ES_1            4
#define IGRAPH_ITERATOR_ES_VECTOR       5

/* -------------------------------------------------- */
/* Constants                                          */
/* -------------------------------------------------- */

typedef enum { IGRAPH_UNDIRECTED=0, IGRAPH_DIRECTED=1 } igraph_i_directed_t;

typedef enum { IGRAPH_NO_LOOPS=0, IGRAPH_LOOPS=1 } igraph_i_loops_t;

typedef enum { IGRAPH_OUT=1, IGRAPH_IN=2, IGRAPH_ALL=3,
	       IGRAPH_TOTAL=3 } igraph_neimode_t;

typedef enum { IGRAPH_WEAK=1, IGRAPH_STRONG=2 } igraph_connectedness_t;

typedef enum { IGRAPH_ADJ_DIRECTED=0, 
	       IGRAPH_ADJ_UNDIRECTED=1, IGRAPH_ADJ_MAX=1,
               IGRAPH_ADJ_UPPER, IGRAPH_ADJ_LOWER, IGRAPH_ADJ_MIN,
	       IGRAPH_ADJ_PLUS } igraph_adjacency_t;

typedef enum { IGRAPH_STAR_OUT=0, IGRAPH_STAR_IN,
	       IGRAPH_STAR_UNDIRECTED } igraph_star_mode_t;

typedef enum { IGRAPH_TREE_OUT=0, IGRAPH_TREE_IN,
	       IGRAPH_TREE_UNDIRECTED } igraph_tree_mode_t;

typedef enum { IGRAPH_ERDOS_RENYI_GNP=0, 
	       IGRAPH_ERDOS_RENYI_GNM } igraph_erdos_renyi_t;

typedef enum { IGRAPH_GET_ADJACENCY_UPPER=0,
	       IGRAPH_GET_ADJACENCY_LOWER,
	       IGRAPH_GET_ADJACENCY_BOTH } igraph_get_adjacency_t;

typedef enum { IGRAPH_DEGSEQ_SIMPLE=0 } igraph_degseq_t;

typedef enum { IGRAPH_TRANSITIVITY_UNDIRECTED=0 } igraph_transitivity_type_t;

typedef enum { IGRAPH_FILEFORMAT_EDGELIST=0,
	       IGRAPH_FILEFORMAT_NCOL,
	       IGRAPH_FILEFORMAT_PAJEK } igraph_fileformat_type_t;

/* -------------------------------------------------- */
/* Vertex sequences                                   */
/* -------------------------------------------------- */

#define IGRAPH_I_STDATA_SIZE 10

struct igraph_vs_t;

typedef struct igraph_i_vstable_t {
  void (*next)(igraph_t *graph, struct igraph_vs_t *vs);
  bool_t (*end)(igraph_t *graph, struct igraph_vs_t *vs);
  void (*reset)(igraph_t *graph, struct igraph_vs_t *vs);
  integer_t(*get)(igraph_t *graph, struct igraph_vs_t *vs);
  int (*unfold)(igraph_t *graph, struct igraph_vs_t *vs);
  void (*destroy)(struct igraph_vs_t *vs);
} igraph_i_vstable_t;

typedef struct igraph_vs_t {
  int type;
  real_t stdata[IGRAPH_I_STDATA_SIZE]; /* working storage */
  void *pdata;			   /* other storage   */
  igraph_i_vstable_t *table;	   /* the table of functions */
  vector_t *v;                     /* vector representation (if any) */
  bool_t vdestroy;                 /* whether the vector should be freed */
  bool_t view;
} igraph_vs_t;

extern igraph_i_vstable_t *igraph_i_vstable[8];

/* -------------------------------------------------- */
/* Edge sequences                                   */
/* -------------------------------------------------- */

struct igraph_es_t;

typedef struct igraph_i_estable_t {
  void (*next)(igraph_t *graph, struct igraph_es_t *es);
  bool_t (*end)(igraph_t *graph, struct igraph_es_t *es);
  void (*reset)(igraph_t *graph, struct igraph_es_t *es);
  integer_t(*get)(igraph_t *graph, struct igraph_es_t *es);
  integer_t (*from)(igraph_t *graph, struct igraph_es_t *es);
  integer_t (*to)(igraph_t *graph, struct igraph_es_t *es);
  int (*unfold)(igraph_t *graph, struct igraph_es_t *es);
  void (*destroy)(struct igraph_es_t *es);
} igraph_i_estable_t;

typedef struct igraph_es_t {
  int type;
  real_t stdata[IGRAPH_I_STDATA_SIZE]; /* working storage */
  void *pdata;			   /* other storage   */
  igraph_i_estable_t *table;	   /* the table of functions */
  vector_t *v;                     /* vector representation (if any) */
  bool_t vdestroy;                 /* whether the vector should be freed */
  bool_t view;
} igraph_es_t;

extern igraph_i_estable_t *igraph_i_estable[6];

/* -------------------------------------------------- */
/* Interface                                          */
/* -------------------------------------------------- */

struct igraph_eit_t;

int igraph_empty(igraph_t *graph, integer_t n, bool_t directed);
int igraph_destroy(igraph_t *graph);
int igraph_copy(igraph_t *to, igraph_t *from);
int igraph_add_edges(igraph_t *graph, vector_t *edges);
int igraph_add_vertices(igraph_t *graph, integer_t nv);
int igraph_delete_edges(igraph_t *graph, vector_t *edges); /* eee */
int igraph_delete_vertices(igraph_t *graph, igraph_vs_t vertices);
integer_t igraph_vcount(igraph_t *graph);
integer_t igraph_ecount(igraph_t *graph);
int igraph_neighbors(igraph_t *graph, vector_t *neis, integer_t vid, 
		     igraph_neimode_t mode); 
bool_t igraph_is_directed(igraph_t *graph);
int igraph_degree(igraph_t *graph, vector_t *res, igraph_vs_t vids, 
		  igraph_neimode_t mode, bool_t loops);

/* -------------------------------------------------- */
/* Vertex sequences, contd.                           */
/* -------------------------------------------------- */

/* Iterator view, for shorthands */
int igraph_vs_create_view_as_vector(igraph_t *graph, igraph_vs_t *vs,
				    igraph_vs_t *newvs);
  
/* Generics */
void igraph_vs_next(igraph_t *graph, igraph_vs_t *vs);
bool_t igraph_vs_end(igraph_t *graph, igraph_vs_t *vs);
void igraph_vs_reset(igraph_t *graph, igraph_vs_t *vs);
integer_t igraph_vs_get(igraph_t *graph, igraph_vs_t *vs);
int igraph_vs_unfold(igraph_t *graph, igraph_vs_t *vs);
void igraph_vs_destroy(igraph_vs_t *vs);

/* Simple vertex iterator, all vertices */
igraph_vs_t igraph_vs_all(igraph_t *graph);

/* Adjacent vertices of a vertex */
igraph_vs_t igraph_vs_adj(igraph_t *graph,
			  integer_t vid, igraph_neimode_t mode);
void igraph_vs_adj_set(igraph_t *graph, igraph_vs_t *vs,
		       integer_t vid, igraph_neimode_t mode);

/* Random walker */
igraph_vs_t igraph_vs_rw(igraph_t *graph,
			 integer_t vid, igraph_neimode_t mode);
long int igraph_vs_rw_length(igraph_t *graph, igraph_vs_t *vs);

/* Random walker with one unit memory */
igraph_vs_t igraph_vs_rw1(igraph_t *graph,
			  integer_t vid, igraph_neimode_t mode);
long int igraph_vs_rw1_length(igraph_t *graph, igraph_vs_t *vs);

/* Empty iterator */
igraph_vs_t igraph_vs_none(igraph_t *graph);

/* Single vertex iterator */
igraph_vs_t igraph_vs_1(igraph_t *igraph, integer_t vid);

/* A sequence of vertices */
igraph_vs_t igraph_vs_seq(igraph_t *igraph, integer_t from, 
			  integer_t to);

/* Iterates over a vector */
igraph_vs_t igraph_vs_vector(igraph_t *igraph, vector_t *vids);

/* Shorthands, a bit dirty */
#define IGRAPH_VS_ALL             (igraph_vs_all(0))
#define IGRAPH_VS_NONE            (igraph_vs_none(0))
#define IGRAPH_VS_1(id)           (igraph_vs_1(0, (id)))
#define IGRAPH_VS_VECTOR(vector)  (igraph_vs_vector(0,(vector)))

/* -------------------------------------------------- */
/* Edge sequences contd.                              */
/* -------------------------------------------------- */

/* Iterator view, for shorthands */
int igraph_es_create_view_as_vector(igraph_t *graph, igraph_es_t *es,
				    igraph_es_t *newes);

/* Generics */
void igraph_es_next(igraph_t *graph, igraph_es_t *es);
bool_t igraph_es_end(igraph_t *graph, igraph_es_t *es);
void igraph_es_reset(igraph_t *graph, igraph_es_t *es);
integer_t igraph_es_get(igraph_t *graph, igraph_es_t *es);
integer_t igraph_es_from(igraph_t *graph, igraph_es_t *es);
integer_t igraph_es_to(igraph_t *graph, igraph_es_t *es);
int igraph_es_unfold(igraph_t *graph, igraph_es_t *es);
void igraph_es_destroy(igraph_es_t *es);

/* Simple edge iterator */
igraph_es_t igraph_es_all(igraph_t *graph);

/* Empty edge iterator */
igraph_es_t igraph_es_none(igraph_t *graph);

/* Single edge edge iterator */
igraph_es_t igraph_es_1(igraph_t *graph, integer_t id);

/* Iterators between vertex sets */
igraph_es_t igraph_es_fromto(igraph_t *igraph, igraph_vs_t from,
			     igraph_vs_t to, igraph_neimode_t mode);

/* Sorted edge iterator */
igraph_es_t igraph_es_fromorder(igraph_t *graph);

/* Adjacent edges of a vertex */
igraph_es_t igraph_es_adj(igraph_t *graph,
			  integer_t vid, igraph_neimode_t mode);
void igraph_es_adj_set(igraph_t *graph, igraph_es_t *es,
			integer_t vid, igraph_neimode_t mode);
integer_t igraph_es_adj_vertex(igraph_t *graph, igraph_es_t *es);

/* View of a vector */
igraph_es_t igraph_es_vector(igraph_t *graph, vector_t *eids);

/* Shorthands, a bit dirty */
#define IGRAPH_ES_ALL             (igraph_es_all(0))
#define IGRAPH_ES_NONE            (igraph_es_none(0))
#define IGRAPH_ES_1(id)           (igraph_es_1(0, (id)))
#define IGRAPH_ES_VECTOR(vector)  (igraph_es_vector(0,(vector)))

/* -------------------------------------------------- */
/* Attributes                                         */
/* -------------------------------------------------- */

int igraph_add_graph_attribute(igraph_t *graph, const char *name,
			       igraph_attribute_type_t type);
int igraph_remove_graph_attribute(igraph_t *graph, const char *name);
int igraph_get_graph_attribute(igraph_t *graph, const char *name, 
			       void **value,
			       igraph_attribute_type_t *type);
int igraph_set_graph_attribute(igraph_t *graph, const char *name, 
			       void *value);
int igraph_list_graph_attributes(igraph_t *graph, igraph_strvector_t *l,
				 vector_t *types);
int igraph_get_graph_attribute_type(igraph_t *graph, const char *name, 
				    igraph_attribute_type_t *type);
bool_t igraph_has_graph_attribute(igraph_t *graph, const char *name);
				  

int igraph_add_vertex_attribute(igraph_t *graph, const char *name,
				igraph_attribute_type_t type);
int igraph_remove_vertex_attribute(igraph_t *graph, const char *name);
int igraph_get_vertex_attribute(igraph_t *graph, const char *name, 
				long int v, void **value,
				igraph_attribute_type_t *type);
int igraph_set_vertex_attribute(igraph_t *graph, const char *name, 
				long int v, void *value);
int igraph_get_vertex_attributes(igraph_t *graph, const char *name, 
				 igraph_vs_t v, void **value);
int igraph_set_vertex_attributes(igraph_t *graph, const char *name, 
				 igraph_vs_t v, void *value);
int igraph_list_vertex_attributes(igraph_t *graph, igraph_strvector_t *l,
				  vector_t *types);
int igraph_get_vertex_attribute_type(igraph_t *graph, const char *name, 
				     igraph_attribute_type_t *type);
bool_t igraph_has_vertex_attribute(igraph_t *graph, const char *name);

int igraph_add_edge_attribute(igraph_t *graph, const char *name, 
			      igraph_attribute_type_t type);
int igraph_remove_edge_attribute(igraph_t *graph, const char *name);
int igraph_get_edge_attribute(igraph_t *graph, const char *name, 
			      long int e, void **value,
			      igraph_attribute_type_t *type);
int igraph_set_edge_attribute(igraph_t *graph, const char *name, 
			      long int e, void *value);
int igraph_get_edge_attributes(igraph_t *graph, const char *name, 
			       vector_t *e, void **value); /* eee */
int igraph_set_edge_attributes(igraph_t *graph, const char *name, 
			       vector_t *e, void *value); /* eee */
int igraph_list_edge_attributes(igraph_t *graph, igraph_strvector_t *l,
				vector_t *types);
int igraph_get_edge_attribute_type(igraph_t *graph, const char *name, 
				   igraph_attribute_type_t *type);
bool_t igraph_has_edge_attribute(igraph_t *graph, const char *name);

/* -------------------------------------------------- */
/* Constructors, deterministic                        */
/* -------------------------------------------------- */

int igraph_create(igraph_t *graph, vector_t *edges, integer_t n, 
		  bool_t directed);
int igraph_adjacency(igraph_t *graph, matrix_t *adjmatrix,
		     igraph_adjacency_t mode);
int igraph_star(igraph_t *graph, integer_t n, igraph_star_mode_t mode, 
		integer_t center);
int igraph_lattice(igraph_t *graph, vector_t *dimvector, integer_t nei, 
		   bool_t directed, bool_t mutual, bool_t circular);
int igraph_ring(igraph_t *graph, integer_t n, bool_t directed, 
		bool_t mutual, bool_t circular);
int igraph_tree(igraph_t *graph, integer_t n, integer_t children, 
		igraph_tree_mode_t type);
int igraph_full(igraph_t *graph, integer_t n, bool_t directed, bool_t loops);

/* -------------------------------------------------- */
/* Constructors, games (=stochastic)                  */
/* -------------------------------------------------- */

int igraph_barabasi_game(igraph_t *graph, integer_t n, integer_t m, 
			 vector_t *outseq, bool_t outpref, bool_t directed);
int igraph_erdos_renyi_game(igraph_t *graph, igraph_erdos_renyi_t type,
			    integer_t n, real_t p,
			    bool_t directed, bool_t loops);
int igraph_erdos_renyi_game_gnp(igraph_t *graph, integer_t n, real_t p,
				bool_t directed, bool_t loops);
int igraph_degree_sequence_game(igraph_t *graph, vector_t *out_deg,
				vector_t *in_deg, igraph_degseq_t method);
int igraph_growing_random_game(igraph_t *graph, integer_t n, 
			       integer_t m, bool_t directed, bool_t citation);
int igraph_aging_prefatt_game(igraph_t *graph, integer_t n, integer_t m,
			      integer_t aging_type, real_t aging_exp);

/* -------------------------------------------------- */
/* Basic query functions                              */
/* -------------------------------------------------- */

bool_t igraph_are_connected(igraph_t *graph, integer_t v1, integer_t v2);

/* -------------------------------------------------- */
/* Structural properties                              */
/* -------------------------------------------------- */

int igraph_diameter(igraph_t *graph, integer_t *res, 
		    bool_t directed, bool_t unconn);
int igraph_minimum_spanning_tree_unweighted(igraph_t *graph, igraph_t *mst);
int igraph_minimum_spanning_tree_prim(igraph_t *graph, igraph_t *mst,
				      vector_t *weights);
int igraph_closeness(igraph_t *graph, vector_t *res, igraph_vs_t vids, 
		     igraph_neimode_t mode);
int igraph_shortest_paths(igraph_t *graph, matrix_t *res, 
			  igraph_vs_t from, igraph_neimode_t mode);
int igraph_get_shortest_paths(igraph_t *graph, vector_t *res,
			      integer_t from, igraph_neimode_t mode);
int igraph_subcomponent(igraph_t *graph, vector_t *res, real_t vid, 
			igraph_neimode_t mode);	
int igraph_betweenness (igraph_t *graph, vector_t *res, igraph_vs_t vids, 
			bool_t directed);
int igraph_edge_betweenness (igraph_t *graph, vector_t *result, 
			     bool_t directed); /* eee + add */
int igraph_subgraph(igraph_t *graph, igraph_t *res, igraph_vs_t vids);
int igraph_average_path_length(igraph_t *graph, real_t *res,
			       bool_t directed, bool_t unconn);
int igraph_simplify(igraph_t *graph, bool_t multiple, bool_t loops);
int igraph_transitivity(igraph_t *graph, vector_t *res, 
			igraph_transitivity_type_t type); /* vvv + add */

/* TODO: degree.distribution (?) */

/* -------------------------------------------------- */
/* Components                                         */
/* -------------------------------------------------- */

int igraph_clusters(igraph_t *graph, vector_t *membership, vector_t *csize, 
		    igraph_connectedness_t mode);
int igraph_is_connected(igraph_t *graph, bool_t *res, 
			igraph_connectedness_t mode);
/* TODO: cluster.distribution (?) */

/* -------------------------------------------------- */
/* Layouts                                            */
/* -------------------------------------------------- */

int igraph_layout_random(igraph_t *graph, matrix_t *res);
int igraph_layout_circle(igraph_t *graph, matrix_t *res);
int igraph_layout_fruchterman_reingold(igraph_t *graph, matrix_t *res,
				       integer_t niter, real_t maxdelta,
				       real_t area, real_t coolexp, 
				       real_t repulserad, bool_t use_seed);
int igraph_layout_kamada_kawai(igraph_t *graph, matrix_t *res,
			       integer_t niter, real_t sigma, 
			       real_t initemp, real_t coolexp,
			       real_t kkconst);
int igraph_layout_springs(igraph_t *graph, matrix_t *res,
			  real_t mass, real_t equil, real_t k,
			  real_t repeqdis, real_t kfr, bool_t repulse);

/* -------------------------------------------------- */
/* Centrality                                         */
/* -------------------------------------------------- */

/* TODO: evcent */

/* -------------------------------------------------- */
/* Cocitation                                         */
/* -------------------------------------------------- */

int igraph_cocitation(igraph_t *graph, matrix_t *res, igraph_vs_t vids);
int igraph_bibcoupling(igraph_t *graph, matrix_t *res, igraph_vs_t vids);

/* -------------------------------------------------- */
/* Community Structure                                */
/* -------------------------------------------------- */

/* TODO: eb.community */
/* TODO: cut.community */
/* TODO: edge.type.matrix */
/* TODO: modularity */
/* TODO:  */

/* -------------------------------------------------- */
/* Conversion                                         */
/* -------------------------------------------------- */

int igraph_get_adjacency(igraph_t *graph, matrix_t *res,
			 igraph_get_adjacency_t type);
int igraph_get_edgelist(igraph_t *graph, vector_t *res, bool_t bycol);

/* -------------------------------------------------- */
/* Read and write foreign formats                     */
/* -------------------------------------------------- */

int igraph_read_graph_edgelist(igraph_t *graph, FILE *instream, 
			       integer_t n, bool_t directed);
int igraph_read_graph_ncol(igraph_t *graph, FILE *instream, bool_t names, 
			  bool_t weights);
int igraph_read_graph_lgl(igraph_t *graph, FILE *instream,
			  bool_t names, bool_t weights);

int igraph_write_graph_edgelist(igraph_t *graph, FILE *outstream);
int igraph_write_graph_ncol(igraph_t *graph, FILE *outstream,
			    const char *names, const char *weights);
int igraph_write_graph_lgl(igraph_t *graph, FILE *outstream,
			   const char *names, const char *weights,
			   bool_t isolates);

/* -------------------------------------------------- */
/* Dynamics measurement                               */
/* -------------------------------------------------- */

int igraph_measure_dynamics_idage(igraph_t *graph, matrix_t *akl, 
				  matrix_t *sd,
				  vector_t *st, integer_t agebins,
				  integer_t maxind, bool_t lsd);
int igraph_measure_dynamics_idage_st(igraph_t *graph, vector_t *res,
				     matrix_t *akl);
int igraph_measure_dynamics_idage_debug(igraph_t *graph, matrix_t *akl,
					matrix_t *sd,
					vector_t *st, integer_t pagebins,
					integer_t pmaxind, bool_t lsd,
					vector_t *estimates, 
					integer_t est_ind, integer_t est_age);

/* -------------------------------------------------- */
/* Other, not graph related                           */
/* -------------------------------------------------- */

int igraph_running_mean(vector_t *data, vector_t *res, integer_t binwidth);
int igraph_random_sample(vector_t *res, integer_t l, integer_t h, 
			 integer_t length);

/* -------------------------------------------------- */
/* For internal use only, should move to other header */
/* -------------------------------------------------- */

typedef struct igraph_i_adjlist_t { 
  integer_t length;
  vector_t *adjs;
} igraph_i_adjlist_t;

int igraph_i_adjlist_init(igraph_t *graph, igraph_i_adjlist_t *al, 
			  igraph_neimode_t mode);
void igraph_i_adjlist_destroy(igraph_i_adjlist_t *al);
/* vector_t *igraph_i_adjlist_get(igraph_i_adjlist_t *al, integer_t no); */
#define igraph_i_adjlist_get(al, no) (&(al)->adjs[(long int)(no)])

__END_DECLS
  
#endif
