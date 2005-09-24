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

#include "types.h"

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
 * \defgroup iterators Vertex and edge iterators
 * \brief Iterators provide a method to walk through some edges or
 * vertices of the graph.
 * 
 * This documentation is not yet written, as it is possible that the concept 
 * of iterators will be reinterpreted soon.
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
} igraph_t;

/* -------------------------------------------------- */
/* Constants                                          */
/* -------------------------------------------------- */

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

/* -------------------------------------------------- */
/* Interface                                          */
/* -------------------------------------------------- */

int igraph_empty(igraph_t *graph, integer_t n, bool_t directed);
int igraph_destroy(igraph_t *graph);
int igraph_copy(igraph_t *to, igraph_t *from);
int igraph_add_edges(igraph_t *graph, vector_t *edges);
int igraph_add_vertices(igraph_t *graph, integer_t nv);
int igraph_delete_edges(igraph_t *graph, vector_t *edges);
int igraph_delete_vertices(igraph_t *graph, vector_t *vertices);
integer_t igraph_vcount(igraph_t *graph);
integer_t igraph_ecount(igraph_t *graph);
int igraph_neighbors(igraph_t *graph, vector_t *neis, integer_t vid, 
		     igraph_neimode_t mode); 
bool_t igraph_is_directed(igraph_t *graph);
int igraph_degree(igraph_t *graph, vector_t *res, vector_t *vids, 
		  igraph_neimode_t mode, bool_t loops);

/* -------------------------------------------------- */
/* Iterators                                          */
/* -------------------------------------------------- */

struct igraph_iterator_t;

typedef struct igraph_iterator_t {
  integer_t type;
  void *data;
  int (*next)(igraph_t *graph, struct igraph_iterator_t *it);
  int (*prev)(igraph_t *graph, struct igraph_iterator_t *it);
  bool_t (*end)(igraph_t *graph, struct igraph_iterator_t *it);
  integer_t (*getvertex)(igraph_t *graph, struct igraph_iterator_t *it);
  integer_t (*getvertexfrom)(igraph_t *graph, struct igraph_iterator_t *it);
  integer_t (*getvertexto)(igraph_t *graph, struct igraph_iterator_t *it);
  integer_t (*getvertexnei)(igraph_t *graph, struct igraph_iterator_t *it);
  integer_t (*getedge)(igraph_t *graph, struct igraph_iterator_t *it);
} igraph_iterator_t;

/* The constructors & destructor */
int igraph_iterator_vid(igraph_t *graph, igraph_iterator_t *it);
int igraph_iterator_vneis(igraph_t *graph, igraph_iterator_t *it, 
			  integer_t vid, igraph_neimode_t mode);
int igraph_iterator_eid(igraph_t *graph, igraph_iterator_t *it);
int igraph_iterator_efromorder(igraph_t *graph, igraph_iterator_t *it);
int igraph_iterator_eneis(igraph_t *graph, igraph_iterator_t *it, 
			  integer_t vid, igraph_neimode_t mode);

int igraph_iterator_destroy(igraph_t *graph, igraph_iterator_t *it);

/* Generics, common functions */
int igraph_next(igraph_t *graph, igraph_iterator_t *it);
int igraph_prev(igraph_t *graph, igraph_iterator_t *it);
bool_t igraph_end(igraph_t *graph, igraph_iterator_t *it);
integer_t igraph_get_vertex_nei(igraph_t *graph, igraph_iterator_t *it);

/* Generics, vertices */
integer_t igraph_get_vertex(igraph_t *graph, igraph_iterator_t *it);

/* Generics, edges */
integer_t igraph_get_vertex_from(igraph_t *graph, igraph_iterator_t *it);
integer_t igraph_get_vertex_to(igraph_t *graph, igraph_iterator_t *it);
integer_t igraph_get_edge(igraph_t *graph, igraph_iterator_t *it);

/* Specifics, simple vertex iterator */
int igraph_next_vid(igraph_t *graph, igraph_iterator_t *it);
int igraph_prev_vid(igraph_t *graph, igraph_iterator_t *it);
bool_t igraph_end_vid(igraph_t *graph, igraph_iterator_t *it);
integer_t igraph_get_vertex_vid(igraph_t *graph, igraph_iterator_t *it);

/* Iterates over the neighbors of a vertex */
int igraph_next_vneis(igraph_t *graph, igraph_iterator_t *it);
int igraph_end_vneis(igraph_t *graph, igraph_iterator_t *it);
integer_t igraph_get_vertex_vneis(igraph_t *graph, igraph_iterator_t *it);

/* Specifics, simple edge iterator */
int igraph_next_eid(igraph_t *graph, igraph_iterator_t *it);
int igraph_prev_eid(igraph_t *graph, igraph_iterator_t *it);
bool_t igraph_end_eid(igraph_t *graph, igraph_iterator_t *it);
integer_t igraph_get_vertex_from_eid(igraph_t *graph, igraph_iterator_t *it);
integer_t igraph_get_vertex_to_eid(igraph_t *graph, igraph_iterator_t *it);
integer_t igraph_get_edge_eid(igraph_t *graph, igraph_iterator_t *it); 

/* Specifics, edge iterator according to the 'from' order */
int igraph_next_efromorder(igraph_t *graph, igraph_iterator_t *it);
int igraph_prev_efromorder(igraph_t *graph, igraph_iterator_t *it);
bool_t igraph_end_efromorder(igraph_t *graph, igraph_iterator_t *it);
integer_t igraph_get_vertex_from_efromorder(igraph_t *graph, 
					    igraph_iterator_t *it);
integer_t igraph_get_vertex_to_efromorder(igraph_t *graph, 
					  igraph_iterator_t *it);
integer_t igraph_get_edge_efromorder(igraph_t *graph, igraph_iterator_t *it);


/* Iterates over the edges to and/or from a vertex */
int igraph_next_eneis(igraph_t *graph, igraph_iterator_t *it);
int igraph_prev_eneis(igraph_t *graph, igraph_iterator_t *it);
bool_t igraph_end_eneis(igraph_t *graph, igraph_iterator_t *it);
integer_t igraph_get_vertex_from_eneis(igraph_t *graph, igraph_iterator_t *it);
integer_t igraph_get_vertex_to_eneis(igraph_t *graph, igraph_iterator_t *it);
integer_t igraph_get_edge_eneis(igraph_t *graph, igraph_iterator_t *it); 
int igraph_iterator_eneis_set(igraph_t *graph, igraph_iterator_t *it, 
			      integer_t vid, igraph_neimode_t mode);
integer_t igraph_get_vertex_nei_eneis(igraph_t *graph, igraph_iterator_t *it);

/* TODO: attributes */

/* -------------------------------------------------- */
/* Error handling                                     */
/* -------------------------------------------------- */

int igraph_error(const char *msg);

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
				vector_t *in_deg, integer_t method);
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
int igraph_closeness(igraph_t *graph, vector_t *res, vector_t *vids, 
		     igraph_neimode_t mode);
int igraph_shortest_paths(igraph_t *graph, matrix_t *res, 
			  vector_t *from, igraph_neimode_t mode);
int igraph_get_shortest_paths(igraph_t *graph, vector_t *res,
			      integer_t from, igraph_neimode_t mode);
int igraph_subcomponent(igraph_t *graph, vector_t *res, real_t vid, 
			igraph_neimode_t mode);
int igraph_betweenness (igraph_t *graph, vector_t *res, vector_t *vids, 
			bool_t directed);
int igraph_edge_betweenness (igraph_t *graph, vector_t *result, 
			     bool_t directed);
int igraph_subgraph(igraph_t *graph, igraph_t *res, vector_t *vids);
int igraph_average_path_length(igraph_t *graph, real_t *res,
			       bool_t directed, bool_t unconn);

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
				       integer_t niter, real_t coolexp,
				       integer_t frame, vector_t *initial,
				       real_t initemp);
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

int igraph_cocitation(igraph_t *graph, matrix_t *res, vector_t *vids);
int igraph_bibcoupling(igraph_t *graph, matrix_t *res, vector_t *vids);

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

int igraph_get_adjacency(igraph_t *graph, vector_t *res);
int igraph_get_edgelist(igraph_t *graph, vector_t *res);

/* -------------------------------------------------- */
/* Read and write foreign formats                     */
/* -------------------------------------------------- */

int igraph_read_graph(const char *filename, integer_t format);
int igraph_write_graph(igraph_t *graph, const char *filename, 
		       integer_t format);

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

__END_DECLS
  
#endif
