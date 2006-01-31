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
  igraph_vector_t from;
  igraph_vector_t to;
  igraph_vector_t oi;
  igraph_vector_t ii;
  igraph_vector_t os;
  igraph_vector_t is;
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
#define IGRAPH_ITERATOR_ES_SEQ          6

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
  void (*next)(const igraph_t *graph, struct igraph_vs_t *vs);
  bool_t (*end)(const igraph_t *graph, const struct igraph_vs_t *vs);
  void (*reset)(const igraph_t *graph, struct igraph_vs_t *vs);
  integer_t(*get)(const igraph_t *graph, const struct igraph_vs_t *vs);
  int (*unfold)(const igraph_t *graph, const struct igraph_vs_t *vs,
		igraph_vector_t *v);
  void (*destroy)(struct igraph_vs_t *vs);
} igraph_i_vstable_t;

typedef struct igraph_vs_t {
  int type;
  real_t stdata[IGRAPH_I_STDATA_SIZE]; /* working storage */
  void *pdata;			       /* additional storage */
  igraph_i_vstable_t *table;	       /* the table of functions */
  bool_t shorthand;
} igraph_vs_t;

/* -------------------------------------------------- */
/* Edge sequences                                   */
/* -------------------------------------------------- */

struct igraph_es_t;

typedef struct igraph_i_estable_t {
  void (*next)(const igraph_t *graph, struct igraph_es_t *es);
  bool_t (*end)(const igraph_t *graph, const struct igraph_es_t *es);
  void (*reset)(const igraph_t *graph, struct igraph_es_t *es);
  integer_t(*get)(const igraph_t *graph, const struct igraph_es_t *es);
  integer_t (*from)(const igraph_t *graph, const struct igraph_es_t *es);
  integer_t (*to)(const igraph_t *graph, const struct igraph_es_t *es);
  int (*unfold)(const igraph_t *graph, const struct igraph_es_t *es,
		igraph_vector_t *v);
  void (*destroy)(struct igraph_es_t *es);
} igraph_i_estable_t;

typedef struct igraph_es_t {
  int type;
  real_t stdata[IGRAPH_I_STDATA_SIZE]; /* working storage */
  void *pdata;			   /* other storage   */
  igraph_i_estable_t *table;	   /* the table of functions */
  bool_t shorthand;
} igraph_es_t;

/* -------------------------------------------------- */
/* Interface                                          */
/* -------------------------------------------------- */

struct igraph_eit_t;

int igraph_empty(igraph_t *graph, integer_t n, bool_t directed);
int igraph_destroy(igraph_t *graph);
int igraph_copy(igraph_t *to, const igraph_t *from);
int igraph_add_edges(igraph_t *graph, const igraph_vector_t *edges);
int igraph_add_vertices(igraph_t *graph, integer_t nv);
int igraph_delete_edges(igraph_t *graph, const igraph_vector_t *edges); /*eee*/
int igraph_delete_vertices(igraph_t *graph, const igraph_vs_t *vertices);
integer_t igraph_vcount(const igraph_t *graph);
integer_t igraph_ecount(const igraph_t *graph);
int igraph_neighbors(const igraph_t *graph, igraph_vector_t *neis, integer_t vid, 
		     igraph_neimode_t mode); 
bool_t igraph_is_directed(const igraph_t *graph);
int igraph_degree(const igraph_t *graph, igraph_vector_t *res, 
		  const igraph_vs_t *vids, igraph_neimode_t mode, 
		  bool_t loops);
int igraph_edge(const igraph_t *graph, integer_t eid, 
		integer_t *from, integer_t *to);		
int igraph_adjacent(const igraph_t *graph, igraph_vector_t *eids, integer_t vid,
		    igraph_neimode_t mode);

/* -------------------------------------------------- */
/* Vertex sequences, contd.                           */
/* -------------------------------------------------- */

/* Generics */
void igraph_vs_next(const igraph_t *graph, igraph_vs_t *vs);
bool_t igraph_vs_end(const igraph_t *graph, const igraph_vs_t *vs);
void igraph_vs_reset(const igraph_t *graph, igraph_vs_t *vs);
integer_t igraph_vs_get(const igraph_t *graph, const igraph_vs_t *vs);
int igraph_vs_unfold(const igraph_t *graph, const igraph_vs_t *vs,
		     igraph_vector_t *v);
void igraph_vs_destroy(igraph_vs_t *vs);

/* Simple vertex iterator, all vertices */
int igraph_vs_all(const igraph_t *graph, igraph_vs_t *it);
const igraph_vs_t *IGRAPH_VS_ALL(const igraph_t *graph);

/* Adjacent vertices of a vertex */
int igraph_vs_adj(const igraph_t *graph, igraph_vs_t *it,
		  integer_t vid, igraph_neimode_t mode);
void igraph_vs_adj_set(const igraph_t *graph, igraph_vs_t *vs,
		       integer_t vid, igraph_neimode_t mode);

/* Random walker */
int igraph_vs_rw(const igraph_t *graph, igraph_vs_t *it,
		 integer_t vid, igraph_neimode_t mode);
long int igraph_vs_rw_length(const igraph_t *graph, const igraph_vs_t *vs);

/* Random walker with one unit memory */
int igraph_vs_rw1(const igraph_t *graph, igraph_vs_t *it,
		  integer_t vid, igraph_neimode_t mode);
long int igraph_vs_rw1_length(const igraph_t *graph, const igraph_vs_t *vs);

/* Empty iterator */
int igraph_vs_none(const igraph_t *graph, igraph_vs_t *it);

/* Single vertex iterator */
int igraph_vs_1(const igraph_t *igraph, igraph_vs_t *it, integer_t vid);
const igraph_vs_t *IGRAPH_VS_1(const igraph_t *graph, integer_t vid);

/* A sequence of vertices */
int igraph_vs_seq(const igraph_t *igraph, igraph_vs_t *it, integer_t from, 
		  integer_t to);

/* Iterates over a vector */
int igraph_vs_vectorview(const igraph_t *igraph, igraph_vs_t *it, 
			 const igraph_vector_t *vids);
int igraph_vs_vectorview_it(const igraph_t *graph, const igraph_vs_t *vs,
			    igraph_vs_t *newvs);
int igraph_vs_vector(const igraph_t *igraph, igraph_vs_t *it, 
		     const igraph_vector_t *vids);
int igraph_vs_vector_small(const igraph_t *igraph, igraph_vs_t *it, ...);
const igraph_vs_t *IGRAPH_VS_VECTOR(const igraph_t *graph, 
				    const igraph_vector_t *vids);
const igraph_vs_t *IGRAPH_VS(const igraph_t *graph, ...);
const igraph_vector_t *igraph_vs_vector_getvector(const igraph_t *graph, 
					   const igraph_vs_t *vs);

/* -------------------------------------------------- */
/* Edge sequences contd.                              */
/* -------------------------------------------------- */

/* Generics */
void igraph_es_next(const igraph_t *graph, igraph_es_t *es);
bool_t igraph_es_end(const igraph_t *graph, const igraph_es_t *es);
void igraph_es_reset(const igraph_t *graph, igraph_es_t *es);
integer_t igraph_es_get(const igraph_t *graph, const igraph_es_t *es);
integer_t igraph_es_from(const igraph_t *graph, const igraph_es_t *es);
integer_t igraph_es_to(const igraph_t *graph, const igraph_es_t *es);
int igraph_es_unfold(const igraph_t *graph, const igraph_es_t *es, 
		     igraph_vector_t *v);
void igraph_es_destroy(igraph_es_t *es);

/* Simple edge iterator */
int igraph_es_all(const igraph_t *graph, igraph_es_t *it);
const igraph_es_t *IGRAPH_ES_ALL(const igraph_t *graph);

/* Empty edge iterator */
int igraph_es_none(const igraph_t *graph, igraph_es_t *it);

/* Single edge edge iterator */
int igraph_es_1(const igraph_t *graph, igraph_es_t *it, integer_t id);
const igraph_es_t *IGRAPH_ES_1(const igraph_t *graph, integer_t eid);

/* A sequence of vertices */
int igraph_es_seq(const igraph_t *igraph, igraph_es_t *it, integer_t from, 
		  integer_t to);

/* Sorted edge iterator */
int igraph_es_fromorder(const igraph_t *graph, igraph_es_t *it);

/* Adjacent edges of a vertex */
int igraph_es_adj(const igraph_t *graph, igraph_es_t *it,
		  integer_t vid, igraph_neimode_t mode);
void igraph_es_adj_set(const igraph_t *graph, igraph_es_t *es,
			integer_t vid, igraph_neimode_t mode);
integer_t igraph_es_adj_vertex(const igraph_t *graph, const igraph_es_t *es);

/* View of a vector */
int igraph_es_vectorview(const igraph_t *graph, igraph_es_t *it, 
			 const igraph_vector_t *eids);
int igraph_es_vectorview_it(const igraph_t *graph, const igraph_es_t *es,
			    igraph_es_t *newes);
int igraph_es_vector(const igraph_t *igraph, igraph_es_t *it, 
		     const igraph_vector_t *eids);
int igraph_es_vector_small(const igraph_t *graph, igraph_es_t *es, ...);
int igraph_es_fromto(const igraph_t *igraph, igraph_es_t *it, 
		     const igraph_vs_t *from,
		     const igraph_vs_t *to, bool_t directed);
const igraph_es_t *IGRAPH_ES_VECTOR(const igraph_t *graph, 
				    const igraph_vector_t *eids);
const igraph_es_t *IGRAPH_ES(const igraph_t *graph, ...);
const igraph_vector_t *igraph_es_vector_getvector(const igraph_t *graph, 
					   const igraph_es_t *vs);

/* -------------------------------------------------- */
/* Attributes                                         */
/* -------------------------------------------------- */

int igraph_add_graph_attribute(igraph_t *graph, const char *name,
			       igraph_attribute_type_t type);
int igraph_remove_graph_attribute(igraph_t *graph, const char *name);
int igraph_get_graph_attribute(const igraph_t *graph, const char *name, 
			       void **value,
			       igraph_attribute_type_t *type);
int igraph_set_graph_attribute(igraph_t *graph, const char *name, 
			       const void *value);
int igraph_list_graph_attributes(const igraph_t *graph, igraph_strvector_t *l,
				 igraph_vector_t *types);
int igraph_get_graph_attribute_type(const igraph_t *graph, const char *name, 
				    igraph_attribute_type_t *type);
bool_t igraph_has_graph_attribute(const igraph_t *graph, const char *name);

int igraph_add_vertex_attribute(igraph_t *graph, const char *name,
				igraph_attribute_type_t type);
int igraph_remove_vertex_attribute(igraph_t *graph, const char *name);
int igraph_get_vertex_attribute(const igraph_t *graph, const char *name, 
				long int v, void **value,
				igraph_attribute_type_t *type);
int igraph_set_vertex_attribute(igraph_t *graph, const char *name, 
				long int v, const void *value);
int igraph_get_vertex_attributes(const igraph_t *graph, const char *name, 
				 const igraph_vs_t *v, void **value);
int igraph_set_vertex_attributes(igraph_t *graph, const char *name, 
				 const igraph_vs_t *v, const void *value);
int igraph_list_vertex_attributes(const igraph_t *graph, igraph_strvector_t *l,
				  igraph_vector_t *types);
int igraph_get_vertex_attribute_type(const igraph_t *graph, const char *name, 
				     igraph_attribute_type_t *type);
bool_t igraph_has_vertex_attribute(const igraph_t *graph, const char *name);

int igraph_add_edge_attribute(igraph_t *graph, const char *name, 
			      igraph_attribute_type_t type);
int igraph_remove_edge_attribute(igraph_t *graph, const char *name);
int igraph_get_edge_attribute(const igraph_t *graph, const char *name, 
			      long int e, void **value,
			      igraph_attribute_type_t *type);
int igraph_set_edge_attribute(igraph_t *graph, const char *name, 
			      long int e, const void *value);
int igraph_get_edge_attributes(const igraph_t *graph, const char *name, 
			       const igraph_es_t *e, void **value);
int igraph_set_edge_attributes(igraph_t *graph, const char *name, 
			       const igraph_es_t *e, const void *value); 
int igraph_list_edge_attributes(const igraph_t *graph, igraph_strvector_t *l,
				igraph_vector_t *types);
int igraph_get_edge_attribute_type(const igraph_t *graph, const char *name, 
				   igraph_attribute_type_t *type);
bool_t igraph_has_edge_attribute(const igraph_t *graph, const char *name);

/* -------------------------------------------------- */
/* Constructors, deterministic                        */
/* -------------------------------------------------- */

int igraph_create(igraph_t *graph, const igraph_vector_t *edges, integer_t n, 
		  bool_t directed);
int igraph_adjacency(igraph_t *graph, igraph_matrix_t *adjmatrix,
		     igraph_adjacency_t mode);
int igraph_star(igraph_t *graph, integer_t n, igraph_star_mode_t mode, 
		integer_t center);
int igraph_lattice(igraph_t *graph, const igraph_vector_t *dimvector, integer_t nei, 
		   bool_t directed, bool_t mutual, bool_t circular);
int igraph_ring(igraph_t *graph, integer_t n, bool_t directed, 
		bool_t mutual, bool_t circular);
int igraph_tree(igraph_t *graph, integer_t n, integer_t children, 
		igraph_tree_mode_t type);
int igraph_full(igraph_t *graph, integer_t n, bool_t directed, bool_t loops);
int igraph_atlas(igraph_t *graph, int number);

/* -------------------------------------------------- */
/* Constructors, games (=stochastic)                  */
/* -------------------------------------------------- */

int igraph_barabasi_game(igraph_t *graph, integer_t n, integer_t m, 
			 const igraph_vector_t *outseq, bool_t outpref, 
			 bool_t directed);
int igraph_erdos_renyi_game(igraph_t *graph, igraph_erdos_renyi_t type,
			    integer_t n, real_t p,
			    bool_t directed, bool_t loops);
int igraph_erdos_renyi_game_gnp(igraph_t *graph, integer_t n, real_t p,
				bool_t directed, bool_t loops);
int igraph_degree_sequence_game(igraph_t *graph, const igraph_vector_t *out_deg,
				const igraph_vector_t *in_deg, 
				igraph_degseq_t method);
int igraph_growing_random_game(igraph_t *graph, integer_t n, 
			       integer_t m, bool_t directed, bool_t citation);
int igraph_aging_prefatt_game(igraph_t *graph, integer_t n, integer_t m,
			      integer_t aging_type, real_t aging_exp);

/* -------------------------------------------------- */
/* Basic query functions                              */
/* -------------------------------------------------- */

bool_t igraph_are_connected(const igraph_t *graph, integer_t v1, integer_t v2);

/* -------------------------------------------------- */
/* Structural properties                              */
/* -------------------------------------------------- */

int igraph_diameter(const igraph_t *graph, integer_t *res, 
		    bool_t directed, bool_t unconn);
int igraph_minimum_spanning_tree_unweighted(const igraph_t *graph, 
					    igraph_t *mst);
int igraph_minimum_spanning_tree_prim(const igraph_t *graph, igraph_t *mst,
				      const igraph_vector_t *weights);
int igraph_closeness(const igraph_t *graph, igraph_vector_t *res, 
		     const igraph_vs_t *vids, igraph_neimode_t mode);
int igraph_shortest_paths(const igraph_t *graph, igraph_matrix_t *res, 
			  const igraph_vs_t *from, igraph_neimode_t mode);
int igraph_get_shortest_paths(const igraph_t *graph, igraph_vector_t *res,
			      integer_t from, igraph_neimode_t mode);
int igraph_subcomponent(const igraph_t *graph, igraph_vector_t *res, real_t vid, 
			igraph_neimode_t mode);	
int igraph_betweenness (const igraph_t *graph, igraph_vector_t *res, 
			const igraph_vs_t *vids, bool_t directed);
int igraph_edge_betweenness (const igraph_t *graph, igraph_vector_t *result,
			     bool_t directed); /* eee + add */
int igraph_pagerank(const igraph_t *graph, igraph_vector_t *res, 
		    const igraph_vs_t *vids, bool_t directed, integer_t niter, 
		    real_t eps, real_t damping);
int igraph_subgraph(const igraph_t *graph, igraph_t *res, 
		    const igraph_vs_t *vids);
int igraph_average_path_length(const igraph_t *graph, real_t *res,
			       bool_t directed, bool_t unconn);
int igraph_simplify(igraph_t *graph, bool_t multiple, bool_t loops);
int igraph_transitivity(const igraph_t *graph, igraph_vector_t *res, 
			igraph_transitivity_type_t type); /* vvv + add */

/* TODO: degree.distribution (?) */

/* -------------------------------------------------- */
/* Components                                         */
/* -------------------------------------------------- */

int igraph_clusters(const igraph_t *graph, igraph_vector_t *membership, 
		    igraph_vector_t *csize, igraph_connectedness_t mode);
int igraph_is_connected(const igraph_t *graph, bool_t *res, 
			igraph_connectedness_t mode);
int igraph_decompose(const igraph_t *graph, igraph_vector_ptr_t *components, 
		     igraph_connectedness_t mode, 
		     long int maxcompno, long int minelements);

/* TODO: cluster.distribution (?) */

/* -------------------------------------------------- */
/* Layouts                                            */
/* -------------------------------------------------- */

int igraph_layout_random(const igraph_t *graph, igraph_matrix_t *res);
int igraph_layout_circle(const igraph_t *graph, igraph_matrix_t *res);
int igraph_layout_fruchterman_reingold(const igraph_t *graph, igraph_matrix_t *res,
				       integer_t niter, real_t maxdelta,
				       real_t area, real_t coolexp, 
				       real_t repulserad, bool_t use_seed);
int igraph_layout_kamada_kawai(const igraph_t *graph, igraph_matrix_t *res,
			       integer_t niter, real_t sigma, 
			       real_t initemp, real_t coolexp,
			       real_t kkconst);
int igraph_layout_springs(const igraph_t *graph, igraph_matrix_t *res,
			  real_t mass, real_t equil, real_t k,
			  real_t repeqdis, real_t kfr, bool_t repulse);
int igraph_layout_lgl(const igraph_t *graph, igraph_matrix_t *res,
		      integer_t maxiter, real_t maxdelta, 
		      real_t area, real_t coolexp,
		      real_t repulserad, real_t cellsize, integer_t root);

/* -------------------------------------------------- */
/* Visitor-like functions                             */
/* -------------------------------------------------- */

int igraph_bfs(igraph_t *graph, integer_t vid, igraph_neimode_t mode,
	       igraph_vector_t *vids, igraph_vector_t *layers,
	       igraph_vector_t *parents);

/* -------------------------------------------------- */
/* Centrality                                         */
/* -------------------------------------------------- */

/* TODO: evcent */

/* -------------------------------------------------- */
/* Cocitation                                         */
/* -------------------------------------------------- */

int igraph_cocitation(const igraph_t *graph, igraph_matrix_t *res, 
		      const igraph_vs_t *vids);
int igraph_bibcoupling(const igraph_t *graph, igraph_matrix_t *res, 
		       const igraph_vs_t *vids);

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

int igraph_get_adjacency(const igraph_t *graph, igraph_matrix_t *res,
			 igraph_get_adjacency_t type);
int igraph_get_edgelist(const igraph_t *graph, igraph_vector_t *res, bool_t bycol);

/* -------------------------------------------------- */
/* Read and write foreign formats                     */
/* -------------------------------------------------- */

int igraph_read_graph_edgelist(igraph_t *graph, FILE *instream, 
			       integer_t n, bool_t directed);
int igraph_read_graph_ncol(igraph_t *graph, FILE *instream, bool_t names, 
			  bool_t weights);
int igraph_read_graph_lgl(igraph_t *graph, FILE *instream,
			  bool_t names, bool_t weights);

int igraph_write_graph_edgelist(const igraph_t *graph, FILE *outstream);
int igraph_write_graph_ncol(const igraph_t *graph, FILE *outstream,
			    const char *names, const char *weights);
int igraph_write_graph_lgl(const igraph_t *graph, FILE *outstream,
			   const char *names, const char *weights,
			   bool_t isolates);

/* -------------------------------------------------- */
/* Dynamics measurement                               */
/* -------------------------------------------------- */

int igraph_measure_dynamics_idage(const igraph_t *graph, igraph_matrix_t *akl, 
				  igraph_matrix_t *sd,
				  const igraph_vector_t *st, integer_t agebins,
				  integer_t maxind, bool_t lsd);
int igraph_measure_dynamics_idage_st(const igraph_t *graph, igraph_vector_t *res,
				     const igraph_matrix_t *akl);
int igraph_measure_dynamics_idage_debug(const igraph_t *graph, igraph_matrix_t *akl,
					igraph_matrix_t *sd,
					const igraph_vector_t *st, integer_t pagebins,
					integer_t pmaxind, bool_t lsd,
					igraph_vector_t *estimates, 
					integer_t est_ind, integer_t est_age);

/* -------------------------------------------------- */
/* Other, not graph related                           */
/* -------------------------------------------------- */

int igraph_running_mean(const igraph_vector_t *data, igraph_vector_t *res, 
			integer_t binwidth);
int igraph_random_sample(igraph_vector_t *res, integer_t l, integer_t h, 
			 integer_t length);

/* -------------------------------------------------- */
/* For internal use only, should move to other header */
/* -------------------------------------------------- */

typedef struct igraph_i_adjlist_t { 
  integer_t length;
  igraph_vector_t *adjs;
} igraph_i_adjlist_t;

int igraph_i_adjlist_init(const igraph_t *graph, igraph_i_adjlist_t *al, 
			  igraph_neimode_t mode);
void igraph_i_adjlist_destroy(igraph_i_adjlist_t *al);
/* igraph_vector_t *igraph_i_adjlist_get(const igraph_i_adjlist_t *al,  */
/* 			       integer_t no); */
#define igraph_i_adjlist_get(al, no) (&(al)->adjs[(long int)(no)])

__END_DECLS
  
#endif
