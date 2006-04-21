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
   Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA 
   02110-1301 USA

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
#include "interrupt.h"

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
  igraph_integer_t n;
  igraph_bool_t directed;
  igraph_vector_t from;
  igraph_vector_t to;
  igraph_vector_t oi;
  igraph_vector_t ii;
  igraph_vector_t os;
  igraph_vector_t is;
  void *attr;
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
	       IGRAPH_FILEFORMAT_PAJEK,
               IGRAPH_FILEFORMAT_LGL,
               IGRAPH_FILEFORMAT_GRAPHML } igraph_fileformat_type_t;

typedef enum { IGRAPH_REWIRING_SIMPLE=0 } igraph_rewiring_t;

/* -------------------------------------------------- */
/* Vertex sequences                                   */
/* -------------------------------------------------- */

#define IGRAPH_I_STDATA_SIZE 10

struct igraph_vs_t;

typedef struct igraph_i_vstable_t {
  void (*next)(const igraph_t *graph, struct igraph_vs_t *vs);
  igraph_bool_t (*end)(const igraph_t *graph, const struct igraph_vs_t *vs);
  void (*reset)(const igraph_t *graph, struct igraph_vs_t *vs);
  igraph_integer_t(*get)(const igraph_t *graph, const struct igraph_vs_t *vs);
  int (*unfold)(const igraph_t *graph, const struct igraph_vs_t *vs,
		igraph_vector_t *v);
  void (*destroy)(struct igraph_vs_t *vs);
} igraph_i_vstable_t;

typedef struct igraph_vs_t {
  int type;
  igraph_real_t stdata[IGRAPH_I_STDATA_SIZE]; /* working storage */
  void *pdata;			       /* additional storage */
  igraph_i_vstable_t *table;	       /* the table of functions */
  igraph_bool_t shorthand;
} igraph_vs_t;

/* -------------------------------------------------- */
/* Edge sequences                                   */
/* -------------------------------------------------- */

struct igraph_es_t;

typedef struct igraph_i_estable_t {
  void (*next)(const igraph_t *graph, struct igraph_es_t *es);
  igraph_bool_t (*end)(const igraph_t *graph, const struct igraph_es_t *es);
  void (*reset)(const igraph_t *graph, struct igraph_es_t *es);
  igraph_integer_t(*get)(const igraph_t *graph, const struct igraph_es_t *es);
  igraph_integer_t (*from)(const igraph_t *graph, const struct igraph_es_t *es);
  igraph_integer_t (*to)(const igraph_t *graph, const struct igraph_es_t *es);
  int (*unfold)(const igraph_t *graph, const struct igraph_es_t *es,
		igraph_vector_t *v);
  void (*destroy)(struct igraph_es_t *es);
} igraph_i_estable_t;

typedef struct igraph_es_t {
  int type;
  igraph_real_t stdata[IGRAPH_I_STDATA_SIZE]; /* working storage */
  void *pdata;			   /* other storage   */
  igraph_i_estable_t *table;	   /* the table of functions */
  igraph_bool_t shorthand;
} igraph_es_t;

/* -------------------------------------------------- */
/* Interface                                          */
/* -------------------------------------------------- */

struct igraph_eit_t;

int igraph_empty(igraph_t *graph, igraph_integer_t n, igraph_bool_t directed);
int igraph_destroy(igraph_t *graph);
int igraph_copy(igraph_t *to, const igraph_t *from);
int igraph_add_edges(igraph_t *graph, const igraph_vector_t *edges);
int igraph_add_vertices(igraph_t *graph, igraph_integer_t nv);
int igraph_delete_edges(igraph_t *graph, const igraph_vector_t *edges); /*eee*/
int igraph_delete_vertices(igraph_t *graph, const igraph_vs_t *vertices);
igraph_integer_t igraph_vcount(const igraph_t *graph);
igraph_integer_t igraph_ecount(const igraph_t *graph);
int igraph_neighbors(const igraph_t *graph, igraph_vector_t *neis, igraph_integer_t vid, 
		     igraph_neimode_t mode); 
igraph_bool_t igraph_is_directed(const igraph_t *graph);
int igraph_degree(const igraph_t *graph, igraph_vector_t *res, 
		  const igraph_vs_t *vids, igraph_neimode_t mode, 
		  igraph_bool_t loops);
int igraph_edge(const igraph_t *graph, igraph_integer_t eid, 
		igraph_integer_t *from, igraph_integer_t *to);		
int igraph_adjacent(const igraph_t *graph, igraph_vector_t *eids, igraph_integer_t vid,
		    igraph_neimode_t mode);

/* -------------------------------------------------- */
/* Vertex sequences, contd.                           */
/* -------------------------------------------------- */

/* Generics */
void igraph_vs_next(const igraph_t *graph, igraph_vs_t *vs);
igraph_bool_t igraph_vs_end(const igraph_t *graph, const igraph_vs_t *vs);
void igraph_vs_reset(const igraph_t *graph, igraph_vs_t *vs);
igraph_integer_t igraph_vs_get(const igraph_t *graph, const igraph_vs_t *vs);
int igraph_vs_unfold(const igraph_t *graph, const igraph_vs_t *vs,
		     igraph_vector_t *v);
void igraph_vs_destroy(igraph_vs_t *vs);

/* Simple vertex iterator, all vertices */
int igraph_vs_all(const igraph_t *graph, igraph_vs_t *it);
const igraph_vs_t *IGRAPH_VS_ALL(const igraph_t *graph);

/* Adjacent vertices of a vertex */
int igraph_vs_adj(const igraph_t *graph, igraph_vs_t *it,
		  igraph_integer_t vid, igraph_neimode_t mode);
void igraph_vs_adj_set(const igraph_t *graph, igraph_vs_t *vs,
		       igraph_integer_t vid, igraph_neimode_t mode);

/* Random walker */
int igraph_vs_rw(const igraph_t *graph, igraph_vs_t *it,
		 igraph_integer_t vid, igraph_neimode_t mode);
long int igraph_vs_rw_length(const igraph_t *graph, const igraph_vs_t *vs);

/* Random walker with one unit memory */
int igraph_vs_rw1(const igraph_t *graph, igraph_vs_t *it,
		  igraph_integer_t vid, igraph_neimode_t mode);
long int igraph_vs_rw1_length(const igraph_t *graph, const igraph_vs_t *vs);

/* Empty iterator */
int igraph_vs_none(const igraph_t *graph, igraph_vs_t *it);

/* Single vertex iterator */
int igraph_vs_1(const igraph_t *igraph, igraph_vs_t *it, igraph_integer_t vid);
const igraph_vs_t *IGRAPH_VS_1(const igraph_t *graph, igraph_integer_t vid);

/* A sequence of vertices */
int igraph_vs_seq(const igraph_t *igraph, igraph_vs_t *it, igraph_integer_t from, 
		  igraph_integer_t to);

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
igraph_bool_t igraph_es_end(const igraph_t *graph, const igraph_es_t *es);
void igraph_es_reset(const igraph_t *graph, igraph_es_t *es);
igraph_integer_t igraph_es_get(const igraph_t *graph, const igraph_es_t *es);
igraph_integer_t igraph_es_from(const igraph_t *graph, const igraph_es_t *es);
igraph_integer_t igraph_es_to(const igraph_t *graph, const igraph_es_t *es);
int igraph_es_unfold(const igraph_t *graph, const igraph_es_t *es, 
		     igraph_vector_t *v);
void igraph_es_destroy(igraph_es_t *es);

/* Simple edge iterator */
int igraph_es_all(const igraph_t *graph, igraph_es_t *it);
const igraph_es_t *IGRAPH_ES_ALL(const igraph_t *graph);

/* Empty edge iterator */
int igraph_es_none(const igraph_t *graph, igraph_es_t *it);

/* Single edge edge iterator */
int igraph_es_1(const igraph_t *graph, igraph_es_t *it, igraph_integer_t id);
const igraph_es_t *IGRAPH_ES_1(const igraph_t *graph, igraph_integer_t eid);

/* A sequence of vertices */
int igraph_es_seq(const igraph_t *igraph, igraph_es_t *it, igraph_integer_t from, 
		  igraph_integer_t to);

/* Sorted edge iterator */
int igraph_es_fromorder(const igraph_t *graph, igraph_es_t *it);

/* Adjacent edges of a vertex */
int igraph_es_adj(const igraph_t *graph, igraph_es_t *it,
		  igraph_integer_t vid, igraph_neimode_t mode);
void igraph_es_adj_set(const igraph_t *graph, igraph_es_t *es,
			igraph_integer_t vid, igraph_neimode_t mode);
igraph_integer_t igraph_es_adj_vertex(const igraph_t *graph, const igraph_es_t *es);

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
		     const igraph_vs_t *to, igraph_bool_t directed);
const igraph_es_t *IGRAPH_ES_VECTOR(const igraph_t *graph, 
				    const igraph_vector_t *eids);
const igraph_es_t *IGRAPH_ES(const igraph_t *graph, ...);
const igraph_vector_t *igraph_es_vector_getvector(const igraph_t *graph, 
					   const igraph_es_t *vs);

/* -------------------------------------------------- */
/* Constructors, deterministic                        */
/* -------------------------------------------------- */

int igraph_create(igraph_t *graph, const igraph_vector_t *edges, igraph_integer_t n, 
		  igraph_bool_t directed);
int igraph_adjacency(igraph_t *graph, igraph_matrix_t *adjmatrix,
		     igraph_adjacency_t mode);
int igraph_star(igraph_t *graph, igraph_integer_t n, igraph_star_mode_t mode, 
		igraph_integer_t center);
int igraph_lattice(igraph_t *graph, const igraph_vector_t *dimvector, igraph_integer_t nei, 
		   igraph_bool_t directed, igraph_bool_t mutual, igraph_bool_t circular);
int igraph_ring(igraph_t *graph, igraph_integer_t n, igraph_bool_t directed, 
		igraph_bool_t mutual, igraph_bool_t circular);
int igraph_tree(igraph_t *graph, igraph_integer_t n, igraph_integer_t children, 
		igraph_tree_mode_t type);
int igraph_full(igraph_t *graph, igraph_integer_t n, igraph_bool_t directed, igraph_bool_t loops);
int igraph_atlas(igraph_t *graph, int number);

/* -------------------------------------------------- */
/* Constructors, games (=stochastic)                  */
/* -------------------------------------------------- */

int igraph_barabasi_game(igraph_t *graph, igraph_integer_t n, igraph_integer_t m, 
			 const igraph_vector_t *outseq, igraph_bool_t outpref, 
			 igraph_bool_t directed);
int igraph_erdos_renyi_game(igraph_t *graph, igraph_erdos_renyi_t type,
			    igraph_integer_t n, igraph_real_t p,
			    igraph_bool_t directed, igraph_bool_t loops);
int igraph_erdos_renyi_game_gnp(igraph_t *graph, igraph_integer_t n, igraph_real_t p,
				igraph_bool_t directed, igraph_bool_t loops);
int igraph_degree_sequence_game(igraph_t *graph, const igraph_vector_t *out_deg,
				const igraph_vector_t *in_deg, 
				igraph_degseq_t method);
int igraph_growing_random_game(igraph_t *graph, igraph_integer_t n, 
			       igraph_integer_t m, igraph_bool_t directed, igraph_bool_t citation);
int igraph_aging_prefatt_game(igraph_t *graph, igraph_integer_t n, igraph_integer_t m,
			      igraph_integer_t aging_type, igraph_real_t aging_exp);

int igraph_callaway_traits_game (igraph_t *graph, igraph_integer_t nodes, 
				 igraph_integer_t types, igraph_integer_t edges_per_step, 
				 igraph_vector_t *type_dist,
				 igraph_matrix_t *pref_matrix,
				 igraph_bool_t directed);
int igraph_establishment_game(igraph_t *graph, igraph_integer_t nodes,
			      igraph_integer_t types, igraph_integer_t k,
			      igraph_vector_t *type_dist,
			      igraph_matrix_t *pref_matrix,
			      igraph_bool_t directed);

/* -------------------------------------------------- */
/* Basic query functions                              */
/* -------------------------------------------------- */

igraph_bool_t igraph_are_connected(const igraph_t *graph, igraph_integer_t v1, igraph_integer_t v2);

/* -------------------------------------------------- */
/* Structural properties                              */
/* -------------------------------------------------- */

int igraph_diameter(const igraph_t *graph, igraph_integer_t *res, 
		    igraph_integer_t *from, igraph_integer_t *to,
		    igraph_vector_t *path,
		    igraph_bool_t directed, igraph_bool_t unconn);
int igraph_minimum_spanning_tree_unweighted(const igraph_t *graph, 
					    igraph_t *mst);
int igraph_minimum_spanning_tree_prim(const igraph_t *graph, igraph_t *mst,
				      const igraph_vector_t *weights);
int igraph_closeness(const igraph_t *graph, igraph_vector_t *res, 
		     const igraph_vs_t *vids, igraph_neimode_t mode);
int igraph_shortest_paths(const igraph_t *graph, igraph_matrix_t *res, 
			  const igraph_vs_t *from, igraph_neimode_t mode);
int igraph_get_shortest_paths(const igraph_t *graph, igraph_vector_ptr_t *res,
			      igraph_integer_t from, const igraph_vs_t *to, 
			      igraph_neimode_t mode);
int igraph_get_all_shortest_paths(const igraph_t *graph,
				  igraph_vector_ptr_t *res, 
				  igraph_vector_t *nrgeo,
				  igraph_integer_t from, const igraph_vs_t *to, 
				  igraph_neimode_t mode);
int igraph_subcomponent(const igraph_t *graph, igraph_vector_t *res, igraph_real_t vid, 
			igraph_neimode_t mode);	
int igraph_betweenness (const igraph_t *graph, igraph_vector_t *res, 
			const igraph_vs_t *vids, igraph_bool_t directed);
int igraph_edge_betweenness (const igraph_t *graph, igraph_vector_t *result,
			     igraph_bool_t directed); /* eee + add */
int igraph_pagerank(const igraph_t *graph, igraph_vector_t *res, 
		    const igraph_vs_t *vids, igraph_bool_t directed, igraph_integer_t niter, 
		    igraph_real_t eps, igraph_real_t damping);
int igraph_rewire(igraph_t *graph, igraph_integer_t n, igraph_rewiring_t mode);
int igraph_subgraph(const igraph_t *graph, igraph_t *res, 
		    const igraph_vs_t *vids);
int igraph_average_path_length(const igraph_t *graph, igraph_real_t *res,
			       igraph_bool_t directed, igraph_bool_t unconn);
int igraph_simplify(igraph_t *graph, igraph_bool_t multiple, igraph_bool_t loops);
int igraph_transitivity(const igraph_t *graph, igraph_vector_t *res, 
			igraph_transitivity_type_t type); /* vvv + add */

/* TODO: degree.distribution (?) */

/* -------------------------------------------------- */
/* Components                                         */
/* -------------------------------------------------- */

int igraph_clusters(const igraph_t *graph, igraph_vector_t *membership, 
		    igraph_vector_t *csize, igraph_connectedness_t mode);
int igraph_is_connected(const igraph_t *graph, igraph_bool_t *res, 
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
				       igraph_integer_t niter, igraph_real_t maxdelta,
				       igraph_real_t area, igraph_real_t coolexp, 
				       igraph_real_t repulserad, igraph_bool_t use_seed);
int igraph_layout_grid_fruchterman_reingold(const igraph_t *graph, 
					    igraph_matrix_t *res,
					    igraph_integer_t niter, igraph_real_t maxdelta, 
					    igraph_real_t area, igraph_real_t coolexp,
					    igraph_real_t repulserad, 
					    igraph_real_t cellsize, igraph_bool_t use_seed);
int igraph_layout_kamada_kawai(const igraph_t *graph, igraph_matrix_t *res,
			       igraph_integer_t niter, igraph_real_t sigma, 
			       igraph_real_t initemp, igraph_real_t coolexp,
			       igraph_real_t kkconst);
int igraph_layout_springs(const igraph_t *graph, igraph_matrix_t *res,
			  igraph_real_t mass, igraph_real_t equil, igraph_real_t k,
			  igraph_real_t repeqdis, igraph_real_t kfr, igraph_bool_t repulse);
int igraph_layout_lgl(const igraph_t *graph, igraph_matrix_t *res,
		      igraph_integer_t maxiter, igraph_real_t maxdelta, 
		      igraph_real_t area, igraph_real_t coolexp,
		      igraph_real_t repulserad, igraph_real_t cellsize, igraph_integer_t root);
int igraph_layout_reingold_tilford(const igraph_t *graph, igraph_matrix_t *res,
              long int root);

int igraph_layout_random_3d(const igraph_t *graph, igraph_matrix_t *res);
int igraph_layout_sphere(const igraph_t *graph, igraph_matrix_t *res);
int igraph_layout_fruchterman_reingold_3d(const igraph_t *graph, 
					  igraph_matrix_t *res,
					  igraph_integer_t niter, igraph_real_t maxdelta,
					  igraph_real_t volume, igraph_real_t coolexp,
					  igraph_real_t repulserad,
					  igraph_bool_t use_seed);
int igraph_layout_kamada_kawai_3d(const igraph_t *graph, igraph_matrix_t *res,
				  igraph_integer_t niter, igraph_real_t sigma, 
				  igraph_real_t initemp, igraph_real_t coolexp, 
				  igraph_real_t kkconst);

int igraph_layout_merge_dla(igraph_vector_ptr_t *graphs,
			    igraph_vector_ptr_t *coords, 
			    igraph_matrix_t *res);

/* -------------------------------------------------- */
/* Visitor-like functions                             */
/* -------------------------------------------------- */

int igraph_bfs(igraph_t *graph, igraph_integer_t vid, igraph_neimode_t mode,
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
int igraph_get_edgelist(const igraph_t *graph, igraph_vector_t *res, igraph_bool_t bycol);

/* -------------------------------------------------- */
/* Read and write foreign formats                     */
/* -------------------------------------------------- */

int igraph_read_graph_edgelist(igraph_t *graph, FILE *instream, 
			       igraph_integer_t n, igraph_bool_t directed);
int igraph_read_graph_ncol(igraph_t *graph, FILE *instream,
			   igraph_strvector_t *predefnames, igraph_bool_t names, 
			  igraph_bool_t weights, igraph_bool_t directed);
int igraph_read_graph_lgl(igraph_t *graph, FILE *instream,
			  igraph_bool_t names, igraph_bool_t weights);
int igraph_read_graph_pajek(igraph_t *graph, FILE *instream);
int igraph_read_graph_graphml(igraph_t *graph, FILE *instream,
			      igraph_bool_t directed, int index);

int igraph_write_graph_edgelist(const igraph_t *graph, FILE *outstream);
int igraph_write_graph_ncol(const igraph_t *graph, FILE *outstream,
			    const char *names, const char *weights);
int igraph_write_graph_lgl(const igraph_t *graph, FILE *outstream,
			   const char *names, const char *weights,
			   igraph_bool_t isolates);
int igraph_write_graph_graphml(const igraph_t *graph, FILE *outstream);

/* -------------------------------------------------- */
/* Graph isomorphisms                                 */
/* -------------------------------------------------- */

int igraph_isoclass(const igraph_t *graph, int *isoclass);
int igraph_isomorphic(const igraph_t *graph1, const igraph_t *graph2,
		      igraph_bool_t *iso);
int igraph_isoclass_subgraph(const igraph_t *graph, igraph_vector_t *vids,
			     int *isoclass);
int igraph_isoclass_create(igraph_t *graph, igraph_integer_t size,
			   igraph_integer_t number, igraph_bool_t directed);

/* -------------------------------------------------- */
/* Graph motifs                                       */
/* -------------------------------------------------- */

int igraph_motifs_randesu(const igraph_t *graph, igraph_vector_t *hist, 
			  int size, igraph_vector_t *cut_prob);

int igraph_motifs_randesu_estimate(const igraph_t *graph, igraph_integer_t *est,
				   int size, igraph_vector_t *cut_prob, 
				   igraph_integer_t sample_size, 
				   igraph_vector_t *sample);
int igraph_motifs_randesu_no(const igraph_t *graph, igraph_integer_t *no,
			     int size, igraph_vector_t *cut_prob);

/* -------------------------------------------------- */
/* Progress handlers                                  */
/* -------------------------------------------------- */

typedef int igraph_progress_handler_t(const char *message, igraph_real_t percent,
				      void *data);

extern igraph_progress_handler_t igraph_progress_handler_stderr;

igraph_progress_handler_t *
igraph_set_progress_handler(igraph_progress_handler_t new_handler);

/* -------------------------------------------------- */
/* Graph operators                                    */
/* -------------------------------------------------- */

int igraph_disjoint_union(igraph_t *res, igraph_t *left, igraph_t *right);
int igraph_disjoint_union_many(igraph_t *res, igraph_vector_ptr_t *graphs);
int igraph_union(igraph_t *res, igraph_t *left, igraph_t *right);
int igraph_union_many(igraph_t *res, igraph_vector_ptr_t *graphs);
int igraph_intersection(igraph_t *res, igraph_t *left, igraph_t *right);
int igraph_intersection_many(igraph_t *res, igraph_vector_ptr_t *graphs);
int igraph_difference(igraph_t *res, igraph_t *orig, igraph_t *sub);
int igraph_complementer(igraph_t *res, igraph_t *graph, igraph_bool_t loops);
int igraph_compose(igraph_t *res, igraph_t *g1, igraph_t *g2);

/* -------------------------------------------------- */
/* Dynamics measurement                               */
/* -------------------------------------------------- */

int igraph_measure_dynamics_idage(const igraph_t *graph, igraph_integer_t start_vertex,
				  igraph_matrix_t *akl, 
				  igraph_matrix_t *sd, igraph_matrix_t *confint, 
				  igraph_matrix_t *no,
				  const igraph_vector_t *st, igraph_integer_t agebins,
				  igraph_integer_t maxind, igraph_real_t significance,
				  igraph_bool_t lno);
int igraph_measure_dynamics_idage_st(const igraph_t *graph, igraph_vector_t *res,
				     const igraph_matrix_t *akl);
int igraph_measure_dynamics_idage_debug(const igraph_t *graph, igraph_matrix_t *akl,
					igraph_matrix_t *sd, igraph_matrix_t *confint, 
					igraph_matrix_t *no,
					const igraph_vector_t *st, igraph_integer_t pagebins,
					igraph_integer_t pmaxind, igraph_real_t significance,
					igraph_vector_t *estimates, 
					igraph_integer_t est_ind, igraph_integer_t est_age,
					igraph_bool_t lno);

int igraph_measure_dynamics_id(const igraph_t *graph, igraph_integer_t start_vertex,
			       igraph_matrix_t *ak, igraph_matrix_t *sd,
			       igraph_matrix_t *confint, igraph_matrix_t *no,
			       const igraph_vector_t *st, igraph_integer_t pmaxind,
			       igraph_real_t significance, igraph_bool_t lno);
int igraph_measure_dynamics_id_st(const igraph_t *graph, 
				  igraph_vector_t *res, 
				  const igraph_matrix_t *ak);


/* -------------------------------------------------- */
/* Other, not graph related                           */
/* -------------------------------------------------- */

int igraph_running_mean(const igraph_vector_t *data, igraph_vector_t *res, 
			igraph_integer_t binwidth);
int igraph_random_sample(igraph_vector_t *res, igraph_integer_t l, igraph_integer_t h, 
			 igraph_integer_t length);
int igraph_convex_hull(const igraph_matrix_t *data, igraph_vector_t *resverts,
		       igraph_matrix_t *rescoords);

/* -------------------------------------------------- */
/* For internal use only, should move to other header */
/* -------------------------------------------------- */

typedef struct igraph_i_adjlist_t { 
  igraph_integer_t length;
  igraph_vector_t *adjs;
} igraph_i_adjlist_t;

int igraph_i_adjlist_init(const igraph_t *graph, igraph_i_adjlist_t *al, 
			  igraph_neimode_t mode);
void igraph_i_adjlist_destroy(igraph_i_adjlist_t *al);
/* igraph_vector_t *igraph_i_adjlist_get(const igraph_i_adjlist_t *al,  */
/* 			       igraph_integer_t no); */
#define igraph_i_adjlist_get(al, no) (&(al)->adjs[(long int)(no)])

extern unsigned int igraph_i_isoclass_3[];
extern unsigned int igraph_i_isoclass_4[];
extern unsigned int igraph_i_isoclass_3u[];
extern unsigned int igraph_i_isoclass_4u[];
extern unsigned int igraph_i_isoclass2_3[];
extern unsigned int igraph_i_isoclass2_4[];
extern unsigned int igraph_i_isoclass2_3u[];
extern unsigned int igraph_i_isoclass2_4u[];
extern unsigned int igraph_i_isoclass_3_idx[];
extern unsigned int igraph_i_isoclass_4_idx[];
extern unsigned int igraph_i_isoclass_3u_idx[];
extern unsigned int igraph_i_isoclass_4u_idx[];

/* -------------------------------------------------- */
/* Attributes, this is internal                       */
/* -------------------------------------------------- */

typedef struct igraph_attribute_table_t {
  int (*init)(igraph_t *graph);
  void (*destroy)(igraph_t *graph);
  int (*copy)(igraph_t *to, const igraph_t *from);
  int (*add_vertices)(igraph_t *graph, long int nv);
  void (*delete_vertices)(igraph_t *graph, const igraph_vector_t *eidx,
			  const igraph_vector_t *vidx);
  int (*add_edges)(igraph_t *graph, long int ne);
  void (*delete_edges)(igraph_t *graph, const igraph_vector_t *idx);
} igraph_attribute_table_t;

extern igraph_attribute_table_t *igraph_i_attribute_table;

igraph_attribute_table_t *
igraph_i_set_attribute_table(igraph_attribute_table_t * table);

#define IGRAPH_I_ATTRIBUTE_DESTROY(graph) \
        do {if ((graph)->attr) igraph_i_attribute_destroy(graph);} while(0)
#define IGRAPH_I_ATTRIBUTE_DELETE_VERTICES(graph, eidx, vidx) \
        do {if ((graph)->attr) igraph_i_attribute_delete_vertices((graph),(eidx),(vidx));} while(0)
#define IGRAPH_I_ATTRIBUTE_COPY(to,from) do { \
        int igraph_i_ret=0; \
        if (from->attr) { \
          IGRAPH_CHECK(igraph_i_ret=igraph_i_attribute_copy(to, from)); \
        } \
        if (igraph_i_ret != 0) { \
          IGRAPH_ERROR("", igraph_i_ret); \
        } \
   } while(0)        

int igraph_i_attribute_init(igraph_t *graph);
void igraph_i_attribute_destroy(igraph_t *graph);
int igraph_i_attribute_copy(igraph_t *to, const igraph_t *from);
int igraph_i_attribute_add_vertices(igraph_t *graph, long int nv);
void igraph_i_attribute_delete_vertices(igraph_t *graph, 
					const igraph_vector_t *eidx,
					const igraph_vector_t *vidx);
int igraph_i_attribute_add_edges(igraph_t *graph, long int ne);
void igraph_i_attribute_delete_edges(igraph_t *graph, 
				     const igraph_vector_t *idx);

__END_DECLS
  
#endif
