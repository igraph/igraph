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

typedef enum { IGRAPH_FILEFORMAT_EDGELIST=0,
	       IGRAPH_FILEFORMAT_NCOL,
	       IGRAPH_FILEFORMAT_PAJEK,
               IGRAPH_FILEFORMAT_LGL,
               IGRAPH_FILEFORMAT_GRAPHML } igraph_fileformat_type_t;

typedef enum { IGRAPH_REWIRING_SIMPLE=0 } igraph_rewiring_t;

typedef enum { IGRAPH_EDGEORDER_ID=0,
	       IGRAPH_EDGEORDER_FROM,
	       IGRAPH_EDGEORDER_TO } igraph_edgeorder_type_t;

typedef enum { IGRAPH_TO_DIRECTED_ARBITRARY=0,
	       IGRAPH_TO_DIRECTED_MUTUAL } igraph_to_directed_t;

typedef enum { IGRAPH_TO_UNDIRECTED_EACH=0,
	       IGRAPH_TO_UNDIRECTED_COLLAPSE } igraph_to_undirected_t;

typedef enum { IGRAPH_VCONN_NEI_ERROR=0,
	       IGRAPH_VCONN_NEI_INFINITY,
	       IGRAPH_VCONN_NEI_IGNORE }igraph_vconn_nei_t;

typedef enum { IGRAPH_SPINCOMM_UPDATE_SIMPLE=0,
	       IGRAPH_SPINCOMM_UPDATE_CONFIG } igraph_spincomm_update_t; 

typedef enum { IGRAPH_I_DONT_SIMPLIFY=0,
	       IGRAPH_I_SIMPLIFY,
	       IGRAPH_I_SORT_SIMPLIFY,
	       IGRAPH_I_SORT } igraph_i_lazy_adlist_simplify_t;
	       

/* -------------------------------------------------- */
/* Vertex selectors                                   */
/* -------------------------------------------------- */

#define IGRAPH_VS_ALL       0
#define IGRAPH_VS_ADJ       1
#define IGRAPH_VS_NONE      2
#define IGRAPH_VS_1         3
#define IGRAPH_VS_VECTORPTR 4
#define IGRAPH_VS_VECTOR    5
#define IGRAPH_VS_SEQ       6
#define IGRAPH_VS_NONADJ    7

typedef struct igraph_vs_t {
  int type;
  union {
    igraph_integer_t vid;    	        /* single vertex  */
    const igraph_vector_t *vecptr;      /* vector of vertices  */
    struct {
      igraph_integer_t vid;
      igraph_neimode_t mode;
    } adj;			        /* adjacent vertices  */
    struct {                           
      igraph_integer_t from;
      igraph_integer_t to;
    } seq;                              /* sequence of vertices from:to */
  } data;
} igraph_vs_t;
    
int igraph_vs_all(igraph_vs_t *vs);
igraph_vs_t igraph_vss_all(void);

int igraph_vs_adj(igraph_vs_t *vs, 
		  igraph_integer_t vid, igraph_neimode_t mode);
igraph_vs_t igraph_vss_adj(igraph_integer_t vid, igraph_neimode_t mode);

int igraph_vs_nonadj(igraph_vs_t *vs, igraph_integer_t vid, 
		     igraph_neimode_t mode);

int igraph_vs_none(igraph_vs_t *vs);
igraph_vs_t igraph_vss_none(void);

int igraph_vs_1(igraph_vs_t *vs, igraph_integer_t vid);
igraph_vs_t igraph_vss_1(igraph_integer_t vid);

int igraph_vs_vector(igraph_vs_t *vs,
		     const igraph_vector_t *v);
igraph_vs_t igraph_vss_vector(const igraph_vector_t *v);

int igraph_vs_vector_small(igraph_vs_t *vs, ...);			   

int igraph_vs_vector_copy(igraph_vs_t *vs,
			  const igraph_vector_t *v);

int igraph_vs_seq(igraph_vs_t *vs, igraph_integer_t from, igraph_integer_t to);
igraph_vs_t igraph_vss_seq(igraph_integer_t from, igraph_integer_t to);

void igraph_vs_destroy(igraph_vs_t *vs);

igraph_bool_t igraph_vs_is_all(igraph_vs_t *vs);

int igraph_vs_as_vector(const igraph_t *graph, igraph_vs_t vs, 
			igraph_vector_t *v);

/* -------------------------------------------------- */
/* Vertex iterators                                   */
/* -------------------------------------------------- */

#define IGRAPH_VIT_SEQ       0
#define IGRAPH_VIT_VECTOR    1
#define IGRAPH_VIT_VECTORPTR 2

typedef struct igraph_vit_t {
  int type;
  long int pos;
  long int start;
  long int end;
  const igraph_vector_t *vec;
} igraph_vit_t;

/**
 * \section IGRAPH_VIT Stepping over the vertices
 *
 * <para>After creating an iterator with \ref igraph_vit_create(), it
 * points to the first vertex in the vertex determined by the vertex
 * selector (if there is any). The \ref IGRAPH_VIT_NEXT() macro steps
 * to the next vertex, \ref IGRAPH_VIT_END() checks whether there are
 * more vertices to visit, \ref IGRAPH_VIT_SIZE() gives the total size
 * of the vertices visited so far and to be visited. \ref
 * IGRAPH_VIT_RESET() resets the iterator, it will point to the first
 * vertex again. Finally \ref IGRAPH_VIT_GET() gives the current vertex
 * pointed by the iterator (call this only if \ref IGRAPH_VIT_END()
 * is false).
 * </para>
 * <para>
 * Here is an example on how to step over the neighbors of vertex 0:
 * <informalexample><programlisting>
 * igraph_vs_t vs;
 * igraph_vit_t vit;
 * ...
 * igraph_vs_adj(&amp;vs, 0, IGRAPH_ALL);
 * igraph_vit_create(&amp;graph, vs, &amp;vit);
 * while (!IGRAPH_VIT_END(VIT)) {
 *   printf(" %li", (long int) IGRAPH_VIT_GET(vit));
 *   IGRAPH_VIT_NEXT(vit);
 * }
 * printf("\n");
 * ...
 * igraph_vit_destroy(&amp;vit);
 * igraph_vs_destroy(&amp;vs);
 * </programlisting></informalexample>
 * </para>
 */

/**
 * \define IGRAPH_VIT_NEXT
 * 
 * Steps the iterator to the next vertex. Only call this function if
 * \ref IGRAPH_VIT_END() returns false.
 * \param vit The vertex iterator to step.
 * 
 * Time complexity: O(1).
 */
#define IGRAPH_VIT_NEXT(vit)  (++((vit).pos))
/**
 * \define IGRAPH_VIT_END
 * 
 * Checks whether there are more vertices to step to.
 * \param vit The vertex iterator to check.
 * \return Logical value, if true there are no more vertices to step
 * to.
 * 
 * Time complexity: O(1).
 */
#define IGRAPH_VIT_END(vit)   ((vit).pos >= (vit).end)
/**
 * \define IGRAPH_VIT_SIZE
 * 
 * Gives the number of vertices in a vertex iterator.
 * \param vit The vertex iterator.
 * \return The number of vertices.
 * 
 * Time complexity: O(1).
 */
#define IGRAPH_VIT_SIZE(vit)  ((vit).end - (vit).start)
/**
 * \define IGRAPH_VIT_RESET
 * 
 * Resets a vertex iterator. After calling this macro the iterator
 * will point to the first vertex.
 * \param vit The vertex iterator.
 * 
 * Time complexity: O(1).
 */
#define IGRAPH_VIT_RESET(vit) ((vit).pos = (vit).start)
/**
 * \define IGRAPH_VIT_GET
 * 
 * Gives the vertex id of the current vertex poited to by the
 * iterator. 
 * \param vit The vertex iterator.
 * \return The vertex id of the current vertex.
 * 
 * Time complexity: O(1).
 */
#define IGRAPH_VIT_GET(vit)  \
  (((vit).type == IGRAPH_VIT_SEQ) ? (vit).pos : \
  VECTOR(*(vit).vec)[(vit).pos])

int igraph_vit_create(const igraph_t *graph, 
		      igraph_vs_t vs, igraph_vit_t *vit);
void igraph_vit_destroy(const igraph_vit_t *vit); 

int igraph_vit_as_vector(const igraph_vit_t *vit, igraph_vector_t *v);

/* -------------------------------------------------- */
/* Edge Selectors                                     */
/* -------------------------------------------------- */

#define IGRAPH_ES_ALL       0
#define IGRAPH_ES_ALLFROM   1
#define IGRAPH_ES_ALLTO     2
#define IGRAPH_ES_ADJ       3
#define IGRAPH_ES_NONE      4
#define IGRAPH_ES_1         5
#define IGRAPH_ES_VECTORPTR 6
#define IGRAPH_ES_VECTOR    7
#define IGRAPH_ES_SEQ       8
#define IGRAPH_ES_PAIRS     9
#define IGRAPH_ES_PATH      10
#define IGRAPH_ES_MULTIPAIRS 11

typedef struct igraph_es_t {
  int type;
  union {
    igraph_integer_t vid;
    igraph_integer_t eid;
    const igraph_vector_t *vecptr;
    struct {
      igraph_integer_t vid;
      igraph_neimode_t mode;
    } adj;
    struct {
      igraph_integer_t from;
      igraph_integer_t to;
    } seq;
    struct {
      const igraph_vector_t *ptr;
      igraph_bool_t mode;
    } path;
  } data;
} igraph_es_t;

int igraph_es_all(igraph_es_t *es, 
		  igraph_edgeorder_type_t order);
igraph_es_t igraph_ess_all(igraph_edgeorder_type_t order);

int igraph_es_adj(igraph_es_t *es, 
		  igraph_integer_t vid, igraph_neimode_t mode);

int igraph_es_none(igraph_es_t *es);
igraph_es_t igraph_ess_none(void);

int igraph_es_1(igraph_es_t *es, igraph_integer_t eid);
igraph_es_t igraph_ess_1(igraph_integer_t eid);

int igraph_es_vector(igraph_es_t *es,
		     const igraph_vector_t *v);
igraph_es_t igraph_ess_vector(const igraph_vector_t *v);

int igraph_es_fromto(igraph_es_t *es,
		     igraph_vs_t from, igraph_vs_t to);

int igraph_es_seq(igraph_es_t *es, igraph_integer_t from, igraph_integer_t to);
igraph_es_t igraph_ess_seq(igraph_integer_t from, igraph_integer_t to);

int igraph_es_pairs(igraph_es_t *es, const igraph_vector_t *v, 
		    igraph_bool_t directed);
int igraph_es_pairs_small(igraph_es_t *es, igraph_bool_t directed, ...);

int igraph_es_multipairs(igraph_es_t *es, const igraph_vector_t *v,
			 igraph_bool_t directed);

int igraph_es_path(igraph_es_t *es, const igraph_vector_t *v, 
		   igraph_bool_t directed);
int igraph_es_path_small(igraph_es_t *es, igraph_bool_t directed, ...);

void igraph_es_destroy(igraph_es_t *es);

igraph_bool_t igraph_es_is_all(igraph_es_t *es);

int igraph_es_as_vector(const igraph_t *graph, igraph_es_t es, 
			igraph_vector_t *v);

/* -------------------------------------------------- */
/* Edge Iterators                                     */
/* -------------------------------------------------- */

#define IGRAPH_EIT_SEQ       0
#define IGRAPH_EIT_VECTOR    1
#define IGRAPH_EIT_VECTORPTR 2

typedef struct igraph_eit_t {
  int type;
  long int pos;
  long int start;
  long int end;
  const igraph_vector_t *vec;
} igraph_eit_t;

/**
 * \section IGRAPH_EIT Stepping over the edges
 * 
 * <para>Just like for vertex iterators, macros are provided for
 * stepping over a sequence of edges: \ref IGRAPH_EIT_NEXT() goes to
 * the next edge, \ref IGRAPH_EIT_END() checks whether there are more
 * edges to visit, \ref IGRAPH_EIT_SIZE() gives the number of edges in
 * the edge sequence, \ref IGRAPH_EIT_RESET() resets the iterator to
 * the first edge and \ref IGRAPH_EIT_GET() returns the id of the
 * current edge.</para>
 */

/**
 * \define IGRAPH_EIT_NEXT
 * 
 * Steps the iterator to the next edge. Call this function only if
 * \ref IGRAPH_EIT_END() returns false.
 * \param eit The edge iterator to step.
 * 
 * Time complecity: O(1).
 */
#define IGRAPH_EIT_NEXT(eit) (++((eit).pos))
/**
 * \define IGRAPH_EIT_END
 * 
 * Checks whether there are more edges to step to.
 * \param wit The edge iterator to check.
 * \return Logical value, if true there are no more edges
 * to step to.
 *
 * Time complexity: O(1).
 */
#define IGRAPH_EIT_END(eit)   ((eit).pos >= (eit).end)
/**
 * \define IGRAPH_EIT_SIZE
 * 
 * Gives the number of edges in an edge iterator.
 * \param eit The edge iterator.
 * \return The number of edges.
 * 
 * Time complexity: O(1).
 */
#define IGRAPH_EIT_SIZE(eit)  ((eit).end - (eit).start)
/**
 * \define IGRAPH_EIT_RESET
 * 
 * Resets an ege iterator. After calling this macro the iterator will
 * point to the first edge.
 * \param eit The edge iterator.
 * 
 * Time complexity: O(1).
 */
#define IGRAPH_EIT_RESET(eit) ((eit).pos = (eit).start)
/**
 * \define IGRAPH_EIT_GET
 * 
 * Gives the edge id of the current edge pointed to by an iterator.
 * \param eit The edge iterator.
 * \return The id of the current edge.
 * 
 * Time complexity: O(1).
 */
#define IGRAPH_EIT_GET(eit)  \
  (((eit).type == IGRAPH_EIT_SEQ) ? (eit).pos : \
  VECTOR(*(eit).vec)[(eit).pos])

int igraph_eit_create(const igraph_t *graph, 
		      igraph_es_t es, igraph_eit_t *eit);
void igraph_eit_destroy(const igraph_eit_t *eit); 

int igraph_eit_as_vector(const igraph_eit_t *eit, igraph_vector_t *v);

/* -------------------------------------------------- */
/* Interface                                          */
/* -------------------------------------------------- */

int igraph_empty(igraph_t *graph, igraph_integer_t n, igraph_bool_t directed);
int igraph_empty_attrs(igraph_t *graph, igraph_integer_t n, igraph_bool_t directed, void *attr);
int igraph_destroy(igraph_t *graph);
int igraph_copy(igraph_t *to, const igraph_t *from);
int igraph_add_edges(igraph_t *graph, const igraph_vector_t *edges, 
		     void *attr);
int igraph_add_vertices(igraph_t *graph, igraph_integer_t nv, 
			void *attr);
int igraph_delete_edges(igraph_t *graph, igraph_es_t edges);
int igraph_delete_vertices(igraph_t *graph, const igraph_vs_t vertices);
igraph_integer_t igraph_vcount(const igraph_t *graph);
igraph_integer_t igraph_ecount(const igraph_t *graph);
int igraph_neighbors(const igraph_t *graph, igraph_vector_t *neis, igraph_integer_t vid, 
		     igraph_neimode_t mode); 
igraph_bool_t igraph_is_directed(const igraph_t *graph);
int igraph_degree(const igraph_t *graph, igraph_vector_t *res, 
		  const igraph_vs_t vids, igraph_neimode_t mode, 
		  igraph_bool_t loops);
int igraph_edge(const igraph_t *graph, igraph_integer_t eid, 
		igraph_integer_t *from, igraph_integer_t *to);		
int igraph_get_eid(const igraph_t *graph, igraph_integer_t *eid,
		   igraph_integer_t from, igraph_integer_t to,
		   igraph_bool_t directed);
int igraph_get_eids(const igraph_t *graph, igraph_vector_t *eids,
		    const igraph_vector_t *pairs, igraph_bool_t directed);
int igraph_adjacent(const igraph_t *graph, igraph_vector_t *eids, igraph_integer_t vid,
		    igraph_neimode_t mode);

/* -------------------------------------------------- */
/* Constructors, deterministic                        */
/* -------------------------------------------------- */

int igraph_create(igraph_t *graph, const igraph_vector_t *edges, igraph_integer_t n, 
		  igraph_bool_t directed);
int igraph_small(igraph_t *graph, igraph_integer_t n, igraph_bool_t directed, 
		 ...);
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
int igraph_extended_chordal_ring(igraph_t *graph, igraph_integer_t nodes, 
				 const igraph_matrix_t *W);
int igraph_connect_neighborhood(igraph_t *graph, igraph_integer_t order,
				igraph_neimode_t mode);

/* -------------------------------------------------- */
/* Constructors, games (=stochastic)                  */
/* -------------------------------------------------- */

int igraph_barabasi_game(igraph_t *graph, igraph_integer_t n, igraph_integer_t m, 
			 const igraph_vector_t *outseq, igraph_bool_t outpref, 
			 igraph_bool_t directed);
int igraph_nonlinear_barabasi_game(igraph_t *graph, igraph_integer_t n,
				   igraph_real_t power,
				   igraph_integer_t m,  
				   const igraph_vector_t *outseq,
				   igraph_bool_t outpref,
				   igraph_real_t zeroappeal,
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
int igraph_barabasi_aging_game(igraph_t *graph, 
			       igraph_integer_t nodes,
			       igraph_integer_t m,
			       const igraph_vector_t *outseq,
			       igraph_bool_t outpref,
			       igraph_real_t pa_exp,
			       igraph_real_t aging_exp,
			       igraph_integer_t aging_bin,
			       igraph_real_t zero_deg_appeal,
			       igraph_real_t zero_age_appeal,
			       igraph_real_t deg_coef,
			       igraph_real_t age_coef,
			       igraph_bool_t directed);
int igraph_recent_degree_game(igraph_t *graph, igraph_integer_t n,
			      igraph_real_t power,
			      igraph_integer_t window,
			      igraph_integer_t m,  
			      const igraph_vector_t *outseq,
			      igraph_bool_t outpref,
			      igraph_real_t zero_appeal,
			      igraph_bool_t directed);
int igraph_recent_degree_aging_game(igraph_t *graph,
				    igraph_integer_t nodes,
				    igraph_integer_t m, 
				    const igraph_vector_t *outseq,
				    igraph_bool_t outpref,
				    igraph_real_t pa_exp,
				    igraph_real_t aging_exp,
				    igraph_integer_t aging_bin,
				    igraph_integer_t window,
				    igraph_real_t zero_appeal,
				    igraph_bool_t directed);
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
int igraph_grg_game(igraph_t *graph, igraph_integer_t nodes,
		    igraph_real_t radius, igraph_bool_t torus);
int igraph_preference_game(igraph_t *graph, igraph_integer_t nodes,
			   igraph_integer_t types, igraph_vector_t *type_dist,
			   igraph_matrix_t *pref_matrix,
			   igraph_vector_t *node_type_vec,
			   igraph_bool_t directed, igraph_bool_t loops);
int igraph_asymmetric_preference_game(igraph_t *graph, igraph_integer_t nodes,
				      igraph_integer_t types,
				      igraph_matrix_t *type_dist_matrix,
				      igraph_matrix_t *pref_matrix,
				      igraph_vector_t *node_type_in_vec,
				      igraph_vector_t *node_type_out_vec,
				      igraph_bool_t loops);

int igraph_rewire_edges(igraph_t *graph, igraph_real_t prob);
int igraph_watts_strogatz_game(igraph_t *graph, igraph_integer_t dim,
			       igraph_integer_t size, igraph_integer_t nei,
			       igraph_real_t p);

int igraph_lastcit_game(igraph_t *graph, 
			igraph_integer_t nodes, igraph_integer_t edges_per_node, 
			igraph_integer_t agebins,
			const igraph_vector_t *preference, igraph_bool_t directed);

int igraph_cited_type_game(igraph_t *graph, igraph_integer_t nodes,
			   const igraph_vector_t *types,
			   const igraph_vector_t *pref,
			   igraph_integer_t edges_per_step,
			   igraph_bool_t directed);

int igraph_citing_cited_type_game(igraph_t *graph, igraph_integer_t nodes,
				  const igraph_vector_t *types,
				  const igraph_matrix_t *pref,
				  igraph_integer_t edges_per_step,
				  igraph_bool_t directed);

/* -------------------------------------------------- */
/* Basic query functions                              */
/* -------------------------------------------------- */

int igraph_are_connected(const igraph_t *graph, igraph_integer_t v1, igraph_integer_t v2, igraph_bool_t *res);

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
		     const igraph_vs_t vids, igraph_neimode_t mode);
int igraph_shortest_paths(const igraph_t *graph, igraph_matrix_t *res, 
			  const igraph_vs_t from, igraph_neimode_t mode);
int igraph_get_shortest_paths(const igraph_t *graph, igraph_vector_ptr_t *res,
			      igraph_integer_t from, const igraph_vs_t to, 
			      igraph_neimode_t mode);
int igraph_get_all_shortest_paths(const igraph_t *graph,
				  igraph_vector_ptr_t *res, 
				  igraph_vector_t *nrgeo,
				  igraph_integer_t from, const igraph_vs_t to, 
				  igraph_neimode_t mode);
int igraph_subcomponent(const igraph_t *graph, igraph_vector_t *res, igraph_real_t vid, 
			igraph_neimode_t mode);	
int igraph_betweenness (const igraph_t *graph, igraph_vector_t *res, 
			const igraph_vs_t vids, igraph_bool_t directed);
int igraph_edge_betweenness (const igraph_t *graph, igraph_vector_t *result,
			     igraph_bool_t directed); /* eee + add */
int igraph_pagerank(const igraph_t *graph, igraph_vector_t *res, 
		    const igraph_vs_t vids, igraph_bool_t directed, igraph_integer_t niter, 
		    igraph_real_t eps, igraph_real_t damping);
int igraph_rewire(igraph_t *graph, igraph_integer_t n, igraph_rewiring_t mode);
int igraph_subgraph(const igraph_t *graph, igraph_t *res, 
		    const igraph_vs_t vids);
int igraph_average_path_length(const igraph_t *graph, igraph_real_t *res,
			       igraph_bool_t directed, igraph_bool_t unconn);
int igraph_simplify(igraph_t *graph, igraph_bool_t multiple, igraph_bool_t loops);
int igraph_transitivity_undirected(const igraph_t *graph, 
				   igraph_real_t *res);
int igraph_transitivity_local_undirected(const igraph_t *graph, 
					 igraph_vector_t *res,
					 const igraph_vs_t vids);
int igraph_transitivity_local_undirected1(const igraph_t *graph, 
					  igraph_vector_t *res,
					  const igraph_vs_t vids);
int igraph_transitivity_local_undirected2(const igraph_t *graph, 
					  igraph_vector_t *res,
					  const igraph_vs_t vids);
int igraph_transitivity_local_undirected4(const igraph_t *graph, 
					  igraph_vector_t *res,
					  const igraph_vs_t vids);
int igraph_transitivity_avglocal_undirected(const igraph_t *graph,
					    igraph_real_t *res);
int igraph_reciprocity(const igraph_t *graph, igraph_real_t *res,
		       igraph_bool_t ignore_loops);

int igraph_constraint(const igraph_t *graph, igraph_vector_t *res,
		      igraph_vs_t vids, const igraph_vector_t *weights);
int igraph_maxdegree(const igraph_t *graph, igraph_integer_t *res,
		     igraph_vs_t vids, igraph_neimode_t mode, 
		     igraph_bool_t loops);
int igraph_density(const igraph_t *graph, igraph_real_t *res, 
		   igraph_bool_t loops);

int igraph_neighborhood_size(const igraph_t *graph, igraph_vector_t *res,
			     igraph_vs_t vids, igraph_integer_t order, 
			     igraph_neimode_t mode);
int igraph_neighborhood(const igraph_t *graph, igraph_vector_ptr_t *res,
			igraph_vs_t vids, igraph_integer_t order,
			igraph_neimode_t mode);
int igraph_neighborhood_graphs(const igraph_t *graph, igraph_vector_ptr_t *res,
			       igraph_vs_t vids, igraph_integer_t order,
			       igraph_neimode_t mode);
int igraph_topological_sorting(const igraph_t *graph, igraph_vector_t *res,
			       igraph_neimode_t mode);

/* TODO: degree.distribution (?) */

/* -------------------------------------------------- */
/* Spectral Properties                                */
/* -------------------------------------------------- */

int igraph_laplacian(const igraph_t *graph, igraph_matrix_t *res,
		     igraph_bool_t normalized);

/* -------------------------------------------------- */
/* Components                                         */
/* -------------------------------------------------- */

int igraph_clusters(const igraph_t *graph, igraph_vector_t *membership, 
		    igraph_vector_t *csize, igraph_integer_t *no,
		    igraph_connectedness_t mode);
int igraph_is_connected(const igraph_t *graph, igraph_bool_t *res, 
			igraph_connectedness_t mode);
int igraph_decompose(const igraph_t *graph, igraph_vector_ptr_t *components, 
		     igraph_connectedness_t mode, 
		     long int maxcompno, long int minelements);

/* TODO: cluster.distribution (?) */

/* -------------------------------------------------- */
/* Cliques, maximal independent vertex sets           */
/* -------------------------------------------------- */

int igraph_cliques(const igraph_t *graph, igraph_vector_ptr_t *res,
                   igraph_integer_t min_size, igraph_integer_t max_size);
int igraph_largest_cliques(const igraph_t *graph, 
			   igraph_vector_ptr_t *cliques);
int igraph_maximal_cliques(const igraph_t *graph,
			   igraph_vector_ptr_t *res);
int igraph_clique_number(const igraph_t *graph, igraph_integer_t *no);
int igraph_independent_vertex_sets(const igraph_t *graph,
				   igraph_vector_ptr_t *res,
				   igraph_integer_t min_size,
				   igraph_integer_t max_size);
int igraph_largest_independent_vertex_sets(const igraph_t *graph,
					   igraph_vector_ptr_t *res);
int igraph_maximal_independent_vertex_sets(const igraph_t *graph,
					   igraph_vector_ptr_t *res);
int igraph_independence_number(const igraph_t *graph, igraph_integer_t *no);

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
		      const igraph_vs_t vids);
int igraph_bibcoupling(const igraph_t *graph, igraph_matrix_t *res, 
		       const igraph_vs_t vids);

/* -------------------------------------------------- */
/* Community Structure                                */
/* -------------------------------------------------- */

/* TODO: eb.community */
/* TODO: cut.community */
/* TODO: edge.type.matrix */
/* TODO: modularity */
/* TODO:  */

int igraph_spinglass_community(const igraph_t *graph,
			       const igraph_vector_t *weights,
			       igraph_real_t *modularity,
			       igraph_real_t *temperature,
			       igraph_vector_t *membership, 
			       igraph_vector_t *csize, 
			       igraph_integer_t spins,
			       igraph_bool_t parupdate,
			       igraph_real_t starttemp,
			       igraph_real_t stoptemp,
			       igraph_real_t coolfact,
			       igraph_spincomm_update_t update_rule,
			       igraph_real_t gamma);

int igraph_spinglass_my_community(const igraph_t *graph,
				  const igraph_vector_t *weights,
				  igraph_integer_t vertex,
				  igraph_vector_t *community,
				  igraph_real_t *cohesion,
				  igraph_real_t *adhesion,
				  igraph_integer_t *inner_links,
				  igraph_integer_t *outer_links,
				  igraph_integer_t spins,
				  igraph_spincomm_update_t update_rule,
				  igraph_real_t gamma);

/* -------------------------------------------------- */
/* Conversion                                         */
/* -------------------------------------------------- */

int igraph_get_adjacency(const igraph_t *graph, igraph_matrix_t *res,
			 igraph_get_adjacency_t type);
int igraph_get_edgelist(const igraph_t *graph, igraph_vector_t *res, igraph_bool_t bycol);

int igraph_to_directed(igraph_t *graph, 
		       igraph_to_directed_t flags);
int igraph_to_undirected(igraph_t *graph,
			 igraph_to_undirected_t flags);

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
			      int index);
int igraph_read_graph_dimacs(igraph_t *graph, FILE *instream,
			     igraph_integer_t *source, 
			     igraph_integer_t *target, 
			     igraph_vector_t *capacity, 
			     igraph_bool_t directed);
int igraph_read_graph_graphdb(igraph_t *graph, FILE *instream, 
			      igraph_bool_t directed);

int igraph_write_graph_edgelist(const igraph_t *graph, FILE *outstream);
int igraph_write_graph_ncol(const igraph_t *graph, FILE *outstream,
			    const char *names, const char *weights);
int igraph_write_graph_lgl(const igraph_t *graph, FILE *outstream,
			   const char *names, const char *weights,
			   igraph_bool_t isolates);
int igraph_write_graph_graphml(const igraph_t *graph, FILE *outstream);
int igraph_write_graph_pajek(const igraph_t *graph, FILE *outstream);
int igraph_write_graph_dimacs(const igraph_t *graph, FILE *outstream,
			      long int source, long int target,
			      const igraph_vector_t *capacity);

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
int igraph_isomorphic_vf2(const igraph_t *graph1, const igraph_t *graph2, 
			  igraph_bool_t *iso);

/* -------------------------------------------------- */
/* Graph motifs                                       */
/* -------------------------------------------------- */

int igraph_motifs_randesu(const igraph_t *graph, igraph_vector_t *hist, 
			  int size, const igraph_vector_t *cut_prob);

int igraph_motifs_randesu_estimate(const igraph_t *graph, igraph_integer_t *est,
				   int size, const igraph_vector_t *cut_prob, 
				   igraph_integer_t sample_size, 
				   const igraph_vector_t *sample);
int igraph_motifs_randesu_no(const igraph_t *graph, igraph_integer_t *no,
			     int size, const igraph_vector_t *cut_prob);

/* -------------------------------------------------- */
/* Progress handlers                                  */
/* -------------------------------------------------- */

typedef int igraph_progress_handler_t(const char *message, igraph_real_t percent,
				      void *data);

extern igraph_progress_handler_t igraph_progress_handler_stderr;

igraph_progress_handler_t *
igraph_set_progress_handler(igraph_progress_handler_t new_handler);

int igraph_progress(const char *message, igraph_real_t percent, void *data);

/* -------------------------------------------------- */
/* Graph operators                                    */
/* -------------------------------------------------- */

int igraph_disjoint_union(igraph_t *res, 
			  const igraph_t *left, const igraph_t *right);
int igraph_disjoint_union_many(igraph_t *res, 
			       const igraph_vector_ptr_t *graphs);
int igraph_union(igraph_t *res, const igraph_t *left, const igraph_t *right);
int igraph_union_many(igraph_t *res, const igraph_vector_ptr_t *graphs);
int igraph_intersection(igraph_t *res, 
			const igraph_t *left, const igraph_t *right);
int igraph_intersection_many(igraph_t *res, const igraph_vector_ptr_t *graphs);
int igraph_difference(igraph_t *res, 
		      const igraph_t *orig, const igraph_t *sub);
int igraph_complementer(igraph_t *res, const igraph_t *graph, 
			igraph_bool_t loops);
int igraph_compose(igraph_t *res, const igraph_t *g1, const igraph_t *g2);

/* -------------------------------------------------- */
/* MAximum flows, minimum cuts & such                 */
/* -------------------------------------------------- */

int igraph_maxflow_value(const igraph_t *graph, igraph_real_t *value,
			 igraph_integer_t source, igraph_integer_t target,
			 const igraph_vector_t *capacity);
int igraph_mincut_value(const igraph_t *graph, igraph_real_t *res, 
			const igraph_vector_t *capacity);
int igraph_st_mincut_value(const igraph_t *graph, igraph_real_t *res,
                           igraph_integer_t source, igraph_integer_t target,
			   const igraph_vector_t *capacity);

int igraph_st_vertex_connectivity(const igraph_t *graph, 
				  igraph_integer_t *res,
				  igraph_integer_t source,
				  igraph_integer_t target,
				  igraph_vconn_nei_t neighbors);
int igraph_vertex_connectivity(const igraph_t *graph, igraph_integer_t *res);
int igraph_st_edge_connectivity(const igraph_t *graph, igraph_integer_t *res,
				igraph_integer_t source, 
				igraph_integer_t target);
int igraph_edge_connectivity(const igraph_t *graph, igraph_integer_t *res);
int igraph_edge_disjoint_paths(const igraph_t *graph, igraph_integer_t *res,
			       igraph_integer_t source, 
			       igraph_integer_t target);
int igraph_vertex_disjoint_paths(const igraph_t *graph, igraph_integer_t *res,
				 igraph_integer_t source,
				 igraph_integer_t target);
int igraph_adhesion(const igraph_t *graph, igraph_integer_t *res);
int igraph_cohesion(const igraph_t *graph, igraph_integer_t *res);

/* -------------------------------------------------- */
/* K-Cores                                            */
/* -------------------------------------------------- */

int igraph_coreness(const igraph_t *graph, igraph_vector_t *cores,
		    igraph_neimode_t mode);

/* -------------------------------------------------- */
/* Dynamics measurement                               */
/* -------------------------------------------------- */

int igraph_measure_dynamics_idage(const igraph_t *graph,
				  igraph_matrix_t *akl, 
				  igraph_matrix_t *sd, 
				  igraph_matrix_t *no,
				  igraph_matrix_t *cites,
				  const igraph_vector_t *st, igraph_integer_t agebins,
				  igraph_integer_t maxind);
int igraph_measure_dynamics_idage_st(const igraph_t *graph, igraph_vector_t *res,
				     const igraph_matrix_t *akl);
int igraph_measure_dynamics_idage_expected(const igraph_t *graph,
					   igraph_matrix_t *res,
					   const igraph_matrix_t *akl,
					   const igraph_vector_t *st,
					   igraph_integer_t pmaxind);

int igraph_measure_dynamics_idwindowage(const igraph_t *graph, 
					igraph_matrix_t *akl, 
					igraph_matrix_t *sd, 
					const igraph_vector_t *st, 
					igraph_integer_t pagebins,
					igraph_integer_t pmaxind, 
					igraph_integer_t time_window);
int igraph_measure_dynamics_idwindowage_st(const igraph_t *graph, 
					   igraph_vector_t *res,
					   const igraph_matrix_t *akl,
					   igraph_integer_t time_window);

int igraph_measure_dynamics_citedcat_id_age(const igraph_t *graph,
					    igraph_array3_t *adkl,
					    igraph_array3_t *sd,
					    const igraph_vector_t *st,
					    const igraph_vector_t *cats,
					    igraph_integer_t pno_cats,
					    igraph_integer_t pagebins,
					    igraph_integer_t pmaxind);

int igraph_measure_dynamics_citedcat_id_age_st(const igraph_t *graph,
					       igraph_vector_t *res,
					       const igraph_array3_t *adkl,
					       const igraph_vector_t *cats, 
					       igraph_integer_t pno_cats);

int igraph_measure_dynamics_citingcat_id_age(const igraph_t *graph,
					     igraph_array3_t *adkl,
					     igraph_array3_t *sd,
					     const igraph_vector_t *st,
					     const igraph_vector_t *cats,
					     igraph_integer_t pno_cats,
					     igraph_integer_t pagebins,
					     igraph_integer_t pmaxind);
int igraph_measure_dynamics_citingcat_id_age_st(const igraph_t *graph,
						igraph_vector_t *res,
						const igraph_array3_t *adkl,
						const igraph_vector_t *cats,
						igraph_integer_t pno_cats);

int igraph_measure_dynamics_id(const igraph_t *graph,
			       igraph_matrix_t *ak, igraph_matrix_t *sd,
			       igraph_matrix_t *no, igraph_vector_t *cites,
			       igraph_vector_t *debug,
			       igraph_integer_t debugdeg,
			       const igraph_vector_t *st, igraph_integer_t pmaxind);
int igraph_measure_dynamics_id_st(const igraph_t *graph, 
				  igraph_vector_t *res, 
				  const igraph_matrix_t *ak);
int igraph_measure_dynamics_id_expected(const igraph_t *graph,
					igraph_vector_t *res,
					const igraph_vector_t *ak,
					const igraph_vector_t *st,
					igraph_integer_t pmaxind);
int igraph_measure_dynamics_id_expected2(const igraph_t *graph,
					 igraph_vector_t *res,
					 const igraph_vector_t *ak,
					 const igraph_vector_t *st,
					 igraph_integer_t pmaxind);

int igraph_measure_dynamics_d_d(const igraph_t *graph,
				const igraph_vector_t *ntime,
				const igraph_vector_t *etime,
				igraph_integer_t events,
				igraph_matrix_t *akk,
				igraph_matrix_t *sd,
				const igraph_vector_t *st,
				igraph_integer_t pmadeg);

int igraph_measure_dynamics_d_d_st(const igraph_t *graph,
				   const igraph_vector_t *ntime,
				   const igraph_vector_t *etime,
				   const igraph_matrix_t *akk,
				   igraph_integer_t events,
				   igraph_integer_t maxtotaldeg,
				   igraph_vector_t *st);

int igraph_measure_dynamics_idwindow(const igraph_t *graph, 
				     igraph_matrix_t *ak, 
				     igraph_matrix_t *sd,
				     const igraph_vector_t *st,
				     igraph_integer_t pmaxind,
				     igraph_integer_t time_window);

int igraph_measure_dynamics_idwindow_st(const igraph_t *graph,
					igraph_vector_t *res,
					const igraph_matrix_t *ak,
					igraph_integer_t time_window);

int igraph_measure_dynamics_lastcit(const igraph_t *graph, igraph_vector_t *al,
				    igraph_vector_t *sd,
				    igraph_vector_t *no,
				    const igraph_vector_t *st,
				    igraph_integer_t pagebins);
int igraph_measure_dynamics_lastcit_st(const igraph_t *graph, 
				       igraph_vector_t *res,
				       const igraph_vector_t *al);

int igraph_measure_dynamics_age(const igraph_t *graph, 
				igraph_vector_t *al,
				igraph_vector_t *sd,
				igraph_vector_t *no,
				const igraph_vector_t *st,
				igraph_integer_t pagebins);
int igraph_measure_dynamics_age_st(const igraph_t *graph, 
				   igraph_vector_t *res,
				   const igraph_vector_t *al);

int igraph_measure_dynamics_citedcat(const igraph_t *graph, 
				     const igraph_vector_t *cats,
				     igraph_integer_t pnocats,
				     igraph_vector_t *ak, 
				     igraph_vector_t  *sd,
				     igraph_vector_t *no,
				     const igraph_vector_t *st);
int igraph_measure_dynamics_citedcat_st(const igraph_t *graph,
					igraph_vector_t *res,
					const igraph_vector_t *ak,
					const igraph_vector_t *cats,
					igraph_integer_t pnocats);

int igraph_measure_dynamics_citingcat_citedcat(const igraph_t *graph,
					       igraph_matrix_t *agd,
					       igraph_matrix_t *sd,
					       igraph_matrix_t *no,
					       const igraph_vector_t *st,
					       const igraph_vector_t *cats,
					       igraph_integer_t pnocats);
int igraph_measure_dynamics_citingcat_citedcat_st(const igraph_t *graph,
						  igraph_vector_t *res,
						  const igraph_matrix_t *agd,
						  const igraph_vector_t *cats,
						  igraph_integer_t pnocats);

/* -------------------------------------------------- */
/* Network evolution measurement, new implementation  */
/* -------------------------------------------------- */

int igraph_evolver_d(const igraph_t *graph,
		     igraph_integer_t niter,
		     igraph_vector_t *kernel,		     
		     igraph_vector_t *sd,
		     igraph_vector_t *norm,
		     igraph_vector_t *cites,
		     igraph_vector_t *expected,
		     const igraph_vector_t *debug,
		     igraph_vector_ptr_t *debugres);
int igraph_evolver_mes_d(const igraph_t *graph,
			 igraph_vector_t *kernel,
			 igraph_vector_t *sd,
			 igraph_vector_t *norm,
			 igraph_vector_t *cites,
			 const igraph_vector_t *debug,
			 igraph_vector_ptr_t *debugres,
			 const igraph_vector_t *st,
			 igraph_integer_t pmaxind);
int igraph_evolver_st_d(const igraph_t *graph,
			igraph_vector_t *st,
			const igraph_vector_t *kernel);
int igraph_evolver_exp_d(const igraph_t *graphm,
			 igraph_vector_t *expected,
			 const igraph_vector_t *kernel,
			 const igraph_vector_t *st,
			 igraph_integer_t pmaxind);

int igraph_evolver_ad(const igraph_t *graph,
		      igraph_integer_t niter,
		      igraph_integer_t agebins,
		      igraph_matrix_t *kernel,
		      igraph_matrix_t *sd,
		      igraph_matrix_t *norm,
		      igraph_matrix_t *cites,
		      igraph_matrix_t *expected,
		      const igraph_matrix_t *debug,
		      igraph_vector_ptr_t *debugres);
int igraph_evolver_mes_ad(const igraph_t *graph,
			  igraph_matrix_t *kernel,
			  igraph_matrix_t *sd,
			  igraph_matrix_t *norm,
			  igraph_matrix_t *cites,
			  const igraph_matrix_t *debug,
			  igraph_vector_ptr_t *debugres,
			  const igraph_vector_t *st,
			  igraph_integer_t pmaxind,
			  igraph_integer_t agebins);
int igraph_evolver_ad2(const igraph_t *graph,
		      igraph_integer_t niter,
		      igraph_integer_t agebins,
		      igraph_matrix_t *kernel,
		      igraph_matrix_t *sd,
		      igraph_matrix_t *norm,
		      igraph_matrix_t *cites,
		      igraph_matrix_t *expected,
		      const igraph_matrix_t *debug,
		      igraph_vector_ptr_t *debugres);
int igraph_evolver_mes_ad2(const igraph_t *graph,
			  igraph_matrix_t *kernel,
			  igraph_matrix_t *sd,
			  igraph_matrix_t *norm,
			  igraph_matrix_t *cites,
			  const igraph_matrix_t *debug,
			  igraph_vector_ptr_t *debugres,
			  const igraph_vector_t *st,
			  igraph_integer_t pmaxind,
			  igraph_integer_t agebins);
int igraph_evolver_st_ad(const igraph_t *graph,
			 igraph_vector_t *st,
			 const igraph_matrix_t *kernel);
int igraph_evolver_exp_ad(const igraph_t *graph,
			  igraph_matrix_t *expected,
			  const igraph_matrix_t *kernel,
			  const igraph_vector_t *st,
			  igraph_integer_t pmaxind,
			  igraph_integer_t agebins);

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
int igraph_i_adjlist_init_complementer(const igraph_t *graph,
				       igraph_i_adjlist_t *al, 
				       igraph_neimode_t mode,
				       igraph_bool_t loops);
void igraph_i_adjlist_destroy(igraph_i_adjlist_t *al);
void igraph_i_adjlist_sort(igraph_i_adjlist_t *al);
int igraph_i_adjlist_simplify(igraph_i_adjlist_t *al);
/* igraph_vector_t *igraph_i_adjlist_get(const igraph_i_adjlist_t *al,  */
/* 			       igraph_integer_t no); */
#define igraph_i_adjlist_get(al, no) (&(al)->adjs[(long int)(no)])

typedef struct igraph_i_adjedgelist_t {
  igraph_integer_t length;
  igraph_vector_t *adjs;
} igraph_i_adjedgelist_t;

int igraph_i_adjedgelist_init(const igraph_t *graph, 
			      igraph_i_adjedgelist_t *eal, 
			      igraph_neimode_t mode);
void igraph_i_adjedgelist_destroy(igraph_i_adjedgelist_t *ael);
#define igraph_i_adjedgelist_get(ael, no) (&(ael)->adjs[(long int)(no)])

typedef struct igraph_i_lazy_adjlist_t {
  const igraph_t *graph;
  igraph_integer_t length;
  igraph_vector_t **adjs;
  igraph_neimode_t mode;
  igraph_i_lazy_adlist_simplify_t simplify;
  igraph_vector_t svect;
} igraph_i_lazy_adjlist_t;

int igraph_i_lazy_adjlist_init(const igraph_t *graph,
			       igraph_i_lazy_adjlist_t *al,
			       igraph_neimode_t mode,
			       igraph_i_lazy_adlist_simplify_t simplify);
void igraph_i_lazy_adjlist_destroy(igraph_i_lazy_adjlist_t *al);
/* igraph_vector_t *igraph_i_lazy_adjlist_get(igraph_i_lazy_adjlist_t *al, */
/* 					   igraph_integer_t no); */
#define igraph_i_lazy_adjlist_get(al, no) \
  ((al)->adjs[(long int)(no)] != 0 ? ((al)->adjs[(long int)(no)]) : \
   (igraph_i_lazy_adjlist_get_real(al, no)))
igraph_vector_t *igraph_i_lazy_adjlist_get_real(igraph_i_lazy_adjlist_t *al,
						igraph_integer_t no);

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

#include "attributes.h"

__END_DECLS
  
#endif
