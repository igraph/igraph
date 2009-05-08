/* -*- mode: C -*-  */
/* 
   IGraph library.
   Copyright (C) 2009  Gabor Csardi <csardi@rmki.kfki.hu>
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

#ifndef IGRAPH_ADJLIST_H
#define IGRAPH_ADJLIST_H

#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
# define __BEGIN_DECLS extern "C" {
# define __END_DECLS }
#else
# define __BEGIN_DECLS /* empty */
# define __END_DECLS /* empty */
#endif

#include "igraph_constants.h"
#include "igraph_types.h"
#include "igraph_datatype.h"

__BEGIN_DECLS

typedef struct igraph_adjlist_t { 
  igraph_integer_t length;
  igraph_vector_t *adjs;
} igraph_adjlist_t;

int igraph_adjlist_init(const igraph_t *graph, igraph_adjlist_t *al, 
			  igraph_neimode_t mode);
igraph_integer_t igraph_adjlist_size(const igraph_adjlist_t *al); 
int igraph_adjlist_init_complementer(const igraph_t *graph,
				     igraph_adjlist_t *al, 
				     igraph_neimode_t mode,
				     igraph_bool_t loops);
void igraph_adjlist_destroy(igraph_adjlist_t *al);
void igraph_adjlist_sort(igraph_adjlist_t *al);
int igraph_adjlist_simplify(igraph_adjlist_t *al);
int igraph_adjlist_remove_duplicate(const igraph_t *graph, 
				    igraph_adjlist_t *al);
int igraph_adjlist_print(const igraph_adjlist_t *al, FILE *outfile);
/* igraph_vector_t *igraph_adjlist_get(const igraph_adjlist_t *al,  */
/* 			       igraph_integer_t no); */
/**
 * \define igraph_adjlist_get
 * Query a vector in an adjlist
 * 
 * Returns a pointer to an <type>igraph_vector_t</type> object from an
 * adjacency list. The vector can be modified as desired. 
 * \param al The adjacency list object.
 * \param no The vertex of which the vertex of adjacent vertices are
 *   returned.
 * \return Pointer to the <type>igraph_vector_t</type> object.
 * 
 * Time complexity: O(1).
 */
#define igraph_adjlist_get(al,no) (&(al)->adjs[(long int)(no)])

int igraph_adjlist(igraph_t *graph, const igraph_adjlist_t *adjlist,
		   igraph_bool_t directed, igraph_bool_t duplicate);

typedef struct igraph_adjedgelist_t {
  igraph_integer_t length;
  igraph_vector_t *adjs;
} igraph_adjedgelist_t;

int igraph_adjedgelist_init(const igraph_t *graph, 
			    igraph_adjedgelist_t *eal, 
			    igraph_neimode_t mode);
void igraph_adjedgelist_destroy(igraph_adjedgelist_t *ael);
int igraph_adjedgelist_remove_duplicate(const igraph_t *graph,
					igraph_adjedgelist_t *al);
int igraph_adjedgelist_print(const igraph_adjedgelist_t *al, FILE *outfile);

/**
 * \define igraph_adjedgelist_get
 * Query a vector in an adjedgelist
 *
 * Returns a pointer to an <type>igraph_vector_t</type> object from an
 * adjacency list containing edge ids. The vector can be modified,
 * resized, etc. as desired. 
 * \param graph ael The edge adjacency list.
 * \param no The vertex of which the adjacent edges are returned.
 * \return Pointer to an <type>igraph_vector_t</type> object.
 * 
 * Time complexity: O(1).
 */
#define igraph_adjedgelist_get(ael,no) (&(ael)->adjs[(long int)(no)])

typedef struct igraph_lazy_adjlist_t {
  const igraph_t *graph;
  igraph_integer_t length;
  igraph_vector_t **adjs;
  igraph_neimode_t mode;
  igraph_lazy_adlist_simplify_t simplify;
} igraph_lazy_adjlist_t;

int igraph_lazy_adjlist_init(const igraph_t *graph,
			       igraph_lazy_adjlist_t *al,
			       igraph_neimode_t mode,
			       igraph_lazy_adlist_simplify_t simplify);
void igraph_lazy_adjlist_destroy(igraph_lazy_adjlist_t *al);
/* igraph_vector_t *igraph_lazy_adjlist_get(igraph_lazy_adjlist_t *al, */
/* 					   igraph_integer_t no); */
/**
 * \define igraph_lazy_adjlist_get
 * Query neighbor vertices
 * 
 * If the function is called for the first time for a vertex then the
 * result is stored in the adjacency list and no further query
 * operations are needed when the neighbors of the same vertex are
 * queried again.
 * \param al The lazy adjacency list.
 * \param no The vertex id to query.
 * \return Pointer to a vector. It is allowed to modify it and
 *   modification does not affect the original graph.
 * 
 * Time complexity: O(d), the number of neighbor vertices for the
 * first time, O(1) for subsequent calls.
 */
#define igraph_lazy_adjlist_get(al,no) \
  ((al)->adjs[(long int)(no)] != 0 ? ((al)->adjs[(long int)(no)]) : \
   (igraph_lazy_adjlist_get_real(al, no)))
igraph_vector_t *igraph_lazy_adjlist_get_real(igraph_lazy_adjlist_t *al,
						igraph_integer_t no);

typedef struct igraph_lazy_adjedgelist_t {
  const igraph_t *graph;
  igraph_integer_t length;
  igraph_vector_t **adjs;
  igraph_neimode_t mode;
} igraph_lazy_adjedgelist_t;

int igraph_lazy_adjedgelist_init(const igraph_t *graph,
				   igraph_lazy_adjedgelist_t *al,
				   igraph_neimode_t mode);
void igraph_lazy_adjedgelist_destroy(igraph_lazy_adjedgelist_t *al);
/**
 * \define igraph_lazy_adjedgelist_get
 * Query adjacent edges
 * 
 * If the function is called for the first time for a vertex, then the
 * result is stored in the adjacency list and no further query
 * operations are needed when the adjacent edges of the same vertex are
 * queried again.
 * \param al The lazy adjacency list object.
 * \param no The vertex id to query.
 * \return Pointer to a vector. It is allowed to modify it and
 *   modification does not affect the original graph.
 * 
 * Time complexity: O(d), the number of adjacent edges for the first
 * time, O(1) for subsequent calls with the same \p no argument.
 */
#define igraph_lazy_adjedgelist_get(al,no) \
  ((al)->adjs[(long int)(no)] != 0 ? ((al)->adjs[(long int)(no)]) : \
   (igraph_lazy_adjedgelist_get_real(al, no)))
igraph_vector_t *igraph_lazy_adjedgelist_get_real(igraph_lazy_adjedgelist_t *al,
						    igraph_integer_t no);

__END_DECLS

#endif
