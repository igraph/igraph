/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2009-2012  Gabor Csardi <csardi.gabor@gmail.com>
   334 Harvard street, Cambridge, MA 02139 USA

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

#include "igraph_decls.h"
#include "igraph_constants.h"
#include "igraph_types.h"
#include "igraph_datatype.h"

__BEGIN_DECLS

typedef struct igraph_adjlist_t {
    igraph_integer_t length;
    igraph_vector_int_t *adjs;
} igraph_adjlist_t;

DECLDIR int igraph_adjlist_init(const igraph_t *graph, igraph_adjlist_t *al,
                                igraph_neimode_t mode);
DECLDIR int igraph_adjlist_init_empty(igraph_adjlist_t *al, igraph_integer_t no_of_nodes);
DECLDIR igraph_integer_t igraph_adjlist_size(const igraph_adjlist_t *al);
DECLDIR int igraph_adjlist_init_complementer(const igraph_t *graph,
        igraph_adjlist_t *al,
        igraph_neimode_t mode,
        igraph_bool_t loops);
DECLDIR void igraph_adjlist_destroy(igraph_adjlist_t *al);
DECLDIR void igraph_adjlist_clear(igraph_adjlist_t *al);
DECLDIR void igraph_adjlist_sort(igraph_adjlist_t *al);
DECLDIR int igraph_adjlist_simplify(igraph_adjlist_t *al);
DECLDIR int igraph_adjlist_remove_duplicate(const igraph_t *graph,
        igraph_adjlist_t *al);
DECLDIR int igraph_adjlist_print(const igraph_adjlist_t *al);
DECLDIR int igraph_adjlist_fprint(const igraph_adjlist_t *al, FILE *outfile);
DECLDIR igraph_bool_t igraph_adjlist_has_edge(igraph_adjlist_t* al, igraph_integer_t from, igraph_integer_t to, igraph_bool_t directed);
DECLDIR int igraph_adjlist_replace_edge(igraph_adjlist_t* al, igraph_integer_t from, igraph_integer_t oldto, igraph_integer_t newto, igraph_bool_t directed);

/* igraph_vector_int_t *igraph_adjlist_get(const igraph_adjlist_t *al,  */
/*                 igraph_integer_t no); */
/**
 * \define igraph_adjlist_get
 * Query a vector in an adjlist
 *
 * Returns a pointer to an <type>igraph_vector_int_t</type> object from an
 * adjacency list. The vector can be modified as desired.
 * \param al The adjacency list object.
 * \param no The vertex of which the vertex of adjacent vertices are
 *   returned.
 * \return Pointer to the <type>igraph_vector_int_t</type> object.
 *
 * Time complexity: O(1).
 */
#define igraph_adjlist_get(al,no) (&(al)->adjs[(long int)(no)])

DECLDIR int igraph_adjlist(igraph_t *graph, const igraph_adjlist_t *adjlist,
                           igraph_neimode_t mode, igraph_bool_t duplicate);

typedef struct igraph_inclist_t {
    igraph_integer_t length;
    igraph_vector_int_t *incs;
} igraph_inclist_t;

DECLDIR int igraph_inclist_init(const igraph_t *graph,
                                igraph_inclist_t *il,
                                igraph_neimode_t mode);
DECLDIR int igraph_inclist_init_empty(igraph_inclist_t *il, igraph_integer_t n);
DECLDIR void igraph_inclist_destroy(igraph_inclist_t *il);
DECLDIR void igraph_inclist_clear(igraph_inclist_t *il);
DECLDIR int igraph_inclist_remove_duplicate(const igraph_t *graph,
        igraph_inclist_t *il);
DECLDIR int igraph_inclist_print(const igraph_inclist_t *il);
DECLDIR int igraph_inclist_fprint(const igraph_inclist_t *il, FILE *outfile);

/**
 * \define igraph_inclist_get
 * Query a vector in an incidence list
 *
 * Returns a pointer to an <type>igraph_vector_int_t</type> object from an
 * incidence list containing edge ids. The vector can be modified,
 * resized, etc. as desired.
 * \param il Pointer to the incidence list.
 * \param no The vertex for which the incident edges are returned.
 * \return Pointer to an <type>igraph_vector_int_t</type> object.
 *
 * Time complexity: O(1).
 */
#define igraph_inclist_get(il,no) (&(il)->incs[(long int)(no)])

typedef struct igraph_lazy_adjlist_t {
    const igraph_t *graph;
    igraph_integer_t length;
    igraph_vector_t **adjs;
    igraph_neimode_t mode;
    igraph_lazy_adlist_simplify_t simplify;
} igraph_lazy_adjlist_t;

DECLDIR int igraph_lazy_adjlist_init(const igraph_t *graph,
                                     igraph_lazy_adjlist_t *al,
                                     igraph_neimode_t mode,
                                     igraph_lazy_adlist_simplify_t simplify);
DECLDIR void igraph_lazy_adjlist_destroy(igraph_lazy_adjlist_t *al);
DECLDIR void igraph_lazy_adjlist_clear(igraph_lazy_adjlist_t *al);
/* igraph_vector_t *igraph_lazy_adjlist_get(igraph_lazy_adjlist_t *al, */
/*                     igraph_integer_t no); */
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
DECLDIR igraph_vector_t *igraph_lazy_adjlist_get_real(igraph_lazy_adjlist_t *al,
        igraph_integer_t no);

typedef struct igraph_lazy_inclist_t {
    const igraph_t *graph;
    igraph_integer_t length;
    igraph_vector_t **incs;
    igraph_neimode_t mode;
} igraph_lazy_inclist_t;

DECLDIR int igraph_lazy_inclist_init(const igraph_t *graph,
                                     igraph_lazy_inclist_t *il,
                                     igraph_neimode_t mode);
DECLDIR void igraph_lazy_inclist_destroy(igraph_lazy_inclist_t *il);
DECLDIR void igraph_lazy_inclist_clear(igraph_lazy_inclist_t *il);

/**
 * \define igraph_lazy_inclist_get
 * Query incident edges
 *
 * If the function is called for the first time for a vertex, then the
 * result is stored in the incidence list and no further query
 * operations are needed when the incident edges of the same vertex are
 * queried again.
 * \param al The lazy incidence list object.
 * \param no The vertex id to query.
 * \return Pointer to a vector. It is allowed to modify it and
 *   modification does not affect the original graph.
 *
 * Time complexity: O(d), the number of incident edges for the first
 * time, O(1) for subsequent calls with the same \p no argument.
 */
#define igraph_lazy_inclist_get(al,no) \
    ((al)->incs[(long int)(no)] != 0 ? ((al)->incs[(long int)(no)]) : \
     (igraph_lazy_inclist_get_real(al, no)))
DECLDIR igraph_vector_t *igraph_lazy_inclist_get_real(igraph_lazy_inclist_t *al,
        igraph_integer_t no);

/*************************************************************************
 * DEPRECATED TYPES AND FUNCTIONS
 */

typedef igraph_inclist_t igraph_adjedgelist_t;

DECLDIR int igraph_adjedgelist_init(const igraph_t *graph,
                                    igraph_inclist_t *il,
                                    igraph_neimode_t mode);
DECLDIR void igraph_adjedgelist_destroy(igraph_inclist_t *il);
DECLDIR int igraph_adjedgelist_remove_duplicate(const igraph_t *graph,
        igraph_inclist_t *il);
DECLDIR int igraph_adjedgelist_print(const igraph_inclist_t *il, FILE *outfile);

/**
 * \define igraph_adjedgelist_get
 * Query a vector in an incidence list
 *
 * This macro was superseded by \ref igraph_inclist_get() in igraph 0.6.
 * Please use \ref igraph_inclist_get() instead of this macro.
 *
 * </para><para>
 * Deprecated in version 0.6.
 */
#define igraph_adjedgelist_get(ael,no) (&(ael)->incs[(long int)(no)])

typedef igraph_lazy_inclist_t igraph_lazy_adjedgelist_t;

DECLDIR int igraph_lazy_adjedgelist_init(const igraph_t *graph,
        igraph_lazy_inclist_t *il,
        igraph_neimode_t mode);
DECLDIR void igraph_lazy_adjedgelist_destroy(igraph_lazy_inclist_t *il);

/**
 * \define igraph_lazy_adjedgelist_get
 * Query a vector in a lazy incidence list
 *
 * This macro was superseded by \ref igraph_lazy_inclist_get() in igraph 0.6.
 * Please use \ref igraph_lazy_inclist_get() instead of this macro.
 *
 * </para><para>
 * Deprecated in version 0.6.
 */
#define igraph_lazy_adjedgelist_get(al,no) \
    ((al)->incs[(long int)(no)] != 0 ? ((al)->incs[(long int)(no)]) : \
     (igraph_lazy_adjedgelist_get_real(al, no)))
DECLDIR igraph_vector_t *igraph_lazy_adjedgelist_get_real(igraph_lazy_inclist_t *al,
        igraph_integer_t no);
__END_DECLS

#endif
