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

IGRAPH_EXPORT int igraph_adjlist_init(const igraph_t *graph, igraph_adjlist_t *al,
                                      igraph_neimode_t mode, igraph_loops_t loops,
                                      igraph_multiple_t multiple);
IGRAPH_EXPORT int igraph_adjlist_init_empty(igraph_adjlist_t *al, igraph_integer_t no_of_nodes);
IGRAPH_EXPORT igraph_integer_t igraph_adjlist_size(const igraph_adjlist_t *al);
IGRAPH_EXPORT int igraph_adjlist_init_complementer(const igraph_t *graph,
                                                   igraph_adjlist_t *al,
                                                   igraph_neimode_t mode,
                                                   igraph_bool_t loops);
IGRAPH_EXPORT void igraph_adjlist_destroy(igraph_adjlist_t *al);
IGRAPH_EXPORT void igraph_adjlist_clear(igraph_adjlist_t *al);
IGRAPH_EXPORT void igraph_adjlist_sort(igraph_adjlist_t *al);
IGRAPH_EXPORT int igraph_adjlist_simplify(igraph_adjlist_t *al);
IGRAPH_DEPRECATED IGRAPH_EXPORT int igraph_adjlist_remove_duplicate(const igraph_t *graph,
                                                                    igraph_adjlist_t *al);
IGRAPH_EXPORT int igraph_adjlist_print(const igraph_adjlist_t *al);
IGRAPH_EXPORT int igraph_adjlist_fprint(const igraph_adjlist_t *al, FILE *outfile);
IGRAPH_EXPORT igraph_bool_t igraph_adjlist_has_edge(igraph_adjlist_t* al, igraph_integer_t from, igraph_integer_t to, igraph_bool_t directed);
IGRAPH_EXPORT int igraph_adjlist_replace_edge(igraph_adjlist_t* al, igraph_integer_t from, igraph_integer_t oldto, igraph_integer_t newto, igraph_bool_t directed);

/**
 * \define igraph_adjlist_get
 * \brief Query a vector in an adjacency list.
 *
 * Returns a pointer to an <type>igraph_vector_int_t</type> object from an
 * adjacency list. The vector can be modified as desired.
 * \param al The adjacency list object.
 * \param no The vertex whose adjacent vertices will be returned.
 * \return Pointer to the <type>igraph_vector_int_t</type> object.
 *
 * Time complexity: O(1).
 */
#define igraph_adjlist_get(al,no) (&(al)->adjs[(long int)(no)])

IGRAPH_EXPORT int igraph_adjlist(igraph_t *graph, const igraph_adjlist_t *adjlist,
                                 igraph_neimode_t mode, igraph_bool_t duplicate);

typedef struct igraph_inclist_t {
    igraph_integer_t length;
    igraph_vector_int_t *incs;
} igraph_inclist_t;

IGRAPH_EXPORT int igraph_inclist_init(const igraph_t *graph,
                                      igraph_inclist_t *il,
                                      igraph_neimode_t mode,
                                      igraph_loops_t loops);
IGRAPH_EXPORT int igraph_inclist_init_empty(igraph_inclist_t *il, igraph_integer_t n);
IGRAPH_EXPORT igraph_integer_t igraph_inclist_size(const igraph_inclist_t *al);
IGRAPH_EXPORT void igraph_inclist_destroy(igraph_inclist_t *il);
IGRAPH_EXPORT void igraph_inclist_clear(igraph_inclist_t *il);
IGRAPH_DEPRECATED IGRAPH_EXPORT int igraph_inclist_remove_duplicate(const igraph_t *graph,
                                                                    igraph_inclist_t *il);
IGRAPH_EXPORT int igraph_inclist_print(const igraph_inclist_t *il);
IGRAPH_EXPORT int igraph_inclist_fprint(const igraph_inclist_t *il, FILE *outfile);

/**
 * \define igraph_inclist_get
 * \brief Query a vector in an incidence list.
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
    igraph_vector_int_t **adjs;
    igraph_neimode_t mode;
    igraph_loops_t loops;
    igraph_multiple_t multiple;
    igraph_vector_t dummy;
} igraph_lazy_adjlist_t;

IGRAPH_EXPORT int igraph_lazy_adjlist_init(const igraph_t *graph,
                                           igraph_lazy_adjlist_t *al,
                                           igraph_neimode_t mode,
                                           igraph_loops_t loops,
                                           igraph_multiple_t multiple);
IGRAPH_EXPORT void igraph_lazy_adjlist_destroy(igraph_lazy_adjlist_t *al);
IGRAPH_EXPORT void igraph_lazy_adjlist_clear(igraph_lazy_adjlist_t *al);
IGRAPH_EXPORT igraph_integer_t igraph_lazy_adjlist_size(const igraph_lazy_adjlist_t *al);

/**
 * \define igraph_lazy_adjlist_get
 * \brief Query neighbor vertices.
 *
 * If the function is called for the first time for a vertex then the
 * result is stored in the adjacency list and no further query
 * operations are needed when the neighbors of the same vertex are
 * queried again.
 * \param al The lazy adjacency list.
 * \param no The vertex ID to query.
 * \return Pointer to a vector. It is allowed to modify it and
 *   modification does not affect the original graph.
 *
 * Time complexity: O(d), the number of neighbor vertices for the
 * first time, O(1) for subsequent calls.
 */
#define igraph_lazy_adjlist_get(al,no) \
    ((al)->adjs[(long int)(no)] != 0 ? ((al)->adjs[(long int)(no)]) : \
     (igraph_i_lazy_adjlist_get_real(al, no)))
IGRAPH_EXPORT igraph_vector_int_t *igraph_i_lazy_adjlist_get_real(igraph_lazy_adjlist_t *al, igraph_integer_t no);

typedef struct igraph_lazy_inclist_t {
    const igraph_t *graph;
    igraph_integer_t length;
    igraph_vector_int_t **incs;
    igraph_neimode_t mode;
    igraph_vector_t dummy;
    igraph_loops_t loops;
} igraph_lazy_inclist_t;

IGRAPH_EXPORT int igraph_lazy_inclist_init(const igraph_t *graph,
                                           igraph_lazy_inclist_t *il,
                                           igraph_neimode_t mode,
                                           igraph_loops_t loops);
IGRAPH_EXPORT void igraph_lazy_inclist_destroy(igraph_lazy_inclist_t *il);
IGRAPH_EXPORT void igraph_lazy_inclist_clear(igraph_lazy_inclist_t *il);
IGRAPH_EXPORT igraph_integer_t igraph_lazy_inclist_size(const igraph_lazy_inclist_t *il);

/**
 * \define igraph_lazy_inclist_get
 * \brief Query incident edges.
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
     (igraph_i_lazy_inclist_get_real(al, no)))
IGRAPH_EXPORT igraph_vector_int_t *igraph_i_lazy_inclist_get_real(igraph_lazy_inclist_t *al, igraph_integer_t no);

__END_DECLS

#endif
