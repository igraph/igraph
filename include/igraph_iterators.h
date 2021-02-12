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

#ifndef IGRAPH_ITERATORS_H
#define IGRAPH_ITERATORS_H

#include "igraph_decls.h"
#include "igraph_constants.h"

__BEGIN_DECLS

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
        igraph_integer_t vid;               /* single vertex  */
        const igraph_vector_t *vecptr;      /* vector of vertices  */
        struct {
            igraph_integer_t vid;
            igraph_neimode_t mode;
        } adj;                  /* adjacent vertices  */
        struct {
            igraph_integer_t from;
            igraph_integer_t to;
        } seq;                              /* sequence of vertices from:to */
    } data;
} igraph_vs_t;

IGRAPH_EXPORT int igraph_vs_all(igraph_vs_t *vs);
IGRAPH_EXPORT igraph_vs_t igraph_vss_all(void);

IGRAPH_EXPORT int igraph_vs_adj(igraph_vs_t *vs,
                                igraph_integer_t vid, igraph_neimode_t mode);

IGRAPH_EXPORT int igraph_vs_nonadj(igraph_vs_t *vs, igraph_integer_t vid,
                                   igraph_neimode_t mode);

IGRAPH_EXPORT int igraph_vs_none(igraph_vs_t *vs);
IGRAPH_EXPORT igraph_vs_t igraph_vss_none(void);

IGRAPH_EXPORT int igraph_vs_1(igraph_vs_t *vs, igraph_integer_t vid);
IGRAPH_EXPORT igraph_vs_t igraph_vss_1(igraph_integer_t vid);

IGRAPH_EXPORT int igraph_vs_vector(igraph_vs_t *vs,
                                   const igraph_vector_t *v);
IGRAPH_EXPORT igraph_vs_t igraph_vss_vector(const igraph_vector_t *v);

IGRAPH_EXPORT int igraph_vs_vector_small(igraph_vs_t *vs, ...);

IGRAPH_EXPORT int igraph_vs_vector_copy(igraph_vs_t *vs,
                                        const igraph_vector_t *v);

IGRAPH_EXPORT int igraph_vs_seq(igraph_vs_t *vs, igraph_integer_t from, igraph_integer_t to);
IGRAPH_EXPORT igraph_vs_t igraph_vss_seq(igraph_integer_t from, igraph_integer_t to);

IGRAPH_EXPORT void igraph_vs_destroy(igraph_vs_t *vs);

IGRAPH_EXPORT igraph_bool_t igraph_vs_is_all(const igraph_vs_t *vs);

IGRAPH_EXPORT int igraph_vs_copy(igraph_vs_t* dest, const igraph_vs_t* src);

IGRAPH_EXPORT int igraph_vs_as_vector(const igraph_t *graph, igraph_vs_t vs,
                                      igraph_vector_t *v);
IGRAPH_EXPORT int igraph_vs_size(const igraph_t *graph, const igraph_vs_t *vs,
                                 igraph_integer_t *result);
IGRAPH_EXPORT int igraph_vs_type(const igraph_vs_t *vs);

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
 * pointed to by the iterator (call this only if \ref IGRAPH_VIT_END()
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
 * while (!IGRAPH_VIT_END(vit)) {
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
 * \brief Next vertex.
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
 * \brief Are we at the end?
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
 * \brief Size of a vertex iterator.
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
 * \brief Reset a vertex iterator.
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
 * \brief Query the current position.
 *
 * Gives the vertex id of the current vertex pointed to by the
 * iterator.
 * \param vit The vertex iterator.
 * \return The vertex id of the current vertex.
 *
 * Time complexity: O(1).
 */
#define IGRAPH_VIT_GET(vit)  \
    ((igraph_integer_t)(((vit).type == IGRAPH_VIT_SEQ) ? (vit).pos : \
                        VECTOR(*(vit).vec)[(vit).pos]))

IGRAPH_EXPORT int igraph_vit_create(const igraph_t *graph,
                                    igraph_vs_t vs, igraph_vit_t *vit);
IGRAPH_EXPORT void igraph_vit_destroy(const igraph_vit_t *vit);

IGRAPH_EXPORT int igraph_vit_as_vector(const igraph_vit_t *vit, igraph_vector_t *v);

/* -------------------------------------------------- */
/* Edge Selectors                                     */
/* -------------------------------------------------- */

#define IGRAPH_ES_ALL       0
#define IGRAPH_ES_ALLFROM   1
#define IGRAPH_ES_ALLTO     2
#define IGRAPH_ES_INCIDENT  3
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
        } incident;
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

IGRAPH_EXPORT int igraph_es_all(igraph_es_t *es,
                                igraph_edgeorder_type_t order);
IGRAPH_EXPORT igraph_es_t igraph_ess_all(igraph_edgeorder_type_t order);

IGRAPH_EXPORT int igraph_es_incident(igraph_es_t *es,
                                     igraph_integer_t vid, igraph_neimode_t mode);

IGRAPH_EXPORT int igraph_es_none(igraph_es_t *es);
IGRAPH_EXPORT igraph_es_t igraph_ess_none(void);

IGRAPH_EXPORT int igraph_es_1(igraph_es_t *es, igraph_integer_t eid);
IGRAPH_EXPORT igraph_es_t igraph_ess_1(igraph_integer_t eid);

IGRAPH_EXPORT int igraph_es_vector(igraph_es_t *es,
                                   const igraph_vector_t *v);
IGRAPH_EXPORT igraph_es_t igraph_ess_vector(const igraph_vector_t *v);

IGRAPH_EXPORT int igraph_es_fromto(igraph_es_t *es,
                                   igraph_vs_t from, igraph_vs_t to);

IGRAPH_EXPORT int igraph_es_seq(igraph_es_t *es, igraph_integer_t from, igraph_integer_t to);
IGRAPH_EXPORT igraph_es_t igraph_ess_seq(igraph_integer_t from, igraph_integer_t to);

IGRAPH_EXPORT int igraph_es_vector_copy(igraph_es_t *es, const igraph_vector_t *v);

IGRAPH_EXPORT int igraph_es_pairs(igraph_es_t *es, const igraph_vector_t *v,
                                  igraph_bool_t directed);
IGRAPH_EXPORT int igraph_es_pairs_small(igraph_es_t *es, igraph_bool_t directed, ...);

IGRAPH_EXPORT int igraph_es_multipairs(igraph_es_t *es, const igraph_vector_t *v,
                                       igraph_bool_t directed);

IGRAPH_EXPORT int igraph_es_path(igraph_es_t *es, const igraph_vector_t *v,
                                 igraph_bool_t directed);
IGRAPH_EXPORT int igraph_es_path_small(igraph_es_t *es, igraph_bool_t directed, ...);

IGRAPH_EXPORT void igraph_es_destroy(igraph_es_t *es);

IGRAPH_EXPORT igraph_bool_t igraph_es_is_all(const igraph_es_t *es);

IGRAPH_EXPORT int igraph_es_copy(igraph_es_t* dest, const igraph_es_t* src);

IGRAPH_EXPORT int igraph_es_as_vector(const igraph_t *graph, igraph_es_t es,
                                      igraph_vector_t *v);
IGRAPH_EXPORT int igraph_es_size(const igraph_t *graph, const igraph_es_t *es,
                                 igraph_integer_t *result);
IGRAPH_EXPORT int igraph_es_type(const igraph_es_t *es);


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
 * \brief Next edge.
 *
 * Steps the iterator to the next edge. Call this function only if
 * \ref IGRAPH_EIT_END() returns false.
 * \param eit The edge iterator to step.
 *
 * Time complexity: O(1).
 */
#define IGRAPH_EIT_NEXT(eit) (++((eit).pos))
/**
 * \define IGRAPH_EIT_END
 * \brief Are we at the end?
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
 * \brief Number of edges in the iterator.
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
 * \brief Reset an edge iterator.
 *
 * Resets an edge iterator. After calling this macro the iterator will
 * point to the first edge.
 * \param eit The edge iterator.
 *
 * Time complexity: O(1).
 */
#define IGRAPH_EIT_RESET(eit) ((eit).pos = (eit).start)
/**
 * \define IGRAPH_EIT_GET
 * \brief Query an edge iterator.
 *
 * Gives the edge id of the current edge pointed to by an iterator.
 * \param eit The edge iterator.
 * \return The id of the current edge.
 *
 * Time complexity: O(1).
 */
#define IGRAPH_EIT_GET(eit)  \
    (igraph_integer_t)((((eit).type == IGRAPH_EIT_SEQ) ? (eit).pos : \
                        VECTOR(*(eit).vec)[(eit).pos]))

IGRAPH_EXPORT int igraph_eit_create(const igraph_t *graph,
                                    igraph_es_t es, igraph_eit_t *eit);
IGRAPH_EXPORT void igraph_eit_destroy(const igraph_eit_t *eit);

IGRAPH_EXPORT int igraph_eit_as_vector(const igraph_eit_t *eit, igraph_vector_t *v);

__END_DECLS

#endif
