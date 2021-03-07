/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2005-2012  Gabor Csardi <csardi.gabor@gmail.com>
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
   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
   02110-1301 USA

*/

#include "igraph_iterators.h"
#include "igraph_memory.h"
#include "igraph_interface.h"

#include <string.h>
#include <stdarg.h>

/**
 * \section about_iterators About selectors, iterators
 *
 * <para>Everything about vertices and vertex selectors also applies
 * to edges and edge selectors unless explicitly noted otherwise.</para>
 *
 * <para>The vertex (and edge) selector notion was introduced in igraph 0.2.
 * It is a way to reference a sequence of vertices or edges
 * independently of the graph.</para>
 *
 * <para>While this might sound quite mysterious, it is actually very
 * simple. For example, all vertices of a graph can be selected by
 * \ref igraph_vs_all() and the graph independence means that
 * \ref igraph_vs_all() is not parametrized by a graph object. That is,
 * \ref igraph_vs_all() is the general \em concept of selecting all vertices
 * of a graph. A vertex selector is then a way to specify the class of vertices
 * to be visited. The selector might specify that all vertices of a graph or
 * all the neighbours of a vertex are to be visited. A vertex selector is a
 * way of saying that you want to visit a bunch of vertices, as opposed to a
 * vertex iterator which is a concrete plan for visiting each of the
 * chosen vertices of a specific graph.</para>
 *
 * <para>To determine the actual vertex IDs implied by a vertex selector, you
 * need to apply the concept of selecting vertices to a specific graph object.
 * This can be accomplished by instantiating a vertex iterator using a
 * specific vertex selection concept and a specific graph object. The notion
 * of vertex iterators can be thought of in the following way. Given a
 * specific graph object and the class of vertices to be visited, a vertex
 * iterator is a road map, plan or route for how to visit the chosen
 * vertices.</para>
 *
 * <para>Some vertex selectors have \em immediate versions. These have the
 * prefix \c igraph_vss instead of \c igraph_vs, e.g. \ref igraph_vss_all()
 * instead of \ref igraph_vs_all(). The immediate versions are to be used in
 * the parameter list of the igraph functions, such as \ref igraph_degree().
 * These functions are not associated with any \type igraph_vs_t object, so
 * they have no separate constructors and destructors
 * (destroy functions).</para>
 */

/**
 * \section about_vertex_selectors
 *
 * <para>Vertex selectors are created by vertex selector constructors,
 * can be instantiated with \ref igraph_vit_create(), and are
 * destroyed with \ref igraph_vs_destroy().</para>
 */

/**
 * \function igraph_vs_all
 * \brief Vertex set, all vertices of a graph.
 *
 * \param vs Pointer to an uninitialized \type igraph_vs_t object.
 * \return Error code.
 * \sa \ref igraph_vss_all(), \ref igraph_vs_destroy()
 *
 * This selector includes all vertices of a given graph in
 * increasing vertex id order.
 *
 * </para><para>
 * Time complexity: O(1).
 */

int igraph_vs_all(igraph_vs_t *vs) {
    vs->type = IGRAPH_VS_ALL;
    return 0;
}

/**
 * \function igraph_vss_all
 * \brief All vertices of a graph (immediate version).
 *
 * Immediate vertex selector for all vertices in a graph. It can
 * be used conveniently when some vertex property (e.g. betweenness,
 * degree, etc.) should be calculated for all vertices.
 *
 * \return A vertex selector for all vertices in a graph.
 * \sa \ref igraph_vs_all()
 *
 * Time complexity: O(1).
 */

igraph_vs_t igraph_vss_all(void) {
    igraph_vs_t allvs;
    allvs.type = IGRAPH_VS_ALL;
    return allvs;
}

/**
 * \function igraph_vs_adj
 * \brief Adjacent vertices of a vertex.
 *
 * All neighboring vertices of a given vertex are selected by this
 * selector. The \c mode argument controls the type of the neighboring
 * vertices to be selected. The vertices are visited in increasing vertex
 * ID order, as of igraph version 0.4.
 *
 * \param vs Pointer to an uninitialized vertex selector object.
 * \param vid Vertex ID, the center of the neighborhood.
 * \param mode Decides the type of the neighborhood for directed
 *        graphs. This parameter is ignored for undirected graphs.
 *        Possible values:
 *        \clist
 *        \cli IGRAPH_OUT
 *          All vertices to which there is a directed edge from \c vid. That
 *          is, all the out-neighbors of \c vid.
 *        \cli IGRAPH_IN
 *          All vertices from which there is a directed edge to \c vid. In
 *          other words, all the in-neighbors of \c vid.
 *        \cli IGRAPH_ALL
 *          All vertices to which or from which there is a directed edge
 *          from/to \c vid. That is, all the neighbors of \c vid considered
 *          as if the graph is undirected.
 *        \endclist
 * \return Error code.
 * \sa \ref igraph_vs_destroy()
 *
 * Time complexity: O(1).
 */

int igraph_vs_adj(igraph_vs_t *vs,
                  igraph_integer_t vid, igraph_neimode_t mode) {
    vs->type = IGRAPH_VS_ADJ;
    vs->data.adj.vid = vid;
    vs->data.adj.mode = mode;
    return 0;
}

/**
 * \function igraph_vs_nonadj
 * \brief Non-adjacent vertices of a vertex.
 *
 * All non-neighboring vertices of a given vertex. The \p mode
 * argument controls the type of neighboring vertices \em not to
 * select. Instead of selecting immediate neighbors of \c vid as is done by
 * \ref igraph_vs_adj(), the current function selects vertices that are \em not
 * immediate neighbors of \c vid.
 *
 * \param vs Pointer to an uninitialized vertex selector object.
 * \param vid Vertex ID, the \quote center \endquote of the
 *        non-neighborhood.
 * \param mode The type of neighborhood not to select in directed
 *        graphs. Possible values:
 *        \clist
 *        \cli IGRAPH_OUT
 *          All vertices will be selected except those to which there is a
 *          directed edge from \c vid. That is, we select all vertices
 *          excluding the out-neighbors of \c vid.
 *        \cli IGRAPH_IN
 *          All vertices will be selected except those from which there is a
 *          directed edge to \c vid. In other words, we select all vertices
 *          but the in-neighbors of \c vid.
 *        \cli IGRAPH_ALL
 *          All vertices will be selected except those from or to which there
 *          is a directed edge to or from \c vid. That is, we select all
 *          vertices of \c vid except for its immediate neighbors.
 *        \endclist
 * \return Error code.
 * \sa \ref igraph_vs_destroy()
 *
 * Time complexity: O(1).
 *
 * \example examples/simple/igraph_vs_nonadj.c
 */

int igraph_vs_nonadj(igraph_vs_t *vs, igraph_integer_t vid,
                     igraph_neimode_t mode) {
    vs->type = IGRAPH_VS_NONADJ;
    vs->data.adj.vid = vid;
    vs->data.adj.mode = mode;
    return 0;
}

/**
 * \function igraph_vs_none
 * \brief Empty vertex set.
 *
 * Creates an empty vertex selector.
 *
 * \param vs Pointer to an uninitialized vertex selector object.
 * \return Error code.
 * \sa \ref igraph_vss_none(), \ref igraph_vs_destroy()
 *
 * Time complexity: O(1).
 */

int igraph_vs_none(igraph_vs_t *vs) {
    vs->type = IGRAPH_VS_NONE;
    return 0;
}

/**
 * \function igraph_vss_none
 * \brief Empty vertex set (immediate version).
 *
 * The immediate version of the empty vertex selector.
 *
 * \return An empty vertex selector.
 * \sa \ref igraph_vs_none()
 *
 * Time complexity: O(1).
 */

igraph_vs_t igraph_vss_none(void) {
    igraph_vs_t nonevs;
    nonevs.type = IGRAPH_VS_NONE;
    return nonevs;
}

/**
 * \function igraph_vs_1
 * \brief Vertex set with a single vertex.
 *
 * This vertex selector selects a single vertex.
 *
 * \param vs Pointer to an uninitialized vertex selector object.
 * \param vid The vertex id to be selected.
 * \return Error Code.
 * \sa \ref igraph_vss_1(), \ref igraph_vs_destroy()
 *
 * Time complexity: O(1).
 */

int igraph_vs_1(igraph_vs_t *vs, igraph_integer_t vid) {
    vs->type = IGRAPH_VS_1;
    vs->data.vid = vid;
    return 0;
}

/**
 * \function igraph_vss_1
 * \brief Vertex set with a single vertex (immediate version).
 *
 * The immediate version of the single-vertex selector.
 *
 * \param vid The vertex to be selected.
 * \return A vertex selector containing a single vertex.
 * \sa \ref igraph_vs_1()
 *
 * Time complexity: O(1).
 */

igraph_vs_t igraph_vss_1(igraph_integer_t vid) {
    igraph_vs_t onevs;
    onevs.type = IGRAPH_VS_1;
    onevs.data.vid = vid;
    return onevs;
}

/**
 * \function igraph_vs_vector
 * \brief Vertex set based on a vector.
 *
 * This function makes it possible to handle a \type vector_t
 * temporarily as a vertex selector. The vertex selector should be
 * thought of like a \em view to the vector. If you make changes to
 * the vector that also affects the vertex selector. Destroying the
 * vertex selector does not destroy the vector. (Of course.) Do not
 * destroy the vector before destroying the vertex selector, or you
 * might get strange behavior.
 *
 * \param vs Pointer to an uninitialized vertex selector.
 * \param v Pointer to a \type igraph_vector_t object.
 * \return Error code.
 * \sa \ref igraph_vss_vector(), \ref igraph_vs_destroy()
 *
 * Time complexity: O(1).
 *
 * \example examples/simple/igraph_vs_vector.c
 */

int igraph_vs_vector(igraph_vs_t *vs,
                     const igraph_vector_t *v) {
    vs->type = IGRAPH_VS_VECTORPTR;
    vs->data.vecptr = v;
    return 0;
}

/**
 * \function igraph_vss_vector
 * \brief Vertex set based on a vector (immediate version).
 *
 * This is the immediate version of \ref igraph_vs_vector.
 *
 * \param v Pointer to a \type igraph_vector_t object.
 * \return A vertex selector object containing the vertices in the
 *         vector.
 * \sa \ref igraph_vs_vector()
 *
 * Time complexity: O(1).
 */

igraph_vs_t igraph_vss_vector(const igraph_vector_t *v) {
    igraph_vs_t vecvs;
    vecvs.type = IGRAPH_VS_VECTORPTR;
    vecvs.data.vecptr = v;
    return vecvs;
}

/**
 * \function igraph_vs_vector_small
 * \brief Create a vertex set by giving its elements.
 *
 * This function can be used to create a vertex selector with a couple
 * of vertices. Do not forget to include a <code>-1</code> after the
 * last vertex id. The behavior of the function is undefined if you
 * don't use a <code>-1</code> properly.
 *
 * </para><para>
 * Note that the vertex ids supplied will be parsed as
 * <code>int</code>'s so you cannot supply arbitrarily large (too
 * large for int) vertex ids here.
 *
 * \param vs Pointer to an uninitialized vertex selector object.
 * \param ... Additional parameters, these will be the vertex ids to
 *        be included in the vertex selector. Supply a <code>-1</code>
 *        after the last vertex id.
 * \return Error code.
 * \sa \ref igraph_vs_destroy()
 *
 * Time complexity: O(n), the number of vertex ids supplied.
 */

int igraph_vs_vector_small(igraph_vs_t *vs, ...) {
    va_list ap;
    long int i, n = 0;
    vs->type = IGRAPH_VS_VECTOR;
    vs->data.vecptr = IGRAPH_CALLOC(1, igraph_vector_t);
    if (vs->data.vecptr == 0) {
        IGRAPH_ERROR("Cannot create vertex selector", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, (igraph_vector_t*)vs->data.vecptr);

    va_start(ap, vs);
    while (1) {
        int num = va_arg(ap, int);
        if (num == -1) {
            break;
        }
        n++;
    }
    va_end(ap);

    IGRAPH_VECTOR_INIT_FINALLY((igraph_vector_t*)vs->data.vecptr, n);

    va_start(ap, vs);
    for (i = 0; i < n; i++) {
        VECTOR(*vs->data.vecptr)[i] = (igraph_real_t) va_arg(ap, int);
    }
    va_end(ap);

    IGRAPH_FINALLY_CLEAN(2);
    return 0;
}

/**
 * \function igraph_vs_vector_copy
 * \brief Vertex set based on a vector, with copying.
 *
 * This function makes it possible to handle a \type vector_t
 * permanently as a vertex selector. The vertex selector creates a
 * copy of the original vector, so the vector can safely be destroyed
 * after creating the vertex selector. Changing the original vector
 * will not affect the vertex selector. The vertex selector is
 * responsible for deleting the copy made by itself.
 *
 * \param vs Pointer to an uninitialized vertex selector.
 * \param v Pointer to a \type igraph_vector_t object.
 * \return Error code.
 * \sa \ref igraph_vs_destroy()
 *
 * Time complexity: O(1).
 */

int igraph_vs_vector_copy(igraph_vs_t *vs,
                          const igraph_vector_t *v) {
    vs->type = IGRAPH_VS_VECTOR;
    vs->data.vecptr = IGRAPH_CALLOC(1, igraph_vector_t);
    if (vs->data.vecptr == 0) {
        IGRAPH_ERROR("Cannot create vertex selector", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, (igraph_vector_t*)vs->data.vecptr);
    IGRAPH_CHECK(igraph_vector_copy((igraph_vector_t*)vs->data.vecptr, v));
    IGRAPH_FINALLY_CLEAN(1);
    return 0;
}

/**
 * \function igraph_vs_seq
 * \brief Vertex set, an interval of vertices.
 *
 * Creates a vertex selector containing all vertices with vertex id
 * equal to or bigger than \c from and equal to or smaller than \c
 * to.
 *
 * \param vs Pointer to an uninitialized vertex selector object.
 * \param from The first vertex id to be included in the vertex
 *        selector.
 * \param to The last vertex id to be included in the vertex
 *        selector.
 * \return Error code.
 * \sa \ref igraph_vss_seq(), \ref igraph_vs_destroy()
 *
 * Time complexity: O(1).
 *
 * \example examples/simple/igraph_vs_seq.c
 */

int igraph_vs_seq(igraph_vs_t *vs,
                  igraph_integer_t from, igraph_integer_t to) {
    vs->type = IGRAPH_VS_SEQ;
    vs->data.seq.from = from;
    vs->data.seq.to = to + 1;
    return 0;
}

/**
 * \function igraph_vss_seq
 * \brief An interval of vertices (immediate version).
 *
 * The immediate version of \ref igraph_vs_seq().
 *
 * \param from The first vertex id to be included in the vertex
 *        selector.
 * \param to The last vertex id to be included in the vertex
 *        selector.
 * \return Error code.
 * \sa \ref igraph_vs_seq()
 *
 * Time complexity: O(1).
 */

igraph_vs_t igraph_vss_seq(igraph_integer_t from, igraph_integer_t to) {
    igraph_vs_t vs;
    vs.type = IGRAPH_VS_SEQ;
    vs.data.seq.from = from;
    vs.data.seq.to = to + 1;
    return vs;
}

/**
 * \function igraph_vs_destroy
 * \brief Destroy a vertex set.
 *
 * This function should be called for all vertex selectors when they
 * are not needed. The memory allocated for the vertex selector will
 * be deallocated. Do not call this function on vertex selectors
 * created with the immediate versions of the vertex selector
 * constructors (starting with <code>igraph_vss</code>).
 *
 * \param vs Pointer to a vertex selector object.
 *
 * Time complexity: operating system dependent, usually O(1).
 */

void igraph_vs_destroy(igraph_vs_t *vs) {
    switch (vs->type) {
    case IGRAPH_VS_ALL:
    case IGRAPH_VS_ADJ:
    case IGRAPH_VS_NONE:
    case IGRAPH_VS_1:
    case IGRAPH_VS_VECTORPTR:
    case IGRAPH_VS_SEQ:
    case IGRAPH_VS_NONADJ:
        break;
    case IGRAPH_VS_VECTOR:
        igraph_vector_destroy((igraph_vector_t*)vs->data.vecptr);
        IGRAPH_FREE(vs->data.vecptr);
        break;
    default:
        break;
    }
}

/**
 * \function igraph_vs_is_all
 * \brief Check whether all vertices are included.
 *
 * This function checks whether the vertex selector object was created
 * by \ref igraph_vs_all() or \ref igraph_vss_all(). Note that the
 * vertex selector might contain all vertices in a given graph but if
 * it wasn't created by the two constructors mentioned here the return
 * value will be FALSE.
 *
 * \param vs Pointer to a vertex selector object.
 * \return TRUE (1) if the vertex selector contains all vertices and
 *         FALSE (0) otherwise.
 *
 * Time complexity: O(1).
 */

igraph_bool_t igraph_vs_is_all(const igraph_vs_t *vs) {
    return vs->type == IGRAPH_VS_ALL;
}

int igraph_vs_as_vector(const igraph_t *graph, igraph_vs_t vs,
                        igraph_vector_t *v) {
    igraph_vit_t vit;

    IGRAPH_CHECK(igraph_vit_create(graph, vs, &vit));
    IGRAPH_FINALLY(igraph_vit_destroy, &vit);
    IGRAPH_CHECK(igraph_vit_as_vector(&vit, v));

    igraph_vit_destroy(&vit);
    IGRAPH_FINALLY_CLEAN(1);
    return 0;
}

/**
 * \function igraph_vs_copy
 * \brief Creates a copy of a vertex selector.
 * \param src The selector being copied.
 * \param dest An uninitialized selector that will contain the copy.
 */
int igraph_vs_copy(igraph_vs_t* dest, const igraph_vs_t* src) {
    memcpy(dest, src, sizeof(igraph_vs_t));
    switch (dest->type) {
    case IGRAPH_VS_VECTOR:
        dest->data.vecptr = IGRAPH_CALLOC(1, igraph_vector_t);
        if (!dest->data.vecptr) {
            IGRAPH_ERROR("Cannot copy vertex selector", IGRAPH_ENOMEM);
        }
        IGRAPH_CHECK(igraph_vector_copy((igraph_vector_t*)dest->data.vecptr,
                                        (igraph_vector_t*)src->data.vecptr));
        break;
    }
    return 0;
}

/**
 * \function igraph_vs_type
 * \brief Returns the type of the vertex selector.
 */
int igraph_vs_type(const igraph_vs_t *vs) {
    return vs->type;
}

/**
 * \function igraph_vs_size
 * \brief Returns the size of the vertex selector.
 *
 * The size of the vertex selector is the number of vertices it will
 * yield when it is iterated over.
 *
 * \param graph The graph over which we will iterate.
 * \param result The result will be returned here.
 */
int igraph_vs_size(const igraph_t *graph, const igraph_vs_t *vs,
                   igraph_integer_t *result) {
    igraph_vector_t vec;
    igraph_bool_t *seen;
    long i;

    switch (vs->type) {
    case IGRAPH_VS_NONE:
        *result = 0; return 0;

    case IGRAPH_VS_1:
        *result = 0;
        if (vs->data.vid < igraph_vcount(graph) && vs->data.vid >= 0) {
            *result = 1;
        }
        return 0;

    case IGRAPH_VS_SEQ:
        *result = vs->data.seq.to - vs->data.seq.from;
        return 0;

    case IGRAPH_VS_ALL:
        *result = igraph_vcount(graph); return 0;

    case IGRAPH_VS_ADJ:
        IGRAPH_VECTOR_INIT_FINALLY(&vec, 0);
        IGRAPH_CHECK(igraph_neighbors(graph, &vec, vs->data.adj.vid, vs->data.adj.mode));
        *result = (igraph_integer_t) igraph_vector_size(&vec);
        igraph_vector_destroy(&vec);
        IGRAPH_FINALLY_CLEAN(1);
        return 0;

    case IGRAPH_VS_NONADJ:
        IGRAPH_VECTOR_INIT_FINALLY(&vec, 0);
        IGRAPH_CHECK(igraph_neighbors(graph, &vec, vs->data.adj.vid, vs->data.adj.mode));
        *result = igraph_vcount(graph);
        seen = IGRAPH_CALLOC(*result, igraph_bool_t);
        if (seen == 0) {
            IGRAPH_ERROR("Cannot calculate selector length", IGRAPH_ENOMEM);
        }
        IGRAPH_FINALLY(igraph_free, seen);
        for (i = 0; i < igraph_vector_size(&vec); i++) {
            if (!seen[(long int)VECTOR(vec)[i]]) {
                (*result)--;
                seen[(long int)VECTOR(vec)[i]] = 1;
            }
        }
        igraph_free(seen);
        igraph_vector_destroy(&vec);
        IGRAPH_FINALLY_CLEAN(2);
        return 0;

    case IGRAPH_VS_VECTOR:
    case IGRAPH_VS_VECTORPTR:
        *result = (igraph_integer_t) igraph_vector_size((igraph_vector_t*)vs->data.vecptr);
        return 0;
    }

    IGRAPH_ERROR("Cannot calculate selector length, invalid selector type",
                 IGRAPH_EINVAL);
}

/***************************************************/

/**
 * \function igraph_vit_create
 * \brief Creates a vertex iterator from a vertex selector.
 *
 * This function instantiates a vertex selector object with a given
 * graph. This is the step when the actual vertex ids are created from
 * the \em logical notion of the vertex selector based on the graph.
 * E.g. a vertex selector created with \ref igraph_vs_all() contains
 * knowledge that \em all vertices are included in a (yet indefinite)
 * graph. When instantiating it a vertex iterator object is created,
 * this contains the actual vertex ids in the graph supplied as a
 * parameter.
 *
 * </para><para>
 * The same vertex selector object can be used to instantiate any
 * number vertex iterators.
 *
 * \param graph An \type igraph_t object, a graph.
 * \param vs A vertex selector object.
 * \param vit Pointer to an uninitialized vertex iterator object.
 * \return Error code.
 * \sa \ref igraph_vit_destroy().
 *
 * Time complexity: it depends on the vertex selector type. O(1) for
 * vertex selectors created with \ref igraph_vs_all(), \ref
 * igraph_vs_none(), \ref igraph_vs_1, \ref igraph_vs_vector, \ref
 * igraph_vs_seq(), \ref igraph_vs_vector(), \ref
 * igraph_vs_vector_small(). O(d) for \ref igraph_vs_adj(), d is the
 * number of vertex ids to be included in the iterator. O(|V|) for
 * \ref igraph_vs_nonadj(), |V| is the number of vertices in the graph.
 */

int igraph_vit_create(const igraph_t *graph,
                      igraph_vs_t vs, igraph_vit_t *vit) {
    igraph_vector_t vec;
    igraph_bool_t *seen;
    long int i, j, n;

    switch (vs.type) {
    case IGRAPH_VS_ALL:
        vit->type = IGRAPH_VIT_SEQ;
        vit->pos = 0;
        vit->start = 0;
        vit->end = igraph_vcount(graph);
        break;
    case IGRAPH_VS_ADJ:
        vit->type = IGRAPH_VIT_VECTOR;
        vit->pos = 0;
        vit->start = 0;
        vit->vec = IGRAPH_CALLOC(1, igraph_vector_t);
        if (vit->vec == 0) {
            IGRAPH_ERROR("Cannot create iterator", IGRAPH_ENOMEM);
        }
        IGRAPH_FINALLY(igraph_free, (igraph_vector_t*) vit->vec);
        IGRAPH_VECTOR_INIT_FINALLY((igraph_vector_t*)vit->vec, 0);
        IGRAPH_CHECK(igraph_neighbors(graph, (igraph_vector_t*)vit->vec,
                                      vs.data.adj.vid, vs.data.adj.mode));
        vit->end = igraph_vector_size(vit->vec);
        IGRAPH_FINALLY_CLEAN(2);
        break;
    case IGRAPH_VS_NONADJ:
        vit->type = IGRAPH_VIT_VECTOR;
        vit->pos = 0;
        vit->start = 0;
        vit->vec = IGRAPH_CALLOC(1, igraph_vector_t);
        if (vit->vec == 0) {
            IGRAPH_ERROR("Cannot create iterator", IGRAPH_ENOMEM);
        }
        IGRAPH_FINALLY(igraph_free, (igraph_vector_t*) vit->vec);
        IGRAPH_VECTOR_INIT_FINALLY((igraph_vector_t *) vit->vec, 0);
        IGRAPH_VECTOR_INIT_FINALLY(&vec, 0);
        IGRAPH_CHECK(igraph_neighbors(graph, &vec,
                                      vs.data.adj.vid, vs.data.adj.mode));
        n = igraph_vcount(graph);
        seen = IGRAPH_CALLOC(n, igraph_bool_t);
        if (seen == 0) {
            IGRAPH_ERROR("Cannot create iterator", IGRAPH_ENOMEM);
        }
        IGRAPH_FINALLY(igraph_free, seen);
        for (i = 0; i < igraph_vector_size(&vec); i++) {
            if (! seen [ (long int) VECTOR(vec)[i] ] ) {
                n--;
                seen[ (long int) VECTOR(vec)[i] ] = 1;
            }
        }
        IGRAPH_CHECK(igraph_vector_resize((igraph_vector_t*)vit->vec, n));
        for (i = 0, j = 0; j < n; i++) {
            if (!seen[i]) {
                VECTOR(*vit->vec)[j++] = i;
            }
        }

        IGRAPH_FREE(seen);
        igraph_vector_destroy(&vec);
        vit->end = n;
        IGRAPH_FINALLY_CLEAN(4);
        break;
    case IGRAPH_VS_NONE:
        vit->type = IGRAPH_VIT_SEQ;
        vit->pos = 0;
        vit->start = 0;
        vit->end = 0;
        break;
    case IGRAPH_VS_1:
        vit->type = IGRAPH_VIT_SEQ;
        vit->pos = vs.data.vid;
        vit->start = vs.data.vid;
        vit->end = vs.data.vid + 1;
        if (vit->pos >= igraph_vcount(graph)) {
            IGRAPH_ERROR("Cannot create iterator, invalid vertex id", IGRAPH_EINVVID);
        }
        break;
    case IGRAPH_VS_VECTORPTR:
    case IGRAPH_VS_VECTOR:
        vit->type = IGRAPH_VIT_VECTORPTR;
        vit->pos = 0;
        vit->start = 0;
        vit->vec = vs.data.vecptr;
        vit->end = igraph_vector_size(vit->vec);
        if (!igraph_vector_isininterval(vit->vec, 0, igraph_vcount(graph) - 1)) {
            IGRAPH_ERROR("Cannot create iterator, invalid vertex id", IGRAPH_EINVVID);
        }
        break;
    case IGRAPH_VS_SEQ:
        vit->type = IGRAPH_VIT_SEQ;
        vit->pos = vs.data.seq.from;
        vit->start = vs.data.seq.from;
        vit->end = vs.data.seq.to;
        break;
    default:
        IGRAPH_ERROR("Cannot create iterator, invalid selector", IGRAPH_EINVAL);
        break;
    }
    return 0;
}

/**
 * \function igraph_vit_destroy
 * \brief Destroys a vertex iterator.
 *
 * </para><para>
 * Deallocates memory allocated for a vertex iterator.
 *
 * \param vit Pointer to an initialized vertex iterator object.
 * \sa \ref igraph_vit_create()
 *
 * Time complexity: operating system dependent, usually O(1).
 */

void igraph_vit_destroy(const igraph_vit_t *vit) {
    switch (vit->type) {
    case IGRAPH_VIT_SEQ:
    case IGRAPH_VIT_VECTORPTR:
        break;
    case IGRAPH_VIT_VECTOR:
        igraph_vector_destroy((igraph_vector_t*)vit->vec);
        igraph_free((igraph_vector_t*)vit->vec);
        break;
    default:
        /*     IGRAPH_ERROR("Cannot destroy iterator, unknown type", IGRAPH_EINVAL); */
        break;
    }
}

int igraph_vit_as_vector(const igraph_vit_t *vit, igraph_vector_t *v) {

    long int i;

    IGRAPH_CHECK(igraph_vector_resize(v, IGRAPH_VIT_SIZE(*vit)));

    switch (vit->type) {
    case IGRAPH_VIT_SEQ:
        for (i = 0; i < IGRAPH_VIT_SIZE(*vit); i++) {
            VECTOR(*v)[i] = vit->start + i;
        }
        break;
    case IGRAPH_VIT_VECTOR:
    case IGRAPH_VIT_VECTORPTR:
        for (i = 0; i < IGRAPH_VIT_SIZE(*vit); i++) {
            VECTOR(*v)[i] = VECTOR(*vit->vec)[i];
        }
        break;
    default:
        IGRAPH_ERROR("Cannot convert to vector, unknown iterator type",
                     IGRAPH_EINVAL);
        break;
    }

    return 0;
}

/*******************************************************/

/**
 * \function igraph_es_all
 * \brief Edge set, all edges.
 *
 * \param es Pointer to an uninitialized edge selector object.
 * \param order Constant giving the order in which the edges will be
 *        included in the selector. Possible values:
 *        \c IGRAPH_EDGEORDER_ID, edge id order.
 *        \c IGRAPH_EDGEORDER_FROM, vertex id order, the id of the
 *           \em source vertex counts for directed graphs. The order
 *           of the incident edges of a given vertex is arbitrary.
 *        \c IGRAPH_EDGEORDER_TO, vertex id order, the id of the \em
 *           target vertex counts for directed graphs. The order
 *           of the incident edges of a given vertex is arbitrary.
 *        For undirected graph the latter two is the same.
 * \return Error code.
 * \sa \ref igraph_ess_all(), \ref igraph_es_destroy()
 *
 * Time complexity: O(1).
 */

int igraph_es_all(igraph_es_t *es,
                  igraph_edgeorder_type_t order) {
    switch (order) {
    case IGRAPH_EDGEORDER_ID:
        es->type = IGRAPH_ES_ALL;
        break;
    case IGRAPH_EDGEORDER_FROM:
        es->type = IGRAPH_ES_ALLFROM;
        break;
    case IGRAPH_EDGEORDER_TO:
        es->type = IGRAPH_ES_ALLTO;
        break;
    default:
        IGRAPH_ERROR("Invalid edge order, cannot create selector", IGRAPH_EINVAL);
        break;
    }
    return 0;
}

/**
 * \function igraph_ess_all
 * \brief Edge set, all edges (immediate version)
 *
 * The immediate version of the all-edges selector.
 *
 * \param order Constant giving the order of the edges in the edge
 *        selector. See \ref igraph_es_all() for the possible values.
 * \return The edge selector.
 * \sa \ref igraph_es_all()
 *
 * Time complexity: O(1).
 */

igraph_es_t igraph_ess_all(igraph_edgeorder_type_t order) {
    igraph_es_t es;
    igraph_es_all(&es, order); /* cannot fail */
    return es;
}

/**
 * \function igraph_es_incident
 * \brief Edges incident on a given vertex.
 *
 * \param es Pointer to an uninitialized edge selector object.
 * \param vid Vertex id, of which the incident edges will be
 *        selected.
 * \param mode Constant giving the type of the incident edges to
 *        select. This is ignored for undirected graphs. Possible values:
 *        \c IGRAPH_OUT, outgoing edges;
 *        \c IGRAPH_IN, incoming edges;
 *        \c IGRAPH_ALL, all edges.
 * \return Error code.
 * \sa \ref igraph_es_destroy()
 *
 * Time complexity: O(1).
 */

int igraph_es_incident(igraph_es_t *es,
                       igraph_integer_t vid, igraph_neimode_t mode) {
    es->type = IGRAPH_ES_INCIDENT;
    es->data.incident.vid = vid;
    es->data.incident.mode = mode;
    return 0;
}

/**
 * \function igraph_es_none
 * \brief Empty edge selector.
 *
 * \param es Pointer to an uninitialized edge selector object to
 * initialize.
 * \return Error code.
 * \sa \ref igraph_ess_none(), \ref igraph_es_destroy()
 *
 * Time complexity: O(1).
 */

int igraph_es_none(igraph_es_t *es) {
    es->type = IGRAPH_ES_NONE;
    return 0;
}

/**
 * \function igraph_ess_none
 * \brief Immediate empty edge selector.
 *
 * </para><para>
 * Immediate version of the empty edge selector.
 *
 * \return Initialized empty edge selector.
 * \sa \ref igraph_es_none()
 *
 * Time complexity: O(1).
 */

igraph_es_t igraph_ess_none(void) {
    igraph_es_t es;
    es.type = IGRAPH_ES_NONE;
    return es;
}

/**
 * \function igraph_es_1
 * \brief Edge selector containing a single edge.
 *
 * \param es Pointer to an uninitialized edge selector object.
 * \param eid Edge id of the edge to select.
 * \return Error code.
 * \sa \ref igraph_ess_1(), \ref igraph_es_destroy()
 *
 * Time complexity: O(1).
 */

int igraph_es_1(igraph_es_t *es, igraph_integer_t eid) {
    es->type = IGRAPH_ES_1;
    es->data.eid = eid;
    return 0;
}

/**
 * \function igraph_ess_1
 * \brief Immediate version of the single edge edge selector.
 *
 * \param eid The id of the edge.
 * \return The edge selector.
 * \sa \ref igraph_es_1()
 *
 * Time complexity: O(1).
 */

igraph_es_t igraph_ess_1(igraph_integer_t eid) {
    igraph_es_t es;
    es.type = IGRAPH_ES_1;
    es.data.eid = eid;
    return es;
}

/**
 * \function igraph_es_vector
 * \brief Handle a vector as an edge selector.
 *
 * </para><para>
 * Creates an edge selector which serves as a view to a vector
 * containing edge ids. Do not destroy the vector before destroying
 * the view.
 *
 * Many views can be created to the same vector.
 *
 * \param es Pointer to an uninitialized edge selector.
 * \param v Vector containing edge ids.
 * \return Error code.
 * \sa \ref igraph_ess_vector(), \ref igraph_es_destroy()
 *
 * Time complexity: O(1).
 */

int igraph_es_vector(igraph_es_t *es,
                     const igraph_vector_t *v) {
    es->type = IGRAPH_ES_VECTORPTR;
    es->data.vecptr = v;
    return 0;
}

/**
 * \function igraph_es_vector_copy
 * \brief Edge set, based on a vector, with copying.
 *
 *
 * This function makes it possible to handle a \type vector_t
 * permanently as an edge selector. The edge selector creates a
 * copy of the original vector, so the vector can safely be destroyed
 * after creating the edge selector. Changing the original vector
 * will not affect the edge selector. The edge selector is
 * responsible for deleting the copy made by itself.
 *
 * \param es Pointer to an uninitialized edge selector.
 * \param v Pointer to a \type igraph_vector_t object.
 * \return Error code.
 * \sa \ref igraph_es_destroy()
 *
 * Time complexity: O(1).
 */

int igraph_es_vector_copy(igraph_es_t *es, const igraph_vector_t *v) {
    es->type = IGRAPH_ES_VECTOR;
    es->data.vecptr = IGRAPH_CALLOC(1, igraph_vector_t);
    if (es->data.vecptr == 0) {
        IGRAPH_ERROR("Cannot create edge selector", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, (igraph_vector_t*)es->data.vecptr);
    IGRAPH_CHECK(igraph_vector_copy((igraph_vector_t*)es->data.vecptr, v));
    IGRAPH_FINALLY_CLEAN(1);
    return 0;
}

/**
 * \function igraph_ess_vector
 * \brief Immediate vector view edge selector.
 *
 * </para><para>
 * This is the immediate version of the vector of edge ids edge
 * selector.
 *
 * \param v The vector of edge ids.
 * \return Edge selector, initialized.
 * \sa \ref igraph_es_vector()
 *
 * Time complexity: O(1).
 */

igraph_es_t igraph_ess_vector(const igraph_vector_t *v) {
    igraph_es_t es;
    es.type = IGRAPH_ES_VECTORPTR;
    es.data.vecptr = v;
    return es;
}

/**
 * \function igraph_es_fromto
 * \brief Edge selector, all edges between two vertex sets.
 *
 * </para><para>
 * This function is not implemented yet.
 *
 * \param es Pointer to an uninitialized edge selector.
 * \param from Vertex selector, their outgoing edges will be
 *        selected.
 * \param to Vertex selector, their incoming edges will be selected
 *        from the previous selection.
 * \return Error code.
 * \sa \ref igraph_es_destroy()
 *
 * Time complexity: O(1).
 *
 * \example examples/simple/igraph_es_fromto.c
 */

int igraph_es_fromto(igraph_es_t *es,
                     igraph_vs_t from, igraph_vs_t to) {

    IGRAPH_UNUSED(es); IGRAPH_UNUSED(from); IGRAPH_UNUSED(to);
    IGRAPH_ERROR("igraph_es_fromto not implemented yet", IGRAPH_UNIMPLEMENTED);
    /* TODO */
}

/**
 * \function igraph_es_seq
 * \brief Edge selector, a sequence of edge ids.
 *
 * All edge ids between <code>from</code> and <code>to</code> will be
 * included in the edge selection. This includes <code>from</code> and
 * excludes <code>to</code>.
 *
 * \param es Pointer to an uninitialized edge selector object.
 * \param from The first edge id to be included.
 * \param to The last edge id to be included.
 * \return Error code.
 * \sa \ref igraph_ess_seq(), \ref igraph_es_destroy()
 *
 * Time complexity: O(1).
 */

int igraph_es_seq(igraph_es_t *es,
                  igraph_integer_t from, igraph_integer_t to) {
    es->type = IGRAPH_ES_SEQ;
    es->data.seq.from = from;
    es->data.seq.to = to;
    return 0;
}

/**
 * \function igraph_ess_seq
 * \brief Immediate version of the sequence edge selector.
 *
 * \param from The first edge id to include.
 * \param to The last edge id to include.
 * \return The initialized edge selector.
 * \sa \ref igraph_es_seq()
 *
 * Time complexity: O(1).
 */

igraph_es_t igraph_ess_seq(igraph_integer_t from, igraph_integer_t to) {
    igraph_es_t es;
    es.type = IGRAPH_ES_SEQ;
    es.data.seq.from = from;
    es.data.seq.to = to;
    return es;
}

/**
 * \function igraph_es_pairs
 * \brief Edge selector, multiple edges defined by their endpoints in a vector.
 *
 * The edges between the given pairs of vertices will be included in the
 * edge selection. The vertex pairs must be defined in the vector <code>v</code>,
 * the first element of the vector is the first vertex of the first edge
 * to be selected, the second element is the second vertex of the first
 * edge, the third element is the first vertex of the second edge and
 * so on.
 *
 * \param es Pointer to an uninitialized edge selector object.
 * \param v The vector containing the endpoints of the edges.
 * \param directed Whether the graph is directed or not.
 * \return Error code.
 * \sa \ref igraph_es_pairs_small(), \ref igraph_es_destroy()
 *
 * Time complexity: O(n), the number of edges being selected.
 *
 * \example examples/simple/igraph_es_pairs.c
 */

int igraph_es_pairs(igraph_es_t *es, const igraph_vector_t *v,
                    igraph_bool_t directed) {
    es->type = IGRAPH_ES_PAIRS;
    es->data.path.mode = directed;
    es->data.path.ptr = IGRAPH_CALLOC(1, igraph_vector_t);
    if (es->data.path.ptr == 0) {
        IGRAPH_ERROR("Cannot create edge selector", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, (igraph_vector_t*) es->data.path.ptr);

    IGRAPH_CHECK(igraph_vector_copy((igraph_vector_t*) es->data.path.ptr, v));

    IGRAPH_FINALLY_CLEAN(1);
    return 0;
}

/**
 * \function igraph_es_pairs_small
 * \brief Edge selector, multiple edges defined by their endpoints as arguments.
 *
 * The edges between the given pairs of vertices will be included in the
 * edge selection. The vertex pairs must be given as the arguments of the
 * function call, the third argument is the first vertex of the first edge,
 * the fourth argument is the second vertex of the first edge, the fifth
 * is the first vertex of the second edge and so on. The last element of the
 * argument list must be -1 to denote the end of the argument list.
 *
 * \param es Pointer to an uninitialized edge selector object.
 * \param directed Whether the graph is directed or not.
 * \return Error code.
 * \sa \ref igraph_es_pairs(), \ref igraph_es_destroy()
 *
 * Time complexity: O(n), the number of edges being selected.
 */

int igraph_es_pairs_small(igraph_es_t *es, igraph_bool_t directed, ...) {
    va_list ap;
    long int i, n = 0;
    es->type = IGRAPH_ES_PAIRS;
    es->data.path.mode = directed;
    es->data.path.ptr = IGRAPH_CALLOC(1, igraph_vector_t);
    if (es->data.path.ptr == 0) {
        IGRAPH_ERROR("Cannot create edge selector", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, (igraph_vector_t*)es->data.path.ptr);

    va_start(ap, directed);
    while (1) {
        int num = va_arg(ap, int);
        if (num == -1) {
            break;
        }
        n++;
    }
    va_end(ap);

    IGRAPH_VECTOR_INIT_FINALLY( (igraph_vector_t*) es->data.path.ptr, n);

    va_start(ap, directed);
    for (i = 0; i < n; i++) {
        VECTOR(*es->data.path.ptr)[i] = (igraph_real_t) va_arg(ap, int);
    }
    va_end(ap);

    IGRAPH_FINALLY_CLEAN(2);
    return 0;
}

int igraph_es_multipairs(igraph_es_t *es, const igraph_vector_t *v,
                         igraph_bool_t directed) {
    es->type = IGRAPH_ES_MULTIPAIRS;
    es->data.path.mode = directed;
    es->data.path.ptr = IGRAPH_CALLOC(1, igraph_vector_t);
    if (es->data.path.ptr == 0) {
        IGRAPH_ERROR("Cannot create edge selector", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, (igraph_vector_t*) es->data.path.ptr);

    IGRAPH_CHECK(igraph_vector_copy((igraph_vector_t*) es->data.path.ptr, v));

    IGRAPH_FINALLY_CLEAN(1);
    return 0;
}

/**
 * \function igraph_es_path
 * \brief Edge selector, edge ids on a path.
 *
 * This function takes a vector of vertices and creates a selector of
 * edges between those vertices. Vector {0, 3, 4, 7} will select edges
 * (0 -> 3), (3 -> 4), (4 -> 7). If these edges don't exist then trying
 * to create an iterator using this selector will fail.
 *
 * \param es Pointer to an uninitialized edge selector object.
 * \param v Pointer to a vector of vertex id's along the path.
 * \param directed If edge directions should be taken into account. This
 *                 will be ignored if the graph to select from is undirected.
 * \return Error code.
 * \sa \ref igraph_es_destroy()
 *
 * Time complexity: O(n), the number of vertices.
 */
int igraph_es_path(igraph_es_t *es, const igraph_vector_t *v,
                   igraph_bool_t directed) {
    es->type = IGRAPH_ES_PATH;
    es->data.path.mode = directed;
    es->data.path.ptr = IGRAPH_CALLOC(1, igraph_vector_t);
    if (es->data.path.ptr == 0) {
        IGRAPH_ERROR("Cannot create edge selector", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, (igraph_vector_t*) es->data.path.ptr);

    IGRAPH_CHECK(igraph_vector_copy((igraph_vector_t*) es->data.path.ptr, v));

    IGRAPH_FINALLY_CLEAN(1);
    return 0;
}

int igraph_es_path_small(igraph_es_t *es, igraph_bool_t directed, ...) {
    va_list ap;
    long int i, n = 0;
    es->type = IGRAPH_ES_PATH;
    es->data.path.mode = directed;
    es->data.path.ptr = IGRAPH_CALLOC(1, igraph_vector_t);
    if (es->data.path.ptr == 0) {
        IGRAPH_ERROR("Cannot create edge selector", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, (igraph_vector_t*)es->data.path.ptr);

    va_start(ap, directed);
    while (1) {
        int num = va_arg(ap, int);
        if (num == -1) {
            break;
        }
        n++;
    }
    va_end(ap);

    IGRAPH_VECTOR_INIT_FINALLY( (igraph_vector_t*) es->data.path.ptr, n);

    va_start(ap, directed);
    for (i = 0; i < n; i++) {
        VECTOR(*es->data.path.ptr)[i] = (igraph_real_t) va_arg(ap, int);
    }
    va_end(ap);

    IGRAPH_FINALLY_CLEAN(2);
    return 0;
}

/**
 * \function igraph_es_destroy
 * \brief Destroys an edge selector object.
 *
 * </para><para>
 * Call this function on an edge selector when it is not needed any
 * more. Do \em not call this function on edge selectors created by
 * immediate constructors, those don't need to be destroyed.
 *
 * \param es Pointer to an edge selector object.
 *
 * Time complexity: operating system dependent, usually O(1).
 */

void igraph_es_destroy(igraph_es_t *es) {
    switch (es->type) {
    case IGRAPH_ES_ALL:
    case IGRAPH_ES_ALLFROM:
    case IGRAPH_ES_ALLTO:
    case IGRAPH_ES_INCIDENT:
    case IGRAPH_ES_NONE:
    case IGRAPH_ES_1:
    case IGRAPH_ES_VECTORPTR:
    case IGRAPH_ES_SEQ:
        break;
    case IGRAPH_ES_VECTOR:
        igraph_vector_destroy((igraph_vector_t*)es->data.vecptr);
        IGRAPH_FREE(es->data.vecptr);
        break;
    case IGRAPH_ES_PAIRS:
    case IGRAPH_ES_PATH:
    case IGRAPH_ES_MULTIPAIRS:
        igraph_vector_destroy((igraph_vector_t*)es->data.path.ptr);
        IGRAPH_FREE(es->data.path.ptr);
        break;
    default:
        break;
    }
}

/**
 * \function igraph_es_is_all
 * \brief Check whether an edge selector includes all edges.
 *
 * \param es Pointer to an edge selector object.
 * \return TRUE (1) if <code>es</code> was created with \ref
 * igraph_es_all() or \ref igraph_ess_all(), and FALSE (0) otherwise.
 *
 * Time complexity: O(1).
 */

igraph_bool_t igraph_es_is_all(const igraph_es_t *es) {
    return es->type == IGRAPH_ES_ALL;
}

/**
 * \function igraph_es_copy
 * \brief Creates a copy of an edge selector.
 * \param src The selector being copied.
 * \param dest An uninitialized selector that will contain the copy.
 * \sa \ref igraph_es_destroy()
 */
int igraph_es_copy(igraph_es_t* dest, const igraph_es_t* src) {
    memcpy(dest, src, sizeof(igraph_es_t));
    switch (dest->type) {
    case IGRAPH_ES_VECTOR:
        dest->data.vecptr = IGRAPH_CALLOC(1, igraph_vector_t);
        if (!dest->data.vecptr) {
            IGRAPH_ERROR("Cannot copy edge selector", IGRAPH_ENOMEM);
        }
        IGRAPH_CHECK(igraph_vector_copy((igraph_vector_t*)dest->data.vecptr,
                                        (igraph_vector_t*)src->data.vecptr));
        break;
    case IGRAPH_ES_PATH:
    case IGRAPH_ES_PAIRS:
    case IGRAPH_ES_MULTIPAIRS:
        dest->data.path.ptr = IGRAPH_CALLOC(1, igraph_vector_t);
        if (!dest->data.path.ptr) {
            IGRAPH_ERROR("Cannot copy edge selector", IGRAPH_ENOMEM);
        }
        IGRAPH_CHECK(igraph_vector_copy((igraph_vector_t*)dest->data.path.ptr,
                                        (igraph_vector_t*)src->data.path.ptr));
        break;
    }
    return 0;
}

/**
 * \function igraph_es_as_vector
 * \brief Transform edge selector into vector.
 *
 * </para><para>
 * Call this function on an edge selector to transform it into a vector.
 * This is only implemented for sequence and vector selectors. If the
 * edges do not exist in the graph, this will result in an error.
 *
 * \param graph Pointer to a graph to check if the edges in the selector exist.
 * \param es An edge selector object.
 * \param v Pointer to initialized vector. The result will be stored here.
 *
 * Time complexity: O(n), the number of edges in the selector.
 */
int igraph_es_as_vector(const igraph_t *graph, igraph_es_t es,
                        igraph_vector_t *v) {
    igraph_eit_t eit;

    IGRAPH_CHECK(igraph_eit_create(graph, es, &eit));
    IGRAPH_FINALLY(igraph_eit_destroy, &eit);
    IGRAPH_CHECK(igraph_eit_as_vector(&eit, v));

    igraph_eit_destroy(&eit);
    IGRAPH_FINALLY_CLEAN(1);
    return 0;
}

/**
 * \function igraph_es_type
 * \brief Returns the type of the edge selector.
 */
int igraph_es_type(const igraph_es_t *es) {
    return es->type;
}

static int igraph_i_es_pairs_size(const igraph_t *graph,
                                  const igraph_es_t *es, igraph_integer_t *result);
static int igraph_i_es_path_size(const igraph_t *graph,
                                 const igraph_es_t *es, igraph_integer_t *result);
static int igraph_i_es_multipairs_size(const igraph_t *graph,
                                       const igraph_es_t *es, igraph_integer_t *result);

/**
 * \function igraph_es_size
 * \brief Returns the size of the edge selector.
 *
 * The size of the edge selector is the number of edges it will
 * yield when it is iterated over.
 *
 * \param graph The graph over which we will iterate.
 * \param result The result will be returned here.
 */
int igraph_es_size(const igraph_t *graph, const igraph_es_t *es,
                   igraph_integer_t *result) {
    igraph_vector_t v;

    switch (es->type) {
    case IGRAPH_ES_ALL:
        *result = igraph_ecount(graph);
        return 0;

    case IGRAPH_ES_ALLFROM:
        *result = igraph_ecount(graph);
        return 0;

    case IGRAPH_ES_ALLTO:
        *result = igraph_ecount(graph);
        return 0;

    case IGRAPH_ES_INCIDENT:
        IGRAPH_VECTOR_INIT_FINALLY(&v, 0);
        IGRAPH_CHECK(igraph_incident(graph, &v,
                                     es->data.incident.vid, es->data.incident.mode));
        *result = (igraph_integer_t) igraph_vector_size(&v);
        igraph_vector_destroy(&v);
        IGRAPH_FINALLY_CLEAN(1);
        return 0;

    case IGRAPH_ES_NONE:
        *result = 0;
        return 0;

    case IGRAPH_ES_1:
        if (es->data.eid < igraph_ecount(graph) && es->data.eid >= 0) {
            *result = 1;
        } else {
            *result = 0;
        }
        return 0;

    case IGRAPH_ES_VECTOR:
    case IGRAPH_ES_VECTORPTR:
        *result = (igraph_integer_t) igraph_vector_size((igraph_vector_t*)es->data.vecptr);
        return 0;

    case IGRAPH_ES_SEQ:
        *result = es->data.seq.to - es->data.seq.from;
        return 0;

    case IGRAPH_ES_PAIRS:
        IGRAPH_CHECK(igraph_i_es_pairs_size(graph, es, result));
        return 0;

    case IGRAPH_ES_PATH:
        IGRAPH_CHECK(igraph_i_es_path_size(graph, es, result));
        return 0;

    case IGRAPH_ES_MULTIPAIRS:
        IGRAPH_CHECK(igraph_i_es_multipairs_size(graph, es, result));
        return 0;

    default:
        IGRAPH_ERROR("Cannot calculate selector length, invalid selector type",
                     IGRAPH_EINVAL);
    }
}

static int igraph_i_es_pairs_size(const igraph_t *graph,
                                  const igraph_es_t *es, igraph_integer_t *result) {
    long int n = igraph_vector_size(es->data.path.ptr);
    long int no_of_nodes = igraph_vcount(graph);
    long int i;

    if (n % 2 != 0) {
        IGRAPH_ERROR("Cannot calculate edge selector length from odd number of vertices",
                     IGRAPH_EINVAL);
    }
    if (!igraph_vector_isininterval(es->data.path.ptr, 0, no_of_nodes - 1)) {
        IGRAPH_ERROR("Cannot calculate edge selector length", IGRAPH_EINVVID);
    }

    *result = (igraph_integer_t) (n / 2);
    /* Check for the existence of all edges */
    for (i = 0; i < *result; i++) {
        long int from = (long int) VECTOR(*es->data.path.ptr)[2 * i];
        long int to = (long int) VECTOR(*es->data.path.ptr)[2 * i + 1];
        igraph_integer_t eid;
        IGRAPH_CHECK(igraph_get_eid(graph, &eid, (igraph_integer_t) from,
                                    (igraph_integer_t) to, es->data.path.mode,
                                    /*error=*/ 1));
    }

    return 0;
}

static int igraph_i_es_path_size(const igraph_t *graph,
                                 const igraph_es_t *es, igraph_integer_t *result) {
    long int n = igraph_vector_size(es->data.path.ptr);
    long int no_of_nodes = igraph_vcount(graph);
    long int i;

    if (!igraph_vector_isininterval(es->data.path.ptr, 0, no_of_nodes - 1)) {
        IGRAPH_ERROR("Cannot calculate selector length", IGRAPH_EINVVID);
    }

    if (n <= 1) {
        *result = 0;
    } else {
        *result = (igraph_integer_t) (n - 1);
    }
    for (i = 0; i < *result; i++) {
        long int from = (long int) VECTOR(*es->data.path.ptr)[i];
        long int to = (long int) VECTOR(*es->data.path.ptr)[i + 1];
        igraph_integer_t eid;
        IGRAPH_CHECK(igraph_get_eid(graph, &eid, (igraph_integer_t) from,
                                    (igraph_integer_t) to, es->data.path.mode,
                                    /*error=*/ 1));
    }

    return 0;
}

static int igraph_i_es_multipairs_size(const igraph_t *graph,
                                       const igraph_es_t *es, igraph_integer_t *result) {
    IGRAPH_UNUSED(graph); IGRAPH_UNUSED(es); IGRAPH_UNUSED(result);
    IGRAPH_ERROR("Cannot calculate edge selector length", IGRAPH_UNIMPLEMENTED);
}

/**************************************************/

static int igraph_i_eit_create_allfromto(const igraph_t *graph,
                                         igraph_eit_t *eit,
                                         igraph_neimode_t mode);
static int igraph_i_eit_pairs(const igraph_t *graph,
                              igraph_es_t es, igraph_eit_t *eit);
static int igraph_i_eit_multipairs(const igraph_t *graph,
                                   igraph_es_t es, igraph_eit_t *eit);
static int igraph_i_eit_path(const igraph_t *graph,
                             igraph_es_t es, igraph_eit_t *eit);

static int igraph_i_eit_create_allfromto(const igraph_t *graph,
                                         igraph_eit_t *eit,
                                         igraph_neimode_t mode) {
    igraph_vector_t *vec;
    long int no_of_nodes = igraph_vcount(graph);
    long int no_of_edges = igraph_ecount(graph);
    long int i;

    vec = IGRAPH_CALLOC(1, igraph_vector_t);
    if (vec == 0) {
        IGRAPH_ERROR("Cannot create edge iterator", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, vec);
    IGRAPH_VECTOR_INIT_FINALLY(vec, 0);
    IGRAPH_CHECK(igraph_vector_reserve(vec, no_of_edges));

    if (igraph_is_directed(graph)) {
        igraph_vector_t adj;
        IGRAPH_VECTOR_INIT_FINALLY(&adj, 0);
        for (i = 0; i < no_of_nodes; i++) {
            igraph_incident(graph, &adj, (igraph_integer_t) i, mode);
            igraph_vector_append(vec, &adj);
        }
        igraph_vector_destroy(&adj);
        IGRAPH_FINALLY_CLEAN(1);

    } else {

        igraph_vector_t adj;
        igraph_bool_t *added;
        long int j;
        IGRAPH_VECTOR_INIT_FINALLY(&adj, 0);
        added = IGRAPH_CALLOC(no_of_edges, igraph_bool_t);
        if (added == 0) {
            IGRAPH_ERROR("Cannot create edge iterator", IGRAPH_ENOMEM);
        }
        IGRAPH_FINALLY(igraph_free, added);
        for (i = 0; i < no_of_nodes; i++) {
            igraph_incident(graph, &adj, (igraph_integer_t) i, IGRAPH_ALL);
            for (j = 0; j < igraph_vector_size(&adj); j++) {
                if (!added[ (long int)VECTOR(adj)[j] ]) {
                    igraph_vector_push_back(vec, VECTOR(adj)[j]);
                    added[ (long int)VECTOR(adj)[j] ] += 1;
                }
            }
        }
        igraph_vector_destroy(&adj);
        IGRAPH_FREE(added);
        IGRAPH_FINALLY_CLEAN(2);
    }

    eit->type = IGRAPH_EIT_VECTOR;
    eit->pos = 0;
    eit->start = 0;
    eit->vec = vec;
    eit->end = igraph_vector_size(eit->vec);

    IGRAPH_FINALLY_CLEAN(2);
    return 0;
}

static int igraph_i_eit_pairs(const igraph_t *graph,
                              igraph_es_t es, igraph_eit_t *eit) {
    long int n = igraph_vector_size(es.data.path.ptr);
    long int no_of_nodes = igraph_vcount(graph);
    long int i;

    if (n % 2 != 0) {
        IGRAPH_ERROR("Cannot create edge iterator from odd number of vertices",
                     IGRAPH_EINVAL);
    }
    if (!igraph_vector_isininterval(es.data.path.ptr, 0, no_of_nodes - 1)) {
        IGRAPH_ERROR("Cannot create edge iterator", IGRAPH_EINVVID);
    }

    eit->type = IGRAPH_EIT_VECTOR;
    eit->pos = 0;
    eit->start = 0;
    eit->end = n / 2;
    eit->vec = IGRAPH_CALLOC(1, igraph_vector_t);
    if (eit->vec == 0) {
        IGRAPH_ERROR("Cannot create edge iterator", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, (igraph_vector_t*)eit->vec);
    IGRAPH_VECTOR_INIT_FINALLY((igraph_vector_t*)eit->vec, n / 2);

    for (i = 0; i < igraph_vector_size(eit->vec); i++) {
        long int from = (long int) VECTOR(*es.data.path.ptr)[2 * i];
        long int to = (long int) VECTOR(*es.data.path.ptr)[2 * i + 1];
        igraph_integer_t eid;
        IGRAPH_CHECK(igraph_get_eid(graph, &eid, (igraph_integer_t) from,
                                    (igraph_integer_t) to, es.data.path.mode,
                                    /*error=*/ 1));
        VECTOR(*eit->vec)[i] = eid;
    }

    IGRAPH_FINALLY_CLEAN(2);
    return 0;
}

static int igraph_i_eit_multipairs(const igraph_t *graph,
                                   igraph_es_t es, igraph_eit_t *eit) {
    long int n = igraph_vector_size(es.data.path.ptr);
    long int no_of_nodes = igraph_vcount(graph);

    if (n % 2 != 0) {
        IGRAPH_ERROR("Cannot create edge iterator from odd number of vertices",
                     IGRAPH_EINVAL);
    }
    if (!igraph_vector_isininterval(es.data.path.ptr, 0, no_of_nodes - 1)) {
        IGRAPH_ERROR("Cannot create edge iterator", IGRAPH_EINVVID);
    }

    eit->type = IGRAPH_EIT_VECTOR;
    eit->pos = 0;
    eit->start = 0;
    eit->end = n / 2;
    eit->vec = IGRAPH_CALLOC(1, igraph_vector_t);
    if (eit->vec == 0) {
        IGRAPH_ERROR("Cannot create edge iterator", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, (igraph_vector_t*)eit->vec);
    IGRAPH_VECTOR_INIT_FINALLY((igraph_vector_t*)eit->vec, n / 2);

    IGRAPH_CHECK(igraph_get_eids_multi(graph, (igraph_vector_t *) eit->vec,
                                       /*pairs=*/ es.data.path.ptr, /*path=*/ 0,
                                       es.data.path.mode, /*error=*/ 1));

    IGRAPH_FINALLY_CLEAN(2);
    return 0;
}

static int igraph_i_eit_path(const igraph_t *graph,
                             igraph_es_t es, igraph_eit_t *eit) {
    long int n = igraph_vector_size(es.data.path.ptr);
    long int no_of_nodes = igraph_vcount(graph);
    long int i, len;

    if (!igraph_vector_isininterval(es.data.path.ptr, 0, no_of_nodes - 1)) {
        IGRAPH_ERROR("Cannot create edge iterator.", IGRAPH_EINVVID);
    }

    if (n <= 1) {
        len = 0;
    } else {
        len = n - 1;
    }

    eit->type = IGRAPH_EIT_VECTOR;
    eit->pos = 0;
    eit->start = 0;
    eit->end = len;
    eit->vec = IGRAPH_CALLOC(1, igraph_vector_t);
    if (eit->vec == 0) {
        IGRAPH_ERROR("Cannot create edge iterator.", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, (igraph_vector_t*)eit->vec);

    IGRAPH_VECTOR_INIT_FINALLY((igraph_vector_t *)eit->vec, len);

    for (i = 0; i < len; i++) {
        long int from = (long int) VECTOR(*es.data.path.ptr)[i];
        long int to = (long int) VECTOR(*es.data.path.ptr)[i + 1];
        igraph_integer_t eid;
        IGRAPH_CHECK(igraph_get_eid(graph, &eid, (igraph_integer_t) from,
                                    (igraph_integer_t) to, es.data.path.mode,
                                    /*error=*/ 1));
        VECTOR(*eit->vec)[i] = eid;
    }

    IGRAPH_FINALLY_CLEAN(2);
    return 0;
}

/**
 * \function igraph_eit_create
 * \brief Creates an edge iterator from an edge selector.
 *
 * </para><para>
 * This function creates an edge iterator based on an edge selector
 * and a graph.
 *
 * </para><para>
 * The same edge selector can be used to create many edge iterators,
 * also for different graphs.
 *
 * \param graph An \type igraph_t object for which the edge selector
 *        will be instantiated.
 * \param es The edge selector to instantiate.
 * \param eit Pointer to an uninitialized edge iterator.
 * \return Error code.
 * \sa \ref igraph_eit_destroy()
 *
 * Time complexity: depends on the type of the edge selector. For edge
 * selectors created by \ref igraph_es_all(), \ref igraph_es_none(),
 * \ref igraph_es_1(), igraph_es_vector(), igraph_es_seq() it is
 * O(1). For \ref igraph_es_incident() it is O(d) where d is the number of
 * incident edges of the vertex.
 */

int igraph_eit_create(const igraph_t *graph,
                      igraph_es_t es, igraph_eit_t *eit) {
    switch (es.type) {
    case IGRAPH_ES_ALL:
        eit->type = IGRAPH_EIT_SEQ;
        eit->pos = 0;
        eit->start = 0;
        eit->end = igraph_ecount(graph);
        break;
    case IGRAPH_ES_ALLFROM:
        IGRAPH_CHECK(igraph_i_eit_create_allfromto(graph, eit, IGRAPH_OUT));
        break;
    case IGRAPH_ES_ALLTO:
        IGRAPH_CHECK(igraph_i_eit_create_allfromto(graph, eit, IGRAPH_IN));
        break;
    case IGRAPH_ES_INCIDENT:
        eit->type = IGRAPH_EIT_VECTOR;
        eit->pos = 0;
        eit->start = 0;
        eit->vec = IGRAPH_CALLOC(1, igraph_vector_t);
        if (eit->vec == 0) {
            IGRAPH_ERROR("Cannot create iterator.", IGRAPH_ENOMEM);
        }
        IGRAPH_FINALLY(igraph_free, (igraph_vector_t*) eit->vec);
        IGRAPH_VECTOR_INIT_FINALLY((igraph_vector_t*)eit->vec, 0);
        IGRAPH_CHECK(igraph_incident(graph, (igraph_vector_t*)eit->vec,
                                     es.data.incident.vid, es.data.incident.mode));
        eit->end = igraph_vector_size(eit->vec);
        IGRAPH_FINALLY_CLEAN(2);
        break;
    case IGRAPH_ES_NONE:
        eit->type = IGRAPH_EIT_SEQ;
        eit->pos = 0;
        eit->start = 0;
        eit->end = 0;
        break;
    case IGRAPH_ES_1:
        eit->type = IGRAPH_EIT_SEQ;
        eit->pos = es.data.eid;
        eit->start = es.data.eid;
        eit->end = es.data.eid + 1;
        if (eit->pos >= igraph_ecount(graph)) {
            IGRAPH_ERROR("Cannot create iterator, invalid edge id.", IGRAPH_EINVAL);
        }
        break;
    case IGRAPH_ES_VECTOR:
    case IGRAPH_ES_VECTORPTR:
        eit->type = IGRAPH_EIT_VECTORPTR;
        eit->pos = 0;
        eit->start = 0;
        eit->vec = es.data.vecptr;
        eit->end = igraph_vector_size(eit->vec);
        if (!igraph_vector_isininterval(eit->vec, 0, igraph_ecount(graph) - 1)) {
            IGRAPH_ERROR("Cannot create iterator, invalid edge id.", IGRAPH_EINVAL);
        }
        break;
    case IGRAPH_ES_SEQ:
        eit->type = IGRAPH_EIT_SEQ;
        eit->pos = es.data.seq.from;
        eit->start = es.data.seq.from;
        eit->end = es.data.seq.to;
        if (eit->start < 0) {
            IGRAPH_ERROR("Cannot create iterator, invalid edge id.", IGRAPH_EINVAL);
        }
        if (eit->end < 0) {
            IGRAPH_ERROR("Cannot create iterator, invalid edge id.", IGRAPH_EINVAL);
        }
        if (eit->start >= igraph_ecount(graph)) {
            IGRAPH_ERROR("Cannot create iterator, starting edge greater than number of edges.", IGRAPH_EINVAL);
        }
        break;
    case IGRAPH_ES_PAIRS:
        IGRAPH_CHECK(igraph_i_eit_pairs(graph, es, eit));
        break;
    case IGRAPH_ES_MULTIPAIRS:
        IGRAPH_CHECK(igraph_i_eit_multipairs(graph, es, eit));
        break;
    case IGRAPH_ES_PATH:
        IGRAPH_CHECK(igraph_i_eit_path(graph, es, eit));
        break;
    default:
        IGRAPH_ERROR("Cannot create iterator, invalid selector.", IGRAPH_EINVAL);
        break;
    }
    return 0;
}

/**
 * \function igraph_eit_destroy
 * \brief Destroys an edge iterator.
 *
 * \param eit Pointer to an edge iterator to destroy.
 * \sa \ref igraph_eit_create()
 *
 * Time complexity: operating system dependent, usually O(1).
 */

void igraph_eit_destroy(const igraph_eit_t *eit) {
    switch (eit->type) {
    case IGRAPH_EIT_SEQ:
    case IGRAPH_EIT_VECTORPTR:
        break;
    case IGRAPH_EIT_VECTOR:
        igraph_vector_destroy((igraph_vector_t*)eit->vec);
        igraph_free((igraph_vector_t*)eit->vec);
        break;
    default:
        /*     IGRAPH_ERROR("Cannot destroy iterator, unknown type", IGRAPH_EINVAL); */
        break;
    }
}

int igraph_eit_as_vector(const igraph_eit_t *eit, igraph_vector_t *v) {

    long int i;

    IGRAPH_CHECK(igraph_vector_resize(v, IGRAPH_EIT_SIZE(*eit)));

    switch (eit->type) {
    case IGRAPH_EIT_SEQ:
        for (i = 0; i < IGRAPH_EIT_SIZE(*eit); i++) {
            VECTOR(*v)[i] = eit->start + i;
        }
        break;
    case IGRAPH_EIT_VECTOR:
    case IGRAPH_EIT_VECTORPTR:
        for (i = 0; i < IGRAPH_EIT_SIZE(*eit); i++) {
            VECTOR(*v)[i] = VECTOR(*eit->vec)[i];
        }
        break;
    default:
        IGRAPH_ERROR("Cannot convert to vector, unknown iterator type",
                     IGRAPH_EINVAL);
        break;
    }

    return 0;
}
