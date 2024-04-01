/* -*- mode: C -*-  */
/* vim:set ts=4 sw=4 sts=4 et: */
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

#include "igraph_structural.h"

#include "igraph_adjlist.h"
#include "igraph_interface.h"

/**
 * \function igraph_is_simple
 * \brief Decides whether the input graph is a simple graph.
 *
 * A graph is a simple graph if it does not contain loop edges and
 * multiple edges.
 *
 * \param graph The input graph.
 * \param res Pointer to a boolean constant, the result
 *     is stored here.
 * \return Error code.
 *
 * \sa \ref igraph_is_loop() and \ref igraph_is_multiple() to
 * find the loops and multiple edges, \ref igraph_simplify() to
 * get rid of them, or \ref igraph_has_multiple() to decide whether
 * there is at least one multiple edge.
 *
 * Time complexity: O(|V|+|E|).
 */
igraph_error_t igraph_is_simple(const igraph_t *graph, igraph_bool_t *res) {
    igraph_integer_t vc = igraph_vcount(graph);
    igraph_integer_t ec = igraph_ecount(graph);

    /* Is it already known whether the graph has loops or multi-edges? */
    igraph_bool_t known_loop, known_multi;

    /* If it is known, does the graph have them? */
    igraph_bool_t has_loop, has_multi;

    known_loop  = igraph_i_property_cache_has(graph, IGRAPH_PROP_HAS_LOOP);
    if (known_loop) {
        has_loop = igraph_i_property_cache_get_bool(graph, IGRAPH_PROP_HAS_LOOP);
        if (has_loop) {
            *res = false;
            return IGRAPH_SUCCESS;
        }
    }

    known_multi = igraph_i_property_cache_has(graph, IGRAPH_PROP_HAS_MULTI);
    if (known_multi) {
        has_multi = igraph_i_property_cache_get_bool(graph, IGRAPH_PROP_HAS_MULTI);
        if (has_multi) {
            *res = false;
            return IGRAPH_SUCCESS;
        }
    }

    if (known_loop && known_multi) {
        if (!has_loop && !has_multi) {
            *res = true;
            return IGRAPH_SUCCESS;
        }
    }

    /* Up to now, these variables were used to store the cache status.
     * From here on, we re-use them to store the outcome of explicit
     * checks. */

    known_loop = known_multi = false;
    has_loop = has_multi = false; /* be optimistic */

    if (vc == 0 || ec == 0) {
        *res = true;
        known_loop = known_multi = true;
    } else {
        igraph_vector_int_t neis;
        IGRAPH_VECTOR_INT_INIT_FINALLY(&neis, 0);
        for (igraph_integer_t i = 0; i < vc; i++) {
            IGRAPH_CHECK(igraph_neighbors(graph, &neis, i, IGRAPH_OUT));
            const igraph_integer_t n = igraph_vector_int_size(&neis);
            for (igraph_integer_t j = 0; j < n; j++) {
                if (VECTOR(neis)[j] == i) {
                    known_loop = true; has_loop = true; break;
                }
                if (j > 0 && VECTOR(neis)[j - 1] == VECTOR(neis)[j]) {
                    known_multi = true; has_multi = true; break;
                }
            }
        }
        *res = !has_loop && !has_multi;
        if (*res) {
            known_multi = known_loop = true;
        }
        igraph_vector_int_destroy(&neis);
        IGRAPH_FINALLY_CLEAN(1);
    }

    if (known_loop) {
        igraph_i_property_cache_set_bool_checked(graph, IGRAPH_PROP_HAS_LOOP, has_loop);
    }

    if (known_multi) {
        igraph_i_property_cache_set_bool_checked(graph, IGRAPH_PROP_HAS_MULTI, has_multi);
    }

    return IGRAPH_SUCCESS;
}


/**
 * \function igraph_has_multiple
 * \brief Check whether the graph has at least one multiple edge.
 *
 * An edge is a multiple edge if there is another
 * edge with the same head and tail vertices in the graph.
 *
 * </para><para>
 * The return value of this function is cached in the graph itself; calling
 * the function multiple times with no modifications to the graph in between
 * will return a cached value in O(1) time.
 *
 * \param graph The input graph.
 * \param res Pointer to a boolean variable, the result will be stored here.
 * \return Error code.
 *
 * \sa \ref igraph_count_multiple(), \ref igraph_is_multiple() and \ref igraph_simplify().
 *
 * Time complexity: O(e*d), e is the number of edges to check and d is the
 * average degree (out-degree in directed graphs) of the vertices at the
 * tail of the edges.
 *
 * \example examples/simple/igraph_has_multiple.c
 */
igraph_error_t igraph_has_multiple(const igraph_t *graph, igraph_bool_t *res) {
    igraph_integer_t vc = igraph_vcount(graph);
    igraph_integer_t ec = igraph_ecount(graph);
    igraph_bool_t directed = igraph_is_directed(graph);

    IGRAPH_RETURN_IF_CACHED_BOOL(graph, IGRAPH_PROP_HAS_MULTI, res);

    if (vc == 0 || ec == 0) {
        *res = false;
    } else {
        igraph_vector_int_t neis;
        igraph_integer_t i, j, n;
        igraph_bool_t found = false;
        IGRAPH_VECTOR_INT_INIT_FINALLY(&neis, 0);
        for (i = 0; i < vc && !found; i++) {
            IGRAPH_CHECK(igraph_neighbors(graph, &neis, i,
                                          IGRAPH_OUT));
            n = igraph_vector_int_size(&neis);
            for (j = 1; j < n; j++) {
                if (VECTOR(neis)[j - 1] == VECTOR(neis)[j]) {
                    /* If the graph is undirected, loop edges appear twice in the neighbor
                     * list, so check the next item as well */
                    if (directed) {
                        /* Directed, so this is a real multiple edge */
                        found = true; break;
                    } else if (VECTOR(neis)[j - 1] != i) {
                        /* Undirected, but not a loop edge */
                        found = true; break;
                    } else if (j < n - 1 && VECTOR(neis)[j] == VECTOR(neis)[j + 1]) {
                        /* Undirected, loop edge, multiple times */
                        found = true; break;
                    }
                }
            }
        }
        *res = found;
        igraph_vector_int_destroy(&neis);
        IGRAPH_FINALLY_CLEAN(1);
    }

    igraph_i_property_cache_set_bool_checked(graph, IGRAPH_PROP_HAS_MULTI, *res);

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_is_multiple
 * \brief Find the multiple edges in a graph.
 *
 * An edge is a multiple edge if there is another
 * edge with the same head and tail vertices in the graph.
 *
 * </para><para>
 * Note that this function returns true only for the second or more
 * appearances of the multiple edges.
 *
 * \param graph The input graph.
 * \param res Pointer to a boolean vector, the result will be stored
 *        here. It will be resized as needed.
 * \param es The edges to check. Supply \ref igraph_ess_all() if you want
 *        to check all edges.
 * \return Error code.
 *
 * \sa \ref igraph_count_multiple(), \ref igraph_has_multiple() and \ref igraph_simplify().
 *
 * Time complexity: O(e*d), e is the number of edges to check and d is the
 * average degree (out-degree in directed graphs) of the vertices at the
 * tail of the edges.
 *
 * \example examples/simple/igraph_is_multiple.c
 */
igraph_error_t igraph_is_multiple(const igraph_t *graph, igraph_vector_bool_t *res,
                       igraph_es_t es) {
    igraph_eit_t eit;
    igraph_integer_t i, j, n;
    igraph_lazy_inclist_t inclist;

    IGRAPH_CHECK(igraph_eit_create(graph, es, &eit));
    IGRAPH_FINALLY(igraph_eit_destroy, &eit);

    IGRAPH_CHECK(igraph_lazy_inclist_init(graph, &inclist, IGRAPH_OUT, IGRAPH_LOOPS_ONCE));
    IGRAPH_FINALLY(igraph_lazy_inclist_destroy, &inclist);

    IGRAPH_CHECK(igraph_vector_bool_resize(res, IGRAPH_EIT_SIZE(eit)));

    for (i = 0; !IGRAPH_EIT_END(eit); i++, IGRAPH_EIT_NEXT(eit)) {
        igraph_integer_t e = IGRAPH_EIT_GET(eit);
        igraph_integer_t from = IGRAPH_FROM(graph, e);
        igraph_integer_t to = IGRAPH_TO(graph, e);
        igraph_vector_int_t *neis = igraph_lazy_inclist_get(&inclist, from);

        IGRAPH_CHECK_OOM(neis, "Failed to query incident edges.");

        VECTOR(*res)[i] = false;

        n = igraph_vector_int_size(neis);
        for (j = 0; j < n; j++) {
            igraph_integer_t e2 = VECTOR(*neis)[j];
            igraph_integer_t to2 = IGRAPH_OTHER(graph, e2, from);
            if (to2 == to && e2 < e) {
                VECTOR(*res)[i] = true;
            }
        }
    }

    igraph_lazy_inclist_destroy(&inclist);
    igraph_eit_destroy(&eit);
    IGRAPH_FINALLY_CLEAN(2);
    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_count_multiple
 * \brief The multiplicity of some edges in a graph.
 *
 * An edge is called a multiple edge when there is one or more other
 * edge between the same two vertices. The multiplicity of an edge
 * is the number of edges between its endpoints.
 *
 * \param graph The input graph.
 * \param res Pointer to a vector, the result will be stored
 *        here. It will be resized as needed.
 * \param es The edges to check. Supply \ref igraph_ess_all() if you want
 *        to check all edges.
 * \return Error code.
 *
 * \sa \ref igraph_count_multiple_1() if you only need the multiplicity of a
 * single edge; \ref igraph_is_multiple() if you are only interested in whether
 * the graph has at least one edge with multiplicity greater than one;
 * \ref igraph_simplify() to ensure that the graph has no multiple edges.
 *
 * Time complexity: O(E d), E is the number of edges to check and d is the
 * average degree (out-degree in directed graphs) of the vertices at the
 * tail of the edges.
 */
igraph_error_t igraph_count_multiple(const igraph_t *graph, igraph_vector_int_t *res, igraph_es_t es) {
    igraph_eit_t eit;
    igraph_integer_t i, j, n;
    igraph_lazy_adjlist_t adjlist;

    IGRAPH_CHECK(igraph_eit_create(graph, es, &eit));
    IGRAPH_FINALLY(igraph_eit_destroy, &eit);
    IGRAPH_CHECK(igraph_lazy_adjlist_init(graph, &adjlist, IGRAPH_OUT, IGRAPH_LOOPS_ONCE, IGRAPH_MULTIPLE));
    IGRAPH_FINALLY(igraph_lazy_adjlist_destroy, &adjlist);

    IGRAPH_CHECK(igraph_vector_int_resize(res, IGRAPH_EIT_SIZE(eit)));

    for (i = 0; !IGRAPH_EIT_END(eit); i++, IGRAPH_EIT_NEXT(eit)) {
        igraph_integer_t e = IGRAPH_EIT_GET(eit);
        igraph_integer_t from = IGRAPH_FROM(graph, e);
        igraph_integer_t to = IGRAPH_TO(graph, e);
        igraph_vector_int_t *neis = igraph_lazy_adjlist_get(&adjlist, from);

        IGRAPH_CHECK_OOM(neis, "Failed to query adjacent vertices.");

        VECTOR(*res)[i] = 0;

        n = igraph_vector_int_size(neis);
        for (j = 0; j < n; j++) {
            if (VECTOR(*neis)[j] == to) {
                VECTOR(*res)[i]++;
            }
        }
    }

    igraph_lazy_adjlist_destroy(&adjlist);
    igraph_eit_destroy(&eit);
    IGRAPH_FINALLY_CLEAN(2);

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_count_multiple_1
 * \brief The multiplicity of a single edge in a graph.
 *
 * \param graph The input graph.
 * \param res Pointer to an iteger, the result will be stored here.
 * \param eid The ID of the edge to check.
 * \return Error code.
 *
 * \sa \ref igraph_count_multiple() if you need the multiplicity of multiple
 * edges; \ref igraph_is_multiple() if you are only interested in whether the
 * graph has at least one edge with multiplicity greater than one;
 * \ref igraph_simplify() to ensure that the graph has no multiple edges.
 *
 * Time complexity: O(d), where d is the out-degree of the tail of the edge.
 */
igraph_error_t igraph_count_multiple_1(const igraph_t *graph, igraph_integer_t *res, igraph_integer_t eid)
{
    igraph_integer_t i, n, count;
    igraph_integer_t from = IGRAPH_FROM(graph, eid);
    igraph_integer_t to = IGRAPH_TO(graph, eid);
    igraph_vector_int_t vids;

    IGRAPH_VECTOR_INT_INIT_FINALLY(&vids, 0);
    IGRAPH_CHECK(igraph_neighbors(graph, &vids, from, IGRAPH_OUT));

    count = 0;
    n = igraph_vector_int_size(&vids);
    for (i = 0; i < n; i++) {
        if (VECTOR(vids)[i] == to) {
            count++;
        }
    }

    igraph_vector_int_destroy(&vids);
    IGRAPH_FINALLY_CLEAN(1);

    *res = count;

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_is_mutual
 * \brief Check whether some edges of a directed graph are mutual.
 *
 * An (A,B) non-loop directed edge is mutual if the graph contains
 * the (B,A) edge too. Whether directed self-loops are considered mutual
 * is controlled by the \p loops parameter.
 *
 * </para><para>
 * An undirected graph only has mutual edges, by definition.
 *
 * </para><para>
 * Edge multiplicity is not considered here, e.g. if there are two
 * (A,B) edges and one (B,A) edge, then all three are considered to be
 * mutual.
 *
 * \param graph The input graph.
 * \param res Pointer to an initialized vector, the result is stored
 *        here.
 * \param es The sequence of edges to check. Supply
 *        \ref igraph_ess_all() to check all edges.
 * \param loops Boolean, whether to consider directed self-loops
 *        to be mutual.
 * \return Error code.
 *
 * Time complexity: O(n log(d)), n is the number of edges supplied, d
 * is the maximum in-degree of the vertices that are targets of the
 * supplied edges. An upper limit of the time complexity is O(n log(|E|)),
 * |E| is the number of edges in the graph.
 */
igraph_error_t igraph_is_mutual(const igraph_t *graph, igraph_vector_bool_t *res, igraph_es_t es, igraph_bool_t loops) {

    igraph_eit_t eit;
    igraph_lazy_adjlist_t adjlist;
    igraph_integer_t i;

    /* How many edges do we have? */
    IGRAPH_CHECK(igraph_eit_create(graph, es, &eit));
    IGRAPH_FINALLY(igraph_eit_destroy, &eit);
    IGRAPH_CHECK(igraph_vector_bool_resize(res, IGRAPH_EIT_SIZE(eit)));

    /* An undirected graph has mutual edges by definition,
       res is already properly resized */
    if (! igraph_is_directed(graph)) {
        igraph_vector_bool_fill(res, true);
        igraph_eit_destroy(&eit);
        IGRAPH_FINALLY_CLEAN(1);
        return IGRAPH_SUCCESS;
    }

    IGRAPH_CHECK(igraph_lazy_adjlist_init(graph, &adjlist, IGRAPH_OUT, IGRAPH_LOOPS_ONCE, IGRAPH_MULTIPLE));
    IGRAPH_FINALLY(igraph_lazy_adjlist_destroy, &adjlist);

    for (i = 0; ! IGRAPH_EIT_END(eit); i++, IGRAPH_EIT_NEXT(eit)) {
        igraph_integer_t edge = IGRAPH_EIT_GET(eit);
        igraph_integer_t from = IGRAPH_FROM(graph, edge);
        igraph_integer_t to = IGRAPH_TO(graph, edge);

        if (from == to) {
            VECTOR(*res)[i] = loops;
            continue; /* no need to do binsearch for self-loops */
        }

        /* Check whether there is a to->from edge, search for from in the
           out-list of to */
        igraph_vector_int_t *neis = igraph_lazy_adjlist_get(&adjlist, to);
        IGRAPH_CHECK_OOM(neis, "Failed to query neighbors.");
        VECTOR(*res)[i] = igraph_vector_int_binsearch2(neis, from);
    }

    igraph_lazy_adjlist_destroy(&adjlist);
    igraph_eit_destroy(&eit);
    IGRAPH_FINALLY_CLEAN(2);

    return IGRAPH_SUCCESS;
}


/**
 * \function igraph_has_mutual
 * \brief Check whether a directed graph has any mutual edges.
 *
 * An (A,B) non-loop directed edge is mutual if the graph contains
 * the (B,A) edge too. Whether directed self-loops are considered mutual
 * is controlled by the \p loops parameter.
 *
 * </para><para>
 * In undirected graphs, all edges are considered mutual by definition.
 * Thus for undirected graph, this function returns false only when there
 * are no edges.
 *
 * </para><para>
 * To check whether a graph is an oriented graph, use this function in
 * conjunction with \ref igraph_is_directed().
 *
 * \param graph The input graph.
 * \param res Pointer to a boolean, the result will be stored here.
 * \param loops Boolean, whether to consider directed self-loops
 *        to be mutual.
 * \return Error code.
 *
 * Time complexity: O(|E| log(d)) where d is the maximum in-degree.
 */
igraph_error_t igraph_has_mutual(const igraph_t *graph, igraph_bool_t *res, igraph_bool_t loops) {
    igraph_integer_t no_of_edges = igraph_ecount(graph);
    igraph_lazy_adjlist_t adjlist;

    if (! igraph_is_directed(graph)) {
        /* In undirected graphs, all edges are considered mutual, so we just check
         * if there are any edges. */
        *res = no_of_edges > 0;
        return IGRAPH_SUCCESS;
    }

    if (igraph_i_property_cache_has(graph, IGRAPH_PROP_HAS_MUTUAL)) {
        if (igraph_i_property_cache_get_bool(graph, IGRAPH_PROP_HAS_MUTUAL)) {
            /* we know that the graph has at least one mutual non-loop edge
             * (because the cache only stores non-loop edges) */
            *res = true;
            return IGRAPH_SUCCESS;
        } else if (loops) {
            /* no non-loop mutual edges, but maybe we have loops? */
            return igraph_has_loop(graph, res);
        } else {
            /* no non-loop mutual edges, and loops are not to be treated as mutual */
            *res = false;
            return IGRAPH_SUCCESS;
        }
    }

    IGRAPH_CHECK(igraph_lazy_adjlist_init(graph, &adjlist, IGRAPH_OUT, IGRAPH_LOOPS_ONCE, IGRAPH_MULTIPLE));
    IGRAPH_FINALLY(igraph_lazy_adjlist_destroy, &adjlist);

    *res = false; /* assume no mutual edges */
    for (igraph_integer_t edge=0; edge < no_of_edges; edge++) {
        igraph_integer_t from = IGRAPH_FROM(graph, edge);
        igraph_integer_t to = IGRAPH_TO(graph, edge);

        if (from == to) {
            if (loops) {
                *res = true;
                break;
            }
            continue; /* no need to do binsearch for self-loops */
        }

        /* Check whether there is a to->from edge, search for from in the
           out-list of to */
        igraph_vector_int_t *neis = igraph_lazy_adjlist_get(&adjlist, to);
        IGRAPH_CHECK_OOM(neis, "Failed to query neighbors.");
        if (igraph_vector_int_binsearch2(neis, from)) {
            *res = true;
            break;
        }
    }

    igraph_lazy_adjlist_destroy(&adjlist);
    IGRAPH_FINALLY_CLEAN(1);

    /* cache the result if loops are not treated as mutual */
    if (!loops) {
        igraph_i_property_cache_set_bool_checked(graph, IGRAPH_PROP_HAS_MUTUAL, *res);
    }

    return IGRAPH_SUCCESS;
}
