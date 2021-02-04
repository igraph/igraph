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
 * </para><para>
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
int igraph_is_simple(const igraph_t *graph, igraph_bool_t *res) {
    long int vc = igraph_vcount(graph);
    long int ec = igraph_ecount(graph);

    if (vc == 0 || ec == 0) {
        *res = 1;
    } else {
        igraph_vector_t neis;
        long int i, j, n;
        igraph_bool_t found = 0;
        IGRAPH_VECTOR_INIT_FINALLY(&neis, 0);
        for (i = 0; i < vc; i++) {
            IGRAPH_CHECK(igraph_neighbors(graph, &neis, (igraph_integer_t) i, IGRAPH_OUT));
            n = igraph_vector_size(&neis);
            for (j = 0; j < n; j++) {
                if (VECTOR(neis)[j] == i) {
                    found = 1; break;
                }
                if (j > 0 && VECTOR(neis)[j - 1] == VECTOR(neis)[j]) {
                    found = 1; break;
                }
            }
        }
        *res = !found;
        igraph_vector_destroy(&neis);
        IGRAPH_FINALLY_CLEAN(1);
    }

    return 0;
}


/**
 * \function igraph_has_multiple
 * \brief Check whether the graph has at least one multiple edge.
 *
 * </para><para>
 * An edge is a multiple edge if there is another
 * edge with the same head and tail vertices in the graph.
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
int igraph_has_multiple(const igraph_t *graph, igraph_bool_t *res) {
    long int vc = igraph_vcount(graph);
    long int ec = igraph_ecount(graph);
    igraph_bool_t directed = igraph_is_directed(graph);

    if (vc == 0 || ec == 0) {
        *res = 0;
    } else {
        igraph_vector_t neis;
        long int i, j, n;
        igraph_bool_t found = 0;
        IGRAPH_VECTOR_INIT_FINALLY(&neis, 0);
        for (i = 0; i < vc && !found; i++) {
            IGRAPH_CHECK(igraph_neighbors(graph, &neis, (igraph_integer_t) i,
                                          IGRAPH_OUT));
            n = igraph_vector_size(&neis);
            for (j = 1; j < n; j++) {
                if (VECTOR(neis)[j - 1] == VECTOR(neis)[j]) {
                    /* If the graph is undirected, loop edges appear twice in the neighbor
                     * list, so check the next item as well */
                    if (directed) {
                        /* Directed, so this is a real multiple edge */
                        found = 1; break;
                    } else if (VECTOR(neis)[j - 1] != i) {
                        /* Undirected, but not a loop edge */
                        found = 1; break;
                    } else if (j < n - 1 && VECTOR(neis)[j] == VECTOR(neis)[j + 1]) {
                        /* Undirected, loop edge, multiple times */
                        found = 1; break;
                    }
                }
            }
        }
        *res = found;
        igraph_vector_destroy(&neis);
        IGRAPH_FINALLY_CLEAN(1);
    }

    return 0;
}

/**
 * \function igraph_is_multiple
 * \brief Find the multiple edges in a graph.
 *
 * </para><para>
 * An edge is a multiple edge if there is another
 * edge with the same head and tail vertices in the graph.
 *
 * </para><para>
 * Note that this function returns true only for the second or more
 * appearances of the multiple edges.
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
int igraph_is_multiple(const igraph_t *graph, igraph_vector_bool_t *res,
                       igraph_es_t es) {
    igraph_eit_t eit;
    long int i, j, n;
    igraph_lazy_inclist_t inclist;

    IGRAPH_CHECK(igraph_eit_create(graph, es, &eit));
    IGRAPH_FINALLY(igraph_eit_destroy, &eit);

    IGRAPH_CHECK(igraph_lazy_inclist_init(graph, &inclist, IGRAPH_OUT, IGRAPH_LOOPS_ONCE));
    IGRAPH_FINALLY(igraph_lazy_inclist_destroy, &inclist);

    IGRAPH_CHECK(igraph_vector_bool_resize(res, IGRAPH_EIT_SIZE(eit)));

    for (i = 0; !IGRAPH_EIT_END(eit); i++, IGRAPH_EIT_NEXT(eit)) {
        long int e = IGRAPH_EIT_GET(eit);
        long int from = IGRAPH_FROM(graph, e);
        long int to = IGRAPH_TO(graph, e);
        igraph_vector_int_t *neis =
            igraph_lazy_inclist_get(&inclist, (igraph_integer_t) from);

        if (neis == 0) {
            /* Most likely out of memory */
            IGRAPH_ERROR("Out of memory while building lazy incidence list", IGRAPH_ENOMEM);
        }

        VECTOR(*res)[i] = 0;

        n = igraph_vector_int_size(neis);
        for (j = 0; j < n; j++) {
            long int e2 = (long int) VECTOR(*neis)[j];
            long int to2 = IGRAPH_OTHER(graph, e2, from);
            if (to2 == to && e2 < e) {
                VECTOR(*res)[i] = 1;
            }
        }
    }

    igraph_lazy_inclist_destroy(&inclist);
    igraph_eit_destroy(&eit);
    IGRAPH_FINALLY_CLEAN(2);
    return 0;
}


/**
 * \function igraph_count_multiple
 * \brief Count the number of appearances of the edges in a graph.
 *
 * </para><para>
 * If the graph has no multiple edges then the result vector will be
 * filled with ones.
 * (An edge is a multiple edge if there is another
 * edge with the same head and tail vertices in the graph.)
 *
 * </para><para>
 * \param graph The input graph.
 * \param res Pointer to a vector, the result will be stored
 *        here. It will be resized as needed.
 * \param es The edges to check. Supply \ref igraph_ess_all() if you want
 *        to check all edges.
 * \return Error code.
 *
 * \sa \ref igraph_is_multiple() and \ref igraph_simplify().
 *
 * Time complexity: O(E d), E is the number of edges to check and d is the
 * average degree (out-degree in directed graphs) of the vertices at the
 * tail of the edges.
 */
int igraph_count_multiple(const igraph_t *graph, igraph_vector_t *res, igraph_es_t es) {
    igraph_eit_t eit;
    long int i, j, n;
    igraph_lazy_inclist_t inclist;

    IGRAPH_CHECK(igraph_eit_create(graph, es, &eit));
    IGRAPH_FINALLY(igraph_eit_destroy, &eit);
    IGRAPH_CHECK(igraph_lazy_inclist_init(graph, &inclist, IGRAPH_OUT, IGRAPH_LOOPS_ONCE));
    IGRAPH_FINALLY(igraph_lazy_inclist_destroy, &inclist);

    IGRAPH_CHECK(igraph_vector_resize(res, IGRAPH_EIT_SIZE(eit)));

    for (i = 0; !IGRAPH_EIT_END(eit); i++, IGRAPH_EIT_NEXT(eit)) {
        long int e = IGRAPH_EIT_GET(eit);
        long int from = IGRAPH_FROM(graph, e);
        long int to = IGRAPH_TO(graph, e);
        igraph_vector_int_t *neis =
            igraph_lazy_inclist_get(&inclist, (igraph_integer_t) from);
        
        if (neis == 0) {
            /* Most likely out of memory */
            IGRAPH_ERROR("Out of memory while building lazy incidence list", IGRAPH_ENOMEM);
        }

        VECTOR(*res)[i] = 0;
        
        n = igraph_vector_int_size(neis);
        for (j = 0; j < n; j++) {
            long int e2 = (long int) VECTOR(*neis)[j];
            long int to2 = IGRAPH_OTHER(graph, e2, from);
            if (to2 == to) {
                VECTOR(*res)[i] += 1;
            }
        }
    }

    igraph_lazy_inclist_destroy(&inclist);
    igraph_eit_destroy(&eit);
    IGRAPH_FINALLY_CLEAN(2);

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_is_mutual
 * Check whether the edges of a directed graph are mutual.
 *
 * An (A,B) edge is mutual if the graph contains the (B,A) edge, too.
 * </para>
 *
 * <para>An undirected graph only has mutual edges, by definition.
 * </para>
 *
 * <para>Edge multiplicity is not considered here, e.g. if there are two
 * (A,B) edges and one (B,A) edge, then all three are considered to be
 * mutual.
 *
 * \param graph The input graph.
 * \param res Pointer to an initialized vector, the result is stored
 *        here.
 * \param es The sequence of edges to check. Supply
 *        <code>igraph_ess_all()</code> for all edges, see \ref
 *        igraph_ess_all().
 * \return Error code.
 *
 * Time complexity: O(n log(d)), n is the number of edges supplied, d
 * is the maximum in-degree of the vertices that are targets of the
 * supplied edges. An upper limit of the time complexity is O(n log(|E|)),
 * |E| is the number of edges in the graph.
 */
int igraph_is_mutual(igraph_t *graph, igraph_vector_bool_t *res, igraph_es_t es) {

    igraph_eit_t eit;
    igraph_lazy_adjlist_t adjlist;
    long int i;

    /* How many edges do we have? */
    IGRAPH_CHECK(igraph_eit_create(graph, es, &eit));
    IGRAPH_FINALLY(igraph_eit_destroy, &eit);
    IGRAPH_CHECK(igraph_vector_bool_resize(res, IGRAPH_EIT_SIZE(eit)));

    /* An undirected graph has mutual edges by definition,
       res is already properly resized */
    if (! igraph_is_directed(graph)) {
        igraph_vector_bool_fill(res, 1);
        igraph_eit_destroy(&eit);
        IGRAPH_FINALLY_CLEAN(1);
        return 0;
    }

    IGRAPH_CHECK(igraph_lazy_adjlist_init(graph, &adjlist, IGRAPH_OUT, IGRAPH_LOOPS_ONCE, IGRAPH_MULTIPLE));
    IGRAPH_FINALLY(igraph_lazy_adjlist_destroy, &adjlist);

    for (i = 0; ! IGRAPH_EIT_END(eit); i++, IGRAPH_EIT_NEXT(eit)) {
        long int edge = IGRAPH_EIT_GET(eit);
        long int from = IGRAPH_FROM(graph, edge);
        long int to = IGRAPH_TO(graph, edge);

        /* Check whether there is a to->from edge, search for from in the
           out-list of to. We don't search an empty vector, because
           vector_binsearch seems to have a bug with this. */
        igraph_vector_int_t *neis = igraph_lazy_adjlist_get(&adjlist,
                                (igraph_integer_t) to);
        if (igraph_vector_int_empty(neis)) {
            VECTOR(*res)[i] = 0;
        } else {
            VECTOR(*res)[i] = igraph_vector_int_binsearch2(neis, from);
        }
    }

    igraph_lazy_adjlist_destroy(&adjlist);
    igraph_eit_destroy(&eit);
    IGRAPH_FINALLY_CLEAN(2);

    return 0;
}
