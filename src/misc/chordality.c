/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2008-2021  The igraph development team <igraph@igraph.org>

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

#include "igraph_structural.h"

#include "igraph_adjlist.h"
#include "igraph_interface.h"

/**
 * \function igraph_maximum_cardinality_search
 * \brief Maximum cardinality search.
 *
 * This function implements the maximum cardinality search algorithm.
 * It computes a rank \p alpha for each vertex, such that visiting
 * vertices in decreasing rank order corresponds to always choosing
 * the vertex with the most already visited neighbors as the next one
 * to visit.
 *
 * </para><para>
 * Maximum cardinality search is useful in deciding the chordality
 * of a graph. A graph is chordal if and only if any two neighbors
 * of a vertex which are higher in rank than it are connected to
 * each other.
 *
 * </para><para>
 * References:
 *
 * </para><para>
 * Robert E Tarjan and Mihalis Yannakakis: Simple linear-time
 * algorithms to test chordality of graphs, test acyclicity of
 * hypergraphs, and selectively reduce acyclic hypergraphs.
 * SIAM Journal of Computation 13, 566--579, 1984.
 * https://doi.org/10.1137/0213035
 *
 * \param graph The input graph. Edge directions will be ignored.
 * \param alpha Pointer to an initialized vector, the result is stored here.
 *   It will be resized, as needed. Upon return it contains
 *   the rank of the each vertex in the range 0 to <code>n - 1</code>,
 *   where \c n is the number of vertices.
 * \param alpham1 Pointer to an initialized vector or a \c NULL
 *   pointer. If not \c NULL, then the inverse of \p alpha is stored
 *   here. In other words, the elements of \p alpham1 are vertex IDs
 *   in reverse maximum cardinality search order.
 * \return Error code.
 *
 * Time complexity: O(|V|+|E|), linear in terms of the number of
 * vertices and edges.
 *
 * \sa \ref igraph_is_chordal().
 */

igraph_error_t igraph_maximum_cardinality_search(const igraph_t *graph,
                                      igraph_vector_int_t *alpha,
                                      igraph_vector_int_t *alpham1) {

    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_vector_int_t size;
    igraph_vector_int_t head, next, prev; /* doubly linked list with head */
    igraph_integer_t i;
    igraph_adjlist_t adjlist;

    /***************/
    /* local j, v; */
    /***************/

    igraph_integer_t j, v;

    if (no_of_nodes == 0) {
        igraph_vector_int_clear(alpha);
        if (alpham1) {
            igraph_vector_int_clear(alpham1);
        }
        return IGRAPH_SUCCESS;
    }

    IGRAPH_VECTOR_INT_INIT_FINALLY(&size, no_of_nodes);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&head, no_of_nodes);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&next, no_of_nodes);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&prev, no_of_nodes);

    IGRAPH_CHECK(igraph_adjlist_init(graph, &adjlist, IGRAPH_ALL, IGRAPH_NO_LOOPS, IGRAPH_NO_MULTIPLE));
    IGRAPH_FINALLY(igraph_adjlist_destroy, &adjlist);

    IGRAPH_CHECK(igraph_vector_int_resize(alpha, no_of_nodes));
    if (alpham1) {
        IGRAPH_CHECK(igraph_vector_int_resize(alpham1, no_of_nodes));
    }

    /***********************************************/
    /* for i in [0,n-1] -> set(i) := emptyset rof; */
    /***********************************************/

    /* nothing to do, 'head' contains all zeros */

    /*********************************************************/
    /* for v in vertices -> size(v):=0; add v to set(0) rof; */
    /*********************************************************/

    VECTOR(head)[0] = 1;
    for (v = 0; v < no_of_nodes; v++) {
        VECTOR(next)[v] = v + 2;
        VECTOR(prev)[v] = v;
    }
    VECTOR(next)[no_of_nodes - 1] = 0;
    /* size is already all zero */

    /***************/
    /* i:=n; j:=0; */
    /***************/

    i = no_of_nodes; j = 0;

    /**************/
    /* do i>=1 -> */
    /**************/

    while (i >= 1) {
        igraph_integer_t x, k, len;
        igraph_vector_int_t *neis;

        /********************************/
        /* v :=  delete any from set(j) */
        /********************************/

        v = VECTOR(head)[j] - 1;
        x = VECTOR(next)[v];
        VECTOR(head)[j] = x;
        if (x != 0) {
            VECTOR(prev)[x - 1] = 0;
        }

        /*************************************************/
        /* alpha(v) := i; alpham1(i) := v; size(v) := -1 */
        /*************************************************/

        VECTOR(*alpha)[v] = i - 1;
        if (alpham1) {
            VECTOR(*alpham1)[i - 1] = v;
        }
        VECTOR(size)[v] = -1;

        /********************************************/
        /* for {v,w} in E such that size(w) >= 0 -> */
        /********************************************/

        neis = igraph_adjlist_get(&adjlist, v);
        len = igraph_vector_int_size(neis);
        for (k = 0; k < len; k++) {
            igraph_integer_t w = VECTOR(*neis)[k];
            igraph_integer_t ws = VECTOR(size)[w];
            if (ws >= 0) {

                /******************************/
                /* delete w from set(size(w)) */
                /******************************/

                igraph_integer_t nw = VECTOR(next)[w];
                igraph_integer_t pw = VECTOR(prev)[w];
                if (nw != 0) {
                    VECTOR(prev)[nw - 1] = pw;
                }
                if (pw != 0) {
                    VECTOR(next)[pw - 1] = nw;
                } else {
                    VECTOR(head)[ws] = nw;
                }

                /******************************/
                /* size(w) := size(w)+1       */
                /******************************/

                VECTOR(size)[w] += 1;

                /******************************/
                /* add w to set(size(w))      */
                /******************************/

                ws = VECTOR(size)[w];
                nw = VECTOR(head)[ws];
                VECTOR(next)[w] = nw;
                VECTOR(prev)[w] = 0;
                if (nw != 0) {
                    VECTOR(prev)[nw - 1] = w + 1;
                }
                VECTOR(head)[ws] = w + 1;

            }
        }

        /***********************/
        /* i := i-1; j := j+1; */
        /***********************/

        i -= 1;
        j += 1;

        /*********************************************/
        /* do j>=0 and set(j)=emptyset -> j:=j-1; od */
        /*********************************************/

        if (j < no_of_nodes) {
            while (j >= 0 && VECTOR(head)[j] == 0) {
                j--;
            }
        }
    }

    igraph_adjlist_destroy(&adjlist);
    igraph_vector_int_destroy(&prev);
    igraph_vector_int_destroy(&next);
    igraph_vector_int_destroy(&head);
    igraph_vector_int_destroy(&size);
    IGRAPH_FINALLY_CLEAN(5);

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_is_chordal
 * \brief Decides whether a graph is chordal.
 *
 * A graph is chordal if each of its cycles of four or more nodes
 * has a chord, i.e. an edge joining two nodes that are not
 * adjacent in the cycle. An equivalent definition is that any
 * chordless cycles have at most three nodes.
 *
 * If either \p alpha or \p alpham1 is given, then the other is
 * calculated by taking simply the inverse. If neither are given,
 * then \ref igraph_maximum_cardinality_search() is called to calculate
 * them.
 *
 * \param graph The input graph. Edge directions will be ignored.
 * \param alpha Either an alpha vector coming from
 *    \ref igraph_maximum_cardinality_search() (on the same graph), or a
 *    \c NULL pointer.
 * \param alpham1 Either an inverse alpha vector coming from \ref
 *    igraph_maximum_cardinality_search() (on the same graph) or a \c NULL
 *    pointer.
 * \param chordal Pointer to a boolean. If not NULL the result is stored here.
 * \param fill_in Pointer to an initialized vector, or a \c NULL
 *    pointer. If not a \c NULL pointer, then the fill-in, also called the
 *    chordal completion of the graph is stored here.
 *    The chordal completion is a set of edges that are needed to
 *    make the graph chordal. The vector is resized as needed.
 *    Note that the chordal completion returned by this function may not
 *    be minimal, i.e. some of the returned fill-in edges may not be needed
 *    to make the graph chordal.
 * \param newgraph Pointer to an uninitialized graph, or a \c NULL
 *   pointer. If not a null pointer, then a new triangulated graph is
 *   created here. This essentially means adding the fill-in edges to
 *   the original graph.
 * \return Error code.
 *
 * Time complexity: O(n).
 *
 * \sa \ref igraph_maximum_cardinality_search().
 */

igraph_error_t igraph_is_chordal(const igraph_t *graph,
                      const igraph_vector_int_t *alpha,
                      const igraph_vector_int_t *alpham1,
                      igraph_bool_t *chordal,
                      igraph_vector_int_t *fill_in,
                      igraph_t *newgraph) {

    const igraph_integer_t no_of_nodes = igraph_vcount(graph);
    const igraph_vector_int_t *my_alpha = alpha, *my_alpham1 = alpham1;
    igraph_vector_int_t v_alpha, v_alpham1;
    igraph_vector_int_t f, index;
    igraph_integer_t i;
    igraph_adjlist_t adjlist;
    igraph_vector_int_t mark;
    igraph_bool_t calc_edges = fill_in || newgraph;
    igraph_vector_int_t *my_fill_in = fill_in, v_fill_in;

    /*****************/
    /* local v, w, x */
    /*****************/

    igraph_integer_t v, w, x;

    if (alpha && (igraph_vector_int_size(alpha) != no_of_nodes)) {
        IGRAPH_ERRORF("Alpha vector size (%" IGRAPH_PRId ") not equal to number of nodes (%" IGRAPH_PRId ").",
                      IGRAPH_EINVAL, igraph_vector_int_size(alpha), no_of_nodes);
    }

    if (alpham1 && (igraph_vector_int_size(alpham1) != no_of_nodes)) {
        IGRAPH_ERRORF("Inverse alpha vector size (%" IGRAPH_PRId ") not equal to number of nodes (%" IGRAPH_PRId ").",
                      IGRAPH_EINVAL, igraph_vector_int_size(alpham1), no_of_nodes);
    }

    if (!chordal && !calc_edges) {
        /* Nothing to calculate */
        return IGRAPH_SUCCESS;
    }

    if (!alpha && !alpham1) {
        IGRAPH_VECTOR_INT_INIT_FINALLY(&v_alpha, no_of_nodes);
        my_alpha = &v_alpha;
        IGRAPH_VECTOR_INT_INIT_FINALLY(&v_alpham1, no_of_nodes);
        my_alpham1 = &v_alpham1;
        IGRAPH_CHECK(igraph_maximum_cardinality_search(graph,
                     (igraph_vector_int_t*) my_alpha,
                     (igraph_vector_int_t*) my_alpham1));
    } else if (alpha && !alpham1) {
        igraph_integer_t v;
        IGRAPH_VECTOR_INT_INIT_FINALLY(&v_alpham1, no_of_nodes);
        my_alpham1 = &v_alpham1;
        for (v = 0; v < no_of_nodes; v++) {
            igraph_integer_t i = VECTOR(*my_alpha)[v];
            VECTOR(*my_alpham1)[i] = v;
        }
    } else if (!alpha && alpham1) {
        igraph_integer_t i;
        IGRAPH_VECTOR_INT_INIT_FINALLY(&v_alpha, no_of_nodes);
        my_alpha = &v_alpha;
        for (i = 0; i < no_of_nodes; i++) {
            igraph_integer_t v = VECTOR(*my_alpham1)[i];
            VECTOR(*my_alpha)[v] = i;
        }
    }

    if (!fill_in && newgraph) {
        IGRAPH_VECTOR_INT_INIT_FINALLY(&v_fill_in, 0);
        my_fill_in = &v_fill_in;
    }

    IGRAPH_VECTOR_INT_INIT_FINALLY(&f, no_of_nodes);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&index, no_of_nodes);
    IGRAPH_CHECK(igraph_adjlist_init(graph, &adjlist, IGRAPH_ALL, IGRAPH_NO_LOOPS, IGRAPH_NO_MULTIPLE));
    IGRAPH_FINALLY(igraph_adjlist_destroy, &adjlist);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&mark, no_of_nodes);
    if (my_fill_in) {
        igraph_vector_int_clear(my_fill_in);
    }

    if (chordal) {
        *chordal = true;
    }

    /*********************/
    /* for i in [1,n] -> */
    /*********************/

    for (i = 0; i < no_of_nodes; i++) {
        igraph_vector_int_t *neis;
        igraph_integer_t j, len;

        /**********************************************/
        /* w := alpham1(i); f(w) := w; index(w) := i; */
        /**********************************************/

        w = VECTOR(*my_alpham1)[i];
        VECTOR(f)[w] = w;
        VECTOR(index)[w] = i;

        /******************************************/
        /* for {v,w} in E such that alpha(v)<i -> */
        /******************************************/

        neis = igraph_adjlist_get(&adjlist, w);
        len = igraph_vector_int_size(neis);
        for (j = 0; j < len; j++) {
            v = VECTOR(*neis)[j];
            VECTOR(mark)[v] = w + 1;
        }

        for (j = 0; j < len; j++) {
            v = VECTOR(*neis)[j];
            if (VECTOR(*my_alpha)[v] >= i) {
                continue;
            }

            /**********/
            /* x := v */
            /**********/

            x = v;

            /********************/
            /* do index(x)<i -> */
            /********************/

            while (VECTOR(index)[x] < i) {

                /******************/
                /* index(x) := i; */
                /******************/

                VECTOR(index)[x] = i;

                /**********************************/
                /* add {x,w} to E union F(alpha); */
                /**********************************/

                if (VECTOR(mark)[x] != w + 1) {

                    if (chordal) {
                        *chordal = false;
                    }

                    if (my_fill_in) {
                        IGRAPH_CHECK(igraph_vector_int_push_back(my_fill_in, x));
                        IGRAPH_CHECK(igraph_vector_int_push_back(my_fill_in, w));
                    }

                    if (!calc_edges) {
                        /* make sure that we exit from all loops */
                        i = no_of_nodes;
                        j = len;
                        break;
                    }
                }

                /*************/
                /* x := f(x) */
                /*************/

                x = VECTOR(f)[x];

            } /* while (VECTOR(index)[x] < i) */

            /*****************************/
            /* if (f(x)=x -> f(x):=w; fi */
            /*****************************/

            if (VECTOR(f)[x] == x) {
                VECTOR(f)[x] = w;
            }
        }
    }

    igraph_vector_int_destroy(&mark);
    igraph_adjlist_destroy(&adjlist);
    igraph_vector_int_destroy(&index);
    igraph_vector_int_destroy(&f);
    IGRAPH_FINALLY_CLEAN(4);

    if (newgraph) {
        IGRAPH_CHECK(igraph_copy(newgraph, graph));
        IGRAPH_FINALLY(igraph_destroy, newgraph);
        IGRAPH_CHECK(igraph_add_edges(newgraph, my_fill_in, 0));
        IGRAPH_FINALLY_CLEAN(1);
    }

    if (!fill_in && newgraph) {
        igraph_vector_int_destroy(&v_fill_in);
        IGRAPH_FINALLY_CLEAN(1);
    }

    if (!alpha && !alpham1) {
        igraph_vector_int_destroy(&v_alpham1);
        igraph_vector_int_destroy(&v_alpha);
        IGRAPH_FINALLY_CLEAN(2);
    } else if (alpha && !alpham1) {
        igraph_vector_int_destroy(&v_alpham1);
        IGRAPH_FINALLY_CLEAN(1);
    } else if (!alpha && alpham1) {
        igraph_vector_int_destroy(&v_alpha);
        IGRAPH_FINALLY_CLEAN(1);
    }

    return IGRAPH_SUCCESS;
}
