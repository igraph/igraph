/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2006-2012  Gabor Csardi <csardi.gabor@gmail.com>
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

#include "igraph_topology.h"

#include "igraph_adjlist.h"
#include "igraph_interface.h"
#include "igraph_memory.h"
#include "igraph_stack.h"

#include "core/interruption.h"

/**
 * \section about_vf2
 *
 * <para>
 * The VF2 algorithm can search for a subgraph in a larger graph, or check if two
 * graphs are isomorphic. See P. Foggia, C. Sansone, M. Vento, An Improved algorithm for
 * matching large graphs, Proc. of the 3rd IAPR-TC-15 International
 * Workshop on Graph-based Representations, Italy, 2001.
 * </para>
 *
 * <para>
 * VF2 supports both vertex and edge-colored graphs, as well as custom vertex or edge
 * compatibility functions.
 * </para>
 *
 * <para>
 * VF2 works with both directed and undirected graphs. Only simple graphs are supported.
 * Self-loops or multi-edges must not be present in the graphs. Currently, the VF2
 * functions do not check that the input graph is simple: it is the responsibility
 * of the user to pass in valid input.
 * </para>
 */

/**
 * \function igraph_isomorphic_function_vf2
 * The generic VF2 interface
 *
 * </para><para>
 * This function is an implementation of the VF2 isomorphism algorithm,
 * see P. Foggia, C. Sansone, M. Vento, An Improved algorithm for
 * matching large graphs, Proc. of the 3rd IAPR-TC-15 International
 * Workshop on Graph-based Representations, Italy, 2001.</para>
 *
 * <para>For using it you need to define a callback function of type
 * \ref igraph_isohandler_t. This function will be called whenever VF2
 * finds an isomorphism between the two graphs. The mapping between
 * the two graphs will be also provided to this function. If the
 * callback returns a nonzero value then the search is continued,
 * otherwise it stops. The callback function must not destroy the
 * mapping vectors that are passed to it.
 * \param graph1 The first input graph.
 * \param graph2 The second input graph.
 * \param vertex_color1 An optional color vector for the first graph. If
 *   color vectors are given for both graphs, then the isomorphism is
 *   calculated on the colored graphs; i.e. two vertices can match
 *   only if their color also matches. Supply a null pointer here if
 *   your graphs are not colored.
 * \param vertex_color2 An optional color vector for the second graph. See
 *   the previous argument for explanation.
 * \param edge_color1 An optional edge color vector for the first
 *   graph. The matching edges in the two graphs must have matching
 *   colors as well. Supply a null pointer here if your graphs are not
 *   edge-colored.
 * \param edge_color2 The edge color vector for the second graph.
 * \param map12 Pointer to an initialized vector or \c NULL. If not \c
 *   NULL and the supplied graphs are isomorphic then the permutation
 *   taking \p graph1 to \p graph is stored here. If not \c NULL and the
 *   graphs are not isomorphic then a zero-length vector is returned.
 * \param map21 This is the same as \p map12, but for the permutation
 *   taking \p graph2 to \p graph1.
 * \param isohandler_fn The callback function to be called if an
 *   isomorphism is found. See also \ref igraph_isohandler_t.
 * \param node_compat_fn A pointer to a function of type \ref
 *   igraph_isocompat_t. This function will be called by the algorithm to
 *   determine whether two nodes are compatible.
 * \param edge_compat_fn A pointer to a function of type \ref
 *   igraph_isocompat_t. This function will be called by the algorithm to
 *   determine whether two edges are compatible.
 * \param arg Extra argument to supply to functions \p isohandler_fn, \p
 *   node_compat_fn and \p edge_compat_fn.
 * \return Error code.
 *
 * Time complexity: exponential.
 */

int igraph_isomorphic_function_vf2(const igraph_t *graph1, const igraph_t *graph2,
                                   const igraph_vector_int_t *vertex_color1,
                                   const igraph_vector_int_t *vertex_color2,
                                   const igraph_vector_int_t *edge_color1,
                                   const igraph_vector_int_t *edge_color2,
                                   igraph_vector_t *map12,
                                   igraph_vector_t *map21,
                                   igraph_isohandler_t *isohandler_fn,
                                   igraph_isocompat_t *node_compat_fn,
                                   igraph_isocompat_t *edge_compat_fn,
                                   void *arg) {

    long int no_of_nodes = igraph_vcount(graph1);
    long int no_of_edges = igraph_ecount(graph1);
    igraph_vector_t mycore_1, mycore_2, *core_1 = &mycore_1, *core_2 = &mycore_2;
    igraph_vector_t in_1, in_2, out_1, out_2;
    long int in_1_size = 0, in_2_size = 0, out_1_size = 0, out_2_size = 0;
    igraph_vector_int_t *inneis_1, *inneis_2, *outneis_1, *outneis_2;
    long int matched_nodes = 0;
    long int depth;
    long int cand1, cand2;
    long int last1, last2;
    igraph_stack_t path;
    igraph_lazy_adjlist_t inadj1, inadj2, outadj1, outadj2;
    igraph_vector_t indeg1, indeg2, outdeg1, outdeg2;
    long int vsize;

    if (igraph_is_directed(graph1) != igraph_is_directed(graph2)) {
        IGRAPH_ERROR("Cannot compare directed and undirected graphs",
                     IGRAPH_EINVAL);
    }

    if ( (vertex_color1 && !vertex_color2) || (!vertex_color1 && vertex_color2) ) {
        IGRAPH_WARNING("Only one graph is vertex-colored, vertex colors will be ignored");
        vertex_color1 = vertex_color2 = 0;
    }

    if ( (edge_color1 && !edge_color2) || (!edge_color1 && edge_color2)) {
        IGRAPH_WARNING("Only one graph is edge-colored, edge colors will be ignored");
        edge_color1 = edge_color2 = 0;
    }

    if (no_of_nodes != igraph_vcount(graph2) ||
        no_of_edges != igraph_ecount(graph2)) {
        return 0;
    }

    if (vertex_color1) {
        if (igraph_vector_int_size(vertex_color1) != no_of_nodes ||
            igraph_vector_int_size(vertex_color2) != no_of_nodes) {
            IGRAPH_ERROR("Invalid vertex color vector length", IGRAPH_EINVAL);
        }
    }

    if (edge_color1) {
        if (igraph_vector_int_size(edge_color1) != no_of_edges ||
            igraph_vector_int_size(edge_color2) != no_of_edges) {
            IGRAPH_ERROR("Invalid edge color vector length", IGRAPH_EINVAL);
        }
    }

    /* Check color distribution */
    if (vertex_color1) {
        int ret = 0;
        igraph_vector_int_t tmp1, tmp2;
        IGRAPH_CHECK(igraph_vector_int_copy(&tmp1, vertex_color1));
        IGRAPH_FINALLY(igraph_vector_int_destroy, &tmp1);
        IGRAPH_CHECK(igraph_vector_int_copy(&tmp2, vertex_color2));
        IGRAPH_FINALLY(igraph_vector_int_destroy, &tmp2);
        igraph_vector_int_sort(&tmp1);
        igraph_vector_int_sort(&tmp2);
        ret = !igraph_vector_int_all_e(&tmp1, &tmp2);
        igraph_vector_int_destroy(&tmp1);
        igraph_vector_int_destroy(&tmp2);
        IGRAPH_FINALLY_CLEAN(2);
        if (ret) {
            return 0;
        }
    }

    /* Check edge color distribution */
    if (edge_color1) {
        int ret = 0;
        igraph_vector_int_t tmp1, tmp2;
        IGRAPH_CHECK(igraph_vector_int_copy(&tmp1, edge_color1));
        IGRAPH_FINALLY(igraph_vector_int_destroy, &tmp1);
        IGRAPH_CHECK(igraph_vector_int_copy(&tmp2, edge_color2));
        IGRAPH_FINALLY(igraph_vector_int_destroy, &tmp2);
        igraph_vector_int_sort(&tmp1);
        igraph_vector_int_sort(&tmp2);
        ret = !igraph_vector_int_all_e(&tmp1, &tmp2);
        igraph_vector_int_destroy(&tmp1);
        igraph_vector_int_destroy(&tmp2);
        IGRAPH_FINALLY_CLEAN(2);
        if (ret) {
            return 0;
        }
    }

    if (map12) {
        core_1 = map12;
        IGRAPH_CHECK(igraph_vector_resize(core_1, no_of_nodes));
    } else {
        IGRAPH_VECTOR_INIT_FINALLY(core_1, no_of_nodes);
    }
    igraph_vector_fill(core_1, -1);
    if (map21) {
        core_2 = map21;
        IGRAPH_CHECK(igraph_vector_resize(core_2, no_of_nodes));
        igraph_vector_null(core_2);
    } else {
        IGRAPH_VECTOR_INIT_FINALLY(core_2, no_of_nodes);
    }
    igraph_vector_fill(core_2, -1);

    IGRAPH_VECTOR_INIT_FINALLY(&in_1, no_of_nodes);
    IGRAPH_VECTOR_INIT_FINALLY(&in_2, no_of_nodes);
    IGRAPH_VECTOR_INIT_FINALLY(&out_1, no_of_nodes);
    IGRAPH_VECTOR_INIT_FINALLY(&out_2, no_of_nodes);
    IGRAPH_CHECK(igraph_stack_init(&path, 0));
    IGRAPH_FINALLY(igraph_stack_destroy, &path);
    IGRAPH_CHECK(igraph_lazy_adjlist_init(graph1, &inadj1, IGRAPH_IN, IGRAPH_NO_LOOPS, IGRAPH_NO_MULTIPLE));
    IGRAPH_FINALLY(igraph_lazy_adjlist_destroy, &inadj1);
    IGRAPH_CHECK(igraph_lazy_adjlist_init(graph1, &outadj1, IGRAPH_OUT, IGRAPH_NO_LOOPS, IGRAPH_NO_MULTIPLE));
    IGRAPH_FINALLY(igraph_lazy_adjlist_destroy, &outadj1);
    IGRAPH_CHECK(igraph_lazy_adjlist_init(graph2, &inadj2, IGRAPH_IN, IGRAPH_NO_LOOPS, IGRAPH_NO_MULTIPLE));
    IGRAPH_FINALLY(igraph_lazy_adjlist_destroy, &inadj2);
    IGRAPH_CHECK(igraph_lazy_adjlist_init(graph2, &outadj2, IGRAPH_OUT, IGRAPH_NO_LOOPS, IGRAPH_NO_MULTIPLE));
    IGRAPH_FINALLY(igraph_lazy_adjlist_destroy, &outadj2);
    IGRAPH_VECTOR_INIT_FINALLY(&indeg1, 0);
    IGRAPH_VECTOR_INIT_FINALLY(&indeg2, 0);
    IGRAPH_VECTOR_INIT_FINALLY(&outdeg1, 0);
    IGRAPH_VECTOR_INIT_FINALLY(&outdeg2, 0);

    IGRAPH_CHECK(igraph_stack_reserve(&path, no_of_nodes * 2));
    IGRAPH_CHECK(igraph_degree(graph1, &indeg1, igraph_vss_all(),
                               IGRAPH_IN, IGRAPH_LOOPS));
    IGRAPH_CHECK(igraph_degree(graph2, &indeg2, igraph_vss_all(),
                               IGRAPH_IN, IGRAPH_LOOPS));
    IGRAPH_CHECK(igraph_degree(graph1, &outdeg1, igraph_vss_all(),
                               IGRAPH_OUT, IGRAPH_LOOPS));
    IGRAPH_CHECK(igraph_degree(graph2, &outdeg2, igraph_vss_all(),
                               IGRAPH_OUT, IGRAPH_LOOPS));

    depth = 0; last1 = -1; last2 = -1;
    while (depth >= 0) {
        long int i;

        IGRAPH_ALLOW_INTERRUPTION();

        cand1 = -1; cand2 = -1;
        /* Search for the next pair to try */
        if ((in_1_size != in_2_size) ||
            (out_1_size != out_2_size)) {
            /* step back, nothing to do */
        } else if (out_1_size > 0 && out_2_size > 0) {
            /**************************************************************/
            /* cand2, search not always needed */
            if (last2 >= 0) {
                cand2 = last2;
            } else {
                i = 0;
                while (cand2 < 0 && i < no_of_nodes) {
                    if (VECTOR(out_2)[i] > 0 && VECTOR(*core_2)[i] < 0) {
                        cand2 = i;
                    }
                    i++;
                }
            }
            /* search for cand1 now, it should be bigger than last1 */
            i = last1 + 1;
            while (cand1 < 0 && i < no_of_nodes) {
                if (VECTOR(out_1)[i] > 0 && VECTOR(*core_1)[i] < 0) {
                    cand1 = i;
                }
                i++;
            }
        } else if (in_1_size > 0 && in_2_size > 0) {
            /**************************************************************/
            /* cand2, search not always needed */
            if (last2 >= 0) {
                cand2 = last2;
            } else {
                i = 0;
                while (cand2 < 0 && i < no_of_nodes) {
                    if (VECTOR(in_2)[i] > 0 && VECTOR(*core_2)[i] < 0) {
                        cand2 = i;
                    }
                    i++;
                }
            }
            /* search for cand1 now, should be bigger than last1 */
            i = last1 + 1;
            while (cand1 < 0 && i < no_of_nodes) {
                if (VECTOR(in_1)[i] > 0 && VECTOR(*core_1)[i] < 0) {
                    cand1 = i;
                }
                i++;
            }
        } else {
            /**************************************************************/
            /* cand2, search not always needed */
            if (last2 >= 0) {
                cand2 = last2;
            } else {
                i = 0;
                while (cand2 < 0 && i < no_of_nodes) {
                    if (VECTOR(*core_2)[i] < 0) {
                        cand2 = i;
                    }
                    i++;
                }
            }
            /* search for cand1, should be bigger than last1 */
            i = last1 + 1;
            while (cand1 < 0 && i < no_of_nodes) {
                if (VECTOR(*core_1)[i] < 0) {
                    cand1 = i;
                }
                i++;
            }
        }

        /* Ok, we have cand1, cand2 as candidates. Or not? */
        if (cand1 < 0 || cand2 < 0) {
            /**************************************************************/
            /* dead end, step back, if possible. Otherwise we'll terminate */
            if (depth >= 1) {
                last2 = (long int) igraph_stack_pop(&path);
                last1 = (long int) igraph_stack_pop(&path);
                matched_nodes -= 1;
                VECTOR(*core_1)[last1] = -1;
                VECTOR(*core_2)[last2] = -1;

                if (VECTOR(in_1)[last1] != 0) {
                    in_1_size += 1;
                }
                if (VECTOR(out_1)[last1] != 0) {
                    out_1_size += 1;
                }
                if (VECTOR(in_2)[last2] != 0) {
                    in_2_size += 1;
                }
                if (VECTOR(out_2)[last2] != 0) {
                    out_2_size += 1;
                }

                inneis_1 = igraph_lazy_adjlist_get(&inadj1, (igraph_integer_t) last1);
                vsize = igraph_vector_int_size(inneis_1);
                for (i = 0; i < vsize; i++) {
                    long int node = (long int) VECTOR(*inneis_1)[i];
                    if (VECTOR(in_1)[node] == depth) {
                        VECTOR(in_1)[node] = 0;
                        in_1_size -= 1;
                    }
                }
                outneis_1 = igraph_lazy_adjlist_get(&outadj1, (igraph_integer_t) last1);
                vsize = igraph_vector_int_size(outneis_1);
                for (i = 0; i < vsize; i++) {
                    long int node = (long int) VECTOR(*outneis_1)[i];
                    if (VECTOR(out_1)[node] == depth) {
                        VECTOR(out_1)[node] = 0;
                        out_1_size -= 1;
                    }
                }
                inneis_2 = igraph_lazy_adjlist_get(&inadj2, (igraph_integer_t) last2);
                vsize = igraph_vector_int_size(inneis_2);
                for (i = 0; i < vsize; i++) {
                    long int node = (long int) VECTOR(*inneis_2)[i];
                    if (VECTOR(in_2)[node] == depth) {
                        VECTOR(in_2)[node] = 0;
                        in_2_size -= 1;
                    }
                }
                outneis_2 = igraph_lazy_adjlist_get(&outadj2, (igraph_integer_t) last2);
                vsize = igraph_vector_int_size(outneis_2);
                for (i = 0; i < vsize; i++) {
                    long int node = (long int) VECTOR(*outneis_2)[i];
                    if (VECTOR(out_2)[node] == depth) {
                        VECTOR(out_2)[node] = 0;
                        out_2_size -= 1;
                    }
                }

            } /* end of stepping back */

            depth -= 1;

        } else {
            /**************************************************************/
            /* step forward if worth, check if worth first */
            long int xin1 = 0, xin2 = 0, xout1 = 0, xout2 = 0;
            igraph_bool_t end = 0;
            inneis_1 = igraph_lazy_adjlist_get(&inadj1, (igraph_integer_t) cand1);
            outneis_1 = igraph_lazy_adjlist_get(&outadj1, (igraph_integer_t) cand1);
            inneis_2 = igraph_lazy_adjlist_get(&inadj2, (igraph_integer_t) cand2);
            outneis_2 = igraph_lazy_adjlist_get(&outadj2, (igraph_integer_t) cand2);
            if (VECTOR(indeg1)[cand1] != VECTOR(indeg2)[cand2] ||
                VECTOR(outdeg1)[cand1] != VECTOR(outdeg2)[cand2]) {
                end = 1;
            }
            if (vertex_color1 && VECTOR(*vertex_color1)[cand1] != VECTOR(*vertex_color2)[cand2]) {
                end = 1;
            }
            if (node_compat_fn && !node_compat_fn(graph1, graph2,
                                                  (igraph_integer_t) cand1,
                                                  (igraph_integer_t) cand2, arg)) {
                end = 1;
            }

            vsize = igraph_vector_int_size(inneis_1);
            for (i = 0; !end && i < vsize; i++) {
                long int node = (long int) VECTOR(*inneis_1)[i];
                if (VECTOR(*core_1)[node] >= 0) {
                    long int node2 = (long int) VECTOR(*core_1)[node];
                    /* check if there is a node2->cand2 edge */
                    if (!igraph_vector_int_binsearch2(inneis_2, node2)) {
                        end = 1;
                    } else if (edge_color1 || edge_compat_fn) {
                        igraph_integer_t eid1, eid2;
                        igraph_get_eid(graph1, &eid1, (igraph_integer_t) node,
                                       (igraph_integer_t) cand1, /*directed=*/ 1,
                                       /*error=*/ 1);
                        igraph_get_eid(graph2, &eid2, (igraph_integer_t) node2,
                                       (igraph_integer_t) cand2, /*directed=*/ 1,
                                       /*error=*/ 1);
                        if (edge_color1 && VECTOR(*edge_color1)[(long int)eid1] !=
                            VECTOR(*edge_color2)[(long int)eid2]) {
                            end = 1;
                        }
                        if (edge_compat_fn && !edge_compat_fn(graph1, graph2,
                                                              eid1, eid2, arg)) {
                            end = 1;
                        }
                    }
                } else {
                    if (VECTOR(in_1)[node] != 0) {
                        xin1++;
                    }
                    if (VECTOR(out_1)[node] != 0) {
                        xout1++;
                    }
                }
            }
            vsize = igraph_vector_int_size(outneis_1);
            for (i = 0; !end && i < vsize; i++) {
                long int node = (long int) VECTOR(*outneis_1)[i];
                if (VECTOR(*core_1)[node] >= 0) {
                    long int node2 = (long int) VECTOR(*core_1)[node];
                    /* check if there is a cand2->node2 edge */
                    if (!igraph_vector_int_binsearch2(outneis_2, node2)) {
                        end = 1;
                    } else if (edge_color1 || edge_compat_fn) {
                        igraph_integer_t eid1, eid2;
                        igraph_get_eid(graph1, &eid1, (igraph_integer_t) cand1,
                                       (igraph_integer_t) node, /*directed=*/ 1,
                                       /*error=*/ 1);
                        igraph_get_eid(graph2, &eid2, (igraph_integer_t) cand2,
                                       (igraph_integer_t) node2, /*directed=*/ 1,
                                       /*error=*/ 1);
                        if (edge_color1 && VECTOR(*edge_color1)[(long int)eid1] !=
                            VECTOR(*edge_color2)[(long int)eid2]) {
                            end = 1;
                        }
                        if (edge_compat_fn && !edge_compat_fn(graph1, graph2,
                                                              eid1, eid2, arg)) {
                            end = 1;
                        }
                    }
                } else {
                    if (VECTOR(in_1)[node] != 0) {
                        xin1++;
                    }
                    if (VECTOR(out_1)[node] != 0) {
                        xout1++;
                    }
                }
            }
            vsize = igraph_vector_int_size(inneis_2);
            for (i = 0; !end && i < vsize; i++) {
                long int node = (long int) VECTOR(*inneis_2)[i];
                if (VECTOR(*core_2)[node] >= 0) {
                    long int node2 = (long int) VECTOR(*core_2)[node];
                    /* check if there is a node2->cand1 edge */
                    if (!igraph_vector_int_binsearch2(inneis_1, node2)) {
                        end = 1;
                    } else if (edge_color1 || edge_compat_fn) {
                        igraph_integer_t eid1, eid2;
                        igraph_get_eid(graph1, &eid1, (igraph_integer_t) node2,
                                       (igraph_integer_t) cand1, /*directed=*/ 1,
                                       /*error=*/ 1);
                        igraph_get_eid(graph2, &eid2, (igraph_integer_t) node,
                                       (igraph_integer_t) cand2, /*directed=*/ 1,
                                       /*error=*/ 1);
                        if (edge_color1 && VECTOR(*edge_color1)[(long int)eid1] !=
                            VECTOR(*edge_color2)[(long int)eid2]) {
                            end = 1;
                        }
                        if (edge_compat_fn && !edge_compat_fn(graph1, graph2,
                                                              eid1, eid2, arg)) {
                            end = 1;
                        }
                    }
                } else {
                    if (VECTOR(in_2)[node] != 0) {
                        xin2++;
                    }
                    if (VECTOR(out_2)[node] != 0) {
                        xout2++;
                    }
                }
            }
            vsize = igraph_vector_int_size(outneis_2);
            for (i = 0; !end && i < vsize; i++) {
                long int node = (long int) VECTOR(*outneis_2)[i];
                if (VECTOR(*core_2)[node] >= 0) {
                    long int node2 = (long int) VECTOR(*core_2)[node];
                    /* check if there is a cand1->node2 edge */
                    if (!igraph_vector_int_binsearch2(outneis_1, node2)) {
                        end = 1;
                    } else if (edge_color1 || edge_compat_fn) {
                        igraph_integer_t eid1, eid2;
                        igraph_get_eid(graph1, &eid1, (igraph_integer_t) cand1,
                                       (igraph_integer_t) node2, /*directed=*/ 1,
                                       /*error=*/ 1);
                        igraph_get_eid(graph2, &eid2, (igraph_integer_t) cand2,
                                       (igraph_integer_t) node, /*directed=*/ 1,
                                       /*error=*/ 1);
                        if (edge_color1 && VECTOR(*edge_color1)[(long int)eid1] !=
                            VECTOR(*edge_color2)[(long int)eid2]) {
                            end = 1;
                        }
                        if (edge_compat_fn && !edge_compat_fn(graph1, graph2,
                                                              eid1, eid2, arg)) {
                            end = 1;
                        }
                    }
                } else {
                    if (VECTOR(in_2)[node] != 0) {
                        xin2++;
                    }
                    if (VECTOR(out_2)[node] != 0) {
                        xout2++;
                    }
                }
            }

            if (!end && (xin1 == xin2 && xout1 == xout2)) {
                /* Ok, we add the (cand1, cand2) pair to the mapping */
                depth += 1;
                IGRAPH_CHECK(igraph_stack_push(&path, cand1));
                IGRAPH_CHECK(igraph_stack_push(&path, cand2));
                matched_nodes += 1;
                VECTOR(*core_1)[cand1] = cand2;
                VECTOR(*core_2)[cand2] = cand1;

                /* update in_*, out_* */
                if (VECTOR(in_1)[cand1] != 0) {
                    in_1_size -= 1;
                }
                if (VECTOR(out_1)[cand1] != 0) {
                    out_1_size -= 1;
                }
                if (VECTOR(in_2)[cand2] != 0) {
                    in_2_size -= 1;
                }
                if (VECTOR(out_2)[cand2] != 0) {
                    out_2_size -= 1;
                }

                inneis_1 = igraph_lazy_adjlist_get(&inadj1, (igraph_integer_t) cand1);
                vsize = igraph_vector_int_size(inneis_1);
                for (i = 0; i < vsize; i++) {
                    long int node = (long int) VECTOR(*inneis_1)[i];
                    if (VECTOR(in_1)[node] == 0 && VECTOR(*core_1)[node] < 0) {
                        VECTOR(in_1)[node] = depth;
                        in_1_size += 1;
                    }
                }
                outneis_1 = igraph_lazy_adjlist_get(&outadj1, (igraph_integer_t) cand1);
                vsize = igraph_vector_int_size(outneis_1);
                for (i = 0; i < vsize; i++) {
                    long int node = (long int) VECTOR(*outneis_1)[i];
                    if (VECTOR(out_1)[node] == 0 && VECTOR(*core_1)[node] < 0) {
                        VECTOR(out_1)[node] = depth;
                        out_1_size += 1;
                    }
                }
                inneis_2 = igraph_lazy_adjlist_get(&inadj2, (igraph_integer_t) cand2);
                vsize = igraph_vector_int_size(inneis_2);
                for (i = 0; i < vsize; i++) {
                    long int node = (long int) VECTOR(*inneis_2)[i];
                    if (VECTOR(in_2)[node] == 0 && VECTOR(*core_2)[node] < 0) {
                        VECTOR(in_2)[node] = depth;
                        in_2_size += 1;
                    }
                }
                outneis_2 = igraph_lazy_adjlist_get(&outadj2, (igraph_integer_t) cand2);
                vsize = igraph_vector_int_size(outneis_2);
                for (i = 0; i < vsize; i++) {
                    long int node = (long int) VECTOR(*outneis_2)[i];
                    if (VECTOR(out_2)[node] == 0 && VECTOR(*core_2)[node] < 0) {
                        VECTOR(out_2)[node] = depth;
                        out_2_size += 1;
                    }
                }
                last1 = -1; last2 = -1;       /* this the first time here */
            } else {
                last1 = cand1;
                last2 = cand2;
            }

        }

        if (matched_nodes == no_of_nodes && isohandler_fn) {
            if (!isohandler_fn(core_1, core_2, arg)) {
                break;
            }
        }
    }

    igraph_vector_destroy(&outdeg2);
    igraph_vector_destroy(&outdeg1);
    igraph_vector_destroy(&indeg2);
    igraph_vector_destroy(&indeg1);
    igraph_lazy_adjlist_destroy(&outadj2);
    igraph_lazy_adjlist_destroy(&inadj2);
    igraph_lazy_adjlist_destroy(&outadj1);
    igraph_lazy_adjlist_destroy(&inadj1);
    igraph_stack_destroy(&path);
    igraph_vector_destroy(&out_2);
    igraph_vector_destroy(&out_1);
    igraph_vector_destroy(&in_2);
    igraph_vector_destroy(&in_1);
    IGRAPH_FINALLY_CLEAN(13);
    if (!map21) {
        igraph_vector_destroy(core_2);
        IGRAPH_FINALLY_CLEAN(1);
    }
    if (!map12) {
        igraph_vector_destroy(core_1);
        IGRAPH_FINALLY_CLEAN(1);
    }

    return 0;
}

typedef struct {
    igraph_isocompat_t *node_compat_fn, *edge_compat_fn;
    void *arg, *carg;
} igraph_i_iso_cb_data_t;

static igraph_bool_t igraph_i_isocompat_node_cb(
        const igraph_t *graph1,
        const igraph_t *graph2,
        const igraph_integer_t g1_num,
        const igraph_integer_t g2_num,
        void *arg) {
    igraph_i_iso_cb_data_t *data = arg;
    return data->node_compat_fn(graph1, graph2, g1_num, g2_num, data->carg);
}

static igraph_bool_t igraph_i_isocompat_edge_cb(
        const igraph_t *graph1,
        const igraph_t *graph2,
        const igraph_integer_t g1_num,
        const igraph_integer_t g2_num,
        void *arg) {
    igraph_i_iso_cb_data_t *data = arg;
    return data->edge_compat_fn(graph1, graph2, g1_num, g2_num, data->carg);
}

static igraph_bool_t igraph_i_isomorphic_vf2(igraph_vector_t *map12,
                                             igraph_vector_t *map21,
                                             void *arg) {
    igraph_i_iso_cb_data_t *data = arg;
    igraph_bool_t *iso = data->arg;
    IGRAPH_UNUSED(map12); IGRAPH_UNUSED(map21);
    *iso = 1;
    return 0;         /* don't need to continue */
}

/**
 * \function igraph_isomorphic_vf2
 * \brief Isomorphism via VF2
 *
 * </para><para>
 * This function performs the VF2 algorithm via calling \ref
 * igraph_isomorphic_function_vf2().
 *
 * </para><para> Note that this function cannot be used for
 * deciding subgraph isomorphism, use \ref igraph_subisomorphic_vf2()
 * for that.
 * \param graph1 The first graph, may be directed or undirected.
 * \param graph2 The second graph. It must have the same directedness
 *    as \p graph1, otherwise an error is reported.
 * \param vertex_color1 An optional color vector for the first graph. If
 *   color vectors are given for both graphs, then the isomorphism is
 *   calculated on the colored graphs; i.e. two vertices can match
 *   only if their color also matches. Supply a null pointer here if
 *   your graphs are not colored.
 * \param vertex_color2 An optional color vector for the second graph. See
 *   the previous argument for explanation.
 * \param edge_color1 An optional edge color vector for the first
 *   graph. The matching edges in the two graphs must have matching
 *   colors as well. Supply a null pointer here if your graphs are not
 *   edge-colored.
 * \param edge_color2 The edge color vector for the second graph.
 * \param iso Pointer to a logical constant, the result of the
 *    algorithm will be placed here.
 * \param map12 Pointer to an initialized vector or a NULL pointer. If not
 *    a NULL pointer then the mapping from \p graph1 to \p graph2 is
 *    stored here. If the graphs are not isomorphic then the vector is
 *    cleared (i.e. has zero elements).
 * \param map21 Pointer to an initialized vector or a NULL pointer. If not
 *    a NULL pointer then the mapping from \p graph2 to \p graph1 is
 *    stored here. If the graphs are not isomorphic then the vector is
 *    cleared (i.e. has zero elements).
 * \param node_compat_fn A pointer to a function of type \ref
 *   igraph_isocompat_t. This function will be called by the algorithm to
 *   determine whether two nodes are compatible.
 * \param edge_compat_fn A pointer to a function of type \ref
 *   igraph_isocompat_t. This function will be called by the algorithm to
 *   determine whether two edges are compatible.
 * \param arg Extra argument to supply to functions \p node_compat_fn
 *   and \p edge_compat_fn.
 * \return Error code.
 *
 * \sa \ref igraph_subisomorphic_vf2(),
 * \ref igraph_count_isomorphisms_vf2(),
 * \ref igraph_get_isomorphisms_vf2(),
 *
 * Time complexity: exponential, what did you expect?
 *
 * \example examples/simple/igraph_isomorphic_vf2.c
 */

int igraph_isomorphic_vf2(const igraph_t *graph1, const igraph_t *graph2,
                          const igraph_vector_int_t *vertex_color1,
                          const igraph_vector_int_t *vertex_color2,
                          const igraph_vector_int_t *edge_color1,
                          const igraph_vector_int_t *edge_color2,
                          igraph_bool_t *iso, igraph_vector_t *map12,
                          igraph_vector_t *map21,
                          igraph_isocompat_t *node_compat_fn,
                          igraph_isocompat_t *edge_compat_fn,
                          void *arg) {

    igraph_i_iso_cb_data_t data = { node_compat_fn, edge_compat_fn, iso, arg };
    igraph_isocompat_t *ncb = node_compat_fn ? igraph_i_isocompat_node_cb : 0;
    igraph_isocompat_t *ecb = edge_compat_fn ? igraph_i_isocompat_edge_cb : 0;
    *iso = 0;
    IGRAPH_CHECK(igraph_isomorphic_function_vf2(graph1, graph2,
                 vertex_color1, vertex_color2,
                 edge_color1, edge_color2,
                 map12, map21,
                 (igraph_isohandler_t*)
                 igraph_i_isomorphic_vf2,
                 ncb, ecb, &data));
    if (! *iso) {
        if (map12) {
            igraph_vector_clear(map12);
        }
        if (map21) {
            igraph_vector_clear(map21);
        }
    }
    return 0;
}

static igraph_bool_t igraph_i_count_isomorphisms_vf2(
        const igraph_vector_t *map12,
        const igraph_vector_t *map21,
        void *arg) {
    igraph_i_iso_cb_data_t *data = arg;
    igraph_integer_t *count = data->arg;
    IGRAPH_UNUSED(map12); IGRAPH_UNUSED(map21);
    *count += 1;
    return 1;         /* always continue */
}

/**
 * \function igraph_count_isomorphisms_vf2
 * Number of isomorphisms via VF2
 *
 * This function counts the number of isomorphic mappings between two
 * graphs. It uses the generic \ref igraph_isomorphic_function_vf2()
 * function.
 * \param graph1 The first input graph, may be directed or undirected.
 * \param graph2 The second input graph, it must have the same
 *   directedness as \p graph1, or an error will be reported.
 * \param vertex_color1 An optional color vector for the first graph. If
 *   color vectors are given for both graphs, then the isomorphism is
 *   calculated on the colored graphs; i.e. two vertices can match
 *   only if their color also matches. Supply a null pointer here if
 *   your graphs are not colored.
 * \param vertex_color2 An optional color vector for the second graph. See
 *   the previous argument for explanation.
 * \param edge_color1 An optional edge color vector for the first
 *   graph. The matching edges in the two graphs must have matching
 *   colors as well. Supply a null pointer here if your graphs are not
 *   edge-colored.
 * \param edge_color2 The edge color vector for the second graph.
 * \param count Point to an integer, the result will be stored here.
 * \param node_compat_fn A pointer to a function of type \ref
 *   igraph_isocompat_t. This function will be called by the algorithm to
 *   determine whether two nodes are compatible.
 * \param edge_compat_fn A pointer to a function of type \ref
 *   igraph_isocompat_t. This function will be called by the algorithm to
 *   determine whether two edges are compatible.
 * \param arg Extra argument to supply to functions \p node_compat_fn and
 *   \p edge_compat_fn.
 * \return Error code.
 *
 * Time complexity: exponential.
 */

int igraph_count_isomorphisms_vf2(const igraph_t *graph1, const igraph_t *graph2,
                                  const igraph_vector_int_t *vertex_color1,
                                  const igraph_vector_int_t *vertex_color2,
                                  const igraph_vector_int_t *edge_color1,
                                  const igraph_vector_int_t *edge_color2,
                                  igraph_integer_t *count,
                                  igraph_isocompat_t *node_compat_fn,
                                  igraph_isocompat_t *edge_compat_fn,
                                  void *arg) {

    igraph_i_iso_cb_data_t data = { node_compat_fn, edge_compat_fn,
                                    count, arg
                                  };
    igraph_isocompat_t *ncb = node_compat_fn ? igraph_i_isocompat_node_cb : 0;
    igraph_isocompat_t *ecb = edge_compat_fn ? igraph_i_isocompat_edge_cb : 0;
    *count = 0;
    IGRAPH_CHECK(igraph_isomorphic_function_vf2(graph1, graph2,
                 vertex_color1, vertex_color2,
                 edge_color1, edge_color2,
                 0, 0,
                 (igraph_isohandler_t*)
                 igraph_i_count_isomorphisms_vf2,
                 ncb, ecb, &data));
    return 0;
}

static void igraph_i_get_isomorphisms_free(igraph_vector_ptr_t *data) {
    long int i, n = igraph_vector_ptr_size(data);
    for (i = 0; i < n; i++) {
        igraph_vector_t *vec = VECTOR(*data)[i];
        igraph_vector_destroy(vec);
        igraph_free(vec);
    }
}

static igraph_bool_t igraph_i_get_isomorphisms_vf2(
        const igraph_vector_t *map12,
        const igraph_vector_t *map21,
        void *arg) {

    igraph_i_iso_cb_data_t *data = arg;
    igraph_vector_ptr_t *ptrvector = data->arg;
    igraph_vector_t *newvector = IGRAPH_CALLOC(1, igraph_vector_t);
    IGRAPH_UNUSED(map12);
    if (!newvector) {
        igraph_error("Out of memory", IGRAPH_FILE_BASENAME, __LINE__, IGRAPH_ENOMEM);
        return 0;           /* stop right here */
    }
    IGRAPH_FINALLY(igraph_free, newvector);
    IGRAPH_CHECK(igraph_vector_copy(newvector, map21));
    IGRAPH_FINALLY(igraph_vector_destroy, newvector);
    IGRAPH_CHECK(igraph_vector_ptr_push_back(ptrvector, newvector));
    IGRAPH_FINALLY_CLEAN(2);

    return 1;         /* continue finding subisomorphisms */
}

/**
 * \function igraph_get_isomorphisms_vf2
 * \brief Collect all isomorphic mappings of two graphs.
 *
 * This function finds all the isomorphic mappings between two simple
 * graphs. It uses the \ref igraph_isomorphic_function_vf2()
 * function. Call the function with the same graph as \p graph1 and \p
 * graph2 to get automorphisms.
 * \param graph1 The first input graph, may be directed or undirected.
 * \param graph2 The second input graph, it must have the same
 *   directedness as \p graph1, or an error will be reported.
 * \param vertex_color1 An optional color vector for the first graph. If
 *   color vectors are given for both graphs, then the isomorphism is
 *   calculated on the colored graphs; i.e. two vertices can match
 *   only if their color also matches. Supply a null pointer here if
 *   your graphs are not colored.
 * \param vertex_color2 An optional color vector for the second graph. See
 *   the previous argument for explanation.
 * \param edge_color1 An optional edge color vector for the first
 *   graph. The matching edges in the two graphs must have matching
 *   colors as well. Supply a null pointer here if your graphs are not
 *   edge-colored.
 * \param edge_color2 The edge color vector for the second graph.
 * \param maps Pointer vector. On return it is empty if the input graphs
 *   are not isomorphic. Otherwise it contains pointers to
 *   \ref igraph_vector_t objects, each vector is an
 *   isomorphic mapping of \p graph2 to \p graph1. Please note that
 *   you need to 1) Destroy the vectors via \ref
 *   igraph_vector_destroy(), 2) free them via
 *   \ref igraph_free() and then 3) call \ref
 *   igraph_vector_ptr_destroy() on the pointer vector to deallocate all
 *   memory when \p maps is no longer needed.
 * \param node_compat_fn A pointer to a function of type \ref
 *   igraph_isocompat_t. This function will be called by the algorithm to
 *   determine whether two nodes are compatible.
 * \param edge_compat_fn A pointer to a function of type \ref
 *   igraph_isocompat_t. This function will be called by the algorithm to
 *   determine whether two edges are compatible.
 * \param arg Extra argument to supply to functions \p node_compat_fn
 *   and \p edge_compat_fn.
 * \return Error code.
 *
 * Time complexity: exponential.
 */

int igraph_get_isomorphisms_vf2(const igraph_t *graph1,
                                const igraph_t *graph2,
                                const igraph_vector_int_t *vertex_color1,
                                const igraph_vector_int_t *vertex_color2,
                                const igraph_vector_int_t *edge_color1,
                                const igraph_vector_int_t *edge_color2,
                                igraph_vector_ptr_t *maps,
                                igraph_isocompat_t *node_compat_fn,
                                igraph_isocompat_t *edge_compat_fn,
                                void *arg) {

    igraph_i_iso_cb_data_t data = { node_compat_fn, edge_compat_fn, maps, arg };
    igraph_isocompat_t *ncb = node_compat_fn ? igraph_i_isocompat_node_cb : NULL;
    igraph_isocompat_t *ecb = edge_compat_fn ? igraph_i_isocompat_edge_cb : NULL;

    igraph_vector_ptr_clear(maps);
    IGRAPH_FINALLY(igraph_i_get_isomorphisms_free, maps);
    IGRAPH_CHECK(igraph_isomorphic_function_vf2(graph1, graph2,
                 vertex_color1, vertex_color2,
                 edge_color1, edge_color2,
                 NULL, NULL,
                 (igraph_isohandler_t*)
                 igraph_i_get_isomorphisms_vf2,
                 ncb, ecb, &data));
    IGRAPH_FINALLY_CLEAN(1);
    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_subisomorphic_function_vf2
 * Generic VF2 function for subgraph isomorphism problems
 *
 * This function is the pair of \ref igraph_isomorphic_function_vf2(),
 * for subgraph isomorphism problems. It searches for subgraphs of \p
 * graph1 which are isomorphic to \p graph2. When it founds an
 * isomorphic mapping it calls the supplied callback \p isohandler_fn.
 * The mapping (and its inverse) and the additional \p arg argument
 * are supplied to the callback.
 * \param graph1 The first input graph, may be directed or
 *    undirected. This is supposed to be the larger graph.
 * \param graph2 The second input graph, it must have the same
 *    directedness as \p graph1. This is supposed to be the smaller
 *    graph.
 * \param vertex_color1 An optional color vector for the first graph. If
 *   color vectors are given for both graphs, then the subgraph isomorphism is
 *   calculated on the colored graphs; i.e. two vertices can match
 *   only if their color also matches. Supply a null pointer here if
 *   your graphs are not colored.
 * \param vertex_color2 An optional color vector for the second graph. See
 *   the previous argument for explanation.
 * \param edge_color1 An optional edge color vector for the first
 *   graph. The matching edges in the two graphs must have matching
 *   colors as well. Supply a null pointer here if your graphs are not
 *   edge-colored.
 * \param edge_color2 The edge color vector for the second graph.
 * \param map12 Pointer to a vector or \c NULL. If not \c NULL, then an
 *    isomorphic mapping from \p graph1 to \p graph2 is stored here.
 * \param map21 Pointer to a vector ot \c NULL. If not \c NULL, then
 *    an isomorphic mapping from \p graph2 to \p graph1 is stored
 *    here.
 * \param isohandler_fn A pointer to a function of type \ref
 *   igraph_isohandler_t. This will be called whenever a subgraph
 *   isomorphism is found. If the function returns with a non-zero value
 *   then the search is continued, otherwise it stops and the function
 *   returns.
 * \param node_compat_fn A pointer to a function of type \ref
 *   igraph_isocompat_t. This function will be called by the algorithm to
 *   determine whether two nodes are compatible.
 * \param edge_compat_fn A pointer to a function of type \ref
 *   igraph_isocompat_t. This function will be called by the algorithm to
 *   determine whether two edges are compatible.
 * \param arg Extra argument to supply to functions \p isohandler_fn, \p
 *   node_compat_fn and \p edge_compat_fn.
 * \return Error code.
 *
 * Time complexity: exponential.
 */

int igraph_subisomorphic_function_vf2(const igraph_t *graph1,
                                      const igraph_t *graph2,
                                      const igraph_vector_int_t *vertex_color1,
                                      const igraph_vector_int_t *vertex_color2,
                                      const igraph_vector_int_t *edge_color1,
                                      const igraph_vector_int_t *edge_color2,
                                      igraph_vector_t *map12,
                                      igraph_vector_t *map21,
                                      igraph_isohandler_t *isohandler_fn,
                                      igraph_isocompat_t *node_compat_fn,
                                      igraph_isocompat_t *edge_compat_fn,
                                      void *arg) {

    long int no_of_nodes1 = igraph_vcount(graph1),
             no_of_nodes2 = igraph_vcount(graph2);
    long int no_of_edges1 = igraph_ecount(graph1),
             no_of_edges2 = igraph_ecount(graph2);
    igraph_vector_t mycore_1, mycore_2, *core_1 = &mycore_1, *core_2 = &mycore_2;
    igraph_vector_t in_1, in_2, out_1, out_2;
    long int in_1_size = 0, in_2_size = 0, out_1_size = 0, out_2_size = 0;
    igraph_vector_int_t *inneis_1, *inneis_2, *outneis_1, *outneis_2;
    long int matched_nodes = 0;
    long int depth;
    long int cand1, cand2;
    long int last1, last2;
    igraph_stack_t path;
    igraph_lazy_adjlist_t inadj1, inadj2, outadj1, outadj2;
    igraph_vector_t indeg1, indeg2, outdeg1, outdeg2;
    long int vsize;

    if (igraph_is_directed(graph1) != igraph_is_directed(graph2)) {
        IGRAPH_ERROR("Cannot compare directed and undirected graphs",
                     IGRAPH_EINVAL);
    }

    if (no_of_nodes1 < no_of_nodes2 ||
        no_of_edges1 < no_of_edges2) {
        return 0;
    }

    if ( (vertex_color1 && !vertex_color2) || (!vertex_color1 && vertex_color2) ) {
        IGRAPH_WARNING("Only one graph is vertex colored, colors will be ignored");
        vertex_color1 = vertex_color2 = 0;
    }

    if ( (edge_color1 && !edge_color2) || (!edge_color1 && edge_color2) ) {
        IGRAPH_WARNING("Only one graph is edge colored, colors will be ignored");
        edge_color1 = edge_color2 = 0;
    }

    if (vertex_color1) {
        if (igraph_vector_int_size(vertex_color1) != no_of_nodes1 ||
            igraph_vector_int_size(vertex_color2) != no_of_nodes2) {
            IGRAPH_ERROR("Invalid vertex color vector length", IGRAPH_EINVAL);
        }
    }

    if (edge_color1) {
        if (igraph_vector_int_size(edge_color1) != no_of_edges1 ||
            igraph_vector_int_size(edge_color2) != no_of_edges2) {
            IGRAPH_ERROR("Invalid edge color vector length", IGRAPH_EINVAL);
        }
    }

    /* Check color distribution */
    if (vertex_color1) {
        /* TODO */
    }

    /* Check edge color distribution */
    if (edge_color1) {
        /* TODO */
    }

    if (map12) {
        core_1 = map12;
        IGRAPH_CHECK(igraph_vector_resize(core_1, no_of_nodes1));
    } else {
        IGRAPH_VECTOR_INIT_FINALLY(core_1, no_of_nodes1);
    }
    igraph_vector_fill(core_1, -1);
    if (map21) {
        core_2 = map21;
        IGRAPH_CHECK(igraph_vector_resize(core_2, no_of_nodes2));
    } else {
        IGRAPH_VECTOR_INIT_FINALLY(core_2, no_of_nodes2);
    }
    igraph_vector_fill(core_2, -1);
    IGRAPH_VECTOR_INIT_FINALLY(&in_1, no_of_nodes1);
    IGRAPH_VECTOR_INIT_FINALLY(&in_2, no_of_nodes2);
    IGRAPH_VECTOR_INIT_FINALLY(&out_1, no_of_nodes1);
    IGRAPH_VECTOR_INIT_FINALLY(&out_2, no_of_nodes2);
    IGRAPH_CHECK(igraph_stack_init(&path, 0));
    IGRAPH_FINALLY(igraph_stack_destroy, &path);
    IGRAPH_CHECK(igraph_lazy_adjlist_init(graph1, &inadj1, IGRAPH_IN, IGRAPH_NO_LOOPS, IGRAPH_NO_MULTIPLE));
    IGRAPH_FINALLY(igraph_lazy_adjlist_destroy, &inadj1);
    IGRAPH_CHECK(igraph_lazy_adjlist_init(graph1, &outadj1, IGRAPH_OUT, IGRAPH_NO_LOOPS, IGRAPH_NO_MULTIPLE));
    IGRAPH_FINALLY(igraph_lazy_adjlist_destroy, &outadj1);
    IGRAPH_CHECK(igraph_lazy_adjlist_init(graph2, &inadj2, IGRAPH_IN, IGRAPH_NO_LOOPS, IGRAPH_NO_MULTIPLE));
    IGRAPH_FINALLY(igraph_lazy_adjlist_destroy, &inadj2);
    IGRAPH_CHECK(igraph_lazy_adjlist_init(graph2, &outadj2, IGRAPH_OUT, IGRAPH_NO_LOOPS, IGRAPH_NO_MULTIPLE));
    IGRAPH_FINALLY(igraph_lazy_adjlist_destroy, &outadj2);
    IGRAPH_VECTOR_INIT_FINALLY(&indeg1, 0);
    IGRAPH_VECTOR_INIT_FINALLY(&indeg2, 0);
    IGRAPH_VECTOR_INIT_FINALLY(&outdeg1, 0);
    IGRAPH_VECTOR_INIT_FINALLY(&outdeg2, 0);

    IGRAPH_CHECK(igraph_stack_reserve(&path, no_of_nodes2 * 2));
    IGRAPH_CHECK(igraph_degree(graph1, &indeg1, igraph_vss_all(),
                               IGRAPH_IN, IGRAPH_LOOPS));
    IGRAPH_CHECK(igraph_degree(graph2, &indeg2, igraph_vss_all(),
                               IGRAPH_IN, IGRAPH_LOOPS));
    IGRAPH_CHECK(igraph_degree(graph1, &outdeg1, igraph_vss_all(),
                               IGRAPH_OUT, IGRAPH_LOOPS));
    IGRAPH_CHECK(igraph_degree(graph2, &outdeg2, igraph_vss_all(),
                               IGRAPH_OUT, IGRAPH_LOOPS));

    depth = 0; last1 = -1; last2 = -1;
    while (depth >= 0) {
        long int i;

        IGRAPH_ALLOW_INTERRUPTION();

        cand1 = -1; cand2 = -1;
        /* Search for the next pair to try */
        if ((in_1_size < in_2_size) ||
            (out_1_size < out_2_size)) {
            /* step back, nothing to do */
        } else if (out_1_size > 0 && out_2_size > 0) {
            /**************************************************************/
            /* cand2, search not always needed */
            if (last2 >= 0) {
                cand2 = last2;
            } else {
                i = 0;
                while (cand2 < 0 && i < no_of_nodes2) {
                    if (VECTOR(out_2)[i] > 0 && VECTOR(*core_2)[i] < 0) {
                        cand2 = i;
                    }
                    i++;
                }
            }
            /* search for cand1 now, it should be bigger than last1 */
            i = last1 + 1;
            while (cand1 < 0 && i < no_of_nodes1) {
                if (VECTOR(out_1)[i] > 0 && VECTOR(*core_1)[i] < 0) {
                    cand1 = i;
                }
                i++;
            }
        } else if (in_1_size > 0 && in_2_size > 0) {
            /**************************************************************/
            /* cand2, search not always needed */
            if (last2 >= 0) {
                cand2 = last2;
            } else {
                i = 0;
                while (cand2 < 0 && i < no_of_nodes2) {
                    if (VECTOR(in_2)[i] > 0 && VECTOR(*core_2)[i] < 0) {
                        cand2 = i;
                    }
                    i++;
                }
            }
            /* search for cand1 now, should be bigger than last1 */
            i = last1 + 1;
            while (cand1 < 0 && i < no_of_nodes1) {
                if (VECTOR(in_1)[i] > 0 && VECTOR(*core_1)[i] < 0) {
                    cand1 = i;
                }
                i++;
            }
        } else {
            /**************************************************************/
            /* cand2, search not always needed */
            if (last2 >= 0) {
                cand2 = last2;
            } else {
                i = 0;
                while (cand2 < 0 && i < no_of_nodes2) {
                    if (VECTOR(*core_2)[i] < 0) {
                        cand2 = i;
                    }
                    i++;
                }
            }
            /* search for cand1, should be bigger than last1 */
            i = last1 + 1;
            while (cand1 < 0 && i < no_of_nodes1) {
                if (VECTOR(*core_1)[i] < 0) {
                    cand1 = i;
                }
                i++;
            }
        }

        /* Ok, we have cand1, cand2 as candidates. Or not? */
        if (cand1 < 0 || cand2 < 0) {
            /**************************************************************/
            /* dead end, step back, if possible. Otherwise we'll terminate */
            if (depth >= 1) {
                last2 = (long int) igraph_stack_pop(&path);
                last1 = (long int) igraph_stack_pop(&path);
                matched_nodes -= 1;
                VECTOR(*core_1)[last1] = -1;
                VECTOR(*core_2)[last2] = -1;

                if (VECTOR(in_1)[last1] != 0) {
                    in_1_size += 1;
                }
                if (VECTOR(out_1)[last1] != 0) {
                    out_1_size += 1;
                }
                if (VECTOR(in_2)[last2] != 0) {
                    in_2_size += 1;
                }
                if (VECTOR(out_2)[last2] != 0) {
                    out_2_size += 1;
                }

                inneis_1 = igraph_lazy_adjlist_get(&inadj1, (igraph_integer_t) last1);
                vsize = igraph_vector_int_size(inneis_1);
                for (i = 0; i < vsize; i++) {
                    long int node = (long int) VECTOR(*inneis_1)[i];
                    if (VECTOR(in_1)[node] == depth) {
                        VECTOR(in_1)[node] = 0;
                        in_1_size -= 1;
                    }
                }
                outneis_1 = igraph_lazy_adjlist_get(&outadj1, (igraph_integer_t) last1);
                vsize = igraph_vector_int_size(outneis_1);
                for (i = 0; i < vsize; i++) {
                    long int node = (long int) VECTOR(*outneis_1)[i];
                    if (VECTOR(out_1)[node] == depth) {
                        VECTOR(out_1)[node] = 0;
                        out_1_size -= 1;
                    }
                }
                inneis_2 = igraph_lazy_adjlist_get(&inadj2, (igraph_integer_t) last2);
                vsize = igraph_vector_int_size(inneis_2);
                for (i = 0; i < vsize; i++) {
                    long int node = (long int) VECTOR(*inneis_2)[i];
                    if (VECTOR(in_2)[node] == depth) {
                        VECTOR(in_2)[node] = 0;
                        in_2_size -= 1;
                    }
                }
                outneis_2 = igraph_lazy_adjlist_get(&outadj2, (igraph_integer_t) last2);
                vsize = igraph_vector_int_size(outneis_2);
                for (i = 0; i < vsize; i++) {
                    long int node = (long int) VECTOR(*outneis_2)[i];
                    if (VECTOR(out_2)[node] == depth) {
                        VECTOR(out_2)[node] = 0;
                        out_2_size -= 1;
                    }
                }

            } /* end of stepping back */

            depth -= 1;

        } else {
            /**************************************************************/
            /* step forward if worth, check if worth first */
            long int xin1 = 0, xin2 = 0, xout1 = 0, xout2 = 0;
            igraph_bool_t end = 0;
            inneis_1 = igraph_lazy_adjlist_get(&inadj1, (igraph_integer_t) cand1);
            outneis_1 = igraph_lazy_adjlist_get(&outadj1, (igraph_integer_t) cand1);
            inneis_2 = igraph_lazy_adjlist_get(&inadj2, (igraph_integer_t) cand2);
            outneis_2 = igraph_lazy_adjlist_get(&outadj2, (igraph_integer_t) cand2);
            if (VECTOR(indeg1)[cand1] < VECTOR(indeg2)[cand2] ||
                VECTOR(outdeg1)[cand1] < VECTOR(outdeg2)[cand2]) {
                end = 1;
            }
            if (vertex_color1 && VECTOR(*vertex_color1)[cand1] != VECTOR(*vertex_color2)[cand2]) {
                end = 1;
            }
            if (node_compat_fn && !node_compat_fn(graph1, graph2,
                                                  (igraph_integer_t) cand1,
                                                  (igraph_integer_t) cand2, arg)) {
                end = 1;
            }

            vsize = igraph_vector_int_size(inneis_1);
            for (i = 0; !end && i < vsize; i++) {
                long int node = (long int) VECTOR(*inneis_1)[i];
                if (VECTOR(*core_1)[node] < 0) {
                    if (VECTOR(in_1)[node] != 0) {
                        xin1++;
                    }
                    if (VECTOR(out_1)[node] != 0) {
                        xout1++;
                    }
                }
            }
            vsize = igraph_vector_int_size(outneis_1);
            for (i = 0; !end && i < vsize; i++) {
                long int node = (long int) VECTOR(*outneis_1)[i];
                if (VECTOR(*core_1)[node] < 0) {
                    if (VECTOR(in_1)[node] != 0) {
                        xin1++;
                    }
                    if (VECTOR(out_1)[node] != 0) {
                        xout1++;
                    }
                }
            }
            vsize = igraph_vector_int_size(inneis_2);
            for (i = 0; !end && i < vsize; i++) {
                long int node = (long int) VECTOR(*inneis_2)[i];
                if (VECTOR(*core_2)[node] >= 0) {
                    long int node2 = (long int) VECTOR(*core_2)[node];
                    /* check if there is a node2->cand1 edge */
                    if (!igraph_vector_int_binsearch2(inneis_1, node2)) {
                        end = 1;
                    } else if (edge_color1 || edge_compat_fn) {
                        igraph_integer_t eid1, eid2;
                        igraph_get_eid(graph1, &eid1, (igraph_integer_t) node2,
                                       (igraph_integer_t) cand1, /*directed=*/ 1,
                                       /*error=*/ 1);
                        igraph_get_eid(graph2, &eid2, (igraph_integer_t) node,
                                       (igraph_integer_t) cand2, /*directed=*/ 1,
                                       /*error=*/ 1);
                        if (edge_color1 && VECTOR(*edge_color1)[(long int)eid1] !=
                            VECTOR(*edge_color2)[(long int)eid2]) {
                            end = 1;
                        }
                        if (edge_compat_fn && !edge_compat_fn(graph1, graph2,
                                                              eid1, eid2, arg)) {
                            end = 1;
                        }
                    }
                } else {
                    if (VECTOR(in_2)[node] != 0) {
                        xin2++;
                    }
                    if (VECTOR(out_2)[node] != 0) {
                        xout2++;
                    }
                }
            }
            vsize = igraph_vector_int_size(outneis_2);
            for (i = 0; !end && i < vsize; i++) {
                long int node = (long int) VECTOR(*outneis_2)[i];
                if (VECTOR(*core_2)[node] >= 0) {
                    long int node2 = (long int) VECTOR(*core_2)[node];
                    /* check if there is a cand1->node2 edge */
                    if (!igraph_vector_int_binsearch2(outneis_1, node2)) {
                        end = 1;
                    } else if (edge_color1 || edge_compat_fn) {
                        igraph_integer_t eid1, eid2;
                        igraph_get_eid(graph1, &eid1, (igraph_integer_t) cand1,
                                       (igraph_integer_t) node2, /*directed=*/ 1,
                                       /*error=*/ 1);
                        igraph_get_eid(graph2, &eid2, (igraph_integer_t) cand2,
                                       (igraph_integer_t) node, /*directed=*/ 1,
                                       /*error=*/ 1);
                        if (edge_color1 && VECTOR(*edge_color1)[(long int)eid1] !=
                            VECTOR(*edge_color2)[(long int)eid2]) {
                            end = 1;
                        }
                        if (edge_compat_fn && !edge_compat_fn(graph1, graph2,
                                                              eid1, eid2, arg)) {
                            end = 1;
                        }
                    }
                } else {
                    if (VECTOR(in_2)[node] != 0) {
                        xin2++;
                    }
                    if (VECTOR(out_2)[node] != 0) {
                        xout2++;
                    }
                }
            }

            if (!end && (xin1 >= xin2 && xout1 >= xout2)) {
                /* Ok, we add the (cand1, cand2) pair to the mapping */
                depth += 1;
                IGRAPH_CHECK(igraph_stack_push(&path, cand1));
                IGRAPH_CHECK(igraph_stack_push(&path, cand2));
                matched_nodes += 1;
                VECTOR(*core_1)[cand1] = cand2;
                VECTOR(*core_2)[cand2] = cand1;

                /* update in_*, out_* */
                if (VECTOR(in_1)[cand1] != 0) {
                    in_1_size -= 1;
                }
                if (VECTOR(out_1)[cand1] != 0) {
                    out_1_size -= 1;
                }
                if (VECTOR(in_2)[cand2] != 0) {
                    in_2_size -= 1;
                }
                if (VECTOR(out_2)[cand2] != 0) {
                    out_2_size -= 1;
                }

                inneis_1 = igraph_lazy_adjlist_get(&inadj1, (igraph_integer_t) cand1);
                vsize = igraph_vector_int_size(inneis_1);
                for (i = 0; i < vsize; i++) {
                    long int node = (long int) VECTOR(*inneis_1)[i];
                    if (VECTOR(in_1)[node] == 0 && VECTOR(*core_1)[node] < 0) {
                        VECTOR(in_1)[node] = depth;
                        in_1_size += 1;
                    }
                }
                outneis_1 = igraph_lazy_adjlist_get(&outadj1, (igraph_integer_t) cand1);
                vsize = igraph_vector_int_size(outneis_1);
                for (i = 0; i < vsize; i++) {
                    long int node = (long int) VECTOR(*outneis_1)[i];
                    if (VECTOR(out_1)[node] == 0 && VECTOR(*core_1)[node] < 0) {
                        VECTOR(out_1)[node] = depth;
                        out_1_size += 1;
                    }
                }
                inneis_2 = igraph_lazy_adjlist_get(&inadj2, (igraph_integer_t) cand2);
                vsize = igraph_vector_int_size(inneis_2);
                for (i = 0; i < vsize; i++) {
                    long int node = (long int) VECTOR(*inneis_2)[i];
                    if (VECTOR(in_2)[node] == 0 && VECTOR(*core_2)[node] < 0) {
                        VECTOR(in_2)[node] = depth;
                        in_2_size += 1;
                    }
                }
                outneis_2 = igraph_lazy_adjlist_get(&outadj2, (igraph_integer_t) cand2);
                vsize = igraph_vector_int_size(outneis_2);
                for (i = 0; i < vsize; i++) {
                    long int node = (long int) VECTOR(*outneis_2)[i];
                    if (VECTOR(out_2)[node] == 0 && VECTOR(*core_2)[node] < 0) {
                        VECTOR(out_2)[node] = depth;
                        out_2_size += 1;
                    }
                }
                last1 = -1; last2 = -1;       /* this the first time here */
            } else {
                last1 = cand1;
                last2 = cand2;
            }

        }

        if (matched_nodes == no_of_nodes2 && isohandler_fn) {
            if (!isohandler_fn(core_1, core_2, arg)) {
                break;
            }
        }
    }

    igraph_vector_destroy(&outdeg2);
    igraph_vector_destroy(&outdeg1);
    igraph_vector_destroy(&indeg2);
    igraph_vector_destroy(&indeg1);
    igraph_lazy_adjlist_destroy(&outadj2);
    igraph_lazy_adjlist_destroy(&inadj2);
    igraph_lazy_adjlist_destroy(&outadj1);
    igraph_lazy_adjlist_destroy(&inadj1);
    igraph_stack_destroy(&path);
    igraph_vector_destroy(&out_2);
    igraph_vector_destroy(&out_1);
    igraph_vector_destroy(&in_2);
    igraph_vector_destroy(&in_1);
    IGRAPH_FINALLY_CLEAN(13);
    if (!map21) {
        igraph_vector_destroy(core_2);
        IGRAPH_FINALLY_CLEAN(1);
    }
    if (!map12) {
        igraph_vector_destroy(core_1);
        IGRAPH_FINALLY_CLEAN(1);
    }

    return 0;
}

static igraph_bool_t igraph_i_subisomorphic_vf2(
        const igraph_vector_t *map12,
        const igraph_vector_t *map21,
        void *arg) {
    igraph_i_iso_cb_data_t *data = arg;
    igraph_bool_t *iso = data->arg;
    IGRAPH_UNUSED(map12); IGRAPH_UNUSED(map21);
    *iso = 1;
    return 0; /* stop */
}

/**
 * \function igraph_subisomorphic_vf2
 * Decide subgraph isomorphism using VF2
 *
 * Decides whether a subgraph of \p graph1 is isomorphic to \p
 * graph2. It uses \ref igraph_subisomorphic_function_vf2().
 * \param graph1 The first input graph, may be directed or
 *    undirected. This is supposed to be the larger graph.
 * \param graph2 The second input graph, it must have the same
 *    directedness as \p graph1. This is supposed to be the smaller
 *    graph.
 * \param vertex_color1 An optional color vector for the first graph. If
 *   color vectors are given for both graphs, then the subgraph isomorphism is
 *   calculated on the colored graphs; i.e. two vertices can match
 *   only if their color also matches. Supply a null pointer here if
 *   your graphs are not colored.
 * \param vertex_color2 An optional color vector for the second graph. See
 *   the previous argument for explanation.
 * \param edge_color1 An optional edge color vector for the first
 *   graph. The matching edges in the two graphs must have matching
 *   colors as well. Supply a null pointer here if your graphs are not
 *   edge-colored.
 * \param edge_color2 The edge color vector for the second graph.
 * \param iso Pointer to a boolean. The result of the decision problem
 *    is stored here.
 * \param map12 Pointer to a vector or \c NULL. If not \c NULL, then an
 *    isomorphic mapping from \p graph1 to \p graph2 is stored here.
 * \param map21 Pointer to a vector ot \c NULL. If not \c NULL, then
 *    an isomorphic mapping from \p graph2 to \p graph1 is stored
 *    here.
 * \param node_compat_fn A pointer to a function of type \ref
 *   igraph_isocompat_t. This function will be called by the algorithm to
 *   determine whether two nodes are compatible.
 * \param edge_compat_fn A pointer to a function of type \ref
 *   igraph_isocompat_t. This function will be called by the algorithm to
 *   determine whether two edges are compatible.
 * \param arg Extra argument to supply to functions \p node_compat_fn
 *   and \p edge_compat_fn.
 * \return Error code.
 *
 * Time complexity: exponential.
 */

int igraph_subisomorphic_vf2(const igraph_t *graph1, const igraph_t *graph2,
                             const igraph_vector_int_t *vertex_color1,
                             const igraph_vector_int_t *vertex_color2,
                             const igraph_vector_int_t *edge_color1,
                             const igraph_vector_int_t *edge_color2,
                             igraph_bool_t *iso, igraph_vector_t *map12,
                             igraph_vector_t *map21,
                             igraph_isocompat_t *node_compat_fn,
                             igraph_isocompat_t *edge_compat_fn,
                             void *arg) {

    igraph_i_iso_cb_data_t data = { node_compat_fn, edge_compat_fn, iso, arg };
    igraph_isocompat_t *ncb = node_compat_fn ? igraph_i_isocompat_node_cb : 0;
    igraph_isocompat_t *ecb = edge_compat_fn ? igraph_i_isocompat_edge_cb : 0;

    *iso = 0;
    IGRAPH_CHECK(igraph_subisomorphic_function_vf2(graph1, graph2,
                 vertex_color1, vertex_color2,
                 edge_color1, edge_color2,
                 map12, map21,
                 (igraph_isohandler_t *)
                 igraph_i_subisomorphic_vf2,
                 ncb, ecb, &data));
    if (! *iso) {
        if (map12) {
            igraph_vector_clear(map12);
        }
        if (map21) {
            igraph_vector_clear(map21);
        }
    }
    return 0;
}

static igraph_bool_t igraph_i_count_subisomorphisms_vf2(
        const igraph_vector_t *map12,
        const igraph_vector_t *map21,
        void *arg) {
    igraph_i_iso_cb_data_t *data = arg;
    igraph_integer_t *count = data->arg;
    IGRAPH_UNUSED(map12); IGRAPH_UNUSED(map21);
    *count += 1;
    return 1;         /* always continue */
}

/**
 * \function igraph_count_subisomorphisms_vf2
 * Number of subgraph isomorphisms using VF2
 *
 * Count the number of isomorphisms between subgraphs of \p graph1 and
 * \p graph2. This function uses \ref
 * igraph_subisomorphic_function_vf2().
 * \param graph1 The first input graph, may be directed or
 *    undirected. This is supposed to be the larger graph.
 * \param graph2 The second input graph, it must have the same
 *    directedness as \p graph1. This is supposed to be the smaller
 *    graph.
 * \param vertex_color1 An optional color vector for the first graph. If
 *   color vectors are given for both graphs, then the subgraph isomorphism is
 *   calculated on the colored graphs; i.e. two vertices can match
 *   only if their color also matches. Supply a null pointer here if
 *   your graphs are not colored.
 * \param vertex_color2 An optional color vector for the second graph. See
 *   the previous argument for explanation.
 * \param edge_color1 An optional edge color vector for the first
 *   graph. The matching edges in the two graphs must have matching
 *   colors as well. Supply a null pointer here if your graphs are not
 *   edge-colored.
 * \param edge_color2 The edge color vector for the second graph.
 * \param count Pointer to an integer. The number of subgraph
 *    isomorphisms is stored here.
 * \param node_compat_fn A pointer to a function of type \ref
 *   igraph_isocompat_t. This function will be called by the algorithm to
 *   determine whether two nodes are compatible.
 * \param edge_compat_fn A pointer to a function of type \ref
 *   igraph_isocompat_t. This function will be called by the algorithm to
 *   determine whether two edges are compatible.
 * \param arg Extra argument to supply to functions \p node_compat_fn and
 *   \p edge_compat_fn.
 * \return Error code.
 *
 * Time complexity: exponential.
 */

int igraph_count_subisomorphisms_vf2(const igraph_t *graph1, const igraph_t *graph2,
                                     const igraph_vector_int_t *vertex_color1,
                                     const igraph_vector_int_t *vertex_color2,
                                     const igraph_vector_int_t *edge_color1,
                                     const igraph_vector_int_t *edge_color2,
                                     igraph_integer_t *count,
                                     igraph_isocompat_t *node_compat_fn,
                                     igraph_isocompat_t *edge_compat_fn,
                                     void *arg) {

    igraph_i_iso_cb_data_t data = { node_compat_fn, edge_compat_fn,
                                    count, arg
                                  };
    igraph_isocompat_t *ncb = node_compat_fn ? igraph_i_isocompat_node_cb : 0;
    igraph_isocompat_t *ecb = edge_compat_fn ? igraph_i_isocompat_edge_cb : 0;
    *count = 0;
    IGRAPH_CHECK(igraph_subisomorphic_function_vf2(graph1, graph2,
                 vertex_color1, vertex_color2,
                 edge_color1, edge_color2,
                 0, 0,
                 (igraph_isohandler_t*)
                 igraph_i_count_subisomorphisms_vf2,
                 ncb, ecb, &data));
    return 0;
}

static void igraph_i_get_subisomorphisms_free(igraph_vector_ptr_t *data) {
    long int i, n = igraph_vector_ptr_size(data);
    for (i = 0; i < n; i++) {
        igraph_vector_t *vec = VECTOR(*data)[i];
        igraph_vector_destroy(vec);
        igraph_free(vec);
    }
}

static igraph_bool_t igraph_i_get_subisomorphisms_vf2(
        const igraph_vector_t *map12,
        const igraph_vector_t *map21,
        void *arg) {

    igraph_i_iso_cb_data_t *data = arg;
    igraph_vector_ptr_t *vector = data->arg;
    igraph_vector_t *newvector = IGRAPH_CALLOC(1, igraph_vector_t);
    IGRAPH_UNUSED(map12);
    if (!newvector) {
        igraph_error("Out of memory", IGRAPH_FILE_BASENAME, __LINE__, IGRAPH_ENOMEM);
        return 0;           /* stop right here */
    }
    IGRAPH_FINALLY(igraph_free, newvector);
    IGRAPH_CHECK(igraph_vector_copy(newvector, map21));
    IGRAPH_FINALLY(igraph_vector_destroy, newvector);
    IGRAPH_CHECK(igraph_vector_ptr_push_back(vector, newvector));
    IGRAPH_FINALLY_CLEAN(2);

    return 1;         /* continue finding subisomorphisms */
}

/**
 * \function igraph_get_subisomorphisms_vf2
 * \brief Return all subgraph isomorphic mappings.
 *
 * This function collects all isomorphic mappings of \p graph2 to a
 * subgraph of \p graph1. It uses the \ref
 * igraph_subisomorphic_function_vf2() function. The graphs should be simple.
 * \param graph1 The first input graph, may be directed or
 *    undirected. This is supposed to be the larger graph.
 * \param graph2 The second input graph, it must have the same
 *    directedness as \p graph1. This is supposed to be the smaller
 *    graph.
 * \param vertex_color1 An optional color vector for the first graph. If
 *   color vectors are given for both graphs, then the subgraph isomorphism is
 *   calculated on the colored graphs; i.e. two vertices can match
 *   only if their color also matches. Supply a null pointer here if
 *   your graphs are not colored.
 * \param vertex_color2 An optional color vector for the second graph. See
 *   the previous argument for explanation.
 * \param edge_color1 An optional edge color vector for the first
 *   graph. The matching edges in the two graphs must have matching
 *   colors as well. Supply a null pointer here if your graphs are not
 *   edge-colored.
 * \param edge_color2 The edge color vector for the second graph.
 * \param maps Pointer vector. On return it contains pointers to
 *   \ref igraph_vector_t objects, each vector is an
 *   isomorphic mapping of \p graph2 to a subgraph of \p graph1. Please note that
 *   you need to 1) Destroy the vectors via \ref
 *   igraph_vector_destroy(), 2) free them via
 *   \ref igraph_free() and then 3) call \ref
 *   igraph_vector_ptr_destroy() on the pointer vector to deallocate all
 *   memory when \p maps is no longer needed.
 * \param node_compat_fn A pointer to a function of type \ref
 *   igraph_isocompat_t. This function will be called by the algorithm to
 *   determine whether two nodes are compatible.
 * \param edge_compat_fn A pointer to a function of type \ref
 *   igraph_isocompat_t. This function will be called by the algorithm to
 *   determine whether two edges are compatible.
 * \param arg Extra argument to supply to functions \p node_compat_fn
 *   and \p edge_compat_fn.
 * \return Error code.
 *
 * Time complexity: exponential.
 */

int igraph_get_subisomorphisms_vf2(const igraph_t *graph1,
                                   const igraph_t *graph2,
                                   const igraph_vector_int_t *vertex_color1,
                                   const igraph_vector_int_t *vertex_color2,
                                   const igraph_vector_int_t *edge_color1,
                                   const igraph_vector_int_t *edge_color2,
                                   igraph_vector_ptr_t *maps,
                                   igraph_isocompat_t *node_compat_fn,
                                   igraph_isocompat_t *edge_compat_fn,
                                   void *arg) {

    igraph_i_iso_cb_data_t data = { node_compat_fn, edge_compat_fn, maps, arg };
    igraph_isocompat_t *ncb = node_compat_fn ? igraph_i_isocompat_node_cb : NULL;
    igraph_isocompat_t *ecb = edge_compat_fn ? igraph_i_isocompat_edge_cb : NULL;

    igraph_vector_ptr_clear(maps);
    IGRAPH_FINALLY(igraph_i_get_subisomorphisms_free, maps);
    IGRAPH_CHECK(igraph_subisomorphic_function_vf2(graph1, graph2,
                 vertex_color1, vertex_color2,
                 edge_color1, edge_color2,
                 NULL, NULL,
                 (igraph_isohandler_t*)
                 igraph_i_get_subisomorphisms_vf2,
                 ncb, ecb, &data));
    IGRAPH_FINALLY_CLEAN(1);
    return IGRAPH_SUCCESS;
}
