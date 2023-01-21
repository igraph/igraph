/* -*- mode: C -*-  */
/* vim:set ts=4 sw=4 sts=4 et: */
/*
   IGraph library.
   Copyright (C) 2007-2020 The igraph development team

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

#include "igraph_community.h"

#include "igraph_adjlist.h"
#include "igraph_components.h"
#include "igraph_dqueue.h"
#include "igraph_interface.h"
#include "igraph_memory.h"
#include "igraph_nongraph.h"
#include "igraph_progress.h"
#include "igraph_stack.h"

#include "core/indheap.h"
#include "core/interruption.h"

#include <string.h>

static igraph_error_t igraph_i_rewrite_membership_vector(igraph_vector_int_t *membership) {
    const igraph_integer_t no = igraph_vector_int_max(membership) + 1;
    igraph_vector_t idx;
    igraph_integer_t realno = 0;
    const igraph_integer_t len = igraph_vector_int_size(membership);

    IGRAPH_VECTOR_INIT_FINALLY(&idx, no);
    for (igraph_integer_t i = 0; i < len; i++) {
        const igraph_integer_t t = VECTOR(*membership)[i];
        if (VECTOR(idx)[t]) {
            VECTOR(*membership)[i] = VECTOR(idx)[t] - 1;
        } else {
            VECTOR(idx)[t] = ++realno;
            VECTOR(*membership)[i] = VECTOR(idx)[t] - 1;
        }
    }
    igraph_vector_destroy(&idx);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_i_community_eb_get_merges2(const igraph_t *graph,
                                             const igraph_bool_t directed,
                                             const igraph_vector_int_t *edges,
                                             const igraph_vector_t *weights,
                                             igraph_matrix_int_t *res,
                                             igraph_vector_int_t *bridges,
                                             igraph_vector_t *modularity,
                                             igraph_vector_int_t *membership) {

    igraph_vector_int_t mymembership;
    const igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_real_t maxmod = -1;
    igraph_integer_t midx = 0;
    igraph_integer_t no_comps;
    const igraph_bool_t use_directed = directed && igraph_is_directed(graph);
    igraph_integer_t max_merges;

    if (membership) {
        IGRAPH_CHECK(igraph_vector_int_resize(membership, no_of_nodes));
    }
    if (modularity || res || bridges) {
        IGRAPH_CHECK(igraph_connected_components(graph, NULL, NULL, &no_comps, IGRAPH_WEAK));
        max_merges = no_of_nodes - no_comps;

        if (modularity) {
            IGRAPH_CHECK(igraph_vector_resize(modularity,
                                              max_merges + 1));
        }
        if (res) {
            IGRAPH_CHECK(igraph_matrix_int_resize(res, max_merges,
                                              2));
        }
        if (bridges) {
            IGRAPH_CHECK(igraph_vector_int_resize(bridges, max_merges));
        }
    }

    IGRAPH_CHECK(igraph_vector_int_init_range(&mymembership, 0, no_of_nodes));
    IGRAPH_FINALLY(igraph_vector_int_destroy, &mymembership);

    if (membership) {
        IGRAPH_CHECK(igraph_vector_int_update(membership, &mymembership));
    }

    IGRAPH_CHECK(igraph_modularity(graph, &mymembership, weights,
                                   /* resolution */ 1,
                                   use_directed, &maxmod));
    if (modularity) {
        VECTOR(*modularity)[0] = maxmod;
    }

    for (igraph_integer_t i = igraph_vector_int_size(edges) - 1; i >= 0; i--) {
        igraph_integer_t edge = VECTOR(*edges)[i];
        igraph_integer_t from = IGRAPH_FROM(graph, edge);
        igraph_integer_t to = IGRAPH_TO(graph, edge);
        igraph_integer_t c1 = VECTOR(mymembership)[from];
        igraph_integer_t c2 = VECTOR(mymembership)[to];
        igraph_real_t actmod;

        if (c1 != c2) {     /* this is a merge */
            if (res) {
                MATRIX(*res, midx, 0) = c1;
                MATRIX(*res, midx, 1) = c2;
            }
            if (bridges) {
                VECTOR(*bridges)[midx] = i;
            }

            /* The new cluster has id no_of_nodes+midx+1 */
            for (igraph_integer_t j = 0; j < no_of_nodes; j++) {
                if (VECTOR(mymembership)[j] == c1 ||
                    VECTOR(mymembership)[j] == c2) {
                    VECTOR(mymembership)[j] = no_of_nodes + midx;
                }
            }

            IGRAPH_CHECK(igraph_modularity(graph, &mymembership, weights,
                                           /* resolution */ 1,
                                           use_directed, &actmod));
            if (modularity) {
                VECTOR(*modularity)[midx + 1] = actmod;
                if (actmod > maxmod) {
                    maxmod = actmod;
                    if (membership) {
                        IGRAPH_CHECK(igraph_vector_int_update(membership, &mymembership));
                    }
                }
            }

            midx++;
        }
    }

    if (membership) {
        IGRAPH_CHECK(igraph_i_rewrite_membership_vector(membership));
    }

    igraph_vector_int_destroy(&mymembership);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}


/**
 * \function igraph_community_eb_get_merges
 * \brief Calculating the merges, i.e. the dendrogram for an edge betweenness community structure.
 *
 * This function is handy if you have a sequence of edges which are
 * gradually removed from the network and you would like to know how
 * the network falls apart into separate components. The edge sequence
 * may come from the \ref igraph_community_edge_betweenness()
 * function, but this is not necessary. Note that \ref
 * igraph_community_edge_betweenness() can also calculate the
 * dendrogram, via its \p merges argument. Merges happen when the
 * edge removal process is run backwards and two components become
 * connected.
 *
 * \param graph The input graph.
 * \param edges Vector containing the edges to be removed from the
 *    network, all edges are expected to appear exactly once in the
 *    vector.
 * \param directed Whether to use the directed or undirected version
 *    of modularity. Will be ignored for undirected graphs.
 * \param weights An optional vector containing edge weights. If null,
 *     the unweighted modularity scores will be calculated. If not null,
 *     the weighted modularity scores will be calculated. Ignored if both
 *     \p modularity and \p membership are \c NULL pointers.
 * \param res Pointer to an initialized matrix, if not \c NULL then the
 *    dendrogram will be stored here, in the same form as for the
 *    \ref igraph_community_walktrap() function: the matrix has two columns
 *    and each line is a merge given by the IDs of the merged
 *    components. The component IDs are numbered from zero and
 *    component IDs smaller than the number of vertices in the graph
 *    belong to individual vertices. The non-trivial components
 *    containing at least two vertices are numbered from \c n, where \c n is
 *    the number of vertices in the graph. So if the first line
 *    contains \c a and \c b that means that components \c a and \c b
 *    are merged into component \c n, the second line creates
 *    component <code>n+1</code>, etc. The matrix will be resized as needed.
 * \param bridges Pointer to an initialized vector of \c NULL. If not
 *     \c NULL then the indices into \p edges of all edges which caused
 *     one of the merges will be put here. This is equal to all edge removals
 *     which separated the network into more components, in reverse order.
 * \param modularity If not a null pointer, then the modularity values
 *    for the different divisions, corresponding to the merges matrix,
 *    will be stored here.
 * \param membership If not a null pointer, then the membership vector
 *    for the best division (in terms of modularity) will be stored
 *    here.
 * \return Error code.
 *
 * \sa \ref igraph_community_edge_betweenness().
 *
 * Time complexity: O(|E|+|V|log|V|), |V| is the number of vertices,
 * |E| is the number of edges.
 */
igraph_error_t igraph_community_eb_get_merges(const igraph_t *graph,
                                   const igraph_bool_t directed,
                                   const igraph_vector_int_t *edges,
                                   const igraph_vector_t *weights,
                                   igraph_matrix_int_t *res,
                                   igraph_vector_int_t *bridges,
                                   igraph_vector_t *modularity,
                                   igraph_vector_int_t *membership) {

    const igraph_integer_t no_of_nodes = igraph_vcount(graph);
    const igraph_integer_t no_of_edges = igraph_ecount(graph);
    igraph_vector_int_t ptr;
    igraph_integer_t midx = 0;
    igraph_integer_t no_comps;
    const igraph_integer_t no_removed_edges = igraph_vector_int_size(edges);
    igraph_integer_t max_merges;

    if (! igraph_vector_int_isininterval(edges, 0, no_of_edges-1)) {
        IGRAPH_ERROR("Invalid edge ID.", IGRAPH_EINVAL);
    }
    if (no_removed_edges < no_of_edges) {
            IGRAPH_ERRORF("Number of removed edges (%" IGRAPH_PRId ") should be equal to "
                    "number of edges in graph (%" IGRAPH_PRId ").", IGRAPH_EINVAL,
                    no_removed_edges, no_of_edges);
    }

    /* catch null graph early */
    if (no_of_nodes == 0) {
        if (res) {
            IGRAPH_CHECK(igraph_matrix_int_resize(res, 0, 2));
        }
        if (bridges) {
            igraph_vector_int_clear(bridges);
        }
        if (modularity) {
            IGRAPH_CHECK(igraph_vector_resize(modularity, 1));
            VECTOR(*modularity)[0] = IGRAPH_NAN;
        }
        if (membership) {
            igraph_vector_int_clear(membership);
        }
        return IGRAPH_SUCCESS;
    }

    if (membership || modularity) {
        return igraph_i_community_eb_get_merges2(graph,
                directed && igraph_is_directed(graph),
                edges, weights,
                res, bridges, modularity, membership);
    }

    IGRAPH_CHECK(igraph_connected_components(graph, NULL, NULL, &no_comps, IGRAPH_WEAK));

    max_merges = no_of_nodes - no_comps;
    IGRAPH_VECTOR_INT_INIT_FINALLY(&ptr, no_of_nodes * 2 - 1);
    if (res) {
        IGRAPH_CHECK(igraph_matrix_int_resize(res, max_merges, 2));
    }
    if (bridges) {
        IGRAPH_CHECK(igraph_vector_int_resize(bridges, max_merges));
    }

    for (igraph_integer_t i = igraph_vector_int_size(edges) - 1; i >= 0; i--) {
        igraph_integer_t edge = VECTOR(*edges)[i];
        igraph_integer_t from, to, c1, c2, idx;
        IGRAPH_CHECK(igraph_edge(graph, edge, &from, &to));
        idx = from + 1;
        while (VECTOR(ptr)[idx - 1] != 0) {
            idx = VECTOR(ptr)[idx - 1];
        }
        c1 = idx - 1;
        idx = to + 1;
        while (VECTOR(ptr)[idx - 1] != 0) {
            idx = VECTOR(ptr)[idx - 1];
        }
        c2 = idx - 1;
        if (c1 != c2) {     /* this is a merge */
            if (res) {
                MATRIX(*res, midx, 0) = c1;
                MATRIX(*res, midx, 1) = c2;
            }
            if (bridges) {
                VECTOR(*bridges)[midx] = i;
            }

            VECTOR(ptr)[c1] = no_of_nodes + midx + 1;
            VECTOR(ptr)[c2] = no_of_nodes + midx + 1;
            VECTOR(ptr)[from] = no_of_nodes + midx + 1;
            VECTOR(ptr)[to] = no_of_nodes + midx + 1;

            midx++;
        }
    }

    igraph_vector_int_destroy(&ptr);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}

/* Find the smallest active element in the vector */
static igraph_integer_t igraph_i_vector_which_max_not_null(const igraph_vector_t *v,
                                                   const bool *passive) {
    igraph_integer_t which, i = 0, size = igraph_vector_size(v);
    igraph_real_t max;
    while (passive[i]) {
        i++;
    }
    which = i;
    max = VECTOR(*v)[which];
    for (i++; i < size; i++) {
        igraph_real_t elem = VECTOR(*v)[i];
        if (!passive[i] && elem > max) {
            max = elem;
            which = i;
        }
    }

    return which;
}

/**
 * \function igraph_community_edge_betweenness
 * \brief Community finding based on edge betweenness.
 *
 * Community structure detection based on the betweenness of the edges
 * in the network. The algorithm was invented by M. Girvan and
 * M. Newman, see: M. Girvan and M. E. J. Newman: Community structure in
 * social and biological networks, Proc. Nat. Acad. Sci. USA 99, 7821-7826
 * (2002). https://doi.org/10.1073/pnas.122653799
 *
 * </para><para>
 * The idea is that the betweenness of the edges connecting two
 * communities is typically high, as many of the shortest paths
 * between nodes in separate communities go through them. So we
 * gradually remove the edge with highest betweenness from the
 * network, and recalculate edge betweenness after every removal.
 * This way sooner or later the network splits into two components,
 * then after a while one of these components splits again into two smaller
 * components, and so on until all edges are removed. This is a divisive
 * hierarchical approach, the result of which is a dendrogram.
 *
 * </para><para>
 * In directed graphs, when \p directed is set to true, the directed version
 * of betweenness and modularity are used, however, only splits into
 * \em weakly connected components are detected.
 *
 * \param graph The input graph.
 * \param result Pointer to an initialized vector, the result will be
 *     stored here, the IDs of the removed edges in the order of their
 *     removal. It will be resized as needed. It may be \c NULL if
 *     the edge IDs are not needed by the caller.
 * \param edge_betweenness Pointer to an initialized vector or
 *     \c NULL. In the former case the edge betweenness of the removed
 *     edge is stored here. The vector will be resized as needed.
 * \param merges Pointer to an initialized matrix or \c NULL. If not \c NULL
 *     then merges performed by the algorithm are stored here. Even if
 *     this is a divisive algorithm, we can replay it backwards and
 *     note which two clusters were merged. Clusters are numbered from
 *     zero, see the \p merges argument of \ref igraph_community_walktrap()
 *     for details. The matrix will be resized as needed.
 * \param bridges Pointer to an initialized vector of \c NULL. If not
 *     \c NULL then the indices into \p result of all edges which caused
 *     one of the \p merges will be put here. This is equivalent to all edge removals
 *     which separated the network into more components, in reverse order.
 * \param modularity If not a null pointer, then the modularity values
 *     of the different divisions are stored here, in the order
 *     corresponding to the merge matrix. The modularity values will
 *     take weights into account if \p weights is not null.
 * \param membership If not a null pointer, then the membership vector,
 *     corresponding to the highest modularity value, is stored here.
 * \param directed Logical constant. Controls whether to calculate directed
 *    betweenness (i.e. directed paths) for directed graphs, and whether
 *    to use the directed version of modularity. It is ignored for undirected
 *    graphs.
 * \param weights An optional vector containing edge weights. If null,
 *     the unweighted edge betweenness scores will be calculated and
 *     used. If not null, the weighted edge betweenness scores will be
 *     calculated and used.
 * \return Error code.
 *
 * \sa \ref igraph_community_eb_get_merges(), \ref
 * igraph_community_spinglass(), \ref igraph_community_walktrap().
 *
 * Time complexity: O(|V||E|^2), as the betweenness calculation requires
 * O(|V||E|) and we do it |E|-1 times.
 *
 * \example examples/simple/igraph_community_edge_betweenness.c
 */
igraph_error_t igraph_community_edge_betweenness(const igraph_t *graph,
                                      igraph_vector_int_t *result,
                                      igraph_vector_t *edge_betweenness,
                                      igraph_matrix_int_t *merges,
                                      igraph_vector_int_t *bridges,
                                      igraph_vector_t *modularity,
                                      igraph_vector_int_t *membership,
                                      igraph_bool_t directed,
                                      const igraph_vector_t *weights) {

    const igraph_integer_t no_of_nodes = igraph_vcount(graph);
    const igraph_integer_t no_of_edges = igraph_ecount(graph);
    double *distance, *tmpscore;
    double *nrgeo;

    igraph_inclist_t elist_out, elist_in, parents;
    igraph_inclist_t *elist_out_p, *elist_in_p;
    igraph_vector_int_t *neip;
    igraph_integer_t neino;
    igraph_vector_t eb;
    igraph_integer_t maxedge, pos;
    igraph_integer_t from, to;
    igraph_bool_t result_owned = false;
    igraph_stack_int_t stack;
    igraph_real_t steps, steps_done;

    bool *passive;

    /* Needed only for the unweighted case */
    igraph_dqueue_int_t q;

    /* Needed only for the weighted case */
    igraph_2wheap_t heap;

    if (result == NULL) {
        result = IGRAPH_CALLOC(1, igraph_vector_int_t);
        IGRAPH_CHECK_OOM(result, "Insufficient memory for edge betweenness-based community detection.");
        IGRAPH_FINALLY(igraph_free, result);

        IGRAPH_VECTOR_INT_INIT_FINALLY(result, 0);
        result_owned = true;
    }

    directed = directed && igraph_is_directed(graph);
    if (directed) {
        IGRAPH_CHECK(igraph_inclist_init(graph, &elist_out, IGRAPH_OUT, IGRAPH_LOOPS_ONCE));
        IGRAPH_FINALLY(igraph_inclist_destroy, &elist_out);
        IGRAPH_CHECK(igraph_inclist_init(graph, &elist_in, IGRAPH_IN, IGRAPH_LOOPS_ONCE));
        IGRAPH_FINALLY(igraph_inclist_destroy, &elist_in);
        elist_out_p = &elist_out;
        elist_in_p = &elist_in;
    } else {
        IGRAPH_CHECK(igraph_inclist_init(graph, &elist_out, IGRAPH_ALL, IGRAPH_LOOPS_TWICE));
        IGRAPH_FINALLY(igraph_inclist_destroy, &elist_out);
        elist_out_p = elist_in_p = &elist_out;
    }

    distance = IGRAPH_CALLOC(no_of_nodes, double);
    IGRAPH_CHECK_OOM(distance, "Insufficient memory for edge betweenness-based community detection.");
    IGRAPH_FINALLY(igraph_free, distance);

    nrgeo = IGRAPH_CALLOC(no_of_nodes, double);
    IGRAPH_CHECK_OOM(nrgeo, "Insufficient memory for edge betweenness-based community detection.");
    IGRAPH_FINALLY(igraph_free, nrgeo);

    tmpscore = IGRAPH_CALLOC(no_of_nodes, double);
    IGRAPH_CHECK_OOM(tmpscore, "Insufficient memory for edge betweenness-based community detection.");
    IGRAPH_FINALLY(igraph_free, tmpscore);

    if (weights == NULL) {
        IGRAPH_DQUEUE_INT_INIT_FINALLY(&q, 100);
    } else {
        if (igraph_vector_size(weights) != no_of_edges) {
            IGRAPH_ERROR("Weight vector length must agree with number of edges.", IGRAPH_EINVAL);
        }

        if (no_of_edges > 0) {
            /* Must not call vector_min on empty vector */
            igraph_real_t minweight = igraph_vector_min(weights);
            if (minweight <= 0) {
                IGRAPH_ERROR("Weights must be strictly positive.", IGRAPH_EINVAL);
            }

            if (isnan(minweight)) {
                IGRAPH_ERROR("Weights must not be NaN.", IGRAPH_EINVAL);
            }
        }

        if (membership != NULL) {
            IGRAPH_WARNING("Membership vector will be selected based on the highest "
                           "modularity score.");
        }

        if (modularity != NULL || membership != NULL) {
            IGRAPH_WARNING("Modularity calculation with weighted edge betweenness "
                           "community detection might not make sense -- modularity treats edge "
                           "weights as similarities while edge betwenness treats them as "
                           "distances.");
        }

        IGRAPH_CHECK(igraph_2wheap_init(&heap, no_of_nodes));
        IGRAPH_FINALLY(igraph_2wheap_destroy, &heap);
        IGRAPH_CHECK(igraph_inclist_init_empty(&parents, no_of_nodes));
        IGRAPH_FINALLY(igraph_inclist_destroy, &parents);
    }

    IGRAPH_STACK_INT_INIT_FINALLY(&stack, no_of_nodes);

    IGRAPH_CHECK(igraph_vector_int_resize(result, no_of_edges));
    if (edge_betweenness) {
        IGRAPH_CHECK(igraph_vector_resize(edge_betweenness, no_of_edges));
        if (no_of_edges > 0) {
            VECTOR(*edge_betweenness)[no_of_edges - 1] = 0;
        }
    }

    IGRAPH_VECTOR_INIT_FINALLY(&eb, no_of_edges);

    passive = IGRAPH_CALLOC(no_of_edges, bool);
    IGRAPH_CHECK_OOM(passive, "Insufficient memory for edge betweenness-based community detection.");
    IGRAPH_FINALLY(igraph_free, passive);

    /* Estimate the number of steps to be taken.
     * It is assumed that one iteration is O(|E||V|), but |V| is constant
     * anyway, so we will have approximately |E|^2 / 2 steps, and one
     * iteration of the outer loop advances the step counter by the number
     * of remaining edges at that iteration.
     */
    steps = no_of_edges / 2.0 * (no_of_edges + 1);
    steps_done = 0;

    for (igraph_integer_t e = 0; e < no_of_edges; steps_done += no_of_edges - e, e++) {
        IGRAPH_PROGRESS("Edge betweenness community detection: ",
                        100.0 * steps_done / steps, NULL);

        igraph_vector_null(&eb);

        if (weights == NULL) {
            /* Unweighted variant follows */

            /* The following for loop is copied almost intact from
             * igraph_edge_betweenness_cutoff */
            for (igraph_integer_t source = 0; source < no_of_nodes; source++) {

                IGRAPH_ALLOW_INTERRUPTION();

                memset(distance, 0, (size_t) no_of_nodes * sizeof(double));
                memset(nrgeo, 0, (size_t) no_of_nodes * sizeof(double));
                memset(tmpscore, 0, (size_t) no_of_nodes * sizeof(double));
                igraph_stack_int_clear(&stack); /* it should be empty anyway... */

                IGRAPH_CHECK(igraph_dqueue_int_push(&q, source));

                nrgeo[source] = 1;
                distance[source] = 0;

                while (!igraph_dqueue_int_empty(&q)) {
                    igraph_integer_t actnode = igraph_dqueue_int_pop(&q);

                    neip = igraph_inclist_get(elist_out_p, actnode);
                    neino = igraph_vector_int_size(neip);
                    for (igraph_integer_t i = 0; i < neino; i++) {
                        igraph_integer_t edge = VECTOR(*neip)[i];
                        igraph_integer_t neighbor = IGRAPH_OTHER(graph, edge, actnode);
                        if (nrgeo[neighbor] != 0) {
                            /* we've already seen this node, another shortest path? */
                            if (distance[neighbor] == distance[actnode] + 1) {
                                nrgeo[neighbor] += nrgeo[actnode];
                            }
                        } else {
                            /* we haven't seen this node yet */
                            nrgeo[neighbor] += nrgeo[actnode];
                            distance[neighbor] = distance[actnode] + 1;
                            IGRAPH_CHECK(igraph_dqueue_int_push(&q, neighbor));
                            IGRAPH_CHECK(igraph_stack_int_push(&stack, neighbor));
                        }
                    }
                } /* while !igraph_dqueue_int_empty */

                /* Ok, we've the distance of each node and also the number of
                   shortest paths to them. Now we do an inverse search, starting
                   with the farthest nodes. */
                while (!igraph_stack_int_empty(&stack)) {
                    igraph_integer_t actnode = igraph_stack_int_pop(&stack);
                    if (distance[actnode] < 1) {
                        continue;    /* skip source node */
                    }

                    /* set the temporary score of the friends */
                    neip = igraph_inclist_get(elist_in_p, actnode);
                    neino = igraph_vector_int_size(neip);
                    for (igraph_integer_t i = 0; i < neino; i++) {
                        igraph_integer_t edge = VECTOR(*neip)[i];
                        igraph_integer_t neighbor = IGRAPH_OTHER(graph, edge, actnode);
                        if (distance[neighbor] == distance[actnode] - 1 &&
                            nrgeo[neighbor] != 0) {
                            tmpscore[neighbor] +=
                                (tmpscore[actnode] + 1) * nrgeo[neighbor] / nrgeo[actnode];
                            VECTOR(eb)[edge] +=
                                (tmpscore[actnode] + 1) * nrgeo[neighbor] / nrgeo[actnode];
                        }
                    }
                }
                /* Ok, we've the scores for this source */
            } /* for source <= no_of_nodes */
        } else {
            /* Weighted variant follows */

            const igraph_real_t eps = IGRAPH_SHORTEST_PATH_EPSILON;
            int cmp_result;

            /* The following for loop is copied almost intact from
             * igraph_i_edge_betweenness_cutoff_weighted */
            for (igraph_integer_t source = 0; source < no_of_nodes; source++) {
                /* This will contain the edge betweenness in the current step */
                IGRAPH_ALLOW_INTERRUPTION();

                memset(distance, 0, (size_t) no_of_nodes * sizeof(double));
                memset(nrgeo, 0, (size_t) no_of_nodes * sizeof(double));
                memset(tmpscore, 0, (size_t) no_of_nodes * sizeof(double));

                IGRAPH_CHECK(igraph_2wheap_push_with_index(&heap, source, 0));
                distance[source] = 1.0;
                nrgeo[source] = 1;

                while (!igraph_2wheap_empty(&heap)) {
                    igraph_integer_t minnei = igraph_2wheap_max_index(&heap);
                    igraph_real_t mindist = -igraph_2wheap_delete_max(&heap);

                    IGRAPH_CHECK(igraph_stack_int_push(&stack, minnei));

                    neip = igraph_inclist_get(elist_out_p, minnei);
                    neino = igraph_vector_int_size(neip);

                    for (igraph_integer_t i = 0; i < neino; i++) {
                        igraph_integer_t edge = VECTOR(*neip)[i];
                        igraph_integer_t to = IGRAPH_OTHER(graph, edge, minnei);
                        igraph_real_t altdist = mindist + VECTOR(*weights)[edge];
                        igraph_real_t curdist = distance[to];
                        igraph_vector_int_t *v;

                        /* Note: curdist == 0 means infinity, and for this case
                         * cmp_result should be -1. However, this case is handled
                         * specially below, without referring to cmp_result. */
                        cmp_result = igraph_cmp_epsilon(altdist, curdist - 1, eps);

                        if (curdist == 0) {
                            /* This is the first finite distance to 'to' */
                            v = igraph_inclist_get(&parents, to);
                            igraph_vector_int_resize(v, 1);
                            VECTOR(*v)[0] = edge;
                            nrgeo[to] = nrgeo[minnei];
                            distance[to] = altdist + 1.0;
                            IGRAPH_CHECK(igraph_2wheap_push_with_index(&heap, to, -altdist));
                        } else if (cmp_result < 0) {
                            /* This is a shorter path */
                            v = igraph_inclist_get(&parents, to);
                            igraph_vector_int_resize(v, 1);
                            VECTOR(*v)[0] = edge;
                            nrgeo[to] = nrgeo[minnei];
                            distance[to] = altdist + 1.0;
                            igraph_2wheap_modify(&heap, to, -altdist);
                        } else if (cmp_result == 0) {
                            /* Another path with the same length */
                            v = igraph_inclist_get(&parents, to);
                            IGRAPH_CHECK(igraph_vector_int_push_back(v, edge));
                            nrgeo[to] += nrgeo[minnei];
                        }
                    }
                } /* igraph_2wheap_empty(&Q) */

                while (!igraph_stack_int_empty(&stack)) {
                    igraph_integer_t w = igraph_stack_int_pop(&stack);
                    igraph_vector_int_t *parv = igraph_inclist_get(&parents, w);
                    igraph_integer_t parv_len = igraph_vector_int_size(parv);

                    for (igraph_integer_t i = 0; i < parv_len; i++) {
                        igraph_integer_t fedge = VECTOR(*parv)[i];
                        igraph_integer_t neighbor = IGRAPH_OTHER(graph, fedge, w);
                        tmpscore[neighbor] += (tmpscore[w] + 1) * nrgeo[neighbor] / nrgeo[w];
                        VECTOR(eb)[fedge] += (tmpscore[w] + 1) * nrgeo[neighbor] / nrgeo[w];
                    }

                    tmpscore[w] = 0;
                    distance[w] = 0;
                    nrgeo[w] = 0;
                    igraph_vector_int_clear(parv);
                }
            } /* source < no_of_nodes */
        }

        /* Now look for the smallest edge betweenness */
        /* and eliminate that edge from the network */
        maxedge = igraph_i_vector_which_max_not_null(&eb, passive);
        VECTOR(*result)[e] = maxedge;
        if (edge_betweenness) {
            VECTOR(*edge_betweenness)[e] = VECTOR(eb)[maxedge];
            if (!directed) {
                VECTOR(*edge_betweenness)[e] /= 2.0;
            }
        }
        passive[maxedge] = true;
        IGRAPH_CHECK(igraph_edge(graph, maxedge, &from, &to));

        neip = igraph_inclist_get(elist_in_p, to);
        neino = igraph_vector_int_size(neip);
        igraph_vector_int_search(neip, 0, maxedge, &pos);
        VECTOR(*neip)[pos] = VECTOR(*neip)[neino - 1];
        igraph_vector_int_pop_back(neip);

        neip = igraph_inclist_get(elist_out_p, from);
        neino = igraph_vector_int_size(neip);
        igraph_vector_int_search(neip, 0, maxedge, &pos);
        VECTOR(*neip)[pos] = VECTOR(*neip)[neino - 1];
        igraph_vector_int_pop_back(neip);
    }

    IGRAPH_PROGRESS("Edge betweenness community detection: ", 100.0, NULL);

    IGRAPH_FREE(passive);
    igraph_vector_destroy(&eb);
    igraph_stack_int_destroy(&stack);
    IGRAPH_FINALLY_CLEAN(3);

    if (weights == NULL) {
        igraph_dqueue_int_destroy(&q);
        IGRAPH_FINALLY_CLEAN(1);
    } else {
        igraph_2wheap_destroy(&heap);
        igraph_inclist_destroy(&parents);
        IGRAPH_FINALLY_CLEAN(2);
    }
    igraph_free(tmpscore);
    igraph_free(nrgeo);
    igraph_free(distance);
    IGRAPH_FINALLY_CLEAN(3);

    if (directed) {
        igraph_inclist_destroy(&elist_out);
        igraph_inclist_destroy(&elist_in);
        IGRAPH_FINALLY_CLEAN(2);
    } else {
        igraph_inclist_destroy(&elist_out);
        IGRAPH_FINALLY_CLEAN(1);
    }

    if (merges || bridges || modularity || membership) {
        IGRAPH_CHECK(igraph_community_eb_get_merges(graph, directed, result, weights, merges,
                     bridges, modularity,
                     membership));
    }

    if (result_owned) {
        igraph_vector_int_destroy(result);
        IGRAPH_FREE(result);
        IGRAPH_FINALLY_CLEAN(2);
    }

    return IGRAPH_SUCCESS;
}
