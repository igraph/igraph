/* -*- mode: C -*-  */
/* vim:set ts=4 sw=4 sts=4 et: */
/*
   IGraph library.
   Copyright (C) 2007-2012  Gabor Csardi <csardi.gabor@gmail.com>
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

#include "igraph_community.h"
#include "igraph_constructors.h"
#include "igraph_memory.h"
#include "igraph_random.h"
#include "igraph_arpack.h"
#include "igraph_adjlist.h"
#include "igraph_interface.h"
#include "igraph_interrupt_internal.h"
#include "igraph_components.h"
#include "igraph_dqueue.h"
#include "igraph_progress.h"
#include "igraph_stack.h"
#include "igraph_spmatrix.h"
#include "igraph_statusbar.h"
#include "igraph_types_internal.h"
#include "igraph_conversion.h"
#include "igraph_centrality.h"
#include "igraph_structural.h"
#include "config.h"

#include <string.h>
#include <math.h>

#ifdef USING_R
    #include <R.h>
#endif

static int igraph_i_rewrite_membership_vector(igraph_vector_t *membership) {
    long int no = (long int) igraph_vector_max(membership) + 1;
    igraph_vector_t idx;
    long int realno = 0;
    long int i;
    long int len = igraph_vector_size(membership);

    IGRAPH_VECTOR_INIT_FINALLY(&idx, no);
    for (i = 0; i < len; i++) {
        long int t = (long int) VECTOR(*membership)[i];
        if (VECTOR(idx)[t]) {
            VECTOR(*membership)[i] = VECTOR(idx)[t] - 1;
        } else {
            VECTOR(idx)[t] = ++realno;
            VECTOR(*membership)[i] = VECTOR(idx)[t] - 1;
        }
    }
    igraph_vector_destroy(&idx);
    IGRAPH_FINALLY_CLEAN(1);

    return 0;
}

static int igraph_i_community_eb_get_merges2(const igraph_t *graph,
                                             const igraph_vector_t *edges,
                                             const igraph_vector_t *weights,
                                             igraph_matrix_t *res,
                                             igraph_vector_t *bridges,
                                             igraph_vector_t *modularity,
                                             igraph_vector_t *membership) {

    igraph_vector_t mymembership;
    long int no_of_nodes = igraph_vcount(graph);
    long int i;
    igraph_real_t maxmod = -1;
    long int midx = 0;
    igraph_integer_t no_comps;

    IGRAPH_VECTOR_INIT_FINALLY(&mymembership, no_of_nodes);

    if (membership) {
        IGRAPH_CHECK(igraph_vector_resize(membership, no_of_nodes));
    }

    if (modularity || res || bridges) {
        IGRAPH_CHECK(igraph_clusters(graph, 0, 0, &no_comps,
                                     IGRAPH_WEAK));

        if (modularity) {
            IGRAPH_CHECK(igraph_vector_resize(modularity,
                                              no_of_nodes - no_comps + 1));
        }
        if (res) {
            IGRAPH_CHECK(igraph_matrix_resize(res, no_of_nodes - no_comps,
                                              2));
        }
        if (bridges) {
            IGRAPH_CHECK(igraph_vector_resize(bridges,
                                              no_of_nodes - no_comps));
        }
    }

    for (i = 0; i < no_of_nodes; i++) {
        VECTOR(mymembership)[i] = i;
    }
    if (membership) {
        igraph_vector_update(membership, &mymembership);
    }

    IGRAPH_CHECK(igraph_modularity(graph, &mymembership, &maxmod, weights));
    if (modularity) {
        VECTOR(*modularity)[0] = maxmod;
    }

    for (i = igraph_vector_size(edges) - 1; i >= 0; i--) {
        long int edge = (long int) VECTOR(*edges)[i];
        long int from = IGRAPH_FROM(graph, edge);
        long int to = IGRAPH_TO(graph, edge);
        long int c1 = (long int) VECTOR(mymembership)[from];
        long int c2 = (long int) VECTOR(mymembership)[to];
        igraph_real_t actmod;
        long int j;
        if (c1 != c2) {     /* this is a merge */
            if (res) {
                MATRIX(*res, midx, 0) = c1;
                MATRIX(*res, midx, 1) = c2;
            }
            if (bridges) {
                VECTOR(*bridges)[midx] = i + 1;
            }

            /* The new cluster has id no_of_nodes+midx+1 */
            for (j = 0; j < no_of_nodes; j++) {
                if (VECTOR(mymembership)[j] == c1 ||
                    VECTOR(mymembership)[j] == c2) {
                    VECTOR(mymembership)[j] = no_of_nodes + midx;
                }
            }

            IGRAPH_CHECK(igraph_modularity(graph, &mymembership, &actmod, weights));
            if (modularity) {
                VECTOR(*modularity)[midx + 1] = actmod;
                if (actmod > maxmod) {
                    maxmod = actmod;
                    if (membership) {
                        igraph_vector_update(membership, &mymembership);
                    }
                }
            }

            midx++;
        }
    }

    if (membership) {
        IGRAPH_CHECK(igraph_i_rewrite_membership_vector(membership));
    }

    igraph_vector_destroy(&mymembership);
    IGRAPH_FINALLY_CLEAN(1);

    return 0;
}


/**
 * \function igraph_community_eb_get_merges
 * \brief Calculating the merges, i.e. the dendrogram for an edge betweenness community structure
 *
 * </para><para>
 * This function is handy if you have a sequence of edge which are
 * gradually removed from the network and you would like to know how
 * the network falls apart into separate components. The edge sequence
 * may come from the \ref igraph_community_edge_betweenness()
 * function, but this is not necessary. Note that \ref
 * igraph_community_edge_betweenness can also calculate the
 * dendrogram, via its \p merges argument.
 *
 * \param graph The input graph.
 * \param edges Vector containing the edges to be removed from the
 *    network, all edges are expected to appear exactly once in the
 *    vector.
 * \param weights An optional vector containing edge weights. If null,
 *     the unweighted modularity scores will be calculated. If not null,
 *     the weighted modularity scores will be calculated. Ignored if both
 *     \p modularity and \p membership are nulls.
 * \param res Pointer to an initialized matrix, if not NULL then the
 *    dendrogram will be stored here, in the same form as for the \ref
 *    igraph_community_walktrap() function: the matrix has two columns
 *    and each line is a merge given by the ids of the merged
 *    components. The component ids are number from zero and
 *    component ids smaller than the number of vertices in the graph
 *    belong to individual vertices. The non-trivial components
 *    containing at least two vertices are numbered from \c n, \c n is
 *    the number of vertices in the graph. So if the first line
 *    contains \c a and \c b that means that components \c a and \c b
 *    are merged into component \c n, the second line creates
 *    component \c n+1, etc. The matrix will be resized as needed.
 * \param bridges Pointer to an initialized vector or NULL. If not
 *    null then the index of the edge removals which split the network
 *    will be stored here. The vector will be resized as needed.
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

int igraph_community_eb_get_merges(const igraph_t *graph,
                                   const igraph_vector_t *edges,
                                   const igraph_vector_t *weights,
                                   igraph_matrix_t *res,
                                   igraph_vector_t *bridges,
                                   igraph_vector_t *modularity,
                                   igraph_vector_t *membership) {

    long int no_of_nodes = igraph_vcount(graph);
    igraph_vector_t ptr;
    long int i, midx = 0;
    igraph_integer_t no_comps;

    if (membership || modularity) {
        return igraph_i_community_eb_get_merges2(graph, edges, weights, res,
                bridges, modularity,
                membership);
    }

    IGRAPH_CHECK(igraph_clusters(graph, 0, 0, &no_comps, IGRAPH_WEAK));

    IGRAPH_VECTOR_INIT_FINALLY(&ptr, no_of_nodes * 2 - 1);
    if (res) {
        IGRAPH_CHECK(igraph_matrix_resize(res, no_of_nodes - no_comps, 2));
    }
    if (bridges) {
        IGRAPH_CHECK(igraph_vector_resize(bridges, no_of_nodes - no_comps));
    }

    for (i = igraph_vector_size(edges) - 1; i >= 0; i--) {
        igraph_integer_t edge = (igraph_integer_t) VECTOR(*edges)[i];
        igraph_integer_t from, to, c1, c2, idx;
        igraph_edge(graph, edge, &from, &to);
        idx = from + 1;
        while (VECTOR(ptr)[idx - 1] != 0) {
            idx = (igraph_integer_t) VECTOR(ptr)[idx - 1];
        }
        c1 = idx - 1;
        idx = to + 1;
        while (VECTOR(ptr)[idx - 1] != 0) {
            idx = (igraph_integer_t) VECTOR(ptr)[idx - 1];
        }
        c2 = idx - 1;
        if (c1 != c2) {     /* this is a merge */
            if (res) {
                MATRIX(*res, midx, 0) = c1;
                MATRIX(*res, midx, 1) = c2;
            }
            if (bridges) {
                VECTOR(*bridges)[midx] = i + 1;
            }

            VECTOR(ptr)[c1] = no_of_nodes + midx + 1;
            VECTOR(ptr)[c2] = no_of_nodes + midx + 1;
            VECTOR(ptr)[from] = no_of_nodes + midx + 1;
            VECTOR(ptr)[to] = no_of_nodes + midx + 1;

            midx++;
        }
    }

    igraph_vector_destroy(&ptr);
    IGRAPH_FINALLY_CLEAN(1);

    return 0;
}

/* Find the smallest active element in the vector */
static long int igraph_i_vector_which_max_not_null(const igraph_vector_t *v,
                                                   const char *passive) {
    long int which, i = 0, size = igraph_vector_size(v);
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
 * \brief Community finding based on edge betweenness
 *
 * Community structure detection based on the betweenness of the edges
 * in the network. The algorithm was invented by M. Girvan and
 * M. Newman, see: M. Girvan and M. E. J. Newman: Community structure in
 * social and biological networks, Proc. Nat. Acad. Sci. USA 99, 7821-7826
 * (2002).
 *
 * </para><para>
 * The idea is that the betweenness of the edges connecting two
 * communities is typically high, as many of the shortest paths
 * between nodes in separate communities go through them. So we
 * gradually remove the edge with highest betweenness from the
 * network, and recalculate edge betweenness after every removal.
 * This way sooner or later the network falls off to two components,
 * then after a while one of these components falls off to two smaller
 * components, etc. until all edges are removed. This is a divisive
 * hierarchical approach, the result is a dendrogram.
 * \param graph The input graph.
 * \param result Pointer to an initialized vector, the result will be
 *     stored here, the ids of the removed edges in the order of their
 *     removal. It will be resized as needed. It may be NULL if
 *     the edge IDs are not needed by the caller.
 * \param edge_betweenness Pointer to an initialized vector or
 *     NULL. In the former case the edge betweenness of the removed
 *     edge is stored here. The vector will be resized as needed.
 * \param merges Pointer to an initialized matrix or NULL. If not NULL
 *     then merges performed by the algorithm are stored here. Even if
 *     this is a divisive algorithm, we can replay it backwards and
 *     note which two clusters were merged. Clusters are numbered from
 *     zero, see the \p merges argument of \ref
 *     igraph_community_walktrap() for details. The matrix will be
 *     resized as needed.
 * \param bridges Pointer to an initialized vector of NULL. If not
 *     NULL then all edge removals which separated the network into
 *     more components are marked here.
 * \param modularity If not a null pointer, then the modularity values
 *     of the different divisions are stored here, in the order
 *     corresponding to the merge matrix. The modularity values will
 *     take weights into account if \p weights is not null.
 * \param membership If not a null pointer, then the membership vector,
 *     corresponding to the highest modularity value, is stored here.
 * \param directed Logical constant, whether to calculate directed
 *    betweenness (i.e. directed paths) for directed graphs. It is
 *    ignored for undirected graphs.
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

int igraph_community_edge_betweenness(const igraph_t *graph,
                                      igraph_vector_t *result,
                                      igraph_vector_t *edge_betweenness,
                                      igraph_matrix_t *merges,
                                      igraph_vector_t *bridges,
                                      igraph_vector_t *modularity,
                                      igraph_vector_t *membership,
                                      igraph_bool_t directed,
                                      const igraph_vector_t *weights) {

    long int no_of_nodes = igraph_vcount(graph);
    long int no_of_edges = igraph_ecount(graph);
    double *distance, *tmpscore;
    unsigned long long int *nrgeo;
    long int source, i, e;

    igraph_inclist_t elist_out, elist_in, fathers;
    igraph_inclist_t *elist_out_p, *elist_in_p;
    igraph_vector_int_t *neip;
    long int neino;
    igraph_vector_t eb;
    long int maxedge, pos;
    igraph_integer_t from, to;
    igraph_bool_t result_owned = 0;
    igraph_stack_t stack = IGRAPH_STACK_NULL;
    igraph_real_t steps, steps_done;

    char *passive;

    /* Needed only for the unweighted case */
    igraph_dqueue_t q = IGRAPH_DQUEUE_NULL;

    /* Needed only for the weighted case */
    igraph_2wheap_t heap;

    if (result == 0) {
        result = igraph_Calloc(1, igraph_vector_t);
        if (result == 0) {
            IGRAPH_ERROR("edge betweenness community structure failed", IGRAPH_ENOMEM);
        }
        IGRAPH_FINALLY(igraph_free, result);
        IGRAPH_VECTOR_INIT_FINALLY(result, 0);
        result_owned = 1;
    }

    directed = directed && igraph_is_directed(graph);
    if (directed) {
        IGRAPH_CHECK(igraph_inclist_init(graph, &elist_out, IGRAPH_OUT));
        IGRAPH_FINALLY(igraph_inclist_destroy, &elist_out);
        IGRAPH_CHECK(igraph_inclist_init(graph, &elist_in, IGRAPH_IN));
        IGRAPH_FINALLY(igraph_inclist_destroy, &elist_in);
        elist_out_p = &elist_out;
        elist_in_p = &elist_in;
    } else {
        IGRAPH_CHECK(igraph_inclist_init(graph, &elist_out, IGRAPH_ALL));
        IGRAPH_FINALLY(igraph_inclist_destroy, &elist_out);
        elist_out_p = elist_in_p = &elist_out;
    }

    distance = igraph_Calloc(no_of_nodes, double);
    if (distance == 0) {
        IGRAPH_ERROR("edge betweenness community structure failed", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, distance);
    nrgeo = igraph_Calloc(no_of_nodes, unsigned long long int);
    if (nrgeo == 0) {
        IGRAPH_ERROR("edge betweenness community structure failed", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, nrgeo);
    tmpscore = igraph_Calloc(no_of_nodes, double);
    if (tmpscore == 0) {
        IGRAPH_ERROR("edge betweenness community structure failed", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, tmpscore);

    if (weights == 0) {
        IGRAPH_DQUEUE_INIT_FINALLY(&q, 100);
    } else {
        if (igraph_vector_min(weights) <= 0) {
            IGRAPH_ERROR("weights must be strictly positive", IGRAPH_EINVAL);
        }

        if (membership != 0) {
            IGRAPH_WARNING("Membership vector will be selected based on the lowest "\
                           "modularity score.");
        }

        if (modularity != 0 || membership != 0) {
            IGRAPH_WARNING("Modularity calculation with weighted edge betweenness "\
                           "community detection might not make sense -- modularity treats edge "\
                           "weights as similarities while edge betwenness treats them as "\
                           "distances");
        }

        IGRAPH_CHECK(igraph_2wheap_init(&heap, no_of_nodes));
        IGRAPH_FINALLY(igraph_2wheap_destroy, &heap);
        IGRAPH_CHECK(igraph_inclist_init_empty(&fathers,
                                               (igraph_integer_t) no_of_nodes));
        IGRAPH_FINALLY(igraph_inclist_destroy, &fathers);
    }

    IGRAPH_CHECK(igraph_stack_init(&stack, no_of_nodes));
    IGRAPH_FINALLY(igraph_stack_destroy, &stack);

    IGRAPH_CHECK(igraph_vector_resize(result, no_of_edges));
    if (edge_betweenness) {
        IGRAPH_CHECK(igraph_vector_resize(edge_betweenness, no_of_edges));
        if (no_of_edges > 0) {
            VECTOR(*edge_betweenness)[no_of_edges - 1] = 0;
        }
    }

    IGRAPH_VECTOR_INIT_FINALLY(&eb, no_of_edges);

    passive = igraph_Calloc(no_of_edges, char);
    if (!passive) {
        IGRAPH_ERROR("edge betweenness community structure failed", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, passive);

    /* Estimate the number of steps to be taken.
     * It is assumed that one iteration is O(|E||V|), but |V| is constant
     * anyway, so we will have approximately |E|^2 / 2 steps, and one
     * iteration of the outer loop advances the step counter by the number
     * of remaining edges at that iteration.
     */
    steps = no_of_edges / 2.0 * (no_of_edges + 1);
    steps_done = 0;

    for (e = 0; e < no_of_edges; steps_done += no_of_edges - e, e++) {
        IGRAPH_PROGRESS("Edge betweenness community detection: ",
                        100.0 * steps_done / steps, NULL);

        igraph_vector_null(&eb);

        if (weights == 0) {
            /* Unweighted variant follows */

            /* The following for loop is copied almost intact from
             * igraph_edge_betweenness_estimate */
            for (source = 0; source < no_of_nodes; source++) {

                IGRAPH_ALLOW_INTERRUPTION();

                memset(distance, 0, (size_t) no_of_nodes * sizeof(double));
                memset(nrgeo, 0, (size_t) no_of_nodes * sizeof(unsigned long long int));
                memset(tmpscore, 0, (size_t) no_of_nodes * sizeof(double));
                igraph_stack_clear(&stack); /* it should be empty anyway... */

                IGRAPH_CHECK(igraph_dqueue_push(&q, source));

                nrgeo[source] = 1;
                distance[source] = 0;

                while (!igraph_dqueue_empty(&q)) {
                    long int actnode = (long int) igraph_dqueue_pop(&q);

                    neip = igraph_inclist_get(elist_out_p, actnode);
                    neino = igraph_vector_int_size(neip);
                    for (i = 0; i < neino; i++) {
                        igraph_integer_t edge = (igraph_integer_t) VECTOR(*neip)[i], from, to;
                        long int neighbor;
                        igraph_edge(graph, edge, &from, &to);
                        neighbor = actnode != from ? from : to;
                        if (nrgeo[neighbor] != 0) {
                            /* we've already seen this node, another shortest path? */
                            if (distance[neighbor] == distance[actnode] + 1) {
                                nrgeo[neighbor] += nrgeo[actnode];
                            }
                        } else {
                            /* we haven't seen this node yet */
                            nrgeo[neighbor] += nrgeo[actnode];
                            distance[neighbor] = distance[actnode] + 1;
                            IGRAPH_CHECK(igraph_dqueue_push(&q, neighbor));
                            IGRAPH_CHECK(igraph_stack_push(&stack, neighbor));
                        }
                    }
                } /* while !igraph_dqueue_empty */

                /* Ok, we've the distance of each node and also the number of
                   shortest paths to them. Now we do an inverse search, starting
                   with the farthest nodes. */
                while (!igraph_stack_empty(&stack)) {
                    long int actnode = (long int) igraph_stack_pop(&stack);
                    if (distance[actnode] < 1) {
                        continue;    /* skip source node */
                    }

                    /* set the temporary score of the friends */
                    neip = igraph_inclist_get(elist_in_p, actnode);
                    neino = igraph_vector_int_size(neip);
                    for (i = 0; i < neino; i++) {
                        long int edge = (long int) VECTOR(*neip)[i];
                        long int neighbor = IGRAPH_OTHER(graph, edge, actnode);
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

            /* The following for loop is copied almost intact from
             * igraph_i_edge_betweenness_estimate_weighted */
            for (source = 0; source < no_of_nodes; source++) {
                /* This will contain the edge betweenness in the current step */
                IGRAPH_ALLOW_INTERRUPTION();

                memset(distance, 0, (size_t) no_of_nodes * sizeof(double));
                memset(nrgeo, 0, (size_t) no_of_nodes * sizeof(unsigned long long int));
                memset(tmpscore, 0, (size_t) no_of_nodes * sizeof(double));

                igraph_2wheap_push_with_index(&heap, source, 0);
                distance[source] = 1.0;
                nrgeo[source] = 1;

                while (!igraph_2wheap_empty(&heap)) {
                    long int minnei = igraph_2wheap_max_index(&heap);
                    igraph_real_t mindist = -igraph_2wheap_delete_max(&heap);

                    igraph_stack_push(&stack, minnei);

                    neip = igraph_inclist_get(elist_out_p, minnei);
                    neino = igraph_vector_int_size(neip);

                    for (i = 0; i < neino; i++) {
                        long int edge = VECTOR(*neip)[i];
                        long int to = IGRAPH_OTHER(graph, edge, minnei);
                        igraph_real_t altdist = mindist + VECTOR(*weights)[edge];
                        igraph_real_t curdist = distance[to];
                        igraph_vector_int_t *v;

                        if (curdist == 0) {
                            /* This is the first finite distance to 'to' */
                            v = igraph_inclist_get(&fathers, to);
                            igraph_vector_int_resize(v, 1);
                            VECTOR(*v)[0] = edge;
                            nrgeo[to] = nrgeo[minnei];
                            distance[to] = altdist + 1.0;
                            IGRAPH_CHECK(igraph_2wheap_push_with_index(&heap, to, -altdist));
                        } else if (altdist < curdist - 1) {
                            /* This is a shorter path */
                            v = igraph_inclist_get(&fathers, to);
                            igraph_vector_int_resize(v, 1);
                            VECTOR(*v)[0] = edge;
                            nrgeo[to] = nrgeo[minnei];
                            distance[to] = altdist + 1.0;
                            IGRAPH_CHECK(igraph_2wheap_modify(&heap, to, -altdist));
                        } else if (altdist == curdist - 1) {
                            /* Another path with the same length */
                            v = igraph_inclist_get(&fathers, to);
                            igraph_vector_int_push_back(v, edge);
                            nrgeo[to] += nrgeo[minnei];
                        }
                    }
                } /* igraph_2wheap_empty(&Q) */

                while (!igraph_stack_empty(&stack)) {
                    long int w = (long int) igraph_stack_pop(&stack);
                    igraph_vector_int_t *fatv = igraph_inclist_get(&fathers, w);
                    long int fatv_len = igraph_vector_int_size(fatv);

                    for (i = 0; i < fatv_len; i++) {
                        long int fedge = (long int) VECTOR(*fatv)[i];
                        long int neighbor = IGRAPH_OTHER(graph, fedge, w);
                        tmpscore[neighbor] += (tmpscore[w] + 1) * nrgeo[neighbor] / nrgeo[w];
                        VECTOR(eb)[fedge] += (tmpscore[w] + 1) * nrgeo[neighbor] / nrgeo[w];
                    }

                    tmpscore[w] = 0;
                    distance[w] = 0;
                    nrgeo[w] = 0;
                    igraph_vector_int_clear(fatv);
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
        passive[maxedge] = 1;
        igraph_edge(graph, (igraph_integer_t) maxedge, &from, &to);

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

    igraph_free(passive);
    igraph_vector_destroy(&eb);
    igraph_stack_destroy(&stack);
    IGRAPH_FINALLY_CLEAN(3);

    if (weights == 0) {
        igraph_dqueue_destroy(&q);
        IGRAPH_FINALLY_CLEAN(1);
    } else {
        igraph_2wheap_destroy(&heap);
        igraph_inclist_destroy(&fathers);
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
        IGRAPH_CHECK(igraph_community_eb_get_merges(graph, result, weights, merges,
                     bridges, modularity,
                     membership));
    }

    if (result_owned) {
        igraph_vector_destroy(result);
        igraph_Free(result);
        IGRAPH_FINALLY_CLEAN(2);
    }

    return 0;
}


/**
 * \function igraph_community_to_membership
 * \brief Create membership vector from community structure dendrogram
 *
 * This function creates a membership vector from a community
 * structure dendrogram. A membership vector contains for each vertex
 * the id of its graph component, the graph components are numbered
 * from zero, see the same argument of \ref igraph_clusters() for an
 * example of a membership vector.
 *
 * </para><para>
 * Many community detection algorithms return with a \em merges
 * matrix, \ref igraph_community_walktrap() and \ref
 * igraph_community_edge_betweenness() are two examples. The matrix
 * contains the merge operations performed while mapping the
 * hierarchical structure of a network. If the matrix has \c n-1 rows,
 * where \c n is the number of vertices in the graph, then it contains
 * the hierarchical structure of the whole network and it is called a
 * dendrogram.
 *
 * </para><para>
 * This function performs \p steps merge operations as prescribed by
 * the \p merges matrix and returns the current state of the network.
 *
 * </para><para>
 * If \p merges is not a complete dendrogram, it is possible to
 * take \p steps steps if \p steps is not bigger than the number
 * lines in \p merges.
 * \param merges The two-column matrix containing the merge
 *    operations. See \ref igraph_community_walktrap() for the
 *    detailed syntax.
 * \param nodes The number of leaf nodes in the dendrogram
 * \param steps Integer constant, the number of steps to take.
 * \param membership Pointer to an initialized vector, the membership
 *    results will be stored here, if not NULL. The vector will be
 *    resized as needed.
 * \param csize Pointer to an initialized vector, or NULL. If not NULL
 *    then the sizes of the components will be stored here, the vector
 *    will be resized as needed.
 *
 * \sa \ref igraph_community_walktrap(), \ref
 * igraph_community_edge_betweenness(), \ref
 * igraph_community_fastgreedy() for community structure detection
 * algorithms.
 *
 * Time complexity: O(|V|), the number of vertices in the graph.
 */

int igraph_community_to_membership(const igraph_matrix_t *merges,
                                   igraph_integer_t nodes,
                                   igraph_integer_t steps,
                                   igraph_vector_t *membership,
                                   igraph_vector_t *csize) {

    long int no_of_nodes = nodes;
    long int components = no_of_nodes - steps;
    long int i, found = 0;
    igraph_vector_t tmp;

    if (steps > igraph_matrix_nrow(merges)) {
        IGRAPH_ERROR("`steps' to big or `merges' matrix too short", IGRAPH_EINVAL);
    }

    if (membership) {
        IGRAPH_CHECK(igraph_vector_resize(membership, no_of_nodes));
        igraph_vector_null(membership);
    }
    if (csize) {
        IGRAPH_CHECK(igraph_vector_resize(csize, components));
        igraph_vector_null(csize);
    }

    IGRAPH_VECTOR_INIT_FINALLY(&tmp, steps);

    for (i = steps - 1; i >= 0; i--) {
        long int c1 = (long int) MATRIX(*merges, i, 0);
        long int c2 = (long int) MATRIX(*merges, i, 1);

        /* new component? */
        if (VECTOR(tmp)[i] == 0) {
            found++;
            VECTOR(tmp)[i] = found;
        }

        if (c1 < no_of_nodes) {
            long int cid = (long int) VECTOR(tmp)[i] - 1;
            if (membership) {
                VECTOR(*membership)[c1] = cid + 1;
            }
            if (csize) {
                VECTOR(*csize)[cid] += 1;
            }
        } else {
            VECTOR(tmp)[c1 - no_of_nodes] = VECTOR(tmp)[i];
        }

        if (c2 < no_of_nodes) {
            long int cid = (long int) VECTOR(tmp)[i] - 1;
            if (membership) {
                VECTOR(*membership)[c2] = cid + 1;
            }
            if (csize) {
                VECTOR(*csize)[cid] += 1;
            }
        } else {
            VECTOR(tmp)[c2 - no_of_nodes] = VECTOR(tmp)[i];
        }

    }

    if (membership || csize) {
        for (i = 0; i < no_of_nodes; i++) {
            long int tmp = (long int) VECTOR(*membership)[i];
            if (tmp != 0) {
                if (membership) {
                    VECTOR(*membership)[i] = tmp - 1;
                }
            } else {
                if (csize) {
                    VECTOR(*csize)[found] += 1;
                }
                if (membership) {
                    VECTOR(*membership)[i] = found;
                }
                found++;
            }
        }
    }

    igraph_vector_destroy(&tmp);
    IGRAPH_FINALLY_CLEAN(1);

    return 0;
}

/**
 * \function igraph_modularity
 * \brief Calculate the modularity of a graph with respect to some vertex types
 *
 * The modularity of a graph with respect to some division (or vertex
 * types) measures how good the division is, or how separated are the
 * different vertex types from each other. It is defined as
 * Q=1/(2m) * sum((Aij - ki*kj / (2m)) delta(ci,cj), i, j), here `m' is the
 * number of edges, `Aij' is the element of the `A' adjacency matrix
 * in row `i' and column `j', `ki' is the degree of `i', `kj' is the
 * degree of `j', `ci' is the type (or component) of `i', `cj' that of
 * `j', the sum goes over all `i' and `j' pairs of vertices, and
 * `delta(x,y)' is one if x=y and zero otherwise.
 *
 * </para><para>
 * Modularity on weighted graphs is also meaningful. When taking edge
 * weights into account, `Aij' becomes the weight of the corresponding
 * edge (or 0 if there is no edge), `ki' is the total weight of edges
 * incident on vertex `i', `kj' is the total weight of edges incident
 * on vertex `j' and `m' is the total weight of all edges.
 *
 * </para><para>
 * See also Clauset, A.; Newman, M. E. J.; Moore, C. Finding
 * community structure in very large networks, Physical Review E,
 * 2004, 70, 066111.
 * \param graph The input graph. It must be undirected; directed graphs are
 *     not supported yet.
 * \param membership Numeric vector which gives the type of each
 *     vertex, i.e. the component to which it belongs.
 *     It does not have to be consecutive, i.e. empty communities are
 *     allowed.
 * \param modularity Pointer to a real number, the result will be
 *     stored here.
 * \param weights Weight vector or NULL if no weights are specified.
 * \return Error code.
 *
 * Time complexity: O(|V|+|E|), the number of vertices plus the number
 * of edges.
 */

int igraph_modularity(const igraph_t *graph,
                      const igraph_vector_t *membership,
                      igraph_real_t *modularity,
                      const igraph_vector_t *weights) {

    igraph_vector_t e, a;
    long int types = (long int) igraph_vector_max(membership) + 1;
    long int no_of_edges = igraph_ecount(graph);
    long int i;
    igraph_integer_t from, to;
    igraph_real_t m;
    long int c1, c2;

    if (igraph_is_directed(graph)) {
#ifndef USING_R
        IGRAPH_ERROR("modularity is implemented for undirected graphs", IGRAPH_EINVAL);
#else
        REprintf("Modularity is implemented for undirected graphs only.\n");
#endif
    }

    if (igraph_vector_size(membership) < igraph_vcount(graph)) {
        IGRAPH_ERROR("cannot calculate modularity, membership vector too short",
                     IGRAPH_EINVAL);
    }
    if (igraph_vector_min(membership) < 0) {
        IGRAPH_ERROR("Invalid membership vector", IGRAPH_EINVAL);
    }

    IGRAPH_VECTOR_INIT_FINALLY(&e, types);
    IGRAPH_VECTOR_INIT_FINALLY(&a, types);

    if (weights) {
        if (igraph_vector_size(weights) < no_of_edges)
            IGRAPH_ERROR("cannot calculate modularity, weight vector too short",
                         IGRAPH_EINVAL);
        m = igraph_vector_sum(weights);
        for (i = 0; i < no_of_edges; i++) {
            igraph_real_t w = VECTOR(*weights)[i];
            if (w < 0) {
                IGRAPH_ERROR("negative weight in weight vector", IGRAPH_EINVAL);
            }
            igraph_edge(graph, (igraph_integer_t) i, &from, &to);
            c1 = (long int) VECTOR(*membership)[from];
            c2 = (long int) VECTOR(*membership)[to];
            if (c1 == c2) {
                VECTOR(e)[c1] += 2 * w;
            }
            VECTOR(a)[c1] += w;
            VECTOR(a)[c2] += w;
        }
    } else {
        m = no_of_edges;
        for (i = 0; i < no_of_edges; i++) {
            igraph_edge(graph, (igraph_integer_t) i, &from, &to);
            c1 = (long int) VECTOR(*membership)[from];
            c2 = (long int) VECTOR(*membership)[to];
            if (c1 == c2) {
                VECTOR(e)[c1] += 2;
            }
            VECTOR(a)[c1] += 1;
            VECTOR(a)[c2] += 1;
        }
    }

    *modularity = 0.0;
    if (m > 0) {
        for (i = 0; i < types; i++) {
            igraph_real_t tmp = VECTOR(a)[i] / 2 / m;
            *modularity += VECTOR(e)[i] / 2 / m;
            *modularity -= tmp * tmp;
        }
    }

    igraph_vector_destroy(&e);
    igraph_vector_destroy(&a);
    IGRAPH_FINALLY_CLEAN(2);

    return 0;
}

/**
 * \function igraph_modularity_matrix
 * \brief Calculate the modularity matrix
 *
 * This function returns the modularity matrix defined as
 * `B_ij = A_ij - k_i k_j * / 2 m`
 * where `A_ij` denotes the adjacency matrix, `k_i` is the degree of node `i`
 * and `m` is the total weight in the graph. Note that self-loops are multiplied
 * by 2 in this implementation. If weights are specified, the weighted
 * counterparts are used.
 *
 * \param graph   The input graph
 * \param modmat  Pointer to an initialized matrix in which the modularity
 *                matrix is stored.
 * \param weights Edge weights, pointer to a vector. If this is a null pointer
 *                then every edge is assumed to have a weight of 1.
 */

int igraph_modularity_matrix(const igraph_t *graph,
                             igraph_matrix_t *modmat,
                             const igraph_vector_t *weights) {

    long int no_of_nodes = igraph_vcount(graph);
    long int no_of_edges = igraph_ecount(graph);
    igraph_real_t sw = weights ? igraph_vector_sum(weights) : no_of_edges;
    igraph_vector_t deg;
    long int i, j;

    if (weights && igraph_vector_size(weights) != no_of_edges) {
        IGRAPH_ERROR("Invalid weight vector length", IGRAPH_EINVAL);
    }

    IGRAPH_VECTOR_INIT_FINALLY(&deg, no_of_nodes);
    if (!weights) {
        IGRAPH_CHECK(igraph_degree(graph, &deg, igraph_vss_all(), IGRAPH_ALL,
                                   IGRAPH_LOOPS));
    } else {
        IGRAPH_CHECK(igraph_strength(graph, &deg, igraph_vss_all(), IGRAPH_ALL,
                                     IGRAPH_LOOPS, weights));
    }
    IGRAPH_CHECK(igraph_get_adjacency(graph, modmat, IGRAPH_GET_ADJACENCY_BOTH,
                                      /*eids=*/ 0));

    for (i = 0; i < no_of_nodes; i++) {
        MATRIX(*modmat, i, i) *= 2;
    }
    for (i = 0; i < no_of_nodes; i++) {
        for (j = 0; j < no_of_nodes; j++) {
            MATRIX(*modmat, i, j) -= VECTOR(deg)[i] * VECTOR(deg)[j] / 2.0 / sw;
        }
    }

    igraph_vector_destroy(&deg);
    IGRAPH_FINALLY_CLEAN(1);

    return 0;
}

/**
 * \function igraph_reindex_membership
 * \brief Makes the IDs in a membership vector continuous
 *
 * This function reindexes component IDs in a membership vector
 * in a way that the new IDs start from zero and go up to C-1,
 * where C is the number of unique component IDs in the original
 * vector. The supplied membership is expected to fall in the
 * range 0, ..., n - 1.
 *
 * \param  membership  Numeric vector which gives the type of each
 *                     vertex, i.e. the component to which it belongs.
 *                     The vector will be altered in-place.
 * \param  new_to_old  Pointer to a vector which will contain the
 *                     old component ID for each new one, or NULL,
 *                     in which case it is not returned. The vector
 *                     will be resized as needed.
 * \param  nb_clusters Pointer to an integer for the number of
 *                     distinct clusters. If not NULL, this will be
 *                     updated to reflect the number of distinct
 *                     clusters found in membership.
 *
 * Time complexity: should be O(n) for n elements.
 */
int igraph_reindex_membership(igraph_vector_t *membership,
                              igraph_vector_t *new_to_old,
                              igraph_integer_t *nb_clusters) {

    long int i, n = igraph_vector_size(membership);
    igraph_vector_t new_cluster;
    igraph_integer_t i_nb_clusters;

    /* We allow original cluster indices in the range 0, ..., n - 1 */
    IGRAPH_CHECK(igraph_vector_init(&new_cluster, n));
    IGRAPH_FINALLY(igraph_vector_destroy, &new_cluster);

    if (new_to_old) {
        igraph_vector_clear(new_to_old);
    }

    /* Clean clusters. We will store the new cluster + 1 so that membership == 0
     * indicates that no cluster was assigned yet. */
    i_nb_clusters = 1;
    for (i = 0; i < n; i++) {
        long int c = (long int)VECTOR(*membership)[i];

        if (c >= n) {
            IGRAPH_ERROR("Cluster out of range", IGRAPH_EINVAL);
        }

        if (VECTOR(new_cluster)[c] == 0) {
            VECTOR(new_cluster)[c] = (igraph_real_t)i_nb_clusters;
            i_nb_clusters += 1;
            if (new_to_old) {
                IGRAPH_CHECK(igraph_vector_push_back(new_to_old, c));
            }
        }
    }

    /* Assign new membership */
    for (i = 0; i < n; i++) {
        long int c = (long int)VECTOR(*membership)[i];
        VECTOR(*membership)[i] = VECTOR(new_cluster)[c] - 1;
    }
    if (nb_clusters) {
        /* We used the cluster + 1, so correct */
        *nb_clusters = i_nb_clusters - 1;
    }

    igraph_vector_destroy(&new_cluster);

    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}

/********************************************************************/

/**
 * \section about_leading_eigenvector_methods
 *
 * <para>
 * The function documented in these section implements the
 * <quote>leading eigenvector</quote> method developed by Mark Newman and
 * published in MEJ Newman: Finding community structure using the
 * eigenvectors of matrices, Phys Rev E 74:036104 (2006).</para>
 *
 * <para>
 * The heart of the method is the definition of the modularity matrix,
 * B, which is B=A-P, A being the adjacency matrix of the (undirected)
 * network, and P contains the probability that certain edges are
 * present according to the <quote>configuration model</quote> In
 * other words, a Pij element of P is the probability that there is an
 * edge between vertices i and j in a random network in which the
 * degrees of all vertices are the same as in the input graph.</para>
 *
 * <para>
 * The leading eigenvector method works by calculating the eigenvector
 * of the modularity matrix for the largest positive eigenvalue and
 * then separating vertices into two community based on the sign of
 * the corresponding element in the eigenvector. If all elements in
 * the eigenvector are of the same sign that means that the network
 * has no underlying community structure.
 * Check Newman's paper to understand why this is a good method for
 * detecting community structure. </para>
 *
 * <para>
 * The leading eigenvector community structure detection method is
 * implemented in \ref igraph_community_leading_eigenvector(). After
 * the initial split, the following splits are done in a way to
 * optimize modularity regarding to the original network. Note that
 * any further refinement, for example using Kernighan-Lin, as
 * proposed in Section V.A of Newman (2006), is not implemented here.
 * </para>
 *
 * <para>
 * \example examples/simple/igraph_community_leading_eigenvector.c
 * </para>
 */

typedef struct igraph_i_community_leading_eigenvector_data_t {
    igraph_vector_t *idx;
    igraph_vector_t *idx2;
    igraph_adjlist_t *adjlist;
    igraph_inclist_t *inclist;
    igraph_vector_t *tmp;
    long int no_of_edges;
    igraph_vector_t *mymembership;
    long int comm;
    const igraph_vector_t *weights;
    const igraph_t *graph;
    igraph_vector_t *strength;
    igraph_real_t sumweights;
} igraph_i_community_leading_eigenvector_data_t;

static int igraph_i_community_leading_eigenvector(igraph_real_t *to,
                                                  const igraph_real_t *from,
                                                  int n, void *extra) {

    igraph_i_community_leading_eigenvector_data_t *data = extra;
    long int j, k, nlen, size = n;
    igraph_vector_t *idx = data->idx;
    igraph_vector_t *idx2 = data->idx2;
    igraph_vector_t *tmp = data->tmp;
    igraph_adjlist_t *adjlist = data->adjlist;
    igraph_real_t ktx, ktx2;
    long int no_of_edges = data->no_of_edges;
    igraph_vector_t *mymembership = data->mymembership;
    long int comm = data->comm;

    /* Ax */
    for (j = 0; j < size; j++) {
        long int oldid = (long int) VECTOR(*idx)[j];
        igraph_vector_int_t *neis = igraph_adjlist_get(adjlist, oldid);
        nlen = igraph_vector_int_size(neis);
        to[j] = 0.0;
        VECTOR(*tmp)[j] = 0.0;
        for (k = 0; k < nlen; k++) {
            long int nei = (long int) VECTOR(*neis)[k];
            long int neimemb = (long int) VECTOR(*mymembership)[nei];
            if (neimemb == comm) {
                to[j] += from[ (long int) VECTOR(*idx2)[nei] ];
                VECTOR(*tmp)[j] += 1;
            }
        }
    }

    /* Now calculate k^Tx/2m */
    ktx = 0.0; ktx2 = 0.0;
    for (j = 0; j < size; j++) {
        long int oldid = (long int) VECTOR(*idx)[j];
        igraph_vector_int_t *neis = igraph_adjlist_get(adjlist, oldid);
        long int degree = igraph_vector_int_size(neis);
        ktx += from[j] * degree;
        ktx2 += degree;
    }
    ktx = ktx / no_of_edges / 2.0;
    ktx2 = ktx2 / no_of_edges / 2.0;

    /* Now calculate Bx */
    for (j = 0; j < size; j++) {
        long int oldid = (long int) VECTOR(*idx)[j];
        igraph_vector_int_t *neis = igraph_adjlist_get(adjlist, oldid);
        igraph_real_t degree = igraph_vector_int_size(neis);
        to[j] = to[j] - ktx * degree;
        VECTOR(*tmp)[j] = VECTOR(*tmp)[j] - ktx2 * degree;
    }

    /* -d_ij summa l in G B_il */
    for (j = 0; j < size; j++) {
        to[j] -= VECTOR(*tmp)[j] * from[j];
    }

    return 0;
}

static int igraph_i_community_leading_eigenvector2(igraph_real_t *to,
                                                   const igraph_real_t *from,
                                                   int n, void *extra) {

    igraph_i_community_leading_eigenvector_data_t *data = extra;
    long int j, k, nlen, size = n;
    igraph_vector_t *idx = data->idx;
    igraph_vector_t *idx2 = data->idx2;
    igraph_vector_t *tmp = data->tmp;
    igraph_adjlist_t *adjlist = data->adjlist;
    igraph_real_t ktx, ktx2;
    long int no_of_edges = data->no_of_edges;
    igraph_vector_t *mymembership = data->mymembership;
    long int comm = data->comm;

    /* Ax */
    for (j = 0; j < size; j++) {
        long int oldid = (long int) VECTOR(*idx)[j];
        igraph_vector_int_t *neis = igraph_adjlist_get(adjlist, oldid);
        nlen = igraph_vector_int_size(neis);
        to[j] = 0.0;
        VECTOR(*tmp)[j] = 0.0;
        for (k = 0; k < nlen; k++) {
            long int nei = (long int) VECTOR(*neis)[k];
            long int neimemb = (long int) VECTOR(*mymembership)[nei];
            if (neimemb == comm) {
                long int fi = (long int) VECTOR(*idx2)[nei];
                if (fi < size) {
                    to[j] += from[fi];
                }
                VECTOR(*tmp)[j] += 1;
            }
        }
    }

    /* Now calculate k^Tx/2m */
    ktx = 0.0; ktx2 = 0.0;
    for (j = 0; j < size + 1; j++) {
        long int oldid = (long int) VECTOR(*idx)[j];
        igraph_vector_int_t *neis = igraph_adjlist_get(adjlist, oldid);
        long int degree = igraph_vector_int_size(neis);
        if (j < size) {
            ktx += from[j] * degree;
        }
        ktx2 += degree;
    }
    ktx = ktx / no_of_edges / 2.0;
    ktx2 = ktx2 / no_of_edges / 2.0;

    /* Now calculate Bx */
    for (j = 0; j < size; j++) {
        long int oldid = (long int) VECTOR(*idx)[j];
        igraph_vector_int_t *neis = igraph_adjlist_get(adjlist, oldid);
        igraph_real_t degree = igraph_vector_int_size(neis);
        to[j] = to[j] - ktx * degree;
        VECTOR(*tmp)[j] = VECTOR(*tmp)[j] - ktx2 * degree;
    }

    /* -d_ij summa l in G B_il */
    for (j = 0; j < size; j++) {
        to[j] -= VECTOR(*tmp)[j] * from[j];
    }

    return 0;
}

static int igraph_i_community_leading_eigenvector_weighted(igraph_real_t *to,
                                                           const igraph_real_t *from,
                                                           int n, void *extra) {

    igraph_i_community_leading_eigenvector_data_t *data = extra;
    long int j, k, nlen, size = n;
    igraph_vector_t *idx = data->idx;
    igraph_vector_t *idx2 = data->idx2;
    igraph_vector_t *tmp = data->tmp;
    igraph_inclist_t *inclist = data->inclist;
    igraph_real_t ktx, ktx2;
    igraph_vector_t *mymembership = data->mymembership;
    long int comm = data->comm;
    const igraph_vector_t *weights = data->weights;
    const igraph_t *graph = data->graph;
    igraph_vector_t *strength = data->strength;
    igraph_real_t sw = data->sumweights;

    /* Ax */
    for (j = 0; j < size; j++) {
        long int oldid = (long int) VECTOR(*idx)[j];
        igraph_vector_int_t *inc = igraph_inclist_get(inclist, oldid);
        nlen = igraph_vector_int_size(inc);
        to[j] = 0.0;
        VECTOR(*tmp)[j] = 0.0;
        for (k = 0; k < nlen; k++) {
            long int edge = (long int) VECTOR(*inc)[k];
            igraph_real_t w = VECTOR(*weights)[edge];
            long int nei = IGRAPH_OTHER(graph, edge, oldid);
            long int neimemb = (long int) VECTOR(*mymembership)[nei];
            if (neimemb == comm) {
                to[j] += from[ (long int) VECTOR(*idx2)[nei] ] * w;
                VECTOR(*tmp)[j] += w;
            }
        }
    }

    /* k^Tx/2m */
    ktx = 0.0; ktx2 = 0.0;
    for (j = 0; j < size; j++) {
        long int oldid = (long int) VECTOR(*idx)[j];
        igraph_real_t str = VECTOR(*strength)[oldid];
        ktx += from[j] * str;
        ktx2 += str;
    }
    ktx = ktx / sw / 2.0;
    ktx2 = ktx2 / sw / 2.0;

    /* Bx */
    for (j = 0; j < size; j++) {
        long int oldid = (long int) VECTOR(*idx)[j];
        igraph_real_t str = VECTOR(*strength)[oldid];
        to[j] = to[j] - ktx * str;
        VECTOR(*tmp)[j] = VECTOR(*tmp)[j] - ktx2 * str;
    }

    /* -d_ij summa l in G B_il */
    for (j = 0; j < size; j++) {
        to[j] -= VECTOR(*tmp)[j] * from[j];
    }

    return 0;
}

static int igraph_i_community_leading_eigenvector2_weighted(igraph_real_t *to,
                                                            const igraph_real_t *from,
                                                            int n, void *extra) {

    igraph_i_community_leading_eigenvector_data_t *data = extra;
    long int j, k, nlen, size = n;
    igraph_vector_t *idx = data->idx;
    igraph_vector_t *idx2 = data->idx2;
    igraph_vector_t *tmp = data->tmp;
    igraph_inclist_t *inclist = data->inclist;
    igraph_real_t ktx, ktx2;
    igraph_vector_t *mymembership = data->mymembership;
    long int comm = data->comm;
    const igraph_vector_t *weights = data->weights;
    const igraph_t *graph = data->graph;
    igraph_vector_t *strength = data->strength;
    igraph_real_t sw = data->sumweights;

    /* Ax */
    for (j = 0; j < size; j++) {
        long int oldid = (long int) VECTOR(*idx)[j];
        igraph_vector_int_t *inc = igraph_inclist_get(inclist, oldid);
        nlen = igraph_vector_int_size(inc);
        to[j] = 0.0;
        VECTOR(*tmp)[j] = 0.0;
        for (k = 0; k < nlen; k++) {
            long int edge = (long int) VECTOR(*inc)[k];
            igraph_real_t w = VECTOR(*weights)[edge];
            long int nei = IGRAPH_OTHER(graph, edge, oldid);
            long int neimemb = (long int) VECTOR(*mymembership)[nei];
            if (neimemb == comm) {
                long int fi = (long int) VECTOR(*idx2)[nei];
                if (fi < size) {
                    to[j] += from[fi] * w;
                }
                VECTOR(*tmp)[j] += w;
            }
        }
    }

    /* k^Tx/2m */
    ktx = 0.0; ktx2 = 0.0;
    for (j = 0; j < size + 1; j++) {
        long int oldid = (long int) VECTOR(*idx)[j];
        igraph_real_t str = VECTOR(*strength)[oldid];
        if (j < size) {
            ktx += from[j] * str;
        }
        ktx2 += str;
    }
    ktx = ktx / sw / 2.0;
    ktx2 = ktx2 / sw / 2.0;

    /* Bx */
    for (j = 0; j < size; j++) {
        long int oldid = (long int) VECTOR(*idx)[j];
        igraph_real_t str = VECTOR(*strength)[oldid];
        to[j] = to[j] - ktx * str;
        VECTOR(*tmp)[j] = VECTOR(*tmp)[j] - ktx2 * str;
    }

    /* -d_ij summa l in G B_il */
    for (j = 0; j < size; j++) {
        to[j] -= VECTOR(*tmp)[j] * from[j];
    }

    return 0;
}

static void igraph_i_levc_free(igraph_vector_ptr_t *ptr) {
    long int i, n = igraph_vector_ptr_size(ptr);
    for (i = 0; i < n; i++) {
        igraph_vector_t *v = VECTOR(*ptr)[i];
        if (v) {
            igraph_vector_destroy(v);
            igraph_free(v);
        }
    }
}

static void igraph_i_error_handler_none(const char *reason, const char *file,
                                        int line, int igraph_errno) {
    IGRAPH_UNUSED(reason);
    IGRAPH_UNUSED(file);
    IGRAPH_UNUSED(line);
    IGRAPH_UNUSED(igraph_errno);
    /* do nothing */
}


/**
 * \ingroup communities
 * \function igraph_community_leading_eigenvector
 * \brief Leading eigenvector community finding (proper version).
 *
 * Newman's leading eigenvector method for detecting community
 * structure. This is the proper implementation of the recursive,
 * divisive algorithm: each split is done by maximizing the modularity
 * regarding the original network, see MEJ Newman: Finding community
 * structure in networks using the eigenvectors of matrices,
 * Phys Rev E 74:036104 (2006).
 *
 * \param graph The undirected input graph.
 * \param weights The weights of the edges, or a null pointer for
 *    unweighted graphs.
 * \param merges The result of the algorithm, a matrix containing the
 *    information about the splits performed. The matrix is built in
 *    the opposite way however, it is like the result of an
 *    agglomerative algorithm. If at the end of the algorithm (after
 *    \p steps steps was done) there are <quote>p</quote> communities,
 *    then these are numbered from zero to <quote>p-1</quote>. The
 *    first line of the matrix contains the first <quote>merge</quote>
 *    (which is in reality the last split) of two communities into
 *    community <quote>p</quote>, the merge in the second line forms
 *    community <quote>p+1</quote>, etc. The matrix should be
 *    initialized before calling and will be resized as needed.
 *    This argument is ignored of it is \c NULL.
 * \param membership The membership of the vertices after all the
 *    splits were performed will be stored here. The vector must be
 *    initialized  before calling and will be resized as needed.
 *    This argument is ignored if it is \c NULL. This argument can
 *    also be used to supply a starting configuration for the community
 *    finding, in the format of a membership vector. In this case the
 *    \p start argument must be set to 1.
 * \param steps The maximum number of steps to perform. It might
 *    happen that some component (or the whole network) has no
 *    underlying community structure and no further steps can be
 *    done. If you want as many steps as possible then supply the
 *    number of vertices in the network here.
 * \param options The options for ARPACK. \c n is always
 *    overwritten. \c ncv is set to at least 4.
 * \param modularity If not a null pointer, then it must be a pointer
 *    to a real number and the modularity score of the final division
 *    is stored here.
 * \param start Boolean, whether to use the community structure given
 *    in the \p membership argument as a starting point.
 * \param eigenvalues Pointer to an initialized vector or a null
 *    pointer. If not a null pointer, then the eigenvalues calculated
 *    along the community structure detection are stored here. The
 *    non-positive eigenvalues, that do not result a split, are stored
 *    as well.
 * \param eigenvectors If not a null pointer, then the eigenvectors
 *    that are calculated in each step of the algorithm, are stored here,
 *    in a pointer vector. Each eigenvector is stored in an
 *    \ref igraph_vector_t object. The user is responsible of
 *    deallocating the memory that belongs to the individual vectors,
 *    by calling first \ref igraph_vector_destroy(), and then
 *    \ref igraph_free() on them.
 * \param history Pointer to an initialized vector or a null pointer.
 *    If not a null pointer, then a trace of the algorithm is stored
 *    here, encoded numerically. The various operations:
 *    \clist
 *    \cli IGRAPH_LEVC_HIST_START_FULL
 *      Start the algorithm from an initial state where each connected
 *      component is a separate community.
 *    \cli IGRAPH_LEVC_HIST_START_GIVEN
 *      Start the algorithm from a given community structure. The next
 *      value in the vector contains the initial number of
 *      communities.
 *    \cli IGRAPH_LEVC_HIST_SPLIT
 *      Split a community into two communities. The id of the splitted
 *      community is given in the next element of the history vector.
 *      The id of the first new community is the same as the id of the
 *      splitted community. The id of the second community equals to
 *      the number of communities before the split.
 *    \cli IGRAPH_LEVC_HIST_FAILED
 *      Tried to split a community, but it was not worth it, as it
 *      does not result in a bigger modularity value. The id of the
 *      community is given in the next element of the vector.
 *    \endclist
 * \param callback A null pointer or a function of type \ref
 *    igraph_community_leading_eigenvector_callback_t. If given, this
 *    callback function is called after each eigenvector/eigenvalue
 *    calculation. If the callback returns a non-zero value, then the
 *    community finding algorithm stops. See the arguments passed to
 *    the callback at the documentation of \ref
 *    igraph_community_leading_eigenvector_callback_t.
 * \param callback_extra Extra argument to pass to the callback
 *    function.
 * \return Error code.
 *
 * \sa \ref igraph_community_walktrap() and \ref
 * igraph_community_spinglass() for other community structure
 * detection methods.
 *
 * Time complexity: O(|E|+|V|^2*steps), |V| is the number of vertices,
 * |E| the number of edges, <quote>steps</quote> the number of splits
 * performed.
 */

int igraph_community_leading_eigenvector(const igraph_t *graph,
        const igraph_vector_t *weights,
        igraph_matrix_t *merges,
        igraph_vector_t *membership,
        igraph_integer_t steps,
        igraph_arpack_options_t *options,
        igraph_real_t *modularity,
        igraph_bool_t start,
        igraph_vector_t *eigenvalues,
        igraph_vector_ptr_t *eigenvectors,
        igraph_vector_t *history,
        igraph_community_leading_eigenvector_callback_t *callback,
        void *callback_extra) {

    long int no_of_nodes = igraph_vcount(graph);
    long int no_of_edges = igraph_ecount(graph);
    igraph_dqueue_t tosplit;
    igraph_vector_t idx, idx2, mymerges;
    igraph_vector_t strength, tmp;
    long int staken = 0;
    igraph_adjlist_t adjlist;
    igraph_inclist_t inclist;
    long int i, j, k, l;
    long int communities;
    igraph_vector_t vmembership, *mymembership = membership;
    igraph_i_community_leading_eigenvector_data_t extra;
    igraph_arpack_storage_t storage;
    igraph_real_t mod = 0;
    igraph_arpack_function_t *arpcb1 =
        weights ? igraph_i_community_leading_eigenvector_weighted :
        igraph_i_community_leading_eigenvector;
    igraph_arpack_function_t *arpcb2 =
        weights ? igraph_i_community_leading_eigenvector2_weighted :
        igraph_i_community_leading_eigenvector2;
    igraph_real_t sumweights = 0.0;

    if (weights && no_of_edges != igraph_vector_size(weights)) {
        IGRAPH_ERROR("Invalid weight vector length", IGRAPH_EINVAL);
    }

    if (start && !membership) {
        IGRAPH_ERROR("Cannot start from given configuration if memberships "
                     "missing", IGRAPH_EINVAL);
    }

    if (start && membership &&
        igraph_vector_size(membership) != no_of_nodes) {
        IGRAPH_ERROR("Wrong length for vector of predefined memberships",
                     IGRAPH_EINVAL);
    }

    if (start && membership && igraph_vector_max(membership) >= no_of_nodes) {
        IGRAPH_WARNING("Too many communities in membership start vector");
    }

    if (igraph_is_directed(graph)) {
        IGRAPH_WARNING("This method was developed for undirected graphs");
    }

    if (steps < 0 || steps > no_of_nodes - 1) {
        steps = (igraph_integer_t) no_of_nodes - 1;
    }

    if (!membership) {
        mymembership = &vmembership;
        IGRAPH_VECTOR_INIT_FINALLY(mymembership, 0);
    }

    IGRAPH_VECTOR_INIT_FINALLY(&mymerges, 0);
    IGRAPH_CHECK(igraph_vector_reserve(&mymerges, steps * 2));
    IGRAPH_VECTOR_INIT_FINALLY(&idx, 0);
    if (eigenvalues)  {
        igraph_vector_clear(eigenvalues);
    }
    if (eigenvectors) {
        igraph_vector_ptr_clear(eigenvectors);
        IGRAPH_FINALLY(igraph_i_levc_free, eigenvectors);
    }

    IGRAPH_STATUS("Starting leading eigenvector method.\n", 0);

    if (!start) {
        /* Calculate the weakly connected components in the graph and use them as
         * an initial split */
        IGRAPH_CHECK(igraph_clusters(graph, mymembership, &idx, 0, IGRAPH_WEAK));
        communities = igraph_vector_size(&idx);
        IGRAPH_STATUSF(("Starting from %li component(s).\n", 0, communities));
        if (history) {
            IGRAPH_CHECK(igraph_vector_push_back(history,
                                                 IGRAPH_LEVC_HIST_START_FULL));
        }
    } else {
        /* Just create the idx vector for the given membership vector */
        communities = (long int) igraph_vector_max(mymembership) + 1;
        IGRAPH_STATUSF(("Starting from given membership vector with %li "
                        "communities.\n", 0, communities));
        if (history) {
            IGRAPH_CHECK(igraph_vector_push_back(history,
                                                 IGRAPH_LEVC_HIST_START_GIVEN));
            IGRAPH_CHECK(igraph_vector_push_back(history, communities));
        }
        IGRAPH_CHECK(igraph_vector_resize(&idx, communities));
        igraph_vector_null(&idx);
        for (i = 0; i < no_of_nodes; i++) {
            int t = (int) VECTOR(*mymembership)[i];
            VECTOR(idx)[t] += 1;
        }
    }

    IGRAPH_DQUEUE_INIT_FINALLY(&tosplit, 100);
    for (i = 0; i < communities; i++) {
        if (VECTOR(idx)[i] > 2) {
            igraph_dqueue_push(&tosplit, i);
        }
    }
    for (i = 1; i < communities; i++) {
        /* Record merge */
        IGRAPH_CHECK(igraph_vector_push_back(&mymerges, i - 1));
        IGRAPH_CHECK(igraph_vector_push_back(&mymerges, i));
        if (eigenvalues) {
            IGRAPH_CHECK(igraph_vector_push_back(eigenvalues, IGRAPH_NAN));
        }
        if (eigenvectors) {
            igraph_vector_t *v = igraph_Calloc(1, igraph_vector_t);
            if (!v) {
                IGRAPH_ERROR("Cannot do leading eigenvector community detection",
                             IGRAPH_ENOMEM);
            }
            IGRAPH_FINALLY(igraph_free, v);
            IGRAPH_VECTOR_INIT_FINALLY(v, 0);
            IGRAPH_CHECK(igraph_vector_ptr_push_back(eigenvectors, v));
            IGRAPH_FINALLY_CLEAN(2);
        }
        if (history) {
            IGRAPH_CHECK(igraph_vector_push_back(history, IGRAPH_LEVC_HIST_SPLIT));
            IGRAPH_CHECK(igraph_vector_push_back(history, i - 1));
        }
    }
    staken = communities - 1;

    IGRAPH_VECTOR_INIT_FINALLY(&tmp, no_of_nodes);
    IGRAPH_CHECK(igraph_vector_resize(&idx, no_of_nodes));
    igraph_vector_null(&idx);
    IGRAPH_VECTOR_INIT_FINALLY(&idx2, no_of_nodes);
    if (!weights) {
        IGRAPH_CHECK(igraph_adjlist_init(graph, &adjlist, IGRAPH_ALL));
        IGRAPH_FINALLY(igraph_adjlist_destroy, &adjlist);
    } else {
        IGRAPH_CHECK(igraph_inclist_init(graph, &inclist, IGRAPH_ALL));
        IGRAPH_FINALLY(igraph_inclist_destroy, &inclist);
        IGRAPH_VECTOR_INIT_FINALLY(&strength, no_of_nodes);
        IGRAPH_CHECK(igraph_strength(graph, &strength, igraph_vss_all(),
                                     IGRAPH_ALL, IGRAPH_LOOPS, weights));
        sumweights = igraph_vector_sum(weights);
    }

    options->ncv = 0;   /* 0 means "automatic" in igraph_arpack_rssolve */
    options->start = 0;
    options->which[0] = 'L'; options->which[1] = 'A';

    /* Memory for ARPACK */
    /* We are allocating memory for 20 eigenvectors since options->ncv won't be
     * larger than 20 when using automatic mode in igraph_arpack_rssolve */
    IGRAPH_CHECK(igraph_arpack_storage_init(&storage, (int) no_of_nodes, 20,
                                            (int) no_of_nodes, 1));
    IGRAPH_FINALLY(igraph_arpack_storage_destroy, &storage);
    extra.idx = &idx;
    extra.idx2 = &idx2;
    extra.tmp = &tmp;
    extra.adjlist = &adjlist;
    extra.inclist = &inclist;
    extra.weights = weights;
    extra.sumweights = sumweights;
    extra.graph = graph;
    extra.strength = &strength;
    extra.no_of_edges = no_of_edges;
    extra.mymembership = mymembership;

    while (!igraph_dqueue_empty(&tosplit) && staken < steps) {
        long int comm = (long int) igraph_dqueue_pop_back(&tosplit);
        /* depth first search */
        long int size = 0;
        igraph_real_t tmpev;

        IGRAPH_STATUSF(("Trying to split community %li... ", 0, comm));
        IGRAPH_ALLOW_INTERRUPTION();

        for (i = 0; i < no_of_nodes; i++) {
            if (VECTOR(*mymembership)[i] == comm) {
                VECTOR(idx)[size] = i;
                VECTOR(idx2)[i] = size++;
            }
        }

        staken++;
        if (size <= 2) {
            continue;
        }

        /* We solve two eigenproblems, one for the original modularity
           matrix, and one for the modularity matrix after deleting the
           last row and last column from it. This is a trick to find
           multiple leading eigenvalues, because ARPACK is sometimes
           unstable when the first two eigenvalues are requested, but it
           does much better for the single principal eigenvalue. */

        /* We start with the smaller eigenproblem. */

        options->n = (int) size - 1;
        options->info = 0;
        options->nev = 1;
        options->ldv = 0;
        options->ncv = 0;   /* 0 means "automatic" in igraph_arpack_rssolve */
        options->nconv = 0;
        options->lworkl = 0;        /* we surely have enough space */
        extra.comm = comm;

        /* We try calling the solver twice, once from a random starting
           point, once from a fixed one. This is because for some hard
           cases it tends to fail. We need to suppress error handling for
           the first call. */
        {
            int i;
            igraph_error_handler_t *errh =
                igraph_set_error_handler(igraph_i_error_handler_none);
            igraph_warning_handler_t *warnh =
                igraph_set_warning_handler(igraph_warning_handler_ignore);
            igraph_arpack_rssolve(arpcb2, &extra, options, &storage,
                                  /*values=*/ 0, /*vectors=*/ 0);
            igraph_set_error_handler(errh);
            igraph_set_warning_handler(warnh);
            if (options->nconv < 1) {
                /* Call again from a fixed starting point. Note that we cannot use a
                 * fixed all-1 starting vector as sometimes ARPACK would return a
                 * 'starting vector is zero' error -- this is of course not true but
                 * it's a result of ARPACK >= 3.6.3 trying to force the starting vector
                 * into the range of OP (i.e. the matrix being solved). The initial
                 * vector we use here seems to work, but I have no theoretical argument
                 * for its usage; it just happens to work. */
                options->start = 1;
                options->info = 0;
                options->ncv = 0;
                options->lworkl = 0;    /* we surely have enough space */
                for (i = 0; i < options->n ; i++) {
                    storage.resid[i] = i % 2 ? 1 : -1;
                }
                IGRAPH_CHECK(igraph_arpack_rssolve(arpcb2, &extra, options, &storage,
                                                   /*values=*/ 0, /*vectors=*/ 0));
                options->start = 0;
            }
        }

        if (options->nconv < 1) {
            IGRAPH_ERROR("ARPACK did not converge", IGRAPH_ARPACK_FAILED);
        }

        tmpev = storage.d[0];

        /* Now we do the original eigenproblem, again, twice if needed */

        options->n = (int) size;
        options->info = 0;
        options->nev = 1;
        options->ldv = 0;
        options->nconv = 0;
        options->lworkl = 0;    /* we surely have enough space */
        options->ncv = 0;   /* 0 means "automatic" in igraph_arpack_rssolve */

        {
            int i;
            igraph_error_handler_t *errh =
                igraph_set_error_handler(igraph_i_error_handler_none);
            igraph_arpack_rssolve(arpcb1, &extra, options, &storage,
                                  /*values=*/ 0, /*vectors=*/ 0);
            igraph_set_error_handler(errh);
            if (options->nconv < 1) {
                /* Call again from a fixed starting point. See the comment a few lines
                 * above about the exact choice of this starting vector */
                options->start = 1;
                options->info = 0;
                options->ncv = 0;
                options->lworkl = 0;    /* we surely have enough space */
                for (i = 0; i < options->n; i++) {
                    storage.resid[i] = i % 2 ? 1 : -1;
                }
                IGRAPH_CHECK(igraph_arpack_rssolve(arpcb1, &extra, options, &storage,
                                                   /*values=*/ 0, /*vectors=*/ 0));
                options->start = 0;
            }
        }

        if (options->nconv < 1) {
            IGRAPH_ERROR("ARPACK did not converge", IGRAPH_ARPACK_FAILED);
        }

        /* Ok, we have the leading eigenvector of the modularity matrix*/

        /* ---------------------------------------------------------------*/
        /* To avoid numeric errors */
        if (fabs(storage.d[0]) < 1e-8) {
            storage.d[0] = 0;
        }

        /* We replace very small (in absolute value) elements of the
           leading eigenvector with zero, to get the same result,
           consistently.*/
        for (i = 0; i < size; i++) {
            if (fabs(storage.v[i]) < 1e-8) {
                storage.v[i] = 0;
            }
        }

        /* Just to have the always the same result, we multiply by -1
           if the first (nonzero) element is not positive. */
        for (i = 0; i < size; i++) {
            if (storage.v[i] != 0) {
                break;
            }
        }
        if (i < size && storage.v[i] < 0) {
            for (i = 0; i < size; i++) {
                storage.v[i] = - storage.v[i];
            }
        }
        /* ---------------------------------------------------------------*/

        if (callback) {
            igraph_vector_t vv;
            int ret;
            igraph_vector_view(&vv, storage.v, size);
            ret = callback(mymembership, comm, storage.d[0], &vv,
                           arpcb1, &extra, callback_extra);
            if (ret) {
                break;
            }
        }

        if (eigenvalues) {
            IGRAPH_CHECK(igraph_vector_push_back(eigenvalues, storage.d[0]));
        }

        if (eigenvectors) {
            igraph_vector_t *v = igraph_Calloc(1, igraph_vector_t);
            if (!v) {
                IGRAPH_ERROR("Cannot do leading eigenvector community detection",
                             IGRAPH_ENOMEM);
            }
            IGRAPH_FINALLY(igraph_free, v);
            IGRAPH_VECTOR_INIT_FINALLY(v, size);
            for (i = 0; i < size; i++) {
                VECTOR(*v)[i] = storage.v[i];
            }
            IGRAPH_CHECK(igraph_vector_ptr_push_back(eigenvectors, v));
            IGRAPH_FINALLY_CLEAN(2);
        }

        if (storage.d[0] <= 0) {
            IGRAPH_STATUS("no split.\n", 0);
            if (history) {
                IGRAPH_CHECK(igraph_vector_push_back(history,
                                                     IGRAPH_LEVC_HIST_FAILED));
                IGRAPH_CHECK(igraph_vector_push_back(history, comm));
            }
            continue;
        }

        /* Check for multiple leading eigenvalues */

        if (fabs(storage.d[0] - tmpev) < 1e-8) {
            IGRAPH_STATUS("multiple principal eigenvalue, no split.\n", 0);
            if (history) {
                IGRAPH_CHECK(igraph_vector_push_back(history,
                                                     IGRAPH_LEVC_HIST_FAILED));
                IGRAPH_CHECK(igraph_vector_push_back(history, comm));
            }
            continue;
        }

        /* Count the number of vertices in each community after the split */
        l = 0;
        for (j = 0; j < size; j++) {
            if (storage.v[j] < 0) {
                storage.v[j] = -1;
                l++;
            } else {
                storage.v[j] = 1;
            }
        }
        if (l == 0 || l == size) {
            IGRAPH_STATUS("no split.\n", 0);
            if (history) {
                IGRAPH_CHECK(igraph_vector_push_back(history,
                                                     IGRAPH_LEVC_HIST_FAILED));
                IGRAPH_CHECK(igraph_vector_push_back(history, comm));
            }
            continue;
        }

        /* Check that Q increases with our choice of split */
        arpcb1(storage.v + size, storage.v, (int) size, &extra);
        mod = 0;
        for (i = 0; i < size; i++) {
            mod += storage.v[size + i] * storage.v[i];
        }
        if (mod <= 1e-8) {
            IGRAPH_STATUS("no modularity increase, no split.\n", 0);
            if (history) {
                IGRAPH_CHECK(igraph_vector_push_back(history,
                                                     IGRAPH_LEVC_HIST_FAILED));
                IGRAPH_CHECK(igraph_vector_push_back(history, comm));
            }
            continue;
        }

        communities++;
        IGRAPH_STATUS("split.\n", 0);

        /* Rewrite the mymembership vector */
        for (j = 0; j < size; j++) {
            if (storage.v[j] < 0) {
                long int oldid = (long int) VECTOR(idx)[j];
                VECTOR(*mymembership)[oldid] = communities - 1;
            }
        }

        /* Record merge */
        IGRAPH_CHECK(igraph_vector_push_back(&mymerges, comm));
        IGRAPH_CHECK(igraph_vector_push_back(&mymerges, communities - 1));
        if (history) {
            IGRAPH_CHECK(igraph_vector_push_back(history, IGRAPH_LEVC_HIST_SPLIT));
            IGRAPH_CHECK(igraph_vector_push_back(history, comm));
        }

        /* Store the resulting communities in the queue if needed */
        if (l > 1) {
            IGRAPH_CHECK(igraph_dqueue_push(&tosplit, communities - 1));
        }
        if (size - l > 1) {
            IGRAPH_CHECK(igraph_dqueue_push(&tosplit, comm));
        }

    }

    igraph_arpack_storage_destroy(&storage);
    IGRAPH_FINALLY_CLEAN(1);
    if (!weights) {
        igraph_adjlist_destroy(&adjlist);
        IGRAPH_FINALLY_CLEAN(1);
    } else {
        igraph_inclist_destroy(&inclist);
        igraph_vector_destroy(&strength);
        IGRAPH_FINALLY_CLEAN(2);
    }
    igraph_dqueue_destroy(&tosplit);
    igraph_vector_destroy(&tmp);
    igraph_vector_destroy(&idx2);
    IGRAPH_FINALLY_CLEAN(3);

    IGRAPH_STATUS("Done.\n", 0);

    /* reform the mymerges vector */
    if (merges) {
        igraph_vector_null(&idx);
        l = igraph_vector_size(&mymerges);
        k = communities;
        j = 0;
        IGRAPH_CHECK(igraph_matrix_resize(merges, l / 2, 2));
        for (i = l; i > 0; i -= 2) {
            long int from = (long int) VECTOR(mymerges)[i - 1];
            long int to = (long int) VECTOR(mymerges)[i - 2];
            MATRIX(*merges, j, 0) = VECTOR(mymerges)[i - 2];
            MATRIX(*merges, j, 1) = VECTOR(mymerges)[i - 1];
            if (VECTOR(idx)[from] != 0) {
                MATRIX(*merges, j, 1) = VECTOR(idx)[from] - 1;
            }
            if (VECTOR(idx)[to] != 0) {
                MATRIX(*merges, j, 0) = VECTOR(idx)[to] - 1;
            }
            VECTOR(idx)[to] = ++k;
            j++;
        }
    }

    if (eigenvectors) {
        IGRAPH_FINALLY_CLEAN(1);
    }
    igraph_vector_destroy(&idx);
    igraph_vector_destroy(&mymerges);
    IGRAPH_FINALLY_CLEAN(2);

    if (modularity) {
        IGRAPH_CHECK(igraph_modularity(graph, mymembership, modularity,
                                       weights));
    }

    if (!membership) {
        igraph_vector_destroy(mymembership);
        IGRAPH_FINALLY_CLEAN(1);
    }

    return 0;
}

/**
 * \function igraph_le_community_to_membership
 * Vertex membership from the leading eigenvector community structure
 *
 * This function creates a membership vector from the
 * result of \ref igraph_community_leading_eigenvector(),
 * It takes \c membership
 * and performs \c steps merges, according to the supplied
 * \c merges matrix.
 * \param merges The matrix defining the merges to make.
 *     This is usually from the output of the leading eigenvector community
 *     structure detection routines.
 * \param steps The number of steps to make according to \c merges.
 * \param membership Initially the starting membership vector,
 *     on output the resulting membership vector, after performing \c steps merges.
 * \param csize Optionally the sizes of the communities is stored here,
 *     if this is not a null pointer, but an initialized vector.
 * \return Error code.
 *
 * Time complexity: O(|V|), the number of vertices.
 */

int igraph_le_community_to_membership(const igraph_matrix_t *merges,
                                      igraph_integer_t steps,
                                      igraph_vector_t *membership,
                                      igraph_vector_t *csize) {

    long int no_of_nodes = igraph_vector_size(membership);
    igraph_vector_t fake_memb;
    long int components, i;

    if (igraph_matrix_nrow(merges) < steps) {
        IGRAPH_ERROR("`steps' to big or `merges' matrix too short", IGRAPH_EINVAL);
    }

    components = (long int) igraph_vector_max(membership) + 1;
    if (components > no_of_nodes) {
        IGRAPH_ERROR("Invalid membership vector, too many components", IGRAPH_EINVAL);
    }
    if (steps >= components) {
        IGRAPH_ERROR("Cannot make `steps' steps from supplied membership vector",
                     IGRAPH_EINVAL);
    }

    IGRAPH_VECTOR_INIT_FINALLY(&fake_memb, components);

    /* Check membership vector */
    for (i = 0; i < no_of_nodes; i++) {
        if (VECTOR(*membership)[i] < 0) {
            IGRAPH_ERROR("Invalid membership vector, negative id", IGRAPH_EINVAL);
        }
        VECTOR(fake_memb)[ (long int) VECTOR(*membership)[i] ] += 1;
    }
    for (i = 0; i < components; i++) {
        if (VECTOR(fake_memb)[i] == 0) {
            IGRAPH_ERROR("Invalid membership vector, empty cluster", IGRAPH_EINVAL);
        }
    }

    IGRAPH_CHECK(igraph_community_to_membership(merges, (igraph_integer_t)
                 components, steps,
                 &fake_memb, 0));

    /* Ok, now we have the membership of the initial components,
       rewrite the original membership vector. */

    if (csize) {
        IGRAPH_CHECK(igraph_vector_resize(csize, components - steps));
        igraph_vector_null(csize);
    }

    for (i = 0; i < no_of_nodes; i++) {
        VECTOR(*membership)[i] = VECTOR(fake_memb)[ (long int) VECTOR(*membership)[i] ];
        if (csize) {
            VECTOR(*csize)[ (long int) VECTOR(*membership)[i] ] += 1;
        }
    }

    igraph_vector_destroy(&fake_memb);
    IGRAPH_FINALLY_CLEAN(1);
    return 0;
}

/********************************************************************/

/**
 * \ingroup communities
 * \function igraph_community_fluid_communities
 * \brief Community detection algorithm based on the simple idea of
 * several fluids interacting in a non-homogeneous environment
 * (the graph topology), expanding and contracting based on their
 * interaction and density.
 *
 * This function implements the community detection method described in:
 * Parés F, Gasulla DG, et. al. (2018) Fluid Communities: A Competitive,
 * Scalable and Diverse Community Detection Algorithm. In: Complex Networks
 * &amp; Their Applications VI: Proceedings of Complex Networks 2017 (The Sixth
 * International Conference on Complex Networks and Their Applications),
 * Springer, vol 689, p 229.
 *
 * \param graph The input graph. The graph must be simple and connected.
 *   Empty graphs are not supported as well as single vertex graphs.
 *   Edge directions are ignored. Weights are not considered.
 * \param no_of_communities The number of communities to be found. Must be
 *   greater than 0 and fewer than number of vertices in the graph.
 * \param membership The result vector mapping vertices to the communities
 * they are assigned to.
 * \param modularity If not a null pointer, then it must be a pointer
 *   to a real number. The modularity score of the detected community
 *   structure is stored here.
 * \return Error code.
 *
 * Time complexity: O(|E|)
 *
 * \example examples/tests/igraph_community_fluid_communities.c
 */
int igraph_community_fluid_communities(const igraph_t *graph,
                                       igraph_integer_t no_of_communities,
                                       igraph_vector_t *membership,
                                       igraph_real_t *modularity) {
    /* Declaration of variables */
    long int no_of_nodes, i, j, k, kv1;
    igraph_adjlist_t al;
    double max_density;
    igraph_bool_t res, running;
    igraph_vector_t node_order, density, label_counters, dominant_labels, nonzero_labels;
    igraph_vector_int_t com_to_numvertices;

    /* Initialization of variables needed for initial checking */
    no_of_nodes = igraph_vcount(graph);

    /* Checking input values */
    if (no_of_nodes < 2) {
        IGRAPH_ERROR("Empty and single vertex graphs are not supported.", IGRAPH_EINVAL);
    }
    if ((long int) no_of_communities < 1) {
        IGRAPH_ERROR("'no_of_communities' must be greater than 0.", IGRAPH_EINVAL);
    }
    if ((long int) no_of_communities > no_of_nodes) {
        IGRAPH_ERROR("'no_of_communities' can not be greater than number of nodes in "
                     "the graph.", IGRAPH_EINVAL);
    }
    igraph_is_simple(graph, &res);
    if (!res) {
        IGRAPH_ERROR("Only simple graphs are supported.", IGRAPH_EINVAL);
    }
    igraph_is_connected(graph, &res, IGRAPH_WEAK);
    if (!res) {
        IGRAPH_ERROR("Disconnected graphs are not supported.", IGRAPH_EINVAL);
    }
    if (igraph_is_directed(graph)) {
        IGRAPH_WARNING("Edge directions are ignored.");
    }

    /* Internal variables initialization */
    max_density = 1.0;
    running = 1;

    /* Resize membership vector (number of nodes) */
    IGRAPH_CHECK(igraph_vector_resize(membership, no_of_nodes));

    /* Initialize density and com_to_numvertices vectors */
    IGRAPH_CHECK(igraph_vector_init(&density, (long int) no_of_communities));
    IGRAPH_FINALLY(igraph_vector_destroy, &density);
    IGRAPH_CHECK(igraph_vector_int_init(&com_to_numvertices, (long int) no_of_communities));
    IGRAPH_FINALLY(igraph_vector_int_destroy, &com_to_numvertices);

    /* Initialize node ordering vector */
    IGRAPH_CHECK(igraph_vector_init_seq(&node_order, 0, no_of_nodes - 1));
    IGRAPH_FINALLY(igraph_vector_destroy, &node_order);

    /* Initialize the membership vector with 0 values */
    igraph_vector_null(membership);
    /* Initialize densities to max_density */
    igraph_vector_fill(&density, max_density);

    RNG_BEGIN();

    /* Initialize com_to_numvertices and initialize communities into membership vector */
    IGRAPH_CHECK(igraph_vector_shuffle(&node_order));
    for (i = 0; i < no_of_communities; i++) {
        /* Initialize membership at initial nodes for each community
         * where 0 refers to have no label*/
        VECTOR(*membership)[(long int)VECTOR(node_order)[i]] = i + 1.0;
        /* Initialize com_to_numvertices list: Number of vertices for each community */
        VECTOR(com_to_numvertices)[i] = 1;
    }

    /* Create an adjacency list representation for efficiency. */
    IGRAPH_CHECK(igraph_adjlist_init(graph, &al, IGRAPH_ALL));
    IGRAPH_FINALLY(igraph_adjlist_destroy, &al);

    /* Create storage space for counting distinct labels and dominant ones */
    IGRAPH_VECTOR_INIT_FINALLY(&dominant_labels, (long int) no_of_communities);
    IGRAPH_VECTOR_INIT_FINALLY(&nonzero_labels, (long int) no_of_communities);

    IGRAPH_CHECK(igraph_vector_init(&label_counters, (long int) no_of_communities));
    IGRAPH_FINALLY(igraph_vector_destroy, &label_counters);

    /* running is the convergence boolean variable */
    running = 1;
    while (running) {
        /* Declarations of varibales used inside main loop */
        long int v1, size, rand_idx;
        igraph_real_t max_count, label_counter_diff;
        igraph_vector_int_t *neis;
        igraph_bool_t same_label_in_dominant;

        running = 0;

        /* Shuffle the node ordering vector */
        IGRAPH_CHECK(igraph_vector_shuffle(&node_order));
        /* In the prescribed order, loop over the vertices and reassign labels */
        for (i = 0; i < no_of_nodes; i++) {
            /* Clear dominant_labels and nonzero_labels vectors */
            igraph_vector_clear(&dominant_labels);
            igraph_vector_null(&label_counters);

            /* Obtain actual node index */
            v1 = (long int) VECTOR(node_order)[i];
            /* Take into account same label in updating rule */
            kv1 = (long int) VECTOR(*membership)[v1];
            max_count = 0.0;
            if (kv1 != 0) {
                VECTOR(label_counters)[kv1 - 1] += VECTOR(density)[kv1 - 1];
                /* Set up max_count */
                max_count = VECTOR(density)[kv1 - 1];
                /* Initialize dominant_labels */
                IGRAPH_CHECK(igraph_vector_resize(&dominant_labels, 1));
                VECTOR(dominant_labels)[0] = kv1;
            }

            /* Count the weights corresponding to different labels */
            neis = igraph_adjlist_get(&al, v1);
            size = igraph_vector_int_size(neis);
            for (j = 0; j < size; j++) {
                k = (long int) VECTOR(*membership)[(long)VECTOR(*neis)[j]];
                /* skip if it has no label yet */
                if (k == 0) {
                    continue;
                }
                /* Update label counter and evaluate diff against max_count*/
                VECTOR(label_counters)[k - 1] += VECTOR(density)[k - 1];
                label_counter_diff = VECTOR(label_counters)[k - 1] - max_count;
                /* Check if this label must be included in dominant_labels vector */
                if (label_counter_diff > 0.0001) {
                    max_count = VECTOR(label_counters)[k - 1];
                    IGRAPH_CHECK(igraph_vector_resize(&dominant_labels, 1));
                    VECTOR(dominant_labels)[0] = k;
                } else if (-0.0001 < label_counter_diff && label_counter_diff < 0.0001) {
                    IGRAPH_CHECK(igraph_vector_push_back(&dominant_labels, k));
                }
            }

            if (!igraph_vector_empty(&dominant_labels)) {
                /* Maintain same label if it exists in dominant_labels */
                same_label_in_dominant = igraph_vector_contains(&dominant_labels, kv1);

                if (!same_label_in_dominant) {
                    /* We need at least one more iteration */
                    running = 1;

                    /* Select randomly from the dominant labels */
                    rand_idx = RNG_INTEGER(0, igraph_vector_size(&dominant_labels) - 1);
                    k = (long int) VECTOR(dominant_labels)[rand_idx];

                    if (kv1 != 0) {
                        /* Subtract 1 vertex from corresponding community in com_to_numvertices */
                        VECTOR(com_to_numvertices)[kv1 - 1] -= 1;
                        /* Re-calculate density for community kv1 */
                        VECTOR(density)[kv1 - 1] = max_density / VECTOR(com_to_numvertices)[kv1 - 1];
                    }

                    /* Update vertex new label */
                    VECTOR(*membership)[v1] = k;

                    /* Add 1 vertex to corresponding new community in com_to_numvertices */
                    VECTOR(com_to_numvertices)[k - 1] += 1;
                    /* Re-calculate density for new community k */
                    VECTOR(density)[k - 1] = max_density / VECTOR(com_to_numvertices)[k - 1];
                }
            }
        }
    }

    RNG_END();


    /* Shift back the membership vector */
    /* There must be no 0 labels in membership vector at this point */
    for (i = 0; i < no_of_nodes; i++) {
        VECTOR(*membership)[i] -= 1;
        /* Something went wrong: At least one vertex has no community assigned */
        if (VECTOR(*membership)[i] < 0) {
            IGRAPH_ERROR("Something went wrong during execution. One or more vertices got "
                         "no community assigned at algorithm convergence.", IGRAPH_EINTERNAL);
        }
    }

    igraph_adjlist_destroy(&al);
    IGRAPH_FINALLY_CLEAN(1);

    if (modularity) {
        IGRAPH_CHECK(igraph_modularity(graph, membership, modularity,
                                       NULL));
    }

    igraph_vector_destroy(&node_order);
    igraph_vector_destroy(&density);
    igraph_vector_int_destroy(&com_to_numvertices);
    igraph_vector_destroy(&label_counters);
    igraph_vector_destroy(&dominant_labels);
    igraph_vector_destroy(&nonzero_labels);
    IGRAPH_FINALLY_CLEAN(6);

    return 0;
}

/********************************************************************/

/**
 * \ingroup communities
 * \function igraph_community_label_propagation
 * \brief Community detection based on label propagation
 *
 * This function implements the community detection method described in:
 * Raghavan, U.N. and Albert, R. and Kumara, S.: Near linear time algorithm
 * to detect community structures in large-scale networks. Phys Rev E
 * 76, 036106. (2007). This version extends the original method by
 * the ability to take edge weights into consideration and also
 * by allowing some labels to be fixed.
 *
 * </para><para>
 * Weights are taken into account as follows: when the new label of node
 * i is determined, the algorithm iterates over all edges incident on
 * node i and calculate the total weight of edges leading to other
 * nodes with label 0, 1, 2, ..., k-1 (where k is the number of possible
 * labels). The new label of node i will then be the label whose edges
 * (among the ones incident on node i) have the highest total weight.
 *
 * \param graph The input graph, should be undirected to make sense.
 * \param membership The membership vector, the result is returned here.
 *    For each vertex it gives the ID of its community (label).
 * \param weights The weight vector, it should contain a positive
 *    weight for all the edges.
 * \param initial The initial state. If NULL, every vertex will have
 *   a different label at the beginning. Otherwise it must be a vector
 *   with an entry for each vertex. Non-negative values denote different
 *   labels, negative entries denote vertices without labels.
 * \param fixed Boolean vector denoting which labels are fixed. Of course
 *   this makes sense only if you provided an initial state, otherwise
 *   this element will be ignored. Also note that vertices without labels
 *   cannot be fixed.
 * \param modularity If not a null pointer, then it must be a pointer
 *   to a real number. The modularity score of the detected community
 *   structure is stored here.
 * \return Error code.
 *
 * Time complexity: O(m+n)
 *
 * \example examples/simple/igraph_community_label_propagation.c
 */
int igraph_community_label_propagation(const igraph_t *graph,
                                       igraph_vector_t *membership,
                                       const igraph_vector_t *weights,
                                       const igraph_vector_t *initial,
                                       igraph_vector_bool_t *fixed,
                                       igraph_real_t *modularity) {
    long int no_of_nodes = igraph_vcount(graph);
    long int no_of_edges = igraph_ecount(graph);
    long int no_of_not_fixed_nodes = no_of_nodes;
    long int i, j, k;
    igraph_adjlist_t al;
    igraph_inclist_t il;
    igraph_bool_t running = 1;

    igraph_vector_t label_counters, dominant_labels, nonzero_labels, node_order;

    /* The implementation uses a trick to avoid negative array indexing:
     * elements of the membership vector are increased by 1 at the start
     * of the algorithm; this to allow us to denote unlabeled vertices
     * (if any) by zeroes. The membership vector is shifted back in the end
     */

    /* Do some initial checks */
    if (fixed && igraph_vector_bool_size(fixed) != no_of_nodes) {
        IGRAPH_ERROR("Invalid fixed labeling vector length", IGRAPH_EINVAL);
    }
    if (weights) {
        if (igraph_vector_size(weights) != no_of_edges) {
            IGRAPH_ERROR("Invalid weight vector length", IGRAPH_EINVAL);
        } else if (igraph_vector_min(weights) < 0) {
            IGRAPH_ERROR("Weights must be non-negative", IGRAPH_EINVAL);
        }
    }
    if (fixed && !initial) {
        IGRAPH_WARNING("Ignoring fixed vertices as no initial labeling given");
    }

    IGRAPH_CHECK(igraph_vector_resize(membership, no_of_nodes));

    if (initial) {
        if (igraph_vector_size(initial) != no_of_nodes) {
            IGRAPH_ERROR("Invalid initial labeling vector length", IGRAPH_EINVAL);
        }
        /* Check if the labels used are valid, initialize membership vector */
        for (i = 0; i < no_of_nodes; i++) {
            if (VECTOR(*initial)[i] < 0) {
                VECTOR(*membership)[i] = 0;
            } else {
                VECTOR(*membership)[i] = floor(VECTOR(*initial)[i]) + 1;
            }
        }
        if (fixed) {
            for (i = 0; i < no_of_nodes; i++) {
                if (VECTOR(*fixed)[i]) {
                    if (VECTOR(*membership)[i] == 0) {
                        IGRAPH_WARNING("Fixed nodes cannot be unlabeled, ignoring them");
                        VECTOR(*fixed)[i] = 0;
                    } else {
                        no_of_not_fixed_nodes--;
                    }
                }
            }
        }

        i = (long int) igraph_vector_max(membership);
        if (i > no_of_nodes) {
            IGRAPH_ERROR("elements of the initial labeling vector must be between 0 and |V|-1", IGRAPH_EINVAL);
        }
        if (i <= 0) {
            IGRAPH_ERROR("at least one vertex must be labeled in the initial labeling", IGRAPH_EINVAL);
        }
    } else {
        for (i = 0; i < no_of_nodes; i++) {
            VECTOR(*membership)[i] = i + 1;
        }
    }

    /* Create an adjacency/incidence list representation for efficiency.
     * For the unweighted case, the adjacency list is enough. For the
     * weighted case, we need the incidence list */
    if (weights) {
        IGRAPH_CHECK(igraph_inclist_init(graph, &il, IGRAPH_IN));
        IGRAPH_FINALLY(igraph_inclist_destroy, &il);
    } else {
        IGRAPH_CHECK(igraph_adjlist_init(graph, &al, IGRAPH_IN));
        IGRAPH_FINALLY(igraph_adjlist_destroy, &al);
    }

    /* Create storage space for counting distinct labels and dominant ones */
    IGRAPH_VECTOR_INIT_FINALLY(&label_counters, no_of_nodes + 1);
    IGRAPH_VECTOR_INIT_FINALLY(&dominant_labels, 0);
    IGRAPH_VECTOR_INIT_FINALLY(&nonzero_labels, 0);
    IGRAPH_CHECK(igraph_vector_reserve(&dominant_labels, 2));

    RNG_BEGIN();

    /* Initialize node ordering vector with only the not fixed nodes */
    if (fixed) {
        IGRAPH_VECTOR_INIT_FINALLY(&node_order, no_of_not_fixed_nodes);
        for (i = 0, j = 0; i < no_of_nodes; i++) {
            if (!VECTOR(*fixed)[i]) {
                VECTOR(node_order)[j] = i;
                j++;
            }
        }
    } else {
        IGRAPH_CHECK(igraph_vector_init_seq(&node_order, 0, no_of_nodes - 1));
        IGRAPH_FINALLY(igraph_vector_destroy, &node_order);
    }

    running = 1;
    while (running) {
        long int v1, num_neis;
        igraph_real_t max_count;
        igraph_vector_int_t *neis;
        igraph_vector_int_t *ineis;
        igraph_bool_t was_zero;

        running = 0;

        /* Shuffle the node ordering vector */
        IGRAPH_CHECK(igraph_vector_shuffle(&node_order));
        /* In the prescribed order, loop over the vertices and reassign labels */
        for (i = 0; i < no_of_not_fixed_nodes; i++) {
            v1 = (long int) VECTOR(node_order)[i];

            /* Count the weights corresponding to different labels */
            igraph_vector_clear(&dominant_labels);
            igraph_vector_clear(&nonzero_labels);
            max_count = 0.0;
            if (weights) {
                ineis = igraph_inclist_get(&il, v1);
                num_neis = igraph_vector_int_size(ineis);
                for (j = 0; j < num_neis; j++) {
                    k = (long int) VECTOR(*membership)[
                    (long)IGRAPH_OTHER(graph, VECTOR(*ineis)[j], v1) ];
                    if (k == 0) {
                        continue;    /* skip if it has no label yet */
                    }
                    was_zero = (VECTOR(label_counters)[k] == 0);
                    VECTOR(label_counters)[k] += VECTOR(*weights)[(long)VECTOR(*ineis)[j]];
                    if (was_zero && VECTOR(label_counters)[k] != 0) {
                        /* counter just became nonzero */
                        IGRAPH_CHECK(igraph_vector_push_back(&nonzero_labels, k));
                    }
                    if (max_count < VECTOR(label_counters)[k]) {
                        max_count = VECTOR(label_counters)[k];
                        IGRAPH_CHECK(igraph_vector_resize(&dominant_labels, 1));
                        VECTOR(dominant_labels)[0] = k;
                    } else if (max_count == VECTOR(label_counters)[k]) {
                        IGRAPH_CHECK(igraph_vector_push_back(&dominant_labels, k));
                    }
                }
            } else {
                neis = igraph_adjlist_get(&al, v1);
                num_neis = igraph_vector_int_size(neis);
                for (j = 0; j < num_neis; j++) {
                    k = (long int) VECTOR(*membership)[(long)VECTOR(*neis)[j]];
                    if (k == 0) {
                        continue;    /* skip if it has no label yet */
                    }
                    VECTOR(label_counters)[k]++;
                    if (VECTOR(label_counters)[k] == 1) {
                        /* counter just became nonzero */
                        IGRAPH_CHECK(igraph_vector_push_back(&nonzero_labels, k));
                    }
                    if (max_count < VECTOR(label_counters)[k]) {
                        max_count = VECTOR(label_counters)[k];
                        IGRAPH_CHECK(igraph_vector_resize(&dominant_labels, 1));
                        VECTOR(dominant_labels)[0] = k;
                    } else if (max_count == VECTOR(label_counters)[k]) {
                        IGRAPH_CHECK(igraph_vector_push_back(&dominant_labels, k));
                    }
                }
            }

            if (igraph_vector_size(&dominant_labels) > 0) {
                /* Select randomly from the dominant labels */
                k = RNG_INTEGER(0, igraph_vector_size(&dominant_labels) - 1);
                k = (long int) VECTOR(dominant_labels)[k];
                /* Check if the _current_ label of the node is also dominant */
                if (VECTOR(label_counters)[(long)VECTOR(*membership)[v1]] != max_count) {
                    /* Nope, we need at least one more iteration */
                    running = 1;
                }
                VECTOR(*membership)[v1] = k;
            }

            /* Clear the nonzero elements in label_counters */
            num_neis = igraph_vector_size(&nonzero_labels);
            for (j = 0; j < num_neis; j++) {
                VECTOR(label_counters)[(long int)VECTOR(nonzero_labels)[j]] = 0;
            }
        }
    }

    RNG_END();

    /* Shift back the membership vector, permute labels in increasing order */
    /* We recycle label_counters here :) */
    igraph_vector_fill(&label_counters, -1);
    j = 0;
    for (i = 0; i < no_of_nodes; i++) {
        k = (long)VECTOR(*membership)[i] - 1;
        if (k >= 0) {
            if (VECTOR(label_counters)[k] == -1) {
                /* We have seen this label for the first time */
                VECTOR(label_counters)[k] = j;
                k = j;
                j++;
            } else {
                k = (long int) VECTOR(label_counters)[k];
            }
        } else {
            /* This is an unlabeled vertex */
        }
        VECTOR(*membership)[i] = k;
    }

    if (weights) {
        igraph_inclist_destroy(&il);
    } else {
        igraph_adjlist_destroy(&al);
    }
    IGRAPH_FINALLY_CLEAN(1);

    if (modularity) {
        IGRAPH_CHECK(igraph_modularity(graph, membership, modularity,
                                       weights));
    }

    igraph_vector_destroy(&node_order);
    igraph_vector_destroy(&label_counters);
    igraph_vector_destroy(&dominant_labels);
    igraph_vector_destroy(&nonzero_labels);
    IGRAPH_FINALLY_CLEAN(4);

    return 0;
}

/********************************************************************/

/* Structure storing a community */
typedef struct {
    igraph_integer_t size;           /* Size of the community */
    igraph_real_t weight_inside;     /* Sum of edge weights inside community */
    igraph_real_t weight_all;        /* Sum of edge weights starting/ending
                                      in the community */
} igraph_i_multilevel_community;

/* Global community list structure */
typedef struct {
    long int communities_no, vertices_no;  /* Number of communities, number of vertices */
    igraph_real_t weight_sum;              /* Sum of edges weight in the whole graph */
    igraph_i_multilevel_community *item;   /* List of communities */
    igraph_vector_t *membership;           /* Community IDs */
    igraph_vector_t *weights;        /* Graph edge weights */
} igraph_i_multilevel_community_list;

/* Computes the modularity of a community partitioning */
static igraph_real_t igraph_i_multilevel_community_modularity(
    const igraph_i_multilevel_community_list *communities) {
    igraph_real_t result = 0;
    long int i;
    igraph_real_t m = communities->weight_sum;

    for (i = 0; i < communities->vertices_no; i++) {
        if (communities->item[i].size > 0) {
            result += (communities->item[i].weight_inside - communities->item[i].weight_all * communities->item[i].weight_all / m) / m;
        }
    }

    return result;
}

typedef struct {
    long int from;
    long int to;
    long int id;
} igraph_i_multilevel_link;

static int igraph_i_multilevel_link_cmp(const void *a, const void *b) {
    long int r = (((igraph_i_multilevel_link*)a)->from -
                  ((igraph_i_multilevel_link*)b)->from);
    if (r != 0) {
        return (int) r;
    }

    return (int) (((igraph_i_multilevel_link*)a)->to -
                  ((igraph_i_multilevel_link*)b)->to);
}

/* removes multiple edges and returns new edge id's for each edge in |E|log|E| */
static int igraph_i_multilevel_simplify_multiple(igraph_t *graph, igraph_vector_t *eids) {
    long int ecount = igraph_ecount(graph);
    long int i, l = -1, last_from = -1, last_to = -1;
    igraph_bool_t directed = igraph_is_directed(graph);
    igraph_integer_t from, to;
    igraph_vector_t edges;
    igraph_i_multilevel_link *links;

    /* Make sure there's enough space in eids to store the new edge IDs */
    IGRAPH_CHECK(igraph_vector_resize(eids, ecount));

    links = igraph_Calloc(ecount, igraph_i_multilevel_link);
    if (links == 0) {
        IGRAPH_ERROR("multi-level community structure detection failed", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, links);

    for (i = 0; i < ecount; i++) {
        igraph_edge(graph, (igraph_integer_t) i, &from, &to);
        links[i].from = from;
        links[i].to = to;
        links[i].id = i;
    }

    qsort((void*)links, (size_t) ecount, sizeof(igraph_i_multilevel_link),
          igraph_i_multilevel_link_cmp);

    IGRAPH_VECTOR_INIT_FINALLY(&edges, 0);
    for (i = 0; i < ecount; i++) {
        if (links[i].from == last_from && links[i].to == last_to) {
            VECTOR(*eids)[links[i].id] = l;
            continue;
        }

        last_from = links[i].from;
        last_to = links[i].to;

        igraph_vector_push_back(&edges, last_from);
        igraph_vector_push_back(&edges, last_to);

        l++;

        VECTOR(*eids)[links[i].id] = l;
    }

    igraph_Free(links);
    IGRAPH_FINALLY_CLEAN(1);

    igraph_destroy(graph);
    IGRAPH_CHECK(igraph_create(graph, &edges, igraph_vcount(graph), directed));

    igraph_vector_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);

    return 0;
}

typedef struct {
    long int community;
    igraph_real_t weight;
} igraph_i_multilevel_community_link;

static int igraph_i_multilevel_community_link_cmp(const void *a, const void *b) {
    return (int) (((igraph_i_multilevel_community_link*)a)->community -
                  ((igraph_i_multilevel_community_link*)b)->community);
}

/**
 * Given a graph, a community structure and a vertex ID, this method
 * calculates:
 *
 * - edges: the list of edge IDs that are incident on the vertex
 * - weight_all: the total weight of these edges
 * - weight_inside: the total weight of edges that stay within the same
 *   community where the given vertex is right now, excluding loop edges
 * - weight_loop: the total weight of loop edges
 * - links_community and links_weight: together these two vectors list the
 *   communities incident on this vertex and the total weight of edges
 *   pointing to these communities
 */
static int igraph_i_multilevel_community_links(
        const igraph_t *graph,
        const igraph_i_multilevel_community_list *communities,
        igraph_integer_t vertex, igraph_vector_t *edges,
        igraph_real_t *weight_all, igraph_real_t *weight_inside, igraph_real_t *weight_loop,
        igraph_vector_t *links_community, igraph_vector_t *links_weight) {

    long int i, n, last = -1, c = -1;
    igraph_real_t weight = 1;
    long int to, to_community;
    long int community = (long int) VECTOR(*(communities->membership))[(long int)vertex];
    igraph_i_multilevel_community_link *links;

    *weight_all = *weight_inside = *weight_loop = 0;

    igraph_vector_clear(links_community);
    igraph_vector_clear(links_weight);

    /* Get the list of incident edges */
    igraph_incident(graph, edges, vertex, IGRAPH_ALL);

    n = igraph_vector_size(edges);
    links = igraph_Calloc(n, igraph_i_multilevel_community_link);
    if (links == 0) {
        IGRAPH_ERROR("multi-level community structure detection failed", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, links);

    for (i = 0; i < n; i++) {
        long int eidx = (long int) VECTOR(*edges)[i];
        weight = VECTOR(*communities->weights)[eidx];

        to = IGRAPH_OTHER(graph, eidx, vertex);

        *weight_all += weight;
        if (to == vertex) {
            *weight_loop += weight;

            links[i].community = community;
            links[i].weight = 0;
            continue;
        }

        to_community = (long int)VECTOR(*(communities->membership))[to];
        if (community == to_community) {
            *weight_inside += weight;
        }

        /* debug("Link %ld (C: %ld) <-> %ld (C: %ld)\n", vertex, community, to, to_community); */

        links[i].community = to_community;
        links[i].weight = weight;
    }

    /* Sort links by community ID and merge the same */
    qsort((void*)links, (size_t) n, sizeof(igraph_i_multilevel_community_link),
          igraph_i_multilevel_community_link_cmp);
    for (i = 0; i < n; i++) {
        to_community = links[i].community;
        if (to_community != last) {
            igraph_vector_push_back(links_community, to_community);
            igraph_vector_push_back(links_weight, links[i].weight);
            last = to_community;
            c++;
        } else {
            VECTOR(*links_weight)[c] += links[i].weight;
        }
    }

    igraph_free(links);
    IGRAPH_FINALLY_CLEAN(1);

    return 0;
}

static igraph_real_t igraph_i_multilevel_community_modularity_gain(
        const igraph_i_multilevel_community_list *communities,
        igraph_integer_t community, igraph_integer_t vertex,
        igraph_real_t weight_all, igraph_real_t weight_inside) {
    IGRAPH_UNUSED(vertex);
    return weight_inside -
           communities->item[(long int)community].weight_all * weight_all / communities->weight_sum;
}

/* Shrinks communities into single vertices, keeping all the edges.
 * This method is internal because it destroys the graph in-place and
 * creates a new one -- this is fine for the multilevel community
 * detection where a copy of the original graph is used anyway.
 * The membership vector will also be rewritten by the underlying
 * igraph_membership_reindex call */
static int igraph_i_multilevel_shrink(igraph_t *graph, igraph_vector_t *membership) {
    igraph_vector_t edges;
    long int no_of_nodes = igraph_vcount(graph);
    long int no_of_edges = igraph_ecount(graph);
    igraph_bool_t directed = igraph_is_directed(graph);

    long int i;
    igraph_eit_t eit;

    if (no_of_nodes == 0) {
        return 0;
    }

    if (igraph_vector_size(membership) < no_of_nodes) {
        IGRAPH_ERROR("cannot shrink graph, membership vector too short",
                     IGRAPH_EINVAL);
    }

    IGRAPH_VECTOR_INIT_FINALLY(&edges, no_of_edges * 2);

    IGRAPH_CHECK(igraph_reindex_membership(membership, 0, NULL));

    /* Create the new edgelist */
    igraph_eit_create(graph, igraph_ess_all(IGRAPH_EDGEORDER_ID), &eit);
    IGRAPH_FINALLY(igraph_eit_destroy, &eit);
    i = 0;
    while (!IGRAPH_EIT_END(eit)) {
        igraph_integer_t from, to;
        IGRAPH_CHECK(igraph_edge(graph, IGRAPH_EIT_GET(eit), &from, &to));
        VECTOR(edges)[i++] = VECTOR(*membership)[(long int) from];
        VECTOR(edges)[i++] = VECTOR(*membership)[(long int) to];
        IGRAPH_EIT_NEXT(eit);
    }
    igraph_eit_destroy(&eit);
    IGRAPH_FINALLY_CLEAN(1);

    /* Create the new graph */
    igraph_destroy(graph);
    no_of_nodes = (long int) igraph_vector_max(membership) + 1;
    IGRAPH_CHECK(igraph_create(graph, &edges, (igraph_integer_t) no_of_nodes,
                               directed));

    igraph_vector_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);

    return 0;
}

/**
 * \ingroup communities
 * \function igraph_i_community_multilevel_step
 * \brief Performs a single step of the multi-level modularity optimization method
 *
 * This function implements a single step of the multi-level modularity optimization
 * algorithm for finding community structure, see VD Blondel, J-L Guillaume,
 * R Lambiotte and E Lefebvre: Fast unfolding of community hierarchies in large
 * networks, http://arxiv.org/abs/0803.0476 for the details.
 *
 * This function was contributed by Tom Gregorovic.
 *
 * \param graph   The input graph. It must be an undirected graph.
 * \param weights Numeric vector containing edge weights. If \c NULL, every edge
 *     has equal weight. The weights are expected to be non-negative.
 * \param membership The membership vector, the result is returned here.
 *     For each vertex it gives the ID of its community.
 * \param modularity The modularity of the partition is returned here.
 *     \c NULL means that the modularity is not needed.
 * \return Error code.
 *
 * Time complexity: in average near linear on sparse graphs.
 */
static int igraph_i_community_multilevel_step(
        igraph_t *graph,
        igraph_vector_t *weights,
        igraph_vector_t *membership,
        igraph_real_t *modularity) {

    long int i, j;
    long int vcount = igraph_vcount(graph);
    long int ecount = igraph_ecount(graph);
    igraph_integer_t ffrom, fto;
    igraph_real_t q, pass_q;
    int pass;
    igraph_bool_t changed = 0;
    igraph_vector_t links_community;
    igraph_vector_t links_weight;
    igraph_vector_t edges;
    igraph_vector_t temp_membership;
    igraph_i_multilevel_community_list communities;

    /* Initial sanity checks on the input parameters */
    if (igraph_is_directed(graph)) {
        IGRAPH_ERROR("multi-level community detection works for undirected graphs only",
                     IGRAPH_UNIMPLEMENTED);
    }
    if (igraph_vector_size(weights) < igraph_ecount(graph)) {
        IGRAPH_ERROR("multi-level community detection: weight vector too short", IGRAPH_EINVAL);
    }
    if (igraph_vector_any_smaller(weights, 0)) {
        IGRAPH_ERROR("weights must be positive", IGRAPH_EINVAL);
    }

    /* Initialize data structures */
    IGRAPH_VECTOR_INIT_FINALLY(&links_community, 0);
    IGRAPH_VECTOR_INIT_FINALLY(&links_weight, 0);
    IGRAPH_VECTOR_INIT_FINALLY(&edges, 0);
    IGRAPH_VECTOR_INIT_FINALLY(&temp_membership, vcount);
    IGRAPH_CHECK(igraph_vector_resize(membership, vcount));

    /* Initialize list of communities from graph vertices */
    communities.vertices_no = vcount;
    communities.communities_no = vcount;
    communities.weights = weights;
    communities.weight_sum = 2 * igraph_vector_sum(weights);
    communities.membership = membership;
    communities.item = igraph_Calloc(vcount, igraph_i_multilevel_community);
    if (communities.item == 0) {
        IGRAPH_ERROR("multi-level community structure detection failed", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, communities.item);

    /* Still initializing the communities data structure */
    for (i = 0; i < vcount; i++) {
        VECTOR(*communities.membership)[i] = i;
        communities.item[i].size = 1;
        communities.item[i].weight_inside = 0;
        communities.item[i].weight_all = 0;
    }

    /* Some more initialization :) */
    for (i = 0; i < ecount; i++) {
        igraph_real_t weight = 1;
        igraph_edge(graph, (igraph_integer_t) i, &ffrom, &fto);

        weight = VECTOR(*weights)[i];
        communities.item[(long int) ffrom].weight_all += weight;
        communities.item[(long int) fto].weight_all += weight;
        if (ffrom == fto) {
            communities.item[(long int) ffrom].weight_inside += 2 * weight;
        }
    }

    q = igraph_i_multilevel_community_modularity(&communities);
    pass = 1;

    do { /* Pass begin */
        long int temp_communities_no = communities.communities_no;

        pass_q = q;
        changed = 0;

        /* Save the current membership, it will be restored in case of worse result */
        IGRAPH_CHECK(igraph_vector_update(&temp_membership, communities.membership));

        for (i = 0; i < vcount; i++) {
            /* Exclude vertex from its current community */
            igraph_real_t weight_all = 0;
            igraph_real_t weight_inside = 0;
            igraph_real_t weight_loop = 0;
            igraph_real_t max_q_gain = 0;
            igraph_real_t max_weight;
            long int old_id, new_id, n;

            igraph_i_multilevel_community_links(graph, &communities,
                                                (igraph_integer_t) i, &edges,
                                                &weight_all, &weight_inside,
                                                &weight_loop, &links_community,
                                                &links_weight);

            old_id = (long int)VECTOR(*(communities.membership))[i];
            new_id = old_id;

            /* Update old community */
            igraph_vector_set(communities.membership, i, -1);
            communities.item[old_id].size--;
            if (communities.item[old_id].size == 0) {
                communities.communities_no--;
            }
            communities.item[old_id].weight_all -= weight_all;
            communities.item[old_id].weight_inside -= 2 * weight_inside + weight_loop;

            /* debug("Remove %ld all: %lf Inside: %lf\n", i, -weight_all, -2*weight_inside + weight_loop); */

            /* Find new community to join with the best modification gain */
            max_q_gain = 0;
            max_weight = weight_inside;
            n = igraph_vector_size(&links_community);

            for (j = 0; j < n; j++) {
                long int c = (long int) VECTOR(links_community)[j];
                igraph_real_t w = VECTOR(links_weight)[j];

                igraph_real_t q_gain =
                    igraph_i_multilevel_community_modularity_gain(&communities,
                            (igraph_integer_t) c,
                            (igraph_integer_t) i,
                            weight_all, w);
                /* debug("Link %ld -> %ld weight: %lf gain: %lf\n", i, c, (double) w, (double) q_gain); */
                if (q_gain > max_q_gain) {
                    new_id = c;
                    max_q_gain = q_gain;
                    max_weight = w;
                }
            }

            /* debug("Added vertex %ld to community %ld (gain %lf).\n", i, new_id, (double) max_q_gain); */

            /* Add vertex to "new" community and update it */
            igraph_vector_set(communities.membership, i, new_id);
            if (communities.item[new_id].size == 0) {
                communities.communities_no++;
            }
            communities.item[new_id].size++;
            communities.item[new_id].weight_all += weight_all;
            communities.item[new_id].weight_inside += 2 * max_weight + weight_loop;

            if (new_id != old_id) {
                changed++;
            }
        }

        q = igraph_i_multilevel_community_modularity(&communities);

        if (changed && (q > pass_q)) {
            /* debug("Pass %d (changed: %d) Communities: %ld Modularity from %lf to %lf\n",
              pass, changed, communities.communities_no, (double) pass_q, (double) q); */
            pass++;
        } else {
            /* No changes or the modularity became worse, restore last membership */
            IGRAPH_CHECK(igraph_vector_update(communities.membership, &temp_membership));
            communities.communities_no = temp_communities_no;
            break;
        }

        IGRAPH_ALLOW_INTERRUPTION();
    } while (changed && (q > pass_q)); /* Pass end */

    if (modularity) {
        *modularity = q;
    }

    /* debug("Result Communities: %ld Modularity: %lf\n",
      communities.communities_no, (double) q); */

    IGRAPH_CHECK(igraph_reindex_membership(membership, 0, NULL));

    /* Shrink the nodes of the graph according to the present community structure
     * and simplify the resulting graph */

    /* TODO: check if we really need to copy temp_membership */
    IGRAPH_CHECK(igraph_vector_update(&temp_membership, membership));
    IGRAPH_CHECK(igraph_i_multilevel_shrink(graph, &temp_membership));
    igraph_vector_destroy(&temp_membership);
    IGRAPH_FINALLY_CLEAN(1);

    /* Update edge weights after shrinking and simplification */
    /* Here we reuse the edges vector as we don't need the previous contents anymore */
    /* TODO: can we use igraph_simplify here? */
    IGRAPH_CHECK(igraph_i_multilevel_simplify_multiple(graph, &edges));

    /* We reuse the links_weight vector to store the old edge weights */
    IGRAPH_CHECK(igraph_vector_update(&links_weight, weights));
    igraph_vector_fill(weights, 0);

    for (i = 0; i < ecount; i++) {
        VECTOR(*weights)[(long int)VECTOR(edges)[i]] += VECTOR(links_weight)[i];
    }

    igraph_free(communities.item);
    igraph_vector_destroy(&links_community);
    igraph_vector_destroy(&links_weight);
    igraph_vector_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(4);

    return 0;
}

/**
 * \ingroup communities
 * \function igraph_community_multilevel
 * \brief Finding community structure by multi-level optimization of modularity
 *
 * This function implements the multi-level modularity optimization
 * algorithm for finding community structure, see
 * VD Blondel, J-L Guillaume, R Lambiotte and E Lefebvre: Fast unfolding of
 * community hierarchies in large networks, J Stat Mech P10008 (2008)
 * for the details (preprint: http://arxiv.org/abs/arXiv:0803.0476).
 *
 * It is based on the modularity measure and a hierarchical approach.
 * Initially, each vertex is assigned to a community on its own. In every step,
 * vertices are re-assigned to communities in a local, greedy way: each vertex
 * is moved to the community with which it achieves the highest contribution to
 * modularity. When no vertices can be reassigned, each community is considered
 * a vertex on its own, and the process starts again with the merged communities.
 * The process stops when there is only a single vertex left or when the modularity
 * cannot be increased any more in a step.
 *
 * This function was contributed by Tom Gregorovic.
 *
 * \param graph The input graph. It must be an undirected graph.
 * \param weights Numeric vector containing edge weights. If \c NULL, every edge
 *    has equal weight. The weights are expected to be non-negative.
 * \param membership The membership vector, the result is returned here.
 *    For each vertex it gives the ID of its community. The vector
 *    must be initialized and it will be resized accordingly.
 * \param memberships Numeric matrix that will contain the membership
 *     vector after each level, if not \c NULL. It must be initialized and
 *     it will be resized accordingly.
 * \param modularity Numeric vector that will contain the modularity score
 *     after each level, if not \c NULL. It must be initialized and it
 *     will be resized accordingly.
 * \return Error code.
 *
 * Time complexity: in average near linear on sparse graphs.
 *
 * \example examples/simple/igraph_community_multilevel.c
 */

int igraph_community_multilevel(const igraph_t *graph,
                                const igraph_vector_t *weights, igraph_vector_t *membership,
                                igraph_matrix_t *memberships, igraph_vector_t *modularity) {

    igraph_t g;
    igraph_vector_t w, m, level_membership;
    igraph_real_t prev_q = -1, q = -1;
    int i, level = 1;
    long int vcount = igraph_vcount(graph);

    /* Make a copy of the original graph, we will do the merges on the copy */
    IGRAPH_CHECK(igraph_copy(&g, graph));
    IGRAPH_FINALLY(igraph_destroy, &g);

    if (weights) {
        IGRAPH_CHECK(igraph_vector_copy(&w, weights));
        IGRAPH_FINALLY(igraph_vector_destroy, &w);
    } else {
        IGRAPH_VECTOR_INIT_FINALLY(&w, igraph_ecount(&g));
        igraph_vector_fill(&w, 1);
    }

    IGRAPH_VECTOR_INIT_FINALLY(&m, vcount);
    IGRAPH_VECTOR_INIT_FINALLY(&level_membership, vcount);

    if (memberships || membership) {
        /* Put each vertex in its own community */
        for (i = 0; i < vcount; i++) {
            VECTOR(level_membership)[i] = i;
        }
    }
    if (memberships) {
        /* Resize the membership matrix to have vcount columns and no rows */
        IGRAPH_CHECK(igraph_matrix_resize(memberships, 0, vcount));
    }
    if (modularity) {
        /* Clear the modularity vector */
        igraph_vector_clear(modularity);
    }

    while (1) {
        /* Remember the previous modularity and vertex count, do a single step */
        igraph_integer_t step_vcount = igraph_vcount(&g);

        prev_q = q;
        IGRAPH_CHECK(igraph_i_community_multilevel_step(&g, &w, &m, &q));

        /* Were there any merges? If not, we have to stop the process */
        if (igraph_vcount(&g) == step_vcount || q < prev_q) {
            break;
        }

        if (memberships || membership) {
            for (i = 0; i < vcount; i++) {
                /* Readjust the membership vector */
                VECTOR(level_membership)[i] = VECTOR(m)[(long int) VECTOR(level_membership)[i]];
            }
        }

        if (modularity) {
            /* If we have to return the modularity scores, add it to the modularity vector */
            IGRAPH_CHECK(igraph_vector_push_back(modularity, q));
        }

        if (memberships) {
            /* If we have to return the membership vectors at each level, store the new
             * membership vector */
            IGRAPH_CHECK(igraph_matrix_add_rows(memberships, 1));
            IGRAPH_CHECK(igraph_matrix_set_row(memberships, &level_membership, level - 1));
        }

        /* debug("Level: %d Communities: %ld Modularity: %f\n", level, (long int) igraph_vcount(&g),
          (double) q); */

        /* Increase the level counter */
        level++;
    }

    /* It might happen that there are no merges, so every vertex is in its
       own community. We still might want the modularity score for that. */
    if (modularity && igraph_vector_size(modularity) == 0) {
        igraph_vector_t tmp;
        igraph_real_t mod;
        int i;
        IGRAPH_VECTOR_INIT_FINALLY(&tmp, vcount);
        for (i = 0; i < vcount; i++) {
            VECTOR(tmp)[i] = i;
        }
        IGRAPH_CHECK(igraph_modularity(graph, &tmp, &mod, weights));
        igraph_vector_destroy(&tmp);
        IGRAPH_FINALLY_CLEAN(1);
        IGRAPH_CHECK(igraph_vector_resize(modularity, 1));
        VECTOR(*modularity)[0] = mod;
    }

    /* If we need the final membership vector, copy it to the output */
    if (membership) {
        IGRAPH_CHECK(igraph_vector_resize(membership, vcount));
        for (i = 0; i < vcount; i++) {
            VECTOR(*membership)[i] = VECTOR(level_membership)[i];
        }
    }

    /* Destroy the copy of the graph */
    igraph_destroy(&g);

    /* Destroy the temporary vectors */
    igraph_vector_destroy(&m);
    igraph_vector_destroy(&w);
    igraph_vector_destroy(&level_membership);
    IGRAPH_FINALLY_CLEAN(4);

    return 0;
}


static int igraph_i_compare_communities_vi(const igraph_vector_t *v1,
                                           const igraph_vector_t *v2, igraph_real_t* result);
static int igraph_i_compare_communities_nmi(const igraph_vector_t *v1,
                                            const igraph_vector_t *v2, igraph_real_t* result);
static int igraph_i_compare_communities_rand(const igraph_vector_t *v1,
                                             const igraph_vector_t *v2, igraph_real_t* result, igraph_bool_t adjust);
static int igraph_i_split_join_distance(const igraph_vector_t *v1,
                                        const igraph_vector_t *v2, igraph_integer_t* distance12,
                                        igraph_integer_t* distance21);

/**
 * \ingroup communities
 * \function igraph_compare_communities
 * \brief Compares community structures using various metrics
 *
 * This function assesses the distance between two community structures
 * using the variation of information (VI) metric of Meila (2003), the
 * normalized mutual information (NMI) of Danon et al (2005), the
 * split-join distance of van Dongen (2000), the Rand index of Rand (1971)
 * or the adjusted Rand index of Hubert and Arabie (1985).
 *
 * </para><para>
 * References:
 *
 * </para><para>
 * Meila M: Comparing clusterings by the variation of information.
 * In: Schölkopf B, Warmuth MK (eds.). Learning Theory and Kernel Machines:
 * 16th Annual Conference on Computational Learning Theory and 7th Kernel
 * Workshop, COLT/Kernel 2003, Washington, DC, USA. Lecture Notes in Computer
 * Science, vol. 2777, Springer, 2003. ISBN: 978-3-540-40720-1.
 *
 * </para><para>
 * Danon L, Diaz-Guilera A, Duch J, Arenas A: Comparing community structure
 * identification. J Stat Mech P09008, 2005.
 *
 * </para><para>
 * van Dongen S: Performance criteria for graph clustering and Markov cluster
 * experiments. Technical Report INS-R0012, National Research Institute for
 * Mathematics and Computer Science in the Netherlands, Amsterdam, May 2000.
 *
 * </para><para>
 * Rand WM: Objective criteria for the evaluation of clustering methods.
 * J Am Stat Assoc 66(336):846-850, 1971.
 *
 * </para><para>
 * Hubert L and Arabie P: Comparing partitions. Journal of Classification
 * 2:193-218, 1985.
 *
 * \param  comm1   the membership vector of the first community structure
 * \param  comm2   the membership vector of the second community structure
 * \param  result  the result is stored here.
 * \param  method  the comparison method to use. \c IGRAPH_COMMCMP_VI
 *                 selects the variation of information (VI) metric of
 *                 Meila (2003), \c IGRAPH_COMMCMP_NMI selects the
 *                 normalized mutual information measure proposed by
 *                 Danon et al (2005), \c IGRAPH_COMMCMP_SPLIT_JOIN
 *                 selects the split-join distance of van Dongen (2000),
 *                 \c IGRAPH_COMMCMP_RAND selects the unadjusted Rand
 *                 index (1971) and \c IGRAPH_COMMCMP_ADJUSTED_RAND
 *                 selects the adjusted Rand index.
 *
 * \return  Error code.
 *
 * Time complexity: O(n log(n)).
 */
int igraph_compare_communities(const igraph_vector_t *comm1,
                               const igraph_vector_t *comm2, igraph_real_t* result,
                               igraph_community_comparison_t method) {
    igraph_vector_t c1, c2;

    if (igraph_vector_size(comm1) != igraph_vector_size(comm2)) {
        IGRAPH_ERROR("community membership vectors have different lengths", IGRAPH_EINVAL);
    }

    /* Copy and reindex membership vectors to make sure they are continuous */
    IGRAPH_CHECK(igraph_vector_copy(&c1, comm1));
    IGRAPH_FINALLY(igraph_vector_destroy, &c1);

    IGRAPH_CHECK(igraph_vector_copy(&c2, comm2));
    IGRAPH_FINALLY(igraph_vector_destroy, &c2);

    IGRAPH_CHECK(igraph_reindex_membership(&c1, 0, NULL));
    IGRAPH_CHECK(igraph_reindex_membership(&c2, 0, NULL));

    switch (method) {
    case IGRAPH_COMMCMP_VI:
        IGRAPH_CHECK(igraph_i_compare_communities_vi(&c1, &c2, result));
        break;

    case IGRAPH_COMMCMP_NMI:
        IGRAPH_CHECK(igraph_i_compare_communities_nmi(&c1, &c2, result));
        break;

    case IGRAPH_COMMCMP_SPLIT_JOIN: {
        igraph_integer_t d12, d21;
        IGRAPH_CHECK(igraph_i_split_join_distance(&c1, &c2, &d12, &d21));
        *result = d12 + d21;
    }
    break;

    case IGRAPH_COMMCMP_RAND:
    case IGRAPH_COMMCMP_ADJUSTED_RAND:
        IGRAPH_CHECK(igraph_i_compare_communities_rand(&c1, &c2, result,
                     method == IGRAPH_COMMCMP_ADJUSTED_RAND));
        break;

    default:
        IGRAPH_ERROR("unknown community comparison method", IGRAPH_EINVAL);
    }

    /* Clean up everything */
    igraph_vector_destroy(&c1);
    igraph_vector_destroy(&c2);
    IGRAPH_FINALLY_CLEAN(2);

    return 0;
}

/**
 * \ingroup communities
 * \function igraph_split_join_distance
 * \brief Calculates the split-join distance of two community structures
 *
 * The split-join distance between partitions A and B is the sum of the
 * projection distance of A from B and the projection distance of B from
 * A. The projection distance is an asymmetric measure and it is defined
 * as follows:
 *
 * </para><para>
 * First, each set in partition A is evaluated against all sets in partition
 * B. For each set in partition A, the best matching set in partition B is
 * found and the overlap size is calculated. (Matching is quantified by the
 * size of the overlap between the two sets). Then, the maximal overlap sizes
 * for each set in A are summed together and subtracted from the number of
 * elements in A.
 *
 * </para><para>
 * The split-join distance will be returned in two arguments, \c distance12
 * will contain the projection distance of the first partition from the
 * second, while \c distance21 will be the projection distance of the second
 * partition from the first. This makes it easier to detect whether a
 * partition is a subpartition of the other, since in this case, the
 * corresponding distance will be zero.
 *
 * </para><para>
 * Reference:
 *
 * </para><para>
 * van Dongen S: Performance criteria for graph clustering and Markov cluster
 * experiments. Technical Report INS-R0012, National Research Institute for
 * Mathematics and Computer Science in the Netherlands, Amsterdam, May 2000.
 *
 * \param  comm1       the membership vector of the first community structure
 * \param  comm2       the membership vector of the second community structure
 * \param  distance12  pointer to an \c igraph_integer_t, the projection distance
 *                     of the first community structure from the second one will be
 *                     returned here.
 * \param  distance21  pointer to an \c igraph_integer_t, the projection distance
 *                     of the second community structure from the first one will be
 *                     returned here.
 * \return  Error code.
 *
 * \see \ref igraph_compare_communities() with the \c IGRAPH_COMMCMP_SPLIT_JOIN
 * method if you are not interested in the individual distances but only the sum
 * of them.
 *
 * Time complexity: O(n log(n)).
 */
int igraph_split_join_distance(const igraph_vector_t *comm1,
                               const igraph_vector_t *comm2, igraph_integer_t *distance12,
                               igraph_integer_t *distance21) {
    igraph_vector_t c1, c2;

    if (igraph_vector_size(comm1) != igraph_vector_size(comm2)) {
        IGRAPH_ERROR("community membership vectors have different lengths", IGRAPH_EINVAL);
    }

    /* Copy and reindex membership vectors to make sure they are continuous */
    IGRAPH_CHECK(igraph_vector_copy(&c1, comm1));
    IGRAPH_FINALLY(igraph_vector_destroy, &c1);

    IGRAPH_CHECK(igraph_vector_copy(&c2, comm2));
    IGRAPH_FINALLY(igraph_vector_destroy, &c2);

    IGRAPH_CHECK(igraph_reindex_membership(&c1, 0, NULL));
    IGRAPH_CHECK(igraph_reindex_membership(&c2, 0, NULL));

    IGRAPH_CHECK(igraph_i_split_join_distance(&c1, &c2, distance12, distance21));

    /* Clean up everything */
    igraph_vector_destroy(&c1);
    igraph_vector_destroy(&c2);
    IGRAPH_FINALLY_CLEAN(2);

    return 0;
}

/**
 * Calculates the entropy and the mutual information for two reindexed community
 * membership vectors v1 and v2. This is needed by both Meila's and Danon's
 * community comparison measure.
 */
static int igraph_i_entropy_and_mutual_information(const igraph_vector_t* v1,
        const igraph_vector_t* v2, double* h1, double* h2, double* mut_inf) {
    long int i, n = igraph_vector_size(v1);
    long int k1 = (long int)igraph_vector_max(v1) + 1;
    long int k2 = (long int)igraph_vector_max(v2) + 1;
    double *p1, *p2;
    igraph_spmatrix_t m;
    igraph_spmatrix_iter_t mit;

    p1 = igraph_Calloc(k1, double);
    if (p1 == 0) {
        IGRAPH_ERROR("igraph_i_entropy_and_mutual_information failed", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, p1);
    p2 = igraph_Calloc(k2, double);
    if (p2 == 0) {
        IGRAPH_ERROR("igraph_i_entropy_and_mutual_information failed", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, p2);

    /* Calculate the entropy of v1 */
    *h1 = 0.0;
    for (i = 0; i < n; i++) {
        p1[(long int)VECTOR(*v1)[i]]++;
    }
    for (i = 0; i < k1; i++) {
        p1[i] /= n;
        *h1 -= p1[i] * log(p1[i]);
    }

    /* Calculate the entropy of v2 */
    *h2 = 0.0;
    for (i = 0; i < n; i++) {
        p2[(long int)VECTOR(*v2)[i]]++;
    }
    for (i = 0; i < k2; i++) {
        p2[i] /= n;
        *h2 -= p2[i] * log(p2[i]);
    }

    /* We will only need the logs of p1 and p2 from now on */
    for (i = 0; i < k1; i++) {
        p1[i] = log(p1[i]);
    }
    for (i = 0; i < k2; i++) {
        p2[i] = log(p2[i]);
    }

    /* Calculate the mutual information of v1 and v2 */
    *mut_inf = 0.0;
    IGRAPH_CHECK(igraph_spmatrix_init(&m, k1, k2));
    IGRAPH_FINALLY(igraph_spmatrix_destroy, &m);
    for (i = 0; i < n; i++) {
        IGRAPH_CHECK(igraph_spmatrix_add_e(&m,
                                           (int)VECTOR(*v1)[i], (int)VECTOR(*v2)[i], 1));
    }
    IGRAPH_CHECK(igraph_spmatrix_iter_create(&mit, &m));
    IGRAPH_FINALLY(igraph_spmatrix_iter_destroy, &mit);
    while (!igraph_spmatrix_iter_end(&mit)) {
        double p = mit.value / n;
        *mut_inf += p * (log(p) - p1[mit.ri] - p2[mit.ci]);
        igraph_spmatrix_iter_next(&mit);
    }

    igraph_spmatrix_iter_destroy(&mit);
    igraph_spmatrix_destroy(&m);
    igraph_Free(p1); igraph_Free(p2);

    IGRAPH_FINALLY_CLEAN(4);

    return 0;
}

/**
 * Implementation of the normalized mutual information (NMI) measure of
 * Danon et al. This function assumes that the community membership
 * vectors have already been normalized using igraph_reindex_communities().
 *
 * </para><para>
 * Reference: Danon L, Diaz-Guilera A, Duch J, Arenas A: Comparing community
 * structure identification. J Stat Mech P09008, 2005.
 *
 * </para><para>
 * Time complexity: O(n log(n))
 */
static int igraph_i_compare_communities_nmi(const igraph_vector_t *v1, const igraph_vector_t *v2,
                                     igraph_real_t* result) {
    double h1, h2, mut_inf;

    IGRAPH_CHECK(igraph_i_entropy_and_mutual_information(v1, v2, &h1, &h2, &mut_inf));

    if (h1 == 0 && h2 == 0) {
        *result = 1;
    } else {
        *result = 2 * mut_inf / (h1 + h2);
    }

    return IGRAPH_SUCCESS;
}

/**
 * Implementation of the variation of information metric (VI) of
 * Meila et al. This function assumes that the community membership
 * vectors have already been normalized using igraph_reindex_communities().
 *
 * </para><para>
 * Reference: Meila M: Comparing clusterings by the variation of information.
 * In: Schölkopf B, Warmuth MK (eds.). Learning Theory and Kernel Machines:
 * 16th Annual Conference on Computational Learning Theory and 7th Kernel
 * Workshop, COLT/Kernel 2003, Washington, DC, USA. Lecture Notes in Computer
 * Science, vol. 2777, Springer, 2003. ISBN: 978-3-540-40720-1.
 *
 * </para><para>
 * Time complexity: O(n log(n))
 */
static int igraph_i_compare_communities_vi(const igraph_vector_t *v1, const igraph_vector_t *v2,
                                    igraph_real_t* result) {
    double h1, h2, mut_inf;

    IGRAPH_CHECK(igraph_i_entropy_and_mutual_information(v1, v2, &h1, &h2, &mut_inf));
    *result = h1 + h2 - 2 * mut_inf;

    return IGRAPH_SUCCESS;
}

/**
 * \brief Calculates the confusion matrix for two clusterings.
 *
 * </para><para>
 * This function assumes that the community membership vectors have already
 * been normalized using igraph_reindex_communities().
 *
 * </para><para>
 * Time complexity: O(n log(max(k1, k2))), where n is the number of vertices, k1
 * and k2 are the number of clusters in each of the clusterings.
 */
static int igraph_i_confusion_matrix(const igraph_vector_t *v1, const igraph_vector_t *v2,
                              igraph_spmatrix_t *m) {
    long int k1 = (long int)igraph_vector_max(v1) + 1;
    long int k2 = (long int)igraph_vector_max(v2) + 1;
    long int i, n = igraph_vector_size(v1);

    IGRAPH_CHECK(igraph_spmatrix_resize(m, k1, k2));
    for (i = 0; i < n; i++) {
        IGRAPH_CHECK(igraph_spmatrix_add_e(m,
                                           (int)VECTOR(*v1)[i], (int)VECTOR(*v2)[i], 1));
    }

    return IGRAPH_SUCCESS;
}

/**
 * Implementation of the split-join distance of van Dongen.
 *
 * </para><para>
 * This function assumes that the community membership vectors have already
 * been normalized using igraph_reindex_communities().
 *
 * </para><para>
 * Reference: van Dongen S: Performance criteria for graph clustering and Markov
 * cluster experiments. Technical Report INS-R0012, National Research Institute
 * for Mathematics and Computer Science in the Netherlands, Amsterdam, May 2000.
 *
 * </para><para>
 * Time complexity: O(n log(max(k1, k2))), where n is the number of vertices, k1
 * and k2 are the number of clusters in each of the clusterings.
 */
static int igraph_i_split_join_distance(const igraph_vector_t *v1, const igraph_vector_t *v2,
                                 igraph_integer_t* distance12, igraph_integer_t* distance21) {
    long int n = igraph_vector_size(v1);
    igraph_vector_t rowmax, colmax;
    igraph_spmatrix_t m;
    igraph_spmatrix_iter_t mit;

    /* Calculate the confusion matrix */
    IGRAPH_CHECK(igraph_spmatrix_init(&m, 1, 1));
    IGRAPH_FINALLY(igraph_spmatrix_destroy, &m);
    IGRAPH_CHECK(igraph_i_confusion_matrix(v1, v2, &m));

    /* Initialize vectors that will store the row/columnwise maxima */
    IGRAPH_VECTOR_INIT_FINALLY(&rowmax, igraph_spmatrix_nrow(&m));
    IGRAPH_VECTOR_INIT_FINALLY(&colmax, igraph_spmatrix_ncol(&m));

    /* Find the row/columnwise maxima */
    IGRAPH_CHECK(igraph_spmatrix_iter_create(&mit, &m));
    IGRAPH_FINALLY(igraph_spmatrix_iter_destroy, &mit);
    while (!igraph_spmatrix_iter_end(&mit)) {
        if (mit.value > VECTOR(rowmax)[mit.ri]) {
            VECTOR(rowmax)[mit.ri] = mit.value;
        }
        if (mit.value > VECTOR(colmax)[mit.ci]) {
            VECTOR(colmax)[mit.ci] = mit.value;
        }
        igraph_spmatrix_iter_next(&mit);
    }
    igraph_spmatrix_iter_destroy(&mit);
    IGRAPH_FINALLY_CLEAN(1);

    /* Calculate the distances */
    *distance12 = (igraph_integer_t) (n - igraph_vector_sum(&rowmax));
    *distance21 = (igraph_integer_t) (n - igraph_vector_sum(&colmax));

    igraph_vector_destroy(&rowmax);
    igraph_vector_destroy(&colmax);
    igraph_spmatrix_destroy(&m);
    IGRAPH_FINALLY_CLEAN(3);

    return IGRAPH_SUCCESS;
}

/**
 * Implementation of the adjusted and unadjusted Rand indices.
 *
 * </para><para>
 * This function assumes that the community membership vectors have already
 * been normalized using igraph_reindex_communities().
 *
 * </para><para>
 * References:
 *
 * </para><para>
 * Rand WM: Objective criteria for the evaluation of clustering methods. J Am
 * Stat Assoc 66(336):846-850, 1971.
 *
 * </para><para>
 * Hubert L and Arabie P: Comparing partitions. Journal of Classification
 * 2:193-218, 1985.
 *
 * </para><para>
 * Time complexity: O(n log(max(k1, k2))), where n is the number of vertices, k1
 * and k2 are the number of clusters in each of the clusterings.
 */
static int igraph_i_compare_communities_rand(
        const igraph_vector_t *v1, const igraph_vector_t *v2,
        igraph_real_t *result, igraph_bool_t adjust) {
    igraph_spmatrix_t m;
    igraph_spmatrix_iter_t mit;
    igraph_vector_t rowsums, colsums;
    long int i, nrow, ncol;
    double rand, n;
    double frac_pairs_in_1, frac_pairs_in_2;

    /* Calculate the confusion matrix */
    IGRAPH_CHECK(igraph_spmatrix_init(&m, 1, 1));
    IGRAPH_FINALLY(igraph_spmatrix_destroy, &m);
    IGRAPH_CHECK(igraph_i_confusion_matrix(v1, v2, &m));

    /* The unadjusted Rand index is defined as (a+d) / (a+b+c+d), where:
     *
     * - a is the number of pairs in the same cluster both in v1 and v2. This
     *   equals the sum of n(i,j) choose 2 for all i and j.
     *
     * - b is the number of pairs in the same cluster in v1 and in different
     *   clusters in v2. This is sum n(i,*) choose 2 for all i minus a.
     *   n(i,*) is the number of elements in cluster i in v1.
     *
     * - c is the number of pairs in the same cluster in v2 and in different
     *   clusters in v1. This is sum n(*,j) choose 2 for all j minus a.
     *   n(*,j) is the number of elements in cluster j in v2.
     *
     * - d is (n choose 2) - a - b - c.
     *
     * Therefore, a+d = (n choose 2) - b - c
     *                = (n choose 2) - sum (n(i,*) choose 2)
     *                               - sum (n(*,j) choose 2)
     *                               + 2 * sum (n(i,j) choose 2).
     *
     * Since a+b+c+d = (n choose 2) and this goes in the denominator, we can
     * just as well start dividing each term in a+d by (n choose 2), which
     * yields:
     *
     * 1 - sum( n(i,*)/n * (n(i,*)-1)/(n-1) )
     *   - sum( n(*,i)/n * (n(*,i)-1)/(n-1) )
     *   + sum( n(i,j)/n * (n(i,j)-1)/(n-1) ) * 2
     */

    /* Calculate row and column sums */
    nrow = igraph_spmatrix_nrow(&m);
    ncol = igraph_spmatrix_ncol(&m);
    n = igraph_vector_size(v1) + 0.0;
    IGRAPH_VECTOR_INIT_FINALLY(&rowsums, nrow);
    IGRAPH_VECTOR_INIT_FINALLY(&colsums, ncol);
    IGRAPH_CHECK(igraph_spmatrix_rowsums(&m, &rowsums));
    IGRAPH_CHECK(igraph_spmatrix_colsums(&m, &colsums));

    /* Start calculating the unadjusted Rand index */
    rand = 0.0;
    IGRAPH_CHECK(igraph_spmatrix_iter_create(&mit, &m));
    IGRAPH_FINALLY(igraph_spmatrix_iter_destroy, &mit);
    while (!igraph_spmatrix_iter_end(&mit)) {
        rand += (mit.value / n) * (mit.value - 1) / (n - 1);
        igraph_spmatrix_iter_next(&mit);
    }
    igraph_spmatrix_iter_destroy(&mit);
    IGRAPH_FINALLY_CLEAN(1);

    frac_pairs_in_1 = frac_pairs_in_2 = 0.0;
    for (i = 0; i < nrow; i++) {
        frac_pairs_in_1 += (VECTOR(rowsums)[i] / n) * (VECTOR(rowsums)[i] - 1) / (n - 1);
    }
    for (i = 0; i < ncol; i++) {
        frac_pairs_in_2 += (VECTOR(colsums)[i] / n) * (VECTOR(colsums)[i] - 1) / (n - 1);
    }

    rand = 1.0 + 2 * rand - frac_pairs_in_1 - frac_pairs_in_2;

    if (adjust) {
        double expected = frac_pairs_in_1 * frac_pairs_in_2 +
                          (1 - frac_pairs_in_1) * (1 - frac_pairs_in_2);
        rand = (rand - expected) / (1 - expected);
    }

    igraph_vector_destroy(&rowsums);
    igraph_vector_destroy(&colsums);
    igraph_spmatrix_destroy(&m);
    IGRAPH_FINALLY_CLEAN(3);

    *result = rand;

    return IGRAPH_SUCCESS;
}
