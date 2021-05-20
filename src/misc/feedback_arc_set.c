/* -*- mode: C -*-  */
/* vim:set ts=4 sw=4 sts=4 et: */
/*
   IGraph library.
   Copyright (C) 2011-2012  Gabor Csardi <csardi.gabor@gmail.com>
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

#include "igraph_components.h"
#include "igraph_dqueue.h"
#include "igraph_interface.h"
#include "igraph_memory.h"
#include "igraph_structural.h"
#include "igraph_visitor.h"

#include "internal/glpk_support.h"
#include "misc/feedback_arc_set.h"

int igraph_i_feedback_arc_set_ip(const igraph_t *graph, igraph_vector_t *result,
                                 const igraph_vector_t *weights);


/**
 * \ingroup structural
 * \function igraph_feedback_arc_set
 * \brief Calculates a feedback arc set of the graph using different
 *        algorithms.
 *
 * </para><para>
 * A feedback arc set is a set of edges whose removal makes the graph acyclic.
 * We are usually interested in \em minimum feedback arc sets, i.e. sets of edges
 * whose total weight is minimal among all the feedback arc sets.
 *
 * </para><para>
 * For undirected graphs, the problem is simple: one has to find a maximum weight
 * spanning tree and then remove all the edges not in the spanning tree. For directed
 * graphs, this is an NP-hard problem, and various heuristics are usually used to
 * find an approximate solution to the problem. This function implements a few of
 * these heuristics.
 *
 * \param graph  The graph object.
 * \param result An initialized vector, the result will be returned here.
 * \param weights Weight vector or NULL if no weights are specified.
 * \param algo   The algorithm to use to solve the problem if the graph is directed.
 *        Possible values:
 *        \clist
 *        \cli IGRAPH_FAS_EXACT_IP
 *          Finds a \em minimum feedback arc set using integer programming (IP).
 *          The complexity of this algorithm is exponential of course.
 *        \cli IGRAPH_FAS_APPROX_EADES
 *          Finds a feedback arc set using the heuristic of Eades, Lin and
 *          Smyth (1993). This is guaranteed to be smaller than |E|/2 - |V|/6,
 *          and it is linear in the number of edges (i.e. O(|E|)).
 *          For more details, see Eades P, Lin X and Smyth WF: A fast and effective
 *          heuristic for the feedback arc set problem. In: Proc Inf Process Lett
 *          319-323, 1993.
 *        \endclist
 *
 * \return Error code:
 *         \c IGRAPH_EINVAL if an unknown method was specified or the weight vector
 *            is invalid.
 *
 * \example examples/simple/igraph_feedback_arc_set.c
 * \example examples/simple/igraph_feedback_arc_set_ip.c
 *
 * Time complexity: depends on \p algo, see the time complexities there.
 */
int igraph_feedback_arc_set(const igraph_t *graph, igraph_vector_t *result,
                            const igraph_vector_t *weights, igraph_fas_algorithm_t algo) {

    if (weights && igraph_vector_size(weights) < igraph_ecount(graph))
        IGRAPH_ERROR("cannot calculate feedback arc set, weight vector too short",
                     IGRAPH_EINVAL);

    if (!igraph_is_directed(graph)) {
        return igraph_i_feedback_arc_set_undirected(graph, result, weights, 0);
    }

    switch (algo) {
    case IGRAPH_FAS_EXACT_IP:
        return igraph_i_feedback_arc_set_ip(graph, result, weights);

    case IGRAPH_FAS_APPROX_EADES:
        return igraph_i_feedback_arc_set_eades(graph, result, weights, 0);

    default:
        IGRAPH_ERROR("Invalid algorithm", IGRAPH_EINVAL);
    }
}

/**
 * Solves the feedback arc set problem for undirected graphs.
 */
int igraph_i_feedback_arc_set_undirected(const igraph_t *graph, igraph_vector_t *result,
        const igraph_vector_t *weights, igraph_vector_t *layering) {
    igraph_vector_t edges;
    long int i, j, n, no_of_nodes = igraph_vcount(graph);

    IGRAPH_VECTOR_INIT_FINALLY(&edges, no_of_nodes - 1);
    if (weights) {
        /* Find a maximum weight spanning tree. igraph has a routine for minimum
         * spanning trees, so we negate the weights */
        igraph_vector_t vcopy;
        IGRAPH_CHECK(igraph_vector_copy(&vcopy, weights));
        IGRAPH_FINALLY(igraph_vector_destroy, &vcopy);
        igraph_vector_scale(&vcopy, -1);
        IGRAPH_CHECK(igraph_minimum_spanning_tree(graph, &edges, &vcopy));
        igraph_vector_destroy(&vcopy);
        IGRAPH_FINALLY_CLEAN(1);
    } else {
        /* Any spanning tree will do */
        IGRAPH_CHECK(igraph_minimum_spanning_tree(graph, &edges, 0));
    }

    /* Now we have a bunch of edges that constitute a spanning forest. We have
     * to come up with a layering, and return those edges that are not in the
     * spanning forest */
    igraph_vector_sort(&edges);
    IGRAPH_CHECK(igraph_vector_push_back(&edges, -1));  /* guard element */

    if (result != 0) {
        igraph_vector_clear(result);
        n = igraph_ecount(graph);
        for (i = 0, j = 0; i < n; i++) {
            if (i == VECTOR(edges)[j]) {
                j++;
                continue;
            }
            IGRAPH_CHECK(igraph_vector_push_back(result, i));
        }
    }

    if (layering != 0) {
        igraph_vector_t degrees;
        igraph_vector_t roots;

        IGRAPH_VECTOR_INIT_FINALLY(&degrees, no_of_nodes);
        IGRAPH_VECTOR_INIT_FINALLY(&roots, no_of_nodes);

        IGRAPH_CHECK(igraph_strength(graph, &degrees, igraph_vss_all(),
                                     IGRAPH_ALL, 0, weights));
        IGRAPH_CHECK((int) igraph_vector_qsort_ind(&degrees, &roots,
                     /* descending = */ 1));
        IGRAPH_CHECK(igraph_bfs(graph,
                                /* root = */ 0,
                                /* roots = */ &roots,
                                /* mode = */ IGRAPH_OUT,
                                /* unreachable = */ 0,
                                /* restricted = */ 0,
                                /* order = */ 0,
                                /* rank = */ 0,
                                /* father = */ 0,
                                /* pred = */ 0,
                                /* succ = */ 0,
                                /* dist = */ layering,
                                /* callback = */ 0,
                                /* extra = */ 0));

        igraph_vector_destroy(&degrees);
        igraph_vector_destroy(&roots);
        IGRAPH_FINALLY_CLEAN(2);
    }

    igraph_vector_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}

/**
 * Solves the feedback arc set problem using the heuristics of Eades et al.
 */
int igraph_i_feedback_arc_set_eades(const igraph_t *graph, igraph_vector_t *result,
                                    const igraph_vector_t *weights, igraph_vector_t *layers) {
    long int i, j, k, v, eid, no_of_nodes = igraph_vcount(graph), nodes_left;
    igraph_dqueue_t sources, sinks;
    igraph_vector_t neis;
    igraph_vector_t indegrees, outdegrees;
    igraph_vector_t instrengths, outstrengths;
    long int* ordering;
    long int order_next_pos = 0, order_next_neg = -1;
    igraph_real_t diff, maxdiff;

    ordering = IGRAPH_CALLOC(no_of_nodes, long int);
    IGRAPH_FINALLY(igraph_free, ordering);

    IGRAPH_VECTOR_INIT_FINALLY(&indegrees, no_of_nodes);
    IGRAPH_VECTOR_INIT_FINALLY(&outdegrees, no_of_nodes);
    IGRAPH_VECTOR_INIT_FINALLY(&instrengths, no_of_nodes);
    IGRAPH_VECTOR_INIT_FINALLY(&outstrengths, no_of_nodes);
    IGRAPH_VECTOR_INIT_FINALLY(&neis, 0);
    IGRAPH_CHECK(igraph_dqueue_init(&sources, 0));
    IGRAPH_FINALLY(igraph_dqueue_destroy, &sources);
    IGRAPH_CHECK(igraph_dqueue_init(&sinks, 0));
    IGRAPH_FINALLY(igraph_dqueue_destroy, &sinks);

    IGRAPH_CHECK(igraph_degree(graph, &indegrees, igraph_vss_all(), IGRAPH_IN, 0));
    IGRAPH_CHECK(igraph_degree(graph, &outdegrees, igraph_vss_all(), IGRAPH_OUT, 0));

    if (weights) {
        IGRAPH_CHECK(igraph_strength(graph, &instrengths, igraph_vss_all(), IGRAPH_IN, 0, weights));
        IGRAPH_CHECK(igraph_strength(graph, &outstrengths, igraph_vss_all(), IGRAPH_OUT, 0, weights));
    } else {
        IGRAPH_CHECK(igraph_vector_update(&instrengths, &indegrees));
        IGRAPH_CHECK(igraph_vector_update(&outstrengths, &outdegrees));
    }

    /* Find initial sources and sinks */
    nodes_left = no_of_nodes;
    for (i = 0; i < no_of_nodes; i++) {
        if (VECTOR(indegrees)[i] == 0) {
            if (VECTOR(outdegrees)[i] == 0) {
                /* Isolated vertex, we simply ignore it */
                nodes_left--;
                ordering[i] = order_next_pos++;
                VECTOR(indegrees)[i] = VECTOR(outdegrees)[i] = -1;
            } else {
                /* This is a source */
                igraph_dqueue_push(&sources, i);
            }
        } else if (VECTOR(outdegrees)[i] == 0) {
            /* This is a sink */
            igraph_dqueue_push(&sinks, i);
        }
    }

    /* While we have any nodes left... */
    while (nodes_left > 0) {
        /* (1) Remove the sources one by one */
        while (!igraph_dqueue_empty(&sources)) {
            i = (long)igraph_dqueue_pop(&sources);
            /* Add the node to the ordering */
            ordering[i] = order_next_pos++;
            /* Exclude the node from further searches */
            VECTOR(indegrees)[i] = VECTOR(outdegrees)[i] = -1;
            /* Get the neighbors and decrease their degrees */
            IGRAPH_CHECK(igraph_incident(graph, &neis, (igraph_integer_t) i,
                                         IGRAPH_OUT));
            j = igraph_vector_size(&neis);
            for (i = 0; i < j; i++) {
                eid = (long int) VECTOR(neis)[i];
                k = IGRAPH_TO(graph, eid);
                if (VECTOR(indegrees)[k] <= 0) {
                    /* Already removed, continue */
                    continue;
                }
                VECTOR(indegrees)[k]--;
                VECTOR(instrengths)[k] -= (weights ? VECTOR(*weights)[eid] : 1.0);
                if (VECTOR(indegrees)[k] == 0) {
                    IGRAPH_CHECK(igraph_dqueue_push(&sources, k));
                }
            }
            nodes_left--;
        }

        /* (2) Remove the sinks one by one */
        while (!igraph_dqueue_empty(&sinks)) {
            i = (long)igraph_dqueue_pop(&sinks);
            /* Maybe the vertex became sink and source at the same time, hence it
             * was already removed in the previous iteration. Check it. */
            if (VECTOR(indegrees)[i] < 0) {
                continue;
            }
            /* Add the node to the ordering */
            ordering[i] = order_next_neg--;
            /* Exclude the node from further searches */
            VECTOR(indegrees)[i] = VECTOR(outdegrees)[i] = -1;
            /* Get the neighbors and decrease their degrees */
            IGRAPH_CHECK(igraph_incident(graph, &neis, (igraph_integer_t) i,
                                         IGRAPH_IN));
            j = igraph_vector_size(&neis);
            for (i = 0; i < j; i++) {
                eid = (long int) VECTOR(neis)[i];
                k = IGRAPH_FROM(graph, eid);
                if (VECTOR(outdegrees)[k] <= 0) {
                    /* Already removed, continue */
                    continue;
                }
                VECTOR(outdegrees)[k]--;
                VECTOR(outstrengths)[k] -= (weights ? VECTOR(*weights)[eid] : 1.0);
                if (VECTOR(outdegrees)[k] == 0) {
                    IGRAPH_CHECK(igraph_dqueue_push(&sinks, k));
                }
            }
            nodes_left--;
        }

        /* (3) No more sources or sinks. Find the node with the largest
         * difference between its out-strength and in-strength */
        v = -1; maxdiff = -IGRAPH_INFINITY;
        for (i = 0; i < no_of_nodes; i++) {
            if (VECTOR(outdegrees)[i] < 0) {
                continue;
            }
            diff = VECTOR(outstrengths)[i] - VECTOR(instrengths)[i];
            if (diff > maxdiff) {
                maxdiff = diff;
                v = i;
            }
        }
        if (v >= 0) {
            /* Remove vertex v */
            ordering[v] = order_next_pos++;
            /* Remove outgoing edges */
            IGRAPH_CHECK(igraph_incident(graph, &neis, (igraph_integer_t) v,
                                         IGRAPH_OUT));
            j = igraph_vector_size(&neis);
            for (i = 0; i < j; i++) {
                eid = (long int) VECTOR(neis)[i];
                k = IGRAPH_TO(graph, eid);
                if (VECTOR(indegrees)[k] <= 0) {
                    /* Already removed, continue */
                    continue;
                }
                VECTOR(indegrees)[k]--;
                VECTOR(instrengths)[k] -= (weights ? VECTOR(*weights)[eid] : 1.0);
                if (VECTOR(indegrees)[k] == 0) {
                    IGRAPH_CHECK(igraph_dqueue_push(&sources, k));
                }
            }
            /* Remove incoming edges */
            IGRAPH_CHECK(igraph_incident(graph, &neis, (igraph_integer_t) v,
                                         IGRAPH_IN));
            j = igraph_vector_size(&neis);
            for (i = 0; i < j; i++) {
                eid = (long int) VECTOR(neis)[i];
                k = IGRAPH_FROM(graph, eid);
                if (VECTOR(outdegrees)[k] <= 0) {
                    /* Already removed, continue */
                    continue;
                }
                VECTOR(outdegrees)[k]--;
                VECTOR(outstrengths)[k] -= (weights ? VECTOR(*weights)[eid] : 1.0);
                if (VECTOR(outdegrees)[k] == 0 && VECTOR(indegrees)[k] > 0) {
                    IGRAPH_CHECK(igraph_dqueue_push(&sinks, k));
                }
            }

            VECTOR(outdegrees)[v] = -1;
            VECTOR(indegrees)[v] = -1;
            nodes_left--;
        }
    }

    igraph_dqueue_destroy(&sinks);
    igraph_dqueue_destroy(&sources);
    igraph_vector_destroy(&neis);
    igraph_vector_destroy(&outstrengths);
    igraph_vector_destroy(&instrengths);
    igraph_vector_destroy(&outdegrees);
    igraph_vector_destroy(&indegrees);
    IGRAPH_FINALLY_CLEAN(7);

    /* Tidy up the ordering */
    for (i = 0; i < no_of_nodes; i++) {
        if (ordering[i] < 0) {
            ordering[i] += no_of_nodes;
        }
    }

    /* Find the feedback edges based on the ordering */
    if (result != 0) {
        igraph_vector_clear(result);
        j = igraph_ecount(graph);
        for (i = 0; i < j; i++) {
            long int from = IGRAPH_FROM(graph, i), to = IGRAPH_TO(graph, i);
            if (from == to || ordering[from] > ordering[to]) {
                IGRAPH_CHECK(igraph_vector_push_back(result, i));
            }
        }
    }

    /* If we have also requested a layering, return that as well */
    if (layers != 0) {
        igraph_vector_t ranks;
        igraph_vector_long_t order_vec;

        IGRAPH_CHECK(igraph_vector_resize(layers, no_of_nodes));
        igraph_vector_null(layers);

        igraph_vector_long_view(&order_vec, ordering, no_of_nodes);

        IGRAPH_VECTOR_INIT_FINALLY(&neis, 0);
        IGRAPH_VECTOR_INIT_FINALLY(&ranks, 0);

        IGRAPH_CHECK((int) igraph_vector_long_qsort_ind(&order_vec, &ranks, 0));

        for (i = 0; i < no_of_nodes; i++) {
            long int from = (long int) VECTOR(ranks)[i];
            IGRAPH_CHECK(igraph_neighbors(graph, &neis, (igraph_integer_t) from,
                                          IGRAPH_OUT));
            k = igraph_vector_size(&neis);
            for (j = 0; j < k; j++) {
                long int to = (long int) VECTOR(neis)[j];
                if (from == to) {
                    continue;
                }
                if (ordering[from] > ordering[to]) {
                    continue;
                }
                if (VECTOR(*layers)[to] < VECTOR(*layers)[from] + 1) {
                    VECTOR(*layers)[to] = VECTOR(*layers)[from] + 1;
                }
            }
        }

        igraph_vector_destroy(&neis);
        igraph_vector_destroy(&ranks);
        IGRAPH_FINALLY_CLEAN(2);
    }

    /* Free the ordering vector */
    igraph_free(ordering);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}

/**
 * Solves the feedback arc set problem using integer programming.
 */
int igraph_i_feedback_arc_set_ip(const igraph_t *graph, igraph_vector_t *result,
                                 const igraph_vector_t *weights) {
#ifndef HAVE_GLPK
    IGRAPH_ERROR("GLPK is not available", IGRAPH_UNIMPLEMENTED);
#else

    igraph_integer_t no_of_components;
    igraph_integer_t no_of_vertices = igraph_vcount(graph);
    igraph_integer_t no_of_edges = igraph_ecount(graph);
    igraph_vector_t membership, ordering, vertex_remapping;
    igraph_vector_ptr_t vertices_by_components, edges_by_components;
    long int i, j, k, l, m, n, from, to;
    igraph_real_t weight;
    glp_prob *ip;
    glp_iocp parm;

    IGRAPH_VECTOR_INIT_FINALLY(&membership, 0);
    IGRAPH_VECTOR_INIT_FINALLY(&ordering, 0);
    IGRAPH_VECTOR_INIT_FINALLY(&vertex_remapping, no_of_vertices);

    igraph_vector_clear(result);

    /* Decompose the graph into connected components */
    IGRAPH_CHECK(igraph_clusters(graph, &membership, 0, &no_of_components,
                                 IGRAPH_WEAK));

    /* Construct vertex and edge lists for each of the components */
    IGRAPH_CHECK(igraph_vector_ptr_init(&vertices_by_components, no_of_components));
    IGRAPH_CHECK(igraph_vector_ptr_init(&edges_by_components, no_of_components));
    IGRAPH_FINALLY(igraph_vector_ptr_destroy_all, &vertices_by_components);
    IGRAPH_FINALLY(igraph_vector_ptr_destroy_all, &edges_by_components);
    for (i = 0; i < no_of_components; i++) {
        igraph_vector_t* vptr;
        vptr = IGRAPH_CALLOC(1, igraph_vector_t);
        if (vptr == 0) {
            IGRAPH_ERROR("cannot calculate feedback arc set using IP", IGRAPH_ENOMEM);
        }
        IGRAPH_FINALLY(igraph_free, vptr);
        IGRAPH_CHECK(igraph_vector_init(vptr, 0));
        IGRAPH_FINALLY_CLEAN(1);
        VECTOR(vertices_by_components)[i] = vptr;
    }
    IGRAPH_VECTOR_PTR_SET_ITEM_DESTRUCTOR(&vertices_by_components, igraph_vector_destroy);
    for (i = 0; i < no_of_components; i++) {
        igraph_vector_t* vptr;
        vptr = IGRAPH_CALLOC(1, igraph_vector_t);
        if (vptr == 0) {
            IGRAPH_ERROR("cannot calculate feedback arc set using IP", IGRAPH_ENOMEM);
        }
        IGRAPH_FINALLY(igraph_free, vptr);
        IGRAPH_CHECK(igraph_vector_init(vptr, 0));
        IGRAPH_FINALLY_CLEAN(1);
        VECTOR(edges_by_components)[i] = vptr;
    }
    IGRAPH_VECTOR_PTR_SET_ITEM_DESTRUCTOR(&edges_by_components, igraph_vector_destroy);
    for (i = 0; i < no_of_vertices; i++) {
        j = (long int) VECTOR(membership)[i];
        IGRAPH_CHECK(igraph_vector_push_back(VECTOR(vertices_by_components)[j], i));
    }
    for (i = 0; i < no_of_edges; i++) {
        j = (long int) VECTOR(membership)[(long)IGRAPH_FROM(graph, i)];
        IGRAPH_CHECK(igraph_vector_push_back(VECTOR(edges_by_components)[j], i));
    }

#define VAR2IDX(i, j) (i*(n-1)+j-(i+1)*i/2)

    /* Configure GLPK */
    IGRAPH_GLPK_SETUP();
    glp_init_iocp(&parm);
    parm.br_tech = GLP_BR_DTH;
    parm.bt_tech = GLP_BT_BLB;
    parm.pp_tech = GLP_PP_ALL;
    parm.presolve = GLP_ON;
    parm.binarize = GLP_OFF;
    parm.cb_func = igraph_i_glpk_interruption_hook;

    /* Solve an IP for feedback arc sets in each of the components */
    for (i = 0; i < no_of_components; i++) {
        igraph_vector_t* vertices_in_comp = (igraph_vector_t*)VECTOR(vertices_by_components)[i];
        igraph_vector_t* edges_in_comp = (igraph_vector_t*)VECTOR(edges_by_components)[i];

        /*
         * Let x_ij denote whether layer(i) < layer(j).
         *
         * The standard formulation of the problem is as follows:
         *
         * max sum_{i,j} w_ij x_ij
         *
         * subject to
         *
         * (1) x_ij + x_ji = 1   (i.e. either layer(i) < layer(j) or layer(i) > layer(j))
         *     for all i < j
         * (2) x_ij + x_jk + x_ki <= 2 for all i < j, i < k, j != k
         *
         * Note that x_ij = 1 implies that x_ji = 0 and vice versa; in other words,
         * x_ij = 1 - x_ji. Thus, we can get rid of the (1) constraints and half of the
         * x_ij variables (where j < i) if we rewrite constraints of type (2) as follows:
         *
         * (2a) x_ij + x_jk - x_ik <= 1 for all i < j, i < k, j < k
         * (2b) x_ij - x_kj - x_ik <= 0 for all i < j, i < k, j > k
         *
         * The goal function then becomes:
         *
         * max sum_{i<j} (w_ij-w_ji) x_ij
         */
        n = igraph_vector_size(vertices_in_comp);
        ip = glp_create_prob();
        IGRAPH_FINALLY(igraph_i_glp_delete_prob, ip);
        glp_set_obj_dir(ip, GLP_MAX);

        /* Construct a mapping from vertex IDs to the [0; n-1] range */
        for (j = 0; j < n; j++) {
            VECTOR(vertex_remapping)[(long)VECTOR(*vertices_in_comp)[j]] = j;
        }

        /* Set up variables */
        k = n * (n - 1) / 2;
        if (k > 0) {
            glp_add_cols(ip, (int) k);
            for (j = 1; j <= k; j++) {
                glp_set_col_kind(ip, (int) j, GLP_BV);
            }
        }

        /* Set up coefficients in the goal function */
        k = igraph_vector_size(edges_in_comp);
        for (j = 0; j < k; j++) {
            l = (long int) VECTOR(*edges_in_comp)[j];
            from = (long int) VECTOR(vertex_remapping)[(long)IGRAPH_FROM(graph, l)];
            to = (long int) VECTOR(vertex_remapping)[(long)IGRAPH_TO(graph, l)];
            if (from == to) {
                continue;
            }

            weight = weights ? VECTOR(*weights)[l] : 1;

            if (from < to) {
                l = VAR2IDX(from, to);
                glp_set_obj_coef(ip, (int) l, glp_get_obj_coef(ip, (int) l) + weight);
            } else {
                l = VAR2IDX(to, from);
                glp_set_obj_coef(ip, (int) l, glp_get_obj_coef(ip, (int) l) - weight);
            }
        }

        /* Add constraints */
        if (n > 1) {
            glp_add_rows(ip, (int)(n * (n - 1) / 2 + n * (n - 1) * (n - 2) / 3));
            m = 1;
            for (j = 0; j < n; j++) {
                int ind[4];
                double val[4] = {0, 1, 1, -1};
                for (k = j + 1; k < n; k++) {
                    ind[1] = (int) VAR2IDX(j, k);
                    /* Type (2a) */
                    val[2] = 1;
                    for (l = k + 1; l < n; l++, m++) {
                        ind[2] = (int) VAR2IDX(k, l);
                        ind[3] = (int) VAR2IDX(j, l);
                        glp_set_row_bnds(ip, (int) m, GLP_UP, 1, 1);
                        glp_set_mat_row(ip, (int) m, 3, ind, val);
                    }
                    /* Type (2b) */
                    val[2] = -1;
                    for (l = j + 1; l < k; l++, m++) {
                        ind[2] = (int) VAR2IDX(l, k);
                        ind[3] = (int) VAR2IDX(j, l);
                        glp_set_row_bnds(ip, (int) m, GLP_UP, 0, 0);
                        glp_set_mat_row(ip, (int) m, 3, ind, val);
                    }
                }
            }
        }

        /* Solve the problem */
        IGRAPH_GLPK_CHECK(glp_intopt(ip, &parm), "Feedback arc set using IP failed");

        /* Find the ordering of the vertices */
        IGRAPH_CHECK(igraph_vector_resize(&ordering, n));
        igraph_vector_null(&ordering);
        m = n * (n - 1) / 2;
        j = 0; k = 1;
        for (l = 1; l <= m; l++) {
            /* variable l always corresponds to the (j, k) vertex pair */
            /* printf("(%ld, %ld) = %g\n", i, j, glp_mip_col_val(ip, l)); */
            if (glp_mip_col_val(ip, (int) l) > 0) {
                /* j comes earlier in the ordering than k */
                VECTOR(ordering)[j]++;
            } else {
                /* k comes earlier in the ordering than j */
                VECTOR(ordering)[k]++;
            }
            k++;
            if (k == n) {
                j++; k = j + 1;
            }
        }

        /* Find the feedback edges */
        k = igraph_vector_size(edges_in_comp);
        for (j = 0; j < k; j++) {
            l = (long int) VECTOR(*edges_in_comp)[j];
            from = (long int) VECTOR(vertex_remapping)[(long)IGRAPH_FROM(graph, l)];
            to = (long int) VECTOR(vertex_remapping)[(long)IGRAPH_TO(graph, l)];
            if (from == to || VECTOR(ordering)[from] < VECTOR(ordering)[to]) {
                IGRAPH_CHECK(igraph_vector_push_back(result, l));
            }
        }

        /* Clean up */
        glp_delete_prob(ip);
        IGRAPH_FINALLY_CLEAN(1);
    }

    igraph_vector_ptr_destroy_all(&vertices_by_components);
    igraph_vector_ptr_destroy_all(&edges_by_components);
    igraph_vector_destroy(&vertex_remapping);
    igraph_vector_destroy(&ordering);
    igraph_vector_destroy(&membership);
    IGRAPH_FINALLY_CLEAN(5);

    return IGRAPH_SUCCESS;
#endif
}
