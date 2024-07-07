/* -*- mode: C -*-  */
/* vim:set ts=4 sw=4 sts=4 et: */
/*
   IGraph library.
   Copyright (C) 2007-2022  The igraph development team <igraph@igraph.org>

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

#include "igraph_community.h"

#include "igraph_adjlist.h"
#include "igraph_dqueue.h"
#include "igraph_interface.h"
#include "igraph_memory.h"
#include "igraph_random.h"

#include "core/interruption.h"

/**
 * \ingroup communities
 * \function igraph_community_label_propagation
 * \brief Community detection based on label propagation.
 *
 * This function implements the label propagation-based community detection
 * algorithm described by Raghavan, Albert and Kumara. This version extends
 * the original method by the ability to take edge weights into consideration
 * and also by allowing some labels to be fixed.
 *
 * </para><para>
 * Weights are taken into account as follows: when the new label of node
 * \c i is determined, the algorithm iterates over all edges incident on
 * node \c i and calculate the total weight of edges leading to other
 * nodes with label 0, 1, 2, ..., \c k - 1 (where \c k is the number of possible
 * labels). The new label of node \c i will then be the label whose edges
 * (among the ones incident on node \c i) have the highest total weight.
 *
 * </para><para>
 * For directed graphs, it is important to know that labels can circulate
 * freely only within the strongly connected components of the graph and
 * may propagate in only one direction (or not at all) \em between strongly
 * connected components. You should treat directed edges as directed only
 * if you are aware of the consequences.
 *
 * </para><para>
 * References:
 *
 * </para><para>
 * Raghavan, U.N. and Albert, R. and Kumara, S.:
 * Near linear time algorithm to detect community structures in large-scale networks.
 * Phys Rev E 76, 036106 (2007).
 * https://doi.org/10.1103/PhysRevE.76.036106
 *
 * </para><para>
 * Šubelj, L.: Label propagation for clustering. Chapter in "Advances in
 * Network Clustering and Blockmodeling" edited by P. Doreian, V. Batagelj
 * &amp; A. Ferligoj (Wiley, New York, 2018).
 * https://doi.org/10.1002/9781119483298.ch5
 * https://arxiv.org/abs/1709.05634
 *
 * \param graph The input graph. Note that the algorithm wsa originally
 *    defined for undirected graphs. You are advised to set \p mode to
 *    \c IGRAPH_ALL if you pass a directed graph here to treat it as
 *    undirected.
 * \param membership The membership vector, the result is returned here.
 *    For each vertex it gives the ID of its community (label).
 * \param mode Whether to consider edge directions for the label propagation,
 *    and if so, which direction the labels should propagate. Ignored for
 *    undirected graphs. \c IGRAPH_ALL means to ignore edge directions (even
 *    in directed graphs). \c IGRAPH_OUT means to propagate labels along the
 *    natural direction of the edges. \c IGRAPH_IN means to propagate labels
 *    \em backwards (i.e. from head to tail). It is advised to set this to
 *    \c IGRAPH_ALL unless you are specifically interested in the effect of
 *    edge directions.
 * \param weights The weight vector, it should contain a positive
 *    weight for all the edges.
 * \param initial The initial state. If \c NULL, every vertex will have
 *   a different label at the beginning. Otherwise it must be a vector
 *   with an entry for each vertex. Non-negative values denote different
 *   labels, negative entries denote vertices without labels. Unlabeled
 *   vertices which are not reachable from any labeled ones will remain
 *   unlabeled at the end of the label propagation process, and will be
 *   labeled in an additional step to avoid returning negative values in
 *   \p membership. In undirected graphs, this happens when entire connected
 *   components are unlabeled. Then, each unlabeled component will receive
 *   its own separate label. In directed graphs, the outcome of the
 *   additional labeling should be considered undefined and may change
 *   in the future; please do not rely on it.
 * \param fixed Boolean vector denoting which labels are fixed. Of course
 *   this makes sense only if you provided an initial state, otherwise
 *   this element will be ignored. Note that vertices without labels
 *   cannot be fixed. The fixed status will be ignored for these with a
 *   warning. Also note that label numbers by themselves have no meaning,
 *   and igraph may renumber labels. However, co-membership constraints
 *   will be respected: two vertices can be fixed to be in the same or in
 *   different communities.
 * \param modularity If not a null pointer, then it must be a pointer
 *   to a real number. The modularity score of the detected community
 *   structure is stored here. Note that igraph will calculate the
 *   \em directed modularity if the input graph is directed, even if
 *   you set \p mode to \c IGRAPH_ALL
 * \return Error code.
 *
 * Time complexity: O(m+n)
 *
 * \example examples/simple/igraph_community_label_propagation.c
 */
igraph_error_t igraph_community_label_propagation(const igraph_t *graph,
                                       igraph_vector_int_t *membership,
                                       igraph_neimode_t mode,
                                       const igraph_vector_t *weights,
                                       const igraph_vector_int_t *initial,
                                       const igraph_vector_bool_t *fixed) {
    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_integer_t no_of_edges = igraph_ecount(graph);
    igraph_integer_t no_of_not_fixed_nodes = no_of_nodes;
    igraph_integer_t i, j, k;
    igraph_adjlist_t al;
    igraph_inclist_t il;
    igraph_bool_t running, control_iteration;
    igraph_bool_t unlabelled_left;
    igraph_neimode_t reversed_mode;
    int iter = 0; /* interruption counter */

    igraph_vector_t label_counters; /* real type, stores weight sums */
    igraph_vector_int_t dominant_labels, nonzero_labels, node_order;

    /* We make a copy of 'fixed' as a pointer into 'fixed_copy' after casting
     * away the constness, and promise ourselves that we will make a proper
     * copy of 'fixed' into 'fixed_copy' as soon as we start mutating it */
    igraph_vector_bool_t *fixed_copy = (igraph_vector_bool_t *) fixed;

    /* The implementation uses a trick to avoid negative array indexing:
     * elements of the membership vector are increased by 1 at the start
     * of the algorithm; this to allow us to denote unlabeled vertices
     * (if any) by zeroes. The membership vector is shifted back in the end
     */

    /* Do some initial checks */
    if (fixed && igraph_vector_bool_size(fixed) != no_of_nodes) {
        IGRAPH_ERROR("Fixed labeling vector length must agree with number of nodes.", IGRAPH_EINVAL);
    }
    if (weights) {
        if (igraph_vector_size(weights) != no_of_edges) {
            IGRAPH_ERROR("Length of weight vector must agree with number of edges.", IGRAPH_EINVAL);
        }
        if (no_of_edges > 0) {
            igraph_real_t minweight = igraph_vector_min(weights);
            if (minweight < 0) {
                IGRAPH_ERROR("Weights must not be negative.", IGRAPH_EINVAL);
            }
            if (isnan(minweight)) {
                IGRAPH_ERROR("Weights must not be NaN.", IGRAPH_EINVAL);
            }
        }
    }
    if (fixed && !initial) {
        IGRAPH_WARNING("Ignoring fixed vertices as no initial labeling given.");
    }

    IGRAPH_CHECK(igraph_vector_int_resize(membership, no_of_nodes));

    if (initial) {
        if (igraph_vector_int_size(initial) != no_of_nodes) {
            IGRAPH_ERROR("Initial labeling vector length must agree with number of nodes.", IGRAPH_EINVAL);
        }
        /* Check if the labels used are valid, initialize membership vector */
        for (i = 0; i < no_of_nodes; i++) {
            if (VECTOR(*initial)[i] < 0) {
                VECTOR(*membership)[i] = 0;
            } else {
                VECTOR(*membership)[i] = VECTOR(*initial)[i] + 1;
            }
        }
        if (fixed) {
            for (i = 0; i < no_of_nodes; i++) {
                if (VECTOR(*fixed)[i]) {
                    if (VECTOR(*membership)[i] == 0) {
                        IGRAPH_WARNING("Fixed nodes cannot be unlabeled, ignoring them.");

                        /* We cannot modify 'fixed' because it is const, so we make a copy and
                         * modify 'fixed_copy' instead */
                        if (fixed_copy == fixed) {
                            fixed_copy = IGRAPH_CALLOC(1, igraph_vector_bool_t);
                            if (fixed_copy == 0) {
                                IGRAPH_ERROR("Failed to copy 'fixed' vector.", IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */
                            }

                            IGRAPH_FINALLY(igraph_free, fixed_copy);
                            IGRAPH_CHECK(igraph_vector_bool_init_copy(fixed_copy, fixed));
                            IGRAPH_FINALLY(igraph_vector_bool_destroy, fixed_copy);
                        }

                        VECTOR(*fixed_copy)[i] = false;
                    } else {
                        no_of_not_fixed_nodes--;
                    }
                }
            }
        }

        i = igraph_vector_int_max(membership);
        if (i > no_of_nodes) {
            IGRAPH_ERROR("Elements of the initial labeling vector must be between 0 and |V|-1.", IGRAPH_EINVAL);
        }
    } else {
        for (i = 0; i < no_of_nodes; i++) {
            VECTOR(*membership)[i] = i + 1;
        }
    }

    reversed_mode = IGRAPH_REVERSE_MODE(mode);

    /* From this point onwards we use 'fixed_copy' instead of 'fixed' */

    /* Create an adjacency/incidence list representation for efficiency.
     * For the unweighted case, the adjacency list is enough. For the
     * weighted case, we need the incidence list */
    if (weights) {
        IGRAPH_CHECK(igraph_inclist_init(graph, &il, reversed_mode, IGRAPH_LOOPS_ONCE));
        IGRAPH_FINALLY(igraph_inclist_destroy, &il);
    } else {
        IGRAPH_CHECK(igraph_adjlist_init(graph, &al, reversed_mode, IGRAPH_LOOPS_ONCE, IGRAPH_MULTIPLE));
        IGRAPH_FINALLY(igraph_adjlist_destroy, &al);
    }

    /* Create storage space for counting distinct labels and dominant ones */
    IGRAPH_VECTOR_INIT_FINALLY(&label_counters, no_of_nodes + 1);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&dominant_labels, 0);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&nonzero_labels, 0);
    IGRAPH_CHECK(igraph_vector_int_reserve(&dominant_labels, 2));

    /* Initialize node ordering vector with only the not fixed nodes */
    if (fixed_copy) {
        IGRAPH_VECTOR_INT_INIT_FINALLY(&node_order, no_of_not_fixed_nodes);
        for (i = 0, j = 0; i < no_of_nodes; i++) {
            if (!VECTOR(*fixed_copy)[i]) {
                VECTOR(node_order)[j] = i;
                j++;
            }
        }
    } else {
        IGRAPH_CHECK(igraph_vector_int_init_range(&node_order, 0, no_of_nodes));
        IGRAPH_FINALLY(igraph_vector_int_destroy, &node_order);
    }

    /* There are two alternating types of iterations, one for changing labels and
    the other one for checking the end condition - every vertex in the graph has
    a label to which the maximum number of its neighbors belongs. If control_iteration
    is true, we are just checking the end condition and not relabeling nodes. */
    control_iteration = true;
    running = true;
    while (running) {
        igraph_integer_t v1, num_neis;
        igraph_real_t max_count;
        igraph_vector_int_t *neis;
        igraph_vector_int_t *ineis;
        igraph_bool_t was_zero;

        IGRAPH_ALLOW_INTERRUPTION_LIMITED(iter, 1 << 8);

        if (control_iteration) {
            /* If we are in the control iteration, we expect in the beginning of
            the iteration that all vertices meet the end condition, so 'running' is false.
            If some of them does not, 'running' is set to true later in the code. */
            running = false;
        } else {
            /* Shuffle the node ordering vector if we are in the label updating iteration */
            IGRAPH_CHECK(igraph_vector_int_shuffle(&node_order));
        }

        RNG_BEGIN();
        /* In the prescribed order, loop over the vertices and reassign labels */
        for (i = 0; i < no_of_not_fixed_nodes; i++) {
            v1 = VECTOR(node_order)[i];

            /* Count the weights corresponding to different labels */
            igraph_vector_int_clear(&dominant_labels);
            igraph_vector_int_clear(&nonzero_labels);
            max_count = 0.0;
            if (weights) {
                ineis = igraph_inclist_get(&il, v1);
                num_neis = igraph_vector_int_size(ineis);
                for (j = 0; j < num_neis; j++) {
                    k = VECTOR(*membership)[IGRAPH_OTHER(graph, VECTOR(*ineis)[j], v1)];
                    if (k == 0) {
                        continue;    /* skip if it has no label yet */
                    }
                    was_zero = (VECTOR(label_counters)[k] == 0);
                    VECTOR(label_counters)[k] += VECTOR(*weights)[VECTOR(*ineis)[j]];
                    if (was_zero && VECTOR(label_counters)[k] != 0) {
                        /* counter just became nonzero */
                        IGRAPH_CHECK(igraph_vector_int_push_back(&nonzero_labels, k));
                    }
                    if (max_count < VECTOR(label_counters)[k]) {
                        max_count = VECTOR(label_counters)[k];
                        IGRAPH_CHECK(igraph_vector_int_resize(&dominant_labels, 1));
                        VECTOR(dominant_labels)[0] = k;
                    } else if (max_count == VECTOR(label_counters)[k]) {
                        IGRAPH_CHECK(igraph_vector_int_push_back(&dominant_labels, k));
                    }
                }
            } else {
                neis = igraph_adjlist_get(&al, v1);
                num_neis = igraph_vector_int_size(neis);
                for (j = 0; j < num_neis; j++) {
                    k = VECTOR(*membership)[VECTOR(*neis)[j]];
                    if (k == 0) {
                        continue;    /* skip if it has no label yet */
                    }
                    VECTOR(label_counters)[k]++;
                    if (VECTOR(label_counters)[k] == 1) {
                        /* counter just became nonzero */
                        IGRAPH_CHECK(igraph_vector_int_push_back(&nonzero_labels, k));
                    }
                    if (max_count < VECTOR(label_counters)[k]) {
                        max_count = VECTOR(label_counters)[k];
                        IGRAPH_CHECK(igraph_vector_int_resize(&dominant_labels, 1));
                        VECTOR(dominant_labels)[0] = k;
                    } else if (max_count == VECTOR(label_counters)[k]) {
                        IGRAPH_CHECK(igraph_vector_int_push_back(&dominant_labels, k));
                    }
                }
            }

            if (igraph_vector_int_size(&dominant_labels) > 0) {
                if (control_iteration) {
                    /* Check if the _current_ label of the node is also dominant */
                    if (VECTOR(label_counters)[VECTOR(*membership)[v1]] != max_count) {
                        /* Nope, we need at least one more iteration */
                        running = true;
                    }
                }
                else {
                    /* Select randomly from the dominant labels */
                    k = RNG_INTEGER(0, igraph_vector_int_size(&dominant_labels) - 1);
                    VECTOR(*membership)[v1] = VECTOR(dominant_labels)[k];
                }
            }

            /* Clear the nonzero elements in label_counters */
            num_neis = igraph_vector_int_size(&nonzero_labels);
            for (j = 0; j < num_neis; j++) {
                VECTOR(label_counters)[VECTOR(nonzero_labels)[j]] = 0;
            }
        }
        RNG_END();

        /* Alternating between control iterations and label updating iterations */
        control_iteration = !control_iteration;
    }

    if (weights) {
        igraph_inclist_destroy(&il);
    } else {
        igraph_adjlist_destroy(&al);
    }
    IGRAPH_FINALLY_CLEAN(1);

    /* Shift back the membership vector, permute labels in increasing order */
    /* We recycle label_counters here :) and use it as an integer vector from now on */
    igraph_vector_fill(&label_counters, -1);
    j = 0;
    unlabelled_left = false;
    for (i = 0; i < no_of_nodes; i++) {
        k = VECTOR(*membership)[i] - 1;
        if (k >= 0) {
            if (VECTOR(label_counters)[k] == -1) {
                /* We have seen this label for the first time */
                VECTOR(label_counters)[k] = j;
                k = j;
                j++;
            } else {
                k = (igraph_integer_t) VECTOR(label_counters)[k];
            }
        } else {
            /* This is an unlabeled vertex */
            unlabelled_left = true;
        }
        VECTOR(*membership)[i] = k;
    }

    /* From this point on, unlabelled nodes are represented with -1 (no longer 0). */
#define IS_UNLABELLED(x) (VECTOR(*membership)[x] < 0)

    /* If any nodes are left unlabelled, we assign the remaining labels to them,
     * as well as to all unlabelled nodes reachable from them.
     *
     * Note that only those nodes could remain unlabelled which were unreachable
     * from any labelled ones. Thus, in the undirected case, fully unlabelled
     * connected components remain unlabelled. Here we label each such component
     * with the same label.
     */
    if (unlabelled_left) {
        igraph_dqueue_int_t q;
        igraph_vector_int_t neis;

        /* In the directed case, the outcome depends on the node ordering, thus we
         * shuffle nodes one more time. */
        IGRAPH_CHECK(igraph_vector_int_shuffle(&node_order));

        IGRAPH_VECTOR_INT_INIT_FINALLY(&neis, 0);

        IGRAPH_CHECK(igraph_dqueue_int_init(&q, 0));
        IGRAPH_FINALLY(igraph_dqueue_int_destroy, &q);

        for (i=0; i < no_of_not_fixed_nodes; ++i) {
            igraph_integer_t v = VECTOR(node_order)[i];

            /* Is this node unlabelled? */
            if (IS_UNLABELLED(v)) {
                /* If yes, we label it, and do a BFS to apply the same label
                 * to all other unlabelled nodes reachable from it */
                igraph_dqueue_int_push(&q, v);
                VECTOR(*membership)[v] = j;
                while (!igraph_dqueue_int_empty(&q)) {
                    igraph_integer_t ni, num_neis;
                    igraph_integer_t actnode = igraph_dqueue_int_pop(&q);

                    IGRAPH_CHECK(igraph_neighbors(graph, &neis, actnode, mode));
                    num_neis = igraph_vector_int_size(&neis);

                    for (ni = 0; ni < num_neis; ++ni) {
                        igraph_integer_t neighbor = VECTOR(neis)[ni];
                        if (IS_UNLABELLED(neighbor)) {
                            VECTOR(*membership)[neighbor] = j;
                            IGRAPH_CHECK(igraph_dqueue_int_push(&q, neighbor));
                        }
                    }
                }
                j++;
            }
        }

        igraph_vector_int_destroy(&neis);
        igraph_dqueue_int_destroy(&q);
        IGRAPH_FINALLY_CLEAN(2);
    }

    igraph_vector_int_destroy(&node_order);
    igraph_vector_destroy(&label_counters);
    igraph_vector_int_destroy(&dominant_labels);
    igraph_vector_int_destroy(&nonzero_labels);
    IGRAPH_FINALLY_CLEAN(4);

    if (fixed != fixed_copy) {
        igraph_vector_bool_destroy(fixed_copy);
        IGRAPH_FREE(fixed_copy);
        IGRAPH_FINALLY_CLEAN(2);
    }

    return IGRAPH_SUCCESS;
}
