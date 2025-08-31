/*
   igraph library.
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

static igraph_error_t community_label_propagation(
        const igraph_t *graph,
        igraph_vector_int_t *membership,
        igraph_neimode_t mode,
        const igraph_vector_t *weights,
        const igraph_vector_bool_t *fixed,
        igraph_bool_t retention) {

    const igraph_int_t no_of_nodes = igraph_vcount(graph);
    igraph_int_t no_of_not_fixed_nodes = 0;
    igraph_int_t i, j, k;
    igraph_adjlist_t al;
    igraph_inclist_t il;
    igraph_bool_t running, control_iteration;

    igraph_vector_t label_weights;
    igraph_vector_int_t dominant_labels, nonzero_labels, node_order;
    igraph_neimode_t reverse_mode;
    int iter = 0; /* interruption counter */

    reverse_mode = IGRAPH_REVERSE_MODE(mode);

    /* Create an adjacency/incidence list representation for efficiency.
    * For the unweighted case, the adjacency list is enough. For the
    * weighted case, we need the incidence list */
    if (weights) {
        IGRAPH_CHECK(igraph_inclist_init(graph, &il, reverse_mode, IGRAPH_LOOPS_ONCE));
        IGRAPH_FINALLY(igraph_inclist_destroy, &il);
    } else {
        IGRAPH_CHECK(igraph_adjlist_init(graph, &al, reverse_mode, IGRAPH_LOOPS_ONCE, IGRAPH_MULTIPLE));
        IGRAPH_FINALLY(igraph_adjlist_destroy, &al);
    }

    /* Create storage space for counting distinct labels and dominant ones */
    IGRAPH_VECTOR_INIT_FINALLY(&label_weights, no_of_nodes);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&dominant_labels, 0);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&nonzero_labels, 0);
    IGRAPH_CHECK(igraph_vector_int_reserve(&dominant_labels, 2));

    /* Initialize node ordering vector with only the not fixed nodes */
    if (fixed) {
        IGRAPH_VECTOR_INT_INIT_FINALLY(&node_order, no_of_nodes);
        for (i = 0; i < no_of_nodes; i++) {
            if (!VECTOR(*fixed)[i]) {
                VECTOR(node_order)[no_of_not_fixed_nodes] = i;
                no_of_not_fixed_nodes++;
            }
        }
        IGRAPH_CHECK(igraph_vector_int_resize(&node_order, no_of_not_fixed_nodes));
    } else {
        IGRAPH_CHECK(igraph_vector_int_init_range(&node_order, 0, no_of_nodes));
        IGRAPH_FINALLY(igraph_vector_int_destroy, &node_order);
        no_of_not_fixed_nodes = no_of_nodes;
    }

    /* There are two modes of operation in this implementation: retention or
     * dominance. When using retention, we prefer to keep the current label of a node.
     * Only if the current label is not among the dominant labels will we
     * update the label. If a label changes, we will continue to iterate
     * over all nodes.
     *
     * When not using retention we check for dominance after each iteration. This
     * is implemented as two alternating types of iterations, one for changing
     * labels and the other one for checking the end condition - every vertex in the
     * graph has a label to which the maximum number of its neighbors belongs. If
     * control_iteration is true, we are just checking the end condition and not
     * relabeling nodes.
     */

    control_iteration = true;
    running = true;
    while (running) {
        igraph_int_t v1, num_neis;
        igraph_real_t max_count;
        igraph_vector_int_t *neis;
        igraph_vector_int_t *ineis;
        igraph_bool_t was_zero;

        IGRAPH_ALLOW_INTERRUPTION_LIMITED(iter, 1 << 8);

        if (retention) {
            /* We stop in this iteration by default, unless a label changes */
            running = false;
            /* Shuffle the node ordering vector */
            igraph_vector_int_shuffle(&node_order);
        } else {
            if (control_iteration) {
                /* If we are in the control iteration, we expect in the beginning of
                the iteration that all vertices meet the end condition, so 'running' is false.
                If some of them does not, 'running' is set to true later in the code. */
                running = false;
            } else {
                /* Shuffle the node ordering vector if we are in the label updating iteration */
                igraph_vector_int_shuffle(&node_order);
            }
        }

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
                    if (k < 0) {
                        continue;    /* skip if it has no label yet */
                    }
                    was_zero = (VECTOR(label_weights)[k] == 0);
                    VECTOR(label_weights)[k] += VECTOR(*weights)[VECTOR(*ineis)[j]];

                    if (was_zero && VECTOR(label_weights)[k] != 0) {
                        /* weights just became nonzero */
                        IGRAPH_CHECK(igraph_vector_int_push_back(&nonzero_labels, k));
                    }

                    if (max_count < VECTOR(label_weights)[k]) {
                        max_count = VECTOR(label_weights)[k];
                        IGRAPH_CHECK(igraph_vector_int_resize(&dominant_labels, 1));
                        VECTOR(dominant_labels)[0] = k;
                    } else if (max_count == VECTOR(label_weights)[k]) {
                        IGRAPH_CHECK(igraph_vector_int_push_back(&dominant_labels, k));
                    }
                }
            } else {

                neis = igraph_adjlist_get(&al, v1);
                num_neis = igraph_vector_int_size(neis);

                for (j = 0; j < num_neis; j++) {

                    k = VECTOR(*membership)[VECTOR(*neis)[j]];
                    if (k < 0) {
                        continue;    /* skip if it has no label yet */
                    }
                    VECTOR(label_weights)[k]++;

                    if (VECTOR(label_weights)[k] == 1) {
                        /* weights just became nonzero */
                        IGRAPH_CHECK(igraph_vector_int_push_back(&nonzero_labels, k));
                    }

                    if (max_count < VECTOR(label_weights)[k]) {
                        max_count = VECTOR(label_weights)[k];
                        IGRAPH_CHECK(igraph_vector_int_resize(&dominant_labels, 1));
                        VECTOR(dominant_labels)[0] = k;
                    } else if (max_count == VECTOR(label_weights)[k]) {
                        IGRAPH_CHECK(igraph_vector_int_push_back(&dominant_labels, k));
                    }
                }
            }

            if (igraph_vector_int_size(&dominant_labels) > 0) {
                if (retention) {
                    /* If we are using retention, we first check if the current label
                    is among the maximum label. */
                    j = (long)VECTOR(*membership)[v1];
                    if (j < 0 || /* Doesn't have a label yet */
                        VECTOR(label_weights)[j] == 0 || /* Label not present in neighbors */
                        VECTOR(label_weights)[j] < max_count /* Label not dominant */) {
                        /* Select randomly from the dominant labels */
                        k = RNG_INTEGER(0, igraph_vector_int_size(&dominant_labels) - 1);
                        k = VECTOR(dominant_labels)[(long int)k];
                        /* If label changes, we will continue running */
                        if (k != j) {
                            running = true;
                        }
                        /* Actually change label */
                        VECTOR(*membership)[v1] = k;
                    }
                } else {
                    /* We are not using retention, so check if we should do a control iteration
                    or an update iteration. */
                    if (control_iteration) {
                        /* Check if the _current_ label of the node is also dominant */
                        k = VECTOR(*membership)[v1];
                        if (k < 0 ||                              /* No label assigned yet or */
                            VECTOR(label_weights)[k] < max_count /* Label is not maximum */
                           ) {
                            /* Nope, we need at least one more iteration */
                            running = true;
                        }
                    } else {
                        /* Select randomly from the dominant labels */
                        k = RNG_INTEGER(0, igraph_vector_int_size(&dominant_labels) - 1);
                        VECTOR(*membership)[v1] = VECTOR(dominant_labels)[k];
                    }
                }
            }

            /* Clear the nonzero elements in label_weights */
            num_neis = igraph_vector_int_size(&nonzero_labels);
            for (j = 0; j < num_neis; j++) {
                VECTOR(label_weights)[VECTOR(nonzero_labels)[j]] = 0;
            }
        }

        /* Alternating between control iterations and label updating iterations */
        if (!retention) {
            control_iteration = !control_iteration;
        }
    }

    if (weights) {
        igraph_inclist_destroy(&il);
    } else {
        igraph_adjlist_destroy(&al);
    }
    IGRAPH_FINALLY_CLEAN(1);

    igraph_vector_int_destroy(&node_order);
    igraph_vector_int_destroy(&nonzero_labels);
    igraph_vector_int_destroy(&dominant_labels);
    igraph_vector_destroy(&label_weights);
    IGRAPH_FINALLY_CLEAN(4);

    return IGRAPH_SUCCESS;
}

static igraph_error_t community_fast_label_propagation(
        const igraph_t *graph,
        igraph_vector_int_t *membership,
        igraph_neimode_t mode,
        const igraph_vector_t *weights,
        const igraph_vector_bool_t *fixed) {

    const igraph_int_t no_of_nodes = igraph_vcount(graph);
    igraph_int_t no_of_not_fixed_nodes = 0;
    igraph_int_t i, j, k;
    igraph_inclist_t il;
    igraph_adjlist_t al;

    igraph_vector_t label_weights;
    igraph_vector_int_t dominant_labels, nonzero_labels, node_order;
    igraph_dqueue_int_t queue;
    igraph_vector_bool_t in_queue;
    igraph_neimode_t reverse_mode;
    int iter = 0; /* interruption counter */

    reverse_mode = IGRAPH_REVERSE_MODE(mode);

    if (weights) {
        IGRAPH_CHECK(igraph_inclist_init(graph, &il, reverse_mode, IGRAPH_LOOPS_ONCE));
        IGRAPH_FINALLY(igraph_inclist_destroy, &il);
    } else {
        IGRAPH_CHECK(igraph_adjlist_init(graph, &al, reverse_mode, IGRAPH_LOOPS_ONCE, IGRAPH_MULTIPLE));
        IGRAPH_FINALLY(igraph_adjlist_destroy, &al);
    }

    /* Create storage space for counting distinct labels and dominant ones */
    IGRAPH_VECTOR_INIT_FINALLY(&label_weights, no_of_nodes);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&dominant_labels, 0);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&nonzero_labels, 0);
    IGRAPH_CHECK(igraph_vector_int_reserve(&dominant_labels, 2));

    /* Initialize node ordering vector with only the not fixed nodes */
    IGRAPH_DQUEUE_INT_INIT_FINALLY(&queue, no_of_nodes);
    IGRAPH_VECTOR_BOOL_INIT_FINALLY(&in_queue, no_of_nodes);

    /* Initialize node ordering vector with only the not fixed nodes */
    if (fixed) {
        IGRAPH_VECTOR_INT_INIT_FINALLY(&node_order, no_of_nodes);
        for (i = 0; i < no_of_nodes; i++) {
            if (!VECTOR(*fixed)[i]) {
                VECTOR(node_order)[no_of_not_fixed_nodes] = i;
                no_of_not_fixed_nodes++;
            }
        }
        IGRAPH_CHECK(igraph_vector_int_resize(&node_order, no_of_not_fixed_nodes));
    } else {
        IGRAPH_CHECK(igraph_vector_int_init_range(&node_order, 0, no_of_nodes));
        IGRAPH_FINALLY(igraph_vector_int_destroy, &node_order);
        no_of_not_fixed_nodes = no_of_nodes;
    }

    for (i = 0; i < no_of_not_fixed_nodes; i++) {
        IGRAPH_CHECK(igraph_dqueue_int_push(&queue, VECTOR(node_order)[i]));
        VECTOR(in_queue)[VECTOR(node_order)[i]] = true;
    }
    igraph_vector_int_destroy(&node_order);
    IGRAPH_FINALLY_CLEAN(1);

    while (!igraph_dqueue_int_empty(&queue)) {
        igraph_int_t v1, v2, e = -1, num_neis;
        igraph_real_t max_count;
        igraph_vector_int_t *neis;
        igraph_bool_t was_zero;

        IGRAPH_ALLOW_INTERRUPTION_LIMITED(iter, 1 << 8);

        v1 = igraph_dqueue_int_pop(&queue);
        VECTOR(in_queue)[v1] = false;

        /* Count the weights corresponding to different labels */
        igraph_vector_int_clear(&dominant_labels);
        igraph_vector_int_clear(&nonzero_labels);
        max_count = 0.0;
        if (weights) {
            neis = igraph_inclist_get(&il, v1);
        } else {
            neis = igraph_adjlist_get(&al, v1);
        }

        num_neis = igraph_vector_int_size(neis);
        for (j = 0; j < num_neis; j++) {
            if (weights) {
                e = VECTOR(*neis)[j];
                v2 = IGRAPH_OTHER(graph, e, v1);
            } else {
                v2 = VECTOR(*neis)[j];
            }
            k = VECTOR(*membership)[v2];
            if (k < 0) {
                continue;    /* skip if it has no label yet */
            }
            was_zero = (VECTOR(label_weights)[k] == 0);
            VECTOR(label_weights)[k] += (weights ? VECTOR(*weights)[e] : 1);

            if (was_zero && VECTOR(label_weights)[k] >= 0) {
                /* counter just became non-negative */
                IGRAPH_CHECK(igraph_vector_int_push_back(&nonzero_labels, k));
            }

            if (max_count < VECTOR(label_weights)[k]) {
                max_count = VECTOR(label_weights)[k];
                IGRAPH_CHECK(igraph_vector_int_resize(&dominant_labels, 1));
                VECTOR(dominant_labels)[0] = k;
            } else if (max_count == VECTOR(label_weights)[k]) {
                IGRAPH_CHECK(igraph_vector_int_push_back(&dominant_labels, k));
            }
        }

        if (igraph_vector_int_size(&dominant_labels) > 0) {
            igraph_int_t current_label = VECTOR(*membership)[v1];

            /* Select randomly from the dominant labels */
            k = RNG_INTEGER(0, igraph_vector_int_size(&dominant_labels) - 1);
            igraph_int_t new_label = VECTOR(dominant_labels)[k]; /* a dominant label */

            /* Check if the _current_ label of the node is not the same */
            if (new_label != current_label) {
                /* We still need to consider its neighbors not in the new community */
                for (j = 0; j < num_neis; j++) {
                    if (weights) {
                        e = VECTOR(*neis)[j];
                        v2 = IGRAPH_OTHER(graph, e, v1);
                    } else {
                        v2 = VECTOR(*neis)[j];
                    }
                    if (!VECTOR(in_queue)[v2]) {
                        igraph_int_t neigh_label = VECTOR(*membership)[v2]; /* neighbor community */
                        if (neigh_label != new_label && /* not in new community */
                            (fixed == NULL || !VECTOR(*fixed)[v2]) ) { /* not fixed */
                            IGRAPH_CHECK(igraph_dqueue_int_push(&queue, v2));
                            VECTOR(in_queue)[v2] = true;
                        }
                    }
                }
            }
            VECTOR(*membership)[v1] = new_label;
        }

        /* Clear the nonzero elements in label_weights */
        num_neis = igraph_vector_int_size(&nonzero_labels);
        for (j = 0; j < num_neis; j++) {
            VECTOR(label_weights)[VECTOR(nonzero_labels)[j]] = 0;
        }
    }

    if (weights) {
        igraph_inclist_destroy(&il);
    } else {
        igraph_adjlist_destroy(&al);
    }
    IGRAPH_FINALLY_CLEAN(1);

    igraph_vector_bool_destroy(&in_queue);
    igraph_dqueue_int_destroy(&queue);
    igraph_vector_destroy(&label_weights);
    igraph_vector_int_destroy(&dominant_labels);
    igraph_vector_int_destroy(&nonzero_labels);
    IGRAPH_FINALLY_CLEAN(5);

    return IGRAPH_SUCCESS;
}

/**
 * \ingroup communities
 * \function igraph_community_label_propagation
 * \brief Community detection based on label propagation.
 *
 * This function implements the label propagation-based community detection
 * algorithm described by Raghavan, Albert and Kumara (2007). This version extends
 * the original method by the ability to take edge weights into consideration
 * and also by allowing some labels to be fixed. In addition, it implements
 * the fast label propagation alternative introduced by Traag and Šubelj (2023).
 *
 * </para><para>
 * The algorithm works by iterating over nodes and updating the label of a node
 * based on the labels of its neighbors. The labels that are most frequent among
 * the neighbors are said to be dominant labels. The label of a node is always
 * updated to a dominant label. The algorithm guarantees that the label for each
 * is dominant when it terminates.
 *
 * </para><para>
 * There are several variants implemented, which work slightly differently with
 * the dominance of labels. Nodes with a dominant label might no longer have a
 * dominant label if one of their neighbors change label. In \c
 * IGRAPH_LPA_DOMINANCE an additional iteration over all nodes is made after
 * updating all labels to double check whether all nodes indeed have a dominant
 * label. When updating the label of a node, labels are always sampled from
 * among all dominant labels. The algorithm stops when all nodes have dominant
 * labels. In \c IGRAPH_LPA_RETENTION instead labels are only updated when they
 * are not dominant. That is, they retain their current label whenever the
 * current label is already dominant. The algorithm does not make an additional
 * iteration to check for dominance. Instead, it simply keeps track whether a
 * label has been updated, and terminates if no updates have been made. In \c
 * IGRAPH_LPA_FAST labels are sampled from among all dominant labels, similar to
 * \c IGRAPH_LPA_DOMINANCE. Instead of iterating over all nodes, it keeps track
 * of a queue of nodes that should be considered. Nodes are popped from the
 * queue when they are considered for update. When the label of a node is
 * updated, the node's neighbors are added to the queue again (if they weren't
 * already in the queue). The algorithm terminates when the queue is empty. All
 * variants guarantee that the labels for all nodes are dominant.
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
 * Raghavan, U.N. and Albert, R. and Kumara, S.: Near linear time algorithm to
 * detect community structures in large-scale networks. Phys Rev E 76, 036106
 * (2007). https://doi.org/10.1103/PhysRevE.76.036106
 *
 * </para><para>
 * Šubelj, L.: Label propagation for clustering. Chapter in "Advances in
 * Network Clustering and Blockmodeling" edited by P. Doreian, V. Batagelj
 * &amp; A. Ferligoj (Wiley, New York, 2018).
 * https://doi.org/10.1002/9781119483298.ch5
 * https://arxiv.org/abs/1709.05634
 *
 * </para><para>
 * Traag, V. A., and Šubelj, L.: Large network community detection by fast
 * label propagation. Scientific Reports, 13:1, (2023).
 * https://doi.org/10.1038/s41598-023-29610-z
 * https://arxiv.org/abs/2209.13338
 *
 * \param graph The input graph. Note that the algorithm was originally
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
 * \param lpa_variant Which variant of the label propagation algorithm to run.
 *   \clist
 *   \cli IGRAPH_LPA_DOMINANCE
 *     check for dominance of all nodes after each iteration.
 *   \cli IGRAPH_LPA_RETENTION
 *     keep current label if among dominant labels, only check if labels changed.
 *   \cli IGRAPH_LPA_FAST
 *     sample from dominant labels, only check neighbors.
 *   \endclist
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
        const igraph_vector_bool_t *fixed,
        igraph_lpa_variant_t lpa_variant) {

    const igraph_int_t no_of_nodes = igraph_vcount(graph);
    const igraph_int_t no_of_edges = igraph_ecount(graph);
    igraph_int_t no_of_not_fixed_nodes = no_of_nodes;
    igraph_int_t i, j, k;
    igraph_bool_t unlabelled_left;

    /* We make a copy of 'fixed' as a pointer into 'fixed_copy' after casting
     * away the constness, and promise ourselves that we will make a proper
     * copy of 'fixed' into 'fixed_copy' as soon as we start mutating it */
    igraph_vector_bool_t *fixed_copy = (igraph_vector_bool_t *) fixed;

    /* Unlabelled nodes are represented with -1. */
#define IS_UNLABELLED(x) (VECTOR(*membership)[x] < 0)

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
                VECTOR(*membership)[i] = -1;
            } else {
                VECTOR(*membership)[i] = VECTOR(*initial)[i];
            }
        }
        if (fixed) {
            for (i = 0; i < no_of_nodes; i++) {
                if (VECTOR(*fixed)[i]) {
                    if (IS_UNLABELLED(i)) {
                        IGRAPH_WARNING("Fixed nodes cannot be unlabeled, ignoring them.");

                        /* We cannot modify 'fixed' because it is const, so we make a copy and
                         * modify 'fixed_copy' instead */
                        if (fixed_copy == fixed) {
                            fixed_copy = IGRAPH_CALLOC(1, igraph_vector_bool_t);
                            IGRAPH_CHECK_OOM(fixed_copy, "Insufficient memory for label propagation.");
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
            VECTOR(*membership)[i] = i;
        }
    }

    /* From this point onwards we use 'fixed_copy' instead of 'fixed' */
    switch (lpa_variant) {
    case IGRAPH_LPA_FAST:
        IGRAPH_CHECK(community_fast_label_propagation(graph, membership, mode, weights, fixed_copy));
        break;

    case IGRAPH_LPA_RETENTION:
        IGRAPH_CHECK(community_label_propagation(graph, membership, mode, weights, fixed_copy, /* retention */ true ));
        break;

    case IGRAPH_LPA_DOMINANCE:
        IGRAPH_CHECK(community_label_propagation(graph, membership, mode, weights, fixed_copy, /* retention */ false));
        break;

    default:
        IGRAPH_ERROR("Invalid igraph_lpa_variant_t.", IGRAPH_EINVAL);
    }

    /* Permute labels in increasing order */
    igraph_vector_int_t relabel_label;
    IGRAPH_CHECK(igraph_vector_int_init(&relabel_label, no_of_nodes));
    igraph_vector_int_fill(&relabel_label, -1);
    IGRAPH_FINALLY(igraph_vector_int_destroy, &relabel_label);

    j = 0;
    unlabelled_left = false;
    for (i = 0; i < no_of_nodes; i++) {
        k = VECTOR(*membership)[i];
        if (k >= 0) {
            if (VECTOR(relabel_label)[k] == -1) {
                /* We have seen this label for the first time */
                VECTOR(relabel_label)[k] = j;
                k = j;
                j++;
            } else {
                k = VECTOR(relabel_label)[k];
            }
        } else {
            /* This is an unlabeled vertex */
            unlabelled_left = true;
        }
        VECTOR(*membership)[i] = k;
    }

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

        igraph_vector_int_t node_order;
        /* Initialize node ordering vector with only the not fixed nodes */
        if (fixed) {
            no_of_not_fixed_nodes = 0;
            IGRAPH_VECTOR_INT_INIT_FINALLY(&node_order, no_of_nodes);
            for (i = 0; i < no_of_nodes; i++) {
                if (!VECTOR(*fixed)[i]) {
                    VECTOR(node_order)[no_of_not_fixed_nodes] = i;
                    no_of_not_fixed_nodes++;
                }
            }
            IGRAPH_CHECK(igraph_vector_int_resize(&node_order, no_of_not_fixed_nodes));
        } else {
            IGRAPH_CHECK(igraph_vector_int_init_range(&node_order, 0, no_of_nodes));
            IGRAPH_FINALLY(igraph_vector_int_destroy, &node_order);
            no_of_not_fixed_nodes = no_of_nodes;
        }
        /* Shuffle the node ordering vector */
        igraph_vector_int_shuffle(&node_order);

        IGRAPH_VECTOR_INT_INIT_FINALLY(&neis, 0);

        IGRAPH_CHECK(igraph_dqueue_int_init(&q, 0));
        IGRAPH_FINALLY(igraph_dqueue_int_destroy, &q);

        for (i=0; i < no_of_not_fixed_nodes; ++i) {
            igraph_int_t v = VECTOR(node_order)[i];

            /* Is this node unlabelled? */
            if (IS_UNLABELLED(v)) {
                /* If yes, we label it, and do a BFS to apply the same label
                 * to all other unlabelled nodes reachable from it */
                IGRAPH_CHECK(igraph_dqueue_int_push(&q, v));
                VECTOR(*membership)[v] = j;
                while (!igraph_dqueue_int_empty(&q)) {
                    igraph_int_t ni, num_neis;
                    igraph_int_t actnode = igraph_dqueue_int_pop(&q);

                    IGRAPH_CHECK(igraph_neighbors(graph, &neis, actnode, mode, IGRAPH_LOOPS, IGRAPH_MULTIPLE));
                    num_neis = igraph_vector_int_size(&neis);

                    for (ni = 0; ni < num_neis; ++ni) {
                        igraph_int_t neighbor = VECTOR(neis)[ni];
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
        igraph_vector_int_destroy(&node_order);
        IGRAPH_FINALLY_CLEAN(3);
    }

    igraph_vector_int_destroy(&relabel_label);
    IGRAPH_FINALLY_CLEAN(1);

    if (fixed != fixed_copy) {
        igraph_vector_bool_destroy(fixed_copy);
        IGRAPH_FREE(fixed_copy);
        IGRAPH_FINALLY_CLEAN(2);
    }

    return IGRAPH_SUCCESS;
}
