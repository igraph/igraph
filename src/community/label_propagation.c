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
#include "igraph_interface.h"
#include "igraph_random.h"

/**
 * \ingroup communities
 * \function igraph_community_label_propagation
 * \brief Community detection based on label propagation.
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
 * \c i is determined, the algorithm iterates over all edges incident on
 * node \c i and calculate the total weight of edges leading to other
 * nodes with label 0, 1, 2, ..., \c k - 1 (where \c k is the number of possible
 * labels). The new label of node \c i will then be the label whose edges
 * (among the ones incident on node \c i) have the highest total weight.
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
 *   cannot be fixed. If they are, this vector will be modified to
 *   make it consistent with \p initial.
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
            if (igraph_is_nan(minweight)) {
                IGRAPH_ERROR("Weights must not be NaN.", IGRAPH_EINVAL);
            }
        }
    }
    if (fixed && !initial) {
        IGRAPH_WARNING("Ignoring fixed vertices as no initial labeling given.");
    }

    IGRAPH_CHECK(igraph_vector_resize(membership, no_of_nodes));

    if (initial) {
        if (igraph_vector_size(initial) != no_of_nodes) {
            IGRAPH_ERROR("Initial labeling vector length must agree with number of nodes.", IGRAPH_EINVAL);
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
                        IGRAPH_WARNING("Fixed nodes cannot be unlabeled, ignoring them.");
                        VECTOR(*fixed)[i] = 0;
                    } else {
                        no_of_not_fixed_nodes--;
                    }
                }
            }
        }

        i = (long int) igraph_vector_max(membership);
        if (i > no_of_nodes) {
            IGRAPH_ERROR("Elements of the initial labeling vector must be between 0 and |V|-1.", IGRAPH_EINVAL);
        }
        if (i <= 0) {
            IGRAPH_ERROR("At least one vertex must be labeled in the initial labeling.", IGRAPH_EINVAL);
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
        IGRAPH_CHECK(igraph_inclist_init(graph, &il, IGRAPH_IN, IGRAPH_LOOPS_ONCE));
        IGRAPH_FINALLY(igraph_inclist_destroy, &il);
    } else {
        IGRAPH_CHECK(igraph_adjlist_init(graph, &al, IGRAPH_IN, IGRAPH_LOOPS_ONCE, IGRAPH_MULTIPLE));
        IGRAPH_FINALLY(igraph_adjlist_destroy, &al);
    }

    /* Create storage space for counting distinct labels and dominant ones */
    IGRAPH_VECTOR_INIT_FINALLY(&label_counters, no_of_nodes + 1);
    IGRAPH_VECTOR_INIT_FINALLY(&dominant_labels, 0);
    IGRAPH_VECTOR_INIT_FINALLY(&nonzero_labels, 0);
    IGRAPH_CHECK(igraph_vector_reserve(&dominant_labels, 2));

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

        RNG_BEGIN();
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
        RNG_END();
    }

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
      IGRAPH_CHECK(igraph_modularity(graph, membership, weights,
                                     /* resolution */ 1,
                                     /* directed */ 1, modularity));
    }

    igraph_vector_destroy(&node_order);
    igraph_vector_destroy(&label_counters);
    igraph_vector_destroy(&dominant_labels);
    igraph_vector_destroy(&nonzero_labels);
    IGRAPH_FINALLY_CLEAN(4);

    return IGRAPH_SUCCESS;
}
