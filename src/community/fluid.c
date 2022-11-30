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
#include "igraph_interface.h"
#include "igraph_random.h"
#include "igraph_structural.h"

/**
 * \ingroup communities
 * \function igraph_community_fluid_communities
 * \brief Community detection based on fluids interacting on the graph.
 *
 * The algorithm is based on the simple idea of
 * several fluids interacting in a non-homogeneous environment
 * (the graph topology), expanding and contracting based on their
 * interaction and density. Weighted graphs are not supported.
 *
 * </para><para>
 * This function implements the community detection method described in:
 * Parés F, Gasulla DG, et. al. (2018) Fluid Communities: A Competitive,
 * Scalable and Diverse Community Detection Algorithm. In: Complex Networks
 * &amp; Their Applications VI: Proceedings of Complex Networks 2017 (The Sixth
 * International Conference on Complex Networks and Their Applications),
 * Springer, vol 689, p 229. https://doi.org/10.1007/978-3-319-72150-7_19
 *
 * \param graph The input graph. The graph must be simple and connected.
 *   Edge directions will be ignored.
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
 */
igraph_error_t igraph_community_fluid_communities(const igraph_t *graph,
                                       igraph_integer_t no_of_communities,
                                       igraph_vector_int_t *membership) {
    /* Declaration of variables */
    igraph_integer_t no_of_nodes, i, j, k, kv1;
    igraph_adjlist_t al;
    igraph_real_t max_density;
    igraph_bool_t is_simple, is_connected, running;
    igraph_vector_t density, label_counters;
    igraph_vector_int_t dominant_labels, node_order, com_to_numvertices;

    /* Initialization of variables needed for initial checking */
    no_of_nodes = igraph_vcount(graph);

    /* Checking input values */
    if (no_of_nodes < 2) {
        if (membership) {
            IGRAPH_CHECK(igraph_vector_int_resize(membership, no_of_nodes));
            igraph_vector_int_fill(membership, 0);
        }
        return IGRAPH_SUCCESS;
    }
    if (no_of_communities < 1) {
        IGRAPH_ERROR("Number of requested communities must be greater than zero.", IGRAPH_EINVAL);
    }
    if (no_of_communities > no_of_nodes) {
        IGRAPH_ERROR("Number of requested communities must not be greater than the number of nodes.",
                     IGRAPH_EINVAL);
    }
    IGRAPH_CHECK(igraph_is_simple(graph, &is_simple));
    if (!is_simple) {
        IGRAPH_ERROR("Fluid community detection supports only simple graphs.", IGRAPH_EINVAL);
    }
    if (igraph_is_directed(graph)) {
        /* When the graph is directed, mutual edges are effectively multi-edges as we
         * are ignoring edge directions. */
        igraph_bool_t has_mutual;
        IGRAPH_CHECK(igraph_has_mutual(graph, &has_mutual, false));
        if (has_mutual) {
            IGRAPH_ERROR("Fluid community detection supports only simple graphs.", IGRAPH_EINVAL);
        }
    }
    IGRAPH_CHECK(igraph_is_connected(graph, &is_connected, IGRAPH_WEAK));
    if (!is_connected) {
        IGRAPH_ERROR("Fluid community detection supports only connected graphs.", IGRAPH_EINVAL);
    }
    if (igraph_is_directed(graph)) {
        IGRAPH_WARNING("Edge directions are ignored by fluid community detection.");
    }

    /* Internal variables initialization */
    max_density = 1.0;

    /* Resize membership vector (number of nodes) */
    IGRAPH_CHECK(igraph_vector_int_resize(membership, no_of_nodes));

    /* Initialize density and com_to_numvertices vectors */
    IGRAPH_CHECK(igraph_vector_init(&density, no_of_communities));
    IGRAPH_FINALLY(igraph_vector_destroy, &density);
    IGRAPH_CHECK(igraph_vector_int_init(&com_to_numvertices, no_of_communities));
    IGRAPH_FINALLY(igraph_vector_int_destroy, &com_to_numvertices);

    /* Initialize node ordering vector */
    IGRAPH_CHECK(igraph_vector_int_init_range(&node_order, 0, no_of_nodes));
    IGRAPH_FINALLY(igraph_vector_int_destroy, &node_order);

    /* Initialize the membership vector with 0 values */
    igraph_vector_int_null(membership);
    /* Initialize densities to max_density */
    igraph_vector_fill(&density, max_density);

    /* Initialize com_to_numvertices and initialize communities into membership vector */
    IGRAPH_CHECK(igraph_vector_int_shuffle(&node_order));
    for (i = 0; i < no_of_communities; i++) {
        /* Initialize membership at initial nodes for each community
         * where 0 refers to have no label*/
        VECTOR(*membership)[VECTOR(node_order)[i]] = i + 1;
        /* Initialize com_to_numvertices list: Number of vertices for each community */
        VECTOR(com_to_numvertices)[i] = 1;
    }

    /* Create an adjacency list representation for efficiency. */
    IGRAPH_CHECK(igraph_adjlist_init(graph, &al, IGRAPH_ALL, IGRAPH_LOOPS_TWICE, IGRAPH_MULTIPLE));
    IGRAPH_FINALLY(igraph_adjlist_destroy, &al);

    /* Create storage space for counting distinct labels and dominant ones */
    IGRAPH_VECTOR_INT_INIT_FINALLY(&dominant_labels, no_of_communities);

    IGRAPH_CHECK(igraph_vector_init(&label_counters, no_of_communities));
    IGRAPH_FINALLY(igraph_vector_destroy, &label_counters);

    RNG_BEGIN();

    /* running is the convergence boolean variable */
    running = true;
    while (running) {
        /* Declarations of variables used inside main loop */
        igraph_integer_t v1, size, rand_idx;
        igraph_real_t max_count, label_counter_diff;
        igraph_vector_int_t *neis;
        igraph_bool_t same_label_in_dominant;

        running = false;

        /* Shuffle the node ordering vector */
        IGRAPH_CHECK(igraph_vector_int_shuffle(&node_order));
        /* In the prescribed order, loop over the vertices and reassign labels */
        for (i = 0; i < no_of_nodes; i++) {
            /* Clear dominant_labels and nonzero_labels vectors */
            igraph_vector_int_clear(&dominant_labels);
            igraph_vector_null(&label_counters);

            /* Obtain actual node index */
            v1 = VECTOR(node_order)[i];
            /* Take into account same label in updating rule */
            kv1 = VECTOR(*membership)[v1];
            max_count = 0.0;
            if (kv1 != 0) {
                VECTOR(label_counters)[kv1 - 1] += VECTOR(density)[kv1 - 1];
                /* Set up max_count */
                max_count = VECTOR(density)[kv1 - 1];
                /* Initialize dominant_labels */
                IGRAPH_CHECK(igraph_vector_int_resize(&dominant_labels, 1));
                VECTOR(dominant_labels)[0] = kv1;
            }

            /* Count the weights corresponding to different labels */
            neis = igraph_adjlist_get(&al, v1);
            size = igraph_vector_int_size(neis);
            for (j = 0; j < size; j++) {
                k = VECTOR(*membership)[VECTOR(*neis)[j]];
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
                    IGRAPH_CHECK(igraph_vector_int_resize(&dominant_labels, 1));
                    VECTOR(dominant_labels)[0] = k;
                } else if (-0.0001 < label_counter_diff && label_counter_diff < 0.0001) {
                    IGRAPH_CHECK(igraph_vector_int_push_back(&dominant_labels, k));
                }
            }

            if (!igraph_vector_int_empty(&dominant_labels)) {
                /* Maintain same label if it exists in dominant_labels */
                same_label_in_dominant = igraph_vector_int_contains(&dominant_labels, kv1);

                if (!same_label_in_dominant) {
                    /* We need at least one more iteration */
                    running = true;

                    /* Select randomly from the dominant labels */
                    rand_idx = RNG_INTEGER(0, igraph_vector_int_size(&dominant_labels) - 1);
                    k = VECTOR(dominant_labels)[rand_idx];

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
        IGRAPH_ASSERT(VECTOR(*membership)[i] >= 0); /* all vertices must have a community assigned */
    }

    igraph_adjlist_destroy(&al);
    IGRAPH_FINALLY_CLEAN(1);

    igraph_vector_int_destroy(&node_order);
    igraph_vector_destroy(&density);
    igraph_vector_int_destroy(&com_to_numvertices);
    igraph_vector_destroy(&label_counters);
    igraph_vector_int_destroy(&dominant_labels);
    IGRAPH_FINALLY_CLEAN(5);

    return IGRAPH_SUCCESS;
}
