/*
   igraph library.
   Copyright (C) 2025  The igraph development team <igraph@igraph.org>

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

#include "igraph_operators.h"
#include "igraph_vector_list.h"

/**
 * \function igraph_groups_to_membership
 * \brief Convert a list of groups to a membership vector.
 * \experimental
 *
 * This function creates a membership vector from a list of vertex groups.
 * Each vertex is assigned to a group index. Vertices that appear in the
 * groups list are assigned to their respective group indices (0, 1, 2, ...).
 *
 * </para><para>
 * Vertices that do not appear in any group are assigned to singleton groups,
 * with group indices starting from the number of explicit groups.
 *
 * </para><para>
 * The function validates that no vertex appears in multiple groups and that
 * all vertex indices are within the valid range [0, vcount).
 *
 * \param vcount The total number of vertices. Must be non-negative.
 * \param groups A list of vectors, where each vector contains the vertex
 *        indices belonging to that group.
 * \param membership Pointer to an initialized vector. The membership vector
 *        will be resized to \p vcount and filled with group indices.
 * \return Error code:
 *         \clist
 *           \cli IGRAPH_EINVAL
 *                Invalid vertex index, negative vertex count, or a vertex
 *                appears in multiple groups.
 *         \endclist
 *
 * Time complexity: O(|V| + sum of group sizes), linear in the number of
 * vertices plus the total number of elements in all groups.
 */
igraph_error_t igraph_groups_to_membership(
    igraph_integer_t vcount,
    const igraph_vector_int_list_t *groups,
    igraph_vector_int_t *membership) {

    /* Input validation */
    if (vcount < 0) {
        IGRAPH_ERRORF("Vertex count (%" IGRAPH_PRId ") must be non-negative.",
                      IGRAPH_EINVAL, vcount);
    }
    if (groups == NULL) {
        IGRAPH_ERROR("Groups list must not be NULL.", IGRAPH_EINVAL);
    }
    if (membership == NULL) {
        IGRAPH_ERROR("Membership vector must not be NULL.", IGRAPH_EINVAL);
    }

    IGRAPH_CHECK(igraph_vector_int_resize(membership, vcount));
    igraph_vector_int_fill(membership, -1); /* Initialize with -1 to indicate no group */

    /* Track seen vertices */
    igraph_vector_bool_t seen;
    IGRAPH_VECTOR_BOOL_INIT_FINALLY(&seen, vcount);
    igraph_vector_bool_fill(&seen, false);

    /* Process each group */
    igraph_int_t num_groups = igraph_vector_int_list_size(groups);
    igraph_int_t next_group_idx = num_groups;
    for (igraph_int_t i = 0; i < num_groups; i++) {
        const igraph_vector_int_t *current_group = igraph_vector_int_list_get_ptr(groups, i);

        /* Process each vertex in the group */
        for (igraph_int_t j = 0; j < igraph_vector_int_size(current_group); j++) {
            igraph_int_t vertex = VECTOR(*current_group)[j];

            /* Validate vertex index and check for duplicates */
            if (vertex < 0 || vertex >= vcount) {
                IGRAPH_ERRORF("Invalid vertex index (%" IGRAPH_PRId ") in group %" IGRAPH_PRId ".",
                              IGRAPH_EINVAL, vertex, i);
            }
            if (VECTOR(seen)[vertex]) {
                IGRAPH_ERRORF("Vertex %" IGRAPH_PRId " appears multiple times in groups.",
                              IGRAPH_EINVAL, vertex);
            }
            VECTOR(*membership)[vertex] = i;
            VECTOR(seen)[vertex] = true;
        }
    }

    /* Process vertices that were not in any group (singletons) */
    for (igraph_int_t i = 0; i < vcount; i++) {
        if (VECTOR(*membership)[i] == -1) {
            VECTOR(*membership)[i] = next_group_idx;
            next_group_idx++;
        }
    }

    igraph_vector_bool_destroy(&seen);
    IGRAPH_FINALLY_CLEAN(1);
    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_membership_to_groups
 * \brief Convert a membership vector to a list of groups.
 * \experimental
 *
 * This function creates a list of groups from a membership vector.
 * Each group will contain the vertex indices that have the same
 * membership value.
 *
 * </para><para>
 * The groups are created in order of increasing membership value,
 * and each vertex index appears in exactly one group corresponding
 * to its membership value.
 *
 * \param membership A membership vector, where each element
 *        indicates the group index of the corresponding vertex.
 * \param groups Pointer to an initialized vector list. The groups
 *        will be stored here, with each vector containing the
 *        vertex indices belonging to that group.
 * \return Error code:
 *         \clist
 *           \cli IGRAPH_EINVAL
 *                Invalid membership value (negative or out of bounds).
 *         \endclist
 *
 * Time complexity: O(|V| + |G|), where |V| is the number of vertices
 * and |G| is the number of distinct groups (bounded by |V|).
 *
 * \sa \ref igraph_groups_to_membership() for the inverse operation.
 */
igraph_error_t igraph_membership_to_groups(
    const igraph_vector_int_t *membership,
    igraph_vector_int_list_t *groups) {

    /* Input validation */
    if (membership == NULL) {
        IGRAPH_ERROR("Membership vector must not be NULL.", IGRAPH_EINVAL);
    }
    if (groups == NULL) {
        IGRAPH_ERROR("Groups list must not be NULL.", IGRAPH_EINVAL);
    }

    /* Return empty list for groups if vertex count is zero */
    igraph_int_t vcount = igraph_vector_int_size(membership);
    if (vcount == 0) {
        IGRAPH_CHECK(igraph_vector_int_list_resize(groups, 0));
        return IGRAPH_SUCCESS;
    }

    /* Initialize groups list */
    igraph_int_t num_groups = igraph_vector_int_max(membership) + 1;
    IGRAPH_CHECK(igraph_vector_int_list_resize(groups, num_groups));

    /* Count and allocate groups */
    igraph_vector_int_t group_sizes;
    IGRAPH_VECTOR_INT_INIT_FINALLY(&group_sizes, num_groups);
    igraph_vector_int_fill(&group_sizes, 0);

    for (igraph_int_t i = 0; i < vcount; i++) {
        igraph_int_t group_idx = VECTOR(*membership)[i];
        if (group_idx < 0 || group_idx >= num_groups) {
            IGRAPH_ERRORF("Invalid membership value (%" IGRAPH_PRId ") for vertex %" IGRAPH_PRId ".",
                          IGRAPH_EINVAL, group_idx, i);
        }
        VECTOR(group_sizes)[group_idx]++;
    }

    /* Pre-allocate space for each group */
    for (igraph_int_t i = 0; i < num_groups; i++) {
        igraph_vector_int_t *group = igraph_vector_int_list_get_ptr(groups, i);
        IGRAPH_CHECK(igraph_vector_int_reserve(group, VECTOR(group_sizes)[i]));
    }

    /* Populate each group */
    for (igraph_int_t i = 0; i < vcount; i++) {
        igraph_int_t group_idx = VECTOR(*membership)[i];
        igraph_vector_int_t *group = igraph_vector_int_list_get_ptr(groups, group_idx);
        IGRAPH_CHECK(igraph_vector_int_push_back(group, i));
    }

    IGRAPH_VECTOR_INT_DESTROY(&group_sizes);
    IGRAPH_FINALLY_CLEAN(1);
    return IGRAPH_SUCCESS;
}
