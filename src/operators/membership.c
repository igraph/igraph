#include "igraph_operators.h"
#include "igraph_vector_list.h"
// ... other includes ...

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

    // Input validation
    if (vcount < 0){
        IGRAPH_ERRORF("Vertex count (%" IGRAPH_PRId ") must be non-negative.", 
            IGRAPH_EINVAL, vcount);
    }
    if (groups == NULL){
        IGRAPH_ERROR("Groups list must not be NULL.", IGRAPH_EINVAL);
    }
    if (membership == NULL){
        IGRAPH_ERROR("Membership vector must not be NULL.", IGRAPH_EINVAL);
    }

    IGRAPH_CHECK(igraph_vector_int_resize(membership, vcount));
    igraph_vector_int_fill(membership, -1); // Initialize with -1 to indicate no group

    // Track seen vertices
    igraph_vector_bool_t seen;
    IGRAPH_VECTOR_BOOL_INIT_FINALLY(&seen, vcount);
    igraph_vector_bool_fill(&seen, false);

    // Process each group
    igraph_int_t num_groups = igraph_vector_int_list_size(groups);
    igraph_int_t next_group_idx = 0;
    for (igraph_int_t i = 0; i < num_groups; i++){
        const igraph_vector_int_t *current_group = igraph_vector_int_list_get_ptr(groups, i);

        // Process each vertex in the group
        for (igraph_int_t j = 0; j < igraph_vector_int_size(current_group); j++){
            igraph_int_t vertex = VECTOR(*current_group)[j];

            // Validate vertex index and check for duplicates
            if (vertex < 0 || vertex >= vcount){
                IGRAPH_ERRORF("Invalid vertex index (%" IGRAPH_PRId ") in group %" IGRAPH_PRId ".", 
                    IGRAPH_EINVAL, vertex, i);
            }
            if (VECTOR(seen)[vertex]){
                IGRAPH_ERRORF("Vertex %" IGRAPH_PRId " appears multiple times in groups.", 
                    IGRAPH_EINVAL, vertex);
            }
            VECTOR(*membership)[vertex] = i;
            VECTOR(seen)[vertex] = true;
        }

        // Update next group index
        next_group_idx = i + 1;
    }

    // Process vertices that were not in any group (singletons)
    for (igraph_int_t i = 0; i < vcount; i++){
        if (VECTOR(*membership)[i] == -1){
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
 * \brief Convert membership vector to groups.
 * ...
 */
igraph_error_t igraph_membership_to_groups(
    const igraph_vector_int_t *membership,
    igraph_vector_int_list_t *groups) {
    IGRAPH_ERROR("Function not yet implemented.", IGRAPH_UNIMPLEMENTED);
}
