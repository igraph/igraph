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
#include "igraph_memory.h"
#include "igraph_sparsemat.h"

#include <string.h>
#include <math.h>

/**
 * \function igraph_community_to_membership
 * \brief Creates a membership vector from a community structure dendrogram.
 *
 * This function creates a membership vector from a community
 * structure dendrogram. A membership vector contains for each vertex
 * the id of its graph component, the graph components are numbered
 * from zero, see the same argument of \ref igraph_connected_components()
 * for an example of a membership vector.
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
 *
 * \param merges The two-column matrix containing the merge
 *    operations. See \ref igraph_community_walktrap() for the
 *    detailed syntax.
 * \param nodes The number of leaf nodes in the dendrogram.
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
igraph_error_t igraph_community_to_membership(const igraph_matrix_int_t *merges,
                                   igraph_integer_t nodes,
                                   igraph_integer_t steps,
                                   igraph_vector_int_t *membership,
                                   igraph_vector_int_t *csize) {

    igraph_integer_t no_of_nodes = nodes;
    igraph_integer_t components = no_of_nodes - steps;
    igraph_integer_t i, found = 0;
    igraph_vector_int_t tmp;
    igraph_vector_bool_t already_merged;
    igraph_vector_int_t own_membership;
    igraph_bool_t using_own_membership = false;

    if (steps > igraph_matrix_int_nrow(merges)) {
        IGRAPH_ERRORF("Number of steps is greater than number of rows in merges matrix: found %"
                      IGRAPH_PRId " steps, %" IGRAPH_PRId " rows.", IGRAPH_EINVAL, steps, igraph_matrix_int_nrow(merges));
    }

    if (igraph_matrix_int_ncol(merges) != 2) {
        IGRAPH_ERRORF("The merges matrix should have two columns, but has %" IGRAPH_PRId ".",
                      IGRAPH_EINVAL, igraph_matrix_int_ncol(merges));
    }
    if (steps < 0) {
        IGRAPH_ERRORF("Number of steps should be non-negative, found %" IGRAPH_PRId ".", IGRAPH_EINVAL, steps);
    }

    if (csize != 0 && membership == 0) {
        /* we need a membership vector to calculate 'csize' but the user did
         * not provide one; let's allocate one ourselves */
        IGRAPH_VECTOR_INT_INIT_FINALLY(&own_membership, no_of_nodes);
        using_own_membership = true;
        membership = &own_membership;
    }

    if (membership) {
        IGRAPH_CHECK(igraph_vector_int_resize(membership, no_of_nodes));
        igraph_vector_int_null(membership);
    }
    if (csize) {
        IGRAPH_CHECK(igraph_vector_int_resize(csize, components));
        igraph_vector_int_null(csize);
    }

    IGRAPH_VECTOR_BOOL_INIT_FINALLY(&already_merged, steps + no_of_nodes);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&tmp, steps);

    for (i = steps - 1; i >= 0; i--) {
        igraph_integer_t c1 = MATRIX(*merges, i, 0);
        igraph_integer_t c2 = MATRIX(*merges, i, 1);

        if (VECTOR(already_merged)[c1] == 0) {
            VECTOR(already_merged)[c1] = true;
        } else {
            IGRAPH_ERRORF("Merges matrix contains multiple merges of cluster %" IGRAPH_PRId ".", IGRAPH_EINVAL, c1);
        }
        if (VECTOR(already_merged)[c2] == 0) {
            VECTOR(already_merged)[c2] = true;
        } else {
            IGRAPH_ERRORF("Merges matrix contains multiple merges of cluster %" IGRAPH_PRId ".", IGRAPH_EINVAL, c2);
        }

        /* new component? */
        if (VECTOR(tmp)[i] == 0) {
            found++;
            VECTOR(tmp)[i] = found;
        }

        if (c1 < no_of_nodes) {
            igraph_integer_t cid = VECTOR(tmp)[i] - 1;
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
            igraph_integer_t cid = VECTOR(tmp)[i] - 1;
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
        /* it can never happen that csize != 0 and membership == 0; we have
         * handled that case above */
        for (i = 0; i < no_of_nodes; i++) {
            igraph_integer_t tmp = VECTOR(*membership)[i];
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

    igraph_vector_int_destroy(&tmp);
    igraph_vector_bool_destroy(&already_merged);
    IGRAPH_FINALLY_CLEAN(2);

    if (using_own_membership) {
        igraph_vector_int_destroy(&own_membership);
        IGRAPH_FINALLY_CLEAN(1);
    }

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_reindex_membership
 * \brief Makes the IDs in a membership vector contiguous.
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
 *                     old component ID for each new one, or \c NULL,
 *                     in which case it is not returned. The vector
 *                     will be resized as needed.
 * \param  nb_clusters Pointer to an integer for the number of
 *                     distinct clusters. If not \c NULL, this will be
 *                     updated to reflect the number of distinct
 *                     clusters found in membership.
 *
 * Time complexity: should be O(n) for n elements.
 */
igraph_error_t igraph_reindex_membership(igraph_vector_int_t *membership,
                              igraph_vector_int_t *new_to_old,
                              igraph_integer_t *nb_clusters) {

    igraph_integer_t i, n = igraph_vector_int_size(membership);
    igraph_vector_t new_cluster;
    igraph_integer_t i_nb_clusters;

    /* We allow original cluster indices in the range 0, ..., n - 1 */
    IGRAPH_CHECK(igraph_vector_init(&new_cluster, n));
    IGRAPH_FINALLY(igraph_vector_destroy, &new_cluster);

    if (new_to_old) {
        igraph_vector_int_clear(new_to_old);
    }

    /* Clean clusters. We will store the new cluster + 1 so that membership == 0
     * indicates that no cluster was assigned yet. */
    i_nb_clusters = 1;
    for (i = 0; i < n; i++) {
        igraph_integer_t c = VECTOR(*membership)[i];

        if (c < 0) {
            IGRAPH_ERRORF("Membership indices should be non-negative. "
            "Found member of cluster %" IGRAPH_PRId ".", IGRAPH_EINVAL, c);
        }

        if (c >= n) {
            IGRAPH_ERRORF("Membership indices should be less than total number of vertices. "
            "Found member of cluster %" IGRAPH_PRId ", but only %" IGRAPH_PRId " vertices.", IGRAPH_EINVAL, c, n);
        }

        if (VECTOR(new_cluster)[c] == 0) {
            VECTOR(new_cluster)[c] = (igraph_real_t)i_nb_clusters;
            i_nb_clusters += 1;
            if (new_to_old) {
                IGRAPH_CHECK(igraph_vector_int_push_back(new_to_old, c));
            }
        }
    }

    /* Assign new membership */
    for (i = 0; i < n; i++) {
        igraph_integer_t c = VECTOR(*membership)[i];
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

static igraph_error_t igraph_i_compare_communities_vi(const igraph_vector_int_t *v1,
                                           const igraph_vector_int_t *v2, igraph_real_t* result);
static igraph_error_t igraph_i_compare_communities_nmi(const igraph_vector_int_t *v1,
                                            const igraph_vector_int_t *v2, igraph_real_t* result);
static igraph_error_t igraph_i_compare_communities_rand(const igraph_vector_int_t *v1,
                                             const igraph_vector_int_t *v2, igraph_real_t* result, igraph_bool_t adjust);
static igraph_error_t igraph_i_split_join_distance(const igraph_vector_int_t *v1,
                                        const igraph_vector_int_t *v2, igraph_integer_t* distance12,
                                        igraph_integer_t* distance21);

/**
 * \ingroup communities
 * \function igraph_compare_communities
 * \brief Compares community structures using various metrics.
 *
 * This function assesses the distance between two community structures
 * using the variation of information (VI) metric of Meila (2003), the
 * normalized mutual information (NMI) of Danon et al (2005), the
 * split-join distance of van Dongen (2000), the Rand index of Rand (1971)
 * or the adjusted Rand index of Hubert and Arabie (1985).
 *
 * </para><para>
 * Some of these measures are defined based on the entropy of a discrete
 * random variable associated with a given clustering \c C of vertices.
 * Let \c p_i be the probability that a randomly picked vertex would be part
 * of cluster \c i. Then the entropy of the clustering is
 *
 * </para><para>
 * <code>H(C) = - \sum_i p_i log p_i</code>
 *
 * </para><para>
 * Similarly, we can define the joint entropy of two clusterings \c C_1 and \c C_2
 * based on the probability \c p_ij that a random vertex is part of cluster \c i
 * in the first clustering and cluster \c j in the second one:
 *
 * </para><para>
 * <code>H(C_1, C_2) = - \sum_ii p_ij log p_ij</code>
 *
 * </para><para>
 * The mutual information of \c C_1 and \c C_2 is then
 * <code>MI(C_1, C_2) = H(C_1) + H(C_2) - H(C_1, C_2) >= 0 </code>.
 * A large mutual information indicates a high overlap between the two clusterings.
 * The normalized mutual information, as computed by igraph, is
 *
 * </para><para>
 * <code>NMI(C_1, C_2) = 2 MI(C_1, C_2) / (H(C_1) + H(C_2))</code>.
 *
 * </para><para>
 * It takes its value from the interval (0, 1], with 1 achieved when the two clusterings
 * coincide.
 *
 * </para><para>
 * The variation of information is defined as
 * <code>VI(C_1, C_2) = [H(C_1) - MI(C_1, C_2)] + [H(C_2) - MI(C_1, C_2)]</code>.
 * Lower values of the variation of information indicate a smaller difference between
 * the two clusterings, with <code>VI = 0</code> achieved precisely when they coincide.
 * igraph uses natural units for the variation of information, i.e. it uses the
 * natural logarithm when computing entropies.
 *
 * </para><para>
 * The Rand index is defined as the probability that the two clusterings agree
 * about the cluster memberships of a randomly chosen vertex \em pair. All vertex
 * pairs are considered, and the two clusterings are considered to be in agreement
 * about the memberships of a vertex pair if either the two vertices are in the
 * same cluster in both clusterings, or they are in different clusters in both
 * clusterings. The Rand index is then the number of vertex pairs in agreement,
 * divided by the total number of vertex pairs. A Rand index of zero means that
 * the two clusterings disagree about the membership of all vertex pairs, while
 * 1 means that the two clusterings are identical.
 *
 * </para><para>
 * The adjusted Rand index is similar to the Rand index, but it takes into
 * account that agreement between the two clusterings may also occur by chance
 * even if the two clusterings are chosen completely randomly. The adjusted
 * Rand index therefore subtracts the expected fraction of agreements from the
 * value of the Rand index, and divides the result by one minus the expected
 * fraction of agreements. The maximum value of the adjusted Rand index is
 * still 1 (similarly to the Rand index), indicating maximum agreement, but
 * the value may be less than zero if there is \em less agreement between the
 * two clusterings than what would be expected by chance.
 *
 * </para><para>
 * For an explanation of the split-join distance, see \ref igraph_split_join_distance().
 *
 * </para><para>
 * References:
 *
 * </para><para>
 * Meilă M: Comparing clusterings by the variation of information.
 * In: Schölkopf B, Warmuth MK (eds.). Learning Theory and Kernel Machines:
 * 16th Annual Conference on Computational Learning Theory and 7th Kernel
 * Workshop, COLT/Kernel 2003, Washington, DC, USA. Lecture Notes in Computer
 * Science, vol. 2777, Springer, 2003. ISBN: 978-3-540-40720-1.
 * https://doi.org/10.1007/978-3-540-45167-9_14
 *
 * </para><para>
 * Danon L, Diaz-Guilera A, Duch J, Arenas A: Comparing community structure
 * identification. J Stat Mech P09008, 2005.
 * https://doi.org/10.1088/1742-5468/2005/09/P09008
 *
 * </para><para>
 * van Dongen S: Performance criteria for graph clustering and Markov cluster
 * experiments. Technical Report INS-R0012, National Research Institute for
 * Mathematics and Computer Science in the Netherlands, Amsterdam, May 2000.
 * https://ir.cwi.nl/pub/4461
 *
 * </para><para>
 * Rand WM: Objective criteria for the evaluation of clustering methods.
 * J Am Stat Assoc 66(336):846-850, 1971.
 * https://doi.org/10.2307/2284239
 *
 * </para><para>
 * Hubert L and Arabie P: Comparing partitions. Journal of Classification
 * 2:193-218, 1985.
 * https://doi.org/10.1007/BF01908075
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
 * \sa \ref igraph_split_join_distance().
 *
 * Time complexity: O(n log(n)).
 */
igraph_error_t igraph_compare_communities(const igraph_vector_int_t *comm1,
                               const igraph_vector_int_t *comm2, igraph_real_t* result,
                               igraph_community_comparison_t method) {
    igraph_vector_int_t c1, c2;

    if (igraph_vector_int_size(comm1) != igraph_vector_int_size(comm2)) {
        IGRAPH_ERROR("community membership vectors have different lengths", IGRAPH_EINVAL);
    }

    /* Copy and reindex membership vectors to make sure they are continuous */
    IGRAPH_CHECK(igraph_vector_int_init_copy(&c1, comm1));
    IGRAPH_FINALLY(igraph_vector_int_destroy, &c1);

    IGRAPH_CHECK(igraph_vector_int_init_copy(&c2, comm2));
    IGRAPH_FINALLY(igraph_vector_int_destroy, &c2);

    IGRAPH_CHECK(igraph_reindex_membership(&c1, NULL, NULL));
    IGRAPH_CHECK(igraph_reindex_membership(&c2, NULL, NULL));

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
    igraph_vector_int_destroy(&c1);
    igraph_vector_int_destroy(&c2);
    IGRAPH_FINALLY_CLEAN(2);

    return IGRAPH_SUCCESS;
}

/**
 * \ingroup communities
 * \function igraph_split_join_distance
 * \brief Calculates the split-join distance of two community structures.
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
 * \sa \ref igraph_compare_communities() with the \c IGRAPH_COMMCMP_SPLIT_JOIN
 * method if you are not interested in the individual distances but only the sum
 * of them.
 *
 * Time complexity: O(n log(n)).
 */
igraph_error_t igraph_split_join_distance(const igraph_vector_int_t *comm1,
                               const igraph_vector_int_t *comm2, igraph_integer_t *distance12,
                               igraph_integer_t *distance21) {
    igraph_vector_int_t c1, c2;

    if (igraph_vector_int_size(comm1) != igraph_vector_int_size(comm2)) {
        IGRAPH_ERRORF("Community membership vectors have different lengths: %" IGRAPH_PRId " and %" IGRAPH_PRId ".",
                      IGRAPH_EINVAL, igraph_vector_int_size(comm1), igraph_vector_int_size(comm2));
    }

    /* Copy and reindex membership vectors to make sure they are continuous */
    IGRAPH_CHECK(igraph_vector_int_init_copy(&c1, comm1));
    IGRAPH_FINALLY(igraph_vector_int_destroy, &c1);

    IGRAPH_CHECK(igraph_vector_int_init_copy(&c2, comm2));
    IGRAPH_FINALLY(igraph_vector_int_destroy, &c2);

    IGRAPH_CHECK(igraph_reindex_membership(&c1, NULL, NULL));
    IGRAPH_CHECK(igraph_reindex_membership(&c2, NULL, NULL));

    IGRAPH_CHECK(igraph_i_split_join_distance(&c1, &c2, distance12, distance21));

    /* Clean up everything */
    igraph_vector_int_destroy(&c1);
    igraph_vector_int_destroy(&c2);
    IGRAPH_FINALLY_CLEAN(2);

    return IGRAPH_SUCCESS;
}

/**
 * Calculates the entropy and the mutual information for two reindexed community
 * membership vectors v1 and v2. This is needed by both Meila's and Danon's
 * community comparison measure.
 */
static igraph_error_t igraph_i_entropy_and_mutual_information(const igraph_vector_int_t* v1,
        const igraph_vector_int_t* v2, double* h1, double* h2, double* mut_inf) {
    igraph_integer_t i, n;
    igraph_integer_t k1;
    igraph_integer_t k2;
    igraph_real_t *p1, *p2;
    igraph_sparsemat_t m;
    igraph_sparsemat_t mu; /* uncompressed */
    igraph_sparsemat_iterator_t mit;

    n = igraph_vector_int_size(v1);
    if (n == 0) {
        *h1 = 0;
        *h2 = 0;
        *mut_inf = 0;
        return IGRAPH_SUCCESS;
    }
    k1 = igraph_vector_int_max(v1) + 1;
    k2 = igraph_vector_int_max(v2) + 1;
    p1 = IGRAPH_CALLOC(k1, igraph_real_t);
    if (p1 == 0) {
        IGRAPH_ERROR("Insufficient memory for computing community entropy.", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, p1);
    p2 = IGRAPH_CALLOC(k2, igraph_real_t);
    if (p2 == 0) {
        IGRAPH_ERROR("Insufficient memory for computing community entropy.", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, p2);

    /* Calculate the entropy of v1 */
    *h1 = 0.0;
    for (i = 0; i < n; i++) {
        p1[VECTOR(*v1)[i]]++;
    }
    for (i = 0; i < k1; i++) {
        p1[i] /= n;
        *h1 -= p1[i] * log(p1[i]);
    }

    /* Calculate the entropy of v2 */
    *h2 = 0.0;
    for (i = 0; i < n; i++) {
        p2[VECTOR(*v2)[i]]++;
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
    IGRAPH_CHECK(igraph_sparsemat_init(&mu, k1, k2, n));
    IGRAPH_FINALLY(igraph_sparsemat_destroy, &mu);
    for (i = 0; i < n; i++) {
        IGRAPH_CHECK(igraph_sparsemat_entry(
            &mu, VECTOR(*v1)[i],
            VECTOR(*v2)[i], 1
        ));
    }

    IGRAPH_CHECK(igraph_sparsemat_compress(&mu, &m));
    IGRAPH_FINALLY(igraph_sparsemat_destroy, &m);
    IGRAPH_CHECK(igraph_sparsemat_dupl(&m));

    IGRAPH_CHECK(igraph_sparsemat_iterator_init(&mit, &m));
    while (!igraph_sparsemat_iterator_end(&mit)) {
        double p = igraph_sparsemat_iterator_get(&mit)/ n;
        *mut_inf += p * (log(p) - p1[igraph_sparsemat_iterator_row(&mit)] - p2[igraph_sparsemat_iterator_col(&mit)]);
        igraph_sparsemat_iterator_next(&mit);
    }
    igraph_sparsemat_destroy(&m);
    igraph_sparsemat_destroy(&mu);
    IGRAPH_FREE(p1); IGRAPH_FREE(p2);

    IGRAPH_FINALLY_CLEAN(4);

    return IGRAPH_SUCCESS;
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
static igraph_error_t igraph_i_compare_communities_nmi(const igraph_vector_int_t *v1, const igraph_vector_int_t *v2,
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
static igraph_error_t igraph_i_compare_communities_vi(const igraph_vector_int_t *v1, const igraph_vector_int_t *v2,
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
static igraph_error_t igraph_i_confusion_matrix(const igraph_vector_int_t *v1, const igraph_vector_int_t *v2,
                              igraph_sparsemat_t *m) {
    igraph_integer_t k1, k2, i, n;

    n = igraph_vector_int_size(v1);
    if (n == 0) {
        IGRAPH_CHECK(igraph_sparsemat_resize(m, 0, 0, 0));
        return IGRAPH_SUCCESS;
    }

    k1 = igraph_vector_int_max(v1) + 1;
    k2 = igraph_vector_int_max(v2) + 1;
    IGRAPH_CHECK(igraph_sparsemat_resize(m, k1, k2, n));
    for (i = 0; i < n; i++) {
        IGRAPH_CHECK(igraph_sparsemat_entry(
            m, VECTOR(*v1)[i], VECTOR(*v2)[i], 1
        ));
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
static igraph_error_t igraph_i_split_join_distance(const igraph_vector_int_t *v1, const igraph_vector_int_t *v2,
                                 igraph_integer_t* distance12, igraph_integer_t* distance21) {
    igraph_integer_t n = igraph_vector_int_size(v1);
    igraph_vector_t rowmax, colmax;
    igraph_sparsemat_t m;
    igraph_sparsemat_t mu; /* uncompressed */
    igraph_sparsemat_iterator_t mit;

    if (n == 0) {
        *distance12 = 0;
        *distance21 = 0;
        return IGRAPH_SUCCESS;
    }
    /* Calculate the confusion matrix */
    IGRAPH_CHECK(igraph_sparsemat_init(&mu, 1, 1, 0));
    IGRAPH_FINALLY(igraph_sparsemat_destroy, &mu);
    IGRAPH_CHECK(igraph_i_confusion_matrix(v1, v2, &mu));

    /* Initialize vectors that will store the row/columnwise maxima */
    IGRAPH_VECTOR_INIT_FINALLY(&rowmax, igraph_sparsemat_nrow(&mu));
    IGRAPH_VECTOR_INIT_FINALLY(&colmax, igraph_sparsemat_ncol(&mu));

    /* Find the row/columnwise maxima */
    igraph_sparsemat_compress(&mu, &m);
    IGRAPH_FINALLY(igraph_sparsemat_destroy, &m);
    IGRAPH_CHECK(igraph_sparsemat_dupl(&m));
    IGRAPH_CHECK(igraph_sparsemat_iterator_init(&mit, &m));
    while (!igraph_sparsemat_iterator_end(&mit)) {
        igraph_real_t value = igraph_sparsemat_iterator_get(&mit);
        igraph_integer_t row = igraph_sparsemat_iterator_row(&mit);
        igraph_integer_t col = igraph_sparsemat_iterator_col(&mit);
        if (value > VECTOR(rowmax)[row]) {
            VECTOR(rowmax)[row] = value;
        }
        if (value > VECTOR(colmax)[col]) {
            VECTOR(colmax)[col] = value;
        }
        igraph_sparsemat_iterator_next(&mit);
    }

    /* Calculate the distances */
    *distance12 = (igraph_integer_t) (n - igraph_vector_sum(&rowmax));
    *distance21 = (igraph_integer_t) (n - igraph_vector_sum(&colmax));

    igraph_vector_destroy(&rowmax);
    igraph_vector_destroy(&colmax);
    igraph_sparsemat_destroy(&m);
    igraph_sparsemat_destroy(&mu);
    IGRAPH_FINALLY_CLEAN(4);

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
static igraph_error_t igraph_i_compare_communities_rand(
        const igraph_vector_int_t *v1, const igraph_vector_int_t *v2,
        igraph_real_t *result, igraph_bool_t adjust) {
    igraph_sparsemat_t m;
    igraph_sparsemat_t mu; /* uncompressed */
    igraph_sparsemat_iterator_t mit;
    igraph_vector_t rowsums, colsums;
    igraph_integer_t i, nrow, ncol;
    igraph_real_t rand, n;
    igraph_real_t frac_pairs_in_1, frac_pairs_in_2;

    if (igraph_vector_int_size(v1) <= 1) {
        IGRAPH_ERRORF("Rand indices not defined for only zero or one vertices. "
        "Found membership vector of size %" IGRAPH_PRId ".", IGRAPH_EINVAL, igraph_vector_int_size(v1));
    }

    /* Calculate the confusion matrix */
    IGRAPH_CHECK(igraph_sparsemat_init(&mu, 1, 1, 0));
    IGRAPH_FINALLY(igraph_sparsemat_destroy, &mu);
    IGRAPH_CHECK(igraph_i_confusion_matrix(v1, v2, &mu));

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
    nrow = igraph_sparsemat_nrow(&mu);
    ncol = igraph_sparsemat_ncol(&mu);
    n = igraph_vector_int_size(v1);
    IGRAPH_VECTOR_INIT_FINALLY(&rowsums, nrow);
    IGRAPH_VECTOR_INIT_FINALLY(&colsums, ncol);
    IGRAPH_CHECK(igraph_sparsemat_rowsums(&mu, &rowsums));
    IGRAPH_CHECK(igraph_sparsemat_colsums(&mu, &colsums));

    /* Start calculating the unadjusted Rand index */
    rand = 0.0;
    igraph_sparsemat_compress(&mu, &m);
    IGRAPH_FINALLY(igraph_sparsemat_destroy, &m);
    IGRAPH_CHECK(igraph_sparsemat_dupl(&m));

    IGRAPH_CHECK(igraph_sparsemat_iterator_init(&mit, &m));
    while (!igraph_sparsemat_iterator_end(&mit)) {
        igraph_real_t value = igraph_sparsemat_iterator_get(&mit);
        rand += (value / n) * (value - 1) / (n - 1);
        igraph_sparsemat_iterator_next(&mit);
    }

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
    igraph_sparsemat_destroy(&m);
    igraph_sparsemat_destroy(&mu);
    IGRAPH_FINALLY_CLEAN(4);

    *result = rand;

    return IGRAPH_SUCCESS;
}
