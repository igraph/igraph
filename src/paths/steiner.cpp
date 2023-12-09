/*
   IGraph library.
   Copyright (C) 2022  The igraph development team <igraph@igraph.org>

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
#include "igraph_paths.h"

#include "igraph_components.h"
#include "igraph_interface.h"
#include "igraph_matrix.h"
#include "igraph_structural.h"

#include "core/exceptions.h"
#include "core/interruption.h"

#include <cmath>
#include <algorithm>
#include <iterator>
#include <map>
#include <set>
#include <vector>

#include <iostream>
typedef std::set<igraph_integer_t> int_set;
typedef std::map<std::set<igraph_integer_t>, igraph_integer_t> dictionary;

/*
 *  Generates Subsets of length > 1 for list of integer values.
 *  As algorithm computes the minimum distance between set of variables ranging from
 *  length 2 to steiner_terminals - 1, we need these subsets.
 *  for e.g [2,3,4] --> [ [2,3], [2,4], [3,4], [2,3,4] ]
 *  Populates the subsetMap data structure as well.
 *
 */
static igraph_error_t generateSubsets(
        const igraph_vector_int_t *steinerTerminals, igraph_integer_t n, igraph_integer_t graphsize,
        dictionary &subsetMap, std::set<int_set> &allSubsets)
{
    if (n > sizeof(igraph_integer_t) * 8 - 2) {
        IGRAPH_ERROR(
                "igraph_integer_overflow detected. "
                "The given number of terminals is more than what the computer can handle.", IGRAPH_EINVAL);
    }
    igraph_integer_t count = ((igraph_integer_t)1 << n);
    igraph_integer_t subsetIndex = graphsize;

    /*
        As per algorithm the subsets that are explored have size atleast 2.
        The distance information between singleton vertices is already captured by the adjacency matrix.
    */
    for (igraph_integer_t i = 2; i < count; i++) {

        int_set newSubset;
        for (igraph_integer_t j = 0; j < n; j++) {

            if ((i & ((igraph_integer_t)1 << j)) > 0) {
                newSubset.insert(VECTOR(*steinerTerminals)[j]);
            }
        }

        if (newSubset.size() > 1) {
            if (allSubsets.find(newSubset) == allSubsets.end()) {
                allSubsets.insert(newSubset);
                subsetMap.insert(std::make_pair(newSubset, subsetIndex));
                subsetIndex++;
            }
        }
    }
    return IGRAPH_SUCCESS;
}

static void generateD_E(const int_set &D, std::set<int_set> &allSubsets, igraph_integer_t D1 = -1) {
    /*
    D1 is the element that is always going to be there in all the subsets that we generate
    D has elements for which we are creating subsets
    allSubsets is the container to store all
    */
    int n = D.size();
    igraph_integer_t count = ((igraph_integer_t)1 << n);

    for (igraph_integer_t i = 0; i < count - 1; i++) {
        int_set newSubset;
        if (D1 != (igraph_integer_t) -1) {
            newSubset.insert(D1);
        }
        for (igraph_integer_t j = 0; j < n; j++) {
            if ((i & ((igraph_integer_t)1 << j)) > 0) {
                auto it = D.begin();
                advance(it, j);
                auto jth_elem = *it;
                newSubset.insert(jth_elem);
            }
        }
        allSubsets.insert(newSubset);
    }
    return;
}

/*
 * Purpose: Fetching Index of a subset from subsetMap in order to store and look-up
 * the value of subset from DP table.
 */
static igraph_integer_t fetchIndexofMapofSets(const int_set &subset, const dictionary &subsetMap) {

    auto it = subsetMap.find(subset);
    if (it != subsetMap.end()) {
        return it->second;
    }
    IGRAPH_FATAL("The Subset's index that you tried to find doesn't exist. Hence the code won't run.");
}

/**
 * \function  igraph_steiner_dreyfus_wagner
 * \brief Finds a Steiner tree on an undirected graph using the Dreyfus-Wagner algorithm.
 *
 * A Steiner tree is a tree with minimum length connecting a set of terminal vertices
 * in a graph. If every vertex is a terminal then problem is reduced to minimum spanning tree.
 * The Steiner tree problem is NP-hard but since it's FPT, Dreyfus-Wagner algorithm is able to
 * find an exact Steiner tree on small graphs.
 *
 * </para><para>
 * The algorithm uses a dynamic programming approach and it treats the graph as undirected,
 * i.e. edge directions are ignored.
 *
 * </para><para>
 * Reference:
 *
 * </para><para>
 * S.E Dreyfus, R.A Wagner
 * The Steiner problem in graphs
 * Networks 1, 195-207 (1971).
 * https://doi.org/10.1002/net.3230010302
 *
 * \param graph The graph object.
 * \param weights The edge weights. All edge weights must be non-negative.
 *    Additionally, no edge weight may be NaN.
 * \param terminals Integer vector containing the IDs of the terminal vertices.
 * \param res Pointer to a real number, this will contain the result which is distance of Steiner Tree.
 * \param res_tree Pointer to an initialized integer vector, this will contain the IDs of edges
 *    that are part of Steiner tree.
 * \return Error code.
 *
 * Time complexity: O( 3^k ∗ V + 2^k ∗ V^2 + V∗(V+E) ∗ log(V) )
 * where V and E are the number of vertices and edges
 * and k is the number of Steiner terminals.
 * It is recommended that V <= 50 and k <= 11.
 *
 * \sa \ref igraph_minimum_spanning_tree(), \ref igraph_spanner()
 */

igraph_error_t igraph_steiner_dreyfus_wagner(
    const igraph_t *graph, const igraph_vector_int_t *terminals, const igraph_vector_t *weights,
    igraph_real_t *res, igraph_vector_int_t *res_tree) {
    IGRAPH_HANDLE_EXCEPTIONS_BEGIN;

    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_integer_t no_of_edges = igraph_ecount(graph);
    igraph_integer_t no_of_terminals = igraph_vector_int_size(terminals);

    if (!igraph_vector_int_isininterval(terminals, 0, no_of_nodes)) {
        IGRAPH_ERROR("Invalid vertex ID given as Steiner terminal.", IGRAPH_EINVVID);
    }

    if (weights && igraph_vector_size(weights) != no_of_edges) {

        IGRAPH_ERROR("Invalid weight vector length.", IGRAPH_EINVAL);
    }
    if (no_of_edges > 0 && weights) {
        igraph_real_t minweight = igraph_vector_min(weights);
        if (minweight < 0) {
            IGRAPH_ERRORF("Edge weights must be non-negative, got %g.", IGRAPH_EINVAL, minweight);
        } else if (minweight == 0) {
            /* TODO: can we support zero edge weights? */
            IGRAPH_ERROR("Weight vector contains zero weight.", IGRAPH_EINVAL);
        }
    }

    /* Handle the cases of the null graph and no terminals. */
    if (no_of_nodes == 0 || no_of_terminals == 0 || no_of_terminals == 1) {

        igraph_vector_int_clear(res_tree);
        *res = 0.0;
        return IGRAPH_SUCCESS;
    } else if (no_of_terminals == 2) {
        IGRAPH_CHECK(igraph_get_shortest_path_dijkstra(
                    graph, NULL, res_tree, VECTOR(*terminals)[0],
                    VECTOR(*terminals)[1], weights, IGRAPH_ALL));
        igraph_integer_t tree_size = igraph_vector_int_size(res_tree);

        *res = 0.0;
        if (weights) {
            for (igraph_integer_t i = 0; i < tree_size; i++) {
                *res += VECTOR(*weights)[VECTOR(*res_tree)[i]];
            }
        } else {
            *res = tree_size;
        }

        return IGRAPH_SUCCESS;
    }

    /* Check whether all terminals are within the same connected component. */
    {
        igraph_vector_int_t membership;
        igraph_integer_t no_comps;

        IGRAPH_VECTOR_INT_INIT_FINALLY(&membership, 0);
        IGRAPH_CHECK(igraph_connected_components(graph, &membership, NULL, &no_comps, IGRAPH_WEAK));

        if (no_comps > 1) {
            /* The case of zero terminals was already handled above. */
            igraph_integer_t component_id = VECTOR(membership)[VECTOR(*terminals)[0]];
            for (igraph_integer_t i = 1; i < no_of_terminals; i++) {
                if (VECTOR(membership)[VECTOR(*terminals)[i]] != component_id) {
                    IGRAPH_ERROR(
                            "Not all Steiner terminals are in the same connected component.",
                            IGRAPH_EINVAL);
                }
            }
        }

        igraph_vector_int_destroy(&membership);
        IGRAPH_FINALLY_CLEAN(1);
    }
    /* When every vertex is a Steiner terminal, the problem reduces to
     * finding a minimum spanning tree.
     */
    if (no_of_terminals == no_of_nodes) {
        IGRAPH_CHECK(igraph_minimum_spanning_tree(graph, res_tree, weights));

        igraph_real_t tree_weight = 0.0;
        igraph_integer_t tree_size = igraph_vector_int_size(res_tree);
        if (weights) {
            for (igraph_integer_t i = 0; i < tree_size; i++) {
                tree_weight += VECTOR(*weights)[VECTOR(*res_tree)[i]];
            }
        } else {
            tree_weight = tree_size;
        }
        *res = tree_weight;
        return IGRAPH_SUCCESS;
    }

    igraph_vector_int_t steiner_terminals_copy;
    igraph_matrix_t dp_cache; // dynamic programming table
    igraph_integer_t q;
    std::set<int_set> allSubsets = std::set<int_set>();
    igraph_matrix_t distance;

    IGRAPH_CHECK(igraph_matrix_init(&distance, no_of_nodes, no_of_nodes));
    IGRAPH_FINALLY(igraph_matrix_destroy, &distance);

    /*
     * Compute distances between all pairs of vertices.
     * The Dreyfus - Wagner algorithm needs complete graph information
     * hence this step is necessary.
     */
    IGRAPH_CHECK(igraph_distances_dijkstra(
                graph, &distance, igraph_vss_all(), igraph_vss_all(), weights, IGRAPH_ALL));

    IGRAPH_CHECK(igraph_vector_int_init_copy(&steiner_terminals_copy, terminals));
    IGRAPH_FINALLY(igraph_vector_int_destroy, &steiner_terminals_copy);
    igraph_vector_int_sort(&steiner_terminals_copy);

    q = VECTOR(steiner_terminals_copy)[0];

    igraph_vector_int_remove(&steiner_terminals_copy, 0);

    /*
     * DP table with size number of vertices in the graph + number of subsets
     * Number of subsets is 2 ^ (number of steiner_terminals_copy) - (number of steiner_terminals_copy + 1)
     */

    IGRAPH_CHECK(
            igraph_matrix_init(&dp_cache,
                no_of_nodes +
                    pow(2, igraph_vector_int_size(&steiner_terminals_copy))
                    - (igraph_vector_int_size(&steiner_terminals_copy) + 1),
                no_of_nodes));
    IGRAPH_FINALLY(igraph_matrix_destroy, &dp_cache);

    std::vector<std::vector<int_set> > distance_matrix =
        std::vector<std::vector<int_set> >
        (
         no_of_nodes + pow(2, igraph_vector_int_size(&steiner_terminals_copy)) -
         (igraph_vector_int_size(&steiner_terminals_copy) + 1),
         std::vector<int_set>(no_of_nodes,  int_set())
         );
    /* Addition to the dp_cache where at the same indices of the matrix we are storing
       the set of edge ids that are the shortest path to it
    */
    igraph_matrix_fill(&dp_cache, IGRAPH_INFINITY);
    /*
     * for singleton value the distance in dp cahce is just the
     * distance between same vertices in distance matrix
     * and the set of edges are just the edge between the vertices that we get with dijkstra algorithm
     */
    for (igraph_integer_t i = 0; i < no_of_nodes; i++) {
        IGRAPH_ALLOW_INTERRUPTION();
        for (igraph_integer_t j = 0; j < no_of_nodes; j++) {
            MATRIX(dp_cache, i, j) = MATRIX(distance, i, j);
            igraph_vector_int_t edgelist;
            IGRAPH_VECTOR_INT_INIT_FINALLY(&edgelist, 0);
            igraph_get_shortest_path_dijkstra(graph, nullptr, &edgelist, i, j, weights, IGRAPH_ALL);
            for (int k = 0 ; k < igraph_vector_int_size(&edgelist) ; k++) {
                distance_matrix[i][j].insert(igraph_vector_int_get(&edgelist, k));
            }
            igraph_vector_int_destroy(&edgelist);
            IGRAPH_FINALLY_CLEAN(1);
        }
    }

    /*
     * Usage: Data structure to store index value of subsets.
     * The index would be used in DP table.
     */
    dictionary subsetMap;
    IGRAPH_CHECK(
            generateSubsets(
                &steiner_terminals_copy,
                igraph_vector_int_size(&steiner_terminals_copy),
                no_of_nodes, subsetMap, allSubsets));

    int steiner_terminals_copy_size = igraph_vector_int_size(&steiner_terminals_copy);
    for (igraph_integer_t m = 2; m <= steiner_terminals_copy_size; m++) {
        for (igraph_integer_t i = 0; i < (igraph_integer_t)allSubsets.size(); i++) {
            auto it = allSubsets.begin();
            std::advance(it, i);
            int_set D = *it;

            if (D.size() != m) {
                continue;
            }

            igraph_integer_t indexOfSubsetD;
            indexOfSubsetD = fetchIndexofMapofSets(D, subsetMap);

            IGRAPH_ALLOW_INTERRUPTION();
            for (igraph_integer_t j = 0; j < no_of_nodes; j++) {
                MATRIX(dp_cache, indexOfSubsetD, j) = IGRAPH_INFINITY;
            }

            int_set D_prime = D;
            D_prime.erase(*D.begin());
            std::set<int_set> Subsets = std::set<int_set>();
            generateD_E(D_prime, Subsets, *D.begin());
            // E are subsets of D such that D[1] is in E and E is subset of D where E != D

            for (igraph_integer_t j = 0; j < no_of_nodes; j++) {
                igraph_real_t distance1 = IGRAPH_INFINITY;
                int_set set_u;
                for (auto E : Subsets) {
                    igraph_integer_t indexOfSubsetE =
                        (E.size() == 1)
                            ? *E.begin()
                            : fetchIndexofMapofSets(E, subsetMap);

                    igraph_real_t distanceEJ = MATRIX(dp_cache, indexOfSubsetE, j);

                    int_set DMinusE;

                    std::set_difference(
                            D.begin(), D.end(), E.begin(), E.end(),
                            std::inserter(DMinusE, DMinusE.end()));

                    igraph_integer_t indexOfSubsetDMinusE =
                        (DMinusE.size() == 1)
                            ? *DMinusE.begin()
                            : fetchIndexofMapofSets(DMinusE, subsetMap);

                    if ((distanceEJ + MATRIX(dp_cache, indexOfSubsetDMinusE, j)) < distance1) {
                        distance1 = distanceEJ + MATRIX(dp_cache, indexOfSubsetDMinusE, j);

                        //get the reference to this set and combine them
                        auto& setEJ = distance_matrix[indexOfSubsetE][j];
                        auto& setDMinusEJ = distance_matrix[indexOfSubsetDMinusE][j];
                        set_u.clear();
                        std::set_union(
                                setEJ.begin(), setEJ.end(), setDMinusEJ.begin(), setDMinusEJ.end(),
                                std::inserter(set_u, set_u.end()));
                    }
                }

                for (igraph_integer_t k = 0; k < no_of_nodes; k++) {
                    if ((MATRIX(distance, k, j) + distance1) < MATRIX(dp_cache, indexOfSubsetD, k)) {
                        MATRIX(dp_cache, indexOfSubsetD, k) = MATRIX(distance, k, j) + distance1;

                        //get the reference to this set and combine them
                        auto& setKJ = distance_matrix[k][j];
                        distance_matrix[indexOfSubsetD][k].clear();
                        std::set_union(
                                setKJ.begin(), setKJ.end(), set_u.begin(), set_u.end(),
                                std::inserter(
                                    distance_matrix[indexOfSubsetD][k],
                                    distance_matrix[indexOfSubsetD][k].end()));

                    }
                }
            }
        }
    }


    std::set<int_set> E_subsets = std::set<int_set>();
    int_set C_prime;
    int_set C;
    igraph_integer_t C_1 = igraph_vector_int_get(&steiner_terminals_copy, 0);
    for (int i = 0 ; i < igraph_vector_int_size(&steiner_terminals_copy) ; i++) {
        C.insert(igraph_vector_int_get(&steiner_terminals_copy, i));
        if (i != 0) {
            C_prime.insert(igraph_vector_int_get(&steiner_terminals_copy, i));
        }
    }
    // E are subsets of C such that C[1] is in E and E is subset of C where E != C
    generateD_E(C_prime, E_subsets, C_1);


    igraph_real_t distance2 = IGRAPH_INFINITY;
    int_set final_set;
    for (igraph_integer_t j = 0; j < no_of_nodes; j++) {
        igraph_real_t distance1 = IGRAPH_INFINITY;
        IGRAPH_ALLOW_INTERRUPTION();

        int_set set_u;

        for (auto E : E_subsets) {
            igraph_integer_t indexE =
                (E.size() == 1)
                ? *E.begin()
                : fetchIndexofMapofSets(E, subsetMap);

            int_set CMinusE;
            std::set_difference(
                    C.begin(), C.end(), E.begin(), E.end(),
                    std::inserter(CMinusE, CMinusE.end()));

            igraph_integer_t indexC_E =
                (CMinusE.size() == 1)
                ? *CMinusE.begin()
                : fetchIndexofMapofSets(CMinusE, subsetMap);

            igraph_real_t distanceEJ = MATRIX(dp_cache, indexE, j);
            if (((distanceEJ + (MATRIX(dp_cache, indexC_E, j))) < distance1)) {

                distance1 = distanceEJ + (MATRIX(dp_cache, indexC_E, j));

                //get the reference to this set and combine them
                set_u.clear();
                auto& setEJ = distance_matrix[indexE][j];
                auto& setCMinusEJ = distance_matrix[indexC_E][j];
                std::set_union(
                        setEJ.begin(), setEJ.end(), setCMinusEJ.begin(), setCMinusEJ.end(),
                        std::inserter(set_u, set_u.end()));
            }
        }

        if ((MATRIX(distance, q, j) + distance1) < distance2) {
            distance2 = MATRIX(distance, q, j) + distance1;
            //get the reference to this set and combine them
            final_set.clear();
            auto& setQJ = distance_matrix[q][j];
            std::set_union(
                    setQJ.begin(), setQJ.end(), set_u.begin(), set_u.end(),
                    std::inserter(final_set, final_set.end()));
        }
    }
    *res = distance2;
    //put the edges in the res_tree and move the final value to res,
    //which was stored in distance2 for the scope of this function
    igraph_vector_int_clear(res_tree);
    for (auto elem : final_set) {
        igraph_vector_int_push_back(res_tree, elem);
    }
    igraph_matrix_destroy(&distance);
    igraph_vector_int_destroy(&steiner_terminals_copy);
    igraph_matrix_destroy(&dp_cache);

    IGRAPH_FINALLY_CLEAN(3);

    IGRAPH_HANDLE_EXCEPTIONS_END;
    return IGRAPH_SUCCESS;
}
