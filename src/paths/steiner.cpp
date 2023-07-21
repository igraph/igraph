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
static igraph_error_t generateSubsets(const igraph_vector_int_t *steinerTerminals, igraph_integer_t n, igraph_integer_t graphsize, dictionary &subsetMap, std::set<int_set> &allSubsets) {
    if (n > sizeof(igraph_integer_t) * 8 - 2) {
        IGRAPH_ERROR("igraph_integer_overflow detected. The given number of terminals is more than what the computer can handle.", IGRAPH_EINVAL);
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

static void generateD_E(int D1, const int_set &D, std::set<int_set> &allSubsets) {
    int n = D.size();
    igraph_integer_t count = ((igraph_integer_t)1 << n);

    for (igraph_integer_t i = 0; i < count - 1; i++) {
        int_set newSubset;
        newSubset.insert(D1);
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

/*
 * Purpose: Retriving the value of subset from given index.
 */
static int_set fetchSetsBasedonIndex(igraph_integer_t index, const dictionary &subsetMap) {
    for (const auto &kv : subsetMap) {
        if (kv.second == index) {
            return kv.first;
        }
    }
    IGRAPH_FATAL("The index that you tried to find doesn't exist. Hence the code won't run.");
}

/*
 * Calculating factorial of a number.
 */
static igraph_integer_t factorial(igraph_integer_t n) {
    igraph_integer_t answer = 1;
    for (igraph_integer_t i = 1; i <= n; i++) {
        answer *= i;
    }
    return answer;
}

/*
 * Finding number of combinations nCr
 * Used to determine number of elements to be scanned in DP table
 * for a particular subset during generation of Steiner Tree.
 */
static igraph_integer_t Combination(igraph_integer_t n, igraph_integer_t r) {
    return factorial(n) / (factorial(n - r) * factorial(r));
}

/**
 * Iterates through a particular row of the DP table
 * find the value of column where distance (q,k) + distance(Sk(D))
 * where D is a subset. This signifies a vertex that is part
 * of minimum path from vertex q to D
 *
 * \param dp_cache The DP table.
 * \param indexD Index of the subset D.
 * \param q vertex from who minimum path needs tp be calculated
 */
static igraph_integer_t findMinimumK(igraph_matrix_t *dp_cache, igraph_integer_t indexD, igraph_integer_t q, const dictionary &subsetMap) {

    igraph_integer_t min_col_num = -1;
    igraph_real_t min_sum_for_col;

    int_set D = fetchSetsBasedonIndex(indexD, subsetMap);

    for (igraph_integer_t i = 0; i < dp_cache->ncol; i++) {
        if (q != i) {

            if (min_col_num == -1) {
                min_col_num = i;
                min_sum_for_col = MATRIX(*dp_cache, q, i) + MATRIX(*dp_cache, indexD, i);
            } else if (MATRIX(*dp_cache, q, i) + MATRIX(*dp_cache, indexD, i) < min_sum_for_col) {
                min_col_num = i;
                min_sum_for_col = MATRIX(*dp_cache, q, i) + MATRIX(*dp_cache, indexD, i);
            } else if ((MATRIX(*dp_cache, q, i) + MATRIX(*dp_cache, indexD, i) == min_sum_for_col) && (D.find(i) != D.end())) {
                min_sum_for_col = MATRIX(*dp_cache, q, i) + MATRIX(*dp_cache, indexD, i);
                min_col_num = q;
            }
        }
    }

    return min_col_num;
}

/**
 * \function generate_steiner_tree_exact
 *
 * Generation of the Steiner Tree based on calculation of minimum distances in DP table.
 *
 * \param graph The graph object.
 * \param weights The edge weights. All edge weights must be
 *                non-negative. Additionally, no edge weight may be NaN.
 * \param dp_cache The DP table.
 * \param indexD The index of subset D in DP table.
 * \param q The vertex that was remmoved from steiner terminals.
 * \param vertexlist_all The vector to capture vertices in resultant Steiner Tree.
 * \param edgelist_all The vector to capture edges in resultant Steiner Tree.
 */
static igraph_error_t generate_steiner_tree_exact(const igraph_t *graph, const igraph_vector_t *weights,
        igraph_matrix_t *dp_cache, igraph_integer_t indexD, igraph_integer_t q, igraph_vector_int_t *edgelist_all, const dictionary &subsetMap) {

    int_set C = fetchSetsBasedonIndex(indexD, subsetMap);
    // Initially the value of m is the vertex that was removed from Steiner Terminals
    igraph_integer_t m = q;
    int_set D = C;
    std::set<igraph_integer_t> edgelist_all_set;

    while (D.size() > 1) {

        std::set<igraph_integer_t>::iterator itr;


        indexD = fetchIndexofMapofSets(D, subsetMap);

        // Finding the bridge vertex from m to subset D which is part of shortest path.
        igraph_integer_t k = findMinimumK(dp_cache, indexD, m, subsetMap);

        igraph_vector_int_t vectorlist;
        IGRAPH_VECTOR_INT_INIT_FINALLY(&vectorlist, 0);

        igraph_vector_int_t edgelist;
        IGRAPH_VECTOR_INT_INIT_FINALLY(&edgelist, 0);

        IGRAPH_CHECK(igraph_get_shortest_path_dijkstra(graph, &vectorlist, &edgelist, m, k, weights, IGRAPH_ALL));


        for (int i = 0; i < igraph_vector_int_size(&edgelist); i++) {
            edgelist_all_set.insert(VECTOR(edgelist)[i]);
        }

        igraph_integer_t min_E_value = IGRAPH_INTEGER_MAX;
        int_set min_F;
        /*
         *  When the size of subset is > 2 we need to split the subset into E and F
         *  where E is singleton and F is subset of size (D.size() - 1)
         *  and repeat same process
         */
        if (D.size() > 2) {

            /*
                 We iterate through the subset D and split it into E and F where E is singleton set during every iteration.
                 This process leads to find out value where distance (E,k) + (F,k) is minimum.
            */
            igraph_real_t min_value = IGRAPH_INFINITY;
            igraph_integer_t D_size = D.size();
            for (igraph_integer_t i = 0; i < D_size; i++) {
                igraph_integer_t E_raw = *next(D.begin(), i);
                int_set F;
                int_set E;
                E.insert(E_raw);
                std::set_difference(D.begin(), D.end(), E.begin(), E.end(), std::inserter(F, F.end()));
                igraph_integer_t indexF = fetchIndexofMapofSets(D, subsetMap);
                igraph_real_t temp_value = MATRIX(*dp_cache, *E.begin(), k) + MATRIX(*dp_cache, indexF, k);
                if (temp_value < min_value) {
                    min_value = temp_value;
                    min_E_value = *E.begin();
                    min_F = F;
                }
            }

            igraph_vector_int_t vectorlist_1;
            IGRAPH_VECTOR_INT_INIT_FINALLY(&vectorlist_1, 0);

            igraph_vector_int_t edgelist_1;
            IGRAPH_VECTOR_INT_INIT_FINALLY(&edgelist_1, 0);

            IGRAPH_CHECK(igraph_get_shortest_path_dijkstra(graph, &vectorlist_1, &edgelist_1, k, min_E_value, weights, IGRAPH_ALL));

            for (int i = 0; i < igraph_vector_int_size(&edgelist_1); i++) {
                edgelist_all_set.insert(VECTOR(edgelist_1)[i]);
            }

            igraph_vector_int_destroy(&vectorlist_1);
            igraph_vector_int_destroy(&edgelist_1);
            IGRAPH_FINALLY_CLEAN(2);
        } else {
            /*
                If the size of subset is 2 then just shortest path from
                k to first element of subset and k to second element of subset is sufficient.
            */
            igraph_integer_t E1, F1;

            E1 = *D.begin();
            F1 = *next(D.begin(), 1);

            igraph_vector_int_t vectorlist_1;
            IGRAPH_VECTOR_INT_INIT_FINALLY(&vectorlist_1, 0);

            igraph_vector_int_t edgelist_1;
            IGRAPH_VECTOR_INT_INIT_FINALLY(&edgelist_1, 0);

            IGRAPH_CHECK(igraph_get_shortest_path_dijkstra(graph, &vectorlist_1, &edgelist_1, k, E1, weights, IGRAPH_ALL));

            igraph_vector_int_t vectorlist_2;
            IGRAPH_VECTOR_INT_INIT_FINALLY(&vectorlist_2, 0);

            igraph_vector_int_t edgelist_2;
            IGRAPH_VECTOR_INT_INIT_FINALLY(&edgelist_2, 0);

            IGRAPH_CHECK(igraph_get_shortest_path_dijkstra(graph, &vectorlist_2, &edgelist_2, k, F1, weights, IGRAPH_ALL));

            for (int i = 0; i < igraph_vector_int_size(&edgelist_1); i++) {
                edgelist_all_set.insert(VECTOR(edgelist_1)[i]);
            }

            for (int i = 0; i < igraph_vector_int_size(&edgelist_2); i++) {
                edgelist_all_set.insert(VECTOR(edgelist_2)[i]);
            }

            igraph_vector_int_destroy(&vectorlist_2);
            igraph_vector_int_destroy(&vectorlist_1);

            igraph_vector_int_destroy(&edgelist_1);
            igraph_vector_int_destroy(&edgelist_2);

            IGRAPH_FINALLY_CLEAN(4);

            min_F.insert(F1);
        }

        m = k;
        D = min_F;

        igraph_vector_int_destroy(&vectorlist);
        igraph_vector_int_destroy(&edgelist);

        IGRAPH_FINALLY_CLEAN(2);
    }

    for (auto i : edgelist_all_set) {
        igraph_vector_int_push_back(edgelist_all, i);
    }

    return IGRAPH_SUCCESS;
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
 * It is recommended that V &lt;= 50 and k &lt; 11.
 *
 * \sa \ref igraph_minimum_spanning_tree(), \ref igraph_spanner()
 */

igraph_error_t igraph_steiner_dreyfus_wagner(
    const igraph_t *graph, const igraph_vector_int_t *terminals, const igraph_vector_t *weights,
    igraph_real_t *res, igraph_vector_int_t *res_tree) {
    IGRAPH_HANDLE_EXCEPTIONS_BEGIN;

    const igraph_vector_t *pweights;
    igraph_vector_t iweights; // Will be used when weights are NULL.

    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_integer_t no_of_edges = igraph_ecount(graph);
    igraph_integer_t no_of_terminals = igraph_vector_int_size(terminals);

    if (!igraph_vector_int_isininterval(terminals, 0, no_of_nodes)) {
        IGRAPH_ERROR("Invalid vertex ID given as Steiner terminal.", IGRAPH_EINVVID);
    }
    if (!weights) {

        IGRAPH_CHECK(igraph_vector_init(&iweights, no_of_edges));
        igraph_vector_fill(&iweights, 1);
        pweights = &iweights;
        IGRAPH_FINALLY(igraph_vector_destroy, &iweights);
    } else {
        pweights = weights;
    }

    if (igraph_vector_size(pweights) != no_of_edges) {

        IGRAPH_ERROR("Invalid weight vector length.", IGRAPH_EINVAL);
    }
    if (no_of_edges > 0) {
        igraph_real_t minweight = igraph_vector_min(pweights);
        if (minweight < 0) {
            IGRAPH_ERRORF("Edge weights must be non-negative, got %g.", IGRAPH_EINVAL, minweight);
        } else if (minweight == 0) {
            /* TODO: can we support zero edge weights? */
            IGRAPH_ERROR("Weight vector contains zero weight.", IGRAPH_EINVAL);
        }
    }

    igraph_integer_t flag_terminals = 0;
    /* Handle the cases of the null graph and no terminals. */
    if (no_of_nodes == 0 || no_of_terminals == 0 || no_of_terminals == 1) {

        igraph_vector_int_clear(res_tree);
        *res = 0.0;
        flag_terminals = 1;
    } else if (no_of_terminals == 2) {
        igraph_vector_int_t vertices;
        igraph_real_t tree_weight = 0.0;
        igraph_vector_int_t edges_res;

        IGRAPH_VECTOR_INT_INIT_FINALLY(&vertices, 0);
        IGRAPH_VECTOR_INT_INIT_FINALLY(&edges_res, 0);

        IGRAPH_CHECK(igraph_get_shortest_path_dijkstra(graph, &vertices, &edges_res, VECTOR(*terminals)[0], VECTOR(*terminals)[1], pweights, IGRAPH_ALL));
        igraph_integer_t tree_size = igraph_vector_int_size(&edges_res);

        for (igraph_integer_t i = 0; i < tree_size; i++) {
            tree_weight += VECTOR(*pweights)[VECTOR(edges_res)[i]];
        }
        *res = tree_weight;

        IGRAPH_CHECK(igraph_vector_int_append(res_tree, &edges_res));

        igraph_vector_int_destroy(&vertices);
        igraph_vector_int_destroy(&edges_res);
        IGRAPH_FINALLY_CLEAN(2);

        flag_terminals = 1;
    }

    if (flag_terminals == 1) {
        if (!weights) {
            igraph_vector_destroy(&iweights);
            IGRAPH_FINALLY_CLEAN(1);
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
                    IGRAPH_ERROR("Not all Steiner terminals are in the same connected component.", IGRAPH_EINVAL);
                }
            }
        }

        igraph_vector_int_destroy(&membership);
        IGRAPH_FINALLY_CLEAN(1);
    }
    /* When every vertex is a Steiner terminal, the probably reduces to
     * finding a minimum spanning tree.
     */
    if (no_of_terminals == no_of_nodes) {
        IGRAPH_CHECK(igraph_minimum_spanning_tree(graph, res_tree, pweights));

        igraph_real_t tree_weight = 0.0;
        igraph_integer_t tree_size = igraph_vector_int_size(res_tree);
        for (igraph_integer_t i = 0; i < tree_size; i++) {
            tree_weight += VECTOR(*pweights)[VECTOR(*res_tree)[i]];
        }
        *res = tree_weight;
        return IGRAPH_SUCCESS;
    }

    igraph_vector_int_t steiner_terminals_copy;
    igraph_matrix_t dp_cache; // dynamic programming table
    igraph_integer_t q;
    std::set<int_set> allSubsets = std::set<int_set>();
    igraph_matrix_t distance;

    if (igraph_vector_size(pweights) != no_of_edges) {
        IGRAPH_ERRORF("Weight vector length does not match %" IGRAPH_PRId "vec size and %" IGRAPH_PRId "edges.", IGRAPH_FAILURE, igraph_vector_size(weights), no_of_edges);
    }
    IGRAPH_CHECK(igraph_matrix_init(&distance, no_of_nodes, no_of_nodes));
    IGRAPH_FINALLY(igraph_matrix_destroy, &distance);

    /*
     * Compute distances between all pairs of vertices. The Dreyfus - Wagner algorithm needs complete graph information
     * hence this step is necessary.
     */
    IGRAPH_CHECK(igraph_distances_dijkstra(graph, &distance, igraph_vss_all(), igraph_vss_all(), pweights, IGRAPH_ALL));

    IGRAPH_CHECK(igraph_vector_int_init_copy(&steiner_terminals_copy, terminals));
    IGRAPH_FINALLY(igraph_vector_int_destroy, &steiner_terminals_copy);
    igraph_vector_int_sort(&steiner_terminals_copy);

    /*
        Case - Tree when number of terminals are 3.
        Say, q,t2,t3 are terminals then Steiner Tree = min( distance(q,t1)+ distance(q,t2), distance(t1,t2) + distance(q,t1), distance(t1,t2) + distance(q,t2) )

    */
    if (no_of_terminals == 3) {

        igraph_integer_t steiner_terminals_copy_size = igraph_vector_int_size(&steiner_terminals_copy);
        igraph_real_t min_steiner_tree_dist = IGRAPH_INFINITY;
        for (igraph_integer_t i = 0; i < steiner_terminals_copy_size; i++) {

            igraph_integer_t ter_node = VECTOR(steiner_terminals_copy)[i];
            igraph_real_t weight_inter = 0.0;
            igraph_vector_int_t edges_inter;
            IGRAPH_VECTOR_INT_INIT_FINALLY(&edges_inter, 0);

            for (igraph_integer_t j = 0; j < steiner_terminals_copy_size; j++) {
                if (VECTOR(steiner_terminals_copy)[j] != ter_node) {
                    igraph_vector_int_t vertices_2;
                    igraph_real_t tree_weight = 0.0;
                    igraph_vector_int_t edges;

                    IGRAPH_VECTOR_INT_INIT_FINALLY(&vertices_2, 0);
                    IGRAPH_VECTOR_INT_INIT_FINALLY(&edges, 0);

                    IGRAPH_CHECK(igraph_get_shortest_path_dijkstra(graph, &vertices_2, &edges, ter_node, VECTOR(steiner_terminals_copy)[j], pweights, IGRAPH_ALL));
                    igraph_integer_t tree_size = igraph_vector_int_size(&edges);

                    for (igraph_integer_t k = 0; k < tree_size; k++) {
                        tree_weight += VECTOR(*pweights)[VECTOR(vertices_2)[k]];
                    }
                    weight_inter += tree_weight;
                    IGRAPH_CHECK(igraph_vector_int_append(&edges_inter, &edges));

                    igraph_vector_int_destroy(&vertices_2);
                    igraph_vector_int_destroy(&edges);
                    IGRAPH_FINALLY_CLEAN(2);
                }
            }
            if (weight_inter < min_steiner_tree_dist) {
                igraph_vector_int_update(res_tree, &edges_inter);
                *res = weight_inter;
                min_steiner_tree_dist = weight_inter;
            }

            igraph_vector_int_destroy(&edges_inter);

            IGRAPH_FINALLY_CLEAN(1);
        }

        igraph_matrix_destroy(&distance);
        igraph_vector_int_destroy(&steiner_terminals_copy);
        IGRAPH_FINALLY_CLEAN(2);
        flag_terminals = 1;
    }

    if (flag_terminals == 1) {

        if (!weights) {
            igraph_vector_destroy(&iweights);
            IGRAPH_FINALLY_CLEAN(1);
        }

        return IGRAPH_SUCCESS;
    }

    q = VECTOR(steiner_terminals_copy)[0];

    igraph_vector_int_remove(&steiner_terminals_copy, 0);

    /*
     * DP table with size number of vertices in the graph + 2 ^ (number of steiner_terminals_copy) - (number of steiner_terminals_copy + 1)
     * 2 ^ (number of steiner_terminals_copy) - (number of steiner_terminals_copy + 1) is number of subsets.
     * DP table with size number of vertices in the graph + 2 ^ (number of steiner_terminals_copy) - (number of steiner_terminals_copy + 1)
     * 2 ^ (number of steiner_terminals_copy) - (number of steiner_terminals_copy + 1) is number of subsets.
     */

    IGRAPH_CHECK(igraph_matrix_init(&dp_cache, no_of_nodes + pow(2, igraph_vector_int_size(&steiner_terminals_copy)) - (igraph_vector_int_size(&steiner_terminals_copy) + 1), no_of_nodes));
    IGRAPH_FINALLY(igraph_matrix_destroy, &dp_cache);

    igraph_matrix_fill(&dp_cache, IGRAPH_INFINITY);
    /*
     * for singleton value the distance in dp cahce is just the
     * distance between same vertices in distance matrix
     */
    for (igraph_integer_t i = 0; i < no_of_nodes; i++) {
        IGRAPH_ALLOW_INTERRUPTION();
        for (igraph_integer_t j = 0; j < no_of_nodes; j++) {
            MATRIX(dp_cache, i, j) = MATRIX(distance, i, j);
        }
    }

    /*
     * Usage: Data structure to store index value of subsets.
     * The index would be used in DP table.
     */
    dictionary subsetMap;
    IGRAPH_CHECK(generateSubsets(&steiner_terminals_copy, igraph_vector_int_size(&steiner_terminals_copy), no_of_nodes, subsetMap, allSubsets));

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

            for (igraph_integer_t j = 0; j < no_of_nodes; j++) {
                igraph_real_t distance1 = IGRAPH_INFINITY;

                int_set D_prime = D;
                D_prime.erase(*D.begin());

                std::set<int_set> Subsets = std::set<int_set>();;
                generateD_E(*D.begin(), D_prime, Subsets);

                for (auto E : Subsets) {

                    igraph_integer_t indexOfSubsetE;
                    if (E.size() == 1) {
                        indexOfSubsetE = *E.begin();
                    } else {
                        indexOfSubsetE = fetchIndexofMapofSets(E, subsetMap);
                    }

                    igraph_real_t distanceEJ = MATRIX(dp_cache, indexOfSubsetE, j);


                    int_set DMinusE = D;

                    for (auto elem : E) {
                        DMinusE.erase(elem);
                    }


                    igraph_integer_t indexOfSubsetDMinusE;
                    if (DMinusE.size() == 1) {
                        indexOfSubsetDMinusE = *DMinusE.begin();
                    } else {
                        indexOfSubsetDMinusE = fetchIndexofMapofSets(DMinusE, subsetMap);
                    }

                    distance1 = std::min(distance1, distanceEJ + MATRIX(dp_cache, indexOfSubsetDMinusE, j));
                }

                for (igraph_integer_t k = 0; k < no_of_nodes; k++) {
                    MATRIX(dp_cache, indexOfSubsetD, k) = std::min(MATRIX(dp_cache, indexOfSubsetD, k), MATRIX(distance, k, j) + distance1);
                }
            }
        }
    }


    std::set<int_set> E_subsets = std::set<int_set>();
    int_set C_prime;
    int_set C;
    igraph_integer_t C_1 = igraph_vector_int_get(&steiner_terminals_copy,0);
    for (int i = 0 ; i < igraph_vector_int_size(&steiner_terminals_copy) ; i++){
        C.insert(igraph_vector_int_get(&steiner_terminals_copy,i));
        if (i!=0)   
            C_prime.insert(igraph_vector_int_get(&steiner_terminals_copy,i));
    }
    generateD_E(C_1, C_prime, E_subsets);

    igraph_real_t distance2 = IGRAPH_INFINITY;

    for (igraph_integer_t j = 0; j < no_of_nodes; j++) {
        igraph_real_t distance1 = IGRAPH_INFINITY;
        IGRAPH_ALLOW_INTERRUPTION();
        for (auto E : E_subsets) {
            igraph_integer_t indexE = (E.size() == 1) ? *E.begin() : fetchIndexofMapofSets(E,subsetMap);
            int_set CMinusE;
            std::set_difference(C.begin(),C.end(),E.begin(),E.end(),std::inserter(CMinusE,CMinusE.end()));
            igraph_integer_t indexC_E = (CMinusE.size() == 1) ? *CMinusE.begin() : fetchIndexofMapofSets(CMinusE,subsetMap);
            igraph_real_t distanceFJ = MATRIX(dp_cache,indexE,j);



            if (((distanceFJ + (MATRIX(dp_cache, indexC_E, j))) < distance1)) {
                distance1 = distanceFJ + (MATRIX(dp_cache, indexC_E, j));
            }
        }

        if ((MATRIX(distance, q, j) + distance1) < distance2) {
            distance2 = MATRIX(distance, q, j) + distance1;
        }
    }
    *res = distance2;

    int_set newSet;

    for (igraph_integer_t i = 0; i < steiner_terminals_copy_size; i++) {
        newSet.insert(VECTOR(steiner_terminals_copy)[i]);
    }
    igraph_integer_t indexD = fetchIndexofMapofSets(newSet, subsetMap);


    igraph_vector_int_clear(res_tree);
    IGRAPH_CHECK(generate_steiner_tree_exact(graph, pweights, &dp_cache, indexD, q, res_tree, subsetMap));
    igraph_matrix_destroy(&distance);
    igraph_vector_int_destroy(&steiner_terminals_copy);
    igraph_matrix_destroy(&dp_cache);

    IGRAPH_FINALLY_CLEAN(3);

    if (!weights) {
        igraph_vector_destroy(&iweights);
        IGRAPH_FINALLY_CLEAN(1);
        return IGRAPH_SUCCESS;
    }

    IGRAPH_HANDLE_EXCEPTIONS_END;
    return IGRAPH_SUCCESS;
}
