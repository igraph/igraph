
#include "igraph.h"
#include "core/exceptions.h"

#include <cstring>
#include <cmath>
#include <map>
#include <climits>
#include <vector>
#include <set>
#include <algorithm>
#include <iterator>


/*
 *  Generates Subsets of length > 1 for list of integer values.
 *  As algorithm computes the minimum distance between set of variables ranging from
 *  length 2 to steiner_terminals - 1, we need these subsets.
 *   for e.g [2,3,4] --> [ [2,3], [2,4], [3,4], [2,3,4] ]
 *  Populates the subsetMap data structure as well.
 *
 */

std::set<std::set<igraph_integer_t>> generateSubsets(igraph_vector_int_t steinerTerminals, igraph_integer_t n, igraph_integer_t graphsize, std::map<std::set<igraph_integer_t>, igraph_integer_t>* subsetMap) {
    igraph_integer_t count = ((igraph_integer_t)1 << n);
    std::set<std::set<igraph_integer_t>> allSubsets;
    igraph_integer_t subsetIndex = graphsize;

    // The outer for loop will run 2^n times to get all subset .
    // Here variable i will act as a binary counter

    for (igraph_integer_t i = 0; i < count; i++) {
        // The inner for loop will run n times , As the maximum number of elements a set can have is n
        // This loop will generate a subset
        std::set<igraph_integer_t> newSubset;
        for (igraph_integer_t j = 0; j < n; j++) {
            // This if condition will check if jth bit in binary representation of  i  is set or not
            // if the value of (i & (1 << j)) is greater than 0 , include arr[j] in the current subset
            // otherwise exclude arr[j]
            if ((i & ((igraph_integer_t)1 << j)) > 0) {
                newSubset.insert(VECTOR(steinerTerminals)[j]);
            }
        }

        if (newSubset.size() > 1) {
            if (allSubsets.find(newSubset) == allSubsets.end()) {
                allSubsets.insert(newSubset);
                (*subsetMap).insert(std::make_pair(newSubset, subsetIndex));
                subsetIndex++;
            }
        }
    }
    return allSubsets;
}

/*
 *
 * Purpose: Fetching Index of a subset from subsetMap in order to store and look-up
 * the value of subset from DP table.
 */

igraph_integer_t fetchIndexofMapofSets(std::set<igraph_integer_t> subset, std::map<std::set<igraph_integer_t>, igraph_integer_t>* subsetMap) {
    std::map<std::set<igraph_integer_t>, igraph_integer_t>::iterator it;
    for (it = subsetMap->begin(); it != subsetMap->end(); ++it) {
        if (it->first == subset) {
            return it->second;
        }
    }
    IGRAPH_FATAL("The Subset's index that you tried to find doesn't exist. Hence the code won't run.");
}

/*
 * Purpose: Retriving the value of subset from given index.
 *
 */

std::set<igraph_integer_t> fetchSetsBasedonIndex(igraph_integer_t index, std::map<std::set<igraph_integer_t>, igraph_integer_t>* subsetMap) {
    std::map<std::set<igraph_integer_t>, igraph_integer_t>::iterator it;
    for (it = subsetMap->begin(); it != subsetMap->end(); ++it) {
        if (it->second == index) {
            return it->first;
        }
    }
    IGRAPH_FATAL("The index that you tried to find doesn't exist. Hence the code won't run.");
}

/*
 *
 * Calculating factorial of a number.
 *
 */

igraph_integer_t factorial(igraph_integer_t n) {
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
 *
 */

igraph_integer_t Combination(igraph_integer_t n, igraph_integer_t r) {
    return factorial(n) / (factorial(n - r) * factorial(r));
}

/* \function  findMinimumK
 *
 *
 *  Iterate through a particular row of the DP table
 *  find the value of column where distance (q,k) + distance(Sk(D))
 *  where D is a subset. This signifies a vertex that is part
 *  of minimum path from vertex q to D
 *  \param dp_cache The DP table.
 *
 *  \param indexD Index of the subset D.
 *
 *  \param q vertex from who minimum path needs tp be calculated
 */

igraph_integer_t findMinimumK(igraph_matrix_t *dp_cache, igraph_integer_t indexD, igraph_integer_t q) {

    igraph_integer_t min_col_num = -1;
    igraph_integer_t min_sum_for_col;

    for (igraph_integer_t i = 0; i < dp_cache->ncol; i++) {
        if (q != i) {

            if (min_col_num == -1) {
                min_col_num = i;
                min_sum_for_col = MATRIX(*dp_cache, q, i) + MATRIX(*dp_cache, indexD, i);
            } else if (MATRIX(*dp_cache, q, i) + MATRIX(*dp_cache, indexD, i) < min_sum_for_col) {
                min_col_num = i;
                min_sum_for_col = MATRIX(*dp_cache, q, i) + MATRIX(*dp_cache, indexD, i);
            }
        }
    }

    return min_col_num;
}

/* \function generate_steiner_tree_appx()
 *
 *   Generation of the Steiner Tree based on calculation of minimum distances in DP table.
 *
 *  \param graph The graph object.
 *  \param weights The edge weights. All edge weights must be
 *                 non-negative. Additionally, no edge weight may be NaN.
 *
 *  \param dp_cache The DP table.
 *
 *  \param indexD The index of subset D in DP table.
 *
 *  \param q The vertex that was remmoved from steiner terminals.
 *
 * \param mode How to determine the local neighborhood of each vertex
 *  in directed graphs. Ignored in undirected graphs.
 * \clist
 *         \cli IGRAPH_ALL
 *              take both in- and out-neighbours;
 *              this is a reasonable default for high-level interfaces.
 *         \cli IGRAPH_OUT
 *              take only out-neighbours
 *         \cli IGRAPH_IN
 *              take only in-neighbours
 *         \endclist
 *
 *  \param vectorlist_all The vector to capture vertices in resultant Steiner Tree.
 *
 *  \param edgelist_all The vector to capture edges in resultant Steiner Tree.
 *
 */

igraph_error_t generate_steiner_tree_exact(const igraph_t *graph, const igraph_vector_t *weights,
        igraph_matrix_t *dp_cache, igraph_integer_t indexD, igraph_integer_t q, igraph_vector_int_t *vectorlist_all, igraph_vector_int_t *edgelist_all, std::map<std::set<igraph_integer_t>, igraph_integer_t>* subsetMap) {

    std::set<igraph_integer_t> C = fetchSetsBasedonIndex(indexD, subsetMap);

    // Initially the value of m is the vertex that was removed from Steiner Terminals
    igraph_integer_t m = q;

    igraph_integer_t len = C.size();
    std::set<igraph_integer_t> D = C;

    while (D.size() > 1) {

        indexD = fetchIndexofMapofSets(D, subsetMap);

        // Finding the bridge vertex from m to subset D which is part of shortest path.
        igraph_integer_t k = findMinimumK(dp_cache, indexD, m);

        igraph_vector_int_t vectorlist;
        IGRAPH_CHECK(igraph_vector_int_init(&vectorlist, 1));
        IGRAPH_FINALLY(igraph_vector_int_destroy, &vectorlist);

        igraph_vector_int_t edgelist;
        IGRAPH_CHECK(igraph_vector_int_init(&edgelist, 1));
        IGRAPH_FINALLY(igraph_vector_int_destroy, &edgelist);

        igraph_get_shortest_path_dijkstra(graph, &vectorlist, &edgelist, m, k, weights, IGRAPH_ALL);

        igraph_vector_int_append(vectorlist_all, &vectorlist);
        igraph_vector_int_append(edgelist_all, &edgelist);

        igraph_integer_t min_E_value = IGRAPH_INTEGER_MAX;
        std::set<igraph_integer_t> min_F;
        /*
         *  When the size of subset is > 2 we need to split the subset into E and F
         *  where E is singleton and F is subset of size (D.size() - 1)
         *  and repeat same process
         */
        if (D.size() > 2) {

            /*
             *  Purpose is to fetch the count of number of times scan in the DP table is necessary.
             *  The formula |C| combination (D.size - 1) where combination is mathematical combination
             *  and C is steiner terminal set with one element removed.
             */
            igraph_integer_t numElementsScan = Combination(len, D.size() - 1);
            igraph_real_t min_value = IGRAPH_INFINITY;

            igraph_integer_t holder = fetchIndexofMapofSets(D, subsetMap);
            for (igraph_integer_t i = 1; i <= numElementsScan; i++) {

                // retrieving the set associated with index = holder - i from subsetMap

                std::set<igraph_integer_t> F = fetchSetsBasedonIndex(holder - i, subsetMap);

                std::set<igraph_integer_t> E;

                std::set_difference(D.begin(), D.end(), F.begin(), F.end(), std::inserter(E, E.end()));

                igraph_real_t temp_value = MATRIX(*dp_cache, *E.begin(), k) + MATRIX(*dp_cache, holder - i, k);

                if (temp_value < min_value) {
                    min_value = temp_value;
                    min_E_value = *E.begin();
                    min_F = F;
                }
            }

            igraph_vector_int_t vectorlist_1;
            IGRAPH_CHECK(igraph_vector_int_init(&vectorlist_1, 1));
            IGRAPH_FINALLY(igraph_vector_int_destroy, &vectorlist_1);

            igraph_vector_int_t edgelist_1;
            IGRAPH_CHECK(igraph_vector_int_init(&edgelist_1, 1));
            IGRAPH_FINALLY(igraph_vector_int_destroy, &edgelist_1);

            igraph_get_shortest_path_dijkstra(graph, &vectorlist_1, &edgelist_1, k, min_E_value, weights, IGRAPH_ALL);

            igraph_vector_int_append(vectorlist_all, &vectorlist_1);
            igraph_vector_int_append(edgelist_all, &edgelist_1);

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
            ;

            igraph_vector_int_t vectorlist_1;
            IGRAPH_CHECK(igraph_vector_int_init(&vectorlist_1, 1));
            IGRAPH_FINALLY(igraph_vector_int_destroy, &vectorlist_1);

            igraph_vector_int_t edgelist_1;
            IGRAPH_CHECK(igraph_vector_int_init(&edgelist_1, 1));
            IGRAPH_FINALLY(igraph_vector_int_destroy, &edgelist_1);

            igraph_get_shortest_path_dijkstra(graph, &vectorlist_1, &edgelist_1, k, E1, weights, IGRAPH_ALL);

            igraph_vector_int_t vectorlist_2;
            IGRAPH_CHECK(igraph_vector_int_init(&vectorlist_2, 1));
            IGRAPH_FINALLY(igraph_vector_int_destroy, &vectorlist_2);

            igraph_vector_int_t edgelist_2;
            IGRAPH_CHECK(igraph_vector_int_init(&edgelist_2, 1));
            IGRAPH_FINALLY(igraph_vector_int_destroy, &edgelist_2);

            igraph_get_shortest_path_dijkstra(graph, &vectorlist_2, &edgelist_2, k, F1, weights, IGRAPH_ALL);

            igraph_vector_int_append(vectorlist_all, &vectorlist_1);
            igraph_vector_int_append(vectorlist_all, &vectorlist_2);

            igraph_vector_int_append(edgelist_all, &edgelist_1);
            igraph_vector_int_append(edgelist_all, &edgelist_2);

            igraph_vector_int_destroy(&vectorlist_2);
            igraph_vector_int_destroy(&vectorlist_1);

            igraph_vector_int_destroy(&edgelist_1);
            igraph_vector_int_destroy(&edgelist_2);

            IGRAPH_FINALLY_CLEAN(4);

            std::set<igraph_integer_t> min_F;
            min_F.insert(F1);
        }

        m = k;
        D = min_F;

        igraph_vector_int_destroy(&vectorlist);
        igraph_vector_int_destroy(&edgelist);
        IGRAPH_FINALLY_CLEAN(2);
    }

    return IGRAPH_SUCCESS;
}

/* \function  igraph_steiner_dreyfus_wagner
 *  \brief This function calculates Steiner Tree on undirected graph using Dreyfus-Wagner algorithm
 *
 *
 *  </para><para>
 *   Steiner Tree is a tree with minimum length connecting steiner terminals which are subset of vertices
 *   in a graph. If Steiner terminals equal to number of vertices in the graph then problem is reduced
 *   to minimum spanning tree.
 *  The Steiner Tree problem is NP-Hard but since it's FPT, Dreyfus-Wagner algorithm is able to calculate
 *  Exact value of Steiner Tree on smaller graphs.
 *
 *  The algorithm uses Dynamic Programming approach and it assumes the graph is undirected.
 *
 * </para><para>
 * Reference:
 * S.E Dreyfus, R.A Wagner
 * The steiner problem in graphs
 * Networks Journal, 1971
 * https://doi.org/10.1002/net.3230010302
 *
 * \param graph The graph object.
 *
 * \param res Pointer to a real number, this will contain the result which is distance of Steiner Tree.
 *
 * \param res_tree Pointer to a vector, this will contain the Edges that are part of Steiner Tree
 *
 * \param steiner_terminals Pointer to a vector, this will contain vertices
 *
 * \param weights The edge weights. All edge weights must be
 *       non-negative. Additionally, no edge weight may be NaN.
 *
 * \return Error code:IGRAPH_ERROR
 *
 * Time complexity: O( 3^k ∗ V + 2^k ∗V^2 + V∗(V+E) ∗ log(V) )
 * where V is vertices and E is edges in the graph and k is set of steiner terminals.
 * It's recommended that V <= 50 and k < 11
 *
 * \sa \ref igraph_distances_johnson(), generateSubsets(), fetchIndexofMapofSets(), generate_steiner_tree_appx()
 */

igraph_error_t igraph_steiner_dreyfus_wagner(const igraph_t *graph, const igraph_vector_int_t *steiner_terminals,
        const igraph_vector_t *weights, igraph_real_t *res, igraph_vector_int_t *res_tree) {
    IGRAPH_HANDLE_EXCEPTIONS_BEGIN


    igraph_integer_t no_of_vertices = igraph_vcount(graph);
    igraph_integer_t no_of_edges = igraph_ecount(graph);
    if (no_of_vertices == 0 || (no_of_vertices == 1)) { // graph is empty
        IGRAPH_CHECK(igraph_minimum_spanning_tree(graph, res_tree, weights));
        *res = 0.0;
        return IGRAPH_SUCCESS;
    }
    igraph_bool_t IsConnected;
    igraph_is_connected(graph, &IsConnected, IGRAPH_WEAK);
    if (!IsConnected) {

        IGRAPH_ERROR("The graph is disconnected.", IGRAPH_EINVAL);
    }




    /*
     *  If the steiner terminals is number of vertices in graph then problem
     *  is reduced to minimum spanning tree which is tractable.
     */
    igraph_real_t min = igraph_vector_min(weights);
    if (min < 0) {
        IGRAPH_ERRORF("Weight vector must be non-negative, got %g.", IGRAPH_EINVAL, min);
    } else if (min == 0) {
        IGRAPH_ERROR("Weight vector contains zero weight.", IGRAPH_EINVAL);
    }
    if (igraph_vector_int_size(steiner_terminals) == no_of_vertices) {

        IGRAPH_CHECK(igraph_minimum_spanning_tree(graph, res_tree, weights));
        igraph_real_t size_value = 0.0;
        for (igraph_integer_t i = 0; i < igraph_vector_int_size(res_tree); i++) {
            size_value  += VECTOR(*weights)[VECTOR(*res_tree)[i]];
        }
        *res = size_value;
        return IGRAPH_SUCCESS;
    }


    igraph_vector_int_t steiner_terminals_copy;
    igraph_matrix_t dp_cache; // dynamic programming table
    igraph_integer_t q;
    std::set<std::set<igraph_integer_t>> allSubsets;
    igraph_matrix_t distance;

    if (igraph_vector_size(weights) != no_of_edges) {
        IGRAPH_ERRORF("Weight vector length does not match %" IGRAPH_PRId "vec size and %" IGRAPH_PRId "edges.", IGRAPH_FAILURE, igraph_vector_size(weights), no_of_edges);
    }
    IGRAPH_CHECK(igraph_matrix_init(&distance, no_of_vertices, no_of_vertices));
    IGRAPH_FINALLY(igraph_matrix_destroy, &distance);

    /*
     *  Johnson's algorithm calculates all pairs shortest path
     *  and returns distance matrix. The Dreyfus - Wagner algorithm needs complete graph
     *   hence this step is necessary.
     */
    IGRAPH_CHECK(igraph_distances_johnson(graph, &distance, igraph_vss_all(), igraph_vss_all(), weights));

    /*
        Setting distance from vertex to itself as 0.
    */
    for (igraph_integer_t i = 0; i < no_of_vertices; i++) {
        if (MATRIX(distance, i, i) != 0) {
            MATRIX(distance, i, i) = 0;
        }
    }
    IGRAPH_CHECK(igraph_vector_int_init_copy(&steiner_terminals_copy, steiner_terminals));
    IGRAPH_FINALLY(igraph_vector_int_destroy, &steiner_terminals_copy);
    igraph_vector_int_sort(&steiner_terminals_copy);
    q = VECTOR(steiner_terminals_copy)[0];

    igraph_vector_int_remove(&steiner_terminals_copy, 0);

    /*
     *  DP table with size number of vertices in the graph + 2 ^ (number of steiner_terminals_copy - 1)
     *  2 ^ (number of steiner_terminals_copy - 1) is number of subsets.
     */
    IGRAPH_CHECK(igraph_matrix_init(&dp_cache, no_of_vertices + pow(2, igraph_vector_int_size(&steiner_terminals_copy) - 1), no_of_vertices));
    IGRAPH_FINALLY(igraph_matrix_destroy, &dp_cache);

    igraph_matrix_fill(&dp_cache, IGRAPH_INFINITY);
    /*
     *  for singleton value the distance in dp cahce is just the
     *  distance between same vertices in distance matrix
     */
    for (igraph_integer_t i = 0; i < no_of_vertices; i++) {
        for (igraph_integer_t j = 0; j < no_of_vertices; j++) {
            MATRIX(dp_cache, i, j) = MATRIX(distance, i, j);
        }
    }

    /*
    *
    * Usage: Data structure to store index value of subsets.
    * The index would be used in DP table
    *
    */
    std::map<std::set<igraph_integer_t>, igraph_integer_t> subsetMap;
    allSubsets = generateSubsets(steiner_terminals_copy, igraph_vector_int_size(&steiner_terminals_copy), no_of_vertices, &subsetMap);

    for (igraph_integer_t m = 2; m <= igraph_vector_int_size(&steiner_terminals_copy); m++) {
        for (igraph_integer_t i = 0; i < (igraph_integer_t)allSubsets.size(); i++) {
            auto it = allSubsets.begin();
            std::advance(it, i);
            std::set<igraph_integer_t> D = *it;
            igraph_integer_t indexOfSubsetD;
            indexOfSubsetD = fetchIndexofMapofSets(D, &subsetMap);

            for (igraph_integer_t j = 0; j < no_of_vertices; j++) {
                MATRIX(dp_cache, indexOfSubsetD, j) = IGRAPH_INFINITY;
            }

            for (igraph_integer_t j = 0; j < no_of_vertices; j++) {
                igraph_real_t distance1 = IGRAPH_INFINITY;
                std::set<igraph_integer_t>::iterator subset_D_iterator;

                for (subset_D_iterator = D.begin(); subset_D_iterator != D.end(); subset_D_iterator++) {
                    igraph_integer_t E = *subset_D_iterator;
                    if (E != j) {
                        igraph_integer_t distanceEJ = MATRIX(distance, E, j);

                        std::set<igraph_integer_t> DMinusE = D;

                        /*
                         *  A set with Singleton value E removed from subset D
                         */
                        for (std::set<igraph_integer_t>::iterator iter = DMinusE.begin(); iter != DMinusE.end();) {
                            if (*iter == E) {
                                iter = DMinusE.erase(iter);
                                break;
                            }
                            ++iter;
                        }

                        igraph_integer_t indexOfSubsetDMinusE;
                        if (DMinusE.size() == 1) {
                            std::set<igraph_integer_t>::iterator node = DMinusE.begin();
                            indexOfSubsetDMinusE = *node;
                        } else {
                            indexOfSubsetDMinusE = fetchIndexofMapofSets(DMinusE, &subsetMap);
                        }

                        if ((distanceEJ + MATRIX(dp_cache, indexOfSubsetDMinusE, j)) < distance1) {
                            distance1 = distanceEJ + (MATRIX(dp_cache, indexOfSubsetDMinusE, j));
                        }
                    }
                }

                for (igraph_integer_t k = 0; k < no_of_vertices; k++) {
                    MATRIX(dp_cache, indexOfSubsetD, k) = std::min(MATRIX(dp_cache, indexOfSubsetD, k), MATRIX(distance, k, j) + distance1);
                }
            }
        }
    }

    igraph_real_t distance2 = IGRAPH_INFINITY;

    for (igraph_integer_t j = 0; j < no_of_vertices; j++) {
        igraph_real_t distance1 = IGRAPH_INFINITY;
        for (igraph_integer_t subset_C_iterator = 0; subset_C_iterator < igraph_vector_int_size(steiner_terminals); subset_C_iterator++) {
            igraph_integer_t F = VECTOR(steiner_terminals_copy)[subset_C_iterator];
            igraph_integer_t distanceFJ = MATRIX(distance, F, j);

            std::set<igraph_integer_t> CMinusF;

            for (igraph_integer_t k = 0; k < igraph_vector_int_size(steiner_terminals); k++) {

                if (VECTOR(steiner_terminals_copy)[k] != F) {
                    CMinusF.insert(VECTOR(steiner_terminals_copy)[k]);
                }
            }

            igraph_integer_t indexOfSubsetCMinusF = fetchIndexofMapofSets(CMinusF, &subsetMap);

            if (distanceFJ != 0 && (distanceFJ + (MATRIX(dp_cache, indexOfSubsetCMinusF, j)) < distance1)) {
                distance1 = distanceFJ + (MATRIX(dp_cache, indexOfSubsetCMinusF, j));
            }
        }

        if (q != j && MATRIX(distance, q, j) + distance1 < distance2) {
            distance2 = MATRIX(distance, q, j) + distance1;
        }
    }
    *res = distance2;

    std::set<igraph_integer_t> newSet;
    for (igraph_integer_t i = 0; i < igraph_vector_int_size(&steiner_terminals_copy); i++) {
        newSet.insert(VECTOR(steiner_terminals_copy)[i]);
    }
    igraph_integer_t indexD = fetchIndexofMapofSets(newSet, &subsetMap);

    igraph_vector_int_t vectorlist_all;

    IGRAPH_CHECK(igraph_vector_int_init(&vectorlist_all, 1));

    IGRAPH_CHECK(generate_steiner_tree_exact(graph, weights, &dp_cache, indexD, q, &vectorlist_all, res_tree, &subsetMap));

    // Default vector initiation add 0 and hence it's removed.
    igraph_vector_int_remove(res_tree, 0);

    igraph_matrix_destroy(&distance);
    igraph_vector_int_destroy(&steiner_terminals_copy);

    igraph_matrix_destroy(&dp_cache);
    igraph_vector_int_destroy(&vectorlist_all);

    IGRAPH_FINALLY_CLEAN(3);
    IGRAPH_HANDLE_EXCEPTIONS_END
    return IGRAPH_SUCCESS;
}
