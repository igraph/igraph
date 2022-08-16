#include <igraph.h>

#include "test_utilities.h"

/*
 * This test tries to ensure that property caching works correctly in the most
 * "dangerous" scenarios (goijg from disconnected to connected graphs and
 * back). Feel free to extend the file later on with regression tests for any
 * issues that we might find related to caching.
 *
 * There are no direct APIs to reach into the cache internals (because the
 * presence of the cache is an implementation detail). Therefore, we simply do
 * the following:
 *
 * - perform some operation
 * - query the value of certain properties that we know that are typically
 *   cached (e.g., connectedness, presence of multi- and loop edges) etc
 * - invalidate the cache
 * - perform the query again, check whether the results match
 *
 * The assumption is that functions like igraph_is_simple(), igraph_has_loops()
 * etc work correctly if the cache is empty; correctness of these functions
 * without a cache should be tested in other unit tests, not here.
 */

igraph_error_t has_mutual_nonloop_edge(const igraph_t* graph, igraph_bool_t* result) {
    return igraph_has_mutual(graph, result, /* loops = */ 0);
}

igraph_error_t has_mutual_edge(const igraph_t* graph, igraph_bool_t* result) {
    return igraph_has_mutual(graph, result, /* loops = */ 1);
}

igraph_error_t is_weakly_connected(const igraph_t* graph, igraph_bool_t* result) {
    return igraph_is_connected(graph, result, IGRAPH_WEAK);
}

igraph_error_t is_strongly_connected(const igraph_t* graph, igraph_bool_t* result) {
    return igraph_is_connected(graph, result, IGRAPH_STRONG);
}

igraph_error_t is_forest(const igraph_t* graph, igraph_bool_t* result) {
    return igraph_is_forest(graph, result, /* roots = */ 0, IGRAPH_ALL);
}

void validate_properties(const igraph_t* graph) {
    igraph_bool_t result, cached_result, recalculated_result;

#define CHECK(func) { \
    igraph_invalidate_cache(graph); \
    func(graph, &result); \
    func(graph, &cached_result); \
    igraph_invalidate_cache(graph); \
    func(graph, &recalculated_result); \
    IGRAPH_ASSERT(result == cached_result); \
    IGRAPH_ASSERT(recalculated_result == cached_result); \
}

    CHECK(igraph_is_simple);
    CHECK(igraph_has_loop);
    CHECK(igraph_has_multiple);
    CHECK(igraph_is_dag);
    CHECK(is_forest);
    CHECK(is_weakly_connected);
    CHECK(is_strongly_connected);
    CHECK(has_mutual_edge);
    CHECK(has_mutual_nonloop_edge);
}

void test_basic_operations(igraph_t* graph) {
    validate_properties(graph);
    igraph_add_vertices(graph, 1, /* attr = */ NULL);
    validate_properties(graph);
    igraph_add_vertices(graph, 2, /* attr = */ NULL);
    validate_properties(graph);
    igraph_add_edge(graph, 0, 1);
    validate_properties(graph);
    igraph_add_edge(graph, 0, 2);
    validate_properties(graph);
    igraph_add_edge(graph, 1, 2);
    validate_properties(graph);
    igraph_add_edge(graph, 2, 0);
    validate_properties(graph);
    igraph_delete_edges(graph, igraph_ess_all(IGRAPH_EDGEORDER_ID));
    validate_properties(graph);
    igraph_add_edge(graph, 0, 2);
    validate_properties(graph);
    igraph_delete_vertices(graph, igraph_vss_1(1));
    validate_properties(graph);
    igraph_delete_vertices(graph, igraph_vss_all());
    validate_properties(graph);
}

int test_basic_operations_directed() {
    igraph_t g;

    igraph_empty(&g, 0, IGRAPH_DIRECTED);
    test_basic_operations(&g);
    igraph_destroy(&g);

    return 0;
}

int test_basic_operations_undirected() {
    igraph_t g;

    igraph_empty(&g, 0, IGRAPH_UNDIRECTED);
    test_basic_operations(&g);
    igraph_destroy(&g);

    return 0;
}

int main() {

    RUN_TEST(test_basic_operations_directed);
    RUN_TEST(test_basic_operations_undirected);

    return 0;
}
