
#include <igraph.h>
#include <cstdlib>

inline void check_err(int err) {
    if (err)
        abort();
}

extern "C" int LLVMFuzzerTestOneInput(const uint8_t *Data, size_t Size) {
    igraph_t graph;
    igraph_vector_t edges;

    igraph_set_error_handler(&igraph_error_handler_ignore);

    if (Size % 2 == 1 || Size > 65280) {
        return 0;
    }

    check_err(igraph_vector_init(&edges, Size));
    for (size_t i=0; i < Size; ++i) {
        VECTOR(edges)[i] = Data[i];
    }

    if (! igraph_create(&graph, &edges, 0, IGRAPH_DIRECTED)) {
        igraph_integer_t conn;

        check_err(igraph_vertex_connectivity(&graph, &conn, 0));
        check_err(igraph_to_undirected(&graph, IGRAPH_TO_UNDIRECTED_COLLAPSE, nullptr));
        check_err(igraph_vertex_connectivity(&graph, &conn, 0));

        igraph_destroy(&graph);
    }

    igraph_vector_destroy(&edges);

    return 0;  // Non-zero return values are reserved for future use.
}
