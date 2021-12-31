
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

    if (Size % 2 == 1 || Size > 32640) {
        return 0;
    }

    check_err(igraph_vector_init(&edges, Size));
    for (size_t i=0; i < Size; ++i) {
        VECTOR(edges)[i] = Data[i];
    }

    if (! igraph_create(&graph, &edges, 0, IGRAPH_UNDIRECTED)) {
        {
            igraph_vector_ptr_t separators;
            check_err(igraph_vector_ptr_init(&separators, 0));
            check_err(igraph_all_minimal_st_separators(&graph, &separators));
            IGRAPH_VECTOR_PTR_SET_ITEM_DESTRUCTOR(&separators, igraph_vector_destroy);
            igraph_vector_ptr_destroy_all(&separators);
        }

        {
            igraph_vector_ptr_t separators;
            check_err(igraph_vector_ptr_init(&separators, 0));
            check_err(igraph_minimum_size_separators(&graph, &separators));
            IGRAPH_VECTOR_PTR_SET_ITEM_DESTRUCTOR(&separators, igraph_vector_destroy);
            igraph_vector_ptr_destroy_all(&separators);
        }

        igraph_destroy(&graph);
    }

    igraph_vector_destroy(&edges);

    return 0;  // Non-zero return values are reserved for future use.
}
