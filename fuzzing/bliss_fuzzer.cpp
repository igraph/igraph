
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

    if (Size % 2 == 1 || Size > 512) {
        return 0;
    }

    check_err(igraph_vector_init(&edges, Size));
    for (size_t i=0; i < Size; ++i) {
        VECTOR(edges)[i] = Data[i];
    }

    /* Undirected */
    if (! igraph_create(&graph, &edges, 0, IGRAPH_UNDIRECTED)) {
        igraph_bool_t multi;

        check_err(igraph_has_multiple(&graph, &multi));

        /* Bliss does not support multigraphs and the input is currently not checked */
        if (! multi) {
            igraph_bliss_info_t info;
            igraph_vector_ptr_t generators;
            check_err(igraph_vector_ptr_init(&generators, 0));
            check_err(igraph_automorphism_group(&graph, nullptr, &generators, IGRAPH_BLISS_FS, &info));
            IGRAPH_VECTOR_PTR_SET_ITEM_DESTRUCTOR(&generators, igraph_vector_destroy);
            igraph_vector_ptr_destroy_all(&generators);
        }

        igraph_destroy(&graph);
    }

    /* Directed */
    if (! igraph_create(&graph, &edges, 0, IGRAPH_DIRECTED)) {
        igraph_bool_t multi;

        check_err(igraph_has_multiple(&graph, &multi));

        /* Bliss does not support multigraphs and the input is currently not checked */
        if (! multi) {
            igraph_bliss_info_t info;
            igraph_vector_ptr_t generators;
            check_err(igraph_vector_ptr_init(&generators, 0));
            check_err(igraph_automorphism_group(&graph, nullptr, &generators, IGRAPH_BLISS_FS, &info));
            IGRAPH_VECTOR_PTR_SET_ITEM_DESTRUCTOR(&generators, igraph_vector_destroy);
            igraph_vector_ptr_destroy_all(&generators);
        }

        igraph_destroy(&graph);
    }

    igraph_vector_destroy(&edges);

    return 0;  // Non-zero return values are reserved for future use.
}
