#include <igraph.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

// Optimized Havel-Hakimi Algorithm
int compare(const void *a, const void *b) {
    return (*(int*)b - *(int*)a);
}

bool havel_hakimi(int *degrees, int size) {
    while (1) {
        qsort(degrees, size, sizeof(int), compare);

        if (degrees[0] == 0) {
            return true;  // All degrees are zero
        }

        int max_deg = degrees[0];
        if (max_deg >= size) {
            return false;  // Invalid degree sequence
        }

        degrees[0] = 0;

        for (int i = 1; i <= max_deg; i++) {
            if (degrees[i] == 0) {
                return false;  // Not enough nodes
            }
            degrees[i]--;
        }
    }
    return true;
}

void realize(igraph_vector_int_t *ods, igraph_vector_int_t *ids, igraph_edge_type_sw_t method) {
    igraph_t graph;
    int err;

    err = igraph_realize_degree_sequence(&graph, ods, ids, et, method);

    if (err == IGRAPH_SUCCESS) {
        printf("\n");
        print_graph(&graph);
        igraph_destroy(&graph);
    } else if (err == IGRAPH_UNIMPLEMENTED) {
        printf("not implemented\n");
    } else {
        printf("not graphical\n");
    }
}

void realize2(igraph_vector_int_t *ods, igraph_vector_int_t *ids) {
    printf("Largest:");
    realize(ods, ids, et, IGRAPH_REALIZE_DEGSEQ_LARGEST);
    printf("Smallest:");
    realize(ods, ids, et, IGRAPH_REALIZE_DEGSEQ_SMALLEST);
    printf("Index:");
    realize(ods, ids, et, IGRAPH_REALIZE_DEGSEQ_INDEX);
}
