
#include <igraph.h>
#include <stdio.h>

void print_edges(const igraph_t *graph) {
    long ecount = igraph_ecount(graph);
    long i;

    for (i = 0; i < ecount; ++i) {
        printf("%d %d\n", IGRAPH_FROM(graph, i), IGRAPH_TO(graph, i));
    }
    printf("\n");
}

void print_vector(igraph_vector_t *v) {
    long int i, n;

    n = igraph_vector_size(v);
    for (i = 0; i < n; i++)
        if (i != n - 1 ) {
            printf("%li ", (long int) VECTOR(*v)[i]);
        } else {
            printf("%li", (long int) VECTOR(*v)[i]);
        }
    printf("\n");
}

int main() {
    igraph_t graph;

    {
        igraph_vector_t ds;
        const igraph_real_t rawds[] = { 3, 2, 2, 1 };
        igraph_vector_view(&ds, &rawds[0], sizeof(rawds) / sizeof(igraph_real_t));

        print_vector(&ds);
        printf("\n");

        igraph_realize_degree_sequence(&graph, &ds, NULL, IGRAPH_REALIZE_DEGSEQ_LARGEST);
        print_edges(&graph);
        igraph_destroy(&graph);

        igraph_realize_degree_sequence(&graph, &ds, NULL, IGRAPH_REALIZE_DEGSEQ_SMALLEST);
        print_edges(&graph);
        igraph_destroy(&graph);

        igraph_realize_degree_sequence(&graph, &ds, NULL, IGRAPH_REALIZE_DEGSEQ_INDEX);
        print_edges(&graph);
        igraph_destroy(&graph);
    }

    {
        igraph_vector_t ds;
        const igraph_real_t rawds[] = {1, 3, 3, 4, 1, 2, 1, 1, 1, 3};
        igraph_vector_view(&ds, &rawds[0], sizeof(rawds) / sizeof(igraph_real_t));

        print_vector(&ds);
        printf("\n");

        igraph_realize_degree_sequence(&graph, &ds, NULL, IGRAPH_REALIZE_DEGSEQ_LARGEST);
        print_edges(&graph);
        igraph_destroy(&graph);

        igraph_realize_degree_sequence(&graph, &ds, NULL, IGRAPH_REALIZE_DEGSEQ_SMALLEST);
        print_edges(&graph);
        igraph_destroy(&graph);

        igraph_realize_degree_sequence(&graph, &ds, NULL, IGRAPH_REALIZE_DEGSEQ_INDEX);
        print_edges(&graph);
        igraph_destroy(&graph);
    }

    {
        igraph_vector_t ds;
        const igraph_real_t rawds[] = {2, 0, 3, 2, 2, 2, 2, 3};
        igraph_vector_view(&ds, &rawds[0], sizeof(rawds) / sizeof(igraph_real_t));

        print_vector(&ds);
        printf("\n");

        igraph_realize_degree_sequence(&graph, &ds, NULL, IGRAPH_REALIZE_DEGSEQ_LARGEST);
        print_edges(&graph);
        igraph_destroy(&graph);

        igraph_realize_degree_sequence(&graph, &ds, NULL, IGRAPH_REALIZE_DEGSEQ_SMALLEST);
        print_edges(&graph);
        igraph_destroy(&graph);

        igraph_realize_degree_sequence(&graph, &ds, NULL, IGRAPH_REALIZE_DEGSEQ_INDEX);
        print_edges(&graph);
        igraph_destroy(&graph);
    }

    {
        igraph_vector_t ods, ids;
        const igraph_real_t rawods[] = {3, 0, 1, 1, 1, 1, 0, 1};
        const igraph_real_t rawids[] = {2, 1, 0, 2, 2, 1, 0, 0};
        igraph_vector_view(&ods, &rawods[0], sizeof(rawods) / sizeof(igraph_real_t));
        igraph_vector_view(&ids, &rawids[0], sizeof(rawids) / sizeof(igraph_real_t));

        print_vector(&ods);
        print_vector(&ids);
        printf("\n");

        igraph_realize_degree_sequence(&graph, &ods, &ids, IGRAPH_REALIZE_DEGSEQ_LARGEST);
        print_edges(&graph);
        igraph_destroy(&graph);

        igraph_realize_degree_sequence(&graph, &ods, &ids, IGRAPH_REALIZE_DEGSEQ_SMALLEST);
        print_edges(&graph);
        igraph_destroy(&graph);

        igraph_realize_degree_sequence(&graph, &ods, &ids, IGRAPH_REALIZE_DEGSEQ_INDEX);
        print_edges(&graph);
        igraph_destroy(&graph);
    }

    {
        igraph_vector_t ods, ids;
        const igraph_real_t rawods[] = {3, 1, 2, 3, 1, 2, 2};
        const igraph_real_t rawids[] = {2, 2, 1, 2, 3, 2, 2};
        igraph_vector_view(&ods, &rawods[0], sizeof(rawods) / sizeof(igraph_real_t));
        igraph_vector_view(&ids, &rawids[0], sizeof(rawids) / sizeof(igraph_real_t));

        print_vector(&ods);
        print_vector(&ids);
        printf("\n");

        igraph_realize_degree_sequence(&graph, &ods, &ids, IGRAPH_REALIZE_DEGSEQ_LARGEST);
        print_edges(&graph);
        igraph_destroy(&graph);

        igraph_realize_degree_sequence(&graph, &ods, &ids, IGRAPH_REALIZE_DEGSEQ_SMALLEST);
        print_edges(&graph);
        igraph_destroy(&graph);

        igraph_realize_degree_sequence(&graph, &ods, &ids, IGRAPH_REALIZE_DEGSEQ_INDEX);
        print_edges(&graph);
        igraph_destroy(&graph);
    }

    return 0;
}
