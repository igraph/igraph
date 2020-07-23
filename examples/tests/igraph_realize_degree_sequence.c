
#include <igraph.h>
#include <stdio.h>

#include "test_utilities.inc"

int main() {
    igraph_t graph;

    {
        igraph_vector_t ds;
        const igraph_real_t rawds[] = { 3, 2, 2, 1 };
        igraph_vector_view(&ds, &rawds[0], sizeof(rawds) / sizeof(igraph_real_t));

        printf("\n");
        print_vector_round(&ds, stdout);
        printf("\n");

        igraph_realize_degree_sequence(&graph, &ds, NULL, IGRAPH_SIMPLE_SW, IGRAPH_REALIZE_DEGSEQ_LARGEST);
        print_graph(&graph, stdout);
        igraph_destroy(&graph);

        igraph_realize_degree_sequence(&graph, &ds, NULL, IGRAPH_SIMPLE_SW, IGRAPH_REALIZE_DEGSEQ_SMALLEST);
        print_graph(&graph, stdout);
        igraph_destroy(&graph);

        igraph_realize_degree_sequence(&graph, &ds, NULL, IGRAPH_SIMPLE_SW, IGRAPH_REALIZE_DEGSEQ_INDEX);
        print_graph(&graph, stdout);
        igraph_destroy(&graph);
    }

    {
        igraph_vector_t ds;
        const igraph_real_t rawds[] = {1, 3, 3, 4, 1, 2, 1, 1, 1, 3};
        igraph_vector_view(&ds, &rawds[0], sizeof(rawds) / sizeof(igraph_real_t));

        printf("\n");
        print_vector_round(&ds, stdout);
        printf("\n");

        igraph_realize_degree_sequence(&graph, &ds, NULL, IGRAPH_SIMPLE_SW, IGRAPH_REALIZE_DEGSEQ_LARGEST);
        print_graph(&graph, stdout);
        igraph_destroy(&graph);

        igraph_realize_degree_sequence(&graph, &ds, NULL, IGRAPH_SIMPLE_SW, IGRAPH_REALIZE_DEGSEQ_SMALLEST);
        print_graph(&graph, stdout);
        igraph_destroy(&graph);

        igraph_realize_degree_sequence(&graph, &ds, NULL, IGRAPH_SIMPLE_SW, IGRAPH_REALIZE_DEGSEQ_INDEX);
        print_graph(&graph, stdout);
        igraph_destroy(&graph);
    }

    {
        igraph_vector_t ds;
        const igraph_real_t rawds[] = {2, 0, 3, 2, 2, 2, 2, 3};
        igraph_vector_view(&ds, &rawds[0], sizeof(rawds) / sizeof(igraph_real_t));

        printf("\n");
        print_vector_round(&ds, stdout);
        printf("\n");

        igraph_realize_degree_sequence(&graph, &ds, NULL, IGRAPH_SIMPLE_SW, IGRAPH_REALIZE_DEGSEQ_LARGEST);
        print_graph(&graph, stdout);
        igraph_destroy(&graph);

        igraph_realize_degree_sequence(&graph, &ds, NULL, IGRAPH_SIMPLE_SW, IGRAPH_REALIZE_DEGSEQ_SMALLEST);
        print_graph(&graph, stdout);
        igraph_destroy(&graph);

        igraph_realize_degree_sequence(&graph, &ds, NULL, IGRAPH_SIMPLE_SW, IGRAPH_REALIZE_DEGSEQ_INDEX);
        print_graph(&graph, stdout);
        igraph_destroy(&graph);
    }

    {
        igraph_vector_t ods, ids;
        const igraph_real_t rawods[] = {3, 0, 1, 1, 1, 1, 0, 1};
        const igraph_real_t rawids[] = {2, 1, 0, 2, 2, 1, 0, 0};
        igraph_vector_view(&ods, &rawods[0], sizeof(rawods) / sizeof(igraph_real_t));
        igraph_vector_view(&ids, &rawids[0], sizeof(rawids) / sizeof(igraph_real_t));

        printf("\n");
        print_vector_round(&ods, stdout);
        print_vector_round(&ids, stdout);
        printf("\n");

        igraph_realize_degree_sequence(&graph, &ods, &ids, IGRAPH_SIMPLE_SW, IGRAPH_REALIZE_DEGSEQ_LARGEST);
        print_graph(&graph, stdout);
        igraph_destroy(&graph);

        igraph_realize_degree_sequence(&graph, &ods, &ids, IGRAPH_SIMPLE_SW, IGRAPH_REALIZE_DEGSEQ_SMALLEST);
        print_graph(&graph, stdout);
        igraph_destroy(&graph);

        igraph_realize_degree_sequence(&graph, &ods, &ids, IGRAPH_SIMPLE_SW, IGRAPH_REALIZE_DEGSEQ_INDEX);
        print_graph(&graph, stdout);
        igraph_destroy(&graph);
    }

    {
        igraph_vector_t ods, ids;
        const igraph_real_t rawods[] = {3, 1, 2, 3, 1, 2, 2};
        const igraph_real_t rawids[] = {2, 2, 1, 2, 3, 2, 2};
        igraph_vector_view(&ods, &rawods[0], sizeof(rawods) / sizeof(igraph_real_t));
        igraph_vector_view(&ids, &rawids[0], sizeof(rawids) / sizeof(igraph_real_t));

        printf("\n");
        print_vector_round(&ods, stdout);
        print_vector_round(&ids, stdout);
        printf("\n");

        igraph_realize_degree_sequence(&graph, &ods, &ids, IGRAPH_SIMPLE_SW, IGRAPH_REALIZE_DEGSEQ_LARGEST);
        print_graph(&graph, stdout);
        igraph_destroy(&graph);

        igraph_realize_degree_sequence(&graph, &ods, &ids, IGRAPH_SIMPLE_SW, IGRAPH_REALIZE_DEGSEQ_SMALLEST);
        print_graph(&graph, stdout);
        igraph_destroy(&graph);

        igraph_realize_degree_sequence(&graph, &ods, &ids, IGRAPH_SIMPLE_SW, IGRAPH_REALIZE_DEGSEQ_INDEX);
        print_graph(&graph, stdout);
        igraph_destroy(&graph);
    }

    return 0;
}
