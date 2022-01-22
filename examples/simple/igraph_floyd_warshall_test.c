




int main() {

    igraph_t g;
    igraph_vector_t weights;
    igraph_real_t weights_data_0[] = { 0, 2, 1, 0, 5, 2, 1, 1, 0, 2, 2, 8, 1, 1, 3, 1, 1, 4, 2, 1 };
    igraph_real_t weights_data_1[] = { 6, 7, 8, -4, -2, -3, 9, 2, 7 };
    igraph_real_t weights_data_2[] = { 6, 7, 2, -4, -2, -3, 9, 2, 7 };
    igraph_matrix_t res;

    /* Graph with only positive weights */
    igraph_small(&g, 10, IGRAPH_DIRECTED,
                 0, 1, 0, 2, 0, 3,    1, 2, 1, 4, 1, 5,
                 2, 3, 2, 6,         3, 2, 3, 6,
                 4, 5, 4, 7,         5, 6, 5, 8, 5, 9,
                 7, 5, 7, 8,         8, 9,
                 5, 2,
                 2, 1,
                 -1);
    
    igraph_vector_view(&weights, weights_data_0,
                       sizeof(weights_data_0) / sizeof(igraph_real_t));
    igraph_matrix_init(&res, 0, 0);
    
    igraph_shortest_paths_floyd_warshall(&g, &res,&weights, IGRAPH_OUT);
    print_matrix(&res);

    igraph_matrix_destroy(&res);
    igraph_destroy(&g);

    printf("\n");



}
