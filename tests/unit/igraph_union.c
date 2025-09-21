/*
   igraph library.
   Copyright (C) 2006-2023  The igraph development team <igraph@igraph.org>

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

#include <igraph.h>

#include <stdio.h>

#include "test_utilities.h"

void print_and_destroy_graph_and_maps(igraph_t *uni, igraph_vector_int_list_t *list) {
    printf("Union graph:\n");
    print_graph(uni);
    igraph_destroy(uni);

    printf("Edge maps:\n");
    print_vector_int_list(list);
    igraph_vector_int_list_clear(list);
}

int main(void) {
    igraph_t left, right, uni;
    igraph_vector_int_t edges;
    igraph_vector_ptr_t glist;
    igraph_vector_int_t edge_map1, edge_map2;
    igraph_vector_int_list_t edgemaps;

    printf("BINARY VERSION\n");

    igraph_vector_int_init(&edge_map1, 0);
    igraph_vector_int_init(&edge_map2, 0);

    igraph_small(&left, 0, IGRAPH_DIRECTED,
                 0, 1, 1, 2, 2, 2, 2, 3,
                 -1);

    igraph_small(&right, 0, IGRAPH_DIRECTED,
                 0, 1, 1, 2, 2, 2, 2, 4,
                 -1);

    igraph_union(&uni, &left, &right, &edge_map1, &edge_map2);
    printf("Union graph:\n");
    print_graph(&uni);
    printf("Edge maps:\n");
    print_vector_int(&edge_map1);
    print_vector_int(&edge_map2);

    igraph_destroy(&uni);
    igraph_destroy(&left);
    igraph_destroy(&right);

    igraph_vector_int_destroy(&edge_map1);
    igraph_vector_int_destroy(&edge_map2);

    printf("\n\nN-ARY VERSION\n");

    /* Empty graph list */

    printf("\nEmpty graph list:\n");
    igraph_vector_ptr_init(&glist, 0);
    igraph_vector_int_list_init(&edgemaps, 0);
    igraph_union_many(&uni, &glist, &edgemaps);
    print_and_destroy_graph_and_maps(&uni, &edgemaps);
    igraph_vector_ptr_destroy(&glist);

    /* Non-empty graph list */

    printf("\nNon-empty directed graph list 1:\n");
    igraph_vector_ptr_init(&glist, 10);
    for (igraph_int_t i = 0; i < igraph_vector_ptr_size(&glist); i++) {
        VECTOR(glist)[i] = IGRAPH_CALLOC(1, igraph_t);
        igraph_small(VECTOR(glist)[i], 0, IGRAPH_DIRECTED,
                     0, 1, 1, 0, -1);
    }

    igraph_union_many(&uni, &glist, &edgemaps);

    for (igraph_int_t i = 0; i < igraph_vector_ptr_size(&glist); i++) {
        igraph_destroy(VECTOR(glist)[i]);
        IGRAPH_FREE(VECTOR(glist)[i]);
    }

    print_and_destroy_graph_and_maps(&uni, &edgemaps);
    igraph_vector_ptr_destroy(&glist);

    /* Another non-empty graph list */

    printf("\nNon-empty directed graph list 2:\n");

    igraph_vector_ptr_init(&glist, 10);
    for (igraph_int_t i = 0; i < igraph_vector_ptr_size(&glist); i++) {
        VECTOR(glist)[i] = IGRAPH_CALLOC(1, igraph_t);
        igraph_vector_int_init(&edges, 4);
        VECTOR(edges)[0] = i; VECTOR(edges)[1] = i+1;
        VECTOR(edges)[2] = 1; VECTOR(edges)[3] = 0;
        igraph_create(VECTOR(glist)[i], &edges, 0, IGRAPH_DIRECTED);
        igraph_vector_int_destroy(&edges);
    }

    igraph_union_many(&uni, &glist, &edgemaps);
    igraph_write_graph_edgelist(&uni, stdout);

    for (igraph_int_t i = 0; i < igraph_vector_ptr_size(&glist); i++) {
        igraph_destroy(VECTOR(glist)[i]);
        IGRAPH_FREE(VECTOR(glist)[i]);
    }

    print_and_destroy_graph_and_maps(&uni, &edgemaps);
    igraph_vector_ptr_destroy(&glist);

    /* Undirected graph list*/

    printf("\nUndirected graph list:\n");

    igraph_vector_ptr_init(&glist, 10);
    for (igraph_int_t i = 0; i < igraph_vector_ptr_size(&glist); i++) {
        VECTOR(glist)[i] = IGRAPH_CALLOC(1, igraph_t);
        igraph_vector_int_init(&edges, 4);
        VECTOR(edges)[0] = i; VECTOR(edges)[1] = i+1;
        VECTOR(edges)[2] = 1; VECTOR(edges)[3] = 0;
        igraph_create(VECTOR(glist)[i], &edges, 0, IGRAPH_UNDIRECTED);
        igraph_vector_int_destroy(&edges);
    }

    igraph_union_many(&uni, &glist, &edgemaps);
    igraph_write_graph_edgelist(&uni, stdout);

    for (igraph_int_t i = 0; i < igraph_vector_ptr_size(&glist); i++) {
        igraph_destroy(VECTOR(glist)[i]);
        IGRAPH_FREE(VECTOR(glist)[i]);
    }

    print_and_destroy_graph_and_maps(&uni, &edgemaps);
    igraph_vector_ptr_destroy(&glist);

    igraph_vector_int_list_destroy(&edgemaps);

    VERIFY_FINALLY_STACK();

    return 0;
}
