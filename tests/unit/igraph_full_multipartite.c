/*
   IGraph library.
   Copyright (C) 2022  The igraph development team <igraph@igraph.org>

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
#include "test_utilities.h"

int main() {
    igraph_t g;
    igraph_vector_int_t partitions;
    igraph_vector_int_t types;

    printf("empty directed graph, zero vertices:");
    igraph_vector_int_init(&partitions, 0);
    igraph_vector_int_init(&types, 0);
    igraph_full_multipartite(&g, &types, &partitions, IGRAPH_DIRECTED, IGRAPH_ALL);

    printf("\nEdge list:\n");
    igraph_write_graph_edgelist(&g, stdout);

    printf("\nPartition type:\n");
    igraph_vector_int_print(&types);

    igraph_vector_int_destroy(&partitions);
    igraph_vector_int_destroy(&types);
    igraph_destroy(&g);

    printf("\ndirected graph with one partition, 4 vertices:");
    igraph_vector_int_init(&partitions, 1);
    igraph_vector_int_init(&types, 0);

    VECTOR(partitions)[0] = 4;

    igraph_full_multipartite(&g, &types, &partitions, IGRAPH_DIRECTED, IGRAPH_ALL);

    printf("\nEdge list:\n");
    igraph_write_graph_edgelist(&g, stdout);

    printf("\nPartition type:\n");
    igraph_vector_int_print(&types);
    
    igraph_vector_int_destroy(&partitions);
    igraph_vector_int_destroy(&types);
    igraph_destroy(&g);

    printf("\ndirected graph with 3 partitions:");
    igraph_vector_int_init(&partitions, 3);
    igraph_vector_int_init(&types, 0);

    VECTOR(partitions)[0] = 2;
    VECTOR(partitions)[1] = 3;
    VECTOR(partitions)[2] = 3;

    igraph_full_multipartite(&g, &types, &partitions, IGRAPH_DIRECTED, IGRAPH_ALL);

    printf("\nEdge list:\n");
    igraph_write_graph_edgelist(&g, stdout);

    printf("\nPartition type:\n");
    igraph_vector_int_print(&types);
    
    igraph_vector_int_destroy(&partitions);
    igraph_vector_int_destroy(&types);
    igraph_destroy(&g);

    printf("\ndirected graph, 4 partitions, mode=IN:");
    igraph_vector_int_init(&partitions, 4);
    igraph_vector_int_init(&types, 0);

    VECTOR(partitions)[0] = 2;
    VECTOR(partitions)[1] = 3;
    VECTOR(partitions)[2] = 4;
    VECTOR(partitions)[3] = 2;

    igraph_full_multipartite(&g, &types, &partitions, IGRAPH_DIRECTED, IGRAPH_IN);

    printf("\nEdge list:\n");
    igraph_write_graph_edgelist(&g, stdout);

    printf("\nPartition type:\n");
    igraph_vector_int_print(&types);
    
    igraph_vector_int_destroy(&partitions);
    igraph_vector_int_destroy(&types);
    igraph_destroy(&g);

    printf("\nUndirected graph, 4 partitions:");
    igraph_vector_int_init(&partitions, 4);
    igraph_vector_int_init(&types, 0);

    VECTOR(partitions)[0] = 2;
    VECTOR(partitions)[1] = 3;
    VECTOR(partitions)[2] = 4;
    VECTOR(partitions)[3] = 2;

    igraph_full_multipartite(&g, &types, &partitions, IGRAPH_UNDIRECTED, IGRAPH_ALL);

    printf("\nEdge list:\n");
    igraph_write_graph_edgelist(&g, stdout);

    printf("\nPartition type:\n");
    igraph_vector_int_print(&types);

    igraph_vector_int_destroy(&partitions);
    igraph_vector_int_destroy(&types);
    igraph_destroy(&g);

    printf("\ndirected graph with 3 partitions, all partitions with size zero:");
    igraph_vector_int_init(&partitions, 3);
    igraph_vector_int_init(&types, 0);

    VECTOR(partitions)[0] = 0;
    VECTOR(partitions)[1] = 0;
    VECTOR(partitions)[2] = 0;

    igraph_full_multipartite(&g, &types, &partitions, IGRAPH_DIRECTED, IGRAPH_ALL);

    printf("\nEdge list:\n");
    igraph_write_graph_edgelist(&g, stdout);

    printf("\nPartition type:\n");
    igraph_vector_int_print(&types);
    
    igraph_vector_int_destroy(&partitions);
    igraph_vector_int_destroy(&types);
    igraph_destroy(&g);

    printf("\ndirected graph with 3 partitions, one parition with size zero:");
    igraph_vector_int_init(&partitions, 3);
    igraph_vector_int_init(&types, 0);

    VECTOR(partitions)[0] = 2;
    VECTOR(partitions)[1] = 0;
    VECTOR(partitions)[2] = 3;

    igraph_full_multipartite(&g, &types, &partitions, IGRAPH_DIRECTED, IGRAPH_ALL);

    printf("\nEdge list:\n");
    igraph_write_graph_edgelist(&g, stdout);

    printf("\nPartition type:\n");
    igraph_vector_int_print(&types);
    
    igraph_vector_int_destroy(&partitions);
    igraph_vector_int_destroy(&types);
    igraph_destroy(&g);

    return IGRAPH_SUCCESS;
}
