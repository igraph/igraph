/*
   igraph library.
   Copyright (C) 2026  The igraph development team <igraph@igraph.org>

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

int main(void) {
    igraph_t graph;
    igraph_bool_t triangle_free;

    igraph_setup();

    printf("Is the null graph triangle-free?");
    igraph_empty(&graph, 0, IGRAPH_UNDIRECTED);
    igraph_is_triangle_free(&graph, &triangle_free);
    printf(" %s\n", triangle_free ? "Yes." : "No.");
    igraph_destroy(&graph);

    printf("Is C_4 triangle-free?");
    igraph_cycle_graph(&graph, 4, IGRAPH_UNDIRECTED, false);
    igraph_is_triangle_free(&graph, &triangle_free);
    printf(" %s\n", triangle_free ? "Yes." : "No.");
    igraph_destroy(&graph);

    printf("Is K_4 triangle-free?");
    igraph_full(&graph, 4, IGRAPH_UNDIRECTED, false);
    igraph_is_triangle_free(&graph, &triangle_free);
    printf(" %s\n", triangle_free ? "Yes." : "No.");
    igraph_destroy(&graph);

    printf("Is K_{3,5} triangle-free?");
    igraph_full_bipartite(&graph, NULL, 3, 5, IGRAPH_UNDIRECTED, IGRAPH_ALL);
    igraph_is_triangle_free(&graph, &triangle_free);
    printf(" %s\n", triangle_free ? "Yes." : "No.");
    igraph_destroy(&graph);

    printf("Is a tree triangle-free?");
    igraph_kary_tree(&graph, 13, 3, IGRAPH_TREE_UNDIRECTED);
    igraph_is_triangle_free(&graph, &triangle_free);
    printf(" %s\n", triangle_free ? "Yes." : "No.");
    igraph_destroy(&graph);

    printf("Is the M_4 Mycielski graph triangle-free?");
    igraph_mycielski_graph(&graph, 4);
    igraph_is_triangle_free(&graph, &triangle_free);
    printf(" %s\n", triangle_free ? "Yes." : "No.");
    igraph_destroy(&graph);

    printf("Is Frucht graph triangle-free?");
    igraph_famous(&graph, "Frucht");
    igraph_is_triangle_free(&graph, &triangle_free);
    printf(" %s\n", triangle_free ? "Yes." : "No.");
    igraph_destroy(&graph);

    printf("Is Coxeter graph triangle-free?");
    igraph_famous(&graph, "Coxeter");
    igraph_is_triangle_free(&graph, &triangle_free);
    printf(" %s\n", triangle_free ? "Yes." : "No.");
    igraph_destroy(&graph);

    printf("Is the (2,3) Kautz graph triangle-free when ignoring edge directions?");
    igraph_kautz(&graph, 2, 3);
    igraph_is_triangle_free(&graph, &triangle_free);
    printf(" %s\n", triangle_free ? "Yes." : "No.");
    igraph_destroy(&graph);

    printf("Is the Q_4 hypercube triangle-free?");
    igraph_hypercube(&graph, 4, IGRAPH_UNDIRECTED);
    igraph_is_triangle_free(&graph, &triangle_free);
    printf(" %s\n", triangle_free ? "Yes." : "No.");
    igraph_destroy(&graph);

    return 0;
}
