/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2006-2012  Gabor Csardi <csardi.gabor@gmail.com>
   334 Harvard st, Cambridge MA, 02139 USA

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA
   02110-1301 USA

*/
#include <igraph.h>
#include <stdio.h>

#include "test_utilities.inc"

int main() {
    int i;
    igraph_t g;
    igraph_vector_t values;
    igraph_strvector_t strvalues;
    const char* strings[] = {"foo", "bar", "baz", "spam", "eggs", "bacon"};

    /* Setting up attribute handler */
    igraph_set_attribute_table(&igraph_cattribute_table);

    /* Saving directed graph, no attributes */
    igraph_ring(&g, 5, /* directed = */ 1,
                /* mutual   = */ 0,
                /* circular = */ 1);
    igraph_write_graph_leda(&g, stdout, 0, 0);
    printf("===\n");
    igraph_destroy(&g);

    /* Saving undirected graph, no attributes */
    igraph_ring(&g, 5, /* directed = */ 0,
                /* mutual   = */ 0,
                /* circular = */ 1);
    igraph_write_graph_leda(&g, stdout, 0, 0);
    printf("===\n");
    igraph_destroy(&g);

    /* Saving directed graph with vertex attributes */
    igraph_ring(&g, 5, /* directed = */ 1,
                /* mutual   = */ 0,
                /* circular = */ 1);
    igraph_vector_init_seq(&values, 5, 9);
    SETVANV(&g, "name", &values);
    igraph_write_graph_leda(&g, stdout, "name", 0);
    igraph_vector_destroy(&values);
    printf("===\n");
    DELVAS(&g);
    igraph_strvector_init(&strvalues, 5);
    for (i = 0; i < 5; i++) {
        igraph_strvector_set(&strvalues, i, strings[i]);
    }
    SETVASV(&g, "name", &strvalues);
    igraph_write_graph_leda(&g, stdout, "name", 0);
    igraph_strvector_destroy(&strvalues);
    printf("===\n");
    igraph_destroy(&g);

    /* Saving undirected graph with edge attributes */
    igraph_ring(&g, 5, /* directed = */ 0,
                /* mutual   = */ 0,
                /* circular = */ 1);
    igraph_vector_init_seq(&values, 5, 9);
    SETEANV(&g, "weight", &values);
    igraph_write_graph_leda(&g, stdout, 0, "weight");
    igraph_vector_destroy(&values);
    printf("===\n");
    DELEAS(&g);
    igraph_strvector_init(&strvalues, 5);
    for (i = 0; i < 5; i++) {
        igraph_strvector_set(&strvalues, i, strings[i]);
    }
    SETEASV(&g, "weight", &strvalues);
    igraph_write_graph_leda(&g, stdout, 0, "weight");
    igraph_strvector_destroy(&strvalues);
    printf("===\n");
    igraph_destroy(&g);

    /* Saving undirected graph with edge attributes and large weights */
    igraph_ring(&g, 5, /* directed = */ 0,
                /* mutual   = */ 0,
                /* circular = */ 1);
    igraph_vector_init_seq(&values, 123456789, 123456793);
    SETEANV(&g, "weight", &values);
    igraph_write_graph_leda(&g, stdout, 0, "weight");
    igraph_vector_destroy(&values);
    printf("===\n");
    igraph_destroy(&g);

    VERIFY_FINALLY_STACK();

    return 0;
}
