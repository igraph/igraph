/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2012  Gabor Csardi <csardi.gabor@gmail.com>
   334 Harvard street, Cambridge, MA 02139 USA

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
   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
   02110-1301 USA

*/

#include <igraph.h>

#include "test_utilities.inc"

int print_attributes(const igraph_t *g) {

    igraph_vector_t gtypes, vtypes, etypes;
    igraph_strvector_t gnames, vnames, enames;
    long int i;

    igraph_vector_init(&gtypes, 0);
    igraph_vector_init(&vtypes, 0);
    igraph_vector_init(&etypes, 0);
    igraph_strvector_init(&gnames, 0);
    igraph_strvector_init(&vnames, 0);
    igraph_strvector_init(&enames, 0);

    igraph_cattribute_list(g, &gnames, &gtypes, &vnames, &vtypes,
                           &enames, &etypes);

    for (i = 0; i < igraph_vcount(g); i++) {
        long int j;
        printf("Vertex %li: ", i);
        for (j = 0; j < igraph_strvector_size(&vnames); j++) {
            printf("%s=", STR(vnames, j));
            if (VECTOR(vtypes)[j] == IGRAPH_ATTRIBUTE_NUMERIC) {
                igraph_real_printf(VAN(g, STR(vnames, j), i));
                putchar(' ');
            } else {
                printf("\"%s\" ", VAS(g, STR(vnames, j), i));
            }
        }
        printf("\n");
    }

    for (i = 0; i < igraph_ecount(g); i++) {
        long int j;
        int u = IGRAPH_FROM(g, i), v = IGRAPH_TO(g, i);
        if (u < v && !igraph_is_directed(g)) {
            u = IGRAPH_TO(g, i);
            v = IGRAPH_FROM(g, i);
        }
        printf("Edge %li (%i-%i): ", i, u, v);
        for (j = 0; j < igraph_strvector_size(&enames); j++) {
            printf("%s=", STR(enames, j));
            if (VECTOR(etypes)[j] == IGRAPH_ATTRIBUTE_NUMERIC) {
                igraph_real_printf(EAN(g, STR(enames, j), i));
                putchar(' ');
            } else {
                printf("\"%s\" ", EAS(g, STR(enames, j), i));
            }
        }
        printf("\n");
    }

    igraph_strvector_destroy(&enames);
    igraph_strvector_destroy(&vnames);
    igraph_strvector_destroy(&gnames);
    igraph_vector_destroy(&etypes);
    igraph_vector_destroy(&vtypes);
    igraph_vector_destroy(&gtypes);

    return 0;
}

int main() {
    igraph_t graph;
    FILE *input;

    /* turn on attribute handling */
    igraph_set_attribute_table(&igraph_cattribute_table);

    /* first file, without marginals */

    input = fopen("pajek_bip.net", "r");
    if (input == 0) {
        return 1;
    }

    igraph_read_graph_pajek(&graph, input);
    fclose(input);

    print_attributes(&graph);

    igraph_destroy(&graph);

    /* second file, with marginals */

    printf("---\n");

    input = fopen("pajek_bip2.net", "r");
    if (input == 0) {
        return 1;
    }

    igraph_read_graph_pajek(&graph, input);
    fclose(input);

    print_attributes(&graph);

    igraph_destroy(&graph);

    VERIFY_FINALLY_STACK();

    return 0;
}
