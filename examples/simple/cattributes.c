/* igraph library.
   Copyright (C) 2007-2021  The igraph development team <igraph@igraph.org>

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
#include <string.h>
#include <stdlib.h>

/* Prints graph, vertex and edge attributes stored in a graph. */
void print_attributes(const igraph_t *g) {
    igraph_vector_int_t gtypes, vtypes, etypes;
    igraph_strvector_t gnames, vnames, enames;
    igraph_int_t i, j;

    igraph_vector_int_init(&gtypes, 0);
    igraph_vector_int_init(&vtypes, 0);
    igraph_vector_int_init(&etypes, 0);
    igraph_strvector_init(&gnames, 0);
    igraph_strvector_init(&vnames, 0);
    igraph_strvector_init(&enames, 0);

    igraph_cattribute_list(g,
                           &gnames, &gtypes,
                           &vnames, &vtypes,
                           &enames, &etypes);

    /* graph attributes */
    for (i = 0; i < igraph_strvector_size(&gnames); i++) {
        printf("%s=", igraph_strvector_get(&gnames, i));
        if (VECTOR(gtypes)[i] == IGRAPH_ATTRIBUTE_NUMERIC) {
            igraph_real_printf(GAN(g, igraph_strvector_get(&gnames, i)));
            putchar(' ');
        } else {
            printf("\"%s\" ", GAS(g, igraph_strvector_get(&gnames, i)));
        }
    }
    printf("\n");

    /* vertex attributes */
    for (i = 0; i < igraph_vcount(g); i++) {
        printf("Vertex %" IGRAPH_PRId ": ", i);
        for (j = 0; j < igraph_strvector_size(&vnames); j++) {
            printf("%s=", igraph_strvector_get(&vnames, j));
            if (VECTOR(vtypes)[j] == IGRAPH_ATTRIBUTE_NUMERIC) {
                igraph_real_printf(VAN(g, igraph_strvector_get(&vnames, j), i));
                putchar(' ');
            } else {
                printf("\"%s\" ", VAS(g, igraph_strvector_get(&vnames, j), i));
            }
        }
        printf("\n");
    }

    /* edge attributes */
    for (i = 0; i < igraph_ecount(g); i++) {
        printf("Edge %" IGRAPH_PRId " (%" IGRAPH_PRId "-%" IGRAPH_PRId "): ", i, IGRAPH_FROM(g, i), IGRAPH_TO(g, i));
        for (j = 0; j < igraph_strvector_size(&enames); j++) {
            printf("%s=", igraph_strvector_get(&enames, j));
            if (VECTOR(etypes)[j] == IGRAPH_ATTRIBUTE_NUMERIC) {
                igraph_real_printf(EAN(g, igraph_strvector_get(&enames, j), i));
                putchar(' ');
            } else {
                printf("\"%s\" ", EAS(g, igraph_strvector_get(&enames, j), i));
            }
        }
        printf("\n");
    }

    igraph_strvector_destroy(&enames);
    igraph_strvector_destroy(&vnames);
    igraph_strvector_destroy(&gnames);
    igraph_vector_int_destroy(&etypes);
    igraph_vector_int_destroy(&vtypes);
    igraph_vector_int_destroy(&gtypes);
}

int main(void) {
    igraph_t graph;
    igraph_vector_t y;

    /* Initialize the library. */
    igraph_setup();

    /* Turn on attribute handling. */
    igraph_set_attribute_table(&igraph_cattribute_table);

    igraph_small(&graph, 3, IGRAPH_DIRECTED, 0,1, 1,2, -1);

    /* Set graph attributes. */
    /* numeric */
    SETGAN(&graph, "id", 10);
    /* string */
    SETGAS(&graph, "name", "toy");
    /* boolean */
    SETGAB(&graph, "is_regular", false);

    /* Set edge string attribute. */
    SETEAS(&graph, "color", 1, "RED");

    /* Set vertex attributes as vector. */
    igraph_vector_init(&y, igraph_vcount(&graph));
    igraph_vector_fill(&y, 1.23);
    SETVANV(&graph, "y", &y);
    igraph_vector_destroy(&y);

    /* Set single vertex numeric attribute. */
    SETVAN(&graph, "y", 0, -1);

    /* Delete graph attribute. */
    DELGA(&graph, "is_regular");

    /* Print the final result. */
    print_attributes(&graph);

    /* Delete all remaining attributes. */
    DELALL(&graph);

    /* Destroy the graph. */
    igraph_destroy(&graph);

    return 0;
}
