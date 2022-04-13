/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2007-2012  Gabor Csardi <csardi.gabor@gmail.com>
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
   Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA
   02110-1301 USA

*/

#include <igraph.h>
#include <string.h>
#include <stdlib.h>

int print_attributes(const igraph_t *g);

int main() {

    igraph_t g;
    igraph_vector_t y;

    /* turn on attribute handling */
    igraph_set_attribute_table(&igraph_cattribute_table);

    igraph_small(&g, 3, IGRAPH_DIRECTED, 0,1, 1,2, -1);

    /* Set graph attributes */
    /* numeric */
    SETGAN(&g, "id", 10);
    /* string */
    SETGAS(&g, "name", "toy");
    /* boolean */
    SETGAB(&g, "is_regular", 0);

    /* Set edge string attribute */
    SETEAS(&g, "color", 1, "RED");

    /* Set vertex attributes as vector */
    igraph_vector_init(&y, igraph_vcount(&g));
    igraph_vector_fill(&y, 1.23);
    SETVANV(&g, "y", &y);
    igraph_vector_destroy(&y);

    /* Set single vertex numeric attribute */
    SETVAN(&g, "y", 0, -1);

    /* Delete graph attribute */
    DELGA(&g, "is_regular");

    /* Print the final result */
    print_attributes(&g);

    /* Delete all remaining attributes */
    DELALL(&g);

    /* Destroy the graph */
    igraph_destroy(&g);
    return 0;
}

int print_attributes(const igraph_t *g) {

    igraph_vector_int_t gtypes, vtypes, etypes;
    igraph_strvector_t gnames, vnames, enames;
    igraph_integer_t i;

    igraph_integer_t j;

    igraph_vector_int_init(&gtypes, 0);
    igraph_vector_int_init(&vtypes, 0);
    igraph_vector_int_init(&etypes, 0);
    igraph_strvector_init(&gnames, 0);
    igraph_strvector_init(&vnames, 0);
    igraph_strvector_init(&enames, 0);

    igraph_cattribute_list(g, &gnames, &gtypes, &vnames, &vtypes,
                           &enames, &etypes);

    /* Graph attributes */
    for (i = 0; i < igraph_strvector_size(&gnames); i++) {
        printf("%s=", STR(gnames, i));
        if (VECTOR(gtypes)[i] == IGRAPH_ATTRIBUTE_NUMERIC) {
            igraph_real_printf(GAN(g, STR(gnames, i)));
            putchar(' ');
        } else {
            printf("\"%s\" ", GAS(g, STR(gnames, i)));
        }
    }
    printf("\n");

    for (i = 0; i < igraph_vcount(g); i++) {
        printf("Vertex %" IGRAPH_PRId ": ", i);
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
        printf("Edge %" IGRAPH_PRId " (%" IGRAPH_PRId "-%" IGRAPH_PRId "): ", i, IGRAPH_FROM(g, i), IGRAPH_TO(g, i));
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
    igraph_vector_int_destroy(&etypes);
    igraph_vector_int_destroy(&vtypes);
    igraph_vector_int_destroy(&gtypes);

    return 0;
}
