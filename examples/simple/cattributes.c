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

int print_attributes(const igraph_t *g) {

    igraph_vector_t gtypes, vtypes, etypes;
    igraph_strvector_t gnames, vnames, enames;
    long int i;

    igraph_vector_t vec;
    igraph_strvector_t svec;
    long int j;

    igraph_vector_init(&gtypes, 0);
    igraph_vector_init(&vtypes, 0);
    igraph_vector_init(&etypes, 0);
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
        printf("Edge %li (%i-%i): ", i, (int)IGRAPH_FROM(g, i), (int)IGRAPH_TO(g, i));
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

    /* Check vector-based query functions */
    igraph_vector_init(&vec, 0);
    igraph_strvector_init(&svec, 0);

    for (j = 0; j < igraph_strvector_size(&vnames); j++) {
        if (VECTOR(vtypes)[j] == IGRAPH_ATTRIBUTE_NUMERIC) {
            igraph_cattribute_VANV(g, STR(vnames, j), igraph_vss_all(), &vec);
            for (i = 0; i < igraph_vcount(g); i++) {
                igraph_real_t num = VAN(g, STR(vnames, j), i);
                if (num != VECTOR(vec)[i] &&
                    (!isnan(num) || !isnan(VECTOR(vec)[i]))) {
                    exit(51);
                }
            }
        } else {
            igraph_cattribute_VASV(g, STR(vnames, j), igraph_vss_all(), &svec);
            for (i = 0; i < igraph_vcount(g); i++) {
                const char *str = VAS(g, STR(vnames, j), i);
                if (strcmp(str, STR(svec, i))) {
                    exit(52);
                }
            }
        }
    }

    for (j = 0; j < igraph_strvector_size(&enames); j++) {
        if (VECTOR(etypes)[j] == IGRAPH_ATTRIBUTE_NUMERIC) {
            igraph_cattribute_EANV(g, STR(enames, j),
                                   igraph_ess_all(IGRAPH_EDGEORDER_ID), &vec);
            for (i = 0; i < igraph_ecount(g); i++) {
                igraph_real_t num = EAN(g, STR(enames, j), i);
                if (num != VECTOR(vec)[i] &&
                    (!isnan(num) || !isnan(VECTOR(vec)[i]))) {
                    exit(53);
                }
            }
        } else {
            igraph_cattribute_EASV(g, STR(enames, j),
                                   igraph_ess_all(IGRAPH_EDGEORDER_ID), &svec);
            for (i = 0; i < igraph_ecount(g); i++) {
                const char *str = EAS(g, STR(enames, j), i);
                if (strcmp(str, STR(svec, i))) {
                    exit(54);
                }
            }
        }
    }

    igraph_strvector_destroy(&svec);
    igraph_vector_destroy(&vec);

    igraph_strvector_destroy(&enames);
    igraph_strvector_destroy(&vnames);
    igraph_strvector_destroy(&gnames);
    igraph_vector_destroy(&etypes);
    igraph_vector_destroy(&vtypes);
    igraph_vector_destroy(&gtypes);

    return 0;
}

int main() {

    igraph_t g, g2;
    FILE *ifile;
    igraph_vector_t gtypes, vtypes, etypes;
    igraph_strvector_t gnames, vnames, enames;
    long int i;
    igraph_vector_t y;
    igraph_strvector_t id;
    igraph_vector_bool_t type;
    char str[21];

    /* turn on attribute handling */
    igraph_set_attribute_table(&igraph_cattribute_table);

    ifile = fopen("links.net", "r");
    if (ifile == 0) {
        return 10;
    }
    igraph_read_graph_pajek(&g, ifile);
    fclose(ifile);

    igraph_vector_init(&gtypes, 0);
    igraph_vector_init(&vtypes, 0);
    igraph_vector_init(&etypes, 0);
    igraph_strvector_init(&gnames, 0);
    igraph_strvector_init(&vnames, 0);
    igraph_strvector_init(&enames, 0);

    igraph_cattribute_list(&g, &gnames, &gtypes, &vnames, &vtypes,
                           &enames, &etypes);

    /* List attribute names and types */
    printf("Graph attributes: ");
    for (i = 0; i < igraph_strvector_size(&gnames); i++) {
        printf("%s (%i) ", STR(gnames, i), (int)VECTOR(gtypes)[i]);
    }
    printf("\n");
    printf("Vertex attributes: ");
    for (i = 0; i < igraph_strvector_size(&vnames); i++) {
        printf("%s (%i) ", STR(vnames, i), (int)VECTOR(vtypes)[i]);
    }
    printf("\n");
    printf("Edge attributes: ");
    for (i = 0; i < igraph_strvector_size(&enames); i++) {
        printf("%s (%i) ", STR(enames, i), (int)VECTOR(etypes)[i]);
    }
    printf("\n");

    print_attributes(&g);

    /* Copying a graph */
    igraph_copy(&g2, &g);
    print_attributes(&g2);
    igraph_destroy(&g2);

    /* Adding vertices */
    igraph_add_vertices(&g, 3, 0);
    print_attributes(&g);

    /* Adding edges */
    igraph_add_edge(&g, 1, 1);
    igraph_add_edge(&g, 2, 5);
    igraph_add_edge(&g, 3, 6);
    print_attributes(&g);

    /* Deleting vertices */
    igraph_delete_vertices(&g, igraph_vss_1(1));
    igraph_delete_vertices(&g, igraph_vss_1(4));
    print_attributes(&g);

    /* Deleting edges */
    igraph_delete_edges(&g, igraph_ess_1(igraph_ecount(&g) - 1));
    igraph_delete_edges(&g, igraph_ess_1(0));
    print_attributes(&g);

    /* Set graph attributes */
    SETGAN(&g, "id", 10);
    if (GAN(&g, "id") != 10) {
        return 11;
    }
    SETGAS(&g, "name", "toy");
    if (strcmp(GAS(&g, "name"), "toy")) {
        return 12;
    }
    SETGAB(&g, "is_regular", 0);
    if (GAB(&g, "is_regular") != 0) {
        return 13;
    }

    /* Delete graph attributes */
    DELGA(&g, "id");
    DELGA(&g, "name");
    DELGA(&g, "is_regular");
    igraph_cattribute_list(&g, &gnames, 0, 0, 0, 0, 0);
    if (igraph_strvector_size(&gnames) != 0) {
        return 14;
    }

    /* Delete vertex attributes */
    DELVA(&g, "x");
    DELVA(&g, "shape");
    DELVA(&g, "xfact");
    DELVA(&g, "yfact");
    igraph_cattribute_list(&g, 0, 0, &vnames, 0, 0, 0);
    if (igraph_strvector_size(&vnames) != 3) {
        return 15;
    }

    /* Delete edge attributes */
    igraph_cattribute_list(&g, 0, 0, 0, 0, &enames, 0);
    i = igraph_strvector_size(&enames);
    DELEA(&g, "hook1");
    DELEA(&g, "hook2");
    DELEA(&g, "label");
    igraph_cattribute_list(&g, 0, 0, 0, 0, &enames, 0);
    if (igraph_strvector_size(&enames) != i - 3) {
        return 16;
    }

    /* Set vertex attributes */
    SETVAN(&g, "y", 0, -1);
    SETVAN(&g, "y", 1, 2.1);
    if (VAN(&g, "y", 0) != -1 ||
        VAN(&g, "y", 1) != 2.1) {
        return 17;
    }
    SETVAS(&g, "id", 0, "foo");
    SETVAS(&g, "id", 1, "bar");
    if (strcmp(VAS(&g, "id", 0), "foo") ||
        strcmp(VAS(&g, "id", 1), "bar")) {
        return 18;
    }
    SETVAB(&g, "type", 0, 1);
    SETVAB(&g, "type", 1, 0);
    if (!VAB(&g, "type", 0) || VAB(&g, "type", 1)) {
        return 26;
    }

    /* Set edge attributes */
    SETEAN(&g, "weight", 2, 100.0);
    SETEAN(&g, "weight", 0, -100.1);
    if (EAN(&g, "weight", 2) != 100.0 ||
        EAN(&g, "weight", 0) != -100.1) {
        return 19;
    }
    SETEAS(&g, "color", 2, "RED");
    SETEAS(&g, "color", 0, "Blue");
    if (strcmp(EAS(&g, "color", 2), "RED") ||
        strcmp(EAS(&g, "color", 0), "Blue")) {
        return 20;
    }
    SETEAB(&g, "type", 0, 1);
    SETEAB(&g, "type", 2, 0);
    if (!EAB(&g, "type", 0) || EAB(&g, "type", 2)) {
        return 27;
    }

    /* Set vertex attributes as vector */
    igraph_vector_init(&y, igraph_vcount(&g));
    igraph_vector_fill(&y, 1.23);
    SETVANV(&g, "y", &y);
    igraph_vector_destroy(&y);
    for (i = 0; i < igraph_vcount(&g); i++) {
        if (VAN(&g, "y", i) != 1.23) {
            return 21;
        }
    }
    igraph_vector_init_seq(&y, 0, igraph_vcount(&g) - 1);
    SETVANV(&g, "foobar", &y);
    igraph_vector_destroy(&y);
    for (i = 0; i < igraph_vcount(&g); i++) {
        if (VAN(&g, "foobar", i) != i) {
            return 22;
        }
    }

    igraph_vector_bool_init(&type, igraph_vcount(&g));
    for (i = 0; i < igraph_vcount(&g); i++) {
        VECTOR(type)[i] = (i % 2 == 1);
    }
    SETVABV(&g, "type", &type);
    igraph_vector_bool_destroy(&type);
    for (i = 0; i < igraph_vcount(&g); i++) {
        if (VAB(&g, "type", i) != (i % 2 == 1)) {
            return 28;
        }
    }

    igraph_strvector_init(&id, igraph_vcount(&g));
    for (i = 0; i < igraph_vcount(&g); i++) {
        snprintf(str, sizeof(str) - 1, "%li", i);
        igraph_strvector_set(&id, i, str);
    }
    SETVASV(&g, "foo", &id);
    igraph_strvector_destroy(&id);
    for (i = 0; i < igraph_vcount(&g); i++) {
        printf("%s ", VAS(&g, "foo", i));
    }
    printf("\n");
    igraph_strvector_init(&id, igraph_vcount(&g));
    for (i = 0; i < igraph_vcount(&g); i++) {
        snprintf(str, sizeof(str) - 1, "%li", i);
        igraph_strvector_set(&id, i, str);
    }
    SETVASV(&g, "id", &id);
    igraph_strvector_destroy(&id);
    for (i = 0; i < igraph_vcount(&g); i++) {
        printf("%s ", VAS(&g, "id", i));
    }
    printf("\n");

    /* Set edge attributes as vector */
    igraph_vector_init(&y, igraph_ecount(&g));
    igraph_vector_fill(&y, 12.3);
    SETEANV(&g, "weight", &y);
    igraph_vector_destroy(&y);
    for (i = 0; i < igraph_ecount(&g); i++) {
        if (EAN(&g, "weight", i) != 12.3) {
            return 23;
        }
    }
    igraph_vector_init_seq(&y, 0, igraph_ecount(&g) - 1);
    SETEANV(&g, "foobar", &y);
    igraph_vector_destroy(&y);
    for (i = 0; i < igraph_ecount(&g); i++) {
        if (VAN(&g, "foobar", i) != i) {
            return 24;
        }
    }

    igraph_vector_bool_init(&type, igraph_ecount(&g));
    for (i = 0; i < igraph_ecount(&g); i++) {
        VECTOR(type)[i] = (i % 2 == 1);
    }
    SETEABV(&g, "type", &type);
    igraph_vector_bool_destroy(&type);
    for (i = 0; i < igraph_ecount(&g); i++) {
        if (EAB(&g, "type", i) != (i % 2 == 1)) {
            return 29;
        }
    }

    igraph_strvector_init(&id, igraph_ecount(&g));
    for (i = 0; i < igraph_ecount(&g); i++) {
        snprintf(str, sizeof(str) - 1, "%li", i);
        igraph_strvector_set(&id, i, str);
    }
    SETEASV(&g, "foo", &id);
    igraph_strvector_destroy(&id);
    for (i = 0; i < igraph_ecount(&g); i++) {
        printf("%s ", EAS(&g, "foo", i));
    }
    printf("\n");
    igraph_strvector_init(&id, igraph_ecount(&g));
    for (i = 0; i < igraph_ecount(&g); i++) {
        snprintf(str, sizeof(str) - 1, "%li", i);
        igraph_strvector_set(&id, i, str);
    }
    SETEASV(&g, "color", &id);
    igraph_strvector_destroy(&id);
    for (i = 0; i < igraph_ecount(&g); i++) {
        printf("%s ", EAS(&g, "color", i));
    }
    printf("\n");

    /* Delete all remaining attributes */
    DELALL(&g);
    igraph_cattribute_list(&g, &gnames, &gtypes, &vnames, &vtypes, &enames, &etypes);
    if (igraph_strvector_size(&gnames) != 0 ||
        igraph_strvector_size(&vnames) != 0 ||
        igraph_strvector_size(&enames) != 0) {
        return 25;
    }

    /* Destroy */
    igraph_vector_destroy(&gtypes);
    igraph_vector_destroy(&vtypes);
    igraph_vector_destroy(&etypes);
    igraph_strvector_destroy(&gnames);
    igraph_strvector_destroy(&vnames);
    igraph_strvector_destroy(&enames);

    igraph_destroy(&g);

    return 0;
}
