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

#include "test_utilities.h"

#include <string.h>
#include <stdlib.h>


static void check_vector_queries(const igraph_t *g) {
    igraph_vector_t vec;
    igraph_strvector_t svec;
    igraph_vector_bool_t bvec;
    igraph_strvector_t vnames, enames;
    igraph_vector_int_t vtypes, etypes;
    igraph_integer_t i, j;

    igraph_vector_int_init(&vtypes, 0);
    igraph_vector_int_init(&etypes, 0);
    igraph_strvector_init(&vnames, 0);
    igraph_strvector_init(&enames, 0);


    igraph_vector_init(&vec, 0);
    igraph_strvector_init(&svec, 0);
    igraph_vector_bool_init(&bvec, 0);

    igraph_cattribute_list(g, 0, 0, &vnames, &vtypes,
                           &enames, &etypes);
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
        } else if (VECTOR(vtypes)[j] == IGRAPH_ATTRIBUTE_STRING) {
            igraph_cattribute_VASV(g, STR(vnames, j), igraph_vss_all(), &svec);
            for (i = 0; i < igraph_vcount(g); i++) {
                const char *str = VAS(g, STR(vnames, j), i);
                if (strcmp(str, STR(svec, i))) {
                    exit(52);
                }
            }
        } else {
            igraph_cattribute_VABV(g, STR(vnames, j), igraph_vss_all(), &bvec);
            for (i = 0; i < igraph_vcount(g); i++) {
                igraph_bool_t b = VAB(g, STR(vnames, j), i);
                if (b != VECTOR(bvec)[i]) {
                    exit(53);
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
                    exit(54);
                }
            }
        } else if (VECTOR(etypes)[j] == IGRAPH_ATTRIBUTE_STRING) {
            igraph_cattribute_EASV(g, STR(enames, j),
                                   igraph_ess_all(IGRAPH_EDGEORDER_ID), &svec);
            for (i = 0; i < igraph_ecount(g); i++) {
                const char *str = EAS(g, STR(enames, j), i);
                if (strcmp(str, STR(svec, i))) {
                    exit(55);
                }
            }
        } else {
            igraph_cattribute_EABV(g, STR(enames, j),
                                   igraph_ess_all(IGRAPH_EDGEORDER_ID), &bvec);
            for (i = 0; i < igraph_ecount(g); i++) {
                igraph_bool_t b = EAB(g, STR(enames, j), i);
                if (b != VECTOR(bvec)[i]) {
                    exit(56);
                }
            }
        }
    }

    igraph_strvector_destroy(&svec);
    igraph_vector_destroy(&vec);
    igraph_vector_bool_destroy(&bvec);
    igraph_strvector_destroy(&enames);
    igraph_strvector_destroy(&vnames);
    igraph_vector_int_destroy(&etypes);
    igraph_vector_int_destroy(&vtypes);
}

int main(void) {

    igraph_t g, g2;
    FILE *ifile;
    igraph_integer_t i;
    igraph_vector_t y;
    igraph_strvector_t id;
    igraph_vector_bool_t type;
    char str[21];

    /* turn on attribute handling */
    igraph_set_attribute_table(&igraph_cattribute_table);

    ifile = fopen("links.net", "r");
    IGRAPH_ASSERT(ifile != NULL);

    igraph_read_graph_pajek(&g, ifile);
    fclose(ifile);

    SETGAB(&g, "bool_graph_attr", 1);
    SETEAB(&g, "bool_edge_attr", 0, 1);
    SETVAB(&g, "bool_vertex_attr", 0, 1);

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
    SETGAS(&g, "name", "toy");
    SETGAB(&g, "is_regular", 0);

    printf("Before deleting some attributes:\n");
    print_attributes(&g);
    /* Delete graph attributes */
    DELGA(&g, "id");
    DELGA(&g, "name");
    DELGA(&g, "is_regular");

    /* Delete vertex attributes */
    DELVA(&g, "x");
    DELVA(&g, "shape");
    DELVA(&g, "xfact");
    DELVA(&g, "yfact");

    /* Delete edge attributes */
    DELEA(&g, "hook1");
    DELEA(&g, "hook2");
    DELEA(&g, "label");

    printf("After deleting some attributes:\n");
    print_attributes(&g);

    /* Set vertex attributes */
    SETVAN(&g, "y", 0, -1);
    SETVAN(&g, "y", 1, 2.1);

    SETVAS(&g, "id", 0, "foo");
    SETVAS(&g, "id", 1, "bar");

    SETVAB(&g, "type", 0, 1);
    SETVAB(&g, "type", 1, 0);

    /* Set edge attributes */
    SETEAN(&g, "weight", 2, 100.0);
    SETEAN(&g, "weight", 0, -100.1);

    SETEAS(&g, "color", 2, "RED");
    SETEAS(&g, "color", 0, "Blue");

    SETEAB(&g, "type", 0, 1);
    SETEAB(&g, "type", 2, 0);

    printf("After setting vertex and edge attributes:\n");
    print_attributes(&g);
    check_vector_queries(&g);

    /* Set vertex attributes as vector */
    igraph_vector_init(&y, igraph_vcount(&g));
    igraph_vector_fill(&y, 1.23);
    SETVANV(&g, "y", &y);
    igraph_vector_destroy(&y);

    igraph_vector_init_range(&y, 0, igraph_vcount(&g));
    SETVANV(&g, "foobar", &y);
    igraph_vector_destroy(&y);

    igraph_vector_bool_init(&type, igraph_vcount(&g));
    for (i = 0; i < igraph_vcount(&g); i++) {
        VECTOR(type)[i] = (i % 2 == 1);
    }
    SETVABV(&g, "type", &type);
    igraph_vector_bool_destroy(&type);

    igraph_strvector_init(&id, igraph_vcount(&g));
    for (i = 0; i < igraph_vcount(&g); i++) {
        snprintf(str, sizeof(str) - 1, "%" IGRAPH_PRId, i);
        igraph_strvector_set(&id, i, str);
    }
    SETVASV(&g, "foo", &id);
    igraph_strvector_destroy(&id);

    igraph_strvector_init(&id, igraph_vcount(&g));
    for (i = 0; i < igraph_vcount(&g); i++) {
        snprintf(str, sizeof(str) - 1, "%" IGRAPH_PRId, i);
        igraph_strvector_set(&id, i, str);
    }
    SETVASV(&g, "id", &id);
    igraph_strvector_destroy(&id);

    /* Set edge attributes as vector */
    igraph_vector_init(&y, igraph_ecount(&g));
    igraph_vector_fill(&y, 12.3);
    SETEANV(&g, "weight", &y);
    igraph_vector_destroy(&y);

    igraph_vector_init_range(&y, 0, igraph_ecount(&g));
    SETEANV(&g, "foobar", &y);
    igraph_vector_destroy(&y);

    igraph_vector_bool_init(&type, igraph_ecount(&g));
    for (i = 0; i < igraph_ecount(&g); i++) {
        VECTOR(type)[i] = (i % 2 == 1);
    }
    SETEABV(&g, "type", &type);
    igraph_vector_bool_destroy(&type);

    igraph_strvector_init(&id, igraph_ecount(&g));
    for (i = 0; i < igraph_ecount(&g); i++) {
        snprintf(str, sizeof(str) - 1, "%" IGRAPH_PRId, i);
        igraph_strvector_set(&id, i, str);
    }
    SETEASV(&g, "foo", &id);
    igraph_strvector_destroy(&id);

    igraph_strvector_init(&id, igraph_ecount(&g));
    for (i = 0; i < igraph_ecount(&g); i++) {
        snprintf(str, sizeof(str) - 1, "%" IGRAPH_PRId, i);
        igraph_strvector_set(&id, i, str);
    }
    SETEASV(&g, "color", &id);
    igraph_strvector_destroy(&id);

    printf("After setting vertex and edge attributes by vector:\n");
    print_attributes(&g);
    /* Delete all remaining attributes */
    DELALL(&g);

    printf("After deleting all attributes:\n");
    print_attributes(&g);

    igraph_destroy(&g);

    printf("Setting attributes on vertex 0, permuting with 2\n");
    {
        igraph_t g2;
        igraph_vector_int_t per;

        igraph_vector_int_init_int(&per, 3, 2, 1, 0);

        igraph_small(&g, 3, IGRAPH_DIRECTED, 0,1, 1,2, -1);
        SETVAN(&g, "foo", 0, 5);
        SETVAB(&g, "bar", 0, 1);
        SETVAS(&g, "baz", 0, "foobar");
        printf("Permuting to different graph:\n");
        igraph_permute_vertices(&g, &g2, &per);
        print_attributes(&g2);
        printf("Permuting to same graph:\n");
        igraph_cattribute_table.permute_vertices(&g, &g, &per);
        print_attributes(&g);
        igraph_destroy(&g);
        igraph_destroy(&g2);
        igraph_vector_int_destroy(&per);
    }


    VERIFY_FINALLY_STACK();
    return 0;
}
