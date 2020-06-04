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

void igraph_vector_print(const igraph_vector_t *v) {
    long int i;
    for (i = 0; i < igraph_vector_size(v); i++) {
        printf("%li ", (long int)VECTOR(*v)[i]);
    }
    printf("\n");
}

int main() {

    igraph_t g;
    const igraph_vector_t v = IGRAPH_VECTOR_NULL;
    igraph_real_t edges1[] = { 0, 1, 1, 2, 2, 2, 2, 3, 2, 4, 3, 4 };
    igraph_vector_t was;
    igraph_integer_t size;
    igraph_es_t it;
    long int i;

    igraph_vector_view(&v, edges1, sizeof(edges1) / sizeof(igraph_real_t));
    igraph_vector_init(&was, 0);

    /******************************************/
    /* Directed graph                         */
    /******************************************/

    igraph_create(&g, &v, 0, IGRAPH_DIRECTED);

    /* Simple test, all neighbors */
    for (i = 0; i <= igraph_vector_max(&v); i++) {
        igraph_vector_clear(&was);
        igraph_es_adj(&g, &it, i, IGRAPH_ALL);
        igraph_es_size(&g, &it, &size);
        printf("%ld\n", (long)size);
        while (!igraph_es_end(&g, &it)) {
            igraph_vector_push_back(&was, igraph_es_adj_vertex(&g, &it));
            igraph_es_next(&g, &it);
        }
        igraph_es_destroy(&it);
        igraph_vector_sort(&was);
        igraph_vector_print(&was);
    }

    /* Simple test, outgoing neighbors */
    for (i = 0; i <= igraph_vector_max(&v); i++) {
        igraph_vector_clear(&was);
        igraph_es_adj(&g, &it, i, IGRAPH_OUT);
        igraph_es_size(&g, &it, &size);
        printf("%ld\n", (long)size);
        while (!igraph_es_end(&g, &it)) {
            igraph_vector_push_back(&was, igraph_es_adj_vertex(&g, &it));
            igraph_es_next(&g, &it);
        }
        igraph_es_destroy(&it);
        igraph_vector_sort(&was);
        igraph_vector_print(&was);
    }


    /* Simple test, incoming neighbors */
    for (i = 0; i <= igraph_vector_max(&v); i++) {
        igraph_vector_clear(&was);
        igraph_es_adj(&g, &it, i, IGRAPH_IN);
        igraph_es_size(&g, &it, &size);
        printf("%ld\n", (long)size);
        while (!igraph_es_end(&g, &it)) {
            igraph_vector_push_back(&was, igraph_es_adj_vertex(&g, &it));
            igraph_es_next(&g, &it);
        }
        igraph_es_destroy(&it);
        igraph_vector_sort(&was);
        igraph_vector_print(&was);
    }

    igraph_destroy(&g);

    /******************************************/
    /* Undirected graph                       */
    /******************************************/

    igraph_create(&g, &v, 0, IGRAPH_UNDIRECTED);

    /* Simple test, all neighbors */
    for (i = 0; i <= igraph_vector_max(&v); i++) {
        igraph_vector_clear(&was);
        igraph_es_adj(&g, &it, i, IGRAPH_ALL);
        igraph_es_size(&g, &it, &size);
        printf("%ld\n", (long)size);
        while (!igraph_es_end(&g, &it)) {
            igraph_vector_push_back(&was, igraph_es_adj_vertex(&g, &it));
            igraph_es_next(&g, &it);
        }
        igraph_es_destroy(&it);
        igraph_vector_sort(&was);
        igraph_vector_print(&was);
    }

    igraph_destroy(&g);
    igraph_vector_destroy(&was);

    return 0;
}
