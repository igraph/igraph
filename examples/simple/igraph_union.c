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

#include <stdlib.h>
#include <stdio.h>

void print_vector(igraph_vector_t *v) {
    long int i, l = igraph_vector_size(v);
    for (i = 0; i < l; i++) {
        printf(" %li", (long int) VECTOR(*v)[i]);
    }
    printf("\n");
}

int print_free_vector_ptr(igraph_vector_ptr_t *v) {
    long int i, l = igraph_vector_ptr_size(v);
    printf("---\n");
    for (i = 0; i < l; i++) {
        print_vector(VECTOR(*v)[i]);
        igraph_vector_destroy(VECTOR(*v)[i]);
        IGRAPH_FREE(VECTOR(*v)[i]);
    }
    printf("===\n");
    return 0;
}

int main() {

    igraph_t left, right, uni;
    igraph_vector_t v;
    igraph_vector_ptr_t glist;
    igraph_vector_t edge_map1, edge_map2;
    igraph_vector_ptr_t edgemaps;
    long int i;

    igraph_vector_init(&edge_map1, 0);
    igraph_vector_init(&edge_map2, 0);

    igraph_vector_init_int_end(&v, -1, 0, 1, 1, 2, 2, 2, 2, 3, -1);
    igraph_create(&left, &v, 0, IGRAPH_DIRECTED);
    igraph_vector_destroy(&v);

    igraph_vector_init_int_end(&v, -1, 0, 1, 1, 2, 2, 2, 2, 4, -1);
    igraph_create(&right, &v, 0, IGRAPH_DIRECTED);
    igraph_vector_destroy(&v);

    igraph_union(&uni, &left, &right, &edge_map1, &edge_map2);
    igraph_write_graph_edgelist(&uni, stdout);
    igraph_vector_print(&edge_map1);
    igraph_vector_print(&edge_map2);

    igraph_destroy(&uni);
    igraph_destroy(&left);
    igraph_destroy(&right);
    igraph_vector_destroy(&edge_map1);
    igraph_vector_destroy(&edge_map2);

    /* Empty graph list */
    igraph_vector_ptr_init(&glist, 0);
    igraph_vector_ptr_init(&edgemaps, 0);
    igraph_union_many(&uni, &glist, &edgemaps);
    if (!igraph_is_directed(&uni) || igraph_vcount(&uni) != 0) {
        return 1;
    }
    print_free_vector_ptr(&edgemaps);
    igraph_vector_ptr_destroy(&glist);
    igraph_destroy(&uni);

    /* Non-empty graph list */
    igraph_vector_ptr_init(&glist, 10);
    for (i = 0; i < igraph_vector_ptr_size(&glist); i++) {
        VECTOR(glist)[i] = calloc(1, sizeof(igraph_t));
        igraph_vector_init_int_end(&v, -1, 0, 1, 1, 0, -1);
        igraph_create(VECTOR(glist)[i], &v, 0, IGRAPH_DIRECTED);
        igraph_vector_destroy(&v);
    }

    igraph_union_many(&uni, &glist, &edgemaps);
    igraph_write_graph_edgelist(&uni, stdout);

    for (i = 0; i < igraph_vector_ptr_size(&glist); i++) {
        igraph_destroy(VECTOR(glist)[i]);
        free(VECTOR(glist)[i]);
    }
    print_free_vector_ptr(&edgemaps);
    igraph_vector_ptr_destroy(&glist);
    igraph_destroy(&uni);

    /* Another non-empty graph list */
    igraph_vector_ptr_init(&glist, 10);
    for (i = 0; i < igraph_vector_ptr_size(&glist); i++) {
        VECTOR(glist)[i] = calloc(1, sizeof(igraph_t));
        igraph_vector_init_int_end(&v, -1, i, i + 1, 1, 0, -1);
        igraph_create(VECTOR(glist)[i], &v, 0, IGRAPH_DIRECTED);
        igraph_vector_destroy(&v);
    }

    igraph_union_many(&uni, &glist, &edgemaps);
    igraph_write_graph_edgelist(&uni, stdout);

    for (i = 0; i < igraph_vector_ptr_size(&glist); i++) {
        igraph_destroy(VECTOR(glist)[i]);
        free(VECTOR(glist)[i]);
    }
    print_free_vector_ptr(&edgemaps);
    igraph_vector_ptr_destroy(&glist);
    igraph_destroy(&uni);

    /* Undirected graph list*/
    igraph_vector_ptr_init(&glist, 10);
    for (i = 0; i < igraph_vector_ptr_size(&glist); i++) {
        VECTOR(glist)[i] = calloc(1, sizeof(igraph_t));
        igraph_vector_init_int_end(&v, -1, i, i + 1, 1, 0, -1);
        igraph_create(VECTOR(glist)[i], &v, 0, IGRAPH_UNDIRECTED);
        igraph_vector_destroy(&v);
    }

    igraph_union_many(&uni, &glist, &edgemaps);
    igraph_write_graph_edgelist(&uni, stdout);

    for (i = 0; i < igraph_vector_ptr_size(&glist); i++) {
        igraph_destroy(VECTOR(glist)[i]);
        free(VECTOR(glist)[i]);
    }
    print_free_vector_ptr(&edgemaps);
    igraph_vector_ptr_destroy(&glist);
    igraph_destroy(&uni);

    igraph_vector_ptr_destroy(&edgemaps);

    return 0;
}
