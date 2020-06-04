/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2006-2012  Gabor Csardi <csardi.gabor@gmail.com>
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

#include <stdlib.h>

void print_vector(igraph_vector_t *v) {
    long int i, l = igraph_vector_size(v);
    for (i = 0; i < l; i++) {
        printf(" %li", (long int) VECTOR(*v)[i]);
    }
    printf("\n");
}

int main() {

    igraph_t left, right, uni;
    igraph_vector_t v;
    igraph_vector_ptr_t glist;
    long int i;

    igraph_vector_init_int_end(&v, -1, 0, 1, 1, 2, 2, 2, 2, 3, -1);
    igraph_create(&left, &v, 0, IGRAPH_DIRECTED);
    igraph_vector_destroy(&v);

    igraph_vector_init_int_end(&v, -1, 0, 1, 1, 2, 2, 2, 2, 4, -1);
    igraph_create(&right, &v, 0, IGRAPH_DIRECTED);
    igraph_vector_destroy(&v);

    igraph_disjoint_union(&uni, &left, &right);
    igraph_vector_init(&v, 0);
    igraph_get_edgelist(&uni, &v, 0);
    igraph_vector_sort(&v);
    print_vector(&v);

    igraph_vector_destroy(&v);
    igraph_destroy(&left);
    igraph_destroy(&right);
    igraph_destroy(&uni);

    /* Empty graph list */
    igraph_vector_ptr_init(&glist, 0);
    igraph_disjoint_union_many(&uni, &glist);
    if (!igraph_is_directed(&uni) || igraph_vcount(&uni) != 0) {
        return 1;
    }
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

    igraph_disjoint_union_many(&uni, &glist);
    igraph_vector_init(&v, 0);
    igraph_get_edgelist(&uni, &v, 0);
    igraph_vector_sort(&v);
    print_vector(&v);
    igraph_vector_destroy(&v);

    for (i = 0; i < igraph_vector_ptr_size(&glist); i++) {
        igraph_destroy(VECTOR(glist)[i]);
        free(VECTOR(glist)[i]);
    }
    igraph_vector_ptr_destroy(&glist);
    igraph_destroy(&uni);

    return 0;
}
