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

void print_vector(igraph_vector_t *v) {
    long int i, l = igraph_vector_size(v);
    for (i = 0; i < l; i++) {
        printf(" %li", (long int) VECTOR(*v)[i]);
    }
    printf("\n");
}

int main() {

    igraph_t left, right, isec;
    igraph_vector_t v;
    igraph_vector_ptr_t glist;
    igraph_t g1, g2, g3;
    igraph_vector_t edge_map1, edge_map2;

    igraph_vector_init_int_end(&v, -1, 0, 1, 1, 2, 2, 3, -1);
    igraph_create(&left, &v, 0, IGRAPH_DIRECTED);
    igraph_vector_destroy(&v);

    igraph_vector_init_int_end(&v, -1, 1, 0, 5, 4, 1, 2, 3, 2, -1);
    igraph_create(&right, &v, 0, IGRAPH_DIRECTED);
    igraph_vector_destroy(&v);

    igraph_vector_init(&edge_map1, 0);
    igraph_vector_init(&edge_map2, 0);

    igraph_intersection(&isec, &left, &right, &edge_map1, &edge_map2);
    igraph_vector_init(&v, 0);
    igraph_get_edgelist(&isec, &v, 0);
    printf("---\n");
    print_vector(&v);
    print_vector(&edge_map1);
    print_vector(&edge_map2);
    printf("---\n");
    igraph_vector_destroy(&v);
    igraph_destroy(&left);
    igraph_destroy(&right);
    igraph_destroy(&isec);
    igraph_vector_destroy(&edge_map1);
    igraph_vector_destroy(&edge_map2);

    /* empty graph list */
    igraph_vector_ptr_init(&glist, 0);
    igraph_intersection_many(&isec, &glist, 0);
    if (igraph_vcount(&isec) != 0 || !igraph_is_directed(&isec)) {
        return 1;
    }
    igraph_destroy(&isec);
    igraph_vector_ptr_destroy(&glist);

    /* graph list with an empty graph */
    igraph_vector_ptr_init(&glist, 3);
    igraph_vector_init_int_end(&v, -1, 0, 1, 1, 2, 2, 3, -1);
    igraph_create(&g1, &v, 0, IGRAPH_DIRECTED);
    igraph_vector_destroy(&v);
    igraph_vector_init_int_end(&v, -1, 0, 1, 1, 2, 2, 3, -1);
    igraph_create(&g2, &v, 0, IGRAPH_DIRECTED);
    igraph_vector_destroy(&v);
    igraph_empty(&g3, 10, IGRAPH_DIRECTED);

    VECTOR(glist)[0] = &g1;
    VECTOR(glist)[1] = &g2;
    VECTOR(glist)[2] = &g3;
    igraph_intersection_many(&isec, &glist, 0);
    if (igraph_ecount(&isec) != 0 || igraph_vcount(&isec) != 10) {
        return 2;
    }
    igraph_destroy(&g1);
    igraph_destroy(&g2);
    igraph_destroy(&g3);
    igraph_destroy(&isec);
    igraph_vector_ptr_destroy(&glist);

    /* "proper" graph list */
    igraph_vector_ptr_init(&glist, 3);
    igraph_vector_init_int_end(&v, -1, 0, 1, 1, 2, 2, 3, -1);
    igraph_create(&g1, &v, 0, IGRAPH_DIRECTED);
    igraph_vector_destroy(&v);
    igraph_vector_init_int_end(&v, -1, 0, 1, 1, 2, 2, 3, 3, 2, 4, 5, 6, 5, -1);
    igraph_create(&g2, &v, 0, IGRAPH_DIRECTED);
    igraph_vector_destroy(&v);
    igraph_vector_init_int_end(&v, -1, 2, 3, 1, 0, 1, 2, 3, 2, 4, 5, 6, 5, 2, 3, -1);
    igraph_create(&g3, &v, 0, IGRAPH_DIRECTED);
    igraph_vector_destroy(&v);

    VECTOR(glist)[0] = &g1;
    VECTOR(glist)[1] = &g2;
    VECTOR(glist)[2] = &g3;
    igraph_intersection_many(&isec, &glist, 0);
    igraph_write_graph_edgelist(&isec, stdout);
    igraph_destroy(&g1);
    igraph_destroy(&g2);
    igraph_destroy(&g3);
    igraph_destroy(&isec);
    igraph_vector_ptr_destroy(&glist);

    return 0;
}
