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

int main() {

    igraph_t g;
    igraph_vector_t v;
    int ret;
    igraph_es_t es;

    igraph_vector_init(&v, 8);
    VECTOR(v)[0] = 0;
    VECTOR(v)[1] = 1;
    VECTOR(v)[2] = 1;
    VECTOR(v)[3] = 2;
    VECTOR(v)[4] = 2;
    VECTOR(v)[5] = 3;
    VECTOR(v)[6] = 2;
    VECTOR(v)[7] = 2;
    igraph_create(&g, &v, 0, 0);

    igraph_es_pairs_small(&es, IGRAPH_DIRECTED, 3, 2, -1);
    igraph_delete_edges(&g, es);
    if (igraph_ecount(&g) != 3) {
        return 1;
    }

    /* error test, no such edge to delete */
    igraph_set_error_handler(igraph_error_handler_ignore);
    ret = igraph_delete_edges(&g, es);
    if (ret != IGRAPH_EINVAL) {
        printf("Error code: %i\n", ret);
        return 2;
    }
    if (igraph_ecount(&g) != 3) {
        return 3;
    }

    /* error test, invalid vertex id */
    igraph_es_destroy(&es);
    igraph_es_pairs_small(&es, IGRAPH_DIRECTED, 10, 2, -1);
    ret = igraph_delete_edges(&g, es);
    if (ret != IGRAPH_EINVVID) {
        return 4;
    }
    if (igraph_ecount(&g) != 3) {
        return 5;
    }

    /* error test, invalid (odd) length */
    igraph_es_destroy(&es);
    igraph_es_pairs_small(&es, IGRAPH_DIRECTED, 0, 1, 2, -1);
    ret = igraph_delete_edges(&g, es);
    if (ret != IGRAPH_EINVAL) {
        return 6;
    }
    if (igraph_ecount(&g) != 3) {
        return 7;
    }

    igraph_es_destroy(&es);
    igraph_vector_destroy(&v);
    igraph_destroy(&g);

    return 0;
}
