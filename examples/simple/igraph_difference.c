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

    igraph_t orig, sub, diff;
    igraph_vector_t v;

    /* Subtract from itself */
    printf("subtract itself\n");
    igraph_vector_init_int_end(&v, -1, 0, 1, 1, 2, 2, 1, 4, 5, -1);
    igraph_create(&orig, &v, 0, IGRAPH_DIRECTED);
    igraph_vector_destroy(&v);

    igraph_difference(&diff, &orig, &orig);
    igraph_write_graph_edgelist(&diff, stdout);
    if (igraph_ecount(&diff) != 0 ||
        igraph_vcount(&diff) != igraph_vcount(&orig)) {
        return 1;
    }

    igraph_destroy(&orig);
    igraph_destroy(&diff);

    /* Same for undirected graph */
    printf("subtract itself, undirected\n");
    igraph_vector_init_int_end(&v, -1, 0, 1, 1, 2, 2, 1, 4, 5, -1);
    igraph_create(&orig, &v, 0, IGRAPH_UNDIRECTED);
    igraph_vector_destroy(&v);

    igraph_vector_init_int_end(&v, -1, 1, 0, 1, 2, 2, 1, 4, 5, -1);
    igraph_create(&sub, &v, 0, IGRAPH_UNDIRECTED);
    igraph_vector_destroy(&v);

    igraph_difference(&diff, &orig, &sub);
    igraph_write_graph_edgelist(&diff, stdout);
    if (igraph_ecount(&diff) != 0 ||
        igraph_vcount(&diff) != igraph_vcount(&orig)) {
        return 2;
    }

    igraph_destroy(&orig);
    igraph_destroy(&sub);
    igraph_destroy(&diff);

    /* Subtract the empty graph */
    printf("subtract empty\n");
    igraph_vector_init_int_end(&v, -1, 0, 1, 1, 2, 2, 1, 4, 5, -1);
    igraph_create(&orig, &v, 0, IGRAPH_DIRECTED);
    igraph_vector_destroy(&v);

    igraph_empty(&sub, 3, IGRAPH_DIRECTED);
    igraph_difference(&diff, &orig, &sub);
    igraph_write_graph_edgelist(&diff, stdout);
    if (igraph_ecount(&diff) != igraph_ecount(&orig) ||
        igraph_vcount(&diff) != igraph_vcount(&orig)) {
        return 3;
    }

    igraph_destroy(&orig);
    igraph_destroy(&sub);
    igraph_destroy(&diff);

    /* A `real' example */
    printf("real example\n");
    igraph_vector_init_int_end(&v, -1, 0, 1, 1, 2, 2, 1, 4, 5, 8, 9, -1);
    igraph_create(&orig, &v, 0, IGRAPH_DIRECTED);
    igraph_vector_destroy(&v);

    igraph_vector_init_int_end(&v, -1, 0, 1, 5, 4, 2, 1, 6, 7, -1);
    igraph_create(&sub, &v, 0, IGRAPH_DIRECTED);
    igraph_vector_destroy(&v);

    igraph_difference(&diff, &orig, &sub);
    igraph_write_graph_edgelist(&diff, stdout);

    igraph_destroy(&diff);
    igraph_destroy(&orig);
    igraph_destroy(&sub);

    /* undirected version */
    printf("real example, undirected\n");
    igraph_vector_init_int_end(&v, -1, 0, 1, 1, 2, 2, 1, 4, 5, 8, 9, 8, 10, 8, 13, 8, 11, 8, 12, -1);
    igraph_create(&orig, &v, 0, IGRAPH_UNDIRECTED);
    igraph_vector_destroy(&v);

    igraph_vector_init_int_end(&v, -1, 0, 1, 5, 4, 2, 1, 6, 7, 8, 10, 8, 13, -1);
    igraph_create(&sub, &v, 0, IGRAPH_UNDIRECTED);
    igraph_vector_destroy(&v);

    igraph_difference(&diff, &orig, &sub);
    igraph_write_graph_edgelist(&diff, stdout);

    igraph_destroy(&diff);
    igraph_destroy(&orig);
    igraph_destroy(&sub);

    /* undirected version with loop edge, tests Github issue #597 */
    printf("Github issue #597, undirected\n");
    igraph_vector_init_int_end(&v, -1, 0, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 0, -1);
    igraph_create(&orig, &v, 0, IGRAPH_UNDIRECTED);
    igraph_vector_destroy(&v);

    igraph_vector_init_int_end(&v, -1, 0, 1, 1, 2, 2, 3, 3, 4, 4, 0, -1);
    igraph_create(&sub, &v, 0, IGRAPH_UNDIRECTED);
    igraph_vector_destroy(&v);

    igraph_difference(&diff, &orig, &sub);
    igraph_write_graph_edgelist(&diff, stdout);

    igraph_destroy(&diff);
    igraph_destroy(&orig);
    igraph_destroy(&sub);

    return 0;
}
