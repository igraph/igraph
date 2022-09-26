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

void sort_and_print_vector(igraph_vector_int_t *v) {
    igraph_integer_t i, n = igraph_vector_int_size(v);
    igraph_vector_int_sort(v);
    for (i = 0; i < n; i++) {
        printf(" %" IGRAPH_PRId, VECTOR(*v)[i]);
    }
    printf("\n");
}

int main(void) {

    igraph_t g;
    igraph_vector_int_list_t result;
    igraph_integer_t no;
    igraph_integer_t i;

    igraph_set_warning_handler(igraph_warning_handler_ignore);

    igraph_vector_int_list_init(&result, 0);
    igraph_small(&g, 7, 0, 0, 1, 1, 2, 2, 3, 3, 0, 2, 4, 4, 5, 2, 5, -1);

    igraph_biconnected_components(&g, &no, 0, 0, &result, 0);
    if (no != 2 || no != igraph_vector_int_list_size(&result)) {
        return 1;
    }
    for (i = 0; i < no; i++) {
        sort_and_print_vector(igraph_vector_int_list_get_ptr(&result, i));
    }
    igraph_vector_int_list_clear(&result);

    igraph_biconnected_components(&g, &no, 0, &result, 0, 0);
    if (no != 2 || no != igraph_vector_int_list_size(&result)) {
        return 2;
    }
    for (i = 0; i < no; i++) {
        sort_and_print_vector(igraph_vector_int_list_get_ptr(&result, i));
    }
    igraph_vector_int_list_clear(&result);

    igraph_biconnected_components(&g, &no, &result, 0, 0, 0);
    if (no != 2 || no != igraph_vector_int_list_size(&result)) {
        return 3;
    }
    for (i = 0; i < no; i++) {
        sort_and_print_vector(igraph_vector_int_list_get_ptr(&result, i));
    }

    igraph_vector_int_list_destroy(&result);
    igraph_destroy(&g);

    return 0;
}
