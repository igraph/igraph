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

int compare_vectors(const igraph_vector_int_t *v1, const igraph_vector_int_t *v2) {
    igraph_integer_t s1, s2, i;

    s1 = igraph_vector_int_size(v1);
    s2 = igraph_vector_int_size(v2);
    if (s1 < s2) {
        return -1;
    }
    if (s1 > s2) {
        return 1;
    }
    for (i = 0; i < s1; ++i) {
        if (VECTOR(*v1)[i] < VECTOR(*v2)[i]) {
            return -1;
        }
        if (VECTOR(*v1)[i] > VECTOR(*v2)[i]) {
            return 1;
        }
    }
    return 0;
}

/* Takes a pointer vector of vectors. Sorts each vector, then sorts the pointer vector */
void canonicalize_list(igraph_vector_int_list_t *list) {
    igraph_integer_t len = igraph_vector_int_list_size(list);
    for (igraph_integer_t i = 0; i < len; ++i) {
        igraph_vector_int_sort(igraph_vector_int_list_get_ptr(list, i));
    }
    igraph_vector_int_list_sort(list, &compare_vectors);
}

struct userdata {
    igraph_integer_t i;
    igraph_vector_int_list_t *list;
};

igraph_error_t handler(const igraph_vector_int_t *clique, void *arg) {
    struct userdata *ud;
    igraph_bool_t cont;

    ud = (struct userdata *) arg;
    cont = 1; /* true */

    if (compare_vectors(clique, igraph_vector_int_list_get_ptr(ud->list, ud->i)) != 0) {
        printf("igraph_cliques() and igraph_cliques_callback() give different results.\n");
        cont = 0; /* false */
    }

    ud->i += 1;

    return cont ? IGRAPH_SUCCESS : IGRAPH_STOP;
}

void test_callback(const igraph_t *graph) {
    igraph_vector_int_list_t list;
    struct userdata ud;

    igraph_vector_int_list_init(&list, 0);
    igraph_cliques(graph, &list, 0, 0);

    ud.i = 0;
    ud.list = &list;

    igraph_cliques_callback(graph, 0, 0, &handler, (void *) &ud);

    igraph_vector_int_list_destroy(&list);
}


int main(void) {

    igraph_t g;
    igraph_vector_int_list_t result;
    igraph_es_t es;
    igraph_integer_t omega;
    igraph_integer_t i, j, n;
    const int params[] = {4, -1, 2, 2, 0, 0, -1, -1};

    igraph_set_warning_handler(igraph_warning_handler_ignore);

    igraph_vector_int_list_init(&result, 0);
    igraph_full(&g, 6, 0, 0);
    igraph_es_pairs_small(&es, 0, 0, 1, 0, 2, 3, 5, -1);
    igraph_delete_edges(&g, es);
    igraph_es_destroy(&es);

    for (j = 0; j < sizeof(params) / (2 * sizeof(params[0])); j++) {
        if (params[2 * j + 1] != 0) {
            igraph_cliques(&g, &result, params[2 * j], params[2 * j + 1]);
        } else {
            igraph_largest_cliques(&g, &result);
        }
        n = igraph_vector_int_list_size(&result);
        printf("%" IGRAPH_PRId " cliques found\n", n);
        canonicalize_list(&result);
        for (i = 0; i < n; i++) {
            igraph_vector_int_t* v = igraph_vector_int_list_get_ptr(&result, i);
            igraph_vector_int_print(v);
        }
    }

    igraph_clique_number(&g, &omega);
    printf("omega=%" IGRAPH_PRId "\n", omega);

    test_callback(&g);

    igraph_destroy(&g);

    igraph_kary_tree(&g, 5, 2, IGRAPH_TREE_OUT);
    igraph_cliques(&g, &result, 5, 5);
    if (igraph_vector_int_list_size(&result) != 0) {
        return 1;
    }

    igraph_destroy(&g);
    igraph_vector_int_list_destroy(&result);

    return 0;
}
