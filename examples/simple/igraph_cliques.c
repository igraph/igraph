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

int compare_vectors(const void *p1, const void *p2) {
    igraph_vector_t *v1, *v2;
    long s1, s2, i;

    v1 = *((igraph_vector_t **) p1);
    v2 = *((igraph_vector_t **) p2);
    s1 = igraph_vector_size(v1);
    s2 = igraph_vector_size(v2);
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
void canonicalize_list(igraph_vector_ptr_t *list) {
    long i, len;
    len = igraph_vector_ptr_size(list);
    for (i = 0; i < len; ++i) {
        igraph_vector_sort((igraph_vector_t *) VECTOR(*list)[i]);
    }
    qsort(&(VECTOR(*list)[0]), len, sizeof(void *), &compare_vectors);
}

void print_vector(igraph_vector_t *v) {
    long int i, n = igraph_vector_size(v);
    for (i = 0; i < n; i++) {
        printf(" %li", (long int) VECTOR(*v)[i]);
    }
    printf("\n");
}

void warning_handler_ignore(const char* reason, const char* file, int line, int e) {
}


struct userdata {
    int i;
    igraph_vector_ptr_t *list;
};

igraph_bool_t handler(igraph_vector_t *clique, void *arg) {
    struct userdata *ud;
    igraph_bool_t cont;

    ud = (struct userdata *) arg;
    cont = 1; /* true */

    if (compare_vectors(&clique, &(VECTOR(*(ud->list))[ud->i])) != 0) {
        printf("igraph_cliques() and igraph_cliques_callback() give different results.\n");
        cont = 0; /* false */
    }

    igraph_vector_destroy(clique);
    igraph_free(clique);

    ud->i += 1;

    return cont;
}

void test_callback(const igraph_t *graph) {
    igraph_vector_ptr_t list;
    struct userdata ud;

    igraph_vector_ptr_init(&list, 0);
    igraph_cliques(graph, &list, 0, 0);

    ud.i = 0;
    ud.list = &list;

    igraph_cliques_callback(graph, 0, 0, &handler, (void *) &ud);

    IGRAPH_VECTOR_PTR_SET_ITEM_DESTRUCTOR(&list, igraph_vector_destroy);
    igraph_vector_ptr_destroy_all(&list);
}


int main() {

    igraph_t g;
    igraph_vector_ptr_t result;
    igraph_es_t es;
    igraph_integer_t omega;
    long int i, j, n;
    const int params[] = {4, -1, 2, 2, 0, 0, -1, -1};

    igraph_set_warning_handler(warning_handler_ignore);

    igraph_vector_ptr_init(&result, 0);
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
        n = igraph_vector_ptr_size(&result);
        printf("%ld cliques found\n", (long)n);
        canonicalize_list(&result);
        for (i = 0; i < n; i++) {
            igraph_vector_t* v = (igraph_vector_t*) igraph_vector_ptr_e(&result, i);
            print_vector(v);
            igraph_vector_destroy(v);
            igraph_free(v);
        }
    }

    igraph_clique_number(&g, &omega);
    printf("omega=%ld\n", (long)omega);

    test_callback(&g);

    igraph_destroy(&g);

    igraph_tree(&g, 5, 2, IGRAPH_TREE_OUT);
    igraph_cliques(&g, &result, 5, 5);
    if (igraph_vector_ptr_size(&result) != 0) {
        return 1;
    }

    igraph_destroy(&g);
    igraph_vector_ptr_destroy(&result);

    return 0;
}
