/*
   igraph library.
   Copyright (C) 2020-2024  The igraph development team <igraph@igraph.org>

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

#include <igraph.h>
#include <stdlib.h>

#include "test_utilities.h"

struct userdata {
    int i;
    igraph_vector_int_list_t *list;
};

int compare_vectors(const igraph_vector_int_t *v1, const igraph_vector_int_t *v2) {
    igraph_int_t s1, s2, i;

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


igraph_error_t handler(const igraph_vector_int_t *clique, void *arg) {
    struct userdata *ud;
    igraph_bool_t cont;

    ud = (struct userdata *) arg;
    cont = 1; /* true */

    if (compare_vectors(clique, igraph_vector_int_list_get_ptr(ud->list, ud->i)) != 0) {
        printf("igraph_maximal_cliques() and igraph_maximal_cliques_callback() give different results.\n");
        cont = 0; /* false */
    }

    ud->i += 1;

    return cont ? IGRAPH_SUCCESS : IGRAPH_STOP;
}


igraph_error_t handler_stop(const igraph_vector_int_t *clique, void *arg) {
    /* Stop search as soon as a 3-clique is found. */
    /* Since there are two 3-cliques in the test graph, this will stop the search before it is complete. */
    IGRAPH_UNUSED(arg);
    return (igraph_vector_int_size(clique) == 3) ? IGRAPH_STOP : IGRAPH_SUCCESS;
}


int main(void) {
    igraph_t graph;
    igraph_vector_int_list_t list;
    struct userdata ud;

    igraph_small(&graph, 6, 0,
                 1, 2, 2, 3, 3, 4, 4, 5, 5, 2, 2, 4,
                 -1);

    igraph_vector_int_list_init(&list, 0);
    igraph_maximal_cliques(&graph, &list, IGRAPH_UNLIMITED, IGRAPH_UNLIMITED, IGRAPH_UNLIMITED);

    ud.i = 0;
    ud.list = &list;

    /* Check that the callback function finds the same cliques as igraph_maximal_cliques() */
    IGRAPH_ASSERT(igraph_maximal_cliques_callback(&graph, 0, 0, &handler, (void *) &ud) == IGRAPH_SUCCESS);

    /* Check that the search can be stopped correctly */
    IGRAPH_ASSERT(igraph_maximal_cliques_callback(&graph, 0, 0, &handler_stop, NULL) == IGRAPH_SUCCESS);

    igraph_vector_int_list_destroy(&list);

    igraph_destroy(&graph);

    VERIFY_FINALLY_STACK();

    return 0;
}
