/*
   IGraph library.
   Copyright (C) 2021  The igraph development team <igraph@igraph.org>

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
#include "test_utilities.inc"

int main() {
    igraph_vector_t v1, v2, v3, v4, v5, v6, v7, v8;

    igraph_vector_init_real(&v1, 3, 1e30, 2e30, 9e30);
    igraph_vector_init_real(&v2, 3, 1e30, 2e30, 3e30);
    igraph_vector_init_real(&v3, 2, 1e30, 2e30);
    igraph_vector_init_real(&v4, 0);
    igraph_vector_init_real(&v5, 3, 1e30, 2e30, 9e30);
    igraph_vector_init_real(&v6, 0);
    igraph_vector_init_real(&v7, 3, 9e30, 2e30, 1e30);
    igraph_vector_init_real(&v8, 2, 3e30, 3e30);

    igraph_vector_t *vectors[] = {&v1, &v2, &v3, &v4, &v5, &v6, &v7, &v8};
    long n = sizeof(vectors) / sizeof(igraph_vector_t *);

    printf("Lexicographical ordering:\n");
    igraph_qsort(vectors, n, sizeof(igraph_vector_t *), igraph_vector_lex_cmp);

    for (int i = 0; i < n; i++) {
        print_vector(vectors[i]);
    }

    printf("\nColexicographical ordering:\n");
    igraph_qsort(vectors, n, sizeof(igraph_vector_t *), igraph_vector_colex_cmp);

    for (int i = 0; i < n; i++) {
        print_vector(vectors[i]);
        igraph_vector_destroy(vectors[i]);
    }

    VERIFY_FINALLY_STACK();
    return 0;
}
