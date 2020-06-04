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
   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
   02110-1301 USA

*/

#include <igraph.h>
#include <stdio.h>

void print_vector(const igraph_vector_t *v, FILE *f) {
    long int i;
    for (i = 0; i < igraph_vector_size(v); i++) {
        fprintf(f, " %li", (long int) VECTOR(*v)[i]);
    }
    fprintf(f, "\n");
}

void print_array(const igraph_array3_t *a, FILE *f) {
    long int i, j, k;
    for (k = 0; k < igraph_array3_n(a, 3); k++) {
        for (i = 0; i < igraph_array3_n(a, 1); i++) {
            for (j = 0; j < igraph_array3_n(a, 2); j++) {
                fprintf(f, " %li", (long int) ARRAY3(*a, i, j, k));
            }
            fprintf(f, "\n");
        }
        fprintf(f, "\n");
    }
}

int main() {
    igraph_array3_t a;
    long int i, j, k;
    long int s = 1;

    igraph_array3_init(&a, 5, 4, 3);
    igraph_array3_destroy(&a);

    igraph_array3_init(&a, 5, 4, 3);
    print_array(&a, stdout);
    if (igraph_array3_n(&a, 1) != 5) {
        return 1;
    }
    if (igraph_array3_n(&a, 2) != 4) {
        return 1;
    }
    if (igraph_array3_n(&a, 3) != 3) {
        return 1;
    }
    igraph_array3_destroy(&a);

    igraph_array3_init(&a, 5, 4, 3);
    for (k = 0; k < igraph_array3_n(&a, 3); k++) {
        for (j = 0; j < igraph_array3_n(&a, 2); j++) {
            for (i = 0; i < igraph_array3_n(&a, 1); i++) {
                ARRAY3(a, i, j, k) = s++;
            }
        }
    }
    print_array(&a, stdout);
    print_vector(&a.data, stdout);
    igraph_array3_destroy(&a);

    if (!IGRAPH_FINALLY_STACK_EMPTY) {
        return 2;
    }

    return 0;
}

