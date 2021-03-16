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
    long int i, n = igraph_vector_size(v);
    igraph_real_t sum = 0.0;
    for (i = 0; i < n; i++) {
        if (!igraph_is_nan(VECTOR(*v)[i])) {
            sum += VECTOR(*v)[i];
        }
    }
    for (i = 0; i < n; i++) {
        igraph_real_printf(VECTOR(*v)[i] / sum);
        printf(" ");
    }
    printf("\n");
}

igraph_bool_t print_motif(const igraph_t *graph, igraph_vector_t *vids,
                          int isoclass, void* extra) {
    printf("Class %d: ", isoclass);
    igraph_vector_print(vids);
    return 0;
}


int main() {

    igraph_t g;
    igraph_vector_t hist;
    igraph_vector_t cp_3, cp_4;

    igraph_vector_init_real(&cp_3, 3, 0.0, 0.0, 0.0);
    igraph_vector_init_real(&cp_4, 4, 0.0, 0.0, 0.0, 0.0);

    igraph_ring(&g, 1000, IGRAPH_DIRECTED, 1, 1);
    igraph_vector_init(&hist, 0);
    igraph_motifs_randesu(&g, &hist, 3, &cp_3);
    print_vector(&hist);
    igraph_destroy(&g);
    igraph_vector_destroy(&hist);

    igraph_famous(&g, "bull");
    igraph_motifs_randesu_callback(&g, 3, &cp_3, &print_motif, 0);
    igraph_motifs_randesu_callback(&g, 4, &cp_4, &print_motif, 0);
    igraph_destroy(&g);

    igraph_vector_destroy(&cp_3);
    igraph_vector_destroy(&cp_4);
    return 0;
}
