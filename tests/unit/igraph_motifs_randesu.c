/*
   igraph library.
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

#include "test_utilities.h"

igraph_error_t print_motif(const igraph_t *graph, const igraph_vector_int_t *vids,
                          igraph_int_t isoclass, void* extra) {
    IGRAPH_UNUSED(graph);
    IGRAPH_UNUSED(extra);
    printf("Class %" IGRAPH_PRId ": ", isoclass);
    print_vector_int(vids);
    return IGRAPH_SUCCESS;
}

int main(void) {

    igraph_t g;
    igraph_vector_t hist;
    igraph_int_t size;

    igraph_ring(&g, 1000, IGRAPH_DIRECTED, 1, 1);
    igraph_vector_init(&hist, 0);
    igraph_motifs_randesu(&g, &hist, 3, NULL);
    print_vector(&hist);
    igraph_destroy(&g);
    igraph_vector_destroy(&hist);

    igraph_famous(&g, "Octahedral");
    size = 3;
    printf("Motif size: %" IGRAPH_PRId "\n", size);
    igraph_motifs_randesu_callback(&g, size, NULL, &print_motif, NULL);
    size = 4;
    printf("Motif size: %" IGRAPH_PRId "\n", size);
    igraph_motifs_randesu_callback(&g, size, NULL, &print_motif, NULL);
    size = 5;
    printf("Motif size: %" IGRAPH_PRId "\n", size);
    igraph_motifs_randesu_callback(&g, size, NULL, &print_motif, NULL);
    size = 6;
    printf("Motif size: %" IGRAPH_PRId "\n", size);
    igraph_motifs_randesu_callback(&g, size, NULL, &print_motif, NULL);

    igraph_destroy(&g);

    VERIFY_FINALLY_STACK();

    return 0;
}
