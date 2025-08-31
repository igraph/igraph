/*
   igraph library.
   Copyright (C) 2010-2012  Gabor Csardi <csardi.gabor@gmail.com>
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

int main(void) {
    igraph_t graph;
    igraph_vector_t knn, knnk;
    igraph_vector_t weights;

    /* Initialize the library. */
    igraph_setup();

    igraph_famous(&graph, "Zachary");

    igraph_vector_init(&knn, 0);
    igraph_vector_init(&knnk, 0);

    igraph_avg_nearest_neighbor_degree(&graph, igraph_vss_all(),
                                       IGRAPH_ALL, IGRAPH_ALL,
                                       &knn, &knnk, /*weights=*/ NULL);

    printf("knn: ");
    igraph_vector_print(&knn);
    printf("knn(k): ");
    igraph_vector_print(&knnk);

    igraph_vector_init_range(&weights, 0, igraph_ecount(&graph));

    igraph_avg_nearest_neighbor_degree(&graph, igraph_vss_all(),
                                       IGRAPH_ALL, IGRAPH_ALL,
                                       &knn, &knnk, &weights);
    igraph_vector_destroy(&weights);

    printf("knn: ");
    igraph_vector_print(&knn);
    printf("knn(k): ");
    igraph_vector_print(&knnk);

    igraph_vector_destroy(&knn);
    igraph_vector_destroy(&knnk);

    igraph_destroy(&graph);

    return 0;
}
