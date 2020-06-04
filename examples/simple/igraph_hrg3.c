/* -*- mode: C++ -*-  */
/*
   IGraph library.
   Copyright (C) 2011-2012  Gabor Csardi <csardi.gabor@gmail.com>
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
   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
   02110-1301 USA

*/

#include <igraph.h>
#include <stdio.h>

int main() {
    igraph_t karate;
    igraph_vector_t edges, prob;
    long int i, n;

    igraph_rng_seed(igraph_rng_default(), 42);

    igraph_small(&karate, 34, IGRAPH_UNDIRECTED,
                 0, 1, 0, 2, 0, 3, 0, 4, 0, 5, 0, 6, 0, 7, 0, 8, 0, 10, 0, 11, 0, 12, 0, 13,
                 0, 17, 0, 19, 0, 21, 0, 31,
                 1, 2, 1, 3, 1, 7, 1, 13, 1, 17, 1, 19, 1, 21, 1, 30,
                 2, 3, 2, 7, 2, 27, 2, 28, 2, 32, 2, 9, 2, 8, 2, 13,
                 3, 7, 3, 12, 3, 13,
                 4, 6, 4, 10,
                 5, 6, 5, 10, 5, 16,
                 6, 16,
                 8, 30, 8, 32, 8, 33,
                 9, 33,
                 13, 33,
                 14, 32, 14, 33,
                 15, 32, 15, 33,
                 18, 32, 18, 33,
                 19, 33,
                 20, 32, 20, 33,
                 22, 32, 22, 33,
                 23, 25, 23, 27, 23, 32, 23, 33, 23, 29,
                 24, 25, 24, 27, 24, 31,
                 25, 31,
                 26, 29, 26, 33,
                 27, 33,
                 28, 31, 28, 33,
                 29, 32, 29, 33,
                 30, 32, 30, 33,
                 31, 32, 31, 33,
                 32, 33,
                 -1);

    igraph_vector_init(&edges, 0);
    igraph_vector_init(&prob, 0);
    igraph_hrg_predict(&karate, &edges, &prob, /* hrg= */ 0, /* start= */ 0,
                       /* num_samples= */ 100, /* num_bins= */ 25);

    /* We do some simple validity tests on the resolts only; the exact results
     * are different on i386 vs other platforms due to numerical inaccuracies */
    n = igraph_vector_size(&edges);
    for (i = 0; i < n; i++) {
		if (VECTOR(edges)[i] < 0 || VECTOR(edges)[i] >= igraph_vcount(&karate)) {
            printf("Invalid edges vector:\n");
            igraph_vector_print(&edges);
			return 1;
		}
	}
    n = igraph_vector_size(&prob);
    for (i = 0; i < n; i++) {
		if (VECTOR(prob)[i] < 0 || VECTOR(prob)[i] > 1) {
            printf("Invalid prob vector:\n");
            igraph_vector_print(&prob);
			return 2;
		}
	}

    igraph_destroy(&karate);

    return 0;
}
