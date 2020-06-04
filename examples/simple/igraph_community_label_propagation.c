/* -*- mode: C -*-  */
/* vim:set ts=4 sw=4 sts=4 et: */
/*
   IGraph library.
   Copyright (C) 2007-2012  Gabor Csardi <csardi.gabor@gmail.com>
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
#include <stdio.h>

int main() {
    igraph_t g;
    igraph_vector_t membership;
    igraph_real_t modularity;

    igraph_famous(&g, "Zachary"); /* we use Zachary's karate club network */

    /* label propagation is a stochastic method; the result will depend on the random seed */
    igraph_rng_seed(igraph_rng_default(), 123);

    igraph_vector_init(&membership, 0);
    igraph_community_label_propagation(
                &g, &membership,
                /* weights= */ NULL, /* initial= */ NULL, /* fixed= */ NULL,
                &modularity);

    printf("%d communities found; modularity score is %f.\n",
           (int) (igraph_vector_max(&membership) + 1),
           modularity);

    igraph_vector_destroy(&membership);
    igraph_destroy(&g);

    return 0;
}

