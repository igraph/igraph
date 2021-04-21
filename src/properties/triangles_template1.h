/* -*- mode: C -*-  */
/* vim:set ts=4 sw=4 sts=4 et: */
/*
   IGraph library.
   Copyright (C) 2005-2012  Gabor Csardi <csardi.gabor@gmail.com>
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

long int no_of_nodes = igraph_vcount(graph);
igraph_vit_t vit;
long int nodes_to_calc;
igraph_vector_int_t *neis1, *neis2;
igraph_real_t triangles;
long int i, j, k;
long int neilen1, neilen2;
long int *neis;
igraph_lazy_adjlist_t adjlist;

IGRAPH_CHECK(igraph_vit_create(graph, vids, &vit));
IGRAPH_FINALLY(igraph_vit_destroy, &vit);
nodes_to_calc = IGRAPH_VIT_SIZE(vit);

if (nodes_to_calc == 0) {
    igraph_vector_clear(res);
    igraph_vit_destroy(&vit);
    IGRAPH_FINALLY_CLEAN(1);
    return IGRAPH_SUCCESS;
}

neis = IGRAPH_CALLOC(no_of_nodes, long int);
if (neis == 0) {
    IGRAPH_ERROR("local undirected transitivity failed", IGRAPH_ENOMEM);
}
IGRAPH_FINALLY(igraph_free, neis);

IGRAPH_CHECK(igraph_vector_resize(res, nodes_to_calc));

IGRAPH_CHECK(igraph_lazy_adjlist_init(graph, &adjlist, IGRAPH_ALL, IGRAPH_NO_LOOPS, IGRAPH_NO_MULTIPLE));
IGRAPH_FINALLY(igraph_lazy_adjlist_destroy, &adjlist);

for (i = 0; !IGRAPH_VIT_END(vit); IGRAPH_VIT_NEXT(vit), i++) {
    long int node = IGRAPH_VIT_GET(vit);

    IGRAPH_ALLOW_INTERRUPTION();

    neis1 = igraph_lazy_adjlist_get(&adjlist, (igraph_integer_t) node);
    neilen1 = igraph_vector_int_size(neis1);
    for (j = 0; j < neilen1; j++) {
        neis[ (long int)VECTOR(*neis1)[j] ] = i + 1;
    }
    triangles = 0;

    for (j = 0; j < neilen1; j++) {
        long int v = (long int) VECTOR(*neis1)[j];
        neis2 = igraph_lazy_adjlist_get(&adjlist, (igraph_integer_t) v);
        neilen2 = igraph_vector_int_size(neis2);
        for (k = 0; k < neilen2; k++) {
            long int v2 = (long int) VECTOR(*neis2)[k];
            if (neis[v2] == i + 1) {
                triangles += 1.0;
            }
        }
    }

#ifdef TRANSIT
    if (mode == IGRAPH_TRANSITIVITY_ZERO && neilen1 < 2) {
        VECTOR(*res)[i] = 0.0;
    } else {
        VECTOR(*res)[i] = triangles / neilen1 / (neilen1 - 1);
    }
#else
    VECTOR(*res)[i] = triangles / 2;
#endif
}

igraph_lazy_adjlist_destroy(&adjlist);
IGRAPH_FREE(neis);
igraph_vit_destroy(&vit);
IGRAPH_FINALLY_CLEAN(3);
