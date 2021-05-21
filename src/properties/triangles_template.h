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

#ifdef TRANSIT
#define TRANSIT_TRIEDGES
#endif

long int no_of_nodes = igraph_vcount(graph);
long int node, i, j, nn;
igraph_adjlist_t allneis;
igraph_vector_int_t *neis1, *neis2;
long int neilen1, neilen2;
long int *neis;
long int maxdegree;

#ifdef TRANSIT_TRIEDGES
long int deg1;
#endif

igraph_vector_int_t order;
igraph_vector_int_t rank;
igraph_vector_t degree;

if (no_of_nodes == 0) {
#ifndef TRIANGLES
    igraph_vector_clear(res);
#else
    igraph_vector_int_clear(res);
#endif
    return IGRAPH_SUCCESS;
}

igraph_vector_int_init(&order, no_of_nodes);
IGRAPH_FINALLY(igraph_vector_int_destroy, &order);
IGRAPH_VECTOR_INIT_FINALLY(&degree, no_of_nodes);

IGRAPH_CHECK(igraph_adjlist_init(graph, &allneis, IGRAPH_ALL, IGRAPH_NO_LOOPS, IGRAPH_NO_MULTIPLE));
IGRAPH_FINALLY(igraph_adjlist_destroy, &allneis);

for (i = 0; i < no_of_nodes; i++) {
    VECTOR(degree)[i] = igraph_vector_int_size(igraph_adjlist_get(&allneis, i));
}

maxdegree = (long int) igraph_vector_max(&degree) + 1;
igraph_vector_order1_int(&degree, &order, maxdegree);
igraph_vector_int_init(&rank, no_of_nodes);
IGRAPH_FINALLY(igraph_vector_int_destroy, &rank);
for (i = 0; i < no_of_nodes; i++) {
    VECTOR(rank)[ VECTOR(order)[i] ] = no_of_nodes - i - 1;
}

IGRAPH_CHECK(igraph_i_trans4_al_simplify(&allneis, &rank));

neis = IGRAPH_CALLOC(no_of_nodes, long int);
if (neis == 0) {
    IGRAPH_ERROR("undirected local transitivity failed", IGRAPH_ENOMEM);
}
IGRAPH_FINALLY(igraph_free, neis);

#ifndef TRIANGLES
    IGRAPH_CHECK(igraph_vector_resize(res, no_of_nodes));
    igraph_vector_null(res);
#else
    igraph_vector_int_clear(res);
#endif

for (nn = no_of_nodes - 1; nn >= 0; nn--) {
    node = VECTOR(order)[nn];

    IGRAPH_ALLOW_INTERRUPTION();

    neis1 = igraph_adjlist_get(&allneis, node);
    neilen1 = igraph_vector_int_size(neis1);

#ifdef TRANSIT_TRIEDGES
    deg1 = (long int) VECTOR(degree)[node];
#endif

    /* Mark the neighbors of the node */
    for (i = 0; i < neilen1; i++) {
        neis[ (long int) VECTOR(*neis1)[i] ] = node + 1;
    }

    for (i = 0; i < neilen1; i++) {
        long int nei = (long int) VECTOR(*neis1)[i];
        neis2 = igraph_adjlist_get(&allneis, nei);
        neilen2 = igraph_vector_int_size(neis2);
        for (j = 0; j < neilen2; j++) {
            long int nei2 = (long int) VECTOR(*neis2)[j];
            if (neis[nei2] == node + 1) {
#ifndef TRIANGLES
                VECTOR(*res)[nei2] += 1;
                VECTOR(*res)[nei] += 1;
                VECTOR(*res)[node] += 1;
#else
                IGRAPH_CHECK(igraph_vector_int_push_back(res, node));
                IGRAPH_CHECK(igraph_vector_int_push_back(res, nei));
                IGRAPH_CHECK(igraph_vector_int_push_back(res, nei2));
#endif
            }
        }
    }

#ifdef TRANSIT
    if (mode == IGRAPH_TRANSITIVITY_ZERO && deg1 < 2) {
        VECTOR(*res)[node] = 0.0;
    } else {
        VECTOR(*res)[node] = VECTOR(*res)[node] / deg1 / (deg1 - 1) * 2.0;
    }
#endif
}

igraph_free(neis);
igraph_adjlist_destroy(&allneis);
igraph_vector_int_destroy(&rank);
igraph_vector_destroy(&degree);
igraph_vector_int_destroy(&order);
IGRAPH_FINALLY_CLEAN(5);

#ifdef TRANSIT_TRIEDGES
#undef TRANSIT_TRIEDGES
#endif
