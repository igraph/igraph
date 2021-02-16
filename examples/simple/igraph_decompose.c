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
   Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA
   02110-1301 USA

*/

#include <igraph.h>
#include <stdlib.h>

void free_complist(igraph_vector_ptr_t *complist) {
    long int i;
    for (i = 0; i < igraph_vector_ptr_size(complist); i++) {
        igraph_destroy(VECTOR(*complist)[i]);
        igraph_free(VECTOR(*complist)[i]);
    }
}

int main() {

    igraph_t ring, g;
    igraph_vector_ptr_t complist;
    long int i;
    igraph_real_t edges[] = { 0, 1, 1, 2, 2, 0,
                              3, 4, 4, 5, 5, 6,
                              8, 9, 9, 10
                            };
    igraph_vector_t v;

    /* A ring, a single component */
    igraph_ring(&ring, 10, IGRAPH_UNDIRECTED, 0, 1);

    igraph_vector_ptr_init(&complist, 0);
    igraph_decompose(&ring, &complist, IGRAPH_WEAK, -1, 0);
    igraph_write_graph_edgelist(VECTOR(complist)[0], stdout);
    free_complist(&complist);
    igraph_destroy(&ring);

    /* random graph with a giant component */
    igraph_erdos_renyi_game(&g, IGRAPH_ERDOS_RENYI_GNP, 100, 4.0 / 100,
                            IGRAPH_UNDIRECTED, 0);
    igraph_decompose(&g, &complist, IGRAPH_WEAK, -1, 20);
    if (igraph_vector_ptr_size(&complist) != 1) {
        return 1;
    }
    free_complist(&complist);
    igraph_destroy(&g);

    /* a toy graph, three components maximum, with at least 2 vertices each */
    igraph_create(&g,
                  igraph_vector_view(&v, edges, sizeof(edges) / sizeof(igraph_real_t)),
                  0, IGRAPH_DIRECTED);
    igraph_decompose(&g, &complist, IGRAPH_WEAK, 3, 2);
    for (i = 0; i < igraph_vector_ptr_size(&complist); i++) {
        igraph_write_graph_edgelist(VECTOR(complist)[i], stdout);
    }
    free_complist(&complist);
    igraph_destroy(&g);

    /* The same graph, this time with vertex attributes */
    /*   igraph_vector_init_seq(&idvect, 0, igraph_vcount(&g)-1); */
    /*   igraph_add_vertex_attribute(&g, "id", IGRAPH_ATTRIBUTE_NUM); */
    /*   igraph_set_vertex_attributes(&g, "id", IGRAPH_VS_ALL(&g), &idvect); */
    /*   igraph_vector_destroy(&idvect); */

    /*   igraph_decompose(&g, &complist, IGRAPH_WEAK, 3, 2); */
    /*   for (i=0; i<igraph_vector_ptr_size(&complist); i++) { */
    /*     igraph_t *comp=VECTOR(complist)[i]; */
    /*     igraph_es_t es; */
    /*     igraph_es_all(comp, &es); */
    /*     while (!igraph_es_end(comp, &es)) { */
    /*       igraph_real_t *from, *to;  */
    /*       igraph_get_vertex_attribute(comp, "id", igraph_es_from(comp, &es), */
    /*                (void**) &from, 0); */
    /*       igraph_get_vertex_attribute(comp, "id", igraph_es_to(comp, &es), */
    /*                (void**) &to, 0); */
    /*       printf("%li %li\n", (long int) *from, (long int) *to); */
    /*       igraph_es_next(comp, &es); */
    /*     } */
    /*   } */

    /*   free_complist(&complist); */
    /*   igraph_destroy(&g);   */

    igraph_vector_ptr_destroy(&complist);

    return 0;
}
