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

void print_vector(igraph_vector_t *v, FILE *f) {
    long int i;
    for (i = 0; i < igraph_vector_size(v); i++) {
        fprintf(f, " %li", (long int) VECTOR(*v)[i]);
    }
    fprintf(f, "\n");
}

int main() {
    igraph_t g;
    igraph_integer_t eid;
    igraph_vector_t hist;
    long int i;
    int ret;

    /* DIRECTED */

    igraph_star(&g, 10, IGRAPH_STAR_OUT, 0);

    igraph_vector_init(&hist, 9);

    for (i = 1; i < 10; i++) {
        igraph_get_eid(&g, &eid, 0, i, IGRAPH_DIRECTED, /*error=*/ 1);
        VECTOR(hist)[ (long int) eid ] = 1;
    }
    print_vector(&hist, stdout);

    igraph_vector_destroy(&hist);
    igraph_destroy(&g);

    /* UNDIRECTED */

    igraph_star(&g, 10, IGRAPH_STAR_UNDIRECTED, 0);

    igraph_vector_init(&hist, 9);

    for (i = 1; i < 10; i++) {
        igraph_get_eid(&g, &eid, 0, i, IGRAPH_UNDIRECTED, /*error=*/ 1);
        VECTOR(hist)[ (long int) eid ] += 1;
        igraph_get_eid(&g, &eid, i, 0, IGRAPH_DIRECTED, /*error=*/ 1);
        VECTOR(hist)[ (long int) eid ] += 1;
    }
    print_vector(&hist, stdout);

    igraph_vector_destroy(&hist);
    igraph_destroy(&g);

    /* NON-EXISTANT EDGE */

    igraph_star(&g, 10, IGRAPH_STAR_UNDIRECTED, 0);

    igraph_set_error_handler(igraph_error_handler_ignore);

    ret = igraph_get_eid(&g, &eid, 5, 6, IGRAPH_UNDIRECTED, /*error=*/ 1);
    if (ret != IGRAPH_EINVAL) {
        return 1;
    }

    igraph_destroy(&g);

    return 0;
}

/* Stress test */

/* int main() { */

/*   igraph_t g; */
/*   long int i, n; */
/*   igraph_integer_t from, to, eid; */

/*   igraph_barabasi_game(&g, 10000, 100, 0, 0, 1); */
/*   n=igraph_ecount(&g); */
/*   for (i=0; i<n; i++) { */
/*     igraph_edge(&g, i, &from, &to); */
/*     igraph_get_eid(&g, &eid, from, to, 1, 1); */
/*     igraph_get_eid(&g, &eid, to, from, 0, 1); */
/*     igraph_get_eid(&g, &eid, from, to, 0, 1); */
/*   } */
/*   igraph_destroy(&g); */

/*   igraph_barabasi_game(&g, 10000, 100, 0, 0, 0); */
/*   n=igraph_ecount(&g); */
/*   for (i=0; i<n; i++) { */
/*     igraph_edge(&g, i, &from, &to); */
/*     igraph_get_eid(&g, &eid, from, to, 0, 1); */
/*     igraph_get_eid(&g, &eid, to, from, 0, 1); */
/*   } */
/*   igraph_destroy(&g); */

/*   igraph_erdos_renyi_game(&g, IGRAPH_ERDOS_RENYI_GNP, */
/*            2000, 100.0/2000, 0, 0); */
/*   n=igraph_ecount(&g); */
/*   for (i=0; i<n; i++) { */
/*     igraph_edge(&g, i, &from, &to); */
/*     igraph_get_eid(&g, &eid, from, to, 0, 1); */
/*     igraph_get_eid(&g, &eid, to, from, 0, 1); */
/*   } */
/*   igraph_destroy(&g); */

/*   igraph_full(&g, 500, 0, 0); */
/*   n=igraph_ecount(&g); */
/*   for (i=0; i<n; i++) { */
/*     igraph_edge(&g, i, &from, &to); */
/*     igraph_get_eid(&g, &eid, from, to, 0, 1); */
/*   } */
/*   igraph_destroy(&g); */

/*   igraph_star(&g, 20000, IGRAPH_STAR_OUT, 0); */
/*   n=igraph_ecount(&g); */
/*   for (i=0; i<n; i++) { */
/*     igraph_edge(&g, i, &from, &to); */
/*     igraph_get_eid(&g, &eid, from, to, 0, 1); */
/*   } */
/*   igraph_destroy(&g); */

/*   igraph_star(&g, 20000, IGRAPH_STAR_IN, 0); */
/*   n=igraph_ecount(&g); */
/*   for (i=0; i<n; i++) { */
/*     igraph_edge(&g, i, &from, &to); */
/*     igraph_get_eid(&g, &eid, from, to, 0, 1); */
/*   } */
/*   igraph_destroy(&g); */

/*   igraph_star(&g, 2000000, IGRAPH_STAR_UNDIRECTED, 1999999); */
/*   n=igraph_ecount(&g); */
/*   for (i=0; i<n; i++) { */
/*     igraph_edge(&g, i, &from, &to); */
/*     igraph_get_eid(&g, &eid, from, to, 0, 1); */
/*     igraph_get_eid(&g, &eid, to, from, 0, 1); */
/*   } */
/*   igraph_destroy(&g); */

/*   return 0; */
/* } */
