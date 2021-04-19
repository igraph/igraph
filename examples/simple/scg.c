/* -*- mode: C -*-  */
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

int main() {

    igraph_t g;
    igraph_vector_t ev;
    igraph_t scg_graph;
    igraph_matrix_t scg_matrix;
    igraph_sparsemat_t scg_sparsemat;
    igraph_matrix_t L, R;
    igraph_sparsemat_t Lsparse, Rsparse;
    igraph_matrix_t input_matrix;
    igraph_vector_t groups;
    igraph_vector_t eval;
    igraph_matrix_t evec;

    igraph_tree(&g, 10, /* children= */ 3, IGRAPH_TREE_UNDIRECTED);

    igraph_vector_init(&ev, 1);
    igraph_matrix_init(&L, 0, 0);
    igraph_matrix_init(&R, 0, 0);
    igraph_matrix_init(&scg_matrix, 0, 0);
    igraph_vector_init(&groups, 0);
    igraph_vector_init(&eval, 0);
    igraph_matrix_init(&evec, 0, 0);

#define CALLSYM(algo) do {                      \
        igraph_vector_clear(&eval);                     \
        igraph_matrix_resize(&evec, 0, 0);                  \
        igraph_scg_adjacency(&g, /*matrix=*/ 0, /*sparsemat=*/ 0, &ev,  \
                             /* intervals= */ 3, /* intervals_vector= */ 0, \
                             /* algorithm= */ algo, &eval, &evec,       \
                             /* groups= */ &groups, /* use_arpack= */ 0,    \
                             /* maxiter= */ 0, &scg_graph, &scg_matrix, \
                             &scg_sparsemat, &L, &R,            \
                             &Lsparse, &Rsparse); } while(0)


#define PRINTRES()                      \
    do {                              \
        printf("------------------------------------\n");       \
        igraph_write_graph_edgelist(&scg_graph, stdout);        \
        printf("---\n");                        \
        igraph_vector_print(&groups);               \
        printf("---\n");                        \
        igraph_vector_print(&eval);                 \
        igraph_matrix_print(&evec);                 \
        printf("---\n");                        \
        igraph_sparsemat_print(&scg_sparsemat, stdout);     \
        printf("---\n");                        \
        igraph_sparsemat_print(&Lsparse, stdout);           \
        printf("---\n");                        \
        igraph_sparsemat_print(&Rsparse, stdout);           \
        printf("---\n");                        \
    } while (0)

    VECTOR(ev)[0] = 1;
    CALLSYM(IGRAPH_SCG_EXACT);
    PRINTRES();
    igraph_destroy(&scg_graph);
    igraph_sparsemat_destroy(&scg_sparsemat);
    igraph_sparsemat_destroy(&Lsparse);
    igraph_sparsemat_destroy(&Rsparse);

    VECTOR(ev)[0] = 3;
    CALLSYM(IGRAPH_SCG_EXACT);
    PRINTRES();
    igraph_destroy(&scg_graph);
    igraph_sparsemat_destroy(&scg_sparsemat);
    igraph_sparsemat_destroy(&Lsparse);
    igraph_sparsemat_destroy(&Rsparse);

    igraph_vector_resize(&ev, 2);
    VECTOR(ev)[0] = 1;
    VECTOR(ev)[1] = 3;
    CALLSYM(IGRAPH_SCG_EXACT);
    PRINTRES();
    igraph_destroy(&scg_graph);
    igraph_sparsemat_destroy(&scg_sparsemat);
    igraph_sparsemat_destroy(&Lsparse);
    igraph_sparsemat_destroy(&Rsparse);

#define CALLSYM2(algo) do {                     \
        igraph_vector_clear(&eval);                     \
        igraph_matrix_resize(&evec, 0, 0);                  \
        igraph_scg_adjacency(/* graph=*/ 0, &input_matrix, /*sparsemat=*/ 0, \
                                         &ev, /* intervals= */ 3,           \
                                         /* intervals_vector= */ 0,         \
                                         /* algorithm= */ algo, &eval, &evec,       \
                                         /* groups= */ &groups, /* use_arpack= */ 0,    \
                                         /* maxiter= */ 0, &scg_graph, &scg_matrix, \
                                         &scg_sparsemat, &L, &R,            \
                                         &Lsparse, &Rsparse); } while (0)

    igraph_matrix_init(&input_matrix, 0, 0);
    igraph_get_adjacency(&g, &input_matrix, IGRAPH_GET_ADJACENCY_BOTH,
                         /* eids= */ 0);

    igraph_vector_resize(&ev, 1);
    VECTOR(ev)[0] = 1;
    CALLSYM2(IGRAPH_SCG_EXACT);
    PRINTRES();
    igraph_destroy(&scg_graph);
    igraph_sparsemat_destroy(&scg_sparsemat);
    igraph_sparsemat_destroy(&Lsparse);
    igraph_sparsemat_destroy(&Rsparse);

    VECTOR(ev)[0] = 3;
    CALLSYM2(IGRAPH_SCG_EXACT);
    PRINTRES();
    igraph_destroy(&scg_graph);
    igraph_sparsemat_destroy(&scg_sparsemat);
    igraph_sparsemat_destroy(&Lsparse);
    igraph_sparsemat_destroy(&Rsparse);

    igraph_vector_resize(&ev, 2);
    VECTOR(ev)[0] = 1;
    VECTOR(ev)[1] = 3;
    CALLSYM2(IGRAPH_SCG_EXACT);
    PRINTRES();
    igraph_destroy(&scg_graph);
    igraph_sparsemat_destroy(&scg_sparsemat);
    igraph_sparsemat_destroy(&Lsparse);
    igraph_sparsemat_destroy(&Rsparse);

    igraph_matrix_destroy(&evec);
    igraph_vector_destroy(&eval);
    igraph_vector_destroy(&groups);
    igraph_matrix_destroy(&input_matrix);
    igraph_matrix_destroy(&scg_matrix);
    igraph_matrix_destroy(&L);
    igraph_matrix_destroy(&R);
    igraph_vector_destroy(&ev);
    igraph_destroy(&g);

    /* -------------------------------------------------------------------- */

    return 0;
}
