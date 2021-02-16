/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2010-2012  Gabor Csardi <csardi.gabor@gmail.com>
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

#include <igraph.h>

#include "core/marked_queue.h"
#include "core/estack.h"

#include "flow/flow_internal.h"

#include "test_utilities.inc"

int test_all_st_cuts(const igraph_t *graph,
                     long int source,
                     long int target) {
    igraph_vector_ptr_t cuts, partition1s;
    long int n, i;

    igraph_vector_ptr_init(&cuts, 0);
    igraph_vector_ptr_init(&partition1s, 0);
    igraph_all_st_cuts(graph, &cuts, &partition1s,
                       source, target);

    n = igraph_vector_ptr_size(&partition1s);
    printf("Partitions and cuts:\n");
    for (i = 0; i < n; i++) {
        igraph_vector_t *v = VECTOR(partition1s)[i];
        igraph_vector_t *v2 = VECTOR(cuts)[i];
        printf("P: ");
        igraph_vector_print(v);
        igraph_vector_destroy(v);
        igraph_free(v);
        printf("C: ");
        igraph_vector_print(v2);
        igraph_vector_destroy(v2);
        igraph_free(v2);
    }
    igraph_vector_ptr_destroy(&partition1s);
    igraph_vector_ptr_destroy(&cuts);

    return 0;
}

int main() {
    igraph_t g;
    igraph_vector_ptr_t cuts, partition1s;
    long int i, n;

    igraph_marked_queue_t S;
    igraph_estack_t T;
    long int v;
    igraph_vector_t Isv;

    /* ----------------------------------------------------------- */
    /* This is the example from the Provan-Shier paper,
       for calculating the dominator tree and finding the right pivot
       element */

    igraph_small(&g, 12, IGRAPH_DIRECTED,
                 /* a->b */ 0, 1,
                 /* b->t */ 1, 11,
                 /* c->b */ 2, 1,  /* c->d */ 2, 3,
                 /* d->e */ 3, 4,  /* d->i */ 3, 8,
                 /* e->c */ 4, 2,
                 /* f->c */ 5, 2,  /* f->e */ 5, 4,
                 /* g->d */ 6, 3,  /* g->e */ 6, 4,  /* g->f */ 6, 5,
                 /* g->j */ 6, 9,
                 /* h->g */ 7, 6,  /* h->t */ 7, 11,
                 /* i->a */ 8, 0,
                 /* j->i */ 9, 8,
                 /* s->a */ 10, 0, /* s->c */ 10, 2, /* s->h */ 10, 7,
                 -1);

    /* S={s,a} */
    igraph_marked_queue_init(&S, igraph_vcount(&g));
    igraph_marked_queue_start_batch(&S);
    igraph_marked_queue_push(&S, 10);
    igraph_marked_queue_push(&S, 0);

    /* T={t} */
    igraph_estack_init(&T, igraph_vcount(&g), 1);
    igraph_estack_push(&T, 11);

    igraph_vector_init(&Isv, 0);
    igraph_i_all_st_cuts_pivot(&g, &S, &T,
                               /*source=*/ 10, /*target=*/ 11,
                               &v, &Isv, NULL);

    /* Expected result: v=c, Isv={c,d,e,i} */
    printf("%li; ", v);
    igraph_vector_print(&Isv);

    igraph_vector_destroy(&Isv);
    igraph_estack_destroy(&T);
    igraph_marked_queue_destroy(&S);
    igraph_destroy(&g);

    /* ----------------------------------------------------------- */

    igraph_small(&g, 3, IGRAPH_DIRECTED,
                 0, 1, 1, 2,
                 -1);

    /* S={}, T={} */
    igraph_marked_queue_init(&S, igraph_vcount(&g));
    igraph_estack_init(&T, igraph_vcount(&g), 3);

    igraph_vector_init(&Isv, 0);
    igraph_i_all_st_cuts_pivot(&g, &S, &T,
                               /*source=*/ 0, /*target=*/ 2,
                               &v, &Isv, NULL);
    printf("%li; ", v);
    igraph_vector_print(&Isv);

    igraph_vector_destroy(&Isv);
    igraph_estack_destroy(&T);
    igraph_marked_queue_destroy(&S);
    igraph_destroy(&g);

    /* ----------------------------------------------------------- */

    igraph_small(&g, 3, IGRAPH_DIRECTED,
                 0, 1, 1, 2,
                 -1);

    /* S={}, T={0} */
    igraph_marked_queue_init(&S, igraph_vcount(&g));

    igraph_estack_init(&T, igraph_vcount(&g), 3);
    igraph_estack_push(&T, 0);

    igraph_vector_init(&Isv, 0);
    igraph_i_all_st_cuts_pivot(&g, &S, &T,
                               /*source=*/ 0, /*target=*/ 2,
                               &v, &Isv, NULL);
    printf("%li; ", v);
    igraph_vector_print(&Isv);

    igraph_vector_destroy(&Isv);
    igraph_estack_destroy(&T);
    igraph_marked_queue_destroy(&S);
    igraph_destroy(&g);

    /* ----------------------------------------------------------- */

    igraph_small(&g, 3, IGRAPH_DIRECTED,
                 0, 1, 1, 2,
                 -1);

    /* S={0}, T={} */
    igraph_marked_queue_init(&S, igraph_vcount(&g));
    igraph_marked_queue_push(&S, 0);

    igraph_estack_init(&T, igraph_vcount(&g), 3);

    igraph_vector_init(&Isv, 0);
    igraph_i_all_st_cuts_pivot(&g, &S, &T,
                               /*source=*/ 0, /*target=*/ 2,
                               &v, &Isv, NULL);
    printf("%li; ", v);
    igraph_vector_print(&Isv);

    igraph_vector_destroy(&Isv);
    igraph_estack_destroy(&T);
    igraph_marked_queue_destroy(&S);
    igraph_destroy(&g);

    /* ----------------------------------------------------------- */

    igraph_small(&g, 3, IGRAPH_DIRECTED,
                 0, 1, 1, 2,
                 -1);

    /* S={0}, T={1} */
    igraph_marked_queue_init(&S, igraph_vcount(&g));
    igraph_marked_queue_push(&S, 0);

    igraph_estack_init(&T, igraph_vcount(&g), 3);
    igraph_estack_push(&T, 1);

    igraph_vector_init(&Isv, 0);
    igraph_i_all_st_cuts_pivot(&g, &S, &T,
                               /*source=*/ 0, /*target=*/ 2,
                               &v, &Isv, NULL);
    printf("%li; ", v);
    igraph_vector_print(&Isv);

    igraph_vector_destroy(&Isv);
    igraph_estack_destroy(&T);
    igraph_marked_queue_destroy(&S);
    igraph_destroy(&g);

    /* ----------------------------------------------------------- */

    igraph_small(&g, 3, IGRAPH_DIRECTED,
                 0, 1, 1, 2,
                 -1);

    /* S={0,1}, T={} */
    igraph_marked_queue_init(&S, igraph_vcount(&g));
    igraph_marked_queue_push(&S, 0);
    igraph_marked_queue_push(&S, 1);

    igraph_estack_init(&T, igraph_vcount(&g), 3);

    igraph_vector_init(&Isv, 0);
    igraph_i_all_st_cuts_pivot(&g, &S, &T,
                               /*source=*/ 0, /*target=*/ 2,
                               &v, &Isv, NULL);
    printf("%li; ", v);
    igraph_vector_print(&Isv);

    igraph_vector_destroy(&Isv);
    igraph_estack_destroy(&T);
    igraph_marked_queue_destroy(&S);
    igraph_destroy(&g);

    /* ----------------------------------------------------------- */

    igraph_small(&g, 3, IGRAPH_DIRECTED,
                 0, 1, 1, 2,
                 -1);

    igraph_vector_ptr_init(&partition1s, 0);
    igraph_all_st_cuts(&g, /*cuts=*/ 0, &partition1s,
                       /*source=*/ 0, /*target=*/ 2);

    n = igraph_vector_ptr_size(&partition1s);
    for (i = 0; i < n; i++) {
        igraph_vector_t *v = VECTOR(partition1s)[i];
        igraph_vector_print(v);
        igraph_vector_destroy(v);
        igraph_free(v);
    }
    igraph_vector_ptr_destroy(&partition1s);

    igraph_destroy(&g);

    /* ----------------------------------------------------------- */

    igraph_small(&g, 5, IGRAPH_DIRECTED,
                 0, 1, 1, 2, 1, 3, 2, 4, 3, 4,
                 -1);

    igraph_vector_ptr_init(&partition1s, 0);
    igraph_all_st_cuts(&g, /*cuts=*/ 0, &partition1s,
                       /*source=*/ 0, /*target=*/ 4);

    n = igraph_vector_ptr_size(&partition1s);
    for (i = 0; i < n; i++) {
        igraph_vector_t *v = VECTOR(partition1s)[i];
        igraph_vector_print(v);
        igraph_vector_destroy(v);
        igraph_free(v);
    }
    igraph_vector_ptr_destroy(&partition1s);

    igraph_destroy(&g);

    /* ----------------------------------------------------------- */

    igraph_small(&g, 6, IGRAPH_DIRECTED,
                 0, 1, 1, 2, 1, 3, 2, 4, 3, 4, 1, 5, 5, 4,
                 -1);

    igraph_vector_ptr_init(&cuts, 0);
    igraph_vector_ptr_init(&partition1s, 0);
    igraph_all_st_cuts(&g, &cuts, &partition1s,
                       /*source=*/ 0, /*target=*/ 4);

    n = igraph_vector_ptr_size(&partition1s);
    printf("Partitions and cuts:\n");
    for (i = 0; i < n; i++) {
        igraph_vector_t *v = VECTOR(partition1s)[i];
        igraph_vector_t *v2 = VECTOR(cuts)[i];
        printf("P: ");
        igraph_vector_print(v);
        igraph_vector_destroy(v);
        igraph_free(v);
        printf("C: ");
        igraph_vector_print(v2);
        igraph_vector_destroy(v2);
        igraph_free(v2);
    }
    igraph_vector_ptr_destroy(&partition1s);
    igraph_vector_ptr_destroy(&cuts);

    igraph_destroy(&g);

    /* ----------------------------------------------------------- */

    igraph_small(&g, 3, IGRAPH_DIRECTED,
                 0, 2, 1, 2,
                 -1);

    igraph_vector_ptr_init(&cuts, 0);
    igraph_vector_ptr_init(&partition1s, 0);
    igraph_all_st_cuts(&g, &cuts, &partition1s,
                       /*source=*/ 1, /*target=*/ 2);

    n = igraph_vector_ptr_size(&partition1s);
    printf("Partitions and cuts:\n");
    for (i = 0; i < n; i++) {
        igraph_vector_t *v = VECTOR(partition1s)[i];
        igraph_vector_t *v2 = VECTOR(cuts)[i];
        printf("P: ");
        igraph_vector_print(v);
        igraph_vector_destroy(v);
        igraph_free(v);
        printf("C: ");
        igraph_vector_print(v2);
        igraph_vector_destroy(v2);
        igraph_free(v2);
    }
    igraph_vector_ptr_destroy(&partition1s);
    igraph_vector_ptr_destroy(&cuts);

    igraph_destroy(&g);

    /* ----------------------------------------------------------- */

    igraph_small(&g, 5, IGRAPH_DIRECTED,
                 0, 1, 1, 2, 2, 3, 3, 4, 3, 1,
                 -1);

    igraph_vector_ptr_init(&cuts, 0);
    igraph_vector_ptr_init(&partition1s, 0);
    igraph_all_st_cuts(&g, &cuts, &partition1s,
                       /*source=*/ 0, /*target=*/ 4);

    n = igraph_vector_ptr_size(&partition1s);
    printf("Partitions and cuts:\n");
    for (i = 0; i < n; i++) {
        igraph_vector_t *v = VECTOR(partition1s)[i];
        igraph_vector_t *v2 = VECTOR(cuts)[i];
        printf("P: ");
        igraph_vector_print(v);
        igraph_vector_destroy(v);
        igraph_free(v);
        printf("C: ");
        igraph_vector_print(v2);
        igraph_vector_destroy(v2);
        igraph_free(v2);
    }
    igraph_vector_ptr_destroy(&partition1s);
    igraph_vector_ptr_destroy(&cuts);

    igraph_destroy(&g);

    /* ----------------------------------------------------------- */

    igraph_small(&g, 7, IGRAPH_DIRECTED,
                 0, 1, 0, 2, 1, 3, 2, 3,
                 1, 4, 1, 5, 1, 6,
                 4, 2, 5, 2, 6, 2,
                 -1);

    igraph_vector_ptr_init(&cuts, 0);
    igraph_vector_ptr_init(&partition1s, 0);
    igraph_all_st_cuts(&g, &cuts, &partition1s,
                       /*source=*/ 0, /*target=*/ 3);

    n = igraph_vector_ptr_size(&partition1s);
    printf("Partitions and cuts:\n");
    for (i = 0; i < n; i++) {
        igraph_vector_t *v = VECTOR(partition1s)[i];
        igraph_vector_t *v2 = VECTOR(cuts)[i];
        printf("P: ");
        igraph_vector_print(v);
        igraph_vector_destroy(v);
        igraph_free(v);
        printf("C: ");
        igraph_vector_print(v2);
        igraph_vector_destroy(v2);
        igraph_free(v2);
    }
    igraph_vector_ptr_destroy(&partition1s);
    igraph_vector_ptr_destroy(&cuts);

    /* Check whether it also works if we don't provide partition1s */
    igraph_vector_ptr_init(&cuts, 0);
    igraph_vector_ptr_init(&partition1s, 0);
    igraph_all_st_cuts(&g, &cuts, /*partition1s=*/ 0,
                       /*source=*/ 0, /*target=*/ 3);

    n = igraph_vector_ptr_size(&cuts);
    printf("Cuts only (no partitions):\n");
    for (i = 0; i < n; i++) {
        igraph_vector_t *v2 = VECTOR(cuts)[i];
        printf("C: ");
        igraph_vector_print(v2);
        igraph_vector_destroy(v2);
        igraph_free(v2);
    }
    igraph_vector_ptr_destroy(&partition1s);
    igraph_vector_ptr_destroy(&cuts);

    igraph_destroy(&g);

    /* -----------------------------------------------------------
     * Check problematic cases in issue #1102
     * ----------------------------------------------------------- */

    igraph_small(&g, 4, IGRAPH_DIRECTED,
                 0, 1, 1, 2, 2, 3,
                 -1);
    test_all_st_cuts(&g, 0, 2);
    igraph_destroy(&g);

    igraph_small(&g, 5, IGRAPH_DIRECTED,
                 0, 1, 1, 2, 2, 3, 3, 4,
                 -1);
    test_all_st_cuts(&g, 0, 2);
    test_all_st_cuts(&g, 1, 3);
    igraph_destroy(&g);

    VERIFY_FINALLY_STACK();

    return 0;
}
