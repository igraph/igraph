/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2003-2020  The igraph development team

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

#include "igraph_types.h"

#include "core/cutheap.h"

#define PARENT(x)     ((x)/2)
#define LEFTCHILD(x)  ((x)*2+1)
#define RIGHTCHILD(x) ((x)*2)
#define INACTIVE      IGRAPH_INFINITY
#define UNDEFINED     0.0
#define INDEXINC      1

static void igraph_i_cutheap_switch(igraph_i_cutheap_t *ch,
                                    long int hidx1, long int hidx2) {
    if (hidx1 != hidx2) {
        long int idx1 = (long int) VECTOR(ch->index)[hidx1];
        long int idx2 = (long int) VECTOR(ch->index)[hidx2];

        igraph_real_t tmp = VECTOR(ch->heap)[hidx1];
        VECTOR(ch->heap)[hidx1] = VECTOR(ch->heap)[hidx2];
        VECTOR(ch->heap)[hidx2] = tmp;

        VECTOR(ch->index)[hidx1] = idx2;
        VECTOR(ch->index)[hidx2] = idx1;

        VECTOR(ch->hptr)[idx1] = hidx2 + INDEXINC;
        VECTOR(ch->hptr)[idx2] = hidx1 + INDEXINC;
    }
}

static void igraph_i_cutheap_sink(igraph_i_cutheap_t *ch, long int hidx) {
    long int size = igraph_vector_size(&ch->heap);
    if (LEFTCHILD(hidx) >= size) {
        /* leaf node */
    } else if (RIGHTCHILD(hidx) == size ||
               VECTOR(ch->heap)[LEFTCHILD(hidx)] >=
               VECTOR(ch->heap)[RIGHTCHILD(hidx)]) {
        /* sink to the left if needed */
        if (VECTOR(ch->heap)[hidx] < VECTOR(ch->heap)[LEFTCHILD(hidx)]) {
            igraph_i_cutheap_switch(ch, hidx, LEFTCHILD(hidx));
            igraph_i_cutheap_sink(ch, LEFTCHILD(hidx));
        }
    } else {
        /* sink to the right */
        if (VECTOR(ch->heap)[hidx] < VECTOR(ch->heap)[RIGHTCHILD(hidx)]) {
            igraph_i_cutheap_switch(ch, hidx, RIGHTCHILD(hidx));
            igraph_i_cutheap_sink(ch, RIGHTCHILD(hidx));
        }
    }
}

static void igraph_i_cutheap_shift_up(igraph_i_cutheap_t *ch, long int hidx) {
    if (hidx == 0 || VECTOR(ch->heap)[hidx] < VECTOR(ch->heap)[PARENT(hidx)]) {
        /* at the top */
    } else {
        igraph_i_cutheap_switch(ch, hidx, PARENT(hidx));
        igraph_i_cutheap_shift_up(ch, PARENT(hidx));
    }
}

int igraph_i_cutheap_init(igraph_i_cutheap_t *ch, igraph_integer_t nodes) {
    ch->dnodes = nodes;
    IGRAPH_VECTOR_INIT_FINALLY(&ch->heap, nodes); /* all zero */
    IGRAPH_CHECK(igraph_vector_init_seq(&ch->index, 0, nodes - 1));
    IGRAPH_FINALLY(igraph_vector_destroy, &ch->index);
    IGRAPH_CHECK(igraph_vector_init_seq(&ch->hptr, INDEXINC, nodes + INDEXINC - 1));
    IGRAPH_FINALLY_CLEAN(2);
    return 0;
}

void igraph_i_cutheap_destroy(igraph_i_cutheap_t *ch) {
    igraph_vector_destroy(&ch->hptr);
    igraph_vector_destroy(&ch->index);
    igraph_vector_destroy(&ch->heap);
}

igraph_bool_t igraph_i_cutheap_empty(igraph_i_cutheap_t *ch) {
    return igraph_vector_empty(&ch->heap);
}

/* Number of active vertices */

igraph_integer_t igraph_i_cutheap_active_size(igraph_i_cutheap_t *ch) {
    return (igraph_integer_t) igraph_vector_size(&ch->heap);
}

/* Number of all (defined) vertices */

igraph_integer_t igraph_i_cutheap_size(igraph_i_cutheap_t *ch) {
    return (igraph_integer_t) (ch->dnodes);
}

igraph_real_t igraph_i_cutheap_maxvalue(igraph_i_cutheap_t *ch) {
    return VECTOR(ch->heap)[0];
}

igraph_integer_t igraph_i_cutheap_popmax(igraph_i_cutheap_t *ch) {
    long int size = igraph_vector_size(&ch->heap);
    igraph_integer_t maxindex = (igraph_integer_t) VECTOR(ch->index)[0];
    /* put the last element to the top */
    igraph_i_cutheap_switch(ch, 0, size - 1);
    /* remove the last element */
    VECTOR(ch->hptr)[(long int) igraph_vector_tail(&ch->index)] = INACTIVE;
    igraph_vector_pop_back(&ch->heap);
    igraph_vector_pop_back(&ch->index);
    igraph_i_cutheap_sink(ch, 0);

    return maxindex;
}

/* Update the value of an active vertex, if not active it will be ignored */

int igraph_i_cutheap_update(igraph_i_cutheap_t *ch, igraph_integer_t index,
                            igraph_real_t add) {
    igraph_real_t hidx = VECTOR(ch->hptr)[(long int)index];
    if (hidx != INACTIVE && hidx != UNDEFINED) {
        long int hidx2 = (long int) (hidx - INDEXINC);
        /*     printf("updating vertex %li, heap index %li\n", (long int) index, hidx2); */
        VECTOR(ch->heap)[hidx2] += add;
        igraph_i_cutheap_sink(ch, hidx2);
        igraph_i_cutheap_shift_up(ch, hidx2);
    }
    return 0;
}

/* Reset the value of all vertices to zero and make them active */

int igraph_i_cutheap_reset_undefine(igraph_i_cutheap_t *ch, long int vertex) {
    long int i, j, n = igraph_vector_size(&ch->hptr);
    /* undefine */
    VECTOR(ch->hptr)[vertex] = UNDEFINED;
    ch->dnodes -= 1;

    IGRAPH_CHECK(igraph_vector_resize(&ch->heap, ch->dnodes));
    igraph_vector_null(&ch->heap);

    IGRAPH_CHECK(igraph_vector_resize(&ch->index, ch->dnodes));

    j = 0;
    for (i = 0; i < n; i++) {
        if (VECTOR(ch->hptr)[i] != UNDEFINED) {
            VECTOR(ch->index)[j] = i;
            VECTOR(ch->hptr)[i] = j + INDEXINC;
            j++;
        }
    }

    return 0;
}
