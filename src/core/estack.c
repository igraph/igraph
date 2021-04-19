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
   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
   02110-1301 USA

*/

#include "core/estack.h"

int igraph_estack_init(igraph_estack_t *s, long int setsize,
                       long int stacksize) {
    IGRAPH_CHECK(igraph_vector_bool_init(&s->isin, setsize));
    IGRAPH_FINALLY(igraph_vector_bool_destroy, &s->isin);
    IGRAPH_CHECK(igraph_stack_long_init(&s->stack, stacksize));
    IGRAPH_FINALLY_CLEAN(1);
    return 0;
}

void igraph_estack_destroy(igraph_estack_t *s) {
    igraph_stack_long_destroy(&s->stack);
    igraph_vector_bool_destroy(&s->isin);
}

int igraph_estack_push(igraph_estack_t *s,  long int elem) {
    if ( !VECTOR(s->isin)[elem] ) {
        IGRAPH_CHECK(igraph_stack_long_push(&s->stack, elem));
        VECTOR(s->isin)[elem] = 1;
    }
    return 0;
}

long int igraph_estack_pop(igraph_estack_t *s) {
    long int elem = igraph_stack_long_pop(&s->stack);
    VECTOR(s->isin)[elem] = 0;
    return elem;
}

igraph_bool_t igraph_estack_iselement(const igraph_estack_t *s,
                                      long int elem) {
    return VECTOR(s->isin)[elem];
}

long int igraph_estack_size(const igraph_estack_t *s) {
    return igraph_stack_long_size(&s->stack);
}

#ifndef USING_R
int igraph_estack_print(const igraph_estack_t *s) {
    return igraph_stack_long_print(&s->stack);
}
#endif
