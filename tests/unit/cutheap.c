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

#include "core/cutheap.h"

#include "test_utilities.inc"

int main() {
    igraph_i_cutheap_t ch;
    long int i;

    igraph_i_cutheap_init(&ch, 10);

    for (i = 0; i < 10; i++) {
        igraph_i_cutheap_update(&ch, i, i);
    }
    while (!igraph_i_cutheap_empty(&ch)) {
        long int idx = igraph_i_cutheap_popmax(&ch);
        printf("%li ", idx);
    }
    printf("\n");

    igraph_i_cutheap_destroy(&ch);

    VERIFY_FINALLY_STACK();

    return 0;
}
