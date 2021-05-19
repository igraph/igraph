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
   Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA
   02110-1301 USA

*/

#include <igraph.h>
#include <stdio.h>

int main() {

    igraph_t g;
    igraph_bool_t simple;

    igraph_barabasi_game(/* graph=    */ &g,
                                         /* n=        */ 100,
                                         /* power=    */ 1.0,
                                         /* m=        */ 2,
                                         /* outseq=   */ 0,
                                         /* outpref=  */ 0,
                                         /* A=        */ 1.0,
                                         /* directed= */ IGRAPH_DIRECTED,
                                         /* algo=     */ IGRAPH_BARABASI_PSUMTREE,
                                         /* start_from= */ 0);

    if (igraph_ecount(&g) != 197) {
        return 1;
    }
    if (igraph_vcount(&g) != 100) {
        return 2;
    }
    igraph_is_simple(&g, &simple);
    if (!simple) {
        return 3;
    }

    igraph_destroy(&g);

    /* ============================== */

    igraph_barabasi_game(/* graph=    */ &g,
                                         /* n=        */ 100,
                                         /* power=    */ 1.0,
                                         /* m=        */ 2,
                                         /* outseq=   */ 0,
                                         /* outpref=  */ 0,
                                         /* A=        */ 1.0,
                                         /* directed= */ IGRAPH_DIRECTED,
                                         /* algo=     */ IGRAPH_BARABASI_PSUMTREE_MULTIPLE,
                                         /* start_from= */ 0);

    if (igraph_ecount(&g) != 198) {
        return 4;
    }
    if (igraph_vcount(&g) != 100) {
        return 5;
    }
    igraph_is_simple(&g, &simple);
    if (simple) {
        return 6;
    }

    igraph_destroy(&g);

    /* ============================== */

    igraph_barabasi_game(/* graph=    */ &g,
                                         /* n=        */ 100,
                                         /* power=    */ 1.0,
                                         /* m=        */ 2,
                                         /* outseq=   */ 0,
                                         /* outpref=  */ 0,
                                         /* A=        */ 1.0,
                                         /* directed= */ IGRAPH_DIRECTED,
                                         /* algo=     */ IGRAPH_BARABASI_BAG,
                                         /* start_from= */ 0);

    if (igraph_ecount(&g) != 198) {
        return 7;
    }
    if (igraph_vcount(&g) != 100) {
        return 8;
    }
    igraph_is_simple(&g, &simple);
    if (simple) {
        return 9;
    }

    igraph_destroy(&g);

    return 0;
}
