/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2021  The igraph development team <igraph@igraph.org>

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

#include <igraph.h>

#include "test_utilities.h"

int main(void) {
    igraph_t g;
    igraph_real_t result;

    // recreate example graph from paper. A=0, B=1, etc. Note that O is missing.
    // AC, BC, CD, CE, CF, DF, EF, EG, FG, FH, FJ, GI,
    // GJ, HK, JL, LM, LN, LP, NP, PQ, PR, QS, RT
    igraph_small(&g, 19, IGRAPH_UNDIRECTED, 
        0,2, 1,2, 2,3, 2,4, 2,5, 3,5, 4,5, 4,6, 5,6, 5,7, 5,9, 6,8, 6,9, 7,10,
        9,11, 11,12, 11,13, 11,14, 13,14, 14,15, 14,16, 15,17, 16,18, -1
    );
    igraph_diameter_bound(&g, &result, 5, NULL, NULL, IGRAPH_UNDIRECTED, 0);

    igraph_destroy(&g);

    VERIFY_FINALLY_STACK();

    return 0;
}
