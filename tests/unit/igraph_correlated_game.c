/* -*- mode: C -*-  */
/* vim:set ts=4 sw=4 sts=4 et: */
/*
   IGraph R library.
   Copyright (C) 2003-2013  Gabor Csardi <csardi.gabor@gmail.com>
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

#include "test_utilities.inc"

int main() {
    igraph_t g1, g2;

    igraph_rng_seed(igraph_rng_default(), 9275);

    igraph_erdos_renyi_game(&g1, IGRAPH_ERDOS_RENYI_GNP, 10, .3,
                            IGRAPH_UNDIRECTED, IGRAPH_NO_LOOPS);
    igraph_correlated_game(&g1, &g2, .9, .3, /* permutation=*/ 0);

    IGRAPH_ASSERT(igraph_vcount(&g1) == igraph_vcount(&g2));

    igraph_destroy(&g2);
    igraph_destroy(&g1);

    VERIFY_FINALLY_STACK();

    return 0;
}
