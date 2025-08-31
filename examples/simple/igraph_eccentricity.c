/*
   igraph library.
   Copyright (C) 2024  The igraph development team <igraph@igraph.org>

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

int main(void) {

    igraph_t g;
    igraph_vector_t ecc;

    /* Initialize the library. */
    igraph_setup();

    igraph_vector_init(&ecc, 0);

    igraph_star(&g, 10, IGRAPH_STAR_UNDIRECTED, 0);
    igraph_eccentricity(&g, NULL, &ecc, igraph_vss_all(), IGRAPH_OUT);
    igraph_vector_print(&ecc);
    igraph_destroy(&g);

    igraph_star(&g, 10, IGRAPH_STAR_OUT, 0);
    igraph_eccentricity(&g, NULL, &ecc, igraph_vss_all(), IGRAPH_ALL);
    igraph_vector_print(&ecc);
    igraph_destroy(&g);

    igraph_star(&g, 10, IGRAPH_STAR_OUT, 0);
    igraph_eccentricity(&g, NULL, &ecc, igraph_vss_all(), IGRAPH_OUT);
    igraph_vector_print(&ecc);
    igraph_destroy(&g);

    igraph_vector_destroy(&ecc);

    return 0;
}
