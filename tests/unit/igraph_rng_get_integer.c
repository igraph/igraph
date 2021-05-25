/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2011-2021  The igraph development team

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
    int i;

	/* Seed the RNG, generate 10 random integers */
    igraph_rng_seed(igraph_rng_default(), 42);
    for (i = 0; i < 10; i++) {
        printf("%ld\n", igraph_rng_get_integer(igraph_rng_default(), 10, 100));
    }

	printf("========\n");

	/* Seed the RNG again with the same seed, verify that we get the same
	 * numbers */
    igraph_rng_seed(igraph_rng_default(), 42);
    for (i = 0; i < 10; i++) {
        printf("%ld\n", igraph_rng_get_integer(igraph_rng_default(), 10, 100));
    }

	printf("========\n");

	/* Seed the RNG again with a different seed, verify that we get different
	 * numbers */
    igraph_rng_seed(igraph_rng_default(), 84);
    for (i = 0; i < 10; i++) {
        printf("%ld\n", igraph_rng_get_integer(igraph_rng_default(), 10, 100));
    }

	printf("========\n");

    return 0;
}
