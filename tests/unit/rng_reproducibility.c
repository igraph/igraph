/*
   igraph library.
   Copyright (C) 2019-2022  The igraph development team <igraph@igraph.org>

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

/*
 * This test serves to ensure that the same sequence of random numbers are generated for the
 * same seed on all platforms (different operating systems and 32- or 64-bit systems).
 */

int main(void) {
    int i;
    igraph_rng_seed(igraph_rng_default(), 137);

    for (i = 0; i < 32; ++i) {
        printf("%" IGRAPH_PRId "\n", RNG_INTEGER(0, 100));
    }

    for (i = 0; i < 32; ++i) {
        printf("%g\n", RNG_UNIF(0, 1e-6));
    }
    return 0;
}
