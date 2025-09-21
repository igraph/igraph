/*
   igraph library.
   Copyright (C) 2023  The igraph development team <igraph@igraph.org>

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
#include <float.h>

#include "test_utilities.h"

int main(void) {
    /* Test "normal" cases */
    IGRAPH_ASSERT(igraph_cmp_epsilon(1, 2, 0.25) < 0);
    IGRAPH_ASSERT(igraph_cmp_epsilon(1, 2, 0.5) == 0);
    IGRAPH_ASSERT(igraph_cmp_epsilon(2, 1, 0.25) > 0);
    IGRAPH_ASSERT(igraph_cmp_epsilon(2, 1, 0.5) == 0);
    IGRAPH_ASSERT(igraph_cmp_epsilon(10, 11, 0.25) == 0);
    IGRAPH_ASSERT(igraph_cmp_epsilon(10, 11, 0.5) == 0);

    /* Test infinities. Infinities are not equal to finite numbers with any
     * tolerance */
    IGRAPH_ASSERT(igraph_cmp_epsilon(IGRAPH_INFINITY, 2, 0.5) > 0);
    IGRAPH_ASSERT(igraph_cmp_epsilon(IGRAPH_INFINITY, 2, 20) > 0);
    IGRAPH_ASSERT(igraph_cmp_epsilon(IGRAPH_INFINITY, 2, 200) > 0);
    IGRAPH_ASSERT(igraph_cmp_epsilon(-IGRAPH_INFINITY, 2, 0.5) < 0);
    IGRAPH_ASSERT(igraph_cmp_epsilon(-IGRAPH_INFINITY, 2, 20) < 0);
    IGRAPH_ASSERT(igraph_cmp_epsilon(-IGRAPH_INFINITY, 2, 200) < 0);

    /* Infinities may be equal, though */
    IGRAPH_ASSERT(igraph_cmp_epsilon(IGRAPH_INFINITY, IGRAPH_INFINITY, 0.5) == 0);
    IGRAPH_ASSERT(igraph_cmp_epsilon(IGRAPH_INFINITY, IGRAPH_INFINITY, 0.1) == 0);
    IGRAPH_ASSERT(igraph_cmp_epsilon(IGRAPH_INFINITY, IGRAPH_INFINITY, 1e-7) == 0);
    IGRAPH_ASSERT(igraph_cmp_epsilon(IGRAPH_INFINITY, -IGRAPH_INFINITY, 0.5) > 0);
    IGRAPH_ASSERT(igraph_cmp_epsilon(IGRAPH_INFINITY, -IGRAPH_INFINITY, 0.1) > 0);
    IGRAPH_ASSERT(igraph_cmp_epsilon(IGRAPH_INFINITY, -IGRAPH_INFINITY, 1e-7) > 0);
    IGRAPH_ASSERT(igraph_cmp_epsilon(-IGRAPH_INFINITY, IGRAPH_INFINITY, 0.5) < 0);
    IGRAPH_ASSERT(igraph_cmp_epsilon(-IGRAPH_INFINITY, IGRAPH_INFINITY, 0.1) < 0);
    IGRAPH_ASSERT(igraph_cmp_epsilon(-IGRAPH_INFINITY, IGRAPH_INFINITY, 1e-7) < 0);

    /* NaNs are not equal to anything, not even each other. Sign of the result
     * is not guaranteed */
    IGRAPH_ASSERT(igraph_cmp_epsilon(IGRAPH_INFINITY, IGRAPH_NAN, 0.1) != 0);
    IGRAPH_ASSERT(igraph_cmp_epsilon(-IGRAPH_INFINITY, IGRAPH_NAN, 0.1) != 0);
    IGRAPH_ASSERT(igraph_cmp_epsilon(0, IGRAPH_NAN, 1e-7) != 0);
    IGRAPH_ASSERT(igraph_cmp_epsilon(IGRAPH_NAN, 0, 1e-7) != 0);
    IGRAPH_ASSERT(igraph_cmp_epsilon(IGRAPH_NAN, IGRAPH_NAN, 1e-7) != 0);

    /* Comparison with zero tolerance should be a "strict" comparison */
    IGRAPH_ASSERT(igraph_cmp_epsilon(1, 2, 0) < 0);
    IGRAPH_ASSERT(igraph_cmp_epsilon(2, 1, 0) > 0);
    IGRAPH_ASSERT(igraph_cmp_epsilon(2, 2, 0) == 0);
    IGRAPH_ASSERT(igraph_cmp_epsilon(0, DBL_MIN, 0) < 0);
    IGRAPH_ASSERT(igraph_cmp_epsilon(DBL_MIN, 0, 0) > 0);
    IGRAPH_ASSERT(igraph_cmp_epsilon(-IGRAPH_INFINITY, 0, 0) < 0);
    IGRAPH_ASSERT(igraph_cmp_epsilon(-IGRAPH_INFINITY, -IGRAPH_INFINITY, 0) == 0);
    IGRAPH_ASSERT(igraph_cmp_epsilon(IGRAPH_INFINITY, 0, 0) > 0);
    IGRAPH_ASSERT(igraph_cmp_epsilon(IGRAPH_INFINITY, IGRAPH_INFINITY, 0) == 0);
    IGRAPH_ASSERT(igraph_cmp_epsilon(IGRAPH_NAN, 0, 0) != 0);
    IGRAPH_ASSERT(igraph_cmp_epsilon(IGRAPH_NAN, IGRAPH_NAN, 0) != 0);
    IGRAPH_ASSERT(igraph_cmp_epsilon(IGRAPH_NAN, IGRAPH_INFINITY, 0) != 0);
    IGRAPH_ASSERT(igraph_cmp_epsilon(IGRAPH_NAN, -IGRAPH_INFINITY, 0) != 0);

    return 0;
}
