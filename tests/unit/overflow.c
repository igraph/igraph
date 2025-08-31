/*
   igraph library.
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
#include "math/safe_intop.h"

#define EXPECT_SUCC(expr) IGRAPH_ASSERT((expr) == IGRAPH_SUCCESS)
#define EXPECT_FAIL(expr) IGRAPH_ASSERT((expr) != IGRAPH_SUCCESS)

int main(void) {
    igraph_int_t res;

    igraph_set_error_handler(igraph_error_handler_printignore);

    /* Addition */

    EXPECT_SUCC(igraph_i_safe_add(1, 2, &res));
    EXPECT_SUCC(igraph_i_safe_add(-123, 42, &res));
    EXPECT_SUCC(igraph_i_safe_add(-5, -137, &res));
    EXPECT_SUCC(igraph_i_safe_add(IGRAPH_INTEGER_MAX-1, 1, &res));
    EXPECT_SUCC(igraph_i_safe_add(IGRAPH_INTEGER_MIN+1, -1, &res));
    EXPECT_SUCC(igraph_i_safe_add(IGRAPH_INTEGER_MIN, IGRAPH_INTEGER_MAX, &res));
    EXPECT_SUCC(igraph_i_safe_add(IGRAPH_INTEGER_MAX/2, IGRAPH_INTEGER_MAX/2, &res));

    EXPECT_FAIL(igraph_i_safe_add(IGRAPH_INTEGER_MAX, IGRAPH_INTEGER_MAX, &res));
    EXPECT_FAIL(igraph_i_safe_add(IGRAPH_INTEGER_MAX, 1, &res));
    EXPECT_FAIL(igraph_i_safe_add(IGRAPH_INTEGER_MIN, -1, &res));
    EXPECT_FAIL(igraph_i_safe_add(IGRAPH_INTEGER_MAX/2 + 2, IGRAPH_INTEGER_MAX/2, &res));

    /* Multiplication */

    EXPECT_SUCC(igraph_i_safe_mult(3, 4, &res));
    EXPECT_SUCC(igraph_i_safe_mult(-35, 14, &res));
    EXPECT_SUCC(igraph_i_safe_mult(-5, -9, &res));
    EXPECT_SUCC(igraph_i_safe_mult(IGRAPH_INTEGER_MAX, 1, &res));
    EXPECT_SUCC(igraph_i_safe_mult(1, IGRAPH_INTEGER_MAX, &res));
    EXPECT_SUCC(igraph_i_safe_mult(2, IGRAPH_INTEGER_MAX/2, &res));
    EXPECT_SUCC(igraph_i_safe_mult(2, IGRAPH_INTEGER_MIN/2, &res));

    EXPECT_FAIL(igraph_i_safe_mult(2, IGRAPH_INTEGER_MAX, &res));
    EXPECT_FAIL(igraph_i_safe_mult(2, IGRAPH_INTEGER_MIN, &res));
    EXPECT_FAIL(igraph_i_safe_mult(IGRAPH_INTEGER_MAX, IGRAPH_INTEGER_MAX, &res));
    EXPECT_FAIL(igraph_i_safe_mult(IGRAPH_INTEGER_MAX, IGRAPH_INTEGER_MIN, &res));
    EXPECT_FAIL(igraph_i_safe_mult(IGRAPH_INTEGER_MAX/2, 3, &res));
    EXPECT_FAIL(igraph_i_safe_mult(IGRAPH_INTEGER_MIN/2, 3, &res));
    EXPECT_FAIL(igraph_i_safe_mult(IGRAPH_INTEGER_MAX/2, -3, &res));
    EXPECT_FAIL(igraph_i_safe_mult(IGRAPH_INTEGER_MIN/2, -3, &res));

    /* The following tests assume that, mathematically, INTEGER_MIN = -INTEGER_MAX - 1 */

    EXPECT_SUCC(igraph_i_safe_mult(IGRAPH_INTEGER_MAX, -1, &res));
    EXPECT_SUCC(igraph_i_safe_mult(-1, IGRAPH_INTEGER_MAX, &res));

    EXPECT_FAIL(igraph_i_safe_mult(-1, IGRAPH_INTEGER_MIN, &res));
    EXPECT_FAIL(igraph_i_safe_mult(IGRAPH_INTEGER_MIN, -1, &res));

    return 0;
}
