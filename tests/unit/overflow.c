
#include <igraph.h>
#include "math/safe_intop.h"

#define EXPECT_SUCC(expr) IGRAPH_ASSERT((expr) == IGRAPH_SUCCESS)
#define EXPECT_FAIL(expr) IGRAPH_ASSERT((expr) != IGRAPH_SUCCESS)

int main() {
    igraph_integer_t res;

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
