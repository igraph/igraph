
#ifndef IGRAPH_MATH_SAFE_INTOP_H
#define IGRAPH_MATH_SAFE_INTOP_H

#include "igraph_types.h"

/* These macros raise an error if the operation would result in an overflow.
 * They must only be used in functions that return an igraph_error_t.
 *
 * This code is based on the recommendation of
 * https://wiki.sei.cmu.edu/confluence/display/c/SEI+CERT+C+Coding+Standard
 */

#define IGRAPH_SAFE_ADD(ia, ib, res) \
    do { \
        igraph_integer_t a = (ia), b = (ib); \
        igraph_integer_t sum; \
        if (((b > 0) && (a > (IGRAPH_INTEGER_MAX - b))) || \
            ((b < 0) && (a < (IGRAPH_INTEGER_MIN - b)))) { \
            IGRAPH_ERRORF("Overflow when adding %"IGRAPH_PRId" and %"IGRAPH_PRId".", IGRAPH_EOVERFLOW, a, b); \
        } \
        sum = a+b; \
        *(res) = sum; \
    } while (0)

#define IGRAPH_SAFE_MULT(ia, ib, res) \
    do { \
        igraph_integer_t a = (ia), b = (ib); \
        igraph_integer_t prod; \
        int err=0; \
        if (a > 0) {  /* a is positive */ \
            if (b > 0) {  /* a and b are positive */ \
                if (a > (IGRAPH_INTEGER_MAX / b)) { \
                    err=1; \
                } \
            } else { /* a positive, b nonpositive */ \
                if (b < (IGRAPH_INTEGER_MIN / a)) { \
                    err=1; \
                } \
            } /* a positive, b nonpositive */ \
        } else { /* a is nonpositive */ \
            if (b > 0) { /* a is nonpositive, b is positive */ \
                if (a < (IGRAPH_INTEGER_MIN / b)) { \
                    err=1; \
                } \
            } else { /* a and b are nonpositive */ \
                if ( (a != 0) && (b < (IGRAPH_INTEGER_MAX / a))) { \
                    err=1; \
                } \
            } /* End if a and b are nonpositive */ \
        } /* End if a is nonpositive */ \
        if (err) { \
            IGRAPH_ERRORF("Overflow when multiplying %"IGRAPH_PRId" and %"IGRAPH_PRId".", IGRAPH_EOVERFLOW, a, b); \
        } \
        prod = a*b; \
        *(res) = prod; \
    } while (0)

#endif
