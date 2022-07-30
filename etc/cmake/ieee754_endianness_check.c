/* Checks whether the endianness of IEEE754 doubles matches the endianness of
 * uint64_t on the target system. This is needed to ensure that the trick we
 * employ in igraph_rng_get_unif01() works */

#include <stdint.h>
#include <stdio.h>

union {
    uint64_t as_uint64_t;
    double as_double;
} value;

int main() {
    value.as_uint64_t = 4841376218035192321ULL;
    if (value.as_double == 4510218239279617.0) {
        return 0; /* endianness of uint64_t and double match */
    } else {
        return 1; /* endianness of uint64_t and double do not match */
    }
}
