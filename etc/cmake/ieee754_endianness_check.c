/* Checks whether the endianness of IEEE754 doubles matches the endianness of
 * uint64_t on the target system. This is needed to ensure that the trick we
 * employ in igraph_rng_get_unif01() works. */

#include <stdint.h>
#include <stdio.h>

union {
    uint64_t as_uint64_t;
    double as_double;
} value;

int main(void) {
    value.as_uint64_t = 4841376218035192321ULL;
    if (value.as_double == 4510218239279617.0) {
        /* endianness of uint64_t and double match */
        printf("OK\n");
    }
    /* We always return 0, even for a negative result, this is because we
     * need to tell on the CMake side whether a compiler misconfiguration
     * aborted our program, which can then be detected from a nonzero exit
     * code. */
    return 0;
}
