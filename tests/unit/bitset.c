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

#include "test_utilities.h"

int main(void) {
    igraph_bitset_t v1, v2, v3;
    igraph_int_t n;

    printf("Initialise empty bitset\n");
    n = 0;
    igraph_bitset_init(&v1, n);
    IGRAPH_ASSERT(igraph_bitset_size(&v1) == n);
    igraph_bitset_destroy(&v1);
    printf("\n");

    printf("Initialise bitset of length 32\n");
    n = 32;
    igraph_bitset_init(&v1, n);
    print_bitset(&v1);
    igraph_bitset_destroy(&v1);
    printf("\n");

    printf("Initialise bitset of length 64\n");
    n = 64;
    igraph_bitset_init(&v1, n);
    print_bitset(&v1);
    igraph_bitset_destroy(&v1);
    printf("\n");

    printf("Initialise bitset of length 75\n");
    n = 75;
    igraph_bitset_init(&v1, n);
    print_bitset(&v1);
    IGRAPH_ASSERT(igraph_bitset_size(&v1) == n);
    igraph_bitset_destroy(&v1);
    printf("\n");

    printf("Test IGRAPH_BIT_SET\n");
    n = 35;
    igraph_bitset_init(&v1, n);
    IGRAPH_BIT_SET(v1, 0);
    print_bitset(&v1);
    IGRAPH_BIT_SET(v1, 24);
    print_bitset(&v1);
    IGRAPH_BIT_SET(v1, 13);
    print_bitset(&v1);
    IGRAPH_BIT_SET(v1, 17);
    print_bitset(&v1);
    IGRAPH_BIT_SET(v1, 34);
    print_bitset(&v1);
    IGRAPH_BIT_SET(v1, 13);
    print_bitset(&v1);
    printf("\n");

    printf("Test IGRAPH_BIT_CLEAR\n");
    IGRAPH_BIT_CLEAR(v1, 33);
    print_bitset(&v1);
    IGRAPH_BIT_CLEAR(v1, 34);
    print_bitset(&v1);
    IGRAPH_BIT_CLEAR(v1, 17);
    print_bitset(&v1);
    printf("\n");

    printf("Test printing\n");
    igraph_bitset_print(&v1);
    printf("\n\n");

    printf("Test bitset copy constructor\n");
    igraph_bitset_init_copy(&v2, &v1);
    print_bitset(&v2);
    igraph_bitset_destroy(&v1);
    igraph_bitset_destroy(&v2);
    printf("\n");

    printf("Test bitset resize\n");
    igraph_bitset_init(&v1, 0);
    print_bitset(&v1);
    IGRAPH_ASSERT(igraph_bitset_size(&v1) == 0);
    igraph_bitset_resize(&v1, 10);
    print_bitset(&v1);
    IGRAPH_ASSERT(igraph_bitset_size(&v1) == 10);
    IGRAPH_BIT_SET(v1, 3);
    print_bitset(&v1);
    igraph_bitset_resize(&v1, 64);
    print_bitset(&v1);
    IGRAPH_BIT_SET(v1, 63);
    print_bitset(&v1);
    igraph_bitset_resize(&v1, 70);
    print_bitset(&v1);
    IGRAPH_BIT_SET(v1, 64);
    print_bitset(&v1);
    igraph_bitset_resize(&v1, 64);
    print_bitset(&v1);
    igraph_bitset_resize(&v1, 63);
    print_bitset(&v1);
    IGRAPH_ASSERT(igraph_bitset_size(&v1) == 63);
    igraph_bitset_destroy(&v1);
    printf("\n");

    printf("Test OR\n");
    n = 7;
    igraph_bitset_init(&v1, n);
    igraph_bitset_init(&v2, n);
    igraph_bitset_init(&v3, n);
    IGRAPH_BIT_SET(v1, 0);
    print_bitset(&v1);
    IGRAPH_BIT_SET(v2, 1);
    IGRAPH_BIT_SET(v2, 2);
    IGRAPH_BIT_SET(v2, 5);
    print_bitset(&v2);
    IGRAPH_BIT_SET(v3, 2);
    IGRAPH_BIT_SET(v3, 3);
    IGRAPH_BIT_SET(v3, 6);
    print_bitset(&v3);
    igraph_bitset_or(&v1, &v2, &v3);
    print_bitset(&v1);
    igraph_bitset_destroy(&v1);
    igraph_bitset_destroy(&v2);
    igraph_bitset_destroy(&v3);
    printf("\n");

    printf("Test AND\n");
    n = 7;
    igraph_bitset_init(&v1, n);
    igraph_bitset_init(&v2, n);
    igraph_bitset_init(&v3, n);
    IGRAPH_BIT_SET(v1, 0);
    print_bitset(&v1);
    IGRAPH_BIT_SET(v2, 1);
    IGRAPH_BIT_SET(v2, 2);
    IGRAPH_BIT_SET(v2, 5);
    print_bitset(&v2);
    IGRAPH_BIT_SET(v3, 2);
    IGRAPH_BIT_SET(v3, 3);
    IGRAPH_BIT_SET(v3, 6);
    print_bitset(&v3);
    igraph_bitset_and(&v1, &v2, &v3);
    print_bitset(&v1);
    igraph_bitset_destroy(&v1);
    igraph_bitset_destroy(&v2);
    igraph_bitset_destroy(&v3);
    printf("\n");

    printf("Test XOR\n");
    n = 7;
    igraph_bitset_init(&v2, n);
    igraph_bitset_init(&v1, n);
    igraph_bitset_init(&v3, n);
    IGRAPH_BIT_SET(v1, 0);
    print_bitset(&v1);
    IGRAPH_BIT_SET(v2, 1);
    IGRAPH_BIT_SET(v2, 2);
    IGRAPH_BIT_SET(v2, 5);
    print_bitset(&v2);
    IGRAPH_BIT_SET(v3, 2);
    IGRAPH_BIT_SET(v3, 3);
    IGRAPH_BIT_SET(v3, 6);
    print_bitset(&v3);
    igraph_bitset_xor(&v1, &v2, &v3);
    print_bitset(&v1);
    igraph_bitset_destroy(&v1);
    igraph_bitset_destroy(&v2);
    igraph_bitset_destroy(&v3);
    printf("\n");

    printf("Test NOT\n");
    n = 7;
    igraph_bitset_init(&v1, n);
    igraph_bitset_init(&v2, n);
    IGRAPH_BIT_SET(v1, 1);
    print_bitset(&v1);
    IGRAPH_BIT_SET(v2, 1);
    IGRAPH_BIT_SET(v2, 2);
    IGRAPH_BIT_SET(v2, 5);
    print_bitset(&v2);
    igraph_bitset_not(&v1, &v2);
    print_bitset(&v1);
    igraph_bitset_destroy(&v1);
    igraph_bitset_destroy(&v2);
    printf("\n");

    printf("Test popcount\n");
    n = 75;
    igraph_bitset_init(&v1, n);
    print_bitset(&v1);

    printf("Popcount: %" IGRAPH_PRId "\n", igraph_bitset_popcount(&v1));
    for (igraph_int_t i = 2; i < n; i++) {
        if (i == 2 || i == 3 || i == 5 || i == 7 || i == 11 || !(i%2 == 0 || i%3 == 0 || i%5 == 0 || i%7 == 0 || i%11 == 0)) {
            IGRAPH_BIT_SET(v1, i);
        }
    }
    print_bitset(&v1);
    printf("Popcount: %" IGRAPH_PRId "\n", igraph_bitset_popcount(&v1));
    IGRAPH_BIT_CLEAR(v1, 67);
    print_bitset(&v1);
    printf("Popcount: %" IGRAPH_PRId "\n", igraph_bitset_popcount(&v1));
    IGRAPH_BIT_SET(v1, 9);
    print_bitset(&v1);
    printf("Popcount: %" IGRAPH_PRId "\n", igraph_bitset_popcount(&v1));
    igraph_bitset_not(&v1, &v1);
    print_bitset(&v1);
    printf("Popcount: %" IGRAPH_PRId "\n", igraph_bitset_popcount(&v1));
    igraph_bitset_destroy(&v1);
    printf("\n");

    printf("Test count leading zeros\n");
    n = 67;
    igraph_bitset_init(&v1, n);
    igraph_bitset_not(&v1, &v1);
    print_bitset(&v1);
    printf("Leading zeros: %" IGRAPH_PRId "\n", igraph_bitset_countl_zero(&v1));
    igraph_bitset_not(&v1, &v1);
    print_bitset(&v1);
    printf("Leading zeros: %" IGRAPH_PRId "\n", igraph_bitset_countl_zero(&v1));
    IGRAPH_BIT_SET(v1, 0);
    print_bitset(&v1);
    printf("Leading zeros: %" IGRAPH_PRId "\n", igraph_bitset_countl_zero(&v1));
    IGRAPH_BIT_SET(v1, 23);
    print_bitset(&v1);
    printf("Leading zeros: %" IGRAPH_PRId "\n", igraph_bitset_countl_zero(&v1));
    IGRAPH_BIT_CLEAR(v1, 0);
    print_bitset(&v1);
    printf("Leading zeros: %" IGRAPH_PRId "\n", igraph_bitset_countl_zero(&v1));
    IGRAPH_BIT_SET(v1, n - 3);
    print_bitset(&v1);
    printf("Leading zeros: %" IGRAPH_PRId "\n", igraph_bitset_countl_zero(&v1));
    IGRAPH_BIT_SET(v1, n - 1);
    print_bitset(&v1);
    printf("Leading zeros: %" IGRAPH_PRId "\n", igraph_bitset_countl_zero(&v1));
    igraph_bitset_destroy(&v1);
    printf("\n");


    printf("Test count leading ones\n");
    n = 67;
    igraph_bitset_init(&v1, n);
    igraph_bitset_not(&v1, &v1);
    print_bitset(&v1);
    printf("Leading ones: %" IGRAPH_PRId "\n", igraph_bitset_countl_one(&v1));
    IGRAPH_BIT_CLEAR(v1, 0);
    print_bitset(&v1);
    printf("Leading ones: %" IGRAPH_PRId "\n", igraph_bitset_countl_one(&v1));
    IGRAPH_BIT_CLEAR(v1, 23);
    print_bitset(&v1);
    printf("Leading ones: %" IGRAPH_PRId "\n", igraph_bitset_countl_one(&v1));
    IGRAPH_BIT_SET(v1, 0);
    print_bitset(&v1);
    printf("Leading ones: %" IGRAPH_PRId "\n", igraph_bitset_countl_one(&v1));
    IGRAPH_BIT_CLEAR(v1, n - 3);
    print_bitset(&v1);
    printf("Leading ones: %" IGRAPH_PRId "\n", igraph_bitset_countl_one(&v1));
    IGRAPH_BIT_CLEAR(v1, n - 1);
    print_bitset(&v1);
    printf("Leading ones: %" IGRAPH_PRId "\n", igraph_bitset_countl_one(&v1));
    igraph_bitset_destroy(&v1);
    printf("\n");


    printf("Test count trailing zeros\n");
    n = 67;
    igraph_bitset_init(&v1, n);
    igraph_bitset_not(&v1, &v1);
    print_bitset(&v1);
    printf("Trailing zeros: %" IGRAPH_PRId "\n", igraph_bitset_countr_zero(&v1));
    igraph_bitset_not(&v1, &v1);
    print_bitset(&v1);
    printf("Trailing zeros: %" IGRAPH_PRId "\n", igraph_bitset_countr_zero(&v1));
    IGRAPH_BIT_SET(v1, 0);
    print_bitset(&v1);
    printf("Trailing zeros: %" IGRAPH_PRId "\n", igraph_bitset_countr_zero(&v1));
    IGRAPH_BIT_SET(v1, 23);
    print_bitset(&v1);
    printf("Trailing zeros: %" IGRAPH_PRId "\n", igraph_bitset_countr_zero(&v1));
    IGRAPH_BIT_CLEAR(v1, 0);
    print_bitset(&v1);
    printf("Trailing zeros: %" IGRAPH_PRId "\n", igraph_bitset_countr_zero(&v1));
    IGRAPH_BIT_CLEAR(v1, 23);
    IGRAPH_BIT_SET(v1, n - 1);
    print_bitset(&v1);
    printf("Trailing zeros: %" IGRAPH_PRId "\n", igraph_bitset_countr_zero(&v1));
    IGRAPH_BIT_SET(v1, 2);
    print_bitset(&v1);
    printf("Trailing zeros: %" IGRAPH_PRId "\n", igraph_bitset_countr_zero(&v1));
    igraph_bitset_destroy(&v1);
    printf("\n");

    printf("Test count trailing ones\n");
    n = 67;
    igraph_bitset_init(&v1, n);
    igraph_bitset_not(&v1, &v1);
    print_bitset(&v1);
    printf("Trailing ones: %" IGRAPH_PRId "\n", igraph_bitset_countr_one(&v1));
    IGRAPH_BIT_CLEAR(v1, 23);
    print_bitset(&v1);
    printf("Trailing ones: %" IGRAPH_PRId "\n", igraph_bitset_countr_one(&v1));
    IGRAPH_BIT_CLEAR(v1, 0);
    print_bitset(&v1);
    printf("Trailing ones: %" IGRAPH_PRId "\n", igraph_bitset_countr_one(&v1));
    IGRAPH_BIT_SET(v1, 0);
    print_bitset(&v1);
    printf("Trailing ones: %" IGRAPH_PRId "\n", igraph_bitset_countr_one(&v1));
    IGRAPH_BIT_SET(v1, 23);
    IGRAPH_BIT_CLEAR(v1, n - 1);
    print_bitset(&v1);
    printf("Trailing ones: %" IGRAPH_PRId "\n", igraph_bitset_countr_one(&v1));
    IGRAPH_BIT_CLEAR(v1, 2);
    print_bitset(&v1);
    printf("Trailing ones: %" IGRAPH_PRId "\n", igraph_bitset_countr_one(&v1));
    igraph_bitset_destroy(&v1);
    printf("\n");

    printf("Test checking all elements\n");
    /* 0000 */
    igraph_bitset_init(&v1, 4);
    IGRAPH_ASSERT(! igraph_bitset_is_all_one(&v1));
    IGRAPH_ASSERT(! igraph_bitset_is_any_one(&v1));
    IGRAPH_ASSERT(igraph_bitset_is_all_zero(&v1));
    IGRAPH_ASSERT(igraph_bitset_is_any_zero(&v1));
    /* 1000 */
    IGRAPH_BIT_SET(v1, 3);
    IGRAPH_ASSERT(! igraph_bitset_is_all_one(&v1));
    IGRAPH_ASSERT(igraph_bitset_is_any_one(&v1));
    IGRAPH_ASSERT(! igraph_bitset_is_all_zero(&v1));
    IGRAPH_ASSERT(igraph_bitset_is_any_zero(&v1));
    /* 000 */
    /* This resize leaves a set bit in the unused part of the last word of the.
     * bitset. This helps test that masking out the unused part is working. */
    igraph_bitset_resize(&v1, 3);
    IGRAPH_ASSERT(! igraph_bitset_is_all_one(&v1));
    IGRAPH_ASSERT(! igraph_bitset_is_any_one(&v1));
    IGRAPH_ASSERT(igraph_bitset_is_all_zero(&v1));
    IGRAPH_ASSERT(igraph_bitset_is_any_zero(&v1));
    /* 001 */
    IGRAPH_BIT_SET(v1, 0);
    IGRAPH_ASSERT(! igraph_bitset_is_all_one(&v1));
    IGRAPH_ASSERT(igraph_bitset_is_any_one(&v1));
    IGRAPH_ASSERT(! igraph_bitset_is_all_zero(&v1));
    IGRAPH_ASSERT(igraph_bitset_is_any_zero(&v1));
    /* 111 */
    IGRAPH_BIT_SET(v1, 1); IGRAPH_BIT_SET(v1, 2);
    IGRAPH_ASSERT(igraph_bitset_is_all_one(&v1));
    IGRAPH_ASSERT(igraph_bitset_is_any_one(&v1));
    IGRAPH_ASSERT(! igraph_bitset_is_all_zero(&v1));
    IGRAPH_ASSERT(! igraph_bitset_is_any_zero(&v1));
    igraph_bitset_destroy(&v1);

    igraph_bitset_init(&v1, 67);
    /* all zeros */
    IGRAPH_ASSERT(! igraph_bitset_is_all_one(&v1));
    IGRAPH_ASSERT(! igraph_bitset_is_any_one(&v1));
    IGRAPH_ASSERT(igraph_bitset_is_all_zero(&v1));
    IGRAPH_ASSERT(igraph_bitset_is_any_zero(&v1));
    /* some ones */
    IGRAPH_BIT_SET(v1, 10);
    IGRAPH_ASSERT(! igraph_bitset_is_all_one(&v1));
    IGRAPH_ASSERT(igraph_bitset_is_any_one(&v1));
    IGRAPH_ASSERT(! igraph_bitset_is_all_zero(&v1));
    IGRAPH_ASSERT(igraph_bitset_is_any_zero(&v1));
    /* all ones */
    igraph_bitset_fill(&v1, true);
    IGRAPH_ASSERT(igraph_bitset_is_all_one(&v1));
    IGRAPH_ASSERT(igraph_bitset_is_any_one(&v1));
    IGRAPH_ASSERT(! igraph_bitset_is_all_zero(&v1));
    IGRAPH_ASSERT(! igraph_bitset_is_any_zero(&v1));
    /* some zeros */
    IGRAPH_BIT_CLEAR(v1, 20);
    IGRAPH_ASSERT(! igraph_bitset_is_all_one(&v1));
    IGRAPH_ASSERT(igraph_bitset_is_any_one(&v1));
    IGRAPH_ASSERT(! igraph_bitset_is_all_zero(&v1));
    IGRAPH_ASSERT(igraph_bitset_is_any_zero(&v1));
    igraph_bitset_destroy(&v1);

    /* Empty bitset */
    igraph_bitset_init(&v1, 0);
    IGRAPH_ASSERT(igraph_bitset_is_all_one(&v1));
    IGRAPH_ASSERT(! igraph_bitset_is_any_one(&v1));
    IGRAPH_ASSERT(igraph_bitset_is_all_zero(&v1));
    IGRAPH_ASSERT(! igraph_bitset_is_any_zero(&v1));
    igraph_bitset_destroy(&v1);

    VERIFY_FINALLY_STACK();

    return 0;
}
