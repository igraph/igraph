/* -*- mode: C -*-  */
/* vim:set ts=4 sw=4 sts=4 et: */
/*
   IGraph library.
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
#include <stdlib.h>

#include "test_utilities.h"

int main(void) {
    igraph_vector_int_t v, u, w;
    igraph_integer_t i, n;

    printf("Initialise empty bitset\n");
    n = 0;
    igraph_bitset_init(&v, n);
    igraph_bitset_destroy(&v);
    printf("\n");

    printf("Initialise bitset of length 32\n");
    n = 32;
    igraph_bitset_init(&v, n);
    print_bitset(&v, n);
    igraph_bitset_destroy(&v);
    printf("\n");

    printf("Initialise bitset of length 64\n");
    n = 64;
    igraph_bitset_init(&v, n);
    print_bitset(&v, n);
    igraph_bitset_destroy(&v);
    printf("\n");

    printf("Initialise bitset of length 75\n");
    n = 75;
    igraph_bitset_init(&v, n);
    print_bitset(&v, n);
    igraph_bitset_destroy(&v);
    printf("\n");

    printf("Test IGRAPH_BITSET\n");
    n = 35;
    igraph_bitset_init(&v, n);
    IGRAPH_BITSET(VECTOR(v), 0);
    print_bitset(&v, n);
    IGRAPH_BITSET(VECTOR(v), 24);
    print_bitset(&v, n);
    IGRAPH_BITSET(VECTOR(v), 13);
    print_bitset(&v, n);
    IGRAPH_BITSET(VECTOR(v), 17);
    print_bitset(&v, n);
    IGRAPH_BITSET(VECTOR(v), 34);
    print_bitset(&v, n);
    IGRAPH_BITSET(VECTOR(v), 13);
    print_bitset(&v, n);
    printf("\n");

    printf("Test IGRAPH_BITCLEAR\n");
    IGRAPH_BITCLEAR(VECTOR(v), 33);
    print_bitset(&v, n);
    IGRAPH_BITCLEAR(VECTOR(v), 34);
    print_bitset(&v, n);
    IGRAPH_BITCLEAR(VECTOR(v), 17);
    print_bitset(&v, n);
    igraph_bitset_destroy(&v);
    printf("\n");

    printf("Test OR\n");
    n = 7;
    igraph_bitset_init(&v, n);
    igraph_bitset_init(&u, n);
    igraph_bitset_init(&w, n);
    IGRAPH_BITSET(VECTOR(v), 0);
    print_bitset(&v, n);
    IGRAPH_BITSET(VECTOR(u), 1);
    IGRAPH_BITSET(VECTOR(u), 2);
    IGRAPH_BITSET(VECTOR(u), 5);
    print_bitset(&u, n);
    IGRAPH_BITSET(VECTOR(w), 2);
    IGRAPH_BITSET(VECTOR(w), 3);
    IGRAPH_BITSET(VECTOR(w), 6);
    print_bitset(&w, n);
    igraph_bitset_or(&v, &u, &w, n);
    print_bitset(&v, n);
    igraph_bitset_destroy(&v);
    printf("\n");

    printf("Test AND\n");
    n = 7;
    igraph_bitset_init(&v, n);
    igraph_bitset_init(&u, n);
    igraph_bitset_init(&w, n);
    IGRAPH_BITSET(VECTOR(v), 0);
    print_bitset(&v, n);
    IGRAPH_BITSET(VECTOR(u), 1);
    IGRAPH_BITSET(VECTOR(u), 2);
    IGRAPH_BITSET(VECTOR(u), 5);
    print_bitset(&u, n);
    IGRAPH_BITSET(VECTOR(w), 2);
    IGRAPH_BITSET(VECTOR(w), 3);
    IGRAPH_BITSET(VECTOR(w), 6);
    print_bitset(&w, n);
    igraph_bitset_and(&v, &u, &w, n);
    print_bitset(&v, n);
    igraph_bitset_destroy(&v);
    printf("\n");

    printf("Test XOR\n");
    n = 7;
    igraph_bitset_init(&v, n);
    igraph_bitset_init(&u, n);
    igraph_bitset_init(&w, n);
    IGRAPH_BITSET(VECTOR(v), 0);
    print_bitset(&v, n);
    IGRAPH_BITSET(VECTOR(u), 1);
    IGRAPH_BITSET(VECTOR(u), 2);
    IGRAPH_BITSET(VECTOR(u), 5);
    print_bitset(&u, n);
    IGRAPH_BITSET(VECTOR(w), 2);
    IGRAPH_BITSET(VECTOR(w), 3);
    IGRAPH_BITSET(VECTOR(w), 6);
    print_bitset(&w, n);
    igraph_bitset_xor(&v, &u, &w, n);
    print_bitset(&v, n);
    igraph_bitset_destroy(&v);
    printf("\n");

    printf("Test NOT\n");
    n = 7;
    igraph_bitset_init(&v, n);
    igraph_bitset_init(&u, n);
    IGRAPH_BITSET(VECTOR(v), 1);
    print_bitset(&v, n);
    IGRAPH_BITSET(VECTOR(u), 1);
    IGRAPH_BITSET(VECTOR(u), 2);
    IGRAPH_BITSET(VECTOR(u), 5);
    print_bitset(&u, n);
    igraph_bitset_not(&v, &u, n);
    print_bitset(&v, n);
    igraph_bitset_destroy(&v);
    printf("\n");

    printf("Test popcount\n");
    n = 75;
    igraph_bitset_init(&v, n);
    print_bitset(&v, n);
    printf("Popcount: %ld\n", igraph_bitset_popcount(&v, n));
    for (i = 2; i < n; i++) {
        if (i == 2 || i == 3 || i == 5 || i == 7 || i == 11 || !(i%2 == 0 || i%3 == 0 || i%5 == 0 || i%7 == 0 || i%11 == 0)) {
            IGRAPH_BITSET(VECTOR(v), i);
        }
    }
    print_bitset(&v, n);
    printf("Popcount: %ld\n", igraph_bitset_popcount(&v, n));
    IGRAPH_BITCLEAR(VECTOR(v), 67);
    print_bitset(&v, n);
    printf("Popcount: %ld\n", igraph_bitset_popcount(&v, n));
    IGRAPH_BITSET(VECTOR(v), 9);
    print_bitset(&v, n);
    printf("Popcount: %ld\n", igraph_bitset_popcount(&v, n));
    igraph_bitset_not(&v, &v, n);
    print_bitset(&v, n);
    printf("Popcount: %ld\n", igraph_bitset_popcount(&v, n));
    igraph_bitset_destroy(&v);
    printf("\n");

    printf("Test count leading zeros\n");
    n = 67;
    igraph_bitset_init(&v, n);
    igraph_bitset_not(&v, &v, n);
    print_bitset(&v, n);
    printf("Leading zeros: %ld\n", igraph_bitset_countl_zero(&v, n));
    igraph_bitset_not(&v, &v, n);
    print_bitset(&v, n);
    printf("Leading zeros: %ld\n", igraph_bitset_countl_zero(&v, n));
    IGRAPH_BITSET(VECTOR(v), 0);
    print_bitset(&v, n);
    printf("Leading zeros: %ld\n", igraph_bitset_countl_zero(&v, n));
    IGRAPH_BITSET(VECTOR(v), 23);
    print_bitset(&v, n);
    printf("Leading zeros: %ld\n", igraph_bitset_countl_zero(&v, n));
    IGRAPH_BITCLEAR(VECTOR(v), 0);
    print_bitset(&v, n);
    printf("Leading zeros: %ld\n", igraph_bitset_countl_zero(&v, n));
    IGRAPH_BITSET(VECTOR(v), n - 3);
    print_bitset(&v, n);
    printf("Leading zeros: %ld\n", igraph_bitset_countl_zero(&v, n));
    IGRAPH_BITSET(VECTOR(v), n - 1);
    print_bitset(&v, n);
    printf("Leading zeros: %ld\n", igraph_bitset_countl_zero(&v, n));
    igraph_bitset_destroy(&v);
    printf("\n");


    printf("Test count leading ones\n");
    n = 67;
    igraph_bitset_init(&v, n);
    igraph_bitset_not(&v, &v, n);
    print_bitset(&v, n);
    printf("Leading ones: %ld\n", igraph_bitset_countl_one(&v, n));
    IGRAPH_BITCLEAR(VECTOR(v), 0);
    print_bitset(&v, n);
    printf("Leading ones: %ld\n", igraph_bitset_countl_one(&v, n));
    IGRAPH_BITCLEAR(VECTOR(v), 23);
    print_bitset(&v, n);
    printf("Leading ones: %ld\n", igraph_bitset_countl_one(&v, n));
    IGRAPH_BITSET(VECTOR(v), 0);
    print_bitset(&v, n);
    printf("Leading ones: %ld\n", igraph_bitset_countl_one(&v, n));
    IGRAPH_BITCLEAR(VECTOR(v), n - 3);
    print_bitset(&v, n);
    printf("Leading ones: %ld\n", igraph_bitset_countl_one(&v, n));
    IGRAPH_BITCLEAR(VECTOR(v), n - 1);
    print_bitset(&v, n);
    printf("Leading ones: %ld\n", igraph_bitset_countl_one(&v, n));
    igraph_bitset_destroy(&v);
    printf("\n");


    printf("Test count trailing zeros\n");
    n = 67;
    igraph_bitset_init(&v, n);
    igraph_bitset_not(&v, &v, n);
    print_vector_int(&v);
    print_bitset(&v, n);
    printf("Trailing zeros: %ld\n", igraph_bitset_countr_zero(&v, n));
    igraph_bitset_not(&v, &v, n);
    print_vector_int(&v);
    print_bitset(&v, n);
    printf("Trailing zeros: %ld\n", igraph_bitset_countr_zero(&v, n));
    IGRAPH_BITSET(VECTOR(v), 0);
    print_vector_int(&v);
    print_bitset(&v, n);
    printf("Trailing zeros: %ld\n", igraph_bitset_countr_zero(&v, n));
    IGRAPH_BITSET(VECTOR(v), 23);
    print_bitset(&v, n);
    printf("Trailing zeros: %ld\n", igraph_bitset_countr_zero(&v, n));
    IGRAPH_BITCLEAR(VECTOR(v), 0);
    print_bitset(&v, n);
    printf("Trailing zeros: %ld\n", igraph_bitset_countr_zero(&v, n));
    IGRAPH_BITCLEAR(VECTOR(v), 23);
    IGRAPH_BITSET(VECTOR(v), n - 1);
    print_bitset(&v, n);
    printf("Trailing zeros: %ld\n", igraph_bitset_countr_zero(&v, n));
    IGRAPH_BITSET(VECTOR(v), 2);
    print_bitset(&v, n);
    printf("Trailing zeros: %ld\n", igraph_bitset_countr_zero(&v, n));
    igraph_bitset_destroy(&v);
    printf("\n");

    printf("Test count trailing ones\n");
    n = 67;
    igraph_bitset_init(&v, n);
    igraph_bitset_not(&v, &v, n);
    print_bitset(&v, n);
    printf("Trailing ones: %ld\n", igraph_bitset_countr_one(&v, n));
    IGRAPH_BITCLEAR(VECTOR(v), 23);
    print_bitset(&v, n);
    printf("Trailing ones: %ld\n", igraph_bitset_countr_one(&v, n));
    IGRAPH_BITCLEAR(VECTOR(v), 0);
    print_bitset(&v, n);
    printf("Trailing ones: %ld\n", igraph_bitset_countr_one(&v, n));
    IGRAPH_BITSET(VECTOR(v), 0);
    print_bitset(&v, n);
    printf("Trailing ones: %ld\n", igraph_bitset_countr_one(&v, n));
    IGRAPH_BITSET(VECTOR(v), 23);
    IGRAPH_BITCLEAR(VECTOR(v), n - 1);
    print_bitset(&v, n);
    printf("Trailing ones: %ld\n", igraph_bitset_countr_one(&v, n));
    IGRAPH_BITCLEAR(VECTOR(v), 2);
    print_bitset(&v, n);
    printf("Trailing ones: %ld\n", igraph_bitset_countr_one(&v, n));
    igraph_bitset_destroy(&v);
    printf("\n");

    VERIFY_FINALLY_STACK();

    return 0;
}
