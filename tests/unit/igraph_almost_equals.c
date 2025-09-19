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

#include <float.h>

#include "test_utilities.h"

/* This file is ported from Java; the original source is here:
   https://floating-point-gui.de/errors/NearlyEqualsTest.java */

const double EPS = 0.00001;

void assert_almost_equal_with_eps(double a, double b, double eps, int line) {
    if (!igraph_almost_equals(a, b, eps)) {
        igraph_fatalf("Assertion failed: %g == %g with eps = %g", IGRAPH_FILE_BASENAME, line, a, b, eps);
    }
}

void assert_not_equal_with_eps(double a, double b, double eps, int line) {
    if (igraph_almost_equals(a, b, eps)) {
        igraph_fatalf("Assertion failed: %g != %g with eps = %g", IGRAPH_FILE_BASENAME, line, a, b, eps);
    }
}

#define ASSERT_ALMOST_EQUAL_WITH_EPS(a, b, eps) assert_almost_equal_with_eps(a, b, eps, __LINE__)
#define ASSERT_ALMOST_EQUAL(a, b) assert_almost_equal_with_eps(a, b, EPS, __LINE__)
#define ASSERT_NOT_EQUAL_WITH_EPS(a, b, eps) assert_not_equal_with_eps(a, b, eps, __LINE__)
#define ASSERT_NOT_EQUAL(a, b) assert_not_equal_with_eps(a, b, EPS, __LINE__)

void test_large_numbers(void) {
    ASSERT_ALMOST_EQUAL(1000000, 1000001);
    ASSERT_ALMOST_EQUAL(1000001, 1000000);
    ASSERT_NOT_EQUAL(10000, 10001);
    ASSERT_NOT_EQUAL(10001, 10000);
}

void test_large_negative_numbers(void) {
    ASSERT_ALMOST_EQUAL(-1000000, -1000001);
    ASSERT_ALMOST_EQUAL(-1000001, -1000000);
    ASSERT_NOT_EQUAL(-10000, -10001);
    ASSERT_NOT_EQUAL(-10001, -10000);
}

void test_numbers_around_one(void) {
    ASSERT_ALMOST_EQUAL(1.0000001, 1.0000002);
    ASSERT_ALMOST_EQUAL(1.0000002, 1.0000001);
    ASSERT_NOT_EQUAL(1.0002, 1.0001);
    ASSERT_NOT_EQUAL(1.0001, 1.0002);
}

void test_numbers_around_minus_one(void) {
    ASSERT_ALMOST_EQUAL(-1.0000001, -1.0000002);
    ASSERT_ALMOST_EQUAL(-1.0000002, -1.0000001);
    ASSERT_NOT_EQUAL(-1.0002, -1.0001);
    ASSERT_NOT_EQUAL(-1.0001, -1.0002);
}

void test_small_numbers(void) {
    ASSERT_ALMOST_EQUAL(0.000000001000001, 0.000000001000002);
    ASSERT_ALMOST_EQUAL(0.000000001000002, 0.000000001000001);
    ASSERT_NOT_EQUAL(0.000000000001002, 0.000000000001001);
    ASSERT_NOT_EQUAL(0.000000000001001, 0.000000000001002);
}

void test_small_negative_numbers(void) {
    ASSERT_ALMOST_EQUAL(-0.000000001000001, -0.000000001000002);
    ASSERT_ALMOST_EQUAL(-0.000000001000002, -0.000000001000001);
    ASSERT_NOT_EQUAL(-0.000000000001002, -0.000000000001001);
    ASSERT_NOT_EQUAL(-0.000000000001001, -0.000000000001002);
}

void test_small_differences_away_from_zero(void) {
    ASSERT_ALMOST_EQUAL(0.3, 0.30000003);
    ASSERT_ALMOST_EQUAL(-0.3, -0.30000003);
}

void test_comparisons_involving_zero(void) {
    ASSERT_ALMOST_EQUAL(0, 0);
    ASSERT_ALMOST_EQUAL(0.0, -0.0);
    ASSERT_ALMOST_EQUAL(-0.0, -0.0);
    ASSERT_NOT_EQUAL(0.00000001, 0.0);
    ASSERT_NOT_EQUAL(0.0, 0.00000001);
    ASSERT_NOT_EQUAL(-0.00000001, 0.0);
    ASSERT_NOT_EQUAL(0.0, -0.00000001);

    /* original test contained 1e-40 here, which is a denormalized number in
     * single-precision float world. An equivalent value for doubles is ~1e-320,
     * see : https://docs.oracle.com/javase/8/docs/api/constant-values.html#java.lang.Double.MIN_VALUE
     * The value must be between Double.MIN_VALUE and Double.MIN_NORMAL */
    ASSERT_ALMOST_EQUAL_WITH_EPS(0.0, 1e-320, 0.01);
    ASSERT_ALMOST_EQUAL_WITH_EPS(1e-320, 0.0, 0.01);
}

void test_extreme_values(void) {
    ASSERT_ALMOST_EQUAL(DBL_MAX, DBL_MAX);
    ASSERT_NOT_EQUAL(DBL_MAX, -DBL_MAX);
    ASSERT_NOT_EQUAL(-DBL_MAX, DBL_MAX);
    ASSERT_NOT_EQUAL(DBL_MAX, DBL_MAX / 2);
    ASSERT_NOT_EQUAL(DBL_MAX, -DBL_MAX / 2);
    ASSERT_NOT_EQUAL(-DBL_MAX, DBL_MAX / 2);
}

void test_infinities(void) {
    ASSERT_ALMOST_EQUAL(IGRAPH_INFINITY, IGRAPH_INFINITY);
    ASSERT_ALMOST_EQUAL(-IGRAPH_INFINITY, -IGRAPH_INFINITY);
    ASSERT_NOT_EQUAL(-IGRAPH_INFINITY, IGRAPH_INFINITY);
    ASSERT_NOT_EQUAL(IGRAPH_INFINITY, DBL_MAX);
    ASSERT_NOT_EQUAL(-IGRAPH_INFINITY, -DBL_MAX);
}

void test_nans(void) {
    ASSERT_NOT_EQUAL(IGRAPH_NAN, IGRAPH_NAN);
    ASSERT_NOT_EQUAL(IGRAPH_NAN, 0);
    ASSERT_NOT_EQUAL(-0.0, IGRAPH_NAN);
    ASSERT_NOT_EQUAL(IGRAPH_NAN, -0.0);
    ASSERT_NOT_EQUAL(IGRAPH_NAN, IGRAPH_INFINITY);
    ASSERT_NOT_EQUAL(IGRAPH_INFINITY, IGRAPH_NAN);
    ASSERT_NOT_EQUAL(IGRAPH_NAN, -IGRAPH_INFINITY);
    ASSERT_NOT_EQUAL(-IGRAPH_INFINITY, IGRAPH_NAN);
    ASSERT_NOT_EQUAL(IGRAPH_NAN, DBL_MAX);
    ASSERT_NOT_EQUAL(DBL_MAX, IGRAPH_NAN);
    ASSERT_NOT_EQUAL(IGRAPH_NAN, -DBL_MAX);
    ASSERT_NOT_EQUAL(-DBL_MAX, IGRAPH_NAN);
    ASSERT_NOT_EQUAL(IGRAPH_NAN, DBL_MIN);
    ASSERT_NOT_EQUAL(DBL_MIN, IGRAPH_NAN);
    ASSERT_NOT_EQUAL(IGRAPH_NAN, -DBL_MIN);
    ASSERT_NOT_EQUAL(-DBL_MIN, IGRAPH_NAN);
}

void test_opposite_sides_of_zero(void) {
    ASSERT_NOT_EQUAL(1.000000001, -1);
    ASSERT_NOT_EQUAL(-1, 1.000000001);
    ASSERT_NOT_EQUAL(-1.000000001, 1);
    ASSERT_NOT_EQUAL(1, -1.000000001);

    /* These tests from NearlyEqualsTest.java involved denormalized numbers in
     * Java world, and they were defined as floats. We use doubles, so I converted
     * the values manually by looking up Double.MIN_VALUE */
    ASSERT_ALMOST_EQUAL(49e-324, -49e-324);
}

void test_very_close_to_zero(void) {
#define DBL_DENORM_MIN 4.9e-324
    ASSERT_ALMOST_EQUAL(DBL_DENORM_MIN, DBL_DENORM_MIN);
    ASSERT_ALMOST_EQUAL(DBL_DENORM_MIN, -DBL_DENORM_MIN);
    ASSERT_ALMOST_EQUAL(-DBL_DENORM_MIN, DBL_DENORM_MIN);
    ASSERT_ALMOST_EQUAL(DBL_DENORM_MIN, 0);
    ASSERT_ALMOST_EQUAL(0, DBL_DENORM_MIN);
    ASSERT_ALMOST_EQUAL(-DBL_DENORM_MIN, 0);
    ASSERT_ALMOST_EQUAL(0, -DBL_DENORM_MIN);

    ASSERT_NOT_EQUAL(0.000000001, -DBL_DENORM_MIN);
    ASSERT_NOT_EQUAL(0.000000001, DBL_DENORM_MIN);
    ASSERT_NOT_EQUAL(DBL_DENORM_MIN, 0.000000001);
    ASSERT_NOT_EQUAL(-DBL_DENORM_MIN, 0.000000001);
}

int main(void) {
    test_large_numbers();
    test_large_negative_numbers();
    test_numbers_around_one();
    test_numbers_around_minus_one();
    test_small_numbers();
    test_small_negative_numbers();
    test_small_differences_away_from_zero();
    test_comparisons_involving_zero();
    test_extreme_values();
    test_infinities();
    test_nans();
    test_opposite_sides_of_zero();
    test_very_close_to_zero();

    VERIFY_FINALLY_STACK();

    return 0;
}
