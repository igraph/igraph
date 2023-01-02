/*
   IGraph library.
   Copyright (C) 2022  The igraph development team <igraph@igraph.org>

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
    igraph_vector_t v;
    igraph_vector_t v2;
    igraph_real_t result;
    igraph_bool_t result_bool;
    igraph_real_t nan[3] = {1, IGRAPH_NAN, 2};
    igraph_real_t basic[3] = {1, 5, 2};
    igraph_real_t basic_small[3] = {0, -5, -2};

    printf("Taking sum of squares of empty vector:\n");
    igraph_vector_view(&v, NULL, 0);
    result = igraph_vector_sumsq(&v);
    printf("%g\n", result);

    printf("Taking sum of squares of vector with NaN:\n");
    igraph_vector_view(&v, nan, 3);
    result = igraph_vector_sumsq(&v);
    print_real(stdout, result, "%g");

    printf("\nTaking sum of squares of vector:\n");
    igraph_vector_view(&v, basic, 3);
    result = igraph_vector_sumsq(&v);
    printf("%g\n", result);

    printf("Checking if vector is equal to itself:\n");
    igraph_vector_view(&v, basic, 3);
    igraph_vector_view(&v2, basic, 3);
    result_bool = igraph_vector_is_equal(&v, &v2);
    printf("%d\n", result_bool);

    printf("Checking if vector is equal to other vector:\n");
    igraph_vector_view(&v, basic, 3);
    igraph_vector_view(&v2, nan, 3);
    result_bool = igraph_vector_is_equal(&v, &v2);
    printf("%d\n", result_bool);

    printf("Checking if all elements of vector are less than its own elements:\n");
    igraph_vector_view(&v, basic, 3);
    igraph_vector_view(&v2, basic, 3);
    result_bool = igraph_vector_all_l(&v, &v2);
    printf("%d\n", result_bool);

    printf("Checking if all elements of vector are less than vector with nan:\n");
    igraph_vector_view(&v, basic, 3);
    igraph_vector_view(&v2, nan, 3);
    result_bool = igraph_vector_all_l(&v, &v2);
    printf("%d\n", result_bool);

    printf("Checking if all elements of vector are less than vector with smaller elements:\n");
    igraph_vector_view(&v, basic, 3);
    igraph_vector_view(&v2, basic_small, 3);
    result_bool = igraph_vector_all_l(&v, &v2);
    printf("%d\n", result_bool);

    printf("Checking if all elements of vector are less than vector with larger elements:\n");
    igraph_vector_view(&v, basic_small, 3);
    igraph_vector_view(&v2, basic, 3);
    result_bool = igraph_vector_all_l(&v, &v2);
    printf("%d\n", result_bool);

    printf("Checking if all elements of vector are greater than its own elements:\n");
    igraph_vector_view(&v, basic, 3);
    igraph_vector_view(&v2, basic, 3);
    result_bool = igraph_vector_all_g(&v, &v2);
    printf("%d\n", result_bool);

    printf("Checking if all elements of vector are greater than vector with nan:\n");
    igraph_vector_view(&v, basic, 3);
    igraph_vector_view(&v2, nan, 3);
    result_bool = igraph_vector_all_g(&v, &v2);
    printf("%d\n", result_bool);

    printf("Checking if all elements of vector are greater than vector with smaller elements:\n");
    igraph_vector_view(&v, basic, 3);
    igraph_vector_view(&v2, basic_small, 3);
    result_bool = igraph_vector_all_g(&v, &v2);
    printf("%d\n", result_bool);

    printf("Checking if all elements of vector are greater than vector with larger elements:\n");
    igraph_vector_view(&v, basic_small, 3);
    igraph_vector_view(&v2, basic, 3);
    result_bool = igraph_vector_all_g(&v, &v2);
    printf("%d\n", result_bool);

    printf("Checking vector_printf without nans:\n");
    igraph_vector_view(&v, basic, 3);
    igraph_vector_printf(&v, "--%4g--");

    VERIFY_FINALLY_STACK();
    return 0;
}
