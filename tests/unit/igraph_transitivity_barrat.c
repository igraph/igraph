/*
   IGraph library.
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
#include "test_utilities.inc"

void warning_handler_print_stdout(const char *reason, const char *file,
                                  int line, int igraph_errno) {
    IGRAPH_UNUSED(igraph_errno);
    IGRAPH_UNUSED(file);
    IGRAPH_UNUSED(line);
    fprintf(stdout, "Warning: %s\n", reason);
}

int main() {
    igraph_t g_0, g_1, g_simple;
    igraph_vector_t result, weights_none, weights_simple;

    igraph_small(&g_0, 0, 0, -1);
    igraph_small(&g_1, 1, 0, -1);
    igraph_small(&g_simple, 6, 0, 0,1, 0,2, 1,2, 1,3, 2,3, 3,4, -1);

    igraph_vector_init(&result, 0);
    igraph_vector_init(&weights_none, 0);
    igraph_vector_init_int(&weights_simple, 6, -1, 0, 1, 2, 3, 4);

    igraph_set_warning_handler(warning_handler_print_stdout);

    printf("No vertices, transitivity zero, NULL weights:\n");
    IGRAPH_ASSERT(igraph_transitivity_barrat(&g_0, &result, igraph_vss_all(), /*weights*/ NULL, IGRAPH_TRANSITIVITY_ZERO) == IGRAPH_SUCCESS);
    print_vector(&result);
    printf("\n");

    printf("No vertices, transitivity zero, NULL weights, no vs:\n");
    IGRAPH_ASSERT(igraph_transitivity_barrat(&g_0, &result, igraph_vss_none(), /*weights*/ NULL, IGRAPH_TRANSITIVITY_ZERO) == IGRAPH_SUCCESS);
    print_vector(&result);
    printf("\n");

    printf("No vertices, transitivity zero, no weights:\n");
    IGRAPH_ASSERT(igraph_transitivity_barrat(&g_0, &result, igraph_vss_all(), &weights_none, IGRAPH_TRANSITIVITY_ZERO) == IGRAPH_SUCCESS);
    print_vector(&result);
    printf("\n");

    printf("No vertices, transitivity zero, no weights, no vs:\n");
    IGRAPH_ASSERT(igraph_transitivity_barrat(&g_0, &result, igraph_vss_none(), &weights_none, IGRAPH_TRANSITIVITY_ZERO) == IGRAPH_SUCCESS);
    print_vector(&result);
    printf("\n");

    printf("No vertices, transitivity NAN:\n");
    IGRAPH_ASSERT(igraph_transitivity_barrat(&g_0, &result, igraph_vss_all(), NULL, IGRAPH_TRANSITIVITY_NAN) == IGRAPH_SUCCESS);
    print_vector(&result);
    printf("\n");

    printf("No vertices, transitivity NAN, no vs:\n");
    IGRAPH_ASSERT(igraph_transitivity_barrat(&g_0, &result, igraph_vss_none(), NULL, IGRAPH_TRANSITIVITY_NAN) == IGRAPH_SUCCESS);
    print_vector(&result);
    printf("\n");

    printf("One vertex, transitivity zero:\n");
    IGRAPH_ASSERT(igraph_transitivity_barrat(&g_1, &result, igraph_vss_all(), NULL, IGRAPH_TRANSITIVITY_ZERO) == IGRAPH_SUCCESS);
    print_vector(&result);
    printf("\n");

    printf("One vertex, transitivity NaN:\n");
    IGRAPH_ASSERT(igraph_transitivity_barrat(&g_1, &result, igraph_vss_all(), NULL, IGRAPH_TRANSITIVITY_NAN) == IGRAPH_SUCCESS);
    print_vector(&result);
    printf("\n");

    printf("One vertex, transitivity zero, vs one:\n");
    IGRAPH_ASSERT(igraph_transitivity_barrat(&g_1, &result, igraph_vss_1(0), NULL, IGRAPH_TRANSITIVITY_ZERO) == IGRAPH_SUCCESS);
    print_vector(&result);
    printf("\n");

    printf("One vertex, transitivity NaN, vs one:\n");
    IGRAPH_ASSERT(igraph_transitivity_barrat(&g_1, &result, igraph_vss_1(0), NULL, IGRAPH_TRANSITIVITY_NAN) == IGRAPH_SUCCESS);
    print_vector(&result);
    printf("\n");

    printf("Simple graph, NULL weights, transitivity NAN:\n");
    IGRAPH_ASSERT(igraph_transitivity_barrat(&g_simple, &result, igraph_vss_all(), NULL, IGRAPH_TRANSITIVITY_NAN) == IGRAPH_SUCCESS);
    print_vector(&result);
    printf("\n");

    printf("Simple graph, with weights, transitivity NAN:\n");
    IGRAPH_ASSERT(igraph_transitivity_barrat(&g_simple, &result, igraph_vss_all(), &weights_simple, IGRAPH_TRANSITIVITY_NAN) == IGRAPH_SUCCESS);
    print_vector(&result);
    printf("\n");

    printf("Simple graph, with weights, transitivity zero:\n");
    IGRAPH_ASSERT(igraph_transitivity_barrat(&g_simple, &result, igraph_vss_all(), &weights_simple, IGRAPH_TRANSITIVITY_ZERO) == IGRAPH_SUCCESS);
    print_vector(&result);
    printf("\n");

    printf("Simple graph, with weights, transitivity zero, vss none:\n");
    IGRAPH_ASSERT(igraph_transitivity_barrat(&g_simple, &result, igraph_vss_none(), &weights_simple, IGRAPH_TRANSITIVITY_ZERO) == IGRAPH_SUCCESS);
    print_vector(&result);
    printf("\n");

    VERIFY_FINALLY_STACK();
    igraph_set_error_handler(igraph_error_handler_ignore);

    printf("Wrong weight length, vss all:\n");
    IGRAPH_ASSERT(igraph_transitivity_barrat(&g_simple, &result, igraph_vss_all(), &weights_none, IGRAPH_TRANSITIVITY_ZERO) == IGRAPH_EINVAL);

    printf("Wrong weight length, vss none:\n");
    IGRAPH_ASSERT(igraph_transitivity_barrat(&g_simple, &result, igraph_vss_none(), &weights_none, IGRAPH_TRANSITIVITY_ZERO) == IGRAPH_EINVAL);

    igraph_destroy(&g_0);
    igraph_destroy(&g_1);
    igraph_destroy(&g_simple);

    igraph_vector_destroy(&result);
    igraph_vector_destroy(&weights_none);
    igraph_vector_destroy(&weights_simple);

    VERIFY_FINALLY_STACK();
    return 0;
}
