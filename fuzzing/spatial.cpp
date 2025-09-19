/*
   igraph library.
   Copyright (C) 2025  The igraph development team

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
   02110-1301 USA
*/

#include <igraph.h>

#include <limits>
#include <random>

/* Custom mutator and crossover functions are based on
 * https://rigtorp.se/fuzzing-floating-point-code/
 * A more advanced solution will take into account that we
 * are working with n-dimensional coordinates. */

extern "C" size_t LLVMFuzzerCustomMutator(uint8_t *Data, size_t Size,
                                          size_t MaxSize, unsigned int Seed) {
    double *begin = (double *) Data;
    double *end = (double *) Data + Size / sizeof(double);

    std::minstd_rand gen(Seed);

    auto rfp = [&]() {
        if (gen() < std::bernoulli_distribution(0.8)(gen)) {
            std::normal_distribution<> dist(0.0, 10);
            return dist(gen);
        } else {
            switch (std::uniform_int_distribution<>(0, 9)(gen)) {
            case 0:
                return std::numeric_limits<double>::quiet_NaN();
            case 1:
                return std::numeric_limits<double>::min();
            case 2:
                return std::numeric_limits<double>::max();
            case 3:
                return -std::numeric_limits<double>::min();
            case 4:
                return -std::numeric_limits<double>::max();
            case 5:
                return std::numeric_limits<double>::epsilon();
            case 6:
                return -std::numeric_limits<double>::epsilon();
            case 7:
                return std::numeric_limits<double>::infinity();
            case 8:
                return -std::numeric_limits<double>::infinity();
            case 9:
                return 0.0;
            }
        }
        return 0.0; // never reached, prevents compiler warning
    };

    switch (std::uniform_int_distribution<>(0, 3)(gen)) {
    case 0: { // Change element
        if (begin != end) {
            std::uniform_int_distribution<> d(0, end - begin - 1);
            begin[d(gen)] = rfp();
        }
        break;
    }
    case 1: // Add element
        if (Size + sizeof(double) <= MaxSize) {
            *end = rfp();
            ++end;
        }
        break;
    case 2: // Delete element
        if (begin != end) {
            --end;
        }
        break;
    case 3: // Shuffle elements
        std::shuffle(begin, end, gen);
        break;
    }

    return (end - begin) * sizeof(double);
}


extern "C" size_t LLVMFuzzerCustomCrossOver(const uint8_t *Data1, size_t Size1,
                                            const uint8_t *Data2, size_t Size2,
                                            uint8_t *Out, size_t MaxOutSize,
                                            unsigned int Seed) {
    std::minstd_rand gen(Seed);
    std::bernoulli_distribution bd(0.5);
    size_t n = std::min({Size1, Size2, MaxOutSize}) / sizeof(double);
    for (size_t i = 0; i < n; ++i) {
        ((double *)Out)[i] = bd(gen) ? ((double *)Data1)[i] : ((double *)Data2)[i];
    }
    return n * sizeof(double);
}


extern "C" int LLVMFuzzerTestOneInput(const uint8_t *Data, size_t Size) {
    igraph_matrix_t points1d, points2d, points3d, points4d;
    igraph_t graph;

    igraph_set_error_handler(igraph_error_handler_ignore);
    igraph_set_warning_handler(igraph_warning_handler_ignore);

    if (Size % sizeof(igraph_real_t) != 0 || Size > 120 * sizeof(igraph_real_t)) {
        return 0;
    }

    size_t coord_count = Size / sizeof(double);
    igraph_real_t *coords = (igraph_real_t *) Data;

    points1d = igraph_matrix_view(coords, coord_count / 1, 1);
    points2d = igraph_matrix_view(coords, coord_count / 2, 2);
    points3d = igraph_matrix_view(coords, coord_count / 3, 3);
    points4d = igraph_matrix_view(coords, coord_count / 4, 4);

    {
        igraph_vector_int_t iv;
        igraph_vector_int_init(&iv, 0);
        igraph_convex_hull_2d(&points2d, &iv, NULL);
        igraph_vector_int_destroy(&iv);
    }

    if (igraph_nearest_neighbor_graph(&graph, &points1d, IGRAPH_METRIC_EUCLIDEAN, 2, IGRAPH_INFINITY, IGRAPH_UNDIRECTED) == IGRAPH_SUCCESS) {
        igraph_destroy(&graph);
    }
    if (igraph_nearest_neighbor_graph(&graph, &points2d, IGRAPH_METRIC_EUCLIDEAN, 1, IGRAPH_INFINITY, IGRAPH_DIRECTED) == IGRAPH_SUCCESS) {
        igraph_destroy(&graph);
    }
    if (igraph_nearest_neighbor_graph(&graph, &points3d, IGRAPH_METRIC_EUCLIDEAN, 3, IGRAPH_INFINITY, IGRAPH_UNDIRECTED) == IGRAPH_SUCCESS) {
        igraph_destroy(&graph);
    }
    if (igraph_nearest_neighbor_graph(&graph, &points4d, IGRAPH_METRIC_EUCLIDEAN, 2, IGRAPH_INFINITY, IGRAPH_DIRECTED) == IGRAPH_SUCCESS) {
        igraph_destroy(&graph);
    }

    /* Re-enable after https://github.com/qhull/qhull/issues/166 is fixed. */
    /*
    if (igraph_delaunay_graph(&graph, &points1d) == IGRAPH_SUCCESS) {
        igraph_destroy(&graph);
    }
    if (igraph_delaunay_graph(&graph, &points2d) == IGRAPH_SUCCESS) {
        igraph_destroy(&graph);
    }
    if (igraph_delaunay_graph(&graph, &points3d) == IGRAPH_SUCCESS) {
        igraph_destroy(&graph);
    }
    if (igraph_delaunay_graph(&graph, &points4d) == IGRAPH_SUCCESS) {
        igraph_destroy(&graph);
    }
    */

    IGRAPH_ASSERT(IGRAPH_FINALLY_STACK_EMPTY);

    return 0;  // Non-zero return values are reserved for future use.
}
