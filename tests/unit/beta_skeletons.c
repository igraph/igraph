/*
   igraph library.
   Copyright (C) 2025  The igraph development team <igraph@igraph.org>

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

    igraph_t graph, graph2;
    igraph_bool_t is_same;

    igraph_real_t points[] = {

        0.474217, 0.0314797, 0.208089, 0.439308, 0.967367, 0.530466,
        0.177005, 0.426713, 0.568462, 0.57507, 0.441834, 0.284514, 0.479224,
        0.817988, 0.720209, 0.225744, 0.204941, 0.44297, 0.285318, 0.912984,
        0.831097, 0.0176603, 0.827154, 0.472702, 0.173059, 0.561858,
        0.156276, 0.88019, 0.65935, 0.538207, 0.570379, 0.518081, 0.900553,
        0.656416, 0.726631, 0.863709, 0.380264, 0.287159, 0.31098, 0.230773,
        0.243089, 0.164584, 0.967974, 0.524992, 0.726605, 0.0724703,
        0.739752, 0.447069, 0.0443581, 0.444839
    };

    igraph_real_t trig_lattice_points[] = {0.50000000000000000000000000000000000000000000000000, \
                                           2.5980762113533159402911695122588085504142078807156, 0, \
                                           1.7320508075688772935274463415058723669428052538104, \
                                           1.0000000000000000000000000000000000000000000000000, \
                                           1.7320508075688772935274463415058723669428052538104, \
                                           -0.50000000000000000000000000000000000000000000000000, \
                                           0.86602540378443864676372317075293618347140262690519, \
                                           0.50000000000000000000000000000000000000000000000000, \
                                           0.86602540378443864676372317075293618347140262690519, \
                                           1.5000000000000000000000000000000000000000000000000, \
                                           0.86602540378443864676372317075293618347140262690519, \
                                           -1.0000000000000000000000000000000000000000000000000, 0, 0, 0, \
                                           1.0000000000000000000000000000000000000000000000000, 0, \
                                           2.0000000000000000000000000000000000000000000000000, 0
                                          };
    igraph_real_t rotated_square_lattice_points[] = {
        -0.3594924531727418, 1.3677591805986329, -1.2231182700584293,
        1.8718925443115786, -2.0867440869441167, 2.3760259080245243,
        -2.950369903829804, 2.8801592717374698, 0.14464091054020378,
        2.2313849974843203, -0.7189849063454836, 2.7355183611972658,
        -1.582610723231171, 3.2396517249102117, -2.4462365401168586,
        3.743785088623157, 0.6487742742531495, 3.0950108143700077,
        -0.21485154263253792, 3.5991441780829536, -1.0784773595182253,
        4.103277541795899, -1.9421031764039127, 4.607410905508845,
        1.152907637966095, 3.958636631255695, 0.28928182108040756,
        4.462769994968641, -0.5743439958052798, 4.966903358681586,
        -1.4379698126909672, 5.4710367223945315
    };

    igraph_matrix_t point_mat, point_small_mat, point_singleton_mat, point_null_mat, point_3d_mat, trig_lattice, rot_square_lattice;

    igraph_matrix_init_array(&point_mat, points, 25, 2, false);
    igraph_matrix_init_array(&point_small_mat, points, 2, 2, false);
    igraph_matrix_init_array(&point_singleton_mat, points, 1, 2, false);
    igraph_matrix_init_array(&point_null_mat, points, 0, 2, false);

    igraph_matrix_init_array(&rot_square_lattice, rotated_square_lattice_points, 16, 2, false);
    igraph_matrix_init_array(&trig_lattice, trig_lattice_points, 10, 2, false);
    igraph_matrix_init_array(&point_3d_mat, points, 10, 3, false);

    printf("Lune beta skeleton beta = 2, 25 points\n");
    igraph_lune_beta_skeleton(&graph, &point_mat, 2);
    print_graph_canon(&graph);
    igraph_destroy(&graph);

    printf("Gabriel graph, 25 points\n");
    igraph_lune_beta_skeleton(&graph, &point_mat, 1);
    print_graph_canon(&graph);
    igraph_destroy(&graph);

    printf("Beta = 0.5, 25 points\n");
    igraph_lune_beta_skeleton(&graph, &point_mat, 0.5);
    print_graph_canon(&graph);
    igraph_destroy(&graph);

    printf("Beta = 1.1, Circle based 25 points\n");
    igraph_circle_beta_skeleton(&graph, &point_mat, 1.1);
    print_graph_canon(&graph);
    igraph_destroy(&graph);

    printf("Beta = 1.1, Circle based 2 points\n");
    igraph_circle_beta_skeleton(&graph, &point_small_mat, 1.1);
    print_graph_canon(&graph);
    igraph_destroy(&graph);

    printf("Beta = 1.1, Circle based 1 point\n");
    igraph_circle_beta_skeleton(&graph, &point_singleton_mat, 1.1);
    print_graph_canon(&graph);
    igraph_destroy(&graph);

    printf("Beta = 1.1, Circle based 0 points\n");
    igraph_circle_beta_skeleton(&graph, &point_null_mat, 1.1);
    print_graph_canon(&graph);
    igraph_destroy(&graph);


    printf("Relative neighborhood graph, 10 points 3d\n");
    igraph_lune_beta_skeleton(&graph, &point_3d_mat, 2);
    print_graph_canon(&graph);
    igraph_destroy(&graph);

    igraph_vector_t weights;


    printf("Beta weighted gabriel graph, 2d 25 points cutoff = Infinity\n");
    igraph_vector_init(&weights, 0);
    igraph_beta_weighted_gabriel_graph(&graph, &weights, &point_mat, IGRAPH_INFINITY);
    igraph_lune_beta_skeleton(&graph2, &point_mat, 1);
    igraph_is_same_graph(&graph, &graph2, &is_same);
    igraph_destroy(&graph2);
    igraph_destroy(&graph);
    print_vector(&weights);
    IGRAPH_ASSERT(is_same);
    igraph_vector_destroy(&weights);

    printf("Beta weighted gabriel graph, 2d 25 points cutoff = 5\n");
    igraph_vector_init(&weights, 0);
    igraph_beta_weighted_gabriel_graph(&graph, &weights, &point_mat, 5);
    igraph_destroy(&graph);
    print_vector(&weights);
    igraph_vector_destroy(&weights);

    printf("Beta weighted gabriel graph, 3d 10 points cutoff = Infinity\n");
    igraph_vector_init(&weights, 0);
    igraph_beta_weighted_gabriel_graph(&graph, &weights, &point_3d_mat, IGRAPH_INFINITY);
    igraph_lune_beta_skeleton(&graph2, &point_3d_mat, 1);
    igraph_is_same_graph(&graph, &graph2, &is_same);
    igraph_destroy(&graph2);
    igraph_destroy(&graph);
    print_vector(&weights);
    IGRAPH_ASSERT(is_same);
    igraph_vector_destroy(&weights);

    printf("Beta weighted gabriel graph, 3d 10 points cutoff = 5\n");
    igraph_vector_init(&weights, 0);
    igraph_beta_weighted_gabriel_graph(&graph, &weights, &point_3d_mat, 5);
    igraph_destroy(&graph);
    print_vector(&weights);
    igraph_vector_destroy(&weights);

    printf("Relative neighborhood graph, triangular lattice\n");
    igraph_relative_neighborhood_graph(&graph, &trig_lattice);
    print_graph_canon(&graph);
    igraph_destroy(&graph);
    printf("Lune beta skeleton beta = 2, triangular_lattice\n");
    igraph_lune_beta_skeleton(&graph, &trig_lattice, 2);
    print_graph_canon(&graph);
    igraph_destroy(&graph);

    printf("Gabriel graph of rotated square lattice\n");
    igraph_gabriel_graph(&graph, &rot_square_lattice);
    print_graph_canon(&graph);
    igraph_destroy(&graph);

    igraph_matrix_destroy(&rot_square_lattice);
    igraph_matrix_destroy(&trig_lattice);
    igraph_matrix_destroy(&point_3d_mat);
    igraph_matrix_destroy(&point_mat);
    igraph_matrix_destroy(&point_small_mat);
    igraph_matrix_destroy(&point_singleton_mat);
    igraph_matrix_destroy(&point_null_mat);
    VERIFY_FINALLY_STACK();
    return 0;
}
