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

igraph_error_t delaunay(
        igraph_real_t *points, igraph_int_t numpoints,
        igraph_int_t dimension, igraph_bool_t printing) {

    igraph_matrix_t points_mat;
    igraph_t g;

    IGRAPH_CHECK(igraph_matrix_init_array(&points_mat, points, numpoints, dimension, IGRAPH_ROW_MAJOR));
    IGRAPH_FINALLY(igraph_matrix_destroy, &points_mat);
    IGRAPH_CHECK(igraph_delaunay_graph(&g, &points_mat));
    IGRAPH_FINALLY(igraph_destroy, &g);

    if (printing) {
        printf("%" IGRAPH_PRId "\n", igraph_ecount(&g));
        print_graph_canon(&g);
    }

    igraph_matrix_destroy(&points_mat);
    igraph_destroy(&g);
    IGRAPH_FINALLY_CLEAN(2);
    return IGRAPH_SUCCESS;
}

#define PTCOUNT(points, dim) sizeof(points) / sizeof(points[0]) / dim

int main(void) {
    igraph_real_t points_raw[] = {
        0.311138, 0.59391, 0.67498, 0.1694, 0.822953, 0.904456, 0.129484,
        0.350282, 0.756002, 0.373154, 0.763904, 0.870096, 0.786728,
        0.00985041, 0.406027, 0.504641, 0.174089, 0.358568, 0.388813,
        0.744543, 0.0681765, 0.544161, 0.289917, 0.487095, 0.605402,
        0.270968, 0.978097, 0.955441, 0.495331, 0.688885, 0.693659, 0.579305,
        0.169171, 0.623525, 0.725786, 0.250287, 0.00648143, 0.743068,
        0.71552, 0.534675, 0.757967, 0.533673, 0.804339, 0.113337, 0.445796,
        0.410465, 0.38357, 0.918967, 0.906415, 0.893578, 0.501743, 0.914856,
        0.561324, 0.818718, 0.158764, 0.903373, 0.213128, 0.702836, 0.672937,
        0.291723, 0.855333, 0.0981787, 0.804335, 0.702915, 0.824296,
        0.294826, 0.702266, 0.346332, 0.223656, 0.640166, 0.404121, 0.151721,
        0.32827, 0.290019, 0.875388, 0.221269, 0.260946, 0.443243, 0.97165,
        0.11638, 0.246203, 0.447883, 0.987007, 0.0692745, 0.109563, 0.339652,
        0.0666818, 0.318156, 0.92505, 0.439913, 0.575609, 0.577994, 0.494698,
        0.335231, 0.5605, 0.701834, 0.863751, 0.129853, 0.641124, 0.277378,
        0.897763, 0.429955, 0.331589, 0.925289, 0.980783, 0.835259, 0.291818,
        0.693408, 0.631442, 0.53206, 0.159631, 0.350867, 0.699043, 0.993084,
        0.638868, 0.933556, 0.524775, 0.598018, 0.565078, 0.639017, 0.480785,
        0.0698051, 0.60897, 0.350966, 0.411467, 0.0688254, 0.430729,
        0.761323, 0.359096, 0.590563, 0.372144, 0.1754, 0.422717, 0.296634,
        0.874941, 0.32733, 0.616607, 0.547951, 0.902168, 0.369618, 0.226083,
        0.10907, 0.0655338, 0.456096, 0.633249, 0.676374, 0.9357, 0.945376,
        0.581192, 0.95893, 0.655664, 0.316474, 0.752296, 0.14192, 0.0125042,
        0.595269, 0.084296, 0.281774, 0.6996, 0.66766, 0.478567, 0.505426,
        0.716551, 0.232743, 0.191, 0.0114542, 0.609768, 0.903796, 0.772037,
        0.60114, 0.758929, 0.770639, 0.968881, 0.382973, 0.248955, 0.290848,
        0.752831, 0.0109218, 0.903305, 0.762493, 0.696122, 0.91684,
        0.0138796, 0.45168, 0.668332, 0.0793168, 0.499606, 0.321686,
        0.0934459, 0.499221, 0.92654, 0.917383, 0.479325, 0.147302, 0.9623,
        0.124589, 0.752746, 0.127419, 0.211051, 0.410824
    };

    igraph_real_t point_cloud_3d[] = {
        0.125269, 0.195444, 0.285089, 0.992534, 0.250088, 0.551803, 0.104234,
        0.99641, 0.72404, 0.00427073, 0.721804, 0.00514948, 0.0685917,
        0.504841, 0.884932, 0.320637, 0.22425, 0.107845, 0.460224, 0.923877,
        0.480411, 0.139484, 0.687409, 0.871486, 0.0273725, 0.518508,
        0.472275, 0.57603, 0.848422, 0.79629, 0.844017, 0.51724, 0.486733,
        0.886948, 0.41991, 0.626941, 0.575853, 0.399617, 0.415894, 0.478654,
        0.586383, 0.661319, 0.413907, 0.836301, 0.893907, 0.887531, 0.451777,
        0.392684, 0.488796, 0.1865, 0.518422, 0.087241, 0.424709, 0.210399,
        0.930499, 0.503204, 0.44653, 0.0644197, 0.375475, 0.510661, 0.931541,
        0.742055, 0.739896, 0.320497, 0.949405, 0.00282601, 0.160896,
        0.712019, 0.966189, 0.172367, 0.216689, 0.132008, 0.955499, 0.276767,
        0.271158, 0.085173, 0.306008, 0.621809, 0.826937, 0.702707, 0.367152,
        0.0914104, 0.315055, 0.726888, 0.759314, 0.759621, 0.0954736,
        0.35047, 0.322236, 0.676948, 0.70161, 0.967861, 0.80878, 0.345871,
        0.617758, 0.0890548, 0.947943, 0.93674, 0.00637263, 0.798722, 0.1123,
        0.883446, 0.270618, 0.859522, 0.527832, 0.504072, 0.561741, 0.614159,
        0.611909, 0.779264, 0.0651863, 0.760387, 0.934354, 0.00447453,
        0.652217, 0.00252835, 0.13965, 0.281917, 0.408229, 0.140262,
        0.951417, 0.945453, 0.8125, 0.404407, 0.457046, 0.839765, 0.695236,
        0.345371, 0.372619, 0.33457, 0.484926, 0.329866, 0.813254, 0.0181675,
        0.294805, 0.616989, 0.403054, 0.123351, 0.137132, 0.641034, 0.804549,
        0.776808, 0.193548, 0.0441299, 0.823612, 0.919071, 0.559754,
        0.636797, 0.0478928, 0.36961
    };

    igraph_real_t point_cloud_4d[] = {
        0.768464, 0.05359, 0.69102, 0.246049, 0.996004, 0.189829, 0.00302669,
        0.187796, 0.44531, 0.711051, 0.946157, 0.764919, 0.583837, 0.308961,
        0.219333, 0.885349, 0.480686, 0.620953, 0.0128075, 0.650813,
        0.717234, 0.157298, 0.775933, 0.816187, 0.189236, 0.194153, 0.347867,
        0.942655, 0.177974, 0.120843, 0.298, 0.779916, 0.00154216, 0.755666,
        0.537359, 0.00281348, 0.0687358, 0.298396, 0.19734, 0.940497
    };

    igraph_real_t points_lattice[] = {
        0,0,
        0,1,
        0,2,
        1,0,
        1,1,
        1,2,
        2,0,
        2,1,
        2,2
    };

    igraph_real_t points_collinear[] = {
        0,0,
        2,0,
        0,2,
        1,0
    };

    igraph_real_t points_not_quite_collinear[] = {
        0,0,
        2,0,
        0,2,
        1,0.0000000001
    };

    igraph_real_t points_all_collinear[] = {
        1,1,
        2,2,
        3,3,
        4,4,
    };

    igraph_real_t points_all_coplanar[] = {
            0., 0., 0.,
            -0.9852739637886722, 0.00098610486426615, 0.17098024411421095,
            0.14677925443983744, 0.5177782140612408, 0.842829503226861,
            -0.8384947093488347, 0.518764318925507, 1.0138097473410719
    };

    igraph_real_t points_duplicate[] ={
        1,1,
        0,1,
        1,0,
        1,1
    };

    printf("100 point cloud in 2d\n");
    delaunay(points_raw, PTCOUNT(points_raw, 2), 2, true);

    printf("50 point cloud in 3d\n");
    delaunay(point_cloud_3d, PTCOUNT(point_cloud_3d, 3), 3, true);

    printf("10 point cloud in 4d\n");
    delaunay(point_cloud_4d, PTCOUNT(point_cloud_4d, 4), 4, true);

    printf("3x3 square lattice\n");
    delaunay(points_lattice, PTCOUNT(points_lattice, 2), 2, true);

    printf("Triangle with one subdivided edge\n");
    delaunay(points_collinear, PTCOUNT(points_collinear, 2), 2, true);

    printf("Triangle with one subdivided edge, point then slightly perturbed\n");
    delaunay(points_not_quite_collinear, PTCOUNT(points_not_quite_collinear, 2), 2, true);

    /* TODO: Planned to be supported without error in the future. */
    printf("4 points on a line in 2d\n");
    CHECK_ERROR(delaunay(points_all_collinear, PTCOUNT(points_all_collinear, 2), 2, true), IGRAPH_EINVAL);

    /* TODO: Planned to be supported without error in the future. */
    printf("4 points on a plane in 3d\n");
    CHECK_ERROR(delaunay(points_all_coplanar, PTCOUNT(points_all_coplanar, 3), 3, true), IGRAPH_EINVAL);

    printf("Identical points\n");
    CHECK_ERROR(delaunay(points_duplicate, PTCOUNT(points_duplicate, 2), 2,false), IGRAPH_EINVAL);

    /* The below coordinate arrays are used only partially, hence no PTCOUNT */

    printf("One point in 2d\n");
    delaunay(points_raw, 1, 2, true);

    printf("No points in 2d\n");
    delaunay(points_raw, 0, 2, true);

    printf("0x0 matrix\n");
    delaunay(points_raw, 0, 0, true);

    printf("2 points in 0d\n");
    CHECK_ERROR(delaunay(points_raw, 2, 0, true), IGRAPH_EINVAL);

    printf("4 points in 3d\n");
    delaunay(point_cloud_3d, 4, 3, true);

    /* TODO: Planned to be supported without error in the future. */
    printf("3 points in 3d\n");
    CHECK_ERROR(delaunay(point_cloud_3d, 3, 3, true), IGRAPH_EINVAL);

    VERIFY_FINALLY_STACK();

    return 0;
}
