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

/********************** Cartesian Product ************************/
void test_grid_vs_square_lattice(void) {
    // P4 X P3 should be a 4x3 grid (non-periodic)
    igraph_t p4, p3, product, lattice;
    igraph_bool_t is_iso;
    igraph_ring(&p4, 4, IGRAPH_UNDIRECTED, false, false);  // P4
    igraph_ring(&p3, 3, IGRAPH_UNDIRECTED, false, false);  // P3
    igraph_product(&product, &p4, &p3, IGRAPH_PRODUCT_CARTESIAN);

    igraph_vector_int_t dims;
    igraph_vector_int_init(&dims, 2);
    VECTOR(dims)[0] = 4;
    VECTOR(dims)[1] = 3;

    igraph_square_lattice(&lattice, &dims, 1, IGRAPH_UNDIRECTED, false, NULL);

    igraph_isomorphic(&product, &lattice, &is_iso);
    IGRAPH_ASSERT(is_iso);

    igraph_destroy(&p4);
    igraph_destroy(&p3);
    igraph_destroy(&product);
    igraph_destroy(&lattice);
    igraph_vector_int_destroy(&dims);
}

void test_cylinder_vs_cartesian(void) {
    // C4 X P3 should be a cylindrical lattice (periodic in one dimension)
    igraph_t c4, p3, product, cylinder;
    igraph_bool_t is_iso;
    igraph_ring(&c4, 4, IGRAPH_UNDIRECTED, false, true);  // C4
    igraph_ring(&p3, 3, IGRAPH_UNDIRECTED, false, false);  // P3
    igraph_product(&product, &c4, &p3, IGRAPH_PRODUCT_CARTESIAN);

    igraph_vector_int_t dims;
    igraph_vector_bool_t periodic;
    igraph_vector_int_init(&dims, 2);
    igraph_vector_bool_init(&periodic, 2);
    VECTOR(dims)[0] = 4;
    VECTOR(dims)[1] = 3;
    VECTOR(periodic)[0] = true;  // periodic in first dim (C4)
    VECTOR(periodic)[1] = false;  // non-periodic in second dim (P3)

    igraph_square_lattice(&cylinder, &dims, 1, IGRAPH_UNDIRECTED, false, &periodic);
    igraph_isomorphic(&product, &cylinder, &is_iso);
    IGRAPH_ASSERT(is_iso);

    igraph_destroy(&c4);
    igraph_destroy(&p3);
    igraph_destroy(&product);
    igraph_destroy(&cylinder);
    igraph_vector_int_destroy(&dims);
    igraph_vector_bool_destroy(&periodic);
}


void test_cube_vs_cartesian(void) {
    // K2 X K2 X K2 = Q3 (3-dimensional hypercube)
    igraph_t k2, temp, q3, cube;
    igraph_bool_t is_iso;

    igraph_full(&k2, 2, IGRAPH_UNDIRECTED, false); // K2
    igraph_product(&temp, &k2, &k2, IGRAPH_PRODUCT_CARTESIAN);
    igraph_product(&q3, &temp, &k2, IGRAPH_PRODUCT_CARTESIAN);

    igraph_hypercube(&cube, 3, IGRAPH_UNDIRECTED);
    igraph_isomorphic(&q3, &cube, &is_iso);

    IGRAPH_ASSERT(is_iso);

    igraph_destroy(&k2);
    igraph_destroy(&temp);
    igraph_destroy(&q3);
    igraph_destroy(&cube);
}

void test_torus_vs_cartesian(void) {
    // C4 X C4 = Torus
    igraph_t c4a, c4b, product, torus;
    igraph_bool_t is_iso;

    igraph_ring(&c4a, 4, IGRAPH_UNDIRECTED, false, true);
    igraph_ring(&c4b, 4, IGRAPH_UNDIRECTED, false, true);
    igraph_product(&product, &c4a, &c4b, IGRAPH_PRODUCT_CARTESIAN);

    igraph_vector_int_t dims;
    igraph_vector_bool_t periodic;
    igraph_vector_int_init(&dims, 2);
    igraph_vector_bool_init(&periodic, 2);
    VECTOR(dims)[0] = 4;
    VECTOR(dims)[1] = 4;
    VECTOR(periodic)[0] = true;
    VECTOR(periodic)[1] = true;

    igraph_square_lattice(&torus, &dims, 1, IGRAPH_UNDIRECTED, false, &periodic);
    igraph_isomorphic(&product, &torus, &is_iso);

    IGRAPH_ASSERT(is_iso);

    igraph_destroy(&c4a);
    igraph_destroy(&c4b);
    igraph_destroy(&product);
    igraph_destroy(&torus);
    igraph_vector_int_destroy(&dims);
    igraph_vector_bool_destroy(&periodic);
}

void test_self_loop_cartesian(void) {
    // 1-vertex loop X K4 (circular=false) = K4(circular=true)
    igraph_t v1loop, k4a, k4b, product;
    igraph_bool_t is_iso;

    igraph_full(&v1loop, 1, IGRAPH_UNDIRECTED, true);
    igraph_full(&k4a, 4, IGRAPH_UNDIRECTED, false);
    igraph_product(&product, &v1loop, &k4a, IGRAPH_PRODUCT_CARTESIAN);

    igraph_full(&k4b, 4, IGRAPH_UNDIRECTED, true);
    igraph_isomorphic(&product, &k4b, &is_iso);

    IGRAPH_ASSERT(is_iso);

    igraph_destroy(&v1loop);
    igraph_destroy(&k4a);
    igraph_destroy(&k4b);
    igraph_destroy(&product);
}

void test_multigraph_cartesian(void) {
    igraph_t g1, g2, product;

    // g1: a multigraph with 2 vertices and parallel edges + self-loop
    igraph_vector_int_t edges1;
    igraph_vector_int_init(&edges1, 6);
    VECTOR(edges1)[0] = 0; VECTOR(edges1)[1] = 1;
    VECTOR(edges1)[2] = 0; VECTOR(edges1)[3] = 1;  // Parallel edge
    VECTOR(edges1)[4] = 1; VECTOR(edges1)[5] = 1;  // Self-loop
    igraph_create(&g1, &edges1, 2, IGRAPH_DIRECTED);
    igraph_vector_int_destroy(&edges1);

    // g2: 3 vertices with some parallel edges + self-loop
    igraph_vector_int_t edges2;
    igraph_vector_int_init(&edges2, 8);
    VECTOR(edges2)[0] = 0; VECTOR(edges2)[1] = 1;
    VECTOR(edges2)[2] = 0; VECTOR(edges2)[3] = 1;  // Parallel edge
    VECTOR(edges2)[4] = 1; VECTOR(edges2)[5] = 2;
    VECTOR(edges2)[6] = 2; VECTOR(edges2)[7] = 2;  // Self-loop
    igraph_create(&g2, &edges2, 3, IGRAPH_DIRECTED);
    igraph_vector_int_destroy(&edges2);

    igraph_product(&product, &g1, &g2, IGRAPH_PRODUCT_CARTESIAN);

    // Check expected number of vertices: |V1| * |V2| = 2 * 3 = 6
    IGRAPH_ASSERT(igraph_vcount(&product) == 6);

    // verified by v1*e2 + v2*e1, see: https://en.wikipedia.org/wiki/Graph_product
    IGRAPH_ASSERT(igraph_ecount(&product) == 17);

    igraph_destroy(&g1);
    igraph_destroy(&g2);
    igraph_destroy(&product);
}

/********************** Tensor Product ************************/
// K2 X petersen = G(10,3)
void test_petersen_tensor(void) {
    igraph_t k2, petersen, g_10_3, product;
    igraph_bool_t is_iso;

    igraph_full(&k2, 2, IGRAPH_UNDIRECTED, false);
    igraph_famous(&petersen, "petersen");
    igraph_generalized_petersen(&g_10_3, 10, 3);

    igraph_product(&product, &k2, &petersen, IGRAPH_PRODUCT_TENSOR);

    igraph_isomorphic(&product, &g_10_3, &is_iso);

    IGRAPH_ASSERT(is_iso);

    igraph_destroy(&k2);
    igraph_destroy(&petersen);
    igraph_destroy(&g_10_3);
    igraph_destroy(&product);
}

// C2 X C3 = C6
void test_dir_cycle_tensor(void) {
    igraph_t c2, c3, c6, product;
    igraph_bool_t is_iso;

    igraph_ring(&c2, 2, IGRAPH_DIRECTED, false, true);
    igraph_ring(&c3, 3, IGRAPH_DIRECTED, false, true);
    igraph_ring(&c6, 6, IGRAPH_DIRECTED, false, true);

    igraph_product(&product, &c2, &c3, IGRAPH_PRODUCT_TENSOR);

    igraph_isomorphic(&product, &c6, &is_iso);

    IGRAPH_ASSERT(is_iso);

    igraph_destroy(&c2);
    igraph_destroy(&c3);
    igraph_destroy(&c6);
    igraph_destroy(&product);
}

/********************** Lexicographic Product ************************/
// K2 X K4 = K8
void test_k2_k4_lex(void) {
    igraph_t k2, k4, k8, product;
    igraph_bool_t is_iso;

    igraph_full(&k2, 2, IGRAPH_UNDIRECTED, false);
    igraph_full(&k4, 4, IGRAPH_UNDIRECTED, false);
    igraph_full(&k8, 8, IGRAPH_UNDIRECTED, false);

    igraph_product(&product, &k2, &k4, IGRAPH_PRODUCT_LEXICOGRAPHIC);
    igraph_isomorphic(&product, &k8, &is_iso);
    IGRAPH_ASSERT(is_iso);
    igraph_destroy(&k2);
    igraph_destroy(&k4);
    igraph_destroy(&k8);
    igraph_destroy(&product);
}

/********************** Strong Product ************************/
// P8 X P8 = King Graph
void test_p8_p8_strong(void) {
    igraph_t p8, king_graph, product;
    igraph_bool_t is_iso;

    igraph_vector_int_t edges;
    igraph_vector_int_init(&edges, 0);

    igraph_int_t n = 8;  // n x n chessboard

    for (igraph_int_t i = 0; i < n; i++) {
        for (igraph_int_t j = 0; j < n; j++) {
            igraph_int_t u = i * n + j;
            for (igraph_int_t dx = -1; dx <= 1; dx++) {
                for (igraph_int_t dy = -1; dy <= 1; dy++) {
                    if (dx == 0 && dy == 0) continue;
                    igraph_int_t ni = i + dx, nj = j + dy;
                    if (ni >= 0 && ni < n && nj >= 0 && nj < n) {
                        igraph_int_t v = ni * n + nj;
                        // To avoid duplicate edges, only add u < v
                        if (u < v) {
                            igraph_vector_int_push_back(&edges, u);
                            igraph_vector_int_push_back(&edges, v);
                        }
                    }
                }
            }
        }
    }
    igraph_create(&king_graph, &edges, n * n, IGRAPH_UNDIRECTED);
    igraph_ring(&p8, 8, IGRAPH_UNDIRECTED, false, false);

    igraph_product(&product, &p8, &p8, IGRAPH_PRODUCT_STRONG);
    igraph_isomorphic(&product, &king_graph, &is_iso);
    IGRAPH_ASSERT(is_iso);

    igraph_destroy(&king_graph);
    igraph_destroy(&product);
    igraph_destroy(&p8);
    igraph_vector_int_destroy(&edges);
}

/********************** Modular Product ************************/
// C4 X P2
void test_c4_p2_modular(void) {
    igraph_t c4, p2, res, prod;
    igraph_bool_t is_iso;

    igraph_ring(&c4, 4, IGRAPH_UNDIRECTED, false, true);
    igraph_ring(&p2, 2, IGRAPH_UNDIRECTED, false, false);

    igraph_product(&res, &c4, &p2, IGRAPH_PRODUCT_MODULAR);

    igraph_small(&prod, 8, IGRAPH_UNDIRECTED, 0,3,  0,7,  1,2,  1,6,  2,5,  3,4,  4,7,  5,6,  -1);

    igraph_isomorphic(&res, &prod, &is_iso);
    IGRAPH_ASSERT(is_iso);

    igraph_destroy(&c4);
    igraph_destroy(&p2);
    igraph_destroy(&res);
    igraph_destroy(&prod);
}

int main(void) {
    // CARTESIAN PRODUCT TEST
    test_grid_vs_square_lattice();
    test_cylinder_vs_cartesian();
    test_cube_vs_cartesian();
    test_torus_vs_cartesian();
    test_self_loop_cartesian();
    test_multigraph_cartesian();

    // TENSOR PRODUCT TEST
    test_petersen_tensor();
    test_dir_cycle_tensor();

    // LEXICOGRAPHIC PRODUCT TEST
    test_k2_k4_lex();

    // STRONG PRODUCT TEST
    test_p8_p8_strong();

    // MODULAR PRODUCT TEST
    test_c4_p2_modular();

    VERIFY_FINALLY_STACK();

    return 0;
}
