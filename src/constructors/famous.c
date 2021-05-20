/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2005-2021 The igraph development team

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

#include <string.h>

#include "igraph_constructors.h"

#include "internal/hacks.h"

const igraph_real_t igraph_i_famous_bull[] = {
    5, 5, 0,
    0, 1, 0, 2, 1, 2, 1, 3, 2, 4
};

const igraph_real_t igraph_i_famous_chvatal[] = {
    12, 24, 0,
    5, 6, 6, 7, 7, 8, 8, 9, 5, 9, 4, 5, 4, 8, 2, 8, 2, 6, 0, 6, 0, 9, 3, 9, 3, 7,
    1, 7, 1, 5, 1, 10, 4, 10, 4, 11, 2, 11, 0, 10, 0, 11, 3, 11, 3, 10, 1, 2
};

const igraph_real_t igraph_i_famous_coxeter[] = {
    28, 42, 0,
    0, 1, 0, 2, 0, 7, 1, 4, 1, 13, 2, 3, 2, 8, 3, 6, 3, 9, 4, 5, 4, 12, 5, 6, 5,
    11, 6, 10, 7, 19, 7, 24, 8, 20, 8, 23, 9, 14, 9, 22, 10, 15, 10, 21, 11, 16,
    11, 27, 12, 17, 12, 26, 13, 18, 13, 25, 14, 17, 14, 18, 15, 18, 15, 19, 16, 19,
    16, 20, 17, 20, 21, 23, 21, 26, 22, 24, 22, 27, 23, 25, 24, 26, 25, 27
};

const igraph_real_t igraph_i_famous_cubical[] = {
    8, 12, 0,
    0, 1, 1, 2, 2, 3, 0, 3, 4, 5, 5, 6, 6, 7, 4, 7, 0, 4, 1, 5, 2, 6, 3, 7
};

const igraph_real_t igraph_i_famous_diamond[] = {
    4, 5, 0,
    0, 1, 0, 2, 1, 2, 1, 3, 2, 3
};

const igraph_real_t igraph_i_famous_dodecahedron[] = {
    20, 30, 0,
    0, 1, 0, 4, 0, 5, 1, 2, 1, 6, 2, 3, 2, 7, 3, 4, 3, 8, 4, 9, 5, 10, 5, 11, 6,
    10, 6, 14, 7, 13, 7, 14, 8, 12, 8, 13, 9, 11, 9, 12, 10, 15, 11, 16, 12, 17,
    13, 18, 14, 19, 15, 16, 15, 19, 16, 17, 17, 18, 18, 19
};

const igraph_real_t igraph_i_famous_folkman[] = {
    20, 40, 0,
    0, 5, 0, 8, 0, 10, 0, 13, 1, 7, 1, 9, 1, 12, 1, 14, 2, 6, 2, 8, 2, 11, 2, 13,
    3, 5, 3, 7, 3, 10, 3, 12, 4, 6, 4, 9, 4, 11, 4, 14, 5, 15, 5, 19, 6, 15, 6, 16,
    7, 16, 7, 17, 8, 17, 8, 18, 9, 18, 9, 19, 10, 15, 10, 19, 11, 15, 11, 16, 12,
    16, 12, 17, 13, 17, 13, 18, 14, 18, 14, 19
};

const igraph_real_t igraph_i_famous_franklin[] = {
    12, 18, 0,
    0, 1, 0, 2, 0, 6, 1, 3, 1, 7, 2, 4, 2, 10, 3, 5, 3, 11, 4, 5, 4, 6, 5, 7, 6, 8,
    7, 9, 8, 9, 8, 11, 9, 10, 10, 11
};

const igraph_real_t igraph_i_famous_frucht[] = {
    12, 18, 0,
    0, 1, 0, 2, 0, 11, 1, 3, 1, 6, 2, 5, 2, 10, 3, 4, 3, 6, 4, 8, 4, 11, 5, 9, 5,
    10, 6, 7, 7, 8, 7, 9, 8, 9, 10, 11
};

const igraph_real_t igraph_i_famous_grotzsch[] = {
    11, 20, 0,
    0, 1, 0, 2, 0, 7, 0, 10, 1, 3, 1, 6, 1, 9, 2, 4, 2, 6, 2, 8, 3, 4, 3, 8, 3, 10,
    4, 7, 4, 9, 5, 6, 5, 7, 5, 8, 5, 9, 5, 10
};

const igraph_real_t igraph_i_famous_heawood[] = {
    14, 21, 0,
    0, 1, 0, 5, 0, 13, 1, 2, 1, 10, 2, 3, 2, 7, 3, 4, 3, 12, 4, 5, 4, 9, 5, 6, 6,
    7, 6, 11, 7, 8, 8, 9, 8, 13, 9, 10, 10, 11, 11, 12, 12, 13
};

const igraph_real_t igraph_i_famous_herschel[] = {
    11, 18, 0,
    0, 2, 0, 3, 0, 4, 0, 5, 1, 2, 1, 3, 1, 6, 1, 7, 2, 10, 3, 9, 4, 8, 4, 9, 5, 8,
    5, 10, 6, 8, 6, 9, 7, 8, 7, 10
};

const igraph_real_t igraph_i_famous_house[] = {
    5, 6, 0,
    0, 1, 0, 2, 1, 3, 2, 3, 2, 4, 3, 4
};

const igraph_real_t igraph_i_famous_housex[] = {
    5, 8, 0,
    0, 1, 0, 2, 0, 3, 1, 2, 1, 3, 2, 3, 2, 4, 3, 4
};

const igraph_real_t igraph_i_famous_icosahedron[] = {
    12, 30, 0,
    0, 1, 0, 2, 0, 3, 0, 4, 0, 8, 1, 2, 1, 6, 1, 7, 1, 8, 2, 4, 2, 5, 2, 6, 3, 4,
    3, 8, 3, 9, 3, 11, 4, 5, 4, 11, 5, 6, 5, 10, 5, 11, 6, 7, 6, 10, 7, 8, 7, 9, 7,
    10, 8, 9, 9, 10, 9, 11, 10, 11
};

const igraph_real_t igraph_i_famous_krackhardt_kite[] = {
    10, 18, 0,
    0, 1, 0, 2, 0, 3, 0, 5, 1, 3, 1, 4, 1, 6, 2, 3, 2, 5, 3, 4, 3, 5, 3, 6, 4, 6, 5, 6, 5, 7, 6, 7, 7, 8, 8, 9
};

const igraph_real_t igraph_i_famous_levi[] = {
    30, 45, 0,
    0, 1, 0, 7, 0, 29, 1, 2, 1, 24, 2, 3, 2, 11, 3, 4, 3, 16, 4, 5, 4, 21, 5, 6, 5,
    26, 6, 7, 6, 13, 7, 8, 8, 9, 8, 17, 9, 10, 9, 22, 10, 11, 10, 27, 11, 12, 12,
    13, 12, 19, 13, 14, 14, 15, 14, 23, 15, 16, 15, 28, 16, 17, 17, 18, 18, 19, 18,
    25, 19, 20, 20, 21, 20, 29, 21, 22, 22, 23, 23, 24, 24, 25, 25, 26, 26, 27, 27,
    28, 28, 29
};

const igraph_real_t igraph_i_famous_mcgee[] = {
    24, 36, 0,
    0, 1, 0, 7, 0, 23, 1, 2, 1, 18, 2, 3, 2, 14, 3, 4, 3, 10, 4, 5, 4, 21, 5, 6, 5,
    17, 6, 7, 6, 13, 7, 8, 8, 9, 8, 20, 9, 10, 9, 16, 10, 11, 11, 12, 11, 23, 12,
    13, 12, 19, 13, 14, 14, 15, 15, 16, 15, 22, 16, 17, 17, 18, 18, 19, 19, 20, 20,
    21, 21, 22, 22, 23
};

const igraph_real_t igraph_i_famous_meredith[] = {
    70, 140, 0,
    0, 4, 0, 5, 0, 6, 1, 4, 1, 5, 1, 6, 2, 4, 2, 5, 2, 6, 3, 4, 3, 5, 3, 6, 7, 11,
    7, 12, 7, 13, 8, 11, 8, 12, 8, 13, 9, 11, 9, 12, 9, 13, 10, 11, 10, 12, 10, 13,
    14, 18, 14, 19, 14, 20, 15, 18, 15, 19, 15, 20, 16, 18, 16, 19, 16, 20, 17, 18,
    17, 19, 17, 20, 21, 25, 21, 26, 21, 27, 22, 25, 22, 26, 22, 27, 23, 25, 23, 26,
    23, 27, 24, 25, 24, 26, 24, 27, 28, 32, 28, 33, 28, 34, 29, 32, 29, 33, 29, 34,
    30, 32, 30, 33, 30, 34, 31, 32, 31, 33, 31, 34, 35, 39, 35, 40, 35, 41, 36, 39,
    36, 40, 36, 41, 37, 39, 37, 40, 37, 41, 38, 39, 38, 40, 38, 41, 42, 46, 42, 47,
    42, 48, 43, 46, 43, 47, 43, 48, 44, 46, 44, 47, 44, 48, 45, 46, 45, 47, 45, 48,
    49, 53, 49, 54, 49, 55, 50, 53, 50, 54, 50, 55, 51, 53, 51, 54, 51, 55, 52, 53,
    52, 54, 52, 55, 56, 60, 56, 61, 56, 62, 57, 60, 57, 61, 57, 62, 58, 60, 58, 61,
    58, 62, 59, 60, 59, 61, 59, 62, 63, 67, 63, 68, 63, 69, 64, 67, 64, 68, 64, 69,
    65, 67, 65, 68, 65, 69, 66, 67, 66, 68, 66, 69, 2, 50, 1, 51, 9, 57, 8, 58, 16,
    64, 15, 65, 23, 36, 22, 37, 30, 43, 29, 44, 3, 21, 7, 24, 14, 31, 0, 17, 10,
    28, 38, 42, 35, 66, 59, 63, 52, 56, 45, 49
};

const igraph_real_t igraph_i_famous_noperfectmatching[] = {
    16, 27, 0,
    0, 1, 0, 2, 0, 3, 1, 2, 1, 3, 2, 3, 2, 4, 3, 4, 4, 5, 5, 6, 5, 7, 6, 12, 6, 13,
    7, 8, 7, 9, 8, 9, 8, 10, 8, 11, 9, 10, 9, 11, 10, 11, 12, 13, 12, 14, 12, 15,
    13, 14, 13, 15, 14, 15
};

const igraph_real_t igraph_i_famous_nonline[] = {
    50, 72, 0,
    0, 1, 0, 2, 0, 3, 4, 6, 4, 7, 5, 6, 5, 7, 6, 7, 7, 8, 9, 11, 9, 12, 9, 13, 10,
    11, 10, 12, 10, 13, 11, 12, 11, 13, 12, 13, 14, 15, 15, 16, 15, 17, 16, 17, 16,
    18, 17, 18, 18, 19, 20, 21, 20, 22, 20, 23, 21, 22, 21, 23, 21, 24, 22, 23, 22,
    24, 24, 25, 26, 27, 26, 28, 26, 29, 27, 28, 27, 29, 27, 30, 27, 31, 28, 29, 28,
    30, 28, 31, 30, 31, 32, 34, 32, 35, 32, 36, 33, 34, 33, 35, 33, 37, 34, 35, 36,
    37, 38, 39, 38, 40, 38, 43, 39, 40, 39, 41, 39, 42, 39, 43, 40, 41, 41, 42, 42,
    43, 44, 45, 44, 46, 45, 46, 45, 47, 46, 47, 46, 48, 47, 48, 47, 49, 48, 49
};

const igraph_real_t igraph_i_famous_octahedron[] = {
    6, 12, 0,
    0, 1, 0, 2, 1, 2, 3, 4, 3, 5, 4, 5, 0, 3, 0, 5, 1, 3, 1, 4, 2, 4, 2, 5
};

const igraph_real_t igraph_i_famous_petersen[] = {
    10, 15, 0,
    0, 1, 0, 4, 0, 5, 1, 2, 1, 6, 2, 3, 2, 7, 3, 4, 3, 8, 4, 9, 5, 7, 5, 8, 6, 8, 6, 9, 7, 9
};

const igraph_real_t igraph_i_famous_robertson[] = {
    19, 38, 0,
    0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 10, 11, 11, 12,
    12, 13, 13, 14, 14, 15, 15, 16, 16, 17, 17, 18, 0, 18, 0, 4, 4, 9, 9, 13, 13,
    17, 2, 17, 2, 6, 6, 10, 10, 15, 0, 15, 1, 8, 8, 16, 5, 16, 5, 12, 1, 12, 7, 18,
    7, 14, 3, 14, 3, 11, 11, 18
};

const igraph_real_t igraph_i_famous_smallestcyclicgroup[] = {
    9, 15, 0,
    0, 1, 0, 2, 0, 3, 0, 4, 0, 5, 1, 2, 1, 3, 1, 7, 1, 8, 2, 5, 2, 6, 2, 7, 3, 8,
    4, 5, 6, 7
};

const igraph_real_t igraph_i_famous_tetrahedron[] = {
    4, 6, 0,
    0, 3, 1, 3, 2, 3, 0, 1, 1, 2, 0, 2
};

const igraph_real_t igraph_i_famous_thomassen[] = {
    34, 52, 0,
    0, 2, 0, 3, 1, 3, 1, 4, 2, 4, 5, 7, 5, 8, 6, 8, 6, 9, 7, 9, 10, 12, 10, 13, 11,
    13, 11, 14, 12, 14, 15, 17, 15, 18, 16, 18, 16, 19, 17, 19, 9, 19, 4, 14, 24,
    25, 25, 26, 20, 26, 20, 21, 21, 22, 22, 23, 23, 27, 27, 28, 28, 29, 29, 30, 30,
    31, 31, 32, 32, 33, 24, 33, 5, 24, 6, 25, 7, 26, 8, 20, 0, 20, 1, 21, 2, 22, 3,
    23, 10, 27, 11, 28, 12, 29, 13, 30, 15, 30, 16, 31, 17, 32, 18, 33
};

const igraph_real_t igraph_i_famous_tutte[] = {
    46, 69, 0,
    0, 10, 0, 11, 0, 12, 1, 2, 1, 7, 1, 19, 2, 3, 2, 41, 3, 4, 3, 27, 4, 5, 4, 33,
    5, 6, 5, 45, 6, 9, 6, 29, 7, 8, 7, 21, 8, 9, 8, 22, 9, 24, 10, 13, 10, 14, 11,
    26, 11, 28, 12, 30, 12, 31, 13, 15, 13, 21, 14, 15, 14, 18, 15, 16, 16, 17, 16,
    20, 17, 18, 17, 23, 18, 24, 19, 25, 19, 40, 20, 21, 20, 22, 22, 23, 23, 24, 25,
    26, 25, 38, 26, 34, 27, 28, 27, 39, 28, 34, 29, 30, 29, 44, 30, 35, 31, 32, 31,
    35, 32, 33, 32, 42, 33, 43, 34, 36, 35, 37, 36, 38, 36, 39, 37, 42, 37, 44, 38,
    40, 39, 41, 40, 41, 42, 43, 43, 45, 44, 45
};

const igraph_real_t igraph_i_famous_uniquely3colorable[] = {
    12, 22, 0,
    0, 1, 0, 3, 0, 6, 0, 8, 1, 4, 1, 7, 1, 9, 2, 3, 2, 6, 2, 7, 2, 9, 2, 11, 3, 4,
    3, 10, 4, 5, 4, 11, 5, 6, 5, 7, 5, 8, 5, 10, 8, 11, 9, 10
};

const igraph_real_t igraph_i_famous_walther[] = {
    25, 31, 0,
    0, 1, 1, 2, 1, 8, 2, 3, 2, 13, 3, 4, 3, 16, 4, 5, 5, 6, 5, 19, 6, 7, 6, 20, 7,
    21, 8, 9, 8, 13, 9, 10, 9, 22, 10, 11, 10, 20, 11, 12, 13, 14, 14, 15, 14, 23,
    15, 16, 15, 17, 17, 18, 18, 19, 18, 24, 20, 24, 22, 23, 23, 24
};

const igraph_real_t igraph_i_famous_zachary[] = {
    34, 78, 0,
    0, 1, 0, 2, 0, 3, 0, 4, 0, 5, 0, 6, 0, 7, 0, 8,
    0, 10, 0, 11, 0, 12, 0, 13, 0, 17, 0, 19, 0, 21, 0, 31,
    1, 2, 1, 3, 1, 7, 1, 13, 1, 17, 1, 19, 1, 21, 1, 30,
    2, 3, 2, 7, 2, 27, 2, 28, 2, 32, 2, 9, 2, 8, 2, 13,
    3, 7, 3, 12, 3, 13, 4, 6, 4, 10, 5, 6, 5, 10, 5, 16,
    6, 16, 8, 30, 8, 32, 8, 33, 9, 33, 13, 33, 14, 32, 14, 33,
    15, 32, 15, 33, 18, 32, 18, 33, 19, 33, 20, 32, 20, 33,
    22, 32, 22, 33, 23, 25, 23, 27, 23, 32, 23, 33, 23, 29,
    24, 25, 24, 27, 24, 31, 25, 31, 26, 29, 26, 33, 27, 33,
    28, 31, 28, 33, 29, 32, 29, 33, 30, 32, 30, 33, 31, 32, 31, 33,
    32, 33
};

static int igraph_i_famous(igraph_t *graph, const igraph_real_t *data) {
    long int no_of_nodes = (long int) data[0];
    long int no_of_edges = (long int) data[1];
    igraph_bool_t directed = (igraph_bool_t) data[2];
    igraph_vector_t edges;

    igraph_vector_view(&edges, data + 3, 2 * no_of_edges);
    IGRAPH_CHECK(igraph_create(graph, &edges, (igraph_integer_t) no_of_nodes,
                               directed));
    return 0;
}

/**
 * \function igraph_famous
 * \brief Create a famous graph by simply providing its name.
 *
 * </para><para>
 * The name of the graph can be simply supplied as a string.
 * Note that this function creates graphs which don't take any parameters,
 * there are separate functions for graphs with parameters, e.g. \ref
 * igraph_full() for creating a full graph.
 *
 * </para><para>
 * The following graphs are supported:
 * \clist
 *   \cli Bull
 *           The bull graph, 5 vertices, 5 edges, resembles the
 *           head of a bull if drawn properly.
 *   \cli Chvatal
 *           This is the smallest triangle-free graph that is
 *           both 4-chromatic and 4-regular. According to the Grunbaum
 *           conjecture there exists an m-regular, m-chromatic graph
 *           with n vertices for every m>1 and n>2. The Chvatal graph
 *           is an example for m=4 and n=12. It has 24 edges.
 *   \cli Coxeter
 *           A non-Hamiltonian cubic symmetric graph with 28
 *           vertices and 42 edges.
 *   \cli Cubical
 *           The Platonic graph of the cube. A convex regular
 *           polyhedron with 8 vertices and 12 edges.
 *   \cli Diamond
 *           A graph with 4 vertices and 5 edges, resembles a
 *           schematic diamond if drawn properly.
 *   \cli Dodecahedral, Dodecahedron
 *           Another Platonic solid
 *           with 20 vertices and 30 edges.
 *   \cli Folkman
 *           The semisymmetric graph with minimum number of
 *           vertices, 20 and 40 edges. A semisymmetric graph is
 *           regular, edge transitive and not vertex transitive.
 *   \cli Franklin
 *           This is a graph whose embedding to the Klein
 *           bottle can be colored with six colors, it is a
 *           counterexample to the necessity of the Heawood
 *           conjecture on a Klein bottle. It has 12 vertices and 18
 *           edges.
 *   \cli Frucht
 *           The Frucht Graph is the smallest cubical graph
 *           whose automorphism group consists only of the identity
 *           element. It has 12 vertices and 18 edges.
 *   \cli Grotzsch
 *           The Grötzsch graph is a triangle-free graph with
 *           11 vertices, 20 edges, and chromatic number 4. It is named after
 *           German mathematician Herbert Grötzsch, and its existence
 *           demonstrates that the assumption of planarity is necessary in
 *           Grötzsch's theorem that every triangle-free planar
 *           graph is 3-colorable.
 *   \cli Heawood
 *           The Heawood graph is an undirected graph with 14
 *           vertices and 21 edges. The graph is cubic, and all cycles in the
 *           graph have six or more edges. Every smaller cubic graph has shorter
 *           cycles, so this graph is the 6-cage, the smallest cubic graph of
 *           girth 6.
 *   \cli Herschel
 *           The Herschel graph is the smallest
 *           nonhamiltonian polyhedral graph. It is the
 *           unique such graph on 11 nodes, and has 18 edges.
 *   \cli House
 *           The house graph is a 5-vertex, 6-edge graph, the
 *           schematic draw of a house if drawn properly, basically a
 *           triangle on top of a square.
 *   \cli HouseX
 *           The same as the house graph with an X in the square. 5
 *           vertices and 8 edges.
 *   \cli Icosahedral, Icosahedron
 *           A Platonic solid with 12
 *           vertices and 30 edges.
 *   \cli Krackhardt_Kite
 *           A social network with 10 vertices and 18 edges.
 *           Krackhardt, D. Assessing the Political Landscape:
 *           Structure, Cognition, and Power in Organizations.
 *           Admin. Sci. Quart. 35, 342-369, 1990.
 *   \cli Levi
 *           The graph is a 4-arc transitive cubic graph, it has
 *           30 vertices and 45 edges.
 *   \cli McGee
 *           The McGee graph is the unique 3-regular 7-cage
 *           graph, it has 24 vertices and 36 edges.
 *   \cli Meredith
 *           The Meredith graph is a quartic graph on 70
 *           nodes and 140 edges that is a counterexample to the conjecture that
 *           every 4-regular 4-connected graph is Hamiltonian.
 *   \cli Noperfectmatching
 *           A connected graph with 16 vertices and
 *           27 edges containing no perfect matching. A matching in a graph
 *           is a set of pairwise non-incident edges; that is, no two edges
 *           share a common vertex. A perfect matching is a matching
 *           which covers all vertices of the graph.
 *   \cli Nonline
 *           A graph whose connected components are the 9
 *           graphs whose presence as a vertex-induced subgraph in a
 *           graph makes a nonline graph. It has 50 vertices and 72 edges.
 *   \cli Octahedral, Octahedron
 *           Platonic solid with 6
 *           vertices and 12 edges.
 *   \cli Petersen
 *           A 3-regular graph with 10 vertices and 15 edges. It is
 *           the smallest hypohamiltonian graph, i.e. it is
 *           non-hamiltonian but removing any single vertex from it makes it
 *           Hamiltonian.
 *   \cli Robertson
 *           The unique (4,5)-cage graph, i.e. a 4-regular
 *           graph of girth 5. It has 19 vertices and 38 edges.
 *   \cli Smallestcyclicgroup
 *           A smallest nontrivial graph
 *           whose automorphism group is cyclic. It has 9 vertices and
 *           15 edges.
 *   \cli Tetrahedral, Tetrahedron
 *           Platonic solid with 4
 *           vertices and 6 edges.
 *   \cli Thomassen
 *           The smallest hypotraceable graph,
 *           on 34 vertices and 52 edges. A hypotracable graph does
 *           not contain a Hamiltonian path but after removing any
 *           single vertex from it the remainder always contains a
 *           Hamiltonian path. A graph containing a Hamiltonian path
 *           is called traceable.
 *   \cli Tutte
 *           Tait's Hamiltonian graph conjecture states that
 *           every 3-connected 3-regular planar graph is Hamiltonian.
 *           This graph is a counterexample. It has 46 vertices and 69
 *           edges.
 *   \cli Uniquely3colorable
 *           Returns a 12-vertex, triangle-free
 *           graph with chromatic number 3 that is uniquely
 *           3-colorable.
 *   \cli Walther
 *           An identity graph with 25 vertices and 31
 *           edges. An identity graph has a single graph automorphism,
 *           the trivial one.
 *   \cli Zachary
 *           Social network of friendships between 34 members of a
 *           karate club at a US university in the 1970s. See
 *           W. W. Zachary, An information flow model for conflict and
 *           fission in small groups, Journal of Anthropological
 *           Research 33, 452-473 (1977).
 * \endclist
 *
 * \param graph Pointer to an uninitialized graph object.
 * \param name Character constant, the name of the graph to be
 *     created, it is case insensitive.
 * \return Error code, \c IGRAPH_EINVAL if there is no graph with the
 *     given name.
 *
 * \sa Other functions for creating graph structures:
 * \ref igraph_ring(), \ref igraph_tree(), \ref igraph_lattice(), \ref
 * igraph_full().
 *
 * Time complexity: O(|V|+|E|), the number of vertices plus the number
 * of edges in the graph.
 */

int igraph_famous(igraph_t *graph, const char *name) {

    if (!strcasecmp(name, "bull")) {
        return igraph_i_famous(graph, igraph_i_famous_bull);
    } else if (!strcasecmp(name, "chvatal")) {
        return igraph_i_famous(graph, igraph_i_famous_chvatal);
    } else if (!strcasecmp(name, "coxeter")) {
        return igraph_i_famous(graph, igraph_i_famous_coxeter);
    } else if (!strcasecmp(name, "cubical")) {
        return igraph_i_famous(graph, igraph_i_famous_cubical);
    } else if (!strcasecmp(name, "diamond")) {
        return igraph_i_famous(graph, igraph_i_famous_diamond);
    } else if (!strcasecmp(name, "dodecahedral") ||
               !strcasecmp(name, "dodecahedron")) {
        return igraph_i_famous(graph, igraph_i_famous_dodecahedron);
    } else if (!strcasecmp(name, "folkman")) {
        return igraph_i_famous(graph, igraph_i_famous_folkman);
    } else if (!strcasecmp(name, "franklin")) {
        return igraph_i_famous(graph, igraph_i_famous_franklin);
    } else if (!strcasecmp(name, "frucht")) {
        return igraph_i_famous(graph, igraph_i_famous_frucht);
    } else if (!strcasecmp(name, "grotzsch")) {
        return igraph_i_famous(graph, igraph_i_famous_grotzsch);
    } else if (!strcasecmp(name, "heawood")) {
        return igraph_i_famous(graph, igraph_i_famous_heawood);
    } else if (!strcasecmp(name, "herschel")) {
        return igraph_i_famous(graph, igraph_i_famous_herschel);
    } else if (!strcasecmp(name, "house")) {
        return igraph_i_famous(graph, igraph_i_famous_house);
    } else if (!strcasecmp(name, "housex")) {
        return igraph_i_famous(graph, igraph_i_famous_housex);
    } else if (!strcasecmp(name, "icosahedral") ||
               !strcasecmp(name, "icosahedron")) {
        return igraph_i_famous(graph, igraph_i_famous_icosahedron);
    } else if (!strcasecmp(name, "krackhardt_kite")) {
        return igraph_i_famous(graph, igraph_i_famous_krackhardt_kite);
    } else if (!strcasecmp(name, "levi")) {
        return igraph_i_famous(graph, igraph_i_famous_levi);
    } else if (!strcasecmp(name, "mcgee")) {
        return igraph_i_famous(graph, igraph_i_famous_mcgee);
    } else if (!strcasecmp(name, "meredith")) {
        return igraph_i_famous(graph, igraph_i_famous_meredith);
    } else if (!strcasecmp(name, "noperfectmatching")) {
        return igraph_i_famous(graph, igraph_i_famous_noperfectmatching);
    } else if (!strcasecmp(name, "nonline")) {
        return igraph_i_famous(graph, igraph_i_famous_nonline);
    } else if (!strcasecmp(name, "octahedral") ||
               !strcasecmp(name, "octahedron")) {
        return igraph_i_famous(graph, igraph_i_famous_octahedron);
    } else if (!strcasecmp(name, "petersen")) {
        return igraph_i_famous(graph, igraph_i_famous_petersen);
    } else if (!strcasecmp(name, "robertson")) {
        return igraph_i_famous(graph, igraph_i_famous_robertson);
    } else if (!strcasecmp(name, "smallestcyclicgroup")) {
        return igraph_i_famous(graph, igraph_i_famous_smallestcyclicgroup);
    } else if (!strcasecmp(name, "tetrahedral") ||
               !strcasecmp(name, "tetrahedron")) {
        return igraph_i_famous(graph, igraph_i_famous_tetrahedron);
    } else if (!strcasecmp(name, "thomassen")) {
        return igraph_i_famous(graph, igraph_i_famous_thomassen);
    } else if (!strcasecmp(name, "tutte")) {
        return igraph_i_famous(graph, igraph_i_famous_tutte);
    } else if (!strcasecmp(name, "uniquely3colorable")) {
        return igraph_i_famous(graph, igraph_i_famous_uniquely3colorable);
    } else if (!strcasecmp(name, "walther")) {
        return igraph_i_famous(graph, igraph_i_famous_walther);
    } else if (!strcasecmp(name, "zachary")) {
        return igraph_i_famous(graph, igraph_i_famous_zachary);
    }

    IGRAPH_ERRORF("%s is not a known graph. See the documentation for valid graph names.", IGRAPH_EINVAL, name);
}
