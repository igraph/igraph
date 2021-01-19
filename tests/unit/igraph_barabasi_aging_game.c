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
#include <assert.h>
#include "test_utilities.inc"

int main() {
    igraph_t g;
    igraph_rng_seed(igraph_rng_default(), 42);
    igraph_bool_t tree;

//    /* no vertices */
//    assert(igraph_barabasi_aging_game(
//        &g, /*nodes*/ 0, /*m: edges_per_step*/ 1, 
//        /*outseq: edges per step as vector*/ 0, /*outpref*/ 0,
//        /*pa_exp*/ 1, /*aging_exp*/ 1, /*aging_bin*/ 1,
//        /*zero_deg_appeal*/ 0, /*zero_age_appeal*/ 0, /*deg_coef*/ 0,
//        /*age_coef */ 0, /*directed*/ 0) == IGRAPH_SUCCESS);
//    assert(igraph_vcount(&g) == 0);
//    assert(!igraph_is_directed(&g));
//    igraph_destroy(&g);
//
//    /* one edge per step makes a tree */
//    assert(igraph_barabasi_aging_game(
//           &g, /*nodes*/ 10, /*m: edges_per_step*/ 1, 
//           /*outseq: edges per step as vector*/ 0, /*outpref*/ 0,
//           /*pa_exp*/ 0.5, /*aging_exp*/ -0.5, /*aging_bin*/ 2,
//           /*zero_deg_appeal*/ 0.1, /*zero_age_appeal*/ 0, /*deg_coef*/ 0.1,
//           /*age_coef */ 0.1, /*directed*/ 1) == IGRAPH_SUCCESS);
//    assert(igraph_vcount(&g) == 10);
//    assert(igraph_ecount(&g) == 9);
//    assert(igraph_is_directed(&g));
//    igraph_is_tree(&g, &tree, /* root*/ NULL, IGRAPH_IN);
//    assert(tree);
//    igraph_destroy(&g);
//
//    /*Prefer old vertices to make a star of triple edges */
//    assert(igraph_barabasi_aging_game(
//        &g, /*nodes*/ 5, /*m: edges_per_step*/ 3, 
//        /*outseq: edges per step as vector*/ 0, /*outpref*/ 0,
//        /*pa_exp*/ 0.0, /*aging_exp*/ 100, /*aging_bin*/ 1,
//        /*zero_deg_appeal*/ 1.0, /*zero_age_appeal*/ 0, /*deg_coef*/ 0.0,
//        /*age_coef */ 1, /*directed*/ 1) == IGRAPH_SUCCESS);
//    print_graph_canon(&g);
//    igraph_destroy(&g);
//
//    /* Prefer new vertices to make a line of double edges*/
//    assert(igraph_barabasi_aging_game(
//        &g, /*nodes*/ 5, /*m: edges_per_step*/ 2, 
//        /*outseq: edges per step as vector*/ 0, /*outpref*/ 0,
//        /*pa_exp*/ 0.0, /*aging_exp*/ 100, /*aging_bin*/ 1,
//        /*zero_deg_appeal*/ 1.0, /*zero_age_appeal*/ 1, /*deg_coef*/ 0.0,
//        /*age_coef */ -1, /*directed*/ 1) == IGRAPH_SUCCESS);
//    print_graph_canon(&g);
//    igraph_destroy(&g);
//
//    /* Prefer more edges to make a star*/
//    assert(igraph_barabasi_aging_game(
//        &g, /*nodes*/ 5, /*m: edges_per_step*/ 2, 
//        /*outseq: edges per step as vector*/ 0, /*outpref*/ 0,
//        /*pa_exp*/ 100.0, /*aging_exp*/ 0.0, /*aging_bin*/ 1,
//        /*zero_deg_appeal*/ 1.0, /*zero_age_appeal*/ 1.0, /*deg_coef*/ 1.0,
//        /*age_coef */ 0.0, /*directed*/ 1) == IGRAPH_SUCCESS);
//    print_graph_canon(&g);
//    igraph_destroy(&g);
//
    /* Prefer less edges to make a line*/
    assert(igraph_barabasi_aging_game(
        &g, /*nodes*/ 5, /*m: edges_per_step*/ 2, 
        /*outseq: edges per step as vector*/ 0, /*outpref*/ 0,
        /*pa_exp*/ 1.0, /*aging_exp*/ 0.0, /*aging_bin*/ 1,
        /*zero_deg_appeal*/ 1.0, /*zero_age_appeal*/ 1.0, /*deg_coef*/ -1.0,
        /*age_coef */ 0.0, /*directed*/ 1) == IGRAPH_SUCCESS);
    print_graph_canon(&g);
    //igraph_destroy(&g);

    VERIFY_FINALLY_STACK();
    return 0;

}
