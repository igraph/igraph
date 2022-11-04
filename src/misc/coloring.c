/*
  Heuristic graph coloring algorithms.
  Copyright (C) 2017 Szabolcs Horvat <szhorvat@gmail.com>

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

#include "igraph_coloring.h"

#include "igraph_interface.h"
#include "igraph_adjlist.h"

#include "core/indheap.h"
#include "core/interruption.h"

static igraph_error_t igraph_i_vertex_coloring_greedy_cn(const igraph_t *graph, igraph_vector_int_t *colors) {
    igraph_integer_t i, vertex, maxdeg;
    igraph_integer_t vc = igraph_vcount(graph);
    igraph_2wheap_t cn; /* indexed heap storing number of already coloured neighbours */
    igraph_vector_int_t neigh_colors;
    igraph_adjlist_t adjlist;

    IGRAPH_CHECK(igraph_vector_int_resize(colors, vc));
    igraph_vector_int_fill(colors, 0);

    /* Nothing to do for 0 or 1 vertices.
     * Remember that colours are integers starting from 0,
     * and the 'colors' vector is already 0-initialized above.
     */
    if (vc <= 1) {
        return IGRAPH_SUCCESS;
    }

    IGRAPH_CHECK(igraph_adjlist_init(graph, &adjlist, IGRAPH_ALL, IGRAPH_LOOPS_TWICE, IGRAPH_MULTIPLE));
    IGRAPH_FINALLY(igraph_adjlist_destroy, &adjlist);

    /* find maximum degree and a corresponding vertex */
    {
        igraph_vector_int_t degree;

        IGRAPH_CHECK(igraph_vector_int_init(&degree, 0));
        IGRAPH_FINALLY(igraph_vector_int_destroy, &degree);
        IGRAPH_CHECK(igraph_degree(graph, &degree, igraph_vss_all(), IGRAPH_ALL, 0));

        vertex = igraph_vector_int_which_max(&degree);
        maxdeg = VECTOR(degree)[vertex];

        igraph_vector_int_destroy(&degree);
        IGRAPH_FINALLY_CLEAN(1);
    }

    IGRAPH_CHECK(igraph_vector_int_init(&neigh_colors, 0));
    IGRAPH_CHECK(igraph_vector_int_reserve(&neigh_colors, maxdeg));
    IGRAPH_FINALLY(igraph_vector_int_destroy, &neigh_colors);

    IGRAPH_CHECK(igraph_2wheap_init(&cn, vc));
    IGRAPH_FINALLY(igraph_2wheap_destroy, &cn);
    for (i = 0; i < vc; ++i)
        if (i != vertex) {
            igraph_2wheap_push_with_index(&cn, i, 0); /* should not fail since memory was already reserved */
        }

    while (1) {
        igraph_vector_int_t *neighbors = igraph_adjlist_get(&adjlist, vertex);
        igraph_integer_t neigh_count = igraph_vector_int_size(neighbors);

        /* colour current vertex */
        {
            igraph_integer_t col;

            IGRAPH_CHECK(igraph_vector_int_resize(&neigh_colors, neigh_count));
            for (i = 0; i < neigh_count; ++i) {
                VECTOR(neigh_colors)[i] = VECTOR(*colors)[ VECTOR(*neighbors)[i] ];
            }
            igraph_vector_int_sort(&neigh_colors);

            i = 0;
            col = 0;
            do {
                while (i < neigh_count && VECTOR(neigh_colors)[i] == col) {
                    i++;
                }
                col++;
            } while (i < neigh_count && VECTOR(neigh_colors)[i] == col);

            VECTOR(*colors)[vertex] = col;
        }

        /* increment number of coloured neighbours for each neighbour of vertex */
        for (i = 0; i < neigh_count; ++i) {
            igraph_integer_t idx = VECTOR(*neighbors)[i];
            if (igraph_2wheap_has_elem(&cn, idx)) {
                igraph_2wheap_modify(&cn, idx, igraph_2wheap_get(&cn, idx) + 1);
            }
        }

        /* stop if no more vertices left to colour */
        if (igraph_2wheap_empty(&cn)) {
            break;
        }

        igraph_2wheap_delete_max_index(&cn, &vertex);

        IGRAPH_ALLOW_INTERRUPTION();
    }

    /* subtract 1 from each colour value, so that colours start at 0 */
    igraph_vector_int_add_constant(colors, -1);

    /* free data structures */
    igraph_vector_int_destroy(&neigh_colors);
    igraph_adjlist_destroy(&adjlist);
    igraph_2wheap_destroy(&cn);
    IGRAPH_FINALLY_CLEAN(3);

    return IGRAPH_SUCCESS;
}


/**
 * \function igraph_vertex_coloring_greedy
 * \brief Computes a vertex coloring using a greedy algorithm.
 *
 * </para><para>
 * This function assigns a "color"—represented as a non-negative integer—to
 * each vertex of the graph in such a way that neighboring vertices never have
 * the same color. The obtained coloring is not necessarily minimal.
 *
 * </para><para>
 * Vertices are colored one by one, choosing the smallest color index that
 * differs from that of already colored neighbors.
 * Colors are represented with non-negative integers 0, 1, 2, ...
 *
 * \param graph The input graph.
 * \param colors Pointer to an initialized integer vector. The vertex colors will be stored here.
 * \param heuristic The vertex ordering heuristic to use during greedy coloring. See \ref igraph_coloring_greedy_t
 *
 * \return Error code.
 *
 * \example examples/simple/igraph_coloring.c
 */
igraph_error_t igraph_vertex_coloring_greedy(const igraph_t *graph, igraph_vector_int_t *colors, igraph_coloring_greedy_t heuristic) {
    switch (heuristic) {
    case IGRAPH_COLORING_GREEDY_COLORED_NEIGHBORS:
        return igraph_i_vertex_coloring_greedy_cn(graph, colors);
    default:
        IGRAPH_ERROR("Invalid heuristic for greedy vertex coloring.", IGRAPH_EINVAL);
    }
}

static igraph_integer_t DSatur_node_selection(const igraph_t *graph,igraph_vector_int_t *colors, igraph_vector_int_t * saturation_degree, igraph_integer_t vc){
    igraph_integer_t mst_saurated_node = -1, max_saturation = -1, max_degree = -1;
    for( igraph_integer_t node = 0 ; node < vc ; node++ ){
        //finding an uncolored node with max (saturation degree ,  degree)
        if( VECTOR(*colors)[node] == -1 ){
            if( VECTOR(*saturation_degree)[node] > max_saturation ){
                max_saturation = VECTOR(*saturation_degree)[node];
                mst_saurated_node = node;
            }
            else if( VECTOR(*saturation_degree)[node] == max_saturation ){
                igraph_integer_t degree;
                igraph_degree_1(graph, &degree, node, IGRAPH_ALL, true);
                if( degree > max_degree ){
                    mst_saurated_node = node;
                    max_degree = degree;
                }
            }
        }
    }
    return mst_saurated_node;
}

static void update_saturation_degree(const igraph_t *graph, igraph_vector_int_t * saturation_degree,igraph_vector_int_t *neighbours){
    igraph_integer_t nbr_cnt = igraph_vector_int_size(neighbours);
    for( igraph_integer_t nmbr_indx = 0 ; nmbr_indx < nbr_cnt ; nmbr_indx++ ){
        igraph_integer_t nbr = VECTOR(*neighbours)[nmbr_indx];
        VECTOR(*saturation_degree)[ nbr ] += 1;
    }
}

static igraph_error_t DSatur(const igraph_t *graph,igraph_vector_int_t *colors, igraph_vector_int_t * saturation_degree){
    igraph_vector_int_fill(colors, -1); //-1 as a color means uncolored
    igraph_vector_int_fill(saturation_degree, 0);

    igraph_integer_t vc = igraph_vcount(graph);
    igraph_integer_t vertices_colored = 0;
    igraph_vector_int_t neighbours;
    igraph_vector_int_init(&neighbours, 0);
    while( vertices_colored < vc ){
        igraph_integer_t node_to_color = DSatur_node_selection(graph, colors, saturation_degree, vc);

        IGRAPH_CHECK(igraph_neighbors(graph, &neighbours, node_to_color, IGRAPH_ALL ));
        igraph_integer_t nbr_cnt = igraph_vector_int_size(&neighbours);
        for( igraph_integer_t color = 0 ; color < vc ; color++ ){
            igraph_bool_t viable_color = true;
            for( igraph_integer_t nmbr_indx = 0 ; nmbr_indx < nbr_cnt ; nmbr_indx++ ){
                igraph_integer_t nbr = VECTOR(neighbours)[nmbr_indx];
                if( color == VECTOR(*colors)[nbr] ){
                    viable_color = false;
                    break;
                }
            }
            if( viable_color ){
                VECTOR(*colors)[ node_to_color ] = color;
                break;
            }
        }
        update_saturation_degree(graph, saturation_degree, &neighbours);

        vertices_colored++;

    }
    return IGRAPH_SUCCESS;
}

// static igraph_error_t Color(const igraph_t *graph,igraph_vector_int_t *colors, igraph_vector_int_t * saturation_degree, igraph_matrix_int_t* ants,
//     igraph_integer_t nmbr_colors){

//         return IGRAPH_SUCCESS;
// }

/**
 * \function igraph_vertex_coloring_antColony
 * \brief Computes a vertex coloring using ant colony optimisation metaheuristic.
 *
 * </para><para>
 * This function assigns a "color"—represented as a non-negative integer—to
 * each vertex of the graph in such a way that neighboring vertices never have
 * the same color. The obtained coloring is not necessarily minimal.
 *
 * </para><para>
 * The coloring is performed using the ant-colony optimisation metaheuristic as described in the paper with slight variations
 * A New Ant Algorithm for Graph Coloring - Alain Hertz, Nicolas Zufferey
 * </para><para>
 * The variations are
 * 1) The sorting of the move is first done by how much it would reduce conflicts, and then inplace by the fitness formula given in the paper
 * 2) The selection of moves is not deterministic its probabalistic with the probability distribution being a normal distribution
 * Both of these can be turned off (see parameters) to use the original algorithm
 * \param graph The input graph.
 * \param colors Pointer to an initialized integer vector. The vertex colors will be stored here.
 * \param alpha real number which describes the greedy force when calculating move fitness. Paper suggests keeping this 1.0
 * \param beta  real number which describes the trail force when calculating move fitness. Paper suggests keeping this 5.0
 * \param rho evaporation constant which describes the trail being forgotten per iteration. Paper suggests keeping this at or close to 0.9
 * \param maxIters number of iterations the ant colony loop will run for each color. Paper suggest to keep that at 1000. Setting this to 0 will
 * will jus return the graph colored using the D-Satur algorithm and all checks for other parameters being valid will not run
 * \param tabuRangeUp the lower bound on the maximum tabu tenure. Paper suggests keeping this 9
 * \param tabuFactor the formula for tabu tenure of a move is given by random_int(0, tabuRange) + tabuFactor * NCV(s) where s is the current coloring
 * and NVC is the number of conflicting vertices.
 * \param k_end the number of colors which once achieved will make the algorithm return. Setting this to 0 will allow the algorithm to try and find the optimal coloring
 * \param firstVariationFlag setting this to true will enable the first variatio, which has been observed to improve results
 * \param secondVariationFlag setting this to true will enable second variation, which works well with the first variation.
 * \return Error code.
 *
 * \example examples/simple/igraph_coloring.c
 */
igraph_error_t igraph_vertex_coloring_AntColony(const igraph_t *graph, igraph_vector_int_t *colors,igraph_real_t alpha, igraph_real_t beta,
    igraph_real_t rho, igraph_integer_t maxIters, igraph_integer_t tabuRangeUp, igraph_real_t tabuFactor, igraph_integer_t k_end , igraph_bool_t firstVariationFlag,
    igraph_bool_t secondVariationFlag ) {

    igraph_integer_t vc = igraph_vcount(graph);
    IGRAPH_CHECK(igraph_vector_int_resize(colors, vc));
    if (vc <= 1) {
        igraph_vector_int_fill(colors, 0);
        return IGRAPH_SUCCESS;
    }

    igraph_vector_int_t saturation_degree;
    IGRAPH_CHECK(igraph_vector_int_init(&saturation_degree, vc));
    if( maxIters == 0 ){
        IGRAPH_CHECK(DSatur(graph, colors, &saturation_degree));
        return IGRAPH_SUCCESS;
    }

    if( rho < 0 || rho > 1 ){
        IGRAPH_ERROR("Evaporation constant should be b/w 0 and 1", IGRAPH_EINVAL);
    }
    if( maxIters < 0 ){
        IGRAPH_ERROR("maxIters should be greater than or equal to 0", IGRAPH_EINVAL);
    }
    if( tabuRangeUp < 0 ){
        IGRAPH_ERROR("tabuRangeUp should be greater than or equal to 0", IGRAPH_EINVAL);
    }
    if( tabuFactor < 0.0 ){
        IGRAPH_ERROR("tabuFactor should be greater than or equal to 0", IGRAPH_EINVAL);
    }
    IGRAPH_CHECK(DSatur(graph, colors, &saturation_degree));
    return IGRAPH_SUCCESS;
}
