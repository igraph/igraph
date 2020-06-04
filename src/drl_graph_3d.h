/*
 * Copyright 2007 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
 * certain rights in this software.
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 *     * Redistributions of source code must retain the above copyright
 * notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 * notice, this list of conditions and the following disclaimer in the
 * documentation and/or other materials provided with the distribution.
 *     * Neither the name of Sandia National Laboratories nor the names of
 * its contributors may be used to endorse or promote products derived from
 * this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
 * TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */
// The graph class contains the methods necessary to draw the
// graph.  It calls on the density server class to obtain
// position and density information

#include "DensityGrid_3d.h"
#include "igraph_layout.h"

#include <map>
#include <vector>
#include <ctime>

namespace drl3d {

// layout schedule information
struct layout_schedule {
    int iterations;
    float temperature;
    float attraction;
    float damping_mult;
    time_t time_elapsed;
};

class graph {

public:

    // Methods
    void init_parms ( int rand_seed, float edge_cut, float real_parm );
    void init_parms ( const igraph_layout_drl_options_t *options );
    int read_real ( const igraph_matrix_t *real_mat,
                    const igraph_vector_bool_t *fixed);
    int draw_graph (igraph_matrix_t *res);
    float get_tot_energy ( );

    // Con/Decon
    graph( const igraph_t *igraph,
           const igraph_layout_drl_options_t *options,
           const igraph_vector_t *weights);
    ~graph( ) { }

private:

    // Methods
    int ReCompute ( );
    void update_nodes ( );
    float Compute_Node_Energy ( int node_ind );
    void Solve_Analytic ( int node_ind, float &pos_x, float &pos_y, float &pos_z );
    void get_positions ( std::vector<int> &node_indices, float return_positions[3 * MAX_PROCS] );
    void update_density ( std::vector<int> &node_indices,
                          float old_positions[3 * MAX_PROCS],
                          float new_positions[3 * MAX_PROCS] );
    void update_node_pos ( int node_ind,
                           float old_positions[3 * MAX_PROCS],
                           float new_positions[3 * MAX_PROCS] );

    // MPI information
    int myid, num_procs;

    // graph decomposition information
    int num_nodes;                  // number of nodes in graph
    float highest_sim;              // highest sim for normalization
    std::map <int, int> id_catalog;      // id_catalog[file id] = internal id
    std::map <int, std::map <int, float> > neighbors;     // neighbors of nodes on this proc.

    // graph layout information
    std::vector<Node> positions;
    DensityGrid density_server;

    // original VxOrd information
    int STAGE, iterations;
    float temperature, attraction, damping_mult;
    float min_edges, CUT_END, cut_length_end, cut_off_length, cut_rate;
    bool first_add, fine_first_add, fineDensity;

    // scheduling variables
    layout_schedule liquid;
    layout_schedule expansion;
    layout_schedule cooldown;
    layout_schedule crunch;
    layout_schedule simmer;

    // timing statistics
    time_t start_time, stop_time;

    // online clustering information
    int real_iterations;    // number of iterations to hold .real input fixed
    int tot_iterations;
    int tot_expected_iterations; // for progress bar
    bool real_fixed;
};

} // namespace drl3d
