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
// This file contains the member definitions of the master class

#include <map>
#include <vector>
#include <cmath>

using namespace std;

#include "drl_graph_3d.h"
#include "igraph_random.h"
#include "igraph_interface.h"
#include "igraph_progress.h"
#include "core/interruption.h"
#ifdef MUSE_MPI
    #include <mpi.h>
#endif

namespace drl3d {

graph::graph(const igraph_t *igraph,
             const igraph_layout_drl_options_t *options,
             const igraph_vector_t *weights) {
    myid = 0;
    num_procs = 1;

    STAGE = 0;
    iterations = options->init_iterations;
    temperature = options->init_temperature;
    attraction = options->init_attraction;
    damping_mult = options->init_damping_mult;
    min_edges = 20;
    first_add = fine_first_add = true;
    fineDensity = false;

    // Brian's original Vx schedule
    liquid.iterations = options->liquid_iterations;
    liquid.temperature = options->liquid_temperature;
    liquid.attraction = options->liquid_attraction;
    liquid.damping_mult = options->liquid_damping_mult;
    liquid.time_elapsed = 0;

    expansion.iterations = options->expansion_iterations;
    expansion.temperature = options->expansion_temperature;
    expansion.attraction = options->expansion_attraction;
    expansion.damping_mult = options->expansion_damping_mult;
    expansion.time_elapsed = 0;

    cooldown.iterations = options->cooldown_iterations;
    cooldown.temperature = options->cooldown_temperature;
    cooldown.attraction = options->cooldown_attraction;
    cooldown.damping_mult = options->cooldown_damping_mult;
    cooldown.time_elapsed = 0;

    crunch.iterations = options->crunch_iterations;
    crunch.temperature = options->crunch_temperature;
    crunch.attraction = options->crunch_attraction;
    crunch.damping_mult = options->crunch_damping_mult;
    crunch.time_elapsed = 0;

    simmer.iterations = options->simmer_iterations;
    simmer.temperature = options->simmer_temperature;
    simmer.attraction = options->simmer_attraction;
    simmer.damping_mult = options->simmer_damping_mult;
    simmer.time_elapsed = 0;

    // scan .int file for node info
    highest_sim = 1.0;
    num_nodes = igraph_vcount(igraph);
    long int no_of_edges = igraph_ecount(igraph);
    for (long int i = 0; i < num_nodes; i++) {
        id_catalog[i] = 1;
    }
    map< int, int>::iterator cat_iter;
    for ( cat_iter = id_catalog.begin();
          cat_iter != id_catalog.end(); cat_iter++) {
        cat_iter->second = cat_iter->first;
    }

    // populate node positions and ids
    positions.reserve ( num_nodes );
    for ( cat_iter = id_catalog.begin();
          cat_iter != id_catalog.end();
          cat_iter++ ) {
        positions.push_back ( Node( cat_iter->first ) );
    }

    // read .int file for graph info
    long int node_1, node_2;
    double weight;
    for (long int i = 0; i < no_of_edges; i++) {
        node_1 = IGRAPH_FROM(igraph, i);
        node_2 = IGRAPH_TO(igraph, i);
        weight = weights ? VECTOR(*weights)[i] : 1.0 ;
        (neighbors[id_catalog[node_1]])[id_catalog[node_2]] = weight;
        (neighbors[id_catalog[node_2]])[id_catalog[node_1]] = weight;
    }

    // initialize density server
    density_server.Init();

}

// init_parms -- this subroutine initializes the edge_cut variables
// used in the original VxOrd starting with the edge_cut parameter.
// In our version, edge_cut = 0 means no cutting, 1 = maximum cut.
// We also set the random seed here.

void graph::init_parms ( int rand_seed, float edge_cut, float real_parm ) {

    IGRAPH_UNUSED(rand_seed);
    // first we translate edge_cut the former tcl sliding scale
    //CUT_END = cut_length_end = 39000.0 * (1.0 - edge_cut) + 1000.0;
    CUT_END = cut_length_end = 40000.0 * (1.0 - edge_cut);

    // cut_length_end cannot actually be 0
    if ( cut_length_end <= 1.0 ) {
        cut_length_end = 1.0;
    }

    float cut_length_start = 4.0 * cut_length_end;

    // now we set the parameters used by ReCompute
    cut_off_length = cut_length_start;
    cut_rate = ( cut_length_start - cut_length_end ) / 400.0;

    // finally set the number of iterations to leave .real coords fixed
    int full_comp_iters;
    full_comp_iters = liquid.iterations + expansion.iterations +
                      cooldown.iterations + crunch.iterations + 3;

    // adjust real parm to iterations (do not enter simmer halfway)
    if ( real_parm < 0 ) {
        real_iterations = (int)real_parm;
    } else if ( real_parm == 1) {
        real_iterations = full_comp_iters + simmer.iterations + 100;
    } else {
        real_iterations = (int)(real_parm * full_comp_iters);
    }

    tot_iterations = 0;
    if ( real_iterations > 0 ) {
        real_fixed = true;
    } else {
        real_fixed = false;
    }

    // calculate total expected iterations (for progress bar display)
    tot_expected_iterations = liquid.iterations +
                              expansion.iterations + cooldown.iterations +
                              crunch.iterations + simmer.iterations;

    /*
    // output edge_cutting parms (for debugging)
    cout << "Processor " << myid << ": "
         << "cut_length_end = CUT_END = " << cut_length_end
         << ", cut_length_start = " << cut_length_start
         << ", cut_rate = " << cut_rate << endl;
    */

    // set random seed
    // srand ( rand_seed ); // Don't need this in igraph

}

void graph::init_parms(const igraph_layout_drl_options_t *options) {
    double rand_seed = 0.0;
    double real_in = -1.0;
    init_parms(rand_seed, options->edge_cut, real_in);
}

int graph::read_real ( const igraph_matrix_t *real_mat,
                       const igraph_vector_bool_t *fixed) {
    long int n = igraph_matrix_nrow(real_mat);
    for (long int i = 0; i < n; i++) {
        positions[id_catalog[i]].x = MATRIX(*real_mat, i, 0);
        positions[id_catalog[i]].y = MATRIX(*real_mat, i, 1);
        positions[id_catalog[i]].z = MATRIX(*real_mat, i, 2);
        positions[id_catalog[i]].fixed = fixed ? VECTOR(*fixed)[i] : false;

        if ( real_iterations > 0 ) {
            density_server.Add ( positions[id_catalog[i]], fineDensity );
        }
    }

    return 0;
}

/*********************************************
 * Function: ReCompute                       *
 * Description: Compute the graph locations  *
 * Modified from original code by B. Wylie   *
 ********************************************/

int graph::ReCompute( ) {

    // carryover from original VxOrd
    int MIN = 1;

    /*
    // output parameters (for debugging)
    cout << "ReCompute is using the following parameters: "<< endl;
    cout << "STAGE: " << STAGE << ", iter: " << iterations << ", temp = " << temperature
         << ", attract = " << attraction << ", damping_mult = " << damping_mult
       << ", min_edges = " << min_edges << ", cut_off_length = " << cut_off_length
       << ", fineDensity = " << fineDensity << endl;
    */

    /* igraph progress report */
    float progress = (tot_iterations * 100.0 / tot_expected_iterations);

    switch (STAGE) {
    case 0:
        if (iterations == 0) {
            IGRAPH_PROGRESS("DrL layout (initialization stage)", progress, 0);
        } else {
            IGRAPH_PROGRESS("DrL layout (liquid stage)", progress, 0);
        }
        break;
    case 1:
        IGRAPH_PROGRESS("DrL layout (expansion stage)", progress, 0); break;
    case 2:
        IGRAPH_PROGRESS("DrL layout (cooldown and cluster phase)", progress, 0); break;
    case 3:
        IGRAPH_PROGRESS("DrL layout (crunch phase)", progress, 0); break;
    case 5:
        IGRAPH_PROGRESS("DrL layout (simmer phase)", progress, 0); break;
    case 6:
        IGRAPH_PROGRESS("DrL layout (final phase)", 100.0, 0); break;
    default:
        IGRAPH_PROGRESS("DrL layout (unknown phase)", 0.0, 0); break;
    }

    /* Compute Energies for individual nodes */
    update_nodes ();

    // check to see if we need to free fixed nodes
    tot_iterations++;
    if ( tot_iterations >= real_iterations ) {
        real_fixed = false;
    }


    // ****************************************
    // AUTOMATIC CONTROL SECTION
    // ****************************************

    // STAGE 0: LIQUID
    if (STAGE == 0) {

        if ( iterations == 0 ) {
            start_time = time( NULL );
//          if ( myid == 0 )
//              cout << "Entering liquid stage ...";
        }

        if (iterations < liquid.iterations) {
            temperature = liquid.temperature;
            attraction = liquid.attraction;
            damping_mult = liquid.damping_mult;
            iterations++;
//          if ( myid == 0 )
//              cout << "." << flush;

        } else {

            stop_time = time( NULL );
            liquid.time_elapsed = liquid.time_elapsed + (stop_time - start_time);
            temperature = expansion.temperature;
            attraction = expansion.attraction;
            damping_mult = expansion.damping_mult;
            iterations = 0;

            // go to next stage
            STAGE = 1;
            start_time = time( NULL );

//          if ( myid == 0 )
//              cout << "Entering expansion stage ...";
        }
    }

    // STAGE 1: EXPANSION
    if (STAGE == 1) {

        if (iterations < expansion.iterations) {

            // Play with vars
            if (attraction > 1) {
                attraction -= .05f;
            }
            if (min_edges > 12) {
                min_edges -= .05f;
            }
            cut_off_length -= cut_rate;
            if (damping_mult > .1) {
                damping_mult -= .005f;
            }
            iterations++;
//          if ( myid == 0 ) cout << "." << flush;

        } else {

            stop_time = time( NULL );
            expansion.time_elapsed = expansion.time_elapsed + (stop_time - start_time);
            min_edges = 12;
            damping_mult = cooldown.damping_mult;

            STAGE = 2;
            attraction = cooldown.attraction;
            temperature = cooldown.temperature;
            iterations = 0;
            start_time = time( NULL );

//          if ( myid == 0 )
//              cout << "Entering cool-down stage ...";
        }
    }

    // STAGE 2: Cool down and cluster
    else if (STAGE == 2) {

        if (iterations < cooldown.iterations) {

            // Reduce temperature
            if (temperature > 50) {
                temperature -= 10;
            }

            // Reduce cut length
            if (cut_off_length > cut_length_end) {
                cut_off_length -= cut_rate * 2;
            }
            if (min_edges > MIN) {
                min_edges -= .2f;
            }
            //min_edges = 99;
            iterations++;
//          if ( myid == 0 )
//              cout << "." << flush;

        } else {

            stop_time = time( NULL );
            cooldown.time_elapsed = cooldown.time_elapsed + (stop_time - start_time);
            cut_off_length = cut_length_end;
            temperature = crunch.temperature;
            damping_mult = crunch.damping_mult;
            min_edges = MIN;
            //min_edges = 99; // In other words: no more cutting

            STAGE = 3;
            iterations = 0;
            attraction = crunch.attraction;
            start_time = time( NULL );

//          if ( myid == 0 )
//              cout << "Entering crunch stage ...";
        }
    }

    // STAGE 3: Crunch
    else if (STAGE == 3) {

        if (iterations < crunch.iterations) {
            iterations++;
//          if ( myid == 0 ) cout << "." << flush;
        } else {

            stop_time = time( NULL );
            crunch.time_elapsed = crunch.time_elapsed + (stop_time - start_time);
            iterations = 0;
            temperature = simmer.temperature;
            attraction = simmer.attraction;
            damping_mult = simmer.damping_mult;
            min_edges = 99;
            fineDensity = true;

            STAGE = 5;
            start_time = time( NULL );

//          if ( myid == 0 )
//              cout << "Entering simmer stage ...";
        }
    }

    // STAGE 5: Simmer
    else if ( STAGE == 5 ) {

        if (iterations < simmer.iterations) {
            if (temperature > 50) {
                temperature -= 2;
            }
            iterations++;
//          if ( myid == 0 ) cout << "." << flush;
        } else {
            stop_time = time( NULL );
            simmer.time_elapsed = simmer.time_elapsed + (stop_time - start_time);

            STAGE = 6;

//          if ( myid == 0 )
//              cout << "Layout calculation completed in " <<
//                ( liquid.time_elapsed + expansion.time_elapsed +
//                  cooldown.time_elapsed + crunch.time_elapsed +
//                  simmer.time_elapsed )
//                   << " seconds (not including I/O)."
//                   << endl;
        }
    }

    // STAGE 6: All Done!
    else if ( STAGE == 6) {

        /*
        // output parameters (for debugging)
        cout << "ReCompute is using the following parameters: "<< endl;
        cout << "STAGE: " << STAGE << ", iter: " << iterations << ", temp = " << temperature
             << ", attract = " << attraction << ", damping_mult = " << damping_mult
             << ", min_edges = " << min_edges << ", cut_off_length = " << cut_off_length
             << ", fineDensity = " << fineDensity << endl;
        */

        return 0;
    }

    // ****************************************
    // END AUTOMATIC CONTROL SECTION
    // ****************************************

    // Still need more recomputation
    return 1;

}

// update_nodes -- this function will complete the primary node update
// loop in layout's recompute routine.  It follows exactly the same
// sequence to ensure similarity of parallel layout to the standard layout

void graph::update_nodes ( ) {

    vector<int> node_indices;           // node list of nodes currently being updated
    float old_positions[2 * MAX_PROCS]; // positions before update
    float new_positions[2 * MAX_PROCS]; // positions after update

    bool all_fixed;                     // check if all nodes are fixed

    // initial node list consists of 0,1,...,num_procs
    for ( int i = 0; i < num_procs; i++ ) {
        node_indices.push_back( i );
    }

    // next we calculate the number of nodes there would be if the
    // num_nodes by num_procs schedule grid were perfectly square
    int square_num_nodes = (int)(num_procs + num_procs * floor ((float)(num_nodes - 1) / (float)num_procs ));

    for ( int i = myid; i < square_num_nodes; i += num_procs ) {

        // get old positions
        get_positions ( node_indices, old_positions );

        // default new position is old position
        get_positions ( node_indices, new_positions );

        if ( i < num_nodes ) {

            // advance random sequence according to myid
            for ( int j = 0; j < 2 * myid; j++ ) {
                RNG_UNIF01();
            }
            // rand();

            // calculate node energy possibilities
            if ( !(positions[i].fixed && real_fixed) ) {
                update_node_pos ( i, old_positions, new_positions );
            }

            // advance random sequence for next iteration
            for ( unsigned int j = 2 * myid; j < 2 * (node_indices.size() - 1); j++ ) {
                RNG_UNIF01();
            }
            // rand();

        } else {
            // advance random sequence according to use by
            // the other processors
            for ( unsigned int j = 0; j < 2 * (node_indices.size()); j++ ) {
                RNG_UNIF01();
            }
            //rand();
        }

        // check if anything was actually updated (e.g. everything was fixed)
        all_fixed = true;
        for ( unsigned int j = 0; j < node_indices.size (); j++ )
            if ( !(positions [ node_indices[j] ].fixed && real_fixed) ) {
                all_fixed = false;
            }

        // update positions across processors (if not all fixed)
        if ( !all_fixed ) {
#ifdef MUSE_MPI
            MPI_Allgather ( &new_positions[2 * myid], 2, MPI_FLOAT,
                            new_positions, 2, MPI_FLOAT, MPI_COMM_WORLD );
#endif

            // update positions (old to new)
            update_density ( node_indices, old_positions, new_positions );
        }

        /*
        if ( myid == 0 )
          {
            // output node list (for debugging)
            for ( unsigned int j = 0; j < node_indices.size(); j++ )
              cout << node_indices[j] << " ";
            cout << endl;
          }
        */

        // compute node list for next update
        for ( unsigned int j = 0; j < node_indices.size(); j++ ) {
            node_indices [j] += num_procs;
        }

        while ( !node_indices.empty() && node_indices.back() >= num_nodes ) {
            node_indices.pop_back ( );
        }

    }

    // update first_add and fine_first_add
    first_add = false;
    if ( fineDensity ) {
        fine_first_add = false;
    }

}

// The get_positions function takes the node_indices list
// and returns the corresponding positions in an array.

void graph::get_positions ( vector<int> &node_indices,
                            float return_positions[3 * MAX_PROCS]  ) {

    // fill positions
    for (unsigned int i = 0; i < node_indices.size(); i++) {
        return_positions[3 * i] = positions[ node_indices[i] ].x;
        return_positions[3 * i + 1] = positions[ node_indices[i] ].y;
        return_positions[3 * i + 2] = positions[ node_indices[i] ].z;
    }

}

// update_node_pos -- this subroutine does the actual work of computing
// the new position of a given node.  num_act_proc gives the number
// of active processes at this level for use by the random number
// generators.

void graph::update_node_pos ( int node_ind,
                              float old_positions[3 * MAX_PROCS],
                              float new_positions[3 * MAX_PROCS] ) {

    float energies[2];          // node energies for possible positions
    float updated_pos[2][3];    // possible positions
    float pos_x, pos_y, pos_z;

    // old VxOrd parameter
    float jump_length = .010 * temperature;

    // subtract old node
    density_server.Subtract ( positions[node_ind], first_add, fine_first_add, fineDensity );

    // compute node energy for old solution
    energies[0] = Compute_Node_Energy ( node_ind );

    // move node to centroid position
    Solve_Analytic ( node_ind, pos_x, pos_y, pos_z );
    positions[node_ind].x = updated_pos[0][0] = pos_x;
    positions[node_ind].y = updated_pos[0][1] = pos_y;
    positions[node_ind].z = updated_pos[0][2] = pos_z;

    /*
    // ouput random numbers (for debugging)
    int rand_0, rand_1;
    rand_0 = rand();
    rand_1 = rand();
    cout << myid << ": " << rand_0 << ", " << rand_1 << endl;
    */

    // Do random method (RAND_MAX is C++ maximum random number)
    updated_pos[1][0] = updated_pos[0][0] + (.5 - RNG_UNIF01()) * jump_length;
    updated_pos[1][1] = updated_pos[0][1] + (.5 - RNG_UNIF01()) * jump_length;
    updated_pos[1][2] = updated_pos[0][2] + (.5 - RNG_UNIF01()) * jump_length;

    // compute node energy for random position
    positions[node_ind].x = updated_pos[1][0];
    positions[node_ind].y = updated_pos[1][1];
    positions[node_ind].z = updated_pos[1][2];
    energies[1] = Compute_Node_Energy ( node_ind );

    /*
    // output update possiblities (debugging):
    cout << node_ind << ": (" << updated_pos[0][0] << "," << updated_pos[0][1]
         << "), " << energies[0] << "; (" << updated_pos[1][0] << ","
         << updated_pos[1][1] << "), " << energies[1] << endl;
    */

    // add back old position
    positions[node_ind].x = old_positions[3 * myid];
    positions[node_ind].y = old_positions[3 * myid + 1];
    positions[node_ind].z = old_positions[3 * myid + 2];
    if ( !fineDensity && !first_add ) {
        density_server.Add ( positions[node_ind], fineDensity );
    } else if ( !fine_first_add ) {
        density_server.Add ( positions[node_ind], fineDensity );
    }

    // choose updated node position with lowest energy
    if ( energies[0] < energies[1] ) {
        new_positions[3 * myid] = updated_pos[0][0];
        new_positions[3 * myid + 1] = updated_pos[0][1];
        new_positions[3 * myid + 2] = updated_pos[0][2];
        positions[node_ind].energy = energies[0];
    } else {
        new_positions[3 * myid] = updated_pos[1][0];
        new_positions[3 * myid + 1] = updated_pos[1][1];
        new_positions[3 * myid + 2] = updated_pos[1][2];
        positions[node_ind].energy = energies[1];
    }

}

// update_density takes a sequence of node_indices and their positions and
// updates the positions by subtracting the old positions and adding the
// new positions to the density grid.

void graph::update_density ( vector<int> &node_indices,
                             float old_positions[3 * MAX_PROCS],
                             float new_positions[3 * MAX_PROCS] ) {

    // go through each node and subtract old position from
    // density grid before adding new position
    for ( unsigned int i = 0; i < node_indices.size(); i++ ) {
        positions[node_indices[i]].x = old_positions[3 * i];
        positions[node_indices[i]].y = old_positions[3 * i + 1];
        positions[node_indices[i]].z = old_positions[3 * i + 2];
        density_server.Subtract ( positions[node_indices[i]],
                                  first_add, fine_first_add, fineDensity );

        positions[node_indices[i]].x = new_positions[3 * i];
        positions[node_indices[i]].y = new_positions[3 * i + 1];
        positions[node_indices[i]].z = new_positions[3 * i + 2];
        density_server.Add ( positions[node_indices[i]], fineDensity );
    }

}

/********************************************
* Function: Compute_Node_Energy             *
* Description: Compute the node energy      *
* This code has been modified from the      *
* original code by B. Wylie.                *
*********************************************/

float graph::Compute_Node_Energy( int node_ind ) {

    /* Want to expand 4th power range of attraction */
    float attraction_factor = attraction * attraction *
                              attraction * attraction * 2e-2;

    map <int, float>::iterator EI;
    float x_dis, y_dis, z_dis;
    float energy_distance, weight;
    float node_energy = 0;

    // Add up all connection energies
    for (EI = neighbors[node_ind].begin(); EI != neighbors[node_ind].end(); ++EI) {

        // Get edge weight
        weight = EI->second;

        // Compute x,y distance
        x_dis = positions[ node_ind ].x - positions[ EI->first ].x;
        y_dis = positions[ node_ind ].y - positions[ EI->first ].y;
        z_dis = positions[ node_ind ].z - positions[ EI->first ].z;

        // Energy Distance
        energy_distance = x_dis * x_dis + y_dis * y_dis + z_dis * z_dis;
        if (STAGE < 2) {
            energy_distance *= energy_distance;
        }

        // In the liquid phase we want to discourage long link distances
        if (STAGE == 0) {
            energy_distance *= energy_distance;
        }

        node_energy += weight * attraction_factor * energy_distance;
    }

    // output effect of density (debugging)
    //cout << "[before: " << node_energy;

    // add density
    node_energy += density_server.GetDensity ( positions[ node_ind ].x, positions[ node_ind ].y,
                   positions[ node_ind ].z, fineDensity );

    // after calling density server (debugging)
    //cout << ", after: " << node_energy << "]" << endl;

    // return computated energy
    return node_energy;
}


/*********************************************
* Function: Solve_Analytic                   *
* Description: Compute the node position     *
* This is a modified version of the function *
* originally written by B. Wylie             *
*********************************************/

void graph::Solve_Analytic( int node_ind, float &pos_x, float &pos_y,
                            float &pos_z) {

    map <int, float>::iterator EI;
    float total_weight = 0;
    float x_dis, y_dis, z_dis, x_cen = 0, y_cen = 0, z_cen = 0;
    float x = 0, y = 0, z = 0, dis;
    float damping, weight;

    // Sum up all connections
    for (EI = neighbors[node_ind].begin(); EI != neighbors[node_ind].end(); ++EI) {
        weight = EI->second;
        total_weight += weight;
        x +=  weight * positions[ EI->first ].x;
        y +=  weight * positions[ EI->first ].y;
        z +=  weight * positions[ EI->first ].z;
    }

    // Now set node position
    if (total_weight > 0) {

        // Compute centriod
        x_cen = x / total_weight;
        y_cen = y / total_weight;
        z_cen = z / total_weight;
        damping = 1.0 - damping_mult;
        pos_x = damping * positions[ node_ind ].x + (1.0 - damping) * x_cen;
        pos_y = damping * positions[ node_ind ].y + (1.0 - damping) * y_cen;
        pos_z = damping * positions[ node_ind ].z + (1.0 - damping) * z_cen;
    }

    // No cut edge flag (?)
    if (min_edges == 99) {
        return;
    }

    // Don't cut at end of scale
    if ( CUT_END >= 39500 ) {
        return;
    }

    float num_connections = (float)sqrt((float)neighbors[node_ind].size());
    float maxLength = 0;

    map<int, float>::iterator maxIndex;

    // Go through nodes edges... cutting if necessary
    for (EI = maxIndex = neighbors[node_ind].begin();
         EI != neighbors[node_ind].end(); ++EI) {

        // Check for at least min edges
        if (neighbors[node_ind].size() < min_edges) {
            continue;
        }

        x_dis = x_cen - positions[ EI->first ].x;
        y_dis = y_cen - positions[ EI->first ].y;
        z_dis = z_cen - positions[ EI->first ].z;
        dis = x_dis * x_dis + y_dis * y_dis + z_dis * z_dis;
        dis *= num_connections;

        // Store maximum edge
        if (dis > maxLength) {
            maxLength = dis;
            maxIndex = EI;
        }
    }

    // If max length greater than cut_length then cut
    if (maxLength > cut_off_length) {
        neighbors[ node_ind ].erase( maxIndex );
    }

}


// get_tot_energy adds up the energy for each node to give an estimate of the
// quality of the minimization.

float graph::get_tot_energy ( ) {

    float my_tot_energy, tot_energy;
    my_tot_energy = 0;
    for ( int i = myid; i < num_nodes; i += num_procs ) {
        my_tot_energy += positions[i].energy;
    }

    //vector<Node>::iterator i;
    //for ( i = positions.begin(); i != positions.end(); i++ )
    //  tot_energy += i->energy;

#ifdef MUSE_MPI
    MPI_Reduce ( &my_tot_energy, &tot_energy, 1, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD );
#else
    tot_energy = my_tot_energy;
#endif

    return tot_energy;

}


int graph::draw_graph(igraph_matrix_t *res) {
    int count_iter = 0;
    while (ReCompute()) {
        IGRAPH_ALLOW_INTERRUPTION();
        count_iter++;
    }
    long int n = positions.size();
    IGRAPH_CHECK(igraph_matrix_resize(res, n, 3));
    for (long int i = 0; i < n; i++) {
        MATRIX(*res, i, 0) = positions[i].x;
        MATRIX(*res, i, 1) = positions[i].y;
        MATRIX(*res, i, 2) = positions[i].z;
    }
    return 0;
}

} // namespace drl3d
