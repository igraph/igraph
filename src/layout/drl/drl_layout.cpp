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
// Layout
//
// This program implements a parallel force directed graph drawing
// algorithm.  The algorithm used is based upon a random decomposition
// of the graph and simulated shared memory of node position and density.
// In this version, the simulated shared memory is spread among all processors
//
// The structure of the inputs and outputs of this code will be displayed
// if the program is called without parameters, or if an erroneous
// parameter is passed to the program.
//
// S. Martin
// 5/6/2005

// C++ library routines
#include <map>
#include <vector>

using namespace std;

// layout routines and constants
#include "drl_layout.h"
#include "drl_parse.h"
#include "drl_graph.h"

// MPI
#ifdef MUSE_MPI
    #include <mpi.h>
#endif

using namespace drl;
#include "igraph_layout.h"
#include "igraph_random.h"
#include "igraph_interface.h"

#include "core/exceptions.h"

namespace drl {

// int main(int argc, char **argv) {


//   // initialize MPI
//   int myid, num_procs;

//   #ifdef MUSE_MPI
//     MPI_Init ( &argc, &argv );
//     MPI_Comm_size ( MPI_COMM_WORLD, &num_procs );
//     MPI_Comm_rank ( MPI_COMM_WORLD, &myid );
//   #else
//     myid = 0;
//  num_procs = 1;
//   #endif

//   // parameters that must be broadcast to all processors
//   int rand_seed;
//   float edge_cut;

//   char int_file[MAX_FILE_NAME];
//   char coord_file[MAX_FILE_NAME];
//   char real_file[MAX_FILE_NAME];
//   char parms_file[MAX_FILE_NAME];

//   int int_out = 0;
//   int edges_out = 0;
//   int parms_in = 0;
//   float real_in = -1.0;

//   // user interaction is handled by processor 0
//   if ( myid == 0 )
//   {
//     if ( num_procs > MAX_PROCS )
//  {
//      cout << "Error: Maximum number of processors is " << MAX_PROCS << "." << endl;
//      cout << "Adjust compile time parameter." << endl;
//      #ifdef MUSE_MPI
//        MPI_Abort ( MPI_COMM_WORLD, 1 );
//      #else
//        exit (1);
//      #endif
//  }

//  // get user input
//     parse command_line ( argc, argv );
//  rand_seed = command_line.rand_seed;
//  edge_cut = command_line.edge_cut;
//  int_out = command_line.int_out;
//  edges_out = command_line.edges_out;
//  parms_in = command_line.parms_in;
//  real_in = command_line.real_in;
//  strcpy ( coord_file, command_line.coord_file.c_str() );
//  strcpy ( int_file, command_line.sim_file.c_str() );
//  strcpy ( real_file, command_line.real_file.c_str() );
//  strcpy ( parms_file, command_line.parms_file.c_str() );

//   }

//   // now we initialize all processors by reading .int file
//   #ifdef MUSE_MPI
//     MPI_Bcast ( &int_file, MAX_FILE_NAME, MPI_CHAR, 0, MPI_COMM_WORLD );
//   #endif
//   graph neighbors ( myid, num_procs, int_file );

//   // check for user supplied parameters
//   #ifdef MUSE_MPI
//     MPI_Bcast ( &parms_in, 1, MPI_INT, 0, MPI_COMM_WORLD );
//   #endif
//   if ( parms_in )
//   {
//     #ifdef MUSE_MPI
//    MPI_Bcast ( &parms_file, MAX_FILE_NAME, MPI_CHAR, 0, MPI_COMM_WORLD );
//  #endif
//  neighbors.read_parms ( parms_file );
//   }

//   // set random seed, edge cutting, and real iterations parameters
//   #ifdef MUSE_MPI
//     MPI_Bcast ( &rand_seed, 1, MPI_INT, 0, MPI_COMM_WORLD );
//     MPI_Bcast ( &edge_cut, 1, MPI_FLOAT, 0, MPI_COMM_WORLD );
//  MPI_Bcast ( &real_in, 1, MPI_INT, 0, MPI_COMM_WORLD );
//   #endif
//   neighbors.init_parms ( rand_seed, edge_cut, real_in );

//   // check for .real file with existing coordinates
//   if ( real_in >= 0 )
//   {
//     #ifdef MUSE_MPI
//    MPI_Bcast ( &real_file, MAX_FILE_NAME, MPI_CHAR, 0, MPI_COMM_WORLD );
//  #endif
//  neighbors.read_real ( real_file );
//   }

//   neighbors.draw_graph ( int_out, coord_file );

//   // do we have to write out the edges?
//   #ifdef MUSE_MPI
//     MPI_Bcast ( &edges_out, 1, MPI_INT, 0, MPI_COMM_WORLD );
//   #endif
//   if ( edges_out )
//     {
//    #ifdef MUSE_MPI
//         MPI_Bcast ( &coord_file, MAX_FILE_NAME, MPI_CHAR, 0, MPI_COMM_WORLD );
//    #endif
//       for ( int i = 0; i < num_procs; i++ )
//    {
//      if ( myid == i )
//        neighbors.write_sim ( coord_file );
//      #ifdef MUSE_MPI
//            MPI_Barrier ( MPI_COMM_WORLD );
//      #endif
//    }
//     }

//   // finally we output file and quit
//   float tot_energy;
//   tot_energy = neighbors.get_tot_energy ();
//   if ( myid == 0 )
//   {
//  neighbors.write_coord ( coord_file );
//  cout << "Total Energy: " << tot_energy << "." << endl
//       << "Program terminated successfully." << endl;
//   }

//   // MPI finalize
//   #ifdef MUSE_MPI
//     MPI_Finalize ();
//   #endif

//   return 0;
// }

} // namespace drl

/**
 * \section about_drl
 *
 * <para>
 * DrL is a sophisticated layout generator developed and implemented by
 * Shawn Martin et al. As of October 2012 the original DrL homepage is
 * unfortunately not available. You can read more about this algorithm
 * in the following technical report: Martin, S., Brown, W.M.,
 * Klavans, R., Boyack, K.W., DrL: Distributed Recursive (Graph)
 * Layout. SAND Reports, 2008. 2936: p. 1-10.
 * </para>
 *
 * <para>
 * Only a subset of the complete DrL functionality is
 * included in igraph, parallel runs and recursive, multi-level
 * layouting is not supported.
 * </para>
 *
 * <para>
 * The parameters of the layout are stored in an \ref
 * igraph_layout_drl_options_t structure, this can be initialized by
 * calling the function \ref igraph_layout_drl_options_init().
 * The fields of this structure can then be adjusted by hand if needed.
 * The layout is calculated by an \ref igraph_layout_drl() call.
 * </para>
 */

/**
 * \function igraph_layout_drl_options_init
 * Initialize parameters for the DrL layout generator
 *
 * This function can be used to initialize the struct holding the
 * parameters for the DrL layout generator. There are a number of
 * predefined templates available, it is a good idea to start from one
 * of these by modifying some parameters.
 * \param options The struct to initialize.
 * \param templ The template to use. Currently the following templates
 *     are supplied: \c IGRAPH_LAYOUT_DRL_DEFAULT, \c
 *     IGRAPH_LAYOUT_DRL_COARSEN, \c IGRAPH_LAYOUT_DRL_COARSEST,
 *     \c IGRAPH_LAYOUT_DRL_REFINE and \c IGRAPH_LAYOUT_DRL_FINAL.
 * \return Error code.
 *
 * Time complexity: O(1).
 */

int igraph_layout_drl_options_init(igraph_layout_drl_options_t *options,
                                   igraph_layout_drl_default_t templ) {

    options->edge_cut = 32.0 / 40.0;

    switch (templ) {
    case IGRAPH_LAYOUT_DRL_DEFAULT:
        options->init_iterations   = 0;
        options->init_temperature  = 2000;
        options->init_attraction   = 10;
        options->init_damping_mult = 1.0;

        options->liquid_iterations   = 200;
        options->liquid_temperature  = 2000;
        options->liquid_attraction   = 10;
        options->liquid_damping_mult = 1.0;

        options->expansion_iterations   = 200;
        options->expansion_temperature  = 2000;
        options->expansion_attraction   = 2;
        options->expansion_damping_mult = 1.0;

        options->cooldown_iterations   = 200;
        options->cooldown_temperature  = 2000;
        options->cooldown_attraction   = 1;
        options->cooldown_damping_mult = .1;

        options->crunch_iterations   = 50;
        options->crunch_temperature  = 250;
        options->crunch_attraction   = 1;
        options->crunch_damping_mult = 0.25;

        options->simmer_iterations   = 100;
        options->simmer_temperature  = 250;
        options->simmer_attraction   = .5;
        options->simmer_damping_mult = 0;

        break;
    case IGRAPH_LAYOUT_DRL_COARSEN:
        options->init_iterations   = 0;
        options->init_temperature  = 2000;
        options->init_attraction   = 10;
        options->init_damping_mult = 1.0;

        options->liquid_iterations   = 200;
        options->liquid_temperature  = 2000;
        options->liquid_attraction   = 2;
        options->liquid_damping_mult = 1.0;

        options->expansion_iterations   = 200;
        options->expansion_temperature  = 2000;
        options->expansion_attraction   = 10;
        options->expansion_damping_mult = 1.0;

        options->cooldown_iterations   = 200;
        options->cooldown_temperature  = 2000;
        options->cooldown_attraction   = 1;
        options->cooldown_damping_mult = .1;

        options->crunch_iterations   = 50;
        options->crunch_temperature  = 250;
        options->crunch_attraction   = 1;
        options->crunch_damping_mult = 0.25;

        options->simmer_iterations   = 100;
        options->simmer_temperature  = 250;
        options->simmer_attraction   = .5;
        options->simmer_damping_mult = 0;

        break;
    case IGRAPH_LAYOUT_DRL_COARSEST:
        options->init_iterations   = 0;
        options->init_temperature  = 2000;
        options->init_attraction   = 10;
        options->init_damping_mult = 1.0;

        options->liquid_iterations   = 200;
        options->liquid_temperature  = 2000;
        options->liquid_attraction   = 2;
        options->liquid_damping_mult = 1.0;

        options->expansion_iterations   = 200;
        options->expansion_temperature  = 2000;
        options->expansion_attraction   = 10;
        options->expansion_damping_mult = 1.0;

        options->cooldown_iterations   = 200;
        options->cooldown_temperature  = 2000;
        options->cooldown_attraction   = 1;
        options->cooldown_damping_mult = .1;

        options->crunch_iterations   = 200;
        options->crunch_temperature  = 250;
        options->crunch_attraction   = 1;
        options->crunch_damping_mult = 0.25;

        options->simmer_iterations   = 100;
        options->simmer_temperature  = 250;
        options->simmer_attraction   = .5;
        options->simmer_damping_mult = 0;

        break;
    case IGRAPH_LAYOUT_DRL_REFINE:
        options->init_iterations   = 0;
        options->init_temperature  = 50;
        options->init_attraction   = .5;
        options->init_damping_mult = 0;

        options->liquid_iterations   = 0;
        options->liquid_temperature  = 2000;
        options->liquid_attraction   = 2;
        options->liquid_damping_mult = 1.0;

        options->expansion_iterations   = 50;
        options->expansion_temperature  = 500;
        options->expansion_attraction   = .1;
        options->expansion_damping_mult = .25;

        options->cooldown_iterations   = 50;
        options->cooldown_temperature  = 200;
        options->cooldown_attraction   = 1;
        options->cooldown_damping_mult = .1;

        options->crunch_iterations   = 50;
        options->crunch_temperature  = 250;
        options->crunch_attraction   = 1;
        options->crunch_damping_mult = 0.25;

        options->simmer_iterations   = 0;
        options->simmer_temperature  = 250;
        options->simmer_attraction   = .5;
        options->simmer_damping_mult = 0;

        break;
    case IGRAPH_LAYOUT_DRL_FINAL:
        options->init_iterations   = 0;
        options->init_temperature  = 50;
        options->init_attraction   = .5;
        options->init_damping_mult = 0;

        options->liquid_iterations   = 0;
        options->liquid_temperature  = 2000;
        options->liquid_attraction   = 2;
        options->liquid_damping_mult = 1.0;

        options->expansion_iterations   = 50;
        options->expansion_temperature  = 50;
        options->expansion_attraction   = .1;
        options->expansion_damping_mult = .25;

        options->cooldown_iterations   = 50;
        options->cooldown_temperature  = 200;
        options->cooldown_attraction   = 1;
        options->cooldown_damping_mult = .1;

        options->crunch_iterations   = 50;
        options->crunch_temperature  = 250;
        options->crunch_attraction   = 1;
        options->crunch_damping_mult = 0.25;

        options->simmer_iterations   = 25;
        options->simmer_temperature  = 250;
        options->simmer_attraction   = .5;
        options->simmer_damping_mult = 0;

        break;
    default:
        IGRAPH_ERROR("Unknown DrL template", IGRAPH_EINVAL);
        break;
    }

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_layout_drl
 * The DrL layout generator
 *
 * This function implements the force-directed DrL layout generator.
 * Please see more in the following technical report: Martin, S.,
 * Brown, W.M., Klavans, R., Boyack, K.W., DrL: Distributed Recursive
 * (Graph) Layout. SAND Reports, 2008. 2936: p. 1-10.
 * \param graph The input graph.
 * \param use_seed Logical scalar, if true, then the coordinates
 *    supplied in the \p res argument are used as starting points.
 * \param res Pointer to a matrix, the result layout is stored
 *    here. It will be resized as needed.
 * \param options The parameters to pass to the layout generator.
 * \param weights Edge weights, pointer to a vector. If this is a null
 *    pointer then every edge will have the same weight.
 * \param fixed Pointer to a logical vector, or a null pointer. Originally,
 *    this argument was used in the DrL algorithm to keep the nodes marked
 *    with this argument as fixed; fixed nodes would then keep their
 *    positions in the initial stages of the algorithm. However, due to how
 *    the DrL code imported into igraph is organized, it seems that the
 *    argument does not do anything and we are not sure whether this is a
 *    bug or a feature in DrL. We are leaving the argument here in order not
 *    to break the API, but note that at the present stage it has no effect.
 * \return Error code.
 *
 * Time complexity: ???.
 */

int igraph_layout_drl(const igraph_t *graph, igraph_matrix_t *res,
                      igraph_bool_t use_seed,
                      igraph_layout_drl_options_t *options,
                      const igraph_vector_t *weights,
                      const igraph_vector_bool_t *fixed) {
    const char msg[] = "Damping multipliers cannot be negative, got %f.";

    if (options->init_damping_mult < 0) {
        IGRAPH_ERRORF(msg, IGRAPH_EINVAL, options->init_damping_mult);
    }
    if (options->liquid_damping_mult < 0) {
        IGRAPH_ERRORF(msg, IGRAPH_EINVAL, options->liquid_damping_mult);
    }
    if (options->expansion_damping_mult < 0) {
        IGRAPH_ERRORF(msg, IGRAPH_EINVAL, options->expansion_damping_mult);
    }
    if (options->cooldown_damping_mult < 0) {
        IGRAPH_ERRORF(msg, IGRAPH_EINVAL, options->cooldown_damping_mult);
    }
    if (options->crunch_damping_mult < 0) {
        IGRAPH_ERRORF(msg, IGRAPH_EINVAL, options->crunch_damping_mult);
    }
    if (options->simmer_damping_mult < 0) {
        IGRAPH_ERRORF(msg, IGRAPH_EINVAL, options->simmer_damping_mult);
    }

    IGRAPH_HANDLE_EXCEPTIONS(
        RNG_BEGIN();

        drl::graph neighbors(graph, options, weights);
        neighbors.init_parms(options);
        if (use_seed) {
            IGRAPH_CHECK(igraph_matrix_resize(res, igraph_vcount(graph), 2));
            neighbors.read_real(res, fixed);
        }
        neighbors.draw_graph(res);

        RNG_END();
    );

    return IGRAPH_SUCCESS;
}
