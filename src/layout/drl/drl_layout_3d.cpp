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
#include "drl_layout_3d.h"
#include "drl_parse.h"
#include "drl_graph_3d.h"

// MPI
#ifdef MUSE_MPI
    #include <mpi.h>
#endif

using namespace drl3d;
#include "igraph_layout.h"
#include "igraph_random.h"
#include "igraph_interface.h"

#include "core/exceptions.h"

/**
 * \function igraph_layout_drl_3d
 * The DrL layout generator, 3d version.
 *
 * This function implements the force-directed DrL layout generator.
 * Please see more in the technical report: Martin, S., Brown, W.M.,
 * Klavans, R., Boyack, K.W., DrL: Distributed Recursive (Graph)
 * Layout. SAND Reports, 2008. 2936: p. 1-10.
 *
 * </para><para> This function uses a modified DrL generator that does
 * the layout in three dimensions.
 * \param graph The input graph.
 * \param use_seed Logical scalar, if true, then the coordinates
 *    supplied in the \p res argument are used as starting points.
 * \param res Pointer to a matrix, the result layout is stored
 *    here. It will be resized as needed.
 * \param options The parameters to pass to the layout generator.
 * \param weights Edge weights, pointer to a vector. If this is a null
 *    pointer then every edge will have the same weight.
 * \param fixed Pointer to a logical vector, or a null pointer. This
 *    can be used to fix the position of some vertices. Vertices for
 *    which it is true will not be moved, but stay at the coordinates
 *    given in the \p res matrix. This argument is ignored if it is a
 *    null pointer or if use_seed is false.
 * \return Error code.
 *
 * Time complexity: ???.
 *
 * \sa \ref igraph_layout_drl() for the standard 2d version.
 */

int igraph_layout_drl_3d(const igraph_t *graph, igraph_matrix_t *res,
                         igraph_bool_t use_seed,
                         igraph_layout_drl_options_t *options,
                         const igraph_vector_t *weights,
                         const igraph_vector_bool_t *fixed) {
    IGRAPH_HANDLE_EXCEPTIONS(
        RNG_BEGIN();

        drl3d::graph neighbors(graph, options, weights);
        neighbors.init_parms(options);
        if (use_seed) {
            IGRAPH_CHECK(igraph_matrix_resize(res, igraph_vcount(graph), 3));
            neighbors.read_real(res, fixed);
        }
        neighbors.draw_graph(res);

        RNG_END();
    );

    return 0;
}
