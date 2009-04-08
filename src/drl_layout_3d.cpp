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
#include <iostream>
#include <fstream>
#include <map>
#include <set>
#include <string>
#include <deque>
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
#include "igraph.h"
#include "random.h"

/**
 * \function igraph_layout_drl_3d
 * The DrL layout generator, 3d version.
 * 
 * This function implements the force-directed DrL layout generator.
 * Please see more at http://www.cs.sandia.gov/~smartin/software.html
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
  
  RNG_BEGIN();

  drl3d::graph neighbors(graph, options, weights);
  neighbors.init_parms(options);
  if (use_seed) {
    IGRAPH_CHECK(igraph_matrix_resize(res, igraph_vcount(graph), 3));
    neighbors.read_real(res, fixed);
  }
  neighbors.draw_graph(res);

  RNG_END();
  
  return 0;
}	      
