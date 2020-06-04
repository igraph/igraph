/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2007-2012  Gabor Csardi <csardi.gabor@gmail.com>
   334 Harvard street, Cambridge, MA 02139 USA

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

#include "prpack.h"
#include "prpack/prpack_igraph_graph.h"
#include "prpack/prpack_solver.h"
#include "igraph_error.h"

using namespace prpack;
using namespace std;

/*
 * PRPACK-based implementation of \c igraph_personalized_pagerank.
 *
 * See \c igraph_personalized_pagerank for the documentation of the parameters.
 */
int igraph_personalized_pagerank_prpack(const igraph_t *graph, igraph_vector_t *vector,
                                        igraph_real_t *value, const igraph_vs_t vids,
                                        igraph_bool_t directed, igraph_real_t damping,
                                        igraph_vector_t *reset,
                                        const igraph_vector_t *weights) {
    long int i, no_of_nodes = igraph_vcount(graph), nodes_to_calc;
    igraph_vit_t vit;
    double* u = 0;
    double* v = 0;
    const prpack_result* res;

    if (reset) {
        /* Normalize reset vector so the sum is 1 */
        double reset_sum = igraph_vector_sum(reset);
        if (igraph_vector_min(reset) < 0) {
            IGRAPH_ERROR("the reset vector must not contain negative elements", IGRAPH_EINVAL);
        }
        if (reset_sum == 0) {
            IGRAPH_ERROR("the sum of the elements in the reset vector must not be zero", IGRAPH_EINVAL);
        }

        // Construct the personalization vector
        v = new double[no_of_nodes];
        for (i = 0; i < no_of_nodes; i++) {
            v[i] = VECTOR(*reset)[i] / reset_sum;
        }
    }

    // Construct and run the solver
    prpack_igraph_graph prpack_graph(graph, weights, directed);
    prpack_solver solver(&prpack_graph, false);
    res = solver.solve(damping, 1e-10, u, v, "");

    // Delete the personalization vector
    if (v) {
        delete[] v;
    }

    // Check whether the solver converged
    // TODO: this is commented out because some of the solvers do not implement it yet
    /*
    if (!res->converged) {
        IGRAPH_WARNING("PRPACK solver failed to converge. Results may be inaccurate.");
    }
    */

    // Fill the result vector
    IGRAPH_CHECK(igraph_vit_create(graph, vids, &vit));
    IGRAPH_FINALLY(igraph_vit_destroy, &vit);
    nodes_to_calc = IGRAPH_VIT_SIZE(vit);
    IGRAPH_CHECK(igraph_vector_resize(vector, nodes_to_calc));
    for (IGRAPH_VIT_RESET(vit), i = 0; !IGRAPH_VIT_END(vit);
         IGRAPH_VIT_NEXT(vit), i++) {
        VECTOR(*vector)[i] = res->x[(long int)IGRAPH_VIT_GET(vit)];
    }
    igraph_vit_destroy(&vit);
    IGRAPH_FINALLY_CLEAN(1);

    // TODO: can we get the eigenvalue? We'll just fake it until we can.
    if (value) {
        *value = 1.0;
    }
    delete res;

    return IGRAPH_SUCCESS;
}

