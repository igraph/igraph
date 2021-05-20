/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2007-2021  The igraph development team <igraph@igraph.org>

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

#include "igraph_error.h"

#include "centrality/prpack_internal.h"
#include "centrality/prpack/prpack_igraph_graph.h"
#include "centrality/prpack/prpack_solver.h"
#include "core/exceptions.h"

using namespace prpack;
using namespace std;

/*
 * PRPACK-based implementation of \c igraph_personalized_pagerank.
 *
 * See \c igraph_personalized_pagerank for the documentation of the parameters.
 */
int igraph_i_personalized_pagerank_prpack(const igraph_t *graph, igraph_vector_t *vector,
                                          igraph_real_t *value, const igraph_vs_t vids,
                                          igraph_bool_t directed, igraph_real_t damping,
                                          const igraph_vector_t *reset,
                                          const igraph_vector_t *weights) {
    long int i, no_of_nodes = igraph_vcount(graph), nodes_to_calc;
    igraph_vit_t vit;
    double *u = nullptr;
    double *v = nullptr;
    const prpack_result *res;

    IGRAPH_HANDLE_EXCEPTIONS(
        if (reset) {
            if (igraph_vector_size(reset) != no_of_nodes) {
                IGRAPH_ERROR("Invalid length of reset vector when calculating personalized PageRank scores.", IGRAPH_EINVAL);
            }

            /* Normalize reset vector so the sum is 1 */
            double reset_min = igraph_vector_min(reset);
            if (reset_min < 0) {
                IGRAPH_ERROR("The reset vector must not contain negative elements.", IGRAPH_EINVAL);
            }
            if (igraph_is_nan(reset_min)) {
                IGRAPH_ERROR("The reset vector must not contain NaN values.", IGRAPH_EINVAL);
            }

            double reset_sum = igraph_vector_sum(reset);
            if (reset_sum == 0) {
                IGRAPH_ERROR("The sum of the elements in the reset vector must not be zero.", IGRAPH_EINVAL);
            }

            // Construct the personalization vector
            v = new double[no_of_nodes];
            for (i = 0; i < no_of_nodes; i++) {
                v[i] = VECTOR(*reset)[i] / reset_sum;
            }

            // u is the distribution used when restarting the walk due to being stuck in a sink
            // v is the distribution used when restarting due to damping
            // Here we use the same distribution for both
            u = v;
        }

        // Since PRPACK uses the algebraic method to solve PageRank, damping factors very close to 1.0
        // may lead to numerical instability, the apperance of non-finite values, or the iteration
        // never terminating.
        if (damping > 0.999) {
            IGRAPH_WARNINGF(
                    "Damping factor is %g. "
                    "Damping values close to 1 may lead to numerical instability when using PRPACK.",
                    damping);
        }

        // Construct and run the solver
        prpack_igraph_graph prpack_graph(graph, weights, directed);
        prpack_solver solver(&prpack_graph, false);
        res = solver.solve(damping, 1e-10, u, v, "");

        // Delete the personalization vector
        delete [] v;
    );

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
