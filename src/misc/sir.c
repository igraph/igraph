/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2014  Gabor Csardi <csardi.gabor@gmail.com>
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

#include "igraph_epidemics.h"

#include "igraph_random.h"
#include "igraph_adjlist.h"
#include "igraph_interface.h"
#include "igraph_psumtree.h"
#include "igraph_memory.h"
#include "igraph_structural.h"

#include "core/interruption.h"

int igraph_sir_init(igraph_sir_t *sir) {
    IGRAPH_CHECK(igraph_vector_init(&sir->times, 1));
    IGRAPH_FINALLY(igraph_vector_destroy, &sir->times);
    IGRAPH_CHECK(igraph_vector_int_init(&sir->no_s, 1));
    IGRAPH_FINALLY(igraph_vector_int_destroy, &sir->no_s);
    IGRAPH_CHECK(igraph_vector_int_init(&sir->no_i, 1));
    IGRAPH_FINALLY(igraph_vector_int_destroy, &sir->no_i);
    IGRAPH_CHECK(igraph_vector_int_init(&sir->no_r, 1));
    IGRAPH_FINALLY_CLEAN(3);
    return 0;
}

/**
 * \function igraph_sir_destroy
 * \brief Deallocates memory associated with a SIR simulation run.
 *
 * \param sir The \ref igraph_sir_t object storing the simulation.
 */

void igraph_sir_destroy(igraph_sir_t *sir) {
    igraph_vector_destroy(&sir->times);
    igraph_vector_int_destroy(&sir->no_s);
    igraph_vector_int_destroy(&sir->no_i);
    igraph_vector_int_destroy(&sir->no_r);
}

static void igraph_i_sir_destroy(igraph_vector_ptr_t *v) {
    int i, n = igraph_vector_ptr_size(v);
    for (i = 0; i < n; i++) {
        if ( VECTOR(*v)[i] ) {
            igraph_sir_destroy( VECTOR(*v)[i]) ;
            IGRAPH_FREE( VECTOR(*v)[i] ); /* this also sets the vector_ptr element to NULL */
        }
    }
}

#define S_S 0
#define S_I 1
#define S_R 2

/**
 * \function igraph_sir
 * \brief Performs a number of SIR epidemics model runs on a graph.
 *
 * The SIR model is a simple model from epidemiology. The individuals
 * of the population might be in three states: susceptible, infected
 * and recovered. Recovered people are assumed to be immune to the
 * disease. Susceptibles become infected with a rate that depends on
 * their number of infected neigbors. Infected people become recovered
 * with a constant rate. See these parameters below.
 *
 * </para><para>
 * This function runs multiple simulations, all starting with a
 * single uniformly randomly chosen infected individual. A simulation
 * is stopped when no infected individuals are left.
 *
 * \param graph The graph to perform the model on. For directed graphs
 *        edge directions are ignored and a warning is given.
 * \param beta The rate of infection of an individual that is
 *        susceptible and has a single infected neighbor.
 *        The infection rate of a susceptible individual with n
 *        infected neighbors is n times beta. Formally
 *        this is the rate parameter of an exponential distribution.
 * \param gamma The rate of recovery of an infected individual.
 *        Formally, this is the rate parameter of an exponential
 *        distribution.
 * \param no_sim The number of simulation runs to perform.
 * \param result The result of the simulation is stored here,
 *        in a list of \ref igraph_sir_t objects. To deallocate
 *        memory, the user needs to call \ref igraph_sir_destroy on
 *        each element, before destroying the pointer vector itself
 *        using \ref igraph_vector_ptr_destroy_all().
 * \return Error code.
 *
 * Time complexity: O(no_sim * (|V| + |E| log(|V|))).
 */

int igraph_sir(const igraph_t *graph, igraph_real_t beta,
               igraph_real_t gamma, igraph_integer_t no_sim,
               igraph_vector_ptr_t *result) {

    int infected;
    igraph_vector_int_t status;
    igraph_adjlist_t adjlist;
    int no_of_nodes = igraph_vcount(graph);
    int i, j, ns, ni, nr;
    igraph_vector_int_t *neis;
    igraph_psumtree_t tree;
    igraph_real_t psum;
    int neilen;
    igraph_bool_t simple;

    if (no_of_nodes == 0) {
        IGRAPH_ERROR("Cannot run SIR model on empty graph.", IGRAPH_EINVAL);
    }
    if (igraph_is_directed(graph)) {
        IGRAPH_WARNING("Edge directions are ignored in SIR model.");
    }
    if (beta < 0) {
        IGRAPH_ERROR("The infection rate beta must be non-negative in SIR model.", IGRAPH_EINVAL);
    }
    /* With a recovery rate of zero, the simulation would never stop. */
    if (gamma <= 0) {
        IGRAPH_ERROR("The recovery rate gamma must be positive in SIR model.", IGRAPH_EINVAL);
    }
    if (no_sim <= 0) {
        IGRAPH_ERROR("Number of SIR simulations must be positive.", IGRAPH_EINVAL);
    }

    IGRAPH_CHECK(igraph_is_simple(graph, &simple));
    if (!simple) {
        IGRAPH_ERROR("SIR model only works with simple graphs.", IGRAPH_EINVAL);
    }

    IGRAPH_CHECK(igraph_vector_int_init(&status, no_of_nodes));
    IGRAPH_FINALLY(igraph_vector_int_destroy, &status);
    IGRAPH_CHECK(igraph_adjlist_init(graph, &adjlist, IGRAPH_ALL, IGRAPH_LOOPS_TWICE, IGRAPH_MULTIPLE));
    IGRAPH_FINALLY(igraph_adjlist_destroy, &adjlist);
    IGRAPH_CHECK(igraph_psumtree_init(&tree, no_of_nodes));
    IGRAPH_FINALLY(igraph_psumtree_destroy, &tree);

    IGRAPH_CHECK(igraph_vector_ptr_resize(result, no_sim));
    igraph_vector_ptr_null(result);
    IGRAPH_FINALLY(igraph_i_sir_destroy, result);
    for (i = 0; i < no_sim; i++) {
        igraph_sir_t *sir = IGRAPH_CALLOC(1, igraph_sir_t);
        if (!sir) {
            IGRAPH_ERROR("Cannot run SIR model.", IGRAPH_ENOMEM);
        }
        IGRAPH_CHECK(igraph_sir_init(sir));
        VECTOR(*result)[i] = sir;
    }

    RNG_BEGIN();

    for (j = 0; j < no_sim; j++) {

        igraph_sir_t *sir = VECTOR(*result)[j];
        igraph_vector_t *times_v = &sir->times;
        igraph_vector_int_t *no_s_v = &sir->no_s;
        igraph_vector_int_t *no_i_v = &sir->no_i;
        igraph_vector_int_t *no_r_v = &sir->no_r;

        infected = RNG_INTEGER(0, no_of_nodes - 1);

        /* Initially infected */
        igraph_vector_int_null(&status);
        VECTOR(status)[infected] = S_I;
        ns = no_of_nodes - 1;
        ni = 1;
        nr = 0;

        VECTOR(*times_v)[0] = 0.0;
        VECTOR(*no_s_v)[0]  = ns;
        VECTOR(*no_i_v)[0]  = ni;
        VECTOR(*no_r_v)[0]  = nr;

        if (igraph_psumtree_sum(&tree) != 0) {
            igraph_psumtree_reset(&tree);
        }

        /* Rates */
        IGRAPH_CHECK(igraph_psumtree_update(&tree, infected, gamma));
        neis = igraph_adjlist_get(&adjlist, infected);
        neilen = igraph_vector_int_size(neis);
        for (i = 0; i < neilen; i++) {
            int nei = VECTOR(*neis)[i];
            IGRAPH_CHECK(igraph_psumtree_update(&tree, nei, beta));
        }

        while (ni > 0) {
            igraph_real_t tt;
            igraph_real_t r;
            long int vchange;

            IGRAPH_ALLOW_INTERRUPTION();

            psum = igraph_psumtree_sum(&tree);
            tt = igraph_rng_get_exp(igraph_rng_default(), psum);
            r = RNG_UNIF(0, psum);

            igraph_psumtree_search(&tree, &vchange, r);
            neis = igraph_adjlist_get(&adjlist, vchange);
            neilen = igraph_vector_int_size(neis);

            if (VECTOR(status)[vchange] == S_I) {
                VECTOR(status)[vchange] = S_R;
                ni--; nr++;
                IGRAPH_CHECK(igraph_psumtree_update(&tree, vchange, 0.0));
                for (i = 0; i < neilen; i++) {
                    int nei = VECTOR(*neis)[i];
                    if (VECTOR(status)[nei] == S_S) {
                        igraph_real_t rate = igraph_psumtree_get(&tree, nei);
                        IGRAPH_CHECK(igraph_psumtree_update(&tree, nei, rate - beta));
                    }
                }

            } else { /* S_S */
                VECTOR(status)[vchange] = S_I;
                ns--; ni++;
                IGRAPH_CHECK(igraph_psumtree_update(&tree, vchange, gamma));
                for (i = 0; i < neilen; i++) {
                    int nei = VECTOR(*neis)[i];
                    if (VECTOR(status)[nei] == S_S) {
                        igraph_real_t rate = igraph_psumtree_get(&tree, nei);
                        IGRAPH_CHECK(igraph_psumtree_update(&tree, nei, rate + beta));
                    }
                }
            }

            IGRAPH_CHECK(igraph_vector_push_back(times_v, tt + igraph_vector_tail(times_v)));
            IGRAPH_CHECK(igraph_vector_int_push_back(no_s_v, ns));
            IGRAPH_CHECK(igraph_vector_int_push_back(no_i_v, ni));
            IGRAPH_CHECK(igraph_vector_int_push_back(no_r_v, nr));

        } /* psum > 0 */

    } /* j < no_sim */

    RNG_END();

    igraph_psumtree_destroy(&tree);
    igraph_adjlist_destroy(&adjlist);
    igraph_vector_int_destroy(&status);
    IGRAPH_FINALLY_CLEAN(4);  /* + result */

    return IGRAPH_SUCCESS;
}
