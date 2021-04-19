/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2009-2012  Gabor Csardi <csardi.gabor@gmail.com>
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
   Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA
   02110-1301 USA

*/

#ifndef IGRAPH_LAYOUT_H
#define IGRAPH_LAYOUT_H

#include "igraph_decls.h"
#include "igraph_constants.h"
#include "igraph_types.h"
#include "igraph_vector.h"
#include "igraph_vector_ptr.h"
#include "igraph_matrix.h"
#include "igraph_datatype.h"
#include "igraph_arpack.h"
#include "igraph_iterators.h"

__BEGIN_DECLS

/**
 * \section about_layouts
 *
 * <para>Layout generator functions (or at least most of them) try to place the
 * vertices and edges of a graph on a 2D plane or in 3D space in a way
 * which visually pleases the human eye.</para>
 *
 * <para>They take a graph object and a number of parameters as arguments
 * and return an \type igraph_matrix_t, in which each row gives the
 * coordinates of a vertex.</para>
 */

/* -------------------------------------------------- */
/* Layouts                                            */
/* -------------------------------------------------- */

IGRAPH_EXPORT int igraph_layout_random(const igraph_t *graph, igraph_matrix_t *res);
IGRAPH_EXPORT int igraph_layout_circle(const igraph_t *graph, igraph_matrix_t *res,
                                       igraph_vs_t order);
IGRAPH_EXPORT int igraph_layout_star(const igraph_t *graph, igraph_matrix_t *res,
                                     igraph_integer_t center, const igraph_vector_t *order);
IGRAPH_EXPORT int igraph_layout_grid(const igraph_t *graph, igraph_matrix_t *res, long int width);
IGRAPH_EXPORT int igraph_layout_fruchterman_reingold(const igraph_t *graph,
                                                     igraph_matrix_t *res,
                                                     igraph_bool_t use_seed,
                                                     igraph_integer_t niter,
                                                     igraph_real_t start_temp,
                                                     igraph_layout_grid_t grid,
                                                     const igraph_vector_t *weight,
                                                     const igraph_vector_t *minx,
                                                     const igraph_vector_t *maxx,
                                                     const igraph_vector_t *miny,
                                                     const igraph_vector_t *maxy);

IGRAPH_EXPORT int igraph_layout_kamada_kawai(const igraph_t *graph, igraph_matrix_t *res,
                                             igraph_bool_t use_seed, igraph_integer_t maxiter,
                                             igraph_real_t epsilon, igraph_real_t kkconst,
                                             const igraph_vector_t *weights,
                                             const igraph_vector_t *minx, const igraph_vector_t *maxx,
                                             const igraph_vector_t *miny, const igraph_vector_t *maxy);

IGRAPH_EXPORT int igraph_layout_lgl(const igraph_t *graph, igraph_matrix_t *res,
                                    igraph_integer_t maxiter, igraph_real_t maxdelta,
                                    igraph_real_t area, igraph_real_t coolexp,
                                    igraph_real_t repulserad, igraph_real_t cellsize, igraph_integer_t root);
IGRAPH_EXPORT int igraph_layout_reingold_tilford(const igraph_t *graph, igraph_matrix_t *res,
                                                 igraph_neimode_t mode,
                                                 const igraph_vector_t *roots,
                                                 const igraph_vector_t *rootlevel);
IGRAPH_EXPORT int igraph_layout_reingold_tilford_circular(const igraph_t *graph,
                                                          igraph_matrix_t *res,
                                                          igraph_neimode_t mode,
                                                          const igraph_vector_t *roots,
                                                          const igraph_vector_t *rootlevel);
IGRAPH_EXPORT int igraph_layout_sugiyama(const igraph_t *graph, igraph_matrix_t *res,
                                         igraph_t *extd_graph, igraph_vector_t *extd_to_orig_eids,
                                         const igraph_vector_t* layers, igraph_real_t hgap,
                                         igraph_real_t vgap, long int maxiter, const igraph_vector_t *weights);

IGRAPH_EXPORT int igraph_layout_random_3d(const igraph_t *graph, igraph_matrix_t *res);
IGRAPH_EXPORT int igraph_layout_sphere(const igraph_t *graph, igraph_matrix_t *res);
IGRAPH_EXPORT int igraph_layout_grid_3d(const igraph_t *graph, igraph_matrix_t *res,
                                        long int width, long int height);
IGRAPH_EXPORT int igraph_layout_fruchterman_reingold_3d(const igraph_t *graph,
                                                        igraph_matrix_t *res,
                                                        igraph_bool_t use_seed,
                                                        igraph_integer_t niter,
                                                        igraph_real_t start_temp,
                                                        const igraph_vector_t *weight,
                                                        const igraph_vector_t *minx,
                                                        const igraph_vector_t *maxx,
                                                        const igraph_vector_t *miny,
                                                        const igraph_vector_t *maxy,
                                                        const igraph_vector_t *minz,
                                                        const igraph_vector_t *maxz);

IGRAPH_EXPORT int igraph_layout_kamada_kawai_3d(const igraph_t *graph, igraph_matrix_t *res,
                                                igraph_bool_t use_seed, igraph_integer_t maxiter,
                                                igraph_real_t epsilon, igraph_real_t kkconst,
                                                const igraph_vector_t *weights,
                                                const igraph_vector_t *minx, const igraph_vector_t *maxx,
                                                const igraph_vector_t *miny, const igraph_vector_t *maxy,
                                                const igraph_vector_t *minz, const igraph_vector_t *maxz);

IGRAPH_EXPORT int igraph_layout_graphopt(const igraph_t *graph,
                                         igraph_matrix_t *res, igraph_integer_t niter,
                                         igraph_real_t node_charge, igraph_real_t node_mass,
                                         igraph_real_t spring_length,
                                         igraph_real_t spring_constant,
                                         igraph_real_t max_sa_movement,
                                         igraph_bool_t use_seed);

IGRAPH_EXPORT int igraph_layout_mds(const igraph_t *graph, igraph_matrix_t *res,
                                    const igraph_matrix_t *dist, long int dim);

IGRAPH_EXPORT int igraph_layout_bipartite(const igraph_t *graph,
                                          const igraph_vector_bool_t *types,
                                          igraph_matrix_t *res, igraph_real_t hgap,
                                          igraph_real_t vgap, long int maxiter);

/**
 * \struct igraph_layout_drl_options_t
 * Parameters for the DrL layout generator
 *
 * \member edge_cut The edge cutting parameter.
 *    Edge cutting is done in the late stages of the
 *    algorithm in order to achieve less dense layouts.  Edges are cut
 *    if there is a lot of stress on them (a large value in the
 *    objective function sum).  The edge cutting parameter is a value
 *    between 0 and 1 with 0 representing no edge cutting and 1
 *    representing maximal edge cutting. The default value is 32/40.
 * \member init_iterations Number of iterations, initial phase.
 * \member init_temperature Start temperature, initial phase.
 * \member init_attraction Attraction, initial phase.
 * \member init_damping_mult Damping factor, initial phase.
 * \member liquid_iterations Number of iterations in the liquid phase.
 * \member liquid_temperature Start temperature in the liquid phase.
 * \member liquid_attraction Attraction in the liquid phase.
 * \member liquid_damping_mult Multiplicatie damping factor, liquid phase.
 * \member expansion_iterations Number of iterations in the expansion phase.
 * \member expansion_temperature Start temperature in the expansion phase.
 * \member expansion_attraction Attraction, expansion phase.
 * \member expansion_damping_mult Damping factor, expansion phase.
 * \member cooldown_iterations Number of iterations in the cooldown phase.
 * \member cooldown_temperature Start temperature in the cooldown phase.
 * \member cooldown_attraction Attraction in the cooldown phase.
 * \member cooldown_damping_mult Damping fact int the cooldown phase.
 * \member crunch_iterations Number of iterations in the crunch phase.
 * \member crunch_temperature Start temperature in the crunch phase.
 * \member crunch_attraction Attraction in the crunch phase.
 * \member crunch_damping_mult Damping factor in the crunch phase.
 * \member simmer_iterations Number of iterations in the simmer phase.
 * \member simmer_temperature Start temperature in te simmer phase.
 * \member simmer_attraction Attraction in the simmer phase.
 * \member simmer_damping_mult Multiplicative damping factor in the simmer phase.
 */

typedef struct igraph_layout_drl_options_t {
    igraph_real_t    edge_cut;
    igraph_integer_t init_iterations;
    igraph_real_t    init_temperature;
    igraph_real_t    init_attraction;
    igraph_real_t    init_damping_mult;
    igraph_integer_t liquid_iterations;
    igraph_real_t    liquid_temperature;
    igraph_real_t    liquid_attraction;
    igraph_real_t    liquid_damping_mult;
    igraph_integer_t expansion_iterations;
    igraph_real_t    expansion_temperature;
    igraph_real_t    expansion_attraction;
    igraph_real_t    expansion_damping_mult;
    igraph_integer_t cooldown_iterations;
    igraph_real_t    cooldown_temperature;
    igraph_real_t    cooldown_attraction;
    igraph_real_t    cooldown_damping_mult;
    igraph_integer_t crunch_iterations;
    igraph_real_t    crunch_temperature;
    igraph_real_t    crunch_attraction;
    igraph_real_t    crunch_damping_mult;
    igraph_integer_t simmer_iterations;
    igraph_real_t    simmer_temperature;
    igraph_real_t    simmer_attraction;
    igraph_real_t    simmer_damping_mult;
} igraph_layout_drl_options_t;

/**
 * \typedef igraph_layout_drl_default_t
 * Predefined parameter templates for the DrL layout generator
 *
 * These constants can be used to initialize a set of DrL parameters.
 * These can then be modified according to the user's needs.
 * \enumval IGRAPH_LAYOUT_DRL_DEFAULT The deafult parameters.
 * \enumval IGRAPH_LAYOUT_DRL_COARSEN Slightly modified parameters to
 *      get a coarser layout.
 * \enumval IGRAPH_LAYOUT_DRL_COARSEST An even coarser layout.
 * \enumval IGRAPH_LAYOUT_DRL_REFINE Refine an already calculated layout.
 * \enumval IGRAPH_LAYOUT_DRL_FINAL Finalize an already refined layout.
 */

typedef enum { IGRAPH_LAYOUT_DRL_DEFAULT = 0,
               IGRAPH_LAYOUT_DRL_COARSEN,
               IGRAPH_LAYOUT_DRL_COARSEST,
               IGRAPH_LAYOUT_DRL_REFINE,
               IGRAPH_LAYOUT_DRL_FINAL
             } igraph_layout_drl_default_t;

IGRAPH_EXPORT int igraph_layout_drl_options_init(igraph_layout_drl_options_t *options,
                                                 igraph_layout_drl_default_t templ);
IGRAPH_EXPORT int igraph_layout_drl(const igraph_t *graph, igraph_matrix_t *res,
                                    igraph_bool_t use_seed,
                                    igraph_layout_drl_options_t *options,
                                    const igraph_vector_t *weights,
                                    const igraph_vector_bool_t *fixed);

IGRAPH_EXPORT int igraph_layout_drl_3d(const igraph_t *graph, igraph_matrix_t *res,
                                       igraph_bool_t use_seed,
                                       igraph_layout_drl_options_t *options,
                                       const igraph_vector_t *weights,
                                       const igraph_vector_bool_t *fixed);

IGRAPH_EXPORT int igraph_layout_merge_dla(igraph_vector_ptr_t *graphs,
                                          igraph_vector_ptr_t *coords,
                                          igraph_matrix_t *res);

IGRAPH_EXPORT int igraph_layout_gem(const igraph_t *graph, igraph_matrix_t *res,
                                    igraph_bool_t use_seed, igraph_integer_t maxiter,
                                    igraph_real_t temp_max, igraph_real_t temp_min,
                                    igraph_real_t temp_init);

IGRAPH_EXPORT int igraph_layout_davidson_harel(const igraph_t *graph, igraph_matrix_t *res,
                                               igraph_bool_t use_seed, igraph_integer_t maxiter,
                                               igraph_integer_t fineiter, igraph_real_t cool_fact,
                                               igraph_real_t weight_node_dist, igraph_real_t weight_border,
                                               igraph_real_t weight_edge_lengths,
                                               igraph_real_t weight_edge_crossings,
                                               igraph_real_t weight_node_edge_dist);

__END_DECLS

#endif
