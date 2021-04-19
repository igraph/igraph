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

#ifndef IGRAPH_EPIDEMICS_H
#define IGRAPH_EPIDEMICS_H

#include "igraph_decls.h"
#include "igraph_datatype.h"
#include "igraph_vector.h"
#include "igraph_vector_ptr.h"

__BEGIN_DECLS

/**
 * \struct igraph_sir_t
 * \brief The result of one SIR model simulation.
 *
 * Data structure to store the results of one simulation
 * of the SIR (susceptible-infected-recovered) model on a graph.
 *
 * It has the following members. They are all (real or integer)
 * vectors, and they are of the same length.
 *
 * \member times A vector, the times of the events are stored here.
 * \member no_s An integer vector, the number of susceptibles in
 *              each time step is stored here.
 * \member no_i An integer vector, the number of infected individuals
 *              at each time step, is stored here.
 * \member no_r An integer vector, the number of recovered individuals
 *              is stored here at each time step.
 */

typedef struct igraph_sir_t {
    igraph_vector_t times;
    igraph_vector_int_t no_s, no_i, no_r;
} igraph_sir_t;

IGRAPH_EXPORT int igraph_sir_init(igraph_sir_t *sir);
IGRAPH_EXPORT void igraph_sir_destroy(igraph_sir_t *sir);

IGRAPH_EXPORT int igraph_sir(const igraph_t *graph, igraph_real_t beta,
                             igraph_real_t gamma, igraph_integer_t no_sim,
                             igraph_vector_ptr_t *result);

__END_DECLS

#endif
