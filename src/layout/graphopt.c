/* -*- mode: C -*-  */
/* vim:set ts=4 sw=4 sts=4 et: */
/*
   IGraph library.
   Copyright (C) 2003-2020  The igraph development team

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

#include "igraph_layout.h"

#include "igraph_interface.h"
#include "igraph_progress.h"

#include "core/interruption.h"

#define COULOMBS_CONSTANT 8987500000.0


static igraph_real_t igraph_i_distance_between(
        const igraph_matrix_t *c,
        long int a, long int b);

static int igraph_i_determine_electric_axal_forces(
        const igraph_matrix_t *pos,
        igraph_real_t *x,
        igraph_real_t *y,
        igraph_real_t directed_force,
        igraph_real_t distance,
        long int other_node,
        long int this_node);

static int igraph_i_apply_electrical_force(
        const igraph_matrix_t *pos,
        igraph_vector_t *pending_forces_x,
        igraph_vector_t *pending_forces_y,
        long int other_node, long int this_node,
        igraph_real_t node_charge,
        igraph_real_t distance);

static int igraph_i_determine_spring_axal_forces(
        const igraph_matrix_t *pos,
        igraph_real_t *x, igraph_real_t *y,
        igraph_real_t directed_force,
        igraph_real_t distance,
        igraph_real_t spring_length,
        long int other_node,
        long int this_node);

static int igraph_i_apply_spring_force(
        const igraph_matrix_t *pos,
        igraph_vector_t *pending_forces_x,
        igraph_vector_t *pending_forces_y,
        long int other_node,
        long int this_node, igraph_real_t spring_length,
        igraph_real_t spring_constant);

static int igraph_i_move_nodes(
        igraph_matrix_t *pos,
        const igraph_vector_t *pending_forces_x,
        const igraph_vector_t *pending_forces_y,
        igraph_real_t node_mass,
        igraph_real_t max_sa_movement);

static igraph_real_t igraph_i_distance_between(
        const igraph_matrix_t *c,
        long int a, long int b) {
    igraph_real_t diffx = MATRIX(*c, a, 0) - MATRIX(*c, b, 0);
    igraph_real_t diffy = MATRIX(*c, a, 1) - MATRIX(*c, b, 1);
    return sqrt( diffx * diffx + diffy * diffy );
}

static int igraph_i_determine_electric_axal_forces(const igraph_matrix_t *pos,
        igraph_real_t *x,
        igraph_real_t *y,
        igraph_real_t directed_force,
        igraph_real_t distance,
        long int other_node,
        long int this_node) {

    // We know what the directed force is.  We now need to translate it
    // into the appropriate x and y components.
    // First, assume:
    //                 other_node
    //                    /|
    //  directed_force  /  |
    //                /    | y
    //              /______|
    //    this_node     x
    //
    // other_node.x > this_node.x
    // other_node.y > this_node.y
    // the force will be on this_node away from other_node

    // the proportion (distance/y_distance) is equal to the proportion
    // (directed_force/y_force), as the two triangles are similar.
    // therefore, the magnitude of y_force = (directed_force*y_distance)/distance
    // the sign of y_force is negative, away from other_node

    igraph_real_t x_distance, y_distance;
    y_distance = MATRIX(*pos, other_node, 1) - MATRIX(*pos, this_node, 1);
    if (y_distance < 0) {
        y_distance = -y_distance;
    }
    *y = -1 * ((directed_force * y_distance) / distance);

    // the x component works in exactly the same way.
    x_distance = MATRIX(*pos, other_node, 0) - MATRIX(*pos, this_node, 0);
    if (x_distance < 0) {
        x_distance = -x_distance;
    }
    *x = -1 * ((directed_force * x_distance) / distance);

    // Now we need to reverse the polarity of our answers based on the falsness
    // of our assumptions.
    if (MATRIX(*pos, other_node, 0) < MATRIX(*pos, this_node, 0)) {
        *x = *x * -1;
    }
    if (MATRIX(*pos, other_node, 1) < MATRIX(*pos, this_node, 1)) {
        *y = *y * -1;
    }

    return 0;
}

static int igraph_i_apply_electrical_force(
        const igraph_matrix_t *pos,
        igraph_vector_t *pending_forces_x,
        igraph_vector_t *pending_forces_y,
        long int other_node, long int this_node,
        igraph_real_t node_charge,
        igraph_real_t distance) {

    igraph_real_t directed_force = COULOMBS_CONSTANT *
                                   ((node_charge * node_charge) / (distance * distance));

    igraph_real_t x_force, y_force;
    igraph_i_determine_electric_axal_forces(pos, &x_force, &y_force,
                                            directed_force, distance,
                                            other_node, this_node);

    VECTOR(*pending_forces_x)[this_node] += x_force;
    VECTOR(*pending_forces_y)[this_node] += y_force;
    VECTOR(*pending_forces_x)[other_node] -= x_force;
    VECTOR(*pending_forces_y)[other_node] -= y_force;

    return 0;
}

static int igraph_i_determine_spring_axal_forces(
        const igraph_matrix_t *pos,
        igraph_real_t *x, igraph_real_t *y,
        igraph_real_t directed_force,
        igraph_real_t distance,
        igraph_real_t spring_length,
        long int other_node, long int this_node) {

    // if the spring is just the right size, the forces will be 0, so we can
    // skip the computation.
    //
    // if the spring is too long, our forces will be identical to those computed
    // by determine_electrical_axal_forces() (this_node will be pulled toward
    // other_node).
    //
    // if the spring is too short, our forces will be the opposite of those
    // computed by determine_electrical_axal_forces() (this_node will be pushed
    // away from other_node)
    //
    // finally, since both nodes are movable, only one-half of the total force
    // should be applied to each node, so half the forces for our answer.

    if (distance == spring_length) {
        *x = 0.0;
        *y = 0.0;
    } else {
        igraph_i_determine_electric_axal_forces(pos, x, y, directed_force, distance,
                                                other_node, this_node);
        if (distance < spring_length) {
            *x = -1 * *x;
            *y = -1 * *y;
        }
        *x = 0.5 * *x;
        *y = 0.5 * *y;
    }

    return 0;
}

static int igraph_i_apply_spring_force(
        const igraph_matrix_t *pos,
        igraph_vector_t *pending_forces_x,
        igraph_vector_t *pending_forces_y,
        long int other_node,
        long int this_node, igraph_real_t spring_length,
        igraph_real_t spring_constant) {

    // determined using Hooke's Law:
    //   force = -kx
    // where:
    //   k = spring constant
    //   x = displacement from ideal length in meters

    igraph_real_t distance, displacement, directed_force, x_force, y_force;
    distance = igraph_i_distance_between(pos, other_node, this_node);
    // let's protect ourselves from division by zero by ignoring two nodes that
    // happen to be in the same place.  Since we separate all nodes before we
    // work on any of them, this will only happen in extremely rare circumstances,
    // and when it does, electrical force will probably push one or both of them
    // one way or another anyway.
    if (distance == 0.0) {
        return 0;
    }

    displacement = distance - spring_length;
    if (displacement < 0) {
        displacement = -displacement;
    }
    directed_force = -1 * spring_constant * displacement;
    // remember, this is force directed away from the spring;
    // a negative number is back towards the spring (or, in our case, back towards
    // the other node)

    // get the force that should be applied to >this< node
    igraph_i_determine_spring_axal_forces(pos, &x_force, &y_force,
                                          directed_force, distance, spring_length,
                                          other_node, this_node);

    VECTOR(*pending_forces_x)[this_node] += x_force;
    VECTOR(*pending_forces_y)[this_node] += y_force;
    VECTOR(*pending_forces_x)[other_node] -= x_force;
    VECTOR(*pending_forces_y)[other_node] -= y_force;

    return 0;
}

static int igraph_i_move_nodes(
        igraph_matrix_t *pos,
        const igraph_vector_t *pending_forces_x,
        const igraph_vector_t *pending_forces_y,
        igraph_real_t node_mass,
        igraph_real_t max_sa_movement) {

    // Since each iteration is isolated, time is constant at 1.
    // Therefore:
    //   Force effects acceleration.
    //   acceleration (d(velocity)/time) = velocity
    //   velocity (d(displacement)/time) = displacement
    //   displacement = acceleration

    // determined using Newton's second law:
    //   sum(F) = ma
    // therefore:
    //   acceleration = force / mass
    //   velocity     = force / mass
    //   displacement = force / mass

    long int this_node, no_of_nodes = igraph_vector_size(pending_forces_x);

    for (this_node = 0; this_node < no_of_nodes; this_node++) {

        igraph_real_t x_movement, y_movement;

        x_movement = VECTOR(*pending_forces_x)[this_node] / node_mass;
        if (x_movement > max_sa_movement) {
            x_movement = max_sa_movement;
        } else if (x_movement < -max_sa_movement) {
            x_movement = -max_sa_movement;
        }

        y_movement = VECTOR(*pending_forces_y)[this_node] / node_mass;
        if (y_movement > max_sa_movement) {
            y_movement = max_sa_movement;
        } else if (y_movement < -max_sa_movement) {
            y_movement = -max_sa_movement;
        }

        MATRIX(*pos, this_node, 0) += x_movement;
        MATRIX(*pos, this_node, 1) += y_movement;

    }
    return 0;
}

/**
 * \function igraph_layout_graphopt
 * \brief Optimizes vertex layout via the graphopt algorithm.
 *
 * </para><para>
 * This is a port of the graphopt layout algorithm by Michael Schmuhl.
 * graphopt version 0.4.1 was rewritten in C and the support for
 * layers was removed (might be added later) and a code was a bit
 * reorganized to avoid some unnecessary steps is the node charge (see below)
 * is zero.
 *
 * </para><para>
 * Graphopt uses physical analogies for defining attracting and repelling
 * forces among the vertices and then the physical system is simulated
 * until it reaches an equilibrium. (There is no simulated annealing or
 * anything like that, so a stable fixed point is not guaranteed.)
 *
 * </para><para>
 * See also http://www.schmuhl.org/graphopt/ for the original graphopt.
 * \param graph The input graph.
 * \param res Pointer to an initialized matrix, the result will be stored here
 *    and its initial contents are used as the starting point of the simulation
 *    if the \p use_seed argument is true. Note that in this case the
 *    matrix should have the proper size, otherwise a warning is issued and
 *    the supplied values are ignored. If no starting positions are given
 *    (or they are invalid) then a random starting position is used.
 *    The matrix will be resized if needed.
 * \param niter Integer constant, the number of iterations to perform.
 *    Should be a couple of hundred in general. If you have a large graph
 *    then you might want to only do a few iterations and then check the
 *    result. If it is not good enough you can feed it in again in
 *    the \p res argument. The original graphopt default is 500.
 * \param node_charge The charge of the vertices, used to calculate electric
 *    repulsion. The original graphopt default is 0.001.
 * \param node_mass The mass of the vertices, used for the spring forces.
 *    The original graphopt defaults to 30.
 * \param spring_length The length of the springs.
 *    The original graphopt defaults to zero.
 * \param spring_constant The spring constant, the original graphopt defaults
 *    to one.
 * \param max_sa_movement Real constant, it gives the maximum amount of movement
 *    allowed in a single step along a single axis. The original graphopt
 *    default is 5.
 * \param use_seed Logical scalar, whether to use the positions in \p res as
 *    a starting configuration. See also \p res above.
 * \return Error code.
 *
 * Time complexity: O(n (|V|^2+|E|) ), n is the number of iterations,
 * |V| is the number of vertices, |E| the number
 * of edges. If \p node_charge is zero then it is only O(n|E|).
 */
int igraph_layout_graphopt(const igraph_t *graph, igraph_matrix_t *res,
                           igraph_integer_t niter,
                           igraph_real_t node_charge, igraph_real_t node_mass,
                           igraph_real_t spring_length,
                           igraph_real_t spring_constant,
                           igraph_real_t max_sa_movement,
                           igraph_bool_t use_seed) {

    long int no_of_nodes = igraph_vcount(graph);
    long int no_of_edges = igraph_ecount(graph);
    igraph_vector_t pending_forces_x, pending_forces_y;
    /* Set a flag to calculate (or not) the electrical forces that the nodes */
    /* apply on each other based on if both node types' charges are zero. */
    igraph_bool_t apply_electric_charges = (node_charge != 0);

    long int this_node, other_node, edge;
    igraph_real_t distance;
    long int i;

    IGRAPH_VECTOR_INIT_FINALLY(&pending_forces_x, no_of_nodes);
    IGRAPH_VECTOR_INIT_FINALLY(&pending_forces_y, no_of_nodes);

    if (use_seed) {
        if (igraph_matrix_nrow(res) != no_of_nodes ||
            igraph_matrix_ncol(res) != 2) {
            IGRAPH_WARNING("Invalid size for initial matrix, starting from random layout.");
            IGRAPH_CHECK(igraph_layout_random(graph, res));
        }
    } else {
        IGRAPH_CHECK(igraph_layout_random(graph, res));
    }

    IGRAPH_PROGRESS("Graphopt layout", 0, NULL);
    for (i = niter; i > 0; i--) {
        /* Report progress in approx. every 100th step */
        if (i % 10 == 0) {
            IGRAPH_PROGRESS("Graphopt layout", 100.0 - 100.0 * i / niter, NULL);
        }

        /* Clear pending forces on all nodes */
        igraph_vector_null(&pending_forces_x);
        igraph_vector_null(&pending_forces_y);

        // Apply electrical force applied by all other nodes
        if (apply_electric_charges) {
            // Iterate through all nodes
            for (this_node = 0; this_node < no_of_nodes; this_node++) {
                IGRAPH_ALLOW_INTERRUPTION();
                for (other_node = this_node + 1;
                     other_node < no_of_nodes;
                     other_node++) {
                    distance = igraph_i_distance_between(res, this_node, other_node);
                    // let's protect ourselves from division by zero by ignoring
                    // two nodes that happen to be in the same place.  Since we
                    // separate all nodes before we work on any of them, this
                    // will only happen in extremely rare circumstances, and when
                    // it does, springs will probably pull them apart anyway.
                    // also, if we are more than 50 away, the electric force
                    // will be negligible.
                    // ***** may not always be desirable ****
                    if ((distance != 0.0) && (distance < 500.0)) {
                        //    if (distance != 0.0) {
                        // Apply electrical force from node(counter2) on
                        // node(counter)
                        igraph_i_apply_electrical_force(res, &pending_forces_x,
                                                        &pending_forces_y,
                                                        other_node, this_node,
                                                        node_charge,
                                                        distance);
                    }
                }
            }
        }

        // Apply force from springs
        for (edge = 0; edge < no_of_edges; edge++) {
            long int tthis_node = IGRAPH_FROM(graph, edge);
            long int oother_node = IGRAPH_TO(graph, edge);
            // Apply spring force on both nodes
            igraph_i_apply_spring_force(res, &pending_forces_x, &pending_forces_y,
                                        oother_node, tthis_node, spring_length,
                                        spring_constant);
        }

        // Effect the movement of the nodes based on all pending forces
        igraph_i_move_nodes(res, &pending_forces_x, &pending_forces_y, node_mass,
                            max_sa_movement);
    }
    IGRAPH_PROGRESS("Graphopt layout", 100, NULL);

    igraph_vector_destroy(&pending_forces_y);
    igraph_vector_destroy(&pending_forces_x);
    IGRAPH_FINALLY_CLEAN(2);

    return IGRAPH_SUCCESS;
}
