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
#include "igraph_random.h"
#include "igraph_structural.h"
#include "igraph_visitor.h"

#include "core/grid.h"
#include "core/interruption.h"
#include "core/math.h"

static void igraph_i_norm2d(igraph_real_t *x, igraph_real_t *y) {
    igraph_real_t len = sqrt(*x * *x + *y * *y);
    if (len != 0) {
        *x /= len;
        *y /= len;
    }
}

/**
 * \function igraph_layout_lgl
 * \brief Force based layout algorithm for large graphs.
 *
 * </para><para>
 * This is a layout generator similar to the Large Graph Layout
 * algorithm and program (http://lgl.sourceforge.net/). But unlike LGL, this
 * version uses a Fruchterman-Reingold style simulated annealing
 * algorithm for placing the vertices. The speedup is achieved by
 * placing the vertices on a grid and calculating the repulsion only
 * for vertices which are closer to each other than a limit.
 *
 * \param graph The (initialized) graph object to place. It must be connnected;
 *   disconnected graphs are not handled by the algorithm.
 * \param res Pointer to an initialized matrix object to hold the
 *   result. It will be resized if needed.
 * \param maxit The maximum number of cooling iterations to perform
 *   for each layout step. A reasonable default is 150.
 * \param maxdelta The maximum length of the move allowed for a vertex
 *   in a single iteration. A reasonable default is the number of
 *   vertices.
 * \param area This parameter gives the area of the square on which
 *   the vertices will be placed. A reasonable default value is the
 *   number of vertices squared.
 * \param coolexp The cooling exponent. A reasonable default value is
 *   1.5.
 * \param repulserad Determines the radius at which vertex-vertex
 *   repulsion cancels out attraction of adjacent vertices. A
 *   reasonable default value is \p area times the number of vertices.
 * \param cellsize The size of the grid cells, one side of the
 *   square. A reasonable default value is the fourth root of
 *   \p area (or the square root of the number of vertices if \p area
 *   is also left at its default value).
 * \param proot The root vertex, this is placed first, its neighbors
 *   in the first iteration, second neighbors in the second, etc. If
 *   negative then a random vertex is chosen.
 * \return Error code.
 *
 * Added in version 0.2.</para><para>
 *
 * Time complexity: ideally O(dia*maxit*(|V|+|E|)), |V| is the number
 * of vertices,
 * dia is the diameter of the graph, worst case complexity is still
 * O(dia*maxit*(|V|^2+|E|)), this is the case when all vertices happen to be
 * in the same grid cell.
 */

igraph_error_t igraph_layout_lgl(const igraph_t *graph, igraph_matrix_t *res,
                      igraph_integer_t maxit, igraph_real_t maxdelta,
                      igraph_real_t area, igraph_real_t coolexp,
                      igraph_real_t repulserad, igraph_real_t cellsize,
                      igraph_integer_t proot) {


    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_integer_t no_of_edges = igraph_ecount(graph);
    igraph_t mst;
    igraph_integer_t root;
    igraph_integer_t no_of_layers, actlayer = 0;
    igraph_vector_int_t vids;
    igraph_vector_int_t layers;
    igraph_vector_int_t parents;
    igraph_vector_int_t edges;
    igraph_2dgrid_t grid;
    igraph_vector_int_t eids;
    igraph_vector_t forcex;
    igraph_vector_t forcey;

    igraph_real_t frk = sqrt(area / no_of_nodes);
    igraph_real_t H_n = 0;

    if (no_of_nodes == 0) {
        /* We skip parameter checks for the null graph, as following the recommendations
         * for parameter choices in the documentation would lead to zero values that are
         * considered invalid in general, but don't cause problems for the null graph. */
        IGRAPH_CHECK(igraph_matrix_resize(res, 0, 2));
        return IGRAPH_SUCCESS;
    }

    /* TODO: is zero okay? */
    if (maxit < 0) {
        IGRAPH_ERRORF("Maximum number of iterations must not be negative, got %" IGRAPH_PRId ".", IGRAPH_EINVAL, maxit);
    }

    if (maxdelta <= 0) {
        IGRAPH_ERRORF("Maximum delta must be positive, got %g.", IGRAPH_EINVAL, maxdelta);
    }

    if (area <= 0) {
        IGRAPH_ERRORF("Placement area size must be positive, got %g.", IGRAPH_EINVAL, area);
    }

    if (coolexp <= 0) {
        IGRAPH_ERRORF("Cooling exponent must be positive, got %g.", IGRAPH_EINVAL, coolexp);
    }

    if (repulserad <= 0) {
        IGRAPH_ERRORF("Repusion cutoff radius must be positive, got %g.", IGRAPH_EINVAL, repulserad);
    }

    if (cellsize <= 0) {
        IGRAPH_ERRORF("Cell size must be positive, got %g.", IGRAPH_EINVAL, cellsize);
    }

    IGRAPH_CHECK(igraph_minimum_spanning_tree_unweighted(graph, &mst));
    IGRAPH_FINALLY(igraph_destroy, &mst);

    IGRAPH_CHECK(igraph_matrix_resize(res, no_of_nodes, 2));

    /* Determine the root vertex, random pick right now */
    if (proot < 0) {
        root = RNG_INTEGER(0, no_of_nodes - 1);
    } else {
        root = proot;
    }

    /* Assign the layers */
    IGRAPH_VECTOR_INT_INIT_FINALLY(&vids, 0);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&layers, 0);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&parents, 0);
    if (no_of_nodes > 0) {
        IGRAPH_CHECK(igraph_bfs_simple(&mst, root, IGRAPH_ALL, &vids, &layers, &parents));
    }
    no_of_layers = igraph_vector_int_size(&layers) - 1;

    /* Check whether we have reached all the nodes -- if not, the graph is
     * disconnected */
    if (no_of_nodes > 0 && igraph_vector_int_min(&parents) <= -2) {
        IGRAPH_WARNING("LGL layout does not support disconnected graphs yet.");
    }

    /* We don't need the mst any more */
    igraph_destroy(&mst);
    igraph_empty(&mst, 0, IGRAPH_UNDIRECTED); /* to make finalization work */

    IGRAPH_VECTOR_INT_INIT_FINALLY(&edges, 0);
    IGRAPH_CHECK(igraph_vector_int_reserve(&edges, no_of_edges));
    IGRAPH_VECTOR_INT_INIT_FINALLY(&eids, 0);
    IGRAPH_VECTOR_INIT_FINALLY(&forcex, no_of_nodes);
    IGRAPH_VECTOR_INIT_FINALLY(&forcey, no_of_nodes);

    /* Place the vertices randomly */
    IGRAPH_CHECK(igraph_layout_random(graph, res));
    igraph_matrix_scale(res, 1e6);

    /* This is the grid for calculating the vertices near to a given vertex */
    IGRAPH_CHECK(igraph_2dgrid_init(&grid, res,
                                    -sqrt(area / M_PI), sqrt(area / M_PI), cellsize,
                                    -sqrt(area / M_PI), sqrt(area / M_PI), cellsize));
    IGRAPH_FINALLY(igraph_2dgrid_destroy, &grid);

    /* Place the root vertex */
    igraph_2dgrid_add(&grid, root, 0, 0);

    for (actlayer = 1; actlayer < no_of_layers; actlayer++) {
        H_n += 1.0 / actlayer;
    }

    for (actlayer = 1; actlayer < no_of_layers; actlayer++) {

        igraph_real_t c = 1;
        igraph_integer_t i, j;
        igraph_real_t massx, massy;
        igraph_real_t px, py;
        igraph_real_t sx, sy;

        igraph_integer_t it = 0;
        igraph_real_t epsilon = 10e-6;
        igraph_real_t maxchange = epsilon + 1;
        /* igraph_integer_t pairs; */
        igraph_real_t sconst = sqrt(area / M_PI) / H_n;
        igraph_2dgrid_iterator_t vidit;

        /*     printf("Layer %li:\n", actlayer); */

        /*-----------------------------------------*/
        /* Step 1: place the next layer on spheres */
        /*-----------------------------------------*/

        RNG_BEGIN();

        j = VECTOR(layers)[actlayer];
        for (i = VECTOR(layers)[actlayer - 1];
             i < VECTOR(layers)[actlayer]; i++) {

            igraph_integer_t vid = VECTOR(vids)[i];
            igraph_integer_t par = VECTOR(parents)[vid];

            if (par < 0) {
                /* this is either the root vertex or an unreachable node */
                continue;
            }

            IGRAPH_ALLOW_INTERRUPTION();
            igraph_2dgrid_getcenter(&grid, &massx, &massy);
            igraph_i_norm2d(&massx, &massy);
            px = MATRIX(*res, vid, 0) - MATRIX(*res, par, 0);
            py = MATRIX(*res, vid, 1) - MATRIX(*res, par, 1);
            igraph_i_norm2d(&px, &py);
            sx = c * (massx + px) + MATRIX(*res, vid, 0);
            sy = c * (massy + py) + MATRIX(*res, vid, 1);

            /* The neighbors of 'vid' */
            while (j < VECTOR(layers)[actlayer + 1] && VECTOR(parents)[VECTOR(vids)[j]] == vid) {
                igraph_real_t rx, ry;
                if (actlayer == 1) {
                    igraph_real_t phi = 2 * M_PI / (VECTOR(layers)[2] - 1) * (j - 1);
                    rx = cos(phi);
                    ry = sin(phi);
                } else {
                    rx = RNG_UNIF(-1, 1);
                    ry = RNG_UNIF(-1, 1);
                }
                igraph_i_norm2d(&rx, &ry);
                rx = rx / actlayer * sconst;
                ry = ry / actlayer * sconst;
                igraph_2dgrid_add(&grid, VECTOR(vids)[j], sx + rx, sy + ry);
                j++;
            }
        }

        RNG_END();

        /*-----------------------------------------*/
        /* Step 2: add the edges of the next layer */
        /*-----------------------------------------*/

        for (j = VECTOR(layers)[actlayer];
             j < VECTOR(layers)[actlayer + 1]; j++) {
            igraph_integer_t vid = VECTOR(vids)[j];
            igraph_integer_t k;
            IGRAPH_ALLOW_INTERRUPTION();
            IGRAPH_CHECK(igraph_incident(graph, &eids, vid, IGRAPH_ALL));
            for (k = 0; k < igraph_vector_int_size(&eids); k++) {
                igraph_integer_t eid = VECTOR(eids)[k];
                igraph_integer_t from = IGRAPH_FROM(graph, eid), to = IGRAPH_TO(graph, eid);
                if ((from != vid && igraph_2dgrid_in(&grid, from)) ||
                    (to   != vid && igraph_2dgrid_in(&grid, to))) {
                    igraph_vector_int_push_back(&edges, eid);
                }
            }
        }

        /*-----------------------------------------*/
        /* Step 3: let the springs spring          */
        /*-----------------------------------------*/

        maxchange = epsilon + 1;
        while (it < maxit && maxchange > epsilon) {
            igraph_integer_t jj;
            igraph_real_t t = maxdelta * pow((maxit - it) / (igraph_real_t) maxit, coolexp);
            igraph_integer_t vid, nei;

            IGRAPH_PROGRESS("Large graph layout",
                            100.0 * ((actlayer - 1.0) / (no_of_layers - 1.0) + (it) / (maxit * (no_of_layers - 1.0))),
                            0);

            /* init */
            igraph_vector_null(&forcex);
            igraph_vector_null(&forcey);
            maxchange = 0;

            /* attractive "forces" along the edges */
            for (jj = 0; jj < igraph_vector_int_size(&edges); jj++) {
                igraph_integer_t from = IGRAPH_FROM(graph, VECTOR(edges)[jj]);
                igraph_integer_t to = IGRAPH_TO(graph, VECTOR(edges)[jj]);
                igraph_real_t xd, yd, dist, force;
                IGRAPH_ALLOW_INTERRUPTION();
                xd = MATRIX(*res, from, 0) - MATRIX(*res, to, 0);
                yd = MATRIX(*res, from, 1) - MATRIX(*res, to, 1);
                dist = sqrt(xd*xd + yd*yd);
                if (dist != 0) {
                    xd /= dist;
                    yd /= dist;
                }
                force = dist * dist / frk;
                VECTOR(forcex)[from] -= xd * force;
                VECTOR(forcex)[to]   += xd * force;
                VECTOR(forcey)[from] -= yd * force;
                VECTOR(forcey)[to]   += yd * force;
            }

            /* repulsive "forces" of the vertices nearby */
            /* pairs = 0; */
            igraph_2dgrid_reset(&grid, &vidit);
            while ( (vid = igraph_2dgrid_next(&grid, &vidit) - 1) != -1) {
                while ( (nei = igraph_2dgrid_next_nei(&grid, &vidit) - 1) != -1) {
                    igraph_real_t xd = MATRIX(*res, vid, 0) - MATRIX(*res, nei, 0);
                    igraph_real_t yd = MATRIX(*res, vid, 1) - MATRIX(*res, nei, 1);
                    igraph_real_t dist = sqrt(xd*xd + yd*yd);
                    igraph_real_t force;
                    if (dist < cellsize) {
                        /* pairs++; */
                        if (dist == 0) {
                            dist = epsilon;
                        };
                        xd /= dist; yd /= dist;
                        force = frk * frk * (1.0 / dist - dist * dist / repulserad);
                        VECTOR(forcex)[vid] += xd * force;
                        VECTOR(forcex)[nei] -= xd * force;
                        VECTOR(forcey)[vid] += yd * force;
                        VECTOR(forcey)[nei] -= yd * force;
                    }
                }
            }

            /*       printf("verties: %li iterations: %li\n",  */
            /*       VECTOR(layers)[actlayer+1], pairs); */

            /* apply the changes */
            for (jj = 0; jj < VECTOR(layers)[actlayer + 1]; jj++) {
                igraph_integer_t vvid = VECTOR(vids)[jj];
                igraph_real_t fx = VECTOR(forcex)[vvid];
                igraph_real_t fy = VECTOR(forcey)[vvid];
                igraph_real_t ded = sqrt(fx*fx + fy*fy);
                if (ded > t) {
                    ded = t / ded;
                    fx *= ded; fy *= ded;
                }
                igraph_2dgrid_move(&grid, vvid, fx, fy);
                if (fx > maxchange) {
                    maxchange = fx;
                }
                if (fy > maxchange) {
                    maxchange = fy;
                }
            }
            it++;
            /*       printf("%li iterations, maxchange: %f\n", it, (double)maxchange); */
        }
    }

    IGRAPH_PROGRESS("Large graph layout", 100.0, 0);
    igraph_destroy(&mst);
    igraph_vector_int_destroy(&vids);
    igraph_vector_int_destroy(&layers);
    igraph_vector_int_destroy(&parents);
    igraph_vector_int_destroy(&edges);
    igraph_2dgrid_destroy(&grid);
    igraph_vector_int_destroy(&eids);
    igraph_vector_destroy(&forcex);
    igraph_vector_destroy(&forcey);
    IGRAPH_FINALLY_CLEAN(9);
    return IGRAPH_SUCCESS;

}
