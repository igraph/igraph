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
// This file contains the member definitions of the DensityGrid.h class
// This code is modified from the original code by B.N. Wylie

#include "drl_Node_3d.h"
#include "DensityGrid_3d.h"
#include "igraph_error.h"

#include <deque>
#include <cmath>

using namespace std;

#define GET_BIN(z, y, x) (Bins[(z*GRID_SIZE+y)*GRID_SIZE+x])

namespace drl3d {

//*******************************************************
// Density Grid Destructor -- deallocates memory used
// for Density matrix, fall_off matrix, and node deque.

DensityGrid::~DensityGrid () {
    delete[] Density;
    delete[] fall_off;
    delete[] Bins;
}

/*********************************************
* Function: Density_Grid::Reset              *
* Description: Reset the density grid        *
*********************************************/
// changed from reset to init since we will only
// call this once in the parallel version of layout

void DensityGrid::Init() {

    try {
        Density = new float[GRID_SIZE][GRID_SIZE][GRID_SIZE];
        fall_off = new float[RADIUS * 2 + 1][RADIUS * 2 + 1][RADIUS * 2 + 1];
        Bins = new deque<Node>[GRID_SIZE * GRID_SIZE * GRID_SIZE];
    } catch (bad_alloc&) {
        // cout << "Error: Out of memory! Program stopped." << endl;
#ifdef MUSE_MPI
        MPI_Abort ( MPI_COMM_WORLD, 1 );
#else
        igraph_error("DrL is out of memory", IGRAPH_FILE_BASENAME, __LINE__,
                     IGRAPH_ENOMEM);
        return;
#endif
    }

    // Clear Grid
    int i;
    for (i = 0; i < GRID_SIZE; i++)
        for (int j = 0; j < GRID_SIZE; j++)
            for (int k = 0; k < GRID_SIZE; k++) {
                Density[i][j][k] = 0;
                GET_BIN(i, j, k).erase(GET_BIN(i, j, k).begin(), GET_BIN(i, j, k).end());
            }

    // Compute fall off
    for (i = -RADIUS; i <= RADIUS; i++)
        for (int j = -RADIUS; j <= RADIUS; j++)
            for (int k = -RADIUS; k <= RADIUS; k++) {
                fall_off[i + RADIUS][j + RADIUS][k + RADIUS] =
                    (float)((RADIUS - fabs((float)i)) / RADIUS) *
                    (float)((RADIUS - fabs((float)j)) / RADIUS) *
                    (float)((RADIUS - fabs((float)k)) / RADIUS);
            }

}


/***************************************************
 * Function: DensityGrid::GetDensity               *
 * Description: Get_Density from density grid      *
 **************************************************/
float DensityGrid::GetDensity(float Nx, float Ny, float Nz, bool fineDensity) {
    deque<Node>::iterator BI;
    int x_grid, y_grid, z_grid;
    float x_dist, y_dist, z_dist, distance, density = 0;
    int boundary = 10;  // boundary around plane


    /* Where to look */
    x_grid = (int)((Nx + HALF_VIEW + .5) * VIEW_TO_GRID);
    y_grid = (int)((Ny + HALF_VIEW + .5) * VIEW_TO_GRID);
    z_grid = (int)((Nz + HALF_VIEW + .5) * VIEW_TO_GRID);

    // Check for edges of density grid (10000 is arbitrary high density)
    if (x_grid > GRID_SIZE - boundary || x_grid < boundary) {
        return 10000;
    }
    if (y_grid > GRID_SIZE - boundary || y_grid < boundary) {
        return 10000;
    }
    if (z_grid > GRID_SIZE - boundary || z_grid < boundary) {
        return 10000;
    }

    // Fine density?
    if (fineDensity) {

        // Go through nearest bins
        for (int k = z_grid - 1; k <= z_grid + 1; k++)
            for (int i = y_grid - 1; i <= y_grid + 1; i++)
                for (int j = x_grid - 1; j <= x_grid + 1; j++) {

                    // Look through bin and add fine repulsions
                    for (BI = GET_BIN(k, i, j).begin(); BI < GET_BIN(k, i, j).end(); ++BI) {
                        x_dist =  Nx - (BI->x);
                        y_dist =  Ny - (BI->y);
                        z_dist =  Nz - (BI->z);
                        distance = x_dist * x_dist + y_dist * y_dist + z_dist * z_dist;
                        density += 1e-4 / (distance + 1e-50);
                    }
                }

        // Course density
    } else {

        // Add rough estimate
        density = Density[z_grid][y_grid][x_grid];
        density *= density;
    }

    return density;
}

/// Wrapper functions for the Add and subtract methods
/// Nodes should all be passed by constant ref

void DensityGrid::Add(Node &n, bool fineDensity) {
    if (fineDensity) {
        fineAdd(n);
    } else {
        Add(n);
    }
}

void DensityGrid::Subtract( Node &n, bool first_add,
                            bool fine_first_add, bool fineDensity) {
    if ( fineDensity && !fine_first_add ) {
        fineSubtract (n);
    } else if ( !first_add ) {
        Subtract(n);
    }
}


/***************************************************
 * Function: DensityGrid::Subtract                *
 * Description: Subtract a node from density grid  *
 **************************************************/
void DensityGrid::Subtract(Node &N) {
    int x_grid, y_grid, z_grid, diam;
    float *den_ptr, *fall_ptr;

    /* Where to subtract */
    x_grid = (int)((N.sub_x + HALF_VIEW + .5) * VIEW_TO_GRID);
    y_grid = (int)((N.sub_y + HALF_VIEW + .5) * VIEW_TO_GRID);
    z_grid = (int)((N.sub_z + HALF_VIEW + .5) * VIEW_TO_GRID);
    x_grid -= RADIUS;
    y_grid -= RADIUS;
    z_grid -= RADIUS;
    diam = 2 * RADIUS;

    // check to see that we are inside grid
    if ( (x_grid >= GRID_SIZE) || (x_grid < 0) ||
         (y_grid >= GRID_SIZE) || (y_grid < 0) ||
         (z_grid >= GRID_SIZE) || (z_grid < 0) ) {
#ifdef MUSE_MPI
        MPI_Abort ( MPI_COMM_WORLD, 1 );
#else
        igraph_error("Exceeded density grid in DrL", IGRAPH_FILE_BASENAME,
                     __LINE__, IGRAPH_EDRL);
        return;
#endif
    }

    /* Subtract density values */
    den_ptr = &Density[z_grid][y_grid][x_grid];
    fall_ptr = &fall_off[0][0][0];
    for (int i = 0; i <= diam; i++) {
        for (int j = 0; j <= diam; j++)
            for (int k = 0; k <= diam; k++) {
                *den_ptr++ -= *fall_ptr++;
            }
        den_ptr += GRID_SIZE - (diam + 1);
    }
}

/***************************************************
 * Function: DensityGrid::Add                     *
 * Description: Add a node to the density grid     *
 **************************************************/
void DensityGrid::Add(Node &N) {

    int x_grid, y_grid, z_grid, diam;
    float *den_ptr, *fall_ptr;


    /* Where to add */
    x_grid = (int)((N.x + HALF_VIEW + .5) * VIEW_TO_GRID);
    y_grid = (int)((N.y + HALF_VIEW + .5) * VIEW_TO_GRID);
    z_grid = (int)((N.z + HALF_VIEW + .5) * VIEW_TO_GRID);

    N.sub_x = N.x;
    N.sub_y = N.y;
    N.sub_z = N.z;

    x_grid -= RADIUS;
    y_grid -= RADIUS;
    z_grid -= RADIUS;
    diam = 2 * RADIUS;

    // check to see that we are inside grid
    if ( (x_grid >= GRID_SIZE) || (x_grid < 0) ||
         (y_grid >= GRID_SIZE) || (y_grid < 0) ||
         (z_grid >= GRID_SIZE) || (z_grid < 0) ) {
#ifdef MUSE_MPI
        MPI_Abort ( MPI_COMM_WORLD, 1 );
#else
        igraph_error("Exceeded density grid in DrL", IGRAPH_FILE_BASENAME,
                     __LINE__, IGRAPH_EDRL);
        return;
#endif
    }

    /* Add density values */
    den_ptr = &Density[z_grid][y_grid][x_grid];
    fall_ptr = &fall_off[0][0][0];
    for (int i = 0; i <= diam; i++) {
        for (int j = 0; j <= diam; j++)
            for (int k = 0; k <= diam; k++) {
                *den_ptr++ += *fall_ptr++;
            }
        den_ptr += GRID_SIZE - (diam + 1);
    }

}

/***************************************************
 * Function: DensityGrid::fineSubtract             *
 * Description: Subtract a node from bins          *
 **************************************************/
void DensityGrid::fineSubtract(Node &N) {
    int x_grid, y_grid, z_grid;

    /* Where to subtract */
    x_grid = (int)((N.sub_x + HALF_VIEW + .5) * VIEW_TO_GRID);
    y_grid = (int)((N.sub_y + HALF_VIEW + .5) * VIEW_TO_GRID);
    z_grid = (int)((N.sub_z + HALF_VIEW + .5) * VIEW_TO_GRID);
    GET_BIN(z_grid, y_grid, x_grid).pop_front();
}

/***************************************************
 * Function: DensityGrid::fineAdd                  *
 * Description: Add a node to the bins             *
 **************************************************/
void DensityGrid::fineAdd(Node &N) {
    int x_grid, y_grid, z_grid;

    /* Where to add */
    x_grid = (int)((N.x + HALF_VIEW + .5) * VIEW_TO_GRID);
    y_grid = (int)((N.y + HALF_VIEW + .5) * VIEW_TO_GRID);
    z_grid = (int)((N.z + HALF_VIEW + .5) * VIEW_TO_GRID);
    N.sub_x = N.x;
    N.sub_y = N.y;
    N.sub_z = N.z;
    GET_BIN(z_grid, y_grid, x_grid).push_back(N);
}

} // namespace drl3d
