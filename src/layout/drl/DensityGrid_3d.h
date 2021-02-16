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
#ifndef __DENSITY_GRID_H__
#define __DENSITY_GRID_H__


// Compile time adjustable parameters

#include "drl_layout_3d.h"
#include "drl_Node_3d.h"
#ifdef MUSE_MPI
    #include <mpi.h>
#endif

#include <deque>

namespace drl3d {

class DensityGrid {

public:

    // Methods
    void Init();
    void Subtract(Node &n, bool first_add, bool fine_first_add, bool fineDensity);
    void Add(Node &n, bool fineDensity );
    float GetDensity(float Nx, float Ny, float Nz, bool fineDensity);

    // Contructor/Destructor
    DensityGrid() {};
    ~DensityGrid();

private:

    // Private Members
    void Subtract( Node &N );
    void Add( Node &N );
    void fineSubtract( Node &N );
    void fineAdd( Node &N );

    // new dynamic variables -- SBM
    float (*fall_off)[RADIUS * 2 + 1][RADIUS * 2 + 1];
    float (*Density)[GRID_SIZE][GRID_SIZE];
    std::deque<Node>* Bins;

    // old static variables
    //float fall_off[RADIUS*2+1][RADIUS*2+1];
    //float Density[GRID_SIZE][GRID_SIZE];
    //deque<Node *> Bins[GRID_SIZE][GRID_SIZE];
};

} // namespace drl3d

#endif // __DENSITY_GRID_H__

