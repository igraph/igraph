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
// This file contains compile time parameters which affect the entire
// DrL program.

#define DRL_VERSION "3.2 5/5/2006"

// compile time parameters for MPI message passing
#define MAX_PROCS 256      // maximum number of processors
#define MAX_FILE_NAME 250   // max length of filename
#define MAX_INT_LENGTH 4   // max length of integer suffix of intermediate .coord file

// Compile time adjustable parameters for the Density grid

#define GRID_SIZE 100           // size of Density grid
#define VIEW_SIZE 250.0     // actual physical size of layout plane
// these values use more memory but have
// little effect on performance or layout

#define RADIUS 10               // radius for density fall-off:
// larger values tends to slow down
// the program and clump the data

#define HALF_VIEW 125.0         // 1/2 of VIEW_SIZE
#define VIEW_TO_GRID .4         // ratio of GRID_SIZE to VIEW_SIZE

/*
// original values for VxOrd
#define GRID_SIZE 400           // size of VxOrd Density grid
#define VIEW_SIZE 1600.0        // actual physical size of VxOrd plane
#define RADIUS 10               // radius for density fall-off

#define HALF_VIEW 800           // 1/2 of VIEW_SIZE
#define VIEW_TO_GRID .25        // ratio of GRID_SIZE to VIEW_SIZE
*/
