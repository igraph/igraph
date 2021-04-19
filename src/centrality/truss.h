/*

  Copyright 2017 The Johns Hopkins University Applied Physics Laboratory LLC. All Rights Reserved.

  Truss algorithm for cohesive subgroups.

  Author: Alex Perrone
  Date: 2017-08-03

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

#ifndef TRUSS_H
#define TRUSS_H

#include <igraph.h>

void unpack(const igraph_vector_int_t *tri, igraph_vector_t *unpacked_tri);
void compute_support(const igraph_vector_t *eid, igraph_vector_int_t *support);
void trussness(const igraph_t *graph, igraph_vector_int_t *support,
  igraph_vector_int_t *truss);
int igraph_truss(const igraph_t* graph, igraph_vector_int_t* truss);

#endif
