/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2006-2012  Gabor Csardi <csardi.gabor@gmail.com>
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

/* The original version of this file was written by Jörg Reichardt
   The original copyright notice follows here */

/***************************************************************************
                          NetRoutines.h  -  description
                             -------------------
    begin                : Tue Oct 28 2003
    copyright            : (C) 2003 by Joerg Reichardt
    email                : reichardt@mitte
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef NETROUTINES_H
#define NETROUTINES_H

#include "NetDataTypes.h"
#include "igraph_types.h"
#include "igraph_datatype.h"

int igraph_i_read_network(const igraph_t *graph,
                          const igraph_vector_t *weights,
                          network *net, igraph_bool_t use_weights,
                          unsigned int states);

void reduce_cliques(DLList<ClusterList<NNode*>*>*, FILE *file);
void reduce_cliques2(network*, bool,  long );
void clear_all_markers(network *net);

#endif

