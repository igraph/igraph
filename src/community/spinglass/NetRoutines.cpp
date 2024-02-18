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
                          NetRoutines.cpp  -  description
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

#include "NetRoutines.h"
#include "NetDataTypes.h"

#include "igraph_types.h"
#include "igraph_interface.h"

igraph_error_t igraph_i_read_network_spinglass(
    const igraph_t *graph, const igraph_vector_t *weights,
    network *net, igraph_bool_t use_weights) {

    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_integer_t no_of_edges = igraph_ecount(graph);
    double sum_weight;

    for (igraph_integer_t vid = 0; vid < no_of_nodes; vid++) {
        char name[SPINGLASS_MAX_NAME_LEN];
        snprintf(name, sizeof(name) / sizeof(name[0]), "%" IGRAPH_PRId "", vid+1);
        net->node_list.Push(new NNode(vid, 0, &net->link_list, name));
    }

    sum_weight = 0.0;
    for (igraph_integer_t eid = 0; eid < no_of_edges; eid++) {
        igraph_integer_t v1 = IGRAPH_FROM(graph, eid);
        igraph_integer_t v2 = IGRAPH_TO(graph, eid);
        igraph_real_t w = use_weights ? VECTOR(*weights)[eid] : 1.0;

        NNode *node1 = net->node_list.Get(v1);
        NNode *node2 = net->node_list.Get(v2);

        node1->Connect_To(node2, w);

        sum_weight += w;
    }

    net->sum_weights = sum_weight;

    return IGRAPH_SUCCESS;
}
