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
#include "igraph_conversion.h"

#include <climits>

igraph_error_t igraph_i_read_network(
    const igraph_t *graph, const igraph_vector_t *weights,
    network *net, igraph_bool_t use_weights) {

    double sum_weight = 0.0, min_weight = IGRAPH_POSINFINITY, max_weight = IGRAPH_NEGINFINITY;
    unsigned long min_k = ULONG_MAX, max_k = 0;
    char name[255];
    NNode *node1, *node2;
    DLList_Iter<NNode*> iter;
    igraph_vector_int_t edgelist;
    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_integer_t no_of_edges = igraph_ecount(graph);
    const char *empty = "";

    IGRAPH_VECTOR_INT_INIT_FINALLY(&edgelist, no_of_edges * 2);
    IGRAPH_CHECK(igraph_get_edgelist(graph, &edgelist, false /* rowwise */));

    for (igraph_integer_t ii = 0; ii < no_of_nodes; ii++) {
        net->node_list.Push(new NNode(ii, 0, &net->link_list, empty));
    }

    for (igraph_integer_t ii = 0; ii < no_of_edges; ii++) {
        igraph_integer_t i1 = VECTOR(edgelist)[2 * ii];
        igraph_integer_t i2 = VECTOR(edgelist)[2 * ii + 1];
        igraph_real_t Links = use_weights ? VECTOR(*weights)[ii] : 1.0;

        node1 = net->node_list.Get(i1);
        snprintf(name, sizeof(name) / sizeof(name[0]), "%" IGRAPH_PRId "", i1+1);
        node1->Set_Name(name);

        node2 = net->node_list.Get(i2);
        snprintf(name, sizeof(name) / sizeof(name[0]), "%" IGRAPH_PRId "", i2+1);
        node2->Set_Name(name);

        node1->Connect_To(node2, Links);

        if (Links < min_weight) {
            min_weight = Links;
        }
        if (Links > max_weight) {
            max_weight = Links;
        }
        sum_weight += Links;
    }

    IGRAPH_FINALLY_CLEAN(1);
    igraph_vector_int_destroy(&edgelist);

    node1 = iter.First(&net->node_list);
    while (!iter.End()) {
        if (node1->Get_Degree() > max_k) {
            max_k = node1->Get_Degree();
        }
        if (node1->Get_Degree() < min_k) {
            min_k = node1->Get_Degree();
        }
        node1 = iter.Next();
    }
    net->sum_weights = sum_weight;

    return IGRAPH_SUCCESS;
}
