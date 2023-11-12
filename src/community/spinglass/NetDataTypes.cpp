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
                          NetDataTypes.cpp  -  description
                             -------------------
    begin                : Mon Oct 6 2003
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

#include "NetDataTypes.h"

int NNode::Connect_To(NNode* neighbour, double weight_) {
    NLink *link;
    //sollen doppelte Links erlaubt sein??  NEIN
    if (!neighbour) {
        return 0;
    }
    if (!(neighbours.Is_In_List(neighbour)) && (neighbour != this)) {
        neighbours.Push(neighbour);        // nachbar hier eintragen
        neighbour->neighbours.Push(this); // diesen knoten beim nachbarn eintragen

        link = new NLink(this, neighbour, weight_);     //link erzeugen
        global_link_list->Push(link);                        // in globaler liste eintragen
        n_links.Push(link);                                   // bei diesem Knoten eintragen
        neighbour->n_links.Push(link);                  // beim nachbarn eintragen

        return 1;
    }
    return 0;
}

NLink *NNode::Get_LinkToNeighbour(const NNode* neighbour) {
    DLList_Iter<NLink*> iter;
    NLink *l_cur, *link = nullptr;
    bool found = false;
    // finde einen bestimmten Link aus der Liste der links eines Knotens
    l_cur = iter.First(&n_links);
    while (!iter.End() && !found) {
        if (((l_cur->Get_Start() == this) && (l_cur->Get_End() == neighbour)) || ((l_cur->Get_End() == this) && (l_cur->Get_Start() == neighbour))) {
            found = true;
            link = l_cur;
        }
        l_cur = iter.Next();
    }
    if (found) {
        return link;
    } else {
        return nullptr;
    }
}

igraph_integer_t NNode::Disconnect_From(NNode* neighbour) {
    //sollen doppelte Links erlaubt sein??  s.o.
    neighbours.fDelete(neighbour);
    n_links.fDelete(Get_LinkToNeighbour(neighbour));
    neighbour->n_links.fDelete(neighbour->Get_LinkToNeighbour(this));
    neighbour->neighbours.fDelete(this);
    return 1;
}

igraph_integer_t NNode::Disconnect_From_All() {
    igraph_integer_t number_of_neighbours = 0;
    while (neighbours.Size()) {
        Disconnect_From(neighbours.Pop());
        number_of_neighbours++;
    }
    return number_of_neighbours ;
}
