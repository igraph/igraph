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
#ifdef HAVE_CONFIG_H
    #include <config.h>
#endif

#include "NetDataTypes.h"
#include <cstring>

//#################################################################################
//###############################################################################
//Constructor
NNode::NNode(unsigned long ind, unsigned long c_ind, DLList<NLink*> *ll, const char *n, int states) {
    index = ind;
    cluster_index = c_ind;
    neighbours = new DLList<NNode*>();
    n_links = new DLList<NLink*>();
    global_link_list = ll;
    strcpy(name, n);
    color.red = 0;
    color.green = 0;
    color.blue = 0;
    strcpy(color.pajek_c, "Green");
    clustering = 0.0;
    marker = 0;
    affiliations = 0;
    weight = 0.0;
    affinity = 0.0;
    distance = 0;
    max_states = states;
    state_history = new unsigned long[states + 1];
}

//Destructor
NNode::~NNode() {
    Disconnect_From_All();
    delete neighbours;
    delete n_links;
    delete [] state_history;
    neighbours = NULL;
    n_links = NULL;
    state_history = NULL;
}

void NNode::Add_StateHistory(unsigned int state) {
    if (max_states >= state) {
        state_history[state]++;
    }
}

void NNode::Set_Color(RGBcolor c) {
    color.red = c.red; color.blue = c.blue; color.green = c.green;
    strcpy(color.pajek_c, c.pajek_c);
}

int NNode::Connect_To(NNode* neighbour, double weight_) {
    NLink *link;
    //sollen doppelte Links erlaubt sein??  NEIN
    if (!neighbour) {
        return 0;
    }
    if (!(neighbours->Is_In_List(neighbour)) && (neighbour != this)) {
        neighbours->Push(neighbour);        // nachbar hier eintragen
        neighbour->neighbours->Push(this); // diesen knoten beim nachbarn eintragen

        link = new NLink(this, neighbour, weight_);     //link erzeugen
        global_link_list->Push(link);                        // in globaler liste eintragen
        n_links->Push(link);                                   // bei diesem Knoten eintragen
        neighbour->n_links->Push(link);                  // beim nachbarn eintragen

        return (1);
    }
    return (0);
}

NLink *NNode::Get_LinkToNeighbour(NNode* neighbour) {
    DLList_Iter<NLink*> iter;
    NLink *l_cur, *link = NULL;
    bool found = false;
    // finde einen bestimmten Link aus der Liste der links eines Knotens
    l_cur = iter.First(n_links);
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
        return NULL;
    }
}

int NNode::Disconnect_From(NNode* neighbour) {
    //sollen doppelte Links erlaubt sein??  s.o.
    if (!neighbours) {
        return 0;
    }
    neighbours->fDelete(neighbour);
    n_links->fDelete(Get_LinkToNeighbour(neighbour));
    neighbour->n_links->fDelete(neighbour->Get_LinkToNeighbour(this));
    neighbour->neighbours->fDelete(this);
    return 1;
}

int NNode::Disconnect_From_All() {
    int number_of_neighbours = 0;
    while (neighbours->Size()) {
        Disconnect_From(neighbours->Pop());
        number_of_neighbours++;
    }
    return (number_of_neighbours) ;
}

/*
int NNode::Disconnect_From_All_Grandchildren()
{
 int n_l=links->Size();
 unsigned long pos=0;
 while ((n_l--)>1) {  //alle bis auf das erste loeschen
      pos=(links->Get(n_l+1))->links->Is_In_List(this);
     // printf("%d %d\n",n_l,pos);
      (links->Get(n_l+1))->links->Delete(pos);
  }
 return(pos) ;
}
*/

double NNode::Get_Links_Among_Neigbours() {
//  long neighbours1, neighbours2;
    double lam = 0;
    DLList_Iter<NNode*> iter1, iter2;
//  neighbours1=neighbours->Size();        //so viele Nachbarn hat die Betrachtete Node
    NNode *step1, *step2;
    step1 = iter1.First(neighbours);
    while (!iter1.End()) { //  for (int n1=1;n1<=neighbours1; n1++)
        //step1=neighbours->Get(n1);
        //neighbours2=step1->neighbours->Size();  //so viele Nachbarn hat der n1-ste Nachbar
        step2 = iter2.First(step1->Get_Neighbours());
        while (!iter2.End()) { //for (int n2=1;n2<=neighbours2; n2++)
            //step2=step1->neighbours->Get(n2);
            if (step2->Get_Neighbours()->Is_In_List(this)) {
                lam++;
            }
            step2 = iter2.Next();
        }
        step1 = iter1.Next();
    }
    return (lam / 2.0);
}


double NNode::Get_Clustering() {
    double c;
    unsigned long k;
    k = neighbours->Size();
    if (k <= 1) {
        return (0);
    }
    c = 2.0 * Get_Links_Among_Neigbours() / double(k * k - k);
    return (c);
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++

//Constructor
NLink::NLink(NNode *s, NNode *e, double w) {
    start = s;
    end = e;
    weight = w;
    old_weight = 0;
    marker = 0;
}

//Destructor
NLink::~NLink() {
    if (start && end) {
        start->Disconnect_From(end);
    }
}
