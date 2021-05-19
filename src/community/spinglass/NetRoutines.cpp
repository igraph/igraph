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

int igraph_i_read_network(const igraph_t *graph,
                          const igraph_vector_t *weights,
                          network *net, igraph_bool_t use_weights,
                          unsigned int states) {

    double av_k = 0.0, sum_weight = 0.0, min_weight = 1e60, max_weight = -1e60;
    unsigned long min_k = 999999999, max_k = 0;
    char name[255];
    NNode *node1, *node2;
    DLList_Iter<NNode*> iter;
    igraph_vector_t edgelist;
    long int no_of_nodes = (long int) igraph_vcount(graph);
    long int no_of_edges = (long int) igraph_ecount(graph);
    long int ii;
    const char *empty = "";

    IGRAPH_VECTOR_INIT_FINALLY(&edgelist, no_of_edges * 2);
    IGRAPH_CHECK(igraph_get_edgelist(graph, &edgelist, 0 /* rowwise */));

    for (ii = 0; ii < no_of_nodes; ii++) {
        net->node_list->Push(new NNode(ii, 0, net->link_list, empty, states));
    }

    for (ii = 0; ii < no_of_edges; ii++) {
        long int i1 = (long int)VECTOR(edgelist)[2 * ii];
        long int i2 = (long int)VECTOR(edgelist)[2 * ii + 1];
        igraph_real_t Links;
        if (use_weights) {
            Links = VECTOR(*weights)[ii];
        } else {
            Links = 1.0;
        }

        node1 = net->node_list->Get(i1);
        sprintf(name, "%li", i1+1);
        node1->Set_Name(name);

        node2 = net->node_list->Get(i2);
        sprintf(name, "%li", i2+1);
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
    igraph_vector_destroy(&edgelist);

    node1 = iter.First(net->node_list);
    while (!iter.End()) {
        if (node1->Get_Degree() > max_k) {
            max_k = node1->Get_Degree();
        }
        if (node1->Get_Degree() < min_k) {
            min_k = node1->Get_Degree();
        }
        av_k += node1->Get_Degree();
        node1 = iter.Next();
    }
    net->av_k = av_k / double(net->node_list->Size());
    net->sum_weights = sum_weight;
    net->av_weight = sum_weight / double(net->link_list->Size());
    net->min_k = min_k;
    net->max_k = max_k;
    net->min_weight = min_weight;
    net->max_weight = max_weight;
    net->sum_bids = 0;
    net->min_bids = 0;
    net->max_bids = 0;

    return IGRAPH_SUCCESS;
}

//###############################################################################################################
void reduce_cliques(DLList<ClusterList<NNode*>*> *global_cluster_list, FILE *file) {
    unsigned long size;
    ClusterList<NNode*> *c_cur, *largest_c = NULL;
    DLList<ClusterList<NNode*>*> *subsets;
    DLList_Iter<ClusterList<NNode*>*> c_iter, sub_iter;
    DLList_Iter<NNode*> iter;
    NNode *n_cur;

    if (!(global_cluster_list->Size())) {
        return;
    }
    //wir suchen den groessten Cluster

    c_cur = c_iter.First(global_cluster_list);
    size = 0;
    while (!(c_iter.End())) {
        if (c_cur->Size() > size) {
            size = c_cur->Size();
            largest_c = c_cur;
        }
        c_cur = c_iter.Next();
    }
// printf("Groesster Cluster hat %u Elemente.\n",largest_c->Size());

    //Schauen, ob es Teilmengen gibt, die ebenfalls gefunden wurden
    subsets = new DLList<ClusterList<NNode*>*>();
    c_cur = c_iter.First(global_cluster_list);
    while (!(c_iter.End())) {
        if ((*c_cur < *largest_c || *c_cur == *largest_c) && c_cur != largest_c) { //alle echten Teilcluster von largest_c und die doppelten
            subsets->Push(c_cur);
        }
        c_cur = c_iter.Next();
    }
    // die gefundenen Subsets werden aus der cluster_liste geloescht
    while (subsets->Size()) {
        global_cluster_list->fDelete(subsets->Pop());
    }
    delete subsets;
    // Dann schreiben wir den groessten Cluster in das File
    fprintf(file, "Energie: %1.12f   Nodes:%3lu    -   ", largest_c->Get_Energy(), largest_c->Size());

    n_cur = iter.First(largest_c);
    while (!(iter.End())) {
        fprintf(file, "%s", n_cur->Get_Name());
        n_cur = iter.Next();
        if (n_cur) {
            fprintf(file, ", ");
        }
    }
    fprintf(file, "\n");


    //Schliesslich schmeissen wir noch den eben gefundenen groessten Cluster raus
    global_cluster_list->fDelete(largest_c);
    //und dann geht es von vorn mit der Reduzierten ClusterListe los
    reduce_cliques(global_cluster_list, file);

}
//##################################################################################
void reduce_cliques2(network *net, bool only_double, long marker) {
    unsigned long size;
    ClusterList<NNode*> *c_cur, *largest_c = NULL;
    DLList_Iter<ClusterList<NNode*>*> c_iter;
    do {
        //wir suchen den groessten, nicht markierten Cluster
        size = 0;
        c_cur = c_iter.First(net->cluster_list);
        while (!(c_iter.End())) {
            if ((c_cur->Size() > size) && (c_cur->Get_Marker() != marker)) {
                size = c_cur->Size();
                largest_c = c_cur;
            }
            c_cur = c_iter.Next();
        }
        // printf("Groesster Cluster hat %u Elemente.\n",largest_c->Size());
        //Schauen, ob es Teilmengen gibt, die ebenfalls gefunden wurden
        c_cur = c_iter.First(net->cluster_list);
        while (!(c_iter.End())) {
            if (((!only_double && (*c_cur < *largest_c)) || (*c_cur == *largest_c)) && (c_cur != largest_c)) { //alle echten Teilcluster von largest_c und die doppelten
                net->cluster_list->fDelete(c_cur);
                while (c_cur->Get_Candidates()->Size()) {
                    c_cur->Get_Candidates()->Pop();
                }
                while (c_cur->Size()) {
                    c_cur->Pop();    // die knoten aber nicht loeschen!!
                }
                delete c_cur;    // nicht vergessen, die global geloeschte Clusterliste zu loeschen
            }
            c_cur = c_iter.Next();
        }
        //Schliesslich markieren wir noch den eben gefundenen groessten Cluster
        largest_c->Set_Marker(marker);
    } while (size);
}

//##################################################################################################
unsigned long iterate_nsf_hierarchy(NNode *parent, unsigned long depth, FILE *file) {
    NNode* next_node;
    unsigned long newdepth, maxdepth;
    bool first = true;
    DLList_Iter<NNode*> *iter;
    maxdepth = newdepth = depth;
    iter = new DLList_Iter<NNode*>;
    next_node = iter->First(parent->Get_Neighbours());
    while (!(iter->End())) {
        if (next_node->Get_Marker() > parent->Get_Marker()) { // wir gehen nach unten
            if (first) {
                fprintf(file, ",(");    // eine Neue Klammer auf
            }
            if (first) {
                fprintf(file, "%s", next_node->Get_Name());    // nur vor dem ersten kein Komma
            } else {
                fprintf(file, ",%s", next_node->Get_Name());    // sonst immer mit Komma
            }
            first = false;
            newdepth = iterate_nsf_hierarchy(next_node, depth + 1, file);
            if (maxdepth < newdepth) {
                maxdepth = newdepth;
            }
        }
        next_node = iter->Next();
    }
    if (!first) {
        fprintf(file, ")");    //hat es ueberhaupt einen gegeben?
    }
    //dann klamer zu!
    delete iter;
    return maxdepth;
}

//################################################################
void clear_all_markers(network *net) {
    DLList_Iter<NNode*> iter;
    NNode *n_cur;
    n_cur = iter.First(net->node_list);
    while (!iter.End()) {
        n_cur->Set_Marker(0);
        n_cur = iter.Next();
    }
}
