/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2007-2012  Gabor Csardi <csardi.gabor@gmail.com>
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

/* The original version of this file was written by Pascal Pons
   The original copyright notice follows here. The FSF address was
   fixed by Tamas Nepusz */

// File: communities.cpp
//-----------------------------------------------------------------------------
// Walktrap v0.2 -- Finds community structure of networks using random walks
// Copyright (C) 2004-2005 Pascal Pons
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
// 02110-1301 USA
//-----------------------------------------------------------------------------
// Author   : Pascal Pons
// Email    : pascal.pons@gmail.com
// Web page : http://www-rp.lip6.fr/~latapy/PP/walktrap.html
// Location : Paris, France
// Time     : June 2005
//-----------------------------------------------------------------------------
// see readme.txt for more details

#include "walktrap_communities.h"
#include "config.h"
#include <algorithm>
#include <cmath>

using namespace std;

namespace igraph {

namespace walktrap {

IGRAPH_THREAD_LOCAL int Probabilities::length = 0;
IGRAPH_THREAD_LOCAL Communities* Probabilities::C = 0;
IGRAPH_THREAD_LOCAL float* Probabilities::tmp_vector1 = 0;
IGRAPH_THREAD_LOCAL float* Probabilities::tmp_vector2 = 0;
IGRAPH_THREAD_LOCAL int* Probabilities::id = 0;
IGRAPH_THREAD_LOCAL int* Probabilities::vertices1 = 0;
IGRAPH_THREAD_LOCAL int* Probabilities::vertices2 = 0;
IGRAPH_THREAD_LOCAL int Probabilities::current_id = 0;


Neighbor::Neighbor() {
    next_community1 = 0;
    previous_community1 = 0;
    next_community2 = 0;
    previous_community2 = 0;
    heap_index = -1;
}

Probabilities::~Probabilities() {
    C->memory_used -= memory();
    if (P) {
        delete[] P;
    }
    if (vertices) {
        delete[] vertices;
    }
}

Probabilities::Probabilities(int community) {
    Graph* G = C->G;
    int nb_vertices1 = 0;
    int nb_vertices2 = 0;

    float initial_proba = 1. / float(C->communities[community].size);
    int last =  C->members[C->communities[community].last_member];
    for (int m = C->communities[community].first_member; m != last; m = C->members[m]) {
        tmp_vector1[m] = initial_proba;
        vertices1[nb_vertices1++] = m;
    }

    for (int t = 0; t < length; t++) {
        current_id++;
        if (nb_vertices1 > (G->nb_vertices / 2)) {
            nb_vertices2 = G->nb_vertices;
            for (int i = 0; i < G->nb_vertices; i++) {
                tmp_vector2[i] = 0.;
            }
            if (nb_vertices1 == G->nb_vertices) {
                for (int i = 0; i < G->nb_vertices; i++) {
                    float proba = tmp_vector1[i] / G->vertices[i].total_weight;
                    for (int j = 0; j < G->vertices[i].degree; j++) {
                        tmp_vector2[G->vertices[i].edges[j].neighbor] += proba * G->vertices[i].edges[j].weight;
                    }
                }
            } else {
                for (int i = 0; i < nb_vertices1; i++) {
                    int v1 = vertices1[i];
                    float proba = tmp_vector1[v1] / G->vertices[v1].total_weight;
                    for (int j = 0; j < G->vertices[v1].degree; j++) {
                        tmp_vector2[G->vertices[v1].edges[j].neighbor] += proba * G->vertices[v1].edges[j].weight;
                    }
                }
            }
        } else {
            nb_vertices2 = 0;
            for (int i = 0; i < nb_vertices1; i++) {
                int v1 = vertices1[i];
                float proba = tmp_vector1[v1] / G->vertices[v1].total_weight;
                for (int j = 0; j < G->vertices[v1].degree; j++) {
                    int v2 = G->vertices[v1].edges[j].neighbor;
                    if (id[v2] == current_id) {
                        tmp_vector2[v2] += proba * G->vertices[v1].edges[j].weight;
                    } else {
                        tmp_vector2[v2] = proba * G->vertices[v1].edges[j].weight;
                        id[v2] = current_id;
                        vertices2[nb_vertices2++] = v2;
                    }
                }
            }
        }
        float* tmp = tmp_vector2;
        tmp_vector2 = tmp_vector1;
        tmp_vector1 = tmp;

        int* tmp2 = vertices2;
        vertices2 = vertices1;
        vertices1 = tmp2;

        nb_vertices1 = nb_vertices2;
    }

    if (nb_vertices1 > (G->nb_vertices / 2)) {
        P = new float[G->nb_vertices];
        size = G->nb_vertices;
        vertices = 0;
        if (nb_vertices1 == G->nb_vertices) {
            for (int i = 0; i < G->nb_vertices; i++) {
                P[i] = tmp_vector1[i] / sqrt(G->vertices[i].total_weight);
            }
        } else {
            for (int i = 0; i < G->nb_vertices; i++) {
                P[i] = 0.;
            }
            for (int i = 0; i < nb_vertices1; i++) {
                P[vertices1[i]] = tmp_vector1[vertices1[i]] / sqrt(G->vertices[vertices1[i]].total_weight);
            }
        }
    } else {
        P = new float[nb_vertices1];
        size = nb_vertices1;
        vertices = new int[nb_vertices1];
        int j = 0;
        for (int i = 0; i < G->nb_vertices; i++) {
            if (id[i] == current_id) {
                P[j] = tmp_vector1[i] / sqrt(G->vertices[i].total_weight);
                vertices[j] = i;
                j++;
            }
        }
    }
    C->memory_used += memory();
}

Probabilities::Probabilities(int community1, int community2) {
    // The two following probability vectors must exist.
    // Do not call this function if it is not the case.
    Probabilities* P1 = C->communities[community1].P;
    Probabilities* P2 = C->communities[community2].P;

    float w1 = float(C->communities[community1].size) / float(C->communities[community1].size + C->communities[community2].size);
    float w2 = float(C->communities[community2].size) / float(C->communities[community1].size + C->communities[community2].size);


    if (P1->size == C->G->nb_vertices) {
        P = new float[C->G->nb_vertices];
        size = C->G->nb_vertices;
        vertices = 0;

        if (P2->size == C->G->nb_vertices) { // two full vectors
            for (int i = 0; i < C->G->nb_vertices; i++) {
                P[i] = P1->P[i] * w1 + P2->P[i] * w2;
            }
        } else { // P1 full vector, P2 partial vector
            int j = 0;
            for (int i = 0; i < P2->size; i++) {
                for (; j < P2->vertices[i]; j++) {
                    P[j] = P1->P[j] * w1;
                }
                P[j] = P1->P[j] * w1 + P2->P[i] * w2;
                j++;
            }
            for (; j < C->G->nb_vertices; j++) {
                P[j] = P1->P[j] * w1;
            }
        }
    } else {
        if (P2->size == C->G->nb_vertices) { // P1 partial vector, P2 full vector
            P = new float[C->G->nb_vertices];
            size = C->G->nb_vertices;
            vertices = 0;

            int j = 0;
            for (int i = 0; i < P1->size; i++) {
                for (; j < P1->vertices[i]; j++) {
                    P[j] = P2->P[j] * w2;
                }
                P[j] = P1->P[i] * w1 + P2->P[j] * w2;
                j++;
            }
            for (; j < C->G->nb_vertices; j++) {
                P[j] = P2->P[j] * w2;
            }
        } else { // two partial vectors
            int i = 0;
            int j = 0;
            int nb_vertices1 = 0;
            while ((i < P1->size) && (j < P2->size)) {
                if (P1->vertices[i] < P2->vertices[j]) {
                    tmp_vector1[P1->vertices[i]] = P1->P[i] * w1;
                    vertices1[nb_vertices1++] = P1->vertices[i];
                    i++;
                    continue;
                }
                if (P1->vertices[i] > P2->vertices[j]) {
                    tmp_vector1[P2->vertices[j]] = P2->P[j] * w2;
                    vertices1[nb_vertices1++] = P2->vertices[j];
                    j++;
                    continue;
                }
                tmp_vector1[P1->vertices[i]] = P1->P[i] * w1 + P2->P[j] * w2;
                vertices1[nb_vertices1++] = P1->vertices[i];
                i++;
                j++;
            }
            if (i == P1->size) {
                for (; j < P2->size; j++) {
                    tmp_vector1[P2->vertices[j]] = P2->P[j] * w2;
                    vertices1[nb_vertices1++] = P2->vertices[j];
                }
            } else {
                for (; i < P1->size; i++) {
                    tmp_vector1[P1->vertices[i]] = P1->P[i] * w1;
                    vertices1[nb_vertices1++] = P1->vertices[i];
                }
            }

            if (nb_vertices1 > (C->G->nb_vertices / 2)) {
                P = new float[C->G->nb_vertices];
                size = C->G->nb_vertices;
                vertices = 0;
                for (int i = 0; i < C->G->nb_vertices; i++) {
                    P[i] = 0.;
                }
                for (int i = 0; i < nb_vertices1; i++) {
                    P[vertices1[i]] = tmp_vector1[vertices1[i]];
                }
            } else {
                P = new float[nb_vertices1];
                size = nb_vertices1;
                vertices = new int[nb_vertices1];
                for (int i = 0; i < nb_vertices1; i++) {
                    vertices[i] = vertices1[i];
                    P[i] = tmp_vector1[vertices1[i]];
                }
            }
        }
    }

    C->memory_used += memory();
}

double Probabilities::compute_distance(const Probabilities* P2) const {
    double r = 0.;
    if (vertices) {
        if (P2->vertices) { // two partial vectors
            int i = 0;
            int j = 0;
            while ((i < size) && (j < P2->size)) {
                if (vertices[i] < P2->vertices[j]) {
                    r += P[i] * P[i];
                    i++;
                    continue;
                }
                if (vertices[i] > P2->vertices[j]) {
                    r += P2->P[j] * P2->P[j];
                    j++;
                    continue;
                }
                r += (P[i] - P2->P[j]) * (P[i] - P2->P[j]);
                i++;
                j++;
            }
            if (i == size) {
                for (; j < P2->size; j++) {
                    r += P2->P[j] * P2->P[j];
                }
            } else {
                for (; i < size; i++) {
                    r += P[i] * P[i];
                }
            }
        } else { // P1 partial vector, P2 full vector

            int i = 0;
            for (int j = 0; j < size; j++) {
                for (; i < vertices[j]; i++) {
                    r += P2->P[i] * P2->P[i];
                }
                r += (P[j] - P2->P[i]) * (P[j] - P2->P[i]);
                i++;
            }
            for (; i < P2->size; i++) {
                r += P2->P[i] * P2->P[i];
            }
        }
    } else {
        if (P2->vertices) { // P1 full vector, P2 partial vector
            int i = 0;
            for (int j = 0; j < P2->size; j++) {
                for (; i < P2->vertices[j]; i++) {
                    r += P[i] * P[i];
                }
                r += (P[i] - P2->P[j]) * (P[i] - P2->P[j]);
                i++;
            }
            for (; i < size; i++) {
                r += P[i] * P[i];
            }
        } else { // two full vectors
            for (int i = 0; i < size; i++) {
                r += (P[i] - P2->P[i]) * (P[i] - P2->P[i]);
            }
        }
    }
    return r;
}

long Probabilities::memory() {
    if (vertices) {
        return (sizeof(Probabilities) + long(size) * (sizeof(float) + sizeof(int)));
    } else {
        return (sizeof(Probabilities) + long(size) * sizeof(float));
    }
}

Community::Community() {
    P = 0;
    first_neighbor = 0;
    last_neighbor = 0;
    sub_community_of = -1;
    sub_communities[0] = -1;
    sub_communities[1] = -1;
    sigma = 0.;
    internal_weight = 0.;
    total_weight = 0.;
}

Community::~Community() {
    if (P) {
        delete P;
    }
}


Communities::Communities(Graph* graph, int random_walks_length,
                         long m, igraph_matrix_t *pmerges,
                         igraph_vector_t *pmodularity) {
    max_memory = m;
    memory_used = 0;
    G = graph;
    merges = pmerges;
    mergeidx = 0;
    modularity = pmodularity;

    Probabilities::C = this;
    Probabilities::length = random_walks_length;
    Probabilities::tmp_vector1 = new float[G->nb_vertices];
    Probabilities::tmp_vector2 = new float[G->nb_vertices];
    Probabilities::id = new int[G->nb_vertices];
    for (int i = 0; i < G->nb_vertices; i++) {
        Probabilities::id[i] = 0;
    }
    Probabilities::vertices1 = new int[G->nb_vertices];
    Probabilities::vertices2 = new int[G->nb_vertices];
    Probabilities::current_id = 0;


    members = new int[G->nb_vertices];
    for (int i = 0; i < G->nb_vertices; i++) {
        members[i] = -1;
    }

    H = new Neighbor_heap(G->nb_edges);
    communities = new Community[2 * G->nb_vertices];

// init the n single vertex communities

    if (max_memory != -1) {
        min_delta_sigma = new Min_delta_sigma_heap(G->nb_vertices * 2);
    } else {
        min_delta_sigma = 0;
    }

    for (int i = 0; i < G->nb_vertices; i++) {
        communities[i].this_community = i;
        communities[i].first_member = i;
        communities[i].last_member = i;
        communities[i].size = 1;
        communities[i].sub_community_of = 0;
    }

    nb_communities = G->nb_vertices;
    nb_active_communities = G->nb_vertices;

    for (int i = 0; i < G->nb_vertices; i++)
        for (int j = 0; j < G->vertices[i].degree; j++)
            if (i < G->vertices[i].edges[j].neighbor) {
                communities[i].total_weight += G->vertices[i].edges[j].weight / 2.;
                communities[G->vertices[i].edges[j].neighbor].total_weight += G->vertices[i].edges[j].weight / 2.;
                Neighbor* N = new Neighbor;
                N->community1 = i;
                N->community2 = G->vertices[i].edges[j].neighbor;
                N->delta_sigma = -1. / double(min(G->vertices[i].degree,  G->vertices[G->vertices[i].edges[j].neighbor].degree));
                N->weight = G->vertices[i].edges[j].weight;
                N->exact = false;
                add_neighbor(N);
            }

    if (max_memory != -1) {
        memory_used += min_delta_sigma->memory();
        memory_used += 2 * long(G->nb_vertices) * sizeof(Community);
        memory_used += long(G->nb_vertices) * (2 * sizeof(float) + 3 * sizeof(int)); // the static data of Probabilities class
        memory_used += H->memory() + long(G->nb_edges) * sizeof(Neighbor);
        memory_used += G->memory();
    }

    /*   int c = 0; */
    Neighbor* N = H->get_first();
    if (N == 0) {
        return;    /* this can happen if there are no edges */
    }
    while (!N->exact) {
        update_neighbor(N, compute_delta_sigma(N->community1, N->community2));
        N->exact = true;
        N = H->get_first();
        if (max_memory != -1) {
            manage_memory();
        }
        /* TODO: this could use igraph_progress */
        /*     if(!silent) { */
        /*       c++; */
        /*       for(int k = (500*(c-1))/G->nb_edges + 1; k <= (500*c)/G->nb_edges; k++) { */
        /*  if(k % 50 == 1) {cerr.width(2); cerr << endl << k/ 5 << "% ";} */
        /*  cerr << "."; */
        /*       } */
        /*     } */
    }

}

Communities::~Communities() {
    delete[] members;
    delete[] communities;
    delete H;
    if (min_delta_sigma) {
        delete min_delta_sigma;
    }

    delete[] Probabilities::tmp_vector1;
    delete[] Probabilities::tmp_vector2;
    delete[] Probabilities::id;
    delete[] Probabilities::vertices1;
    delete[] Probabilities::vertices2;
}

float Community::min_delta_sigma() {
    float r = 1.;
    for (Neighbor* N = first_neighbor; N != 0;) {
        if (N->delta_sigma < r) {
            r = N->delta_sigma;
        }
        if (N->community1 == this_community) {
            N = N->next_community1;
        } else {
            N = N->next_community2;
        }
    }
    return r;
}


void Community::add_neighbor(Neighbor* N) { // add a new neighbor at the end of the list
    if (last_neighbor) {
        if (last_neighbor->community1 == this_community) {
            last_neighbor->next_community1 = N;
        } else {
            last_neighbor->next_community2 = N;
        }

        if (N->community1 == this_community) {
            N->previous_community1 = last_neighbor;
        } else {
            N->previous_community2 = last_neighbor;
        }
    } else {
        first_neighbor = N;
        if (N->community1 == this_community) {
            N->previous_community1 = 0;
        } else {
            N->previous_community2 = 0;
        }
    }
    last_neighbor = N;
}

void Community::remove_neighbor(Neighbor* N) {  // remove a neighbor from the list
    if (N->community1 == this_community) {
        if (N->next_community1) {
//      if (N->next_community1->community1 == this_community)
            N->next_community1->previous_community1 = N->previous_community1;
//      else
//  N->next_community1->previous_community2 = N->previous_community1;
        } else {
            last_neighbor = N->previous_community1;
        }
        if (N->previous_community1) {
            if (N->previous_community1->community1 == this_community) {
                N->previous_community1->next_community1 = N->next_community1;
            } else {
                N->previous_community1->next_community2 = N->next_community1;
            }
        } else {
            first_neighbor = N->next_community1;
        }
    } else {
        if (N->next_community2) {
            if (N->next_community2->community1 == this_community) {
                N->next_community2->previous_community1 = N->previous_community2;
            } else {
                N->next_community2->previous_community2 = N->previous_community2;
            }
        } else {
            last_neighbor = N->previous_community2;
        }
        if (N->previous_community2) {
//      if (N->previous_community2->community1 == this_community)
//  N->previous_community2->next_community1 = N->next_community2;
//      else
            N->previous_community2->next_community2 = N->next_community2;
        } else {
            first_neighbor = N->next_community2;
        }
    }
}

void Communities::remove_neighbor(Neighbor* N) {
    communities[N->community1].remove_neighbor(N);
    communities[N->community2].remove_neighbor(N);
    H->remove(N);

    if (max_memory != -1) {
        if (N->delta_sigma == min_delta_sigma->delta_sigma[N->community1]) {
            min_delta_sigma->delta_sigma[N->community1] = communities[N->community1].min_delta_sigma();
            if (communities[N->community1].P) {
                min_delta_sigma->update(N->community1);
            }
        }

        if (N->delta_sigma == min_delta_sigma->delta_sigma[N->community2]) {
            min_delta_sigma->delta_sigma[N->community2] = communities[N->community2].min_delta_sigma();
            if (communities[N->community2].P) {
                min_delta_sigma->update(N->community2);
            }
        }
    }
}

void Communities::add_neighbor(Neighbor* N) {
    communities[N->community1].add_neighbor(N);
    communities[N->community2].add_neighbor(N);
    H->add(N);

    if (max_memory != -1) {
        if (N->delta_sigma < min_delta_sigma->delta_sigma[N->community1]) {
            min_delta_sigma->delta_sigma[N->community1] = N->delta_sigma;
            if (communities[N->community1].P) {
                min_delta_sigma->update(N->community1);
            }
        }

        if (N->delta_sigma < min_delta_sigma->delta_sigma[N->community2]) {
            min_delta_sigma->delta_sigma[N->community2] = N->delta_sigma;
            if (communities[N->community2].P) {
                min_delta_sigma->update(N->community2);
            }
        }
    }
}

void Communities::update_neighbor(Neighbor* N, float new_delta_sigma) {
    if (max_memory != -1) {
        if (new_delta_sigma < min_delta_sigma->delta_sigma[N->community1]) {
            min_delta_sigma->delta_sigma[N->community1] = new_delta_sigma;
            if (communities[N->community1].P) {
                min_delta_sigma->update(N->community1);
            }
        }

        if (new_delta_sigma < min_delta_sigma->delta_sigma[N->community2]) {
            min_delta_sigma->delta_sigma[N->community2] = new_delta_sigma;
            if (communities[N->community2].P) {
                min_delta_sigma->update(N->community2);
            }
        }

        float old_delta_sigma = N->delta_sigma;
        N->delta_sigma = new_delta_sigma;
        H->update(N);

        if (old_delta_sigma == min_delta_sigma->delta_sigma[N->community1]) {
            min_delta_sigma->delta_sigma[N->community1] = communities[N->community1].min_delta_sigma();
            if (communities[N->community1].P) {
                min_delta_sigma->update(N->community1);
            }
        }

        if (old_delta_sigma == min_delta_sigma->delta_sigma[N->community2]) {
            min_delta_sigma->delta_sigma[N->community2] = communities[N->community2].min_delta_sigma();
            if (communities[N->community2].P) {
                min_delta_sigma->update(N->community2);
            }
        }
    } else {
        N->delta_sigma = new_delta_sigma;
        H->update(N);
    }
}

void Communities::manage_memory() {
    while ((memory_used > max_memory) && !min_delta_sigma->is_empty()) {
        int c = min_delta_sigma->get_max_community();
        delete communities[c].P;
        communities[c].P = 0;
        min_delta_sigma->remove_community(c);
    }
}



void Communities::merge_communities(Neighbor* merge_N) {
    int c1 = merge_N->community1;
    int c2 = merge_N->community2;

    communities[nb_communities].first_member = communities[c1].first_member;  // merge the
    communities[nb_communities].last_member = communities[c2].last_member;    // two lists
    members[communities[c1].last_member] = communities[c2].first_member;      // of members

    communities[nb_communities].size = communities[c1].size + communities[c2].size;
    communities[nb_communities].this_community = nb_communities;
    communities[nb_communities].sub_community_of = 0;
    communities[nb_communities].sub_communities[0] = c1;
    communities[nb_communities].sub_communities[1] = c2;
    communities[nb_communities].total_weight = communities[c1].total_weight + communities[c2].total_weight;
    communities[nb_communities].internal_weight = communities[c1].internal_weight + communities[c2].internal_weight + merge_N->weight;
    communities[nb_communities].sigma = communities[c1].sigma + communities[c2].sigma + merge_N->delta_sigma;

    communities[c1].sub_community_of = nb_communities;
    communities[c2].sub_community_of = nb_communities;

// update the new probability vector...

    if (communities[c1].P && communities[c2].P) {
        communities[nb_communities].P = new Probabilities(c1, c2);
    }

    if (communities[c1].P) {
        delete communities[c1].P;
        communities[c1].P = 0;
        if (max_memory != -1) {
            min_delta_sigma->remove_community(c1);
        }
    }
    if (communities[c2].P) {
        delete communities[c2].P;
        communities[c2].P = 0;
        if (max_memory != -1) {
            min_delta_sigma->remove_community(c2);
        }
    }

    if (max_memory != -1) {
        min_delta_sigma->delta_sigma[c1] = -1.;         // to avoid to update the min_delta_sigma for these communities
        min_delta_sigma->delta_sigma[c2] = -1.;         //
        min_delta_sigma->delta_sigma[nb_communities] = -1.;
    }

// update the new neighbors
// by enumerating all the neighbors of c1 and c2

    Neighbor* N1 = communities[c1].first_neighbor;
    Neighbor* N2 = communities[c2].first_neighbor;

    while (N1 && N2) {
        int neighbor_community1;
        int neighbor_community2;

        if (N1->community1 == c1) {
            neighbor_community1 = N1->community2;
        } else {
            neighbor_community1 = N1->community1;
        }
        if (N2->community1 == c2) {
            neighbor_community2 = N2->community2;
        } else {
            neighbor_community2 = N2->community1;
        }

        if (neighbor_community1 < neighbor_community2) {
            Neighbor* tmp = N1;
            if (N1->community1 == c1) {
                N1 = N1->next_community1;
            } else {
                N1 = N1->next_community2;
            }
            remove_neighbor(tmp);
            Neighbor* N = new Neighbor;
            N->weight = tmp->weight;
            N->community1 = neighbor_community1;
            N->community2 = nb_communities;
            N->delta_sigma = (double(communities[c1].size + communities[neighbor_community1].size) * tmp->delta_sigma + double(communities[c2].size) * merge_N->delta_sigma) / (double(communities[c1].size + communities[c2].size + communities[neighbor_community1].size)); //compute_delta_sigma(neighbor_community1, nb_communities);
            N->exact = false;
            delete tmp;
            add_neighbor(N);
        }

        if (neighbor_community2 < neighbor_community1) {
            Neighbor* tmp = N2;
            if (N2->community1 == c2) {
                N2 = N2->next_community1;
            } else {
                N2 = N2->next_community2;
            }
            remove_neighbor(tmp);
            Neighbor* N = new Neighbor;
            N->weight = tmp->weight;
            N->community1 = neighbor_community2;
            N->community2 = nb_communities;
            N->delta_sigma = (double(communities[c1].size) * merge_N->delta_sigma + double(communities[c2].size + communities[neighbor_community2].size) * tmp->delta_sigma) / (double(communities[c1].size + communities[c2].size + communities[neighbor_community2].size)); //compute_delta_sigma(neighbor_community2, nb_communities);
            N->exact = false;
            delete tmp;
            add_neighbor(N);
        }

        if (neighbor_community1 == neighbor_community2) {
            Neighbor* tmp1 = N1;
            Neighbor* tmp2 = N2;
            bool exact = N1->exact && N2->exact;
            if (N1->community1 == c1) {
                N1 = N1->next_community1;
            } else {
                N1 = N1->next_community2;
            }
            if (N2->community1 == c2) {
                N2 = N2->next_community1;
            } else {
                N2 = N2->next_community2;
            }
            remove_neighbor(tmp1);
            remove_neighbor(tmp2);
            Neighbor* N = new Neighbor;
            N->weight = tmp1->weight + tmp2->weight;
            N->community1 = neighbor_community1;
            N->community2 = nb_communities;
            N->delta_sigma = (double(communities[c1].size + communities[neighbor_community1].size) * tmp1->delta_sigma + double(communities[c2].size + communities[neighbor_community1].size) * tmp2->delta_sigma - double(communities[neighbor_community1].size) * merge_N->delta_sigma) / (double(communities[c1].size + communities[c2].size + communities[neighbor_community1].size));
            N->exact = exact;
            delete tmp1;
            delete tmp2;
            add_neighbor(N);
        }
    }


    if (!N1) {
        while (N2) {
//      double delta_sigma2 = N2->delta_sigma;
            int neighbor_community;
            if (N2->community1 == c2) {
                neighbor_community = N2->community2;
            } else {
                neighbor_community = N2->community1;
            }
            Neighbor* tmp = N2;
            if (N2->community1 == c2) {
                N2 = N2->next_community1;
            } else {
                N2 = N2->next_community2;
            }
            remove_neighbor(tmp);
            Neighbor* N = new Neighbor;
            N->weight = tmp->weight;
            N->community1 = neighbor_community;
            N->community2 = nb_communities;
            N->delta_sigma = (double(communities[c1].size) * merge_N->delta_sigma + double(communities[c2].size + communities[neighbor_community].size) * tmp->delta_sigma) / (double(communities[c1].size + communities[c2].size + communities[neighbor_community].size)); //compute_delta_sigma(neighbor_community, nb_communities);
            N->exact = false;
            delete tmp;
            add_neighbor(N);
        }
    }
    if (!N2) {
        while (N1) {
//      double delta_sigma1 = N1->delta_sigma;
            int neighbor_community;
            if (N1->community1 == c1) {
                neighbor_community = N1->community2;
            } else {
                neighbor_community = N1->community1;
            }
            Neighbor* tmp = N1;
            if (N1->community1 == c1) {
                N1 = N1->next_community1;
            } else {
                N1 = N1->next_community2;
            }
            remove_neighbor(tmp);
            Neighbor* N = new Neighbor;
            N->weight = tmp->weight;
            N->community1 = neighbor_community;
            N->community2 = nb_communities;
            N->delta_sigma = (double(communities[c1].size + communities[neighbor_community].size) * tmp->delta_sigma + double(communities[c2].size) * merge_N->delta_sigma) / (double(communities[c1].size + communities[c2].size + communities[neighbor_community].size)); //compute_delta_sigma(neighbor_community, nb_communities);
            N->exact = false;
            delete tmp;
            add_neighbor(N);
        }
    }

    if (max_memory != -1) {
        min_delta_sigma->delta_sigma[nb_communities] = communities[nb_communities].min_delta_sigma();
        min_delta_sigma->update(nb_communities);
    }

    nb_communities++;
    nb_active_communities--;
}

double Communities::merge_nearest_communities() {
    Neighbor* N = H->get_first();
    while (!N->exact) {
        update_neighbor(N, compute_delta_sigma(N->community1, N->community2));
        N->exact = true;
        N = H->get_first();
        if (max_memory != -1) {
            manage_memory();
        }
    }

    double d = N->delta_sigma;
    remove_neighbor(N);

    merge_communities(N);
    if (max_memory != -1) {
        manage_memory();
    }

    if (merges) {
        MATRIX(*merges, mergeidx, 0) = N->community1;
        MATRIX(*merges, mergeidx, 1) = N->community2;
        mergeidx++;
    }

    if (modularity) {
        float Q = 0.;
        for (int i = 0; i < nb_communities; i++) {
            if (communities[i].sub_community_of == 0) {
                Q += (communities[i].internal_weight - communities[i].total_weight * communities[i].total_weight / G->total_weight) / G->total_weight;
            }
        }
        VECTOR(*modularity)[mergeidx] = Q;
    }

    delete N;

    /* This could use igraph_progress */
    /*   if(!silent) { */
    /*     for(int k = (500*(G->nb_vertices - nb_active_communities - 1))/(G->nb_vertices-1) + 1; k <= (500*(G->nb_vertices - nb_active_communities))/(G->nb_vertices-1); k++) { */
    /*       if(k % 50 == 1) {cerr.width(2); cerr << endl << k/ 5 << "% ";} */
    /*       cerr << "."; */
    /*     } */
    /*   } */
    return d;
}

double Communities::compute_delta_sigma(int community1, int community2) {
    if (!communities[community1].P) {
        communities[community1].P = new Probabilities(community1);
        if (max_memory != -1) {
            min_delta_sigma->update(community1);
        }
    }
    if (!communities[community2].P) {
        communities[community2].P = new Probabilities(community2);
        if (max_memory != -1) {
            min_delta_sigma->update(community2);
        }
    }

    return communities[community1].P->compute_distance(communities[community2].P) * double(communities[community1].size) * double(communities[community2].size) / double(communities[community1].size + communities[community2].size);
}

}
}    /* end of namespaces */
