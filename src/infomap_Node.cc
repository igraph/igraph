/* -*- mode: C -*-  */
/* vim:set ts=4 sw=4 sts=4 et: */
/*
   IGraph library.
   Copyright (C) 2011-2012  Gabor Csardi <csardi.gabor@gmail.com>
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

#include "infomap_Node.h"

using namespace std;

Node::Node() {
    exit = 0.0;
    size = 0.0;
    selfLink = 0.0;
}

Node::Node(int nodenr, double tpweight) {
    teleportWeight = tpweight;
    exit = 0.0;
    size = 0.0;
    selfLink = 0.0;
    members.push_back(nodenr); // members = [nodenr]
}

void cpyNode(Node *newNode, Node *oldNode) {
    newNode->exit = oldNode->exit;
    newNode->size = oldNode->size;
    newNode->teleportWeight = oldNode->teleportWeight;
    newNode->danglingSize   = oldNode->danglingSize;

    int Nmembers = oldNode->members.size();
    newNode->members = vector<int>(Nmembers);
    for (int i = 0; i < Nmembers; i++) {
        newNode->members[i] = oldNode->members[i];
    }

    newNode->selfLink = oldNode->selfLink;

    int NoutLinks = oldNode->outLinks.size();
    newNode->outLinks = vector<pair<int, double> >(NoutLinks);
    for (int i = 0; i < NoutLinks; i++) {
        newNode->outLinks[i].first = oldNode->outLinks[i].first;
        newNode->outLinks[i].second = oldNode->outLinks[i].second;
    }

    int NinLinks = oldNode->inLinks.size();
    newNode->inLinks = vector<pair<int, double> >(NinLinks);
    for (int i = 0; i < NinLinks; i++) {
        newNode->inLinks[i].first = oldNode->inLinks[i].first;
        newNode->inLinks[i].second = oldNode->inLinks[i].second;
    }

}

