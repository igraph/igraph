/* -*- mode: C -*-  */
/* vim:set ts=2 sw=2 sts=2 et: */
/* 
   IGraph library.
   Copyright (C) 2011  Gabor Csardi <csardi@rmki.kfki.hu>
   MTA RMKI, Konkoly-Thege Miklos st. 29-33, Budapest 1121, Hungary
   
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

#ifndef NODE_H
#define NODE_H

#include <vector>
#include <utility>

#include "igraph_interface.h"

class Node;
using namespace std;

class Node{
 public:
  
  Node();
  Node(int modulenr,double tpweight);

  vector<int> members;
  vector< pair<int,double> > inLinks;
  vector< pair<int,double> > outLinks;
  double selfLink;

  double teleportWeight;
  double danglingSize;
  double exit;
  double size;
};

void cpyNode(Node *newNode, Node *oldNode);

#endif
