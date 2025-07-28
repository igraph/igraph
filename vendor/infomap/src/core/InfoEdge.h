/*******************************************************************************
 Infomap software package for multi-level network clustering
 Copyright (c) 2013, 2014 Daniel Edler, Anton Holmgren, Martin Rosvall

 This file is part of the Infomap software package.
 See file LICENSE_GPLv3.txt for full license details.
 For more information, see <http://www.mapequation.org>
 ******************************************************************************/

#ifndef INFOEDGE_H_
#define INFOEDGE_H_

#include <iostream>

namespace infomap {

struct EdgeData {
public:
  EdgeData() = default;

  EdgeData(double weight, double flow) : weight(weight), flow(flow) { }

  double weight;
  double flow;
};

class InfoNode;

class InfoEdge {
public:
  InfoEdge(InfoNode& source, InfoNode& target, double weight, double flow)
      : data(weight, flow),
        source(&source),
        target(&target) { }

  InfoNode& other(InfoNode& node) const;

  friend std::ostream& operator<<(std::ostream& out, const InfoEdge& edge);

  EdgeData data;
  InfoNode* source;
  InfoNode* target;
};

} // namespace infomap

#endif // INFOEDGE_H_
