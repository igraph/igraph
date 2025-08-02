/*******************************************************************************
 Infomap software package for multi-level network clustering
 Copyright (c) 2013, 2014 Daniel Edler, Anton Holmgren, Martin Rosvall

 This file is part of the Infomap software package.
 See file LICENSE_GPLv3.txt for full license details.
 For more information, see <http://www.mapequation.org>
 ******************************************************************************/

#include "InfoEdge.h"
#include "InfoNode.h"

namespace infomap {

InfoNode& infomap::InfoEdge::other(InfoNode& node) const
{
  return (node == *source) ? *target : *source;
}

std::ostream& operator<<(std::ostream& out, const InfoEdge& edge)
{
  return out << "(" << *edge.source << ") -> (" << *edge.target << "), flow: "
             << edge.data.flow;
}

} // namespace infomap
