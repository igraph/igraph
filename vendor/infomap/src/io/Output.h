/*******************************************************************************
 Infomap software package for multi-level network clustering
 Copyright (c) 2013, 2014 Daniel Edler, Anton Holmgren, Martin Rosvall

 This file is part of the Infomap software package.
 See file LICENSE_GPLv3.txt for full license details.
 For more information, see <http://www.mapequation.org>
 ******************************************************************************/

#ifndef OUTPUT_H_
#define OUTPUT_H_

#include <map>
#include <string>
#include <utility>

namespace infomap {

class InfomapBase;
class StateNetwork;

std::string writeTree(InfomapBase&, const StateNetwork&, const std::string&, bool states);

std::string writeFlowTree(InfomapBase&, const StateNetwork&, const std::string&, bool states);

std::string writeNewickTree(InfomapBase&, const std::string&, bool states);

std::string writeJsonTree(InfomapBase&, const StateNetwork&, const std::string&, bool states, bool writeLinks);

std::string writeCsvTree(InfomapBase&, const StateNetwork&, const std::string&, bool states);

std::string writeClu(InfomapBase&, const StateNetwork&, const std::string&, bool states, int moduleIndexLevel);

} // namespace infomap

#endif // OUTPUT_H_
