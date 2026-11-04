/*******************************************************************************
 Infomap software package for multi-level network clustering
 Copyright (c) 2013, 2014 Daniel Edler, Anton Holmgren, Martin Rosvall

 This file is part of the Infomap software package.
 See file LICENSE_GPLv3.txt for full license details.
 For more information, see <http://www.mapequation.org>
 ******************************************************************************/

#ifndef NETWORK_BUILDER_H_
#define NETWORK_BUILDER_H_

#include <string>

namespace infomap {

class Network;

void buildNetworkFromInput(Network& network, const std::string& filename);
void buildMetaDataFromInput(Network& network, const std::string& filename);

} // namespace infomap

#endif // NETWORK_BUILDER_H_
