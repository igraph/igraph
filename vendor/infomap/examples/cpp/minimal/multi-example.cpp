/**********************************************************************************

 Infomap software package for multi-level network clustering

 Copyright (c) 2013, 2014 Daniel Edler, Martin Rosvall

 For more information, see <http://www.mapequation.org>


 This file is part of Infomap software package.

 Infomap software package is free software: you can redistribute it and/or modify
 it under the terms of the GNU Affero General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 Infomap software package is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU Affero General Public License for more details.

 You should have received a copy of the GNU Affero General Public License
 along with Infomap software package.  If not, see <http://www.gnu.org/licenses/>.

**********************************************************************************/

#include <iostream>
#include <sstream>
#include <string>

#include <Infomap.h>

using std::string;

void printClusters(infomap::HierarchicalNetwork & tree) {
    std::cout << "\nClusters:\n#layer node clusterIndex:\n";
    for (infomap::LeafIterator leafIt(&tree.getRootNode()); !leafIt.isEnd(); ++leafIt) {
        std::cout << leafIt->stateIndex << " " << leafIt->physIndex << " " << leafIt.moduleIndex() << '\n';
    }
}

int main(int argc, char** argv)
{
	infomap::MemInfomap infomapWrapper("--two-level -N2 --expanded");

	// from (layer physical) to (layer physical) weight
	infomapWrapper.addMultiplexLink(2, 1, 1, 2, 1.0);
	infomapWrapper.addMultiplexLink(1, 2, 2, 1, 1.0);
	infomapWrapper.addMultiplexLink(3, 2, 2, 3, 1.0);

	infomapWrapper.run();

	printClusters(infomapWrapper.tree);
}
