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
#include <Infomap.h>

void printClusters(infomap::InfomapWrapper& infomap)
{
  std::cout << "\nClusters:\n#originalIndex clusterIndex:\n";

  for (auto it = infomap.iterTree(); !it.isEnd(); ++it) {
    if (it->isLeaf())
      std::cout << it->physicalId << " " << it.moduleIndex() << '\n';
  }
}

int main(int argc, char** argv)
{
  std::string inputFilename = "../../../ninetriangles.net";
  std::cout << "Cluster '" << inputFilename << "'...\n";

  infomap::InfomapWrapper infomapWrapper("--two-level -N2");

  auto& network = infomapWrapper.network();
  network.addLink(0, 1);
  network.addLink(0, 2);
  network.addLink(0, 3);
  network.addLink(1, 0);
  network.addLink(1, 2);
  network.addLink(2, 1);
  network.addLink(2, 0);
  network.addLink(3, 0);
  network.addLink(3, 4);
  network.addLink(3, 5);
  network.addLink(4, 3);
  network.addLink(4, 5);
  network.addLink(5, 4);
  network.addLink(5, 3);

  infomapWrapper.run();

  std::cout << "Done. Codelength: " << infomapWrapper.codelength() << "\n";
  printClusters(infomapWrapper);
}
