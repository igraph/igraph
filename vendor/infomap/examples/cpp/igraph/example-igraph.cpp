#include <igraph.h>
#include <Infomap.h>
#include <Infomap-igraph-interface.h>

#include <iostream>

void printClusters(infomap::HierarchicalNetwork& tree)
{
	std::cout << "\nClusters:\n#originalIndex clusterIndex:\n";

	for (infomap::LeafIterator leafIt(&tree.getRootNode()); !leafIt.isEnd(); ++leafIt)
		std::cout << leafIt->originalLeafIndex << " " << leafIt->parentNode->parentIndex << '\n';
}

void partitionNetwork(igraph_t* graph)
{
	infomap::Config config = infomap::init("--two-level");

	infomap::Network network(config);

	infomap::igraphToInfomapNetwork(network, graph);

	infomap::HierarchicalNetwork resultNetwork(config);

	infomap::run(network, resultNetwork);

	printClusters(resultNetwork);
}

int main() {
	igraph_t graph;

	/* Zachary Karate club network */
	igraph_small(&graph, 0, IGRAPH_UNDIRECTED,
               0,  1,  0,  2,  0,  3,  0,  4,  0,  5, //0,  5, 0,  5, 0,  5,
               0,  6,  0,  7,  0,  8,  0, 10,  0, 11,
               0, 12,  0, 13,  0, 17,  0, 19,  0, 21,
               0, 31,  1,  2,  1,  3,  1,  7,  1, 13,
               1, 17,  1, 19,  1, 21,  1, 30,  2,  3,
               2,  7,  2,  8,  2,  9,  2, 13,  2, 27,
               2, 28,  2, 32,  3,  7,  3, 12,  3, 13,
               4,  6,  4, 10,  5,  6,  5, 10,  5, 16,
               6, 16,  8, 30,  8, 32,  8, 33,  9, 33,
               13, 33, 14, 32, 14, 33, 15, 32, 15, 33,
               18, 32, 18, 33, 19, 33, 20, 32, 20, 33,
               22, 32, 22, 33, 23, 25, 23, 27, 23, 29,
               23, 32, 23, 33, 24, 25, 24, 27, 24, 31,
               25, 31, 26, 29, 26, 33, 27, 33, 28, 31,
               28, 33, 29, 32, 29, 33, 30, 32, 30, 33,
               31, 32, 31, 33, 32, 33,
               -1);

	igraph_integer_t diameter;
	igraph_diameter(&graph, &diameter, 0, 0, 0, IGRAPH_UNDIRECTED, 1);
	std::cout << "Created the Zachary Karate club network with diameter " << diameter << ".\n";

	std::cout << "Partition the network with Infomap..." << std::endl;
	partitionNetwork(&graph);

	igraph_destroy(&graph);

	return 0;
}