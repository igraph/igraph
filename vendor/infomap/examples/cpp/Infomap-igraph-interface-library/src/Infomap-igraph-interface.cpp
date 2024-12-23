#include "Infomap-igraph-interface.h"

namespace infomap
{
void igraphToInfomapNetwork(infomap::Network& network, const igraph_t* graph,
	const igraph_vector_t *e_weights, const igraph_vector_t *v_weights)
{
	unsigned int numNodes = static_cast<unsigned int>(igraph_vcount(graph));
	unsigned int numLinks = static_cast<unsigned int>(igraph_ecount(graph));
		
	double linkWeight = 1.0;
	igraph_integer_t from, to;
	
	for (unsigned int i = 0; i < numLinks; ++i)
	{
		igraph_edge(graph, i, &from, &to);
		linkWeight = e_weights ? (double)VECTOR(*e_weights)[i] : 1.0;

		network.addLink(static_cast<unsigned int>(from), static_cast<unsigned int>(to), linkWeight);
		
		double nodeWeight = 1.0;
		if (v_weights) {
			for (unsigned int i = 0; i < numNodes; ++i) {
				nodeWeight = (double)VECTOR(*v_weights)[i];
			}
		}
	}
}
}