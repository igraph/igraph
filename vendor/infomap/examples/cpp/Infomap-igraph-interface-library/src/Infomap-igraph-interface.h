#include <io/Network.h>
#include <igraph_interface.h>

namespace infomap
{
	/**
	* Copy an igraph network to an Infomap network
	*
	* @param network The Infomap network destination
	* @param graph The igraph network source
	* e_weights: Numeric vector giving the weights of the edges. If it is a NULL pointer then all edges
	* will have equal weights. The weights are expected to be positive.
	* v_weights: Numeric vector giving the weights of the vertices. If it is a NULL pointer then all
	* vertices will have equal weights. The weights are expected to be positive.
	*/
	void igraphToInfomapNetwork(infomap::Network& network, const igraph_t* graph,
		const igraph_vector_t *e_weights = NULL, const igraph_vector_t *v_weights = NULL);

	/**
	* @param membership Pointer to a vector. The membership vector is stored here.
	*/
	// void getClusterMembership(igraph_vector_t *membership);
}