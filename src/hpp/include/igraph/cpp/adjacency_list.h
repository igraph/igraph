/* vim:set ts=4 sw=4 sts=4 et: */

#ifndef IGRAPHPP_ADJACENCY_LIST_H
#define IGRAPHPP_ADJACENCY_LIST_H

#include <stdexcept>
#include <igraph/igraph_adjlist.h>
#include <igraph/cpp/graph.h>
#include <igraph/cpp/types.h>
#include <igraph/cpp/util.h>

namespace igraph {

class AdjacencyList : private noncopyable {
private:
    /// The encapsulated adjacency list object
    igraph_adjlist_t m_adjList;

public:
    /// Constructs an adjacency list representation for the given graph
    explicit AdjacencyList(const Graph& graph, NeighborMode mode) : m_adjList() {
        IGRAPH_TRY(igraph_adjlist_init(graph.c_graph(), &m_adjList, mode));
    }

    /// Destroys the adjacency list
    ~AdjacencyList() {
        igraph_adjlist_destroy(&m_adjList);
    }

    /// Queries the neighbors of a given vertex in the adjacency list
    Vector get(long index) {
        return Vector(igraph_adjlist_get(&m_adjList, index));
    }

    /// Queries the neighbors of a given vertex in the adjacency list (const)
    const Vector get(long index) const {
        return Vector(igraph_adjlist_get(&m_adjList, index));
    }

    /// Simplifies the adjacency list, i.e. removes loop and multiple edges
    void simplify() {
        igraph_adjlist_simplify(&m_adjList);
    }

    /// Sorts each vector in the adjacency list
    void sort() {
        igraph_adjlist_sort(&m_adjList);
    }
};

}       // end of namespaces

#endif  // IGRAPHPP_ADJACENCY_LIST_H


