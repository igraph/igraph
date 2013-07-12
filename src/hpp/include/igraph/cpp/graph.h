/* vim:set ts=4 sw=4 sts=4 et: */

#ifndef IGRAPHPP_GRAPH_H
#define IGRAPHPP_GRAPH_H

#include <cassert>
#include <memory>
#include <igraph/igraph_constructors.h>
#include <igraph/igraph_datatype.h>
#include <igraph/igraph_foreign.h>
#include <igraph/igraph_interface.h>
#include <igraph/cpp/attributes.h>
#include <stdexcept>
#include <vector>

#include <igraph/cpp/error.h>
#include <igraph/cpp/types.h>
#include <igraph/cpp/vector.h>

namespace igraph {

class Edge;
class EdgeSelector;
class Vertex;
class VertexSelector;

/// C++-style wrapper around an igraph_t object
class Graph {
private:
    /// The igraph_t instance encapsulated by the wrapper
    igraph_t* m_pGraph;

public:
    /*****************************/
    /* Constructors, destructors */
    /*****************************/

    /// Constructs an empty graph
    Graph(long numVertices = 0, bool directed = false) : m_pGraph(new igraph_t) {
        IGRAPH_TRY(igraph_empty(m_pGraph, numVertices, directed));
    }
    
    /// Constructs a wrapper that wraps the given igraph_t instance
    /**
     * The ownership of the wrapped instance is stolen by the wrapper.
     * The caller should not destroy the graph on its own, ever;
     * the wrapper should be destroyed instead.
     *
     * This function never throws an exception.
     */
    Graph(igraph_t* graph) throw() : m_pGraph(graph) {
        assert(graph);
    }

    /// Constructs a wrapper that wraps the given igraph_t instance
    /**
     * The ownership of the wrapped instance is stolen by the wrapper.
     * The caller should not destroy the graph on its own, ever;
     * the wrapper should be destroyed instead.
     *
     * This function never throws an exception.
     */
    Graph(std::auto_ptr<igraph_t> graph) throw() : m_pGraph(graph.release()) {
    }

    /// Copy constructor
    Graph(const Graph& other) : m_pGraph(new igraph_t) {
        IGRAPH_TRY(igraph_copy(m_pGraph, other.m_pGraph));
    }

    /// Destroys the graph
    ~Graph();

    /******************/
    /* Static methods */
    /******************/

    /********************/
    /* Instance methods */
    /********************/

    /// Adds a single edge to the graph
    void addEdge(integer_t source, integer_t target);

    /// Adds a list of edges to the graph
    void addEdges(const Vector& edges);

    /// Adds a single vertex to the graph
    void addVertex();

    /// Adds the given number of vertices to the graph
    void addVertices(long numVertices);

    /// Returns whether the two vertices are connected
    bool areConnected(long u, long v) const;

    /// Extracts a pointer to the encapsulated graph object
    igraph_t* c_graph() { return m_pGraph; }

    /// Extracts a pointer to the encapsulated graph object (const)
    const igraph_t* c_graph() const { return m_pGraph; }

    /// Returns the degrees of some vertices
    Vector degree(const VertexSelector& vids,
            NeighborMode mode = IGRAPH_ALL, bool loops = false) const;
    /// Returns the degrees of some vertices
    void degree(Vector* result, const VertexSelector& vids,
            NeighborMode mode = IGRAPH_ALL, bool loops = false) const;

    /// Deletes some edges from the graph
    void deleteEdges(const EdgeSelector& es);

    /// Returns the head and tail vertices of an edge
    void edge(integer_t eid, integer_t* from, integer_t* to) const;
    /// Returns the edge with the given index
    Edge edge(integer_t eid);

    /// Returns the number of edges in the graph
    integer_t ecount() const { return igraph_ecount(m_pGraph); }

    /// Returns a copy of the value of the given graph attribute
    AttributeValue getAttribute(const std::string& attribute) const;

    /// Returns the edge list of the graph
    Vector getEdgelist(bool bycol=false) const;
    /// Returns the edge list of the graph
    void getEdgelist(Vector* result, bool bycol=false) const;

    /// Returns the ID of an arbitrary edge between the two given nodes
    integer_t getEid(integer_t source, integer_t target, bool directed=true,
            bool error=false) const;

    /// Returns whether the graph has the given graph attribute
    bool hasAttribute(const std::string& attribute) const;

    /// Returns the edges incident on a given vertex
    void incident(Vector* result, long int vertex, NeighborMode mode = IGRAPH_OUT) const;
    /// Returns the edges incident on a given vertex
    Vector incident(long int vertex, NeighborMode mode = IGRAPH_OUT) const;

    /// Returns whether the graph is directed
    bool isDirected() const;

    /// Checks whether the graph is a simple graph
    bool isSimple() const;

    /// Returns the neighbors of a vertex
    void neighbors(Vector* result, long int vertex, NeighborMode mode = IGRAPH_OUT) const;
    /// Returns the neighbors of a vertex
    Vector neighbors(long int vertex, NeighborMode mode = IGRAPH_OUT) const;

    /// Sets the value of the given graph attribute
    void setAttribute(const std::string& attribute, const AttributeValue& value);

    /// Removes loop and/or multiple edges from the graph
    void simplify(bool multiple=true, bool loops=true);

    /// Returns the number of vertices in the graph
    integer_t vcount() const { return igraph_vcount(m_pGraph); }

    /// Returns the vertex with the given index in the graph
    Vertex vertex(integer_t vid);

    /*************/
    /* Operators */
    /*************/

    /// Disjoint union of two graphs
    Graph operator+(const Graph& other) const;

    /// Assignment operator
    Graph& operator=(const Graph& other) {
        if (&other == this)
            return *this;

        igraph_t new_graph;
        IGRAPH_TRY(igraph_copy(&new_graph, other.m_pGraph));

        if (m_pGraph) {
            igraph_destroy(m_pGraph);
        } else {
            m_pGraph = new igraph_t;
        }
        *m_pGraph = new_graph;

        return *this;
    }

    /// Retrieves the value of the given graph attribute
    /**
     * This method works similar to \c std::map<>.operator[]: if the attribute
     * is found, its value is returned; if the attribute is not found, a new
     * attribute will be created with the default constructor of \c AttributeValue 
     * and this will be returned. Therefore, this operator won't work on const
     * graphs.
     */
    AttributeValue& operator[](const std::string& attribute);

private:
    /// Returns a pointer to the attribute holder of the graph
    AttributeHolder* getAttributeHolder();

    /// Returns a pointer to the attribute holder of the graph (const)
    const AttributeHolder* getAttributeHolder() const;

    friend class Edge;
    friend class EdgeSelector;
    friend class Vertex;
    friend class VertexSelector;
};

}       // end of namespaces

#endif  // IGRAPHPP_GRAPH_H
