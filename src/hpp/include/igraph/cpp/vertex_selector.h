/* vim:set ts=4 sw=4 sts=4 et: */

#ifndef IGRAPHPP_VERTEX_SELECTOR_H
#define IGRAPHPP_VERTEX_SELECTOR_H

#include <igraph/cpp/graph.h>
#include <igraph/cpp/types.h>
#include <igraph/cpp/vector.h>

namespace igraph {

/// C++-style wrapper around an igraph_vs_t
class VertexSelector {
private:
    /// The igraph_vs_t instance encapsulated by the wrapper
    igraph_vs_t m_vs;

    /// The graph the vertex selector refers to
    Graph* m_pGraph;

public:
    /*****************************/
    /* Constructors, destructors */
    /*****************************/

    /// Constructs an empty vertex selector
    explicit VertexSelector(Graph* pGraph = 0) {
        setGraph(pGraph);
        IGRAPH_TRY(igraph_vs_none(&m_vs));
    }

    /// Constructs a vertex selector that selects a single vertex
    VertexSelector(integer_t vid, Graph* pGraph = 0) {
        setGraph(pGraph);
        IGRAPH_TRY(igraph_vs_1(&m_vs, vid));
    }

    /// Constructs a vertex selector that selects a range of vertices
    VertexSelector(integer_t from, integer_t to, Graph* pGraph = 0) {
        setGraph(pGraph);
        IGRAPH_TRY(igraph_vs_seq(&m_vs, from, to));
    }

    /// Constructs a vertex selector that handles a vector as a vertex selector
    VertexSelector(Vector* vector, Graph* pGraph = 0) {
        setGraph(pGraph);
        IGRAPH_TRY(igraph_vs_vector(&m_vs, vector->c_vector()));
    }

    /// Constructs a wrapper that wraps the given igraph_vs_t instance
    /**
     * The ownership of the wrapped instance is stolen by the wrapper.
     * The caller should not destroy the vertex selector on its own, ever;
     * the wrapper should be destroyed instead.
     */
    VertexSelector(igraph_vs_t vs, Graph* pGraph = 0) : m_vs(vs), m_pGraph(pGraph) {}

    /// Destroys the vertex selector
    ~VertexSelector() {
        igraph_vs_destroy(&m_vs);
    }

    /******************/
    /* Static methods */
    /******************/

    /// Creates a vertex selector that selects all vertices
    static VertexSelector All(Graph* pGraph = 0) {
        igraph_vs_t vs;
        IGRAPH_TRY(igraph_vs_all(&vs));
        return VertexSelector(vs, pGraph);
    }

    /********************/
    /* Instance methods */
    /********************/

    /// Returns a pointer to the encapsulated igraph_vs_t instance
    igraph_vs_t* c_vs() {
        return &m_vs;
    }

    /// Returns a pointer to the encapsulated igraph_vs_t instance (const)
    const igraph_vs_t* c_vs() const {
        return &m_vs;
    }
    
    /// Returns the value of the given vertex attribute for the vertices selected by the selector
    AttributeValueVector getAttribute(const std::string& attribute) const;

    /// Returns the graph the vertex selector refers to (const version)
    const Graph* getGraph() const {
        return m_pGraph;
    }

    /// Returns the graph the vertex selector refers to
    Graph* getGraph() {
        return m_pGraph;
    }

    /// Returns whether the graph has the given vertex attribute
    bool hasAttribute(const std::string& attribute) const {
        assert(m_pGraph != 0);
        return m_pGraph->getAttributeHolder()->hasVertexAttribute(attribute);
    }

    /// Returns whether the vertex selector selects all vertices exactly once
    bool isAll() const {
        return igraph_vs_is_all(&m_vs);
    }

    /// Sets the values of the given vertex attribute for the vertices selected by the selector
    void setAttribute(const std::string& attribute,
            const AttributeValueVector& values) const;

    /// Sets the graph the vertex selector refers to
    void setGraph(Graph* pGraph) {
        m_pGraph = pGraph;
    }

    /*************/
    /* Operators */
    /*************/

    /*****************/
    /* Private stuff */
    /*****************/

private:
    /// Assignment operator (intentionally unimplemented)
    VertexSelector& operator=(const VertexSelector&);
};

inline VertexSelector V(Graph* graph) {
    return VertexSelector::All(graph);
}

inline VertexSelector V(Graph& graph) {
    return VertexSelector::All(&graph);
}

}       // end of namespaces

#endif  // IGRAPHPP_VERTEX_SELECTOR_H

