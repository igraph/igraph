/* vim:set ts=4 sw=4 sts=4 et: */

#ifndef IGRAPHPP_EDGE_SELECTOR_H
#define IGRAPHPP_EDGE_SELECTOR_H

#include <igraph/cpp/graph.h>
#include <igraph/cpp/types.h>
#include <igraph/cpp/vector.h>

namespace igraph {

/// C++-style wrapper around an igraph_es_t
class EdgeSelector {
private:
    /// The igraph_es_t instance encapsulated by the wrapper
    igraph_es_t m_es;

    /// The graph the edge selector refers to
    Graph* m_pGraph;

public:
    /*****************************/
    /* Constructors, destructors */
    /*****************************/

    /// Constructs an empty edge selector
    explicit EdgeSelector(Graph* pGraph = 0) {
        setGraph(pGraph);
        IGRAPH_TRY(igraph_es_none(&m_es));
    }

    /// Constructs an edge selector that selects a single edge
    EdgeSelector(integer_t eid, Graph* pGraph = 0) {
        setGraph(pGraph);
        IGRAPH_TRY(igraph_es_1(&m_es, eid));
    }

    /// Constructs an edge selector that handles a vector as an edge selector
    EdgeSelector(Vector* vector, Graph* pGraph = 0) {
        setGraph(pGraph);
        IGRAPH_TRY(igraph_es_vector(&m_es, vector->c_vector()));
    }

    /// Constructs a wrapper that wraps the given igraph_es_t instance
    /**
     * The ownership of the wrapped instance is stolen by the wrapper.
     * The caller should not destroy the edge selector on its own, ever;
     * the wrapper should be destroyed instead.
     */
    EdgeSelector(igraph_es_t es, Graph* pGraph = 0) : m_es(es), m_pGraph(pGraph) {}

    /// Destroys the edge selector
    ~EdgeSelector() {
        igraph_es_destroy(&m_es);
    }

    /******************/
    /* Static methods */
    /******************/

    /// Creates an edge selector that selects all edges
    static EdgeSelector All(EdgeOrderType order = IGRAPH_EDGEORDER_ID, Graph* pGraph = 0) {
        igraph_es_t es;
        IGRAPH_TRY(igraph_es_all(&es, order));
        return EdgeSelector(es, pGraph);
    }

    /// Creates an edge selector from multiple edges defined by their endpoints
    static EdgeSelector Pairs(Vector& vector, bool directed=true, Graph* pGraph = 0) {
        igraph_es_t es;
        IGRAPH_TRY(igraph_es_pairs(&es, vector.c_vector(), directed));
        return EdgeSelector(es, pGraph);
    }

    /********************/
    /* Instance methods */
    /********************/

    /// Returns a pointer to the encapsulated igraph_es_t instance
    igraph_es_t* c_es() {
        return &m_es;
    }

    /// Returns a pointer to the encapsulated igraph_es_t instance (const)
    const igraph_es_t* c_es() const {
        return &m_es;
    }

    /// Returns the graph the edge selector refers to (const version)
    const Graph* getGraph() const {
        return m_pGraph;
    }

    /// Returns the graph the edge selector refers to
    Graph* getGraph() {
        return m_pGraph;
    }

    /// Checks whether an edge selector includes all edges
    bool isAll() {
        return igraph_es_is_all(&m_es);
    }

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
    EdgeSelector& operator=(const EdgeSelector&);
};

inline EdgeSelector E(Graph* graph, EdgeOrderType order = IGRAPH_EDGEORDER_ID) {
    return EdgeSelector::All(order, graph);
}

inline EdgeSelector E(Graph& graph, EdgeOrderType order = IGRAPH_EDGEORDER_ID) {
    return EdgeSelector::All(order, &graph);
}

}       // end of namespaces

#endif  // IGRAPHPP_EDGE_SELECTOR_H

