/* vim:set ts=4 sw=4 sts=4 et: */

#ifndef IGRAPHPP_EDGE_ITERATOR_H
#define IGRAPHPP_EDGE_ITERATOR_H

#include <igraph/cpp/edge.h>
#include <igraph/cpp/edge_selector.h>
#include <igraph/cpp/types.h>

namespace igraph {

/// C++-style wrapper around an igraph_eit_t
class EdgeIterator {
private:
    /// The igraph_eit_t instance encapsulated by the wrapper
    igraph_eit_t m_eit;

    /// The graph the edge selector refers to
    Graph* m_pGraph;

public:
    /*****************************/
    /* Constructors, destructors */
    /*****************************/

    /// Constructs an edge iterator from a selector
    explicit EdgeIterator(EdgeSelector& es) : m_pGraph(es.getGraph()) {
        IGRAPH_TRY(igraph_eit_create(m_pGraph->c_graph(), *es.c_es(), &m_eit));
    }

    /// Constructs a wrapper that wraps the given igraph_eit_t instance
    /**
     * The ownership of the wrapped instance is stolen by the wrapper.
     * The caller should not destroy the edge iterator on its own, ever;
     * the wrapper should be destroyed instead.
     */
    EdgeIterator(igraph_eit_t eit) : m_eit(eit) {}

    /// Destroys the edge iterator
    ~EdgeIterator() {
        igraph_eit_destroy(&m_eit);
    }

    /********************/
    /* Instance methods */
    /********************/

    /// Returns a pointer to the encapsulated igraph_eit_t instance
    igraph_eit_t* c_eit() {
        return &m_eit;
    }

    /// Returns a pointer to the encapsulated igraph_eit_t instance (const)
    const igraph_eit_t* c_eit() const {
        return &m_eit;
    }

    /// Checks whether we are at the end of the iterator
    bool end() {
        return IGRAPH_EIT_END(m_eit);
    }

    /// Returns the index of the current edge of the iterator
    integer_t get() const {
        return IGRAPH_EIT_GET(m_eit);
    }

    /// Resets the iterator
    void reset() {
        IGRAPH_EIT_RESET(m_eit);
    }

    /// Returns the number of edges in the iterator
    long int size() {
        return IGRAPH_EIT_SIZE(m_eit);
    }

    /*************/
    /* Operators */
    /*************/

    /**
     * Increments the iterator, i.e. moves it to the next edge.
     *
     * Use this method if and only if the \ref end() method returns false.
     */
    EdgeIterator& operator++() {
        IGRAPH_EIT_NEXT(m_eit);
        return *this;
    }

    /**
     * Dereferences the iterator, i.e. returns the current Edge.
     */
    Edge operator*() const {
        return m_pGraph->edge(IGRAPH_EIT_GET(m_eit));
    }

    /*****************/
    /* Private stuff */
    /*****************/

private:
    /// Assignment operator (intentionally unimplemented)
    EdgeIterator& operator=(const EdgeIterator&);
};

}       // end of namespaces

#endif  // IGRAPHPP_EDGE_SELECTOR_H

