/* vim:set ts=4 sw=4 sts=4 et: */

#include <igraph/igraph_conversion.h>
#include <igraph/igraph_foreign.h>
#include <igraph/igraph_games.h>
#include <igraph/igraph_operators.h>
#include <igraph/igraph_structural.h>
#include <igraph/cpp/edge.h>
#include <igraph/cpp/edge_selector.h>
#include <igraph/cpp/graph.h>
#include <igraph/cpp/vertex.h>
#include <igraph/cpp/vertex_selector.h>
#include <memory>

namespace igraph {

/// Destructor
Graph::~Graph() {
    if (m_pGraph) {
        igraph_destroy(m_pGraph);
        delete m_pGraph;
    }
}

/******************/
/* Static methods */
/******************/

/********************/
/* Instance methods */
/********************/

void Graph::addEdge(igraph_integer_t source, igraph_integer_t target) {
    Vector edge(2);
    edge[0] = source; edge[1] = target;
    addEdges(edge);
}

void Graph::addEdges(const Vector& edges) {
    assert(m_pGraph);
    IGRAPH_TRY(igraph_add_edges(m_pGraph, edges.c_vector(), 0));
}

void Graph::addVertex() { addVertices(1); }

void Graph::addVertices(long int numVertices) {
    assert(m_pGraph);
    IGRAPH_TRY(igraph_add_vertices(m_pGraph, numVertices, 0));
}

bool Graph::areConnected(long int u, long int v) const {
    igraph_bool_t result;

    assert(m_pGraph);
    IGRAPH_TRY(igraph_are_connected(m_pGraph, u, v, &result));

    return result;
}

Vector Graph::degree(const VertexSelector& vids, NeighborMode mode, bool loops) const {
    Vector result;
    degree(&result, vids, mode, loops);
    return result;
}

void Graph::degree(Vector* result, const VertexSelector& vids,
                   NeighborMode mode, bool loops) const {
    assert(m_pGraph);
    IGRAPH_TRY(igraph_degree(m_pGraph, result->c_vector(), *vids.c_vs(),
                mode, loops));
}

void Graph::deleteEdges(const EdgeSelector& es) {
    assert(m_pGraph);
    IGRAPH_TRY(igraph_delete_edges(m_pGraph, *es.c_es()));
}

void Graph::edge(integer_t eid, integer_t* from, integer_t* to) const {
    assert(m_pGraph);
    IGRAPH_TRY(igraph_edge(m_pGraph, eid, from, to));
}

Edge Graph::edge(integer_t eid) {
    assert(m_pGraph);
    return Edge(this, eid);
}

AttributeValue Graph::getAttribute(const std::string& attribute) const {
    return getAttributeHolder()->getGraphAttribute(attribute);
}

void Graph::getEdgelist(Vector* result, bool bycol) const {
    assert(m_pGraph);
    IGRAPH_TRY(igraph_get_edgelist(m_pGraph, result->c_vector(), bycol));
}

Vector Graph::getEdgelist(bool bycol) const {
    Vector result;
    getEdgelist(&result, bycol);
    return result;
}

integer_t Graph::getEid(integer_t source, integer_t target,
        bool directed, bool error) const {
    integer_t eid;
    IGRAPH_TRY(igraph_get_eid(m_pGraph, &eid, source, target, directed, error));
    return eid;
}

bool Graph::isDirected() const {
    return igraph_is_directed(m_pGraph);
}

bool Graph::isSimple() const {
    bool_t result;
    IGRAPH_TRY(igraph_is_simple(m_pGraph, &result));
    return result;
}

bool Graph::hasAttribute(const std::string& attribute) const {
    return getAttributeHolder()->hasGraphAttribute(attribute);
}

void Graph::incident(Vector* result, long int vertex, NeighborMode mode) const {
    assert(m_pGraph);
    IGRAPH_TRY(igraph_incident(m_pGraph, result->c_vector(), vertex, mode));
}

Vector Graph::incident(long int vertex, NeighborMode mode) const {
    Vector result;
    incident(&result, vertex, mode);
    return result;
}

void Graph::neighbors(Vector* result, long int vertex, NeighborMode mode) const {
    assert(m_pGraph);
    IGRAPH_TRY(igraph_neighbors(m_pGraph, result->c_vector(), vertex, mode));
}

Vector Graph::neighbors(long int vertex, NeighborMode mode) const {
    Vector result;
    neighbors(&result, vertex, mode);
    return result;
}

void Graph::setAttribute(const std::string& attribute, const AttributeValue& value) {
    return getAttributeHolder()->setGraphAttribute(attribute, value);
}

void Graph::simplify(bool multiple, bool loops) {
    // TODO: last argument (attribute combination)
    assert(m_pGraph);
    IGRAPH_TRY(igraph_simplify(m_pGraph, multiple, loops, 0));
}

Vertex Graph::vertex(integer_t vid) {
    assert(m_pGraph);
    return Vertex(this, vid);
}

/*************/
/* Operators */
/*************/

Graph Graph::operator+(const Graph& other) const {
    std::auto_ptr<igraph_t> result(new igraph_t);
    IGRAPH_TRY(igraph_disjoint_union(result.get(), m_pGraph, other.m_pGraph));
    return Graph(result.release());
}

AttributeValue& Graph::operator[](const std::string& attribute) {
    return getAttributeHolder()->getGraphAttributeReference(attribute);
}

/***************************************************************************/

AttributeHolder* Graph::getAttributeHolder() {
    return static_cast<AttributeHolder*>(m_pGraph->attr);
}

const AttributeHolder* Graph::getAttributeHolder() const {
    return static_cast<AttributeHolder*>(m_pGraph->attr);
}

}         // end of namespaces

