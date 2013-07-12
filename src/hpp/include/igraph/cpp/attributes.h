/* vim:set ts=4 sw=4 sts=4 et: */

#ifndef IGRAPHPP_ATTRIBUTES_H
#define IGRAPHPP_ATTRIBUTES_H

#include <map>
#include <vector>
#include <igraph/cpp/any.hpp>
#include <igraph/cpp/str_vector.h>
#include <igraph/cpp/types.h>

namespace igraph {

/// Typedef for attribute values
typedef any AttributeValue;

/// Typedef for a vector that stores attributes
typedef std::vector<AttributeValue> AttributeValueVector;

class AttributeHandlerImpl;

/// Attribute holder class for graphs
class AttributeHolder {
private:
    /// Typedef for the graph attribute storage
    typedef std::map<std::string, AttributeValue> GraphAttributeMap;

    /// Typedef for the vertex attribute storage
    typedef std::map<std::string, AttributeValueVector> VertexAttributeMap;

    /// Typedef for the edge attribute storage
    typedef std::map<std::string, AttributeValueVector> EdgeAttributeMap;

    /// Storage for the graph attributes
    GraphAttributeMap m_graphAttributes;
    
    /// Storage for the vertex attributes
    VertexAttributeMap m_vertexAttributes;

    /// Storage for the edge attributes
    EdgeAttributeMap m_edgeAttributes;

public:
    /// Returns a reference to the value of the given edge attribute
    AttributeValueVector& getEdgeAttributeReference(const std::string& attribute,
            integer_t ecount) {
        AttributeValueVector& result = m_edgeAttributes[attribute];
        if (result.size() < static_cast<size_t>(ecount)) {
            result.resize(ecount);
        }
        return result;
    }

    /// Returns a copy of the value of the given edge attribute
    AttributeValueVector getEdgeAttribute(const std::string& attribute) const;

    /// Returns the value of the given edge attribute for a given edge index
    AttributeValue getEdgeAttribute(const std::string& attribute,
        long int index) const;

    /// Stores the list of edge attributes in an StrVector
    void getEdgeAttributeList(StrVector& container) const;

    /// Returns a reference to the value of the given graph attribute
    AttributeValue& getGraphAttributeReference(const std::string& attribute) {
        return m_graphAttributes[attribute];
    }

    /// Returns a copy of the value of the given graph attribute
    AttributeValue getGraphAttribute(const std::string& attribute) const;

    /// Stores the list of graph attributes in an StrVector
    void getGraphAttributeList(StrVector& container) const;

    /// Returns a reference to the value of the given vertex attribute
    AttributeValueVector& getVertexAttributeReference(const std::string& attribute,
            integer_t vcount) {
        AttributeValueVector& result = m_vertexAttributes[attribute];
        if (result.size() < static_cast<size_t>(vcount)) {
            result.resize(vcount);
        }
        return result;
    }

    /// Returns a copy of the value of the given vertex attribute
    AttributeValueVector getVertexAttribute(const std::string& attribute) const;

    /// Returns the value of the given vertex attribute for a given vertex index
    AttributeValue getVertexAttribute(const std::string& attribute,
        long int index) const;

    /// Stores the list of vertex attributes in an StrVector
    void getVertexAttributeList(StrVector& container) const;

    /// Checks whether a given edge attribute exists
    bool hasEdgeAttribute(const std::string& attribute) const;

    /// Checks whether a given graph attribute exists
    bool hasGraphAttribute(const std::string& attribute) const;

    /// Checks whether a given vertex attribute exists
    bool hasVertexAttribute(const std::string& attribute) const;

    /// Sets the value of a given graph attribute
    void setGraphAttribute(const std::string& attribute, const AttributeValue& value);

    /// Sets the value of a given vertex attribute for all vertices
    void setVertexAttribute(const std::string& attribute, const AttributeValueVector& values);

    /// Sets the value of a given edge attribute for all edges
    void setEdgeAttribute(const std::string& attribute, const AttributeValueVector& values);

    friend class AttributeHandlerImpl;

private:
    /// Copies the graph attributes from another holder
    void copyGraphAttributesFrom(const AttributeHolder& other);

    /// Copies the vertex attributes from another holder
    void copyVertexAttributesFrom(const AttributeHolder& other);

    /// Copies the edge attributes from another holder
    void copyEdgeAttributesFrom(const AttributeHolder& other);
};

/// Attribute handler to be attached to igraph's attribute handler interface
struct AttributeHandler {
    /// Attaches the attribute handler to igraph
    static void attach();
};

}               // end of namespaces

#endif          // IGRAPHPP_ATTRIBUTES_H

