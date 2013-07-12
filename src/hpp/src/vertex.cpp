/* vim:set ts=4 sw=4 sts=4 et: */

#include <igraph/cpp/graph.h>
#include <igraph/cpp/types.h>
#include <igraph/cpp/vertex.h>

namespace igraph {

AttributeValue Vertex::getAttribute(const std::string& attribute) const {
    const AttributeHolder& attributeHolder = *m_pGraph->getAttributeHolder();
    return attributeHolder.getVertexAttribute(attribute, m_index);
}

AttributeValue Vertex::getAttribute(const std::string& attribute,
        const AttributeValue& defaultValue) const {
    const AttributeHolder& attributeHolder = *m_pGraph->getAttributeHolder();
    if (!attributeHolder.hasVertexAttribute(attribute))
        return defaultValue;
    return attributeHolder.getVertexAttribute(attribute, m_index);
}

void Vertex::setAttribute(const std::string& attribute, const AttributeValue& value) {
    AttributeHolder& attributeHolder = *m_pGraph->getAttributeHolder();
    AttributeValueVector& vec = attributeHolder.getVertexAttributeReference(attribute,
            m_pGraph->vcount());
    vec[m_index] = value;
}

}         // end of namespace

