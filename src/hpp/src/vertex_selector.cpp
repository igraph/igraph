/* vim:set ts=4 sw=4 sts=4 et: */

#include <igraph/cpp/graph.h>
#include <igraph/cpp/vertex_selector.h>

namespace igraph {

AttributeValueVector VertexSelector::getAttribute(const std::string& attribute) const {
    const AttributeHolder& attributeHolder = *m_pGraph->getAttributeHolder();

    if (isAll())
        return attributeHolder.getVertexAttribute(attribute);

    AttributeValueVector result;
    igraph_vit_t vit;
    IGRAPH_TRY(igraph_vit_create(m_pGraph->c_graph(), m_vs, &vit));
    while (!IGRAPH_VIT_END(vit)) {
        // TODO: vectorize this!
        result.push_back(attributeHolder.getVertexAttribute(attribute, (long int)IGRAPH_VIT_GET(vit)));
        IGRAPH_VIT_NEXT(vit);
    }
    igraph_vit_destroy(&vit);

    return result;
}

void VertexSelector::setAttribute(const std::string& attribute,
        const AttributeValueVector& values) const {
    /*
    if (isAll()) {
        m_pGraph->getAttributeHolder()->setVertexAttribute(attribute, values);
        return;
    }
    */

    AttributeValueVector& allAttrs =
        m_pGraph->getAttributeHolder()->getVertexAttributeReference(attribute,
                m_pGraph->vcount());

    igraph_vit_t vit;
    AttributeValueVector::const_iterator it = values.begin();
    AttributeValueVector::const_iterator endIt = values.end();
    IGRAPH_TRY(igraph_vit_create(m_pGraph->c_graph(), m_vs, &vit));
    while (!IGRAPH_VIT_END(vit)) {
        assert(it != endIt);
        allAttrs[(long int)IGRAPH_VIT_NEXT(vit)] = *it;
        it++;
        IGRAPH_VIT_NEXT(vit);
    }
    igraph_vit_destroy(&vit);
}

}           // end of namespaces
