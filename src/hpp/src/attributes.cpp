/* vim:set ts=4 sw=4 sts=4 et: */

#include <igraph/igraph_attributes.h>
#include <igraph/igraph_interface.h>
#include <igraph/igraph_iterators.h>
#include <igraph/cpp/attributes.h>
#include <igraph/cpp/error.h>
#include <igraph/cpp/ptr_vector.h>
#include <igraph/cpp/vector.h>

namespace igraph {

/***************************************************************************/

struct first_selector {
    template <typename T>
    typename T::first_type operator()(T pair) const {
        return pair.first;
    }
};

/***************************************************************************/

void AttributeHolder::copyGraphAttributesFrom(const AttributeHolder& other) {
    m_graphAttributes = other.m_graphAttributes;
}

void AttributeHolder::copyVertexAttributesFrom(const AttributeHolder& other) {
    m_vertexAttributes = other.m_vertexAttributes;
}

void AttributeHolder::copyEdgeAttributesFrom(const AttributeHolder& other) {
    m_edgeAttributes = other.m_edgeAttributes;
}

AttributeValueVector AttributeHolder::getEdgeAttribute(const std::string& attribute) const {
    EdgeAttributeMap::const_iterator it = m_edgeAttributes.find(attribute);
    if (it == m_edgeAttributes.end())
        return AttributeValueVector();
    return it->second;
}

AttributeValue AttributeHolder::getEdgeAttribute(const std::string& attribute,
        long int index) const {
    EdgeAttributeMap::const_iterator it = m_edgeAttributes.find(attribute);
    if (it == m_edgeAttributes.end())
        return AttributeValue();
    return it->second[index];
}

void AttributeHolder::getEdgeAttributeList(StrVector& container) const {
    EdgeAttributeMap::const_iterator it = m_edgeAttributes.begin();
    while (it != m_edgeAttributes.end()) {
        container.push_back(it->first);
        ++it;
    }
}

AttributeValue AttributeHolder::getGraphAttribute(const std::string& attribute) const {
    GraphAttributeMap::const_iterator it = m_graphAttributes.find(attribute);
    if (it == m_graphAttributes.end())
        return AttributeValue();
    return it->second;
}

void AttributeHolder::getGraphAttributeList(StrVector& container) const {
    GraphAttributeMap::const_iterator it = m_graphAttributes.begin();
    while (it != m_graphAttributes.end()) {
        container.push_back(it->first);
        ++it;
    }
}

AttributeValueVector AttributeHolder::getVertexAttribute(const std::string& attribute) const {
    VertexAttributeMap::const_iterator it = m_vertexAttributes.find(attribute);
    if (it == m_vertexAttributes.end())
        return AttributeValueVector();
    return it->second;
}

AttributeValue AttributeHolder::getVertexAttribute(const std::string& attribute,
        long int index) const {
    VertexAttributeMap::const_iterator it = m_vertexAttributes.find(attribute);
    if (it == m_vertexAttributes.end())
        return AttributeValue();
    return it->second[index];
}

void AttributeHolder::getVertexAttributeList(StrVector& container) const {
    VertexAttributeMap::const_iterator it = m_vertexAttributes.begin();
    while (it != m_vertexAttributes.end()) {
        container.push_back(it->first);
        ++it;
    }
}

bool AttributeHolder::hasEdgeAttribute(const std::string& attribute) const {
    return m_edgeAttributes.find(attribute) != m_edgeAttributes.end();
}

bool AttributeHolder::hasGraphAttribute(const std::string& attribute) const {
    return m_graphAttributes.find(attribute) != m_graphAttributes.end();
}

bool AttributeHolder::hasVertexAttribute(const std::string& attribute) const {
    return m_vertexAttributes.find(attribute) != m_vertexAttributes.end();
}

void AttributeHolder::setGraphAttribute(const std::string& attribute,
        const AttributeValue& value) {
    m_graphAttributes[attribute] = value;
}

void AttributeHolder::setVertexAttribute(const std::string& attribute,
        const AttributeValueVector& values) {
    m_vertexAttributes[attribute] = values;
}

void AttributeHolder::setEdgeAttribute(const std::string& attribute,
        const AttributeValueVector& values) {
    m_edgeAttributes[attribute] = values;
}

/***************************************************************************/

class AttributeHandlerImpl {
public:
    static int init(igraph_t *graph, igraph_vector_ptr_t *attr) {
        IGRAPHPP_TRY_NEW(graph->attr, AttributeHolder);
        return IGRAPH_SUCCESS;
    }

    static void destroy(igraph_t *graph) {
        if (graph->attr)
            delete static_cast<AttributeHolder*>(graph->attr);
        graph->attr = NULL;
    }

    static int copy(igraph_t *to, const igraph_t *from,
            igraph_bool_t ga, igraph_bool_t va, igraph_bool_t ea) {
        IGRAPHPP_TRY_NEW(to->attr, AttributeHolder);

        if (from->attr == NULL)
            return IGRAPH_SUCCESS;

        const AttributeHolder& from0 = *(static_cast<AttributeHolder*>(from->attr));
        AttributeHolder& to0 = *(static_cast<AttributeHolder*>(to->attr));

        if (ga)
            to0.copyGraphAttributesFrom(from0);
        if (va)
            to0.copyVertexAttributesFrom(from0);
        if (ea)
            to0.copyEdgeAttributesFrom(from0);

        return IGRAPH_SUCCESS;
    }

    /***********************************************************************/

    static igraph_attribute_type_t infer_type(const AttributeValue& value) {
        const std::type_info& type = value.type();

        if (type == typeid(bool)
                || type == typeid(short) || type == typeid(unsigned short)
                || type == typeid(int) || type == typeid(unsigned int)
                || type == typeid(long) || type == typeid(unsigned long)
                || type == typeid(long long) || type == typeid(unsigned long long)
                || type == typeid(float) || type == typeid(double)
                || type == typeid(long double))
            return IGRAPH_ATTRIBUTE_NUMERIC;

        if (type == typeid(const char*) || type == typeid(std::string))
            return IGRAPH_ATTRIBUTE_STRING;

        return IGRAPH_ATTRIBUTE_DEFAULT;
    }

    /***********************************************************************/

#define HANDLE_TYPE(T) {            \
    if (type == typeid(T)) {        \
        return *value.as<T>();      \
    }                               \
}

    static real_t as_numeric_attribute_value(const any& value) {
        const std::type_info& type = value.type();

        if (type == typeid(void))
            return 0;

        HANDLE_TYPE(bool);
        HANDLE_TYPE(short int);
        HANDLE_TYPE(unsigned short int);
        HANDLE_TYPE(int);
        HANDLE_TYPE(unsigned int);
        HANDLE_TYPE(long int);
        HANDLE_TYPE(unsigned long int);
        HANDLE_TYPE(long long int);
        HANDLE_TYPE(unsigned long long int);
        HANDLE_TYPE(float);
        HANDLE_TYPE(double);
        HANDLE_TYPE(long double);

        return 0;
    }

    static const char* as_string_attribute_value(const any& value) {
        static const char* empty = "";

        const std::type_info& type = value.type();
        if (type == typeid(void))
            return empty;

        HANDLE_TYPE(const char*);

        if (value.type() == typeid(std::string))
            return value.as<std::string>()->c_str();

        return empty;
    }

#undef HANDLE_TYPE

    /***********************************************************************/

    template <typename AttrMapType>
    static int add_vertices_edges_helper(AttrMapType& attrs, long int lastIndex,
            long int numNewEntries, igraph_vector_ptr_t* attr_vector_ptr) {
        // Extend the vertex attribute vectors by numNewEntries new elements
        for (AttributeHolder::VertexAttributeMap::iterator it = attrs.begin();
                it != attrs.end(); it++) {
            it->second.resize(lastIndex + numNewEntries);
        }

        // Do we have attributes with the newly added vertices? If not, return.
        if (attr_vector_ptr == 0)
            return IGRAPH_SUCCESS;

        PtrVector<igraph_attribute_record_t*> attrRecords(attr_vector_ptr);
        long int i, j;

        // For each attribute record...
        for (PtrVector<igraph_attribute_record_t*>::iterator it = attrRecords.begin();
                it != attrRecords.end(); it++) {
           igraph_attribute_record_t* record = (*it);

            // Do we have an attribute with this name? If not, add it.
            if (attrs.find(record->name) == attrs.end())
                attrs[record->name].resize(lastIndex + numNewEntries);

            AttributeValueVector& vec = attrs[record->name];
            switch (record->type) {
                case IGRAPH_ATTRIBUTE_NUMERIC:
                    {
                        Vector values(static_cast<igraph_vector_t*>(
                                    const_cast<void*>(record->value)
                        ), false);
                        for (i = 0, j = lastIndex; i < numNewEntries; i++, j++)
                            vec[j] = values[i];
                    }
                    break;

                case IGRAPH_ATTRIBUTE_STRING:
                    {
                        StrVector values(static_cast<igraph_strvector_t*>(
                                    const_cast<void*>(record->value)), false);
                        for (i = 0, j = lastIndex; i < numNewEntries; i++, j++)
                            vec[j] = std::string(values[i]);
                    }
                    break;

                default:
                    // Unsupported attribute type, just continue
                    // TODO: show a warning?
                    break;
            }
        }

        return IGRAPH_SUCCESS;
    }

    static int add_vertices(igraph_t *graph, long int nv, igraph_vector_ptr_t *attr) {
        AttributeHolder& holder = *(static_cast<AttributeHolder*>(graph->attr));
        AttributeHolder::VertexAttributeMap& attrs = holder.m_vertexAttributes;
        return add_vertices_edges_helper(attrs, igraph_vcount(graph)-nv, nv, attr);
    }

    /***********************************************************************/

    static int add_edges(igraph_t *graph, const igraph_vector_t *edges,
                  igraph_vector_ptr_t *attr) {
        AttributeHolder& holder = *(static_cast<AttributeHolder*>(graph->attr));
        AttributeHolder::EdgeAttributeMap& attrs = holder.m_edgeAttributes;
        long int ec = igraph_vector_size(edges) / 2;
        return add_vertices_edges_helper(attrs, igraph_ecount(graph)-ec, ec, attr);
    }

    /***********************************************************************/

    static int get_info(const igraph_t *graph,
            igraph_strvector_t *gnames, igraph_vector_t *gtypes,
            igraph_strvector_t *vnames, igraph_vector_t *vtypes,
            igraph_strvector_t *enames, igraph_vector_t *etypes) {
        AttributeHolder& holder = *(static_cast<AttributeHolder*>(graph->attr));
        igraph_attribute_type_t attr_type;
        integer_t vcount = igraph_vcount(graph);
        integer_t ecount = igraph_ecount(graph);

        if (gnames) {
            StrVector gnamesWrapper(gnames);
            holder.getGraphAttributeList(gnamesWrapper);

            // TODO: gtypes
        }

        if (vnames) {
            StrVector vnamesWrapper(vnames);
            Vector vtypesWrapper(vtypes);

            holder.getVertexAttributeList(vnamesWrapper);

            for (StrVector::const_iterator it = vnamesWrapper.begin();
                    it != vnamesWrapper.end(); ++it) {
                AttributeValueVector& ref = holder.getVertexAttributeReference(
                        *it, vcount);
                IGRAPH_CHECK(gettype_helper(graph, &attr_type, ref.begin(), ref.end()));
                vtypesWrapper.push_back(attr_type);
            }
        }

        if (enames) {
            StrVector enamesWrapper(enames);
            Vector etypesWrapper(etypes);

            holder.getEdgeAttributeList(enamesWrapper);

            for (StrVector::const_iterator it = enamesWrapper.begin();
                    it != enamesWrapper.end(); ++it) {
                AttributeValueVector& ref = holder.getEdgeAttributeReference(
                        *it, ecount);
                IGRAPH_CHECK(gettype_helper(graph, &attr_type, ref.begin(), ref.end()));
                etypesWrapper.push_back(attr_type);
            }
        }

        return IGRAPH_SUCCESS;
    }

    /***********************************************************************/

    static int get_numeric_graph_attr(const igraph_t *graph, const char *name,
            igraph_vector_t *value) {
        AttributeHolder& holder = *(static_cast<AttributeHolder*>(graph->attr));
        AttributeHolder::GraphAttributeMap& attrs = holder.m_graphAttributes;
        AttributeHolder::GraphAttributeMap::const_iterator it;

        it = attrs.find(name);
        if (it == attrs.end())
            IGRAPH_ERROR("Unknown attribute", IGRAPH_EINVAL);

        IGRAPH_CHECK(igraph_vector_resize(value, 1));
        VECTOR(*value)[0] = *it->second.as<igraph_real_t>();

        return IGRAPH_SUCCESS;
    }

    /***********************************************************************/

    static int get_numeric_vertex_attr(const igraph_t *graph, const char *name,
            igraph_vs_t vs, igraph_vector_t *value) {
        integer_t vcount = igraph_vcount(graph);
        integer_t i;
        AttributeHolder& holder = *(static_cast<AttributeHolder*>(graph->attr));

        if (!holder.hasVertexAttribute(name))
            IGRAPH_ERROR("Unknown attribute", IGRAPH_EINVAL);

        AttributeValueVector& ref = holder.getVertexAttributeReference(name, vcount);
        Vector valueWrapper(value);

        i = 0;

        if (igraph_vs_is_all(&vs)) {
            // Shortcut
            valueWrapper.resize(vcount);
            for (AttributeValueVector::const_iterator it = ref.begin(); it != ref.end(); ++it, ++i) {
                valueWrapper[i] = as_numeric_attribute_value(*it);
            }
        } else {
            igraph_vit_t vit;
            
            valueWrapper.resize(0);

            IGRAPH_CHECK(igraph_vit_create(graph, vs, &vit));
            IGRAPH_FINALLY(igraph_vit_destroy, &vit);

            while (!IGRAPH_VIT_END(vit)) {
                integer_t vid = IGRAPH_VIT_GET(vit);
                valueWrapper.push_back(as_numeric_attribute_value(ref[vid]));
                IGRAPH_VIT_NEXT(vit);
            }

            igraph_vit_destroy(&vit);
            IGRAPH_FINALLY_CLEAN(1);
        }

        return IGRAPH_SUCCESS;
    }

    /***********************************************************************/

    static int get_numeric_edge_attr(const igraph_t *graph, const char *name,
            igraph_es_t es, igraph_vector_t *value) {
        integer_t ecount = igraph_ecount(graph);
        integer_t i;
        AttributeHolder& holder = *(static_cast<AttributeHolder*>(graph->attr));

        if (!holder.hasEdgeAttribute(name))
            IGRAPH_ERROR("Unknown attribute", IGRAPH_EINVAL);

        AttributeValueVector& ref = holder.getEdgeAttributeReference(name, ecount);
        Vector valueWrapper(value);

        i = 0;

        if (igraph_es_is_all(&es)) {
            // Shortcut
            valueWrapper.resize(ecount);
            for (AttributeValueVector::const_iterator it = ref.begin(); it != ref.end(); ++it, ++i) {
                valueWrapper[i] = as_numeric_attribute_value(*it);
            }
        } else {
            igraph_eit_t eit;

            valueWrapper.resize(0);

            IGRAPH_CHECK(igraph_eit_create(graph, es, &eit));
            IGRAPH_FINALLY(igraph_eit_destroy, &eit);

            while (!IGRAPH_EIT_END(eit)) {
                integer_t eid = IGRAPH_EIT_GET(eit);
                valueWrapper.push_back(as_numeric_attribute_value(ref[eid]));
                IGRAPH_EIT_NEXT(eit);
            }

            igraph_eit_destroy(&eit);
            IGRAPH_FINALLY_CLEAN(1);
        }

        return IGRAPH_SUCCESS;
    }

    /***********************************************************************/

    static int get_string_graph_attr(const igraph_t *graph, const char *name,
            igraph_strvector_t *value) {
        AttributeHolder& holder = *(static_cast<AttributeHolder*>(graph->attr));
        AttributeHolder::GraphAttributeMap& attrs = holder.m_graphAttributes;
        AttributeHolder::GraphAttributeMap::const_iterator it;

        it = attrs.find(name);
        if (it == attrs.end())
            IGRAPH_ERROR("Unknown attribute", IGRAPH_EINVAL);

        IGRAPH_CHECK(igraph_strvector_resize(value, 1));
        IGRAPH_CHECK(igraph_strvector_set(value, 0, *it->second.as<const char*>()));

        return IGRAPH_SUCCESS;
    }

    /***********************************************************************/

    static int get_string_vertex_attr(const igraph_t *graph, const char *name,
            igraph_vs_t vs, igraph_strvector_t *value) {
        integer_t vcount = igraph_vcount(graph);
        integer_t i;
        AttributeHolder& holder = *(static_cast<AttributeHolder*>(graph->attr));

        if (!holder.hasVertexAttribute(name))
            IGRAPH_ERROR("Unknown attribute", IGRAPH_EINVAL);

        AttributeValueVector& ref = holder.getVertexAttributeReference(name, vcount);
        StrVector valueWrapper(value);

        i = 0;

        if (igraph_vs_is_all(&vs)) {
            // Shortcut
            valueWrapper.resize(vcount);
            for (AttributeValueVector::const_iterator it = ref.begin(); it != ref.end(); ++it, ++i) {
                igraph_strvector_set(value, i, *(it->as<const char*>()));
            }
        } else {
            igraph_vit_t vit;
            
            valueWrapper.resize(0);

            IGRAPH_CHECK(igraph_vit_create(graph, vs, &vit));
            IGRAPH_FINALLY(igraph_vit_destroy, &vit);

            while (!IGRAPH_VIT_END(vit)) {
                integer_t vid = IGRAPH_VIT_GET(vit);
                valueWrapper.push_back(as_string_attribute_value(ref[vid]));
                IGRAPH_VIT_NEXT(vit);
            }

            igraph_vit_destroy(&vit);
            IGRAPH_FINALLY_CLEAN(1);
        }

        return IGRAPH_SUCCESS;
    }

    /***********************************************************************/

    static int get_string_edge_attr(const igraph_t *graph, const char *name,
            igraph_es_t es, igraph_strvector_t *value) {
        integer_t ecount = igraph_ecount(graph);
        integer_t i;
        AttributeHolder& holder = *(static_cast<AttributeHolder*>(graph->attr));

        if (!holder.hasEdgeAttribute(name))
            IGRAPH_ERROR("Unknown attribute", IGRAPH_EINVAL);

        AttributeValueVector& ref = holder.getEdgeAttributeReference(name, ecount);
        StrVector valueWrapper(value);

        i = 0;

        if (igraph_es_is_all(&es)) {
            // Shortcut
            valueWrapper.resize(ecount);
            for (AttributeValueVector::const_iterator it = ref.begin(); it != ref.end(); ++it, ++i) {
                igraph_strvector_set(value, i, *(it->as<const char*>()));
            }
        } else {
            igraph_eit_t eit;

            valueWrapper.resize(0);

            IGRAPH_CHECK(igraph_eit_create(graph, es, &eit));
            IGRAPH_FINALLY(igraph_eit_destroy, &eit);

            while (!IGRAPH_EIT_END(eit)) {
                integer_t eid = IGRAPH_EIT_GET(eit);
                valueWrapper.push_back(as_string_attribute_value(ref[eid]));
                IGRAPH_EIT_NEXT(eit);
            }

            igraph_eit_destroy(&eit);
            IGRAPH_FINALLY_CLEAN(1);
        }

        return IGRAPH_SUCCESS;
    }

    /***********************************************************************/

    template <typename It>
    static int gettype_helper(const igraph_t *graph, igraph_attribute_type_t *type,
            const It& begin, const It& end) {
        igraph_attribute_type_t result = IGRAPH_ATTRIBUTE_DEFAULT;

        /* Infer the attribute type by evaluating the types of the items */
        for (It it = begin; it != end; ++it) {
            if (it->type() == typeid(void))
                continue;

            igraph_attribute_type_t current_type = infer_type(*it);

            if (result == IGRAPH_ATTRIBUTE_DEFAULT)
                result = current_type;
            else if (result != current_type) {
                // Mixed attributes are treated as strings
                result = IGRAPH_ATTRIBUTE_STRING;
                break;
            }
        }

        *type = result;
        return IGRAPH_SUCCESS;
    }

    template <typename T>
    static int gettype_helper(const igraph_t *graph, igraph_attribute_type_t *type,
            const typename T::const_iterator& it) {
        const typename T::mapped_type& values = it->second;
        return gettype_helper(graph, type, values.begin(), values.end());
    }

    template <typename T>
    static int gettype_helper(const igraph_t *graph, igraph_attribute_type_t *type,
            const char* name, const T& map) {
        typename T::const_iterator it = map.find(name);
        if (it == map.end())
            IGRAPH_ERROR("Unknown attribute", IGRAPH_EINVAL);

        return gettype_helper<T>(graph, type, it);
    }

    static int gettype(const igraph_t *graph, igraph_attribute_type_t *type,
            igraph_attribute_elemtype_t elemtype, const char *name) {
        AttributeHolder& holder = *(static_cast<AttributeHolder*>(graph->attr));
        switch (elemtype) {
        case IGRAPH_ATTRIBUTE_GRAPH:
            *type = infer_type(holder.getGraphAttribute(name));
            return IGRAPH_SUCCESS;

        case IGRAPH_ATTRIBUTE_VERTEX:
            return gettype_helper(graph, type, name, holder.m_vertexAttributes);

        case IGRAPH_ATTRIBUTE_EDGE:
            return gettype_helper(graph, type, name, holder.m_edgeAttributes);

        default:
            IGRAPH_ERROR("Unknown attribute element type", IGRAPH_EINVAL);
        }
    }

    /***********************************************************************/

    template <typename T>
    static igraph_bool_t has_attr_helper(const T& map, const char* name) {
        return map.find(name) != map.end();
    }

    static igraph_bool_t has_attr(const igraph_t *graph,
            igraph_attribute_elemtype_t type, const char *name) {
        AttributeHolder& holder = *(static_cast<AttributeHolder*>(graph->attr));

        switch (type) {
            case IGRAPH_ATTRIBUTE_GRAPH:
                return has_attr_helper(holder.m_graphAttributes, name);

            case IGRAPH_ATTRIBUTE_VERTEX:
                return has_attr_helper(holder.m_vertexAttributes, name);

            case IGRAPH_ATTRIBUTE_EDGE:
                return has_attr_helper(holder.m_edgeAttributes, name);

            default:
                IGRAPH_ERROR("Unknown attribute element type", IGRAPH_EINVAL);
        }
    }
};


/***************************************************************************/

static igraph_attribute_table_t cpp_attribute_handler = {
    &AttributeHandlerImpl::init,
    &AttributeHandlerImpl::destroy,
    &AttributeHandlerImpl::copy,
    &AttributeHandlerImpl::add_vertices,
    /* permute_vertices = */ 0,
    /* combine_vertices = */ 0,
    &AttributeHandlerImpl::add_edges,
    /* permute_edges = */ 0,
    /* combine_edges = */ 0,
    &AttributeHandlerImpl::get_info,
    &AttributeHandlerImpl::has_attr,
    &AttributeHandlerImpl::gettype,
    &AttributeHandlerImpl::get_numeric_graph_attr,
    &AttributeHandlerImpl::get_string_graph_attr,
		/* get_bool_graph_attr = */ 0,
    &AttributeHandlerImpl::get_numeric_vertex_attr,
    &AttributeHandlerImpl::get_string_vertex_attr,
		/* get_bool_vertex_attr = */ 0,
    &AttributeHandlerImpl::get_numeric_edge_attr,
    &AttributeHandlerImpl::get_string_edge_attr,
		/* get_bool_edge_attr = */ 0
};

void AttributeHandler::attach() {
    igraph_i_set_attribute_table(&cpp_attribute_handler);
}

}            // end of namespace

