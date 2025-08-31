/*
   igraph library.
   Copyright (C) 2005-2012  Gabor Csardi <csardi.gabor@gmail.com>
   334 Harvard street, Cambridge, MA 02139 USA

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
   02110-1301 USA

*/

#include "igraph_attributes.h"
#include "igraph_memory.h"

#include "graph/attributes.h"
#include "internal/hacks.h" /* strdup */

#include <string.h>
#include <stdarg.h>

/**
 * \section about_attributes
 *
 * <para>Attributes are numbers, boolean values or strings associated with
 * the vertices or edges of a graph, or with the graph itself. E.g. you may
 * label vertices with symbolic names or attach numeric weights to the edges
 * of a graph. In addition to these three basic types, a custom object
 * type is supported as well.</para>
 *
 * <para>igraph attributes are designed to be flexible and extensible.
 * In igraph attributes are implemented via an interface abstraction:
 * any type implementing the functions in the interface can be used
 * for storing vertex, edge and graph attributes. This means that
 * different attribute implementations can be used together with
 * igraph. This is reasonable: if igraph is used from Python attributes can be
 * of any Python type, from R all R types are allowed. There is also an
 * experimental attribute implementation to be used when programming
 * in C, but by default it is currently turned off.</para>
 *
 * <para>First we briefly look over how attribute handlers can be
 * implemented. This is not something a user does every day. It is
 * rather typically the job of the high level interface writers. (But
 * it is possible to write an interface without implementing
 * attributes.) Then we show the experimental C attribute handler.</para>
 */

/**
 * \section about_attribute_table
 * <para>It is possible to attach an attribute handling
 * interface to \a igraph. This is simply a table of functions, of
 * type \ref igraph_attribute_table_t. These functions are invoked to
 * notify the attribute handling code about the structural changes in
 * a graph. See the documentation of this type for details.</para>
 *
 * <para>By default there is no attribute interface attached to \a igraph.
 * To attach one, call \ref igraph_set_attribute_table with your new
 * table. This is normally done on program startup, and is kept untouched
 * for the program's lifetime. It must be done before any graph object
 * is created, as graphs created with a given attribute handler
 * cannot be manipulated while a different attribute handler is
 * active.</para>
 */

/**
 * \section about_attribute_record
 *
 * <para>Functions in the attribute handler interface may refer to
 * \em "attribute records" or \em "attribute record lists". An attribute record
 * is simply a triplet consisting of an attribute name, an attribute type and
 * a vector containing the values of the attribute. Attribute record lists are
 * typed containers that contain a sequence of attribute records. Attribute
 * record lists own the attribute records that they contain, and similarly,
 * attribute records own the vectors contained in them. Destroying an attribute
 * record destroys the vector of values inside it, and destroying an attribute
 * record list destroys all attribute records in the list.</para>
 */

/**
 * \section about_attribute_combination
 *
 * <para>Several graph operations may collapse multiple vertices or edges into
 * a single one. Attribute combination lists are used to indicate to the attribute
 * handler how to combine the attributes of the original vertices or edges and
 * how to derive the final attribute value that is to be assigned to the collapsed
 * vertex or edge. For example, \ref igraph_simplify() removes loops and combines
 * multiple edges into a single one; in case of a graph with an edge attribute
 * named \c weight the attribute combination list can tell the attribute handler
 * whether the weight of a collapsed edge should be the sum, the mean or some other
 * function of the weights of the original edges that were collapsed into one.</para>
 *
 * <para>One attribute combination list may contain several attribute combination
 * records, one for each vertex or edge attribute that is to be handled during the
 * operation.</para>
 */

static void igraph_i_attribute_record_set_type(
    igraph_attribute_record_t *attr, igraph_attribute_type_t type, void *ptr
);
static void igraph_i_attribute_record_destroy_values(igraph_attribute_record_t *attr);

/**
 * \function igraph_attribute_record_init
 * \brief Initializes an attribute record with a given name and type.
 *
 * \param attr  the attribute record to initialize
 * \param name  name of the attribute
 * \param type  type of the attribute
 * \return Error code:
 *         \c IGRAPH_ENOMEM if there is not enough memory.
 *
 * Time complexity: O(1).
 */
igraph_error_t igraph_attribute_record_init(
    igraph_attribute_record_t *attr, const char* name, igraph_attribute_type_t type
) {
    attr->name = NULL;
    attr->type = IGRAPH_ATTRIBUTE_UNSPECIFIED;
    attr->value.as_raw = NULL;
    attr->default_value.string = NULL;

    IGRAPH_CHECK(igraph_attribute_record_set_name(attr, name));
    IGRAPH_CHECK(igraph_attribute_record_set_type(attr, type));

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_attribute_record_init_copy
 * \brief Initializes an attribute record by copying another record.
 *
 * </para><para>
 * Copies made by this function are deep copies: a full copy of the value
 * vector contained in the record is placed in the new record so they become
 * independent of each other.
 *
 * \param to    the attribute record to initialize
 * \param from  the attribute record to copy data from
 * \return Error code:
 *         \c IGRAPH_ENOMEM if there is not enough memory.
 *
 * Time complexity: operating system dependent, usually O(n), where n is the
 * size of the value vector in the attribute record.
 */
igraph_error_t igraph_attribute_record_init_copy(
    igraph_attribute_record_t *to, const igraph_attribute_record_t *from
) {
    IGRAPH_CHECK(igraph_attribute_record_init(to, from->name, from->type));

    switch (from->type) {
        case IGRAPH_ATTRIBUTE_NUMERIC:
            IGRAPH_CHECK(igraph_vector_update(to->value.as_vector, from->value.as_vector));
            break;

        case IGRAPH_ATTRIBUTE_STRING:
            IGRAPH_CHECK(igraph_strvector_update(to->value.as_strvector, from->value.as_strvector));
            break;

        case IGRAPH_ATTRIBUTE_BOOLEAN:
            IGRAPH_CHECK(igraph_vector_bool_update(to->value.as_vector_bool, from->value.as_vector_bool));
            break;

        case IGRAPH_ATTRIBUTE_UNSPECIFIED:
            break;

        default:
            break;
    }

    return IGRAPH_SUCCESS;
}

static void igraph_i_attribute_record_destroy_values(igraph_attribute_record_t *attr) {
    IGRAPH_ASSERT(attr != NULL);

    if (attr->value.as_raw) {
        switch (attr->type) {
            case IGRAPH_ATTRIBUTE_NUMERIC:
                igraph_vector_destroy(attr->value.as_vector);
                break;

            case IGRAPH_ATTRIBUTE_STRING:
                igraph_strvector_destroy(attr->value.as_strvector);
                break;

            case IGRAPH_ATTRIBUTE_BOOLEAN:
                igraph_vector_bool_destroy(attr->value.as_vector_bool);
                break;

            default:
                break;
        }

        igraph_free(attr->value.as_raw);
        attr->value.as_raw = NULL;
    }

    switch (attr->type) {
        case IGRAPH_ATTRIBUTE_NUMERIC:
            attr->default_value.numeric = 0;
            break;

        case IGRAPH_ATTRIBUTE_STRING:
            if (attr->default_value.string) {
                igraph_free(attr->default_value.string);
                attr->default_value.string = NULL;
            }
            break;

        case IGRAPH_ATTRIBUTE_BOOLEAN:
            attr->default_value.boolean = 0;
            break;

        default:
            break;
    }

    attr->type = IGRAPH_ATTRIBUTE_UNSPECIFIED;
}

/**
 * \function igraph_attribute_record_destroy
 * \brief Destroys an attribute record.
 *
 * \param attr  the previously initialized attribute record to destroy.
 *
 * Time complexity: operating system dependent.
 */
void igraph_attribute_record_destroy(igraph_attribute_record_t *attr) {
    igraph_i_attribute_record_destroy_values(attr);

    if (attr->name) {
        igraph_free(attr->name);
        attr->name = NULL;
    }
}

/**
 * \function igraph_attribute_record_check_type
 * \brief Checks whether the type of the attribute record is equal to an expected type.
 *
 * \param attr  the attribute record to test
 * \param type  the expected type of the attribute record
 * \return Error code:
 *         \c IGRAPH_EINVAL if the type of the attribute record is not equal to
 *         the expected type
 *
 * Time complexity: O(1).
 */
igraph_error_t igraph_attribute_record_check_type(
    const igraph_attribute_record_t *attr, igraph_attribute_type_t type
) {
    if (type != attr->type) {
        switch (type) {
            case IGRAPH_ATTRIBUTE_STRING:
                IGRAPH_ERROR("String attribute expected.", IGRAPH_EINVAL);
                break;
            case IGRAPH_ATTRIBUTE_NUMERIC:
                IGRAPH_ERROR("Numeric attribute expected.", IGRAPH_EINVAL);
                break;
            case IGRAPH_ATTRIBUTE_BOOLEAN:
                IGRAPH_ERROR("Boolean attribute expected.", IGRAPH_EINVAL);
                break;
            case IGRAPH_ATTRIBUTE_OBJECT:
                IGRAPH_ERROR("Object attribute expected.", IGRAPH_EINVAL);
                break;
            default:
                IGRAPH_ERROR("Attribute with unknown type expected.", IGRAPH_EINVAL);
                break;
        }
    }

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_attribute_record_size
 * \brief Returns the size of the value vector in an attribute record.
 *
 * \param attr      the attribute record to query
 * \return the number of elements in the value vector of the attribute record
 */
igraph_int_t igraph_attribute_record_size(const igraph_attribute_record_t *attr) {
    IGRAPH_ASSERT(attr != NULL);

    switch (attr->type) {
        case IGRAPH_ATTRIBUTE_NUMERIC:
            return igraph_vector_size(attr->value.as_vector);

        case IGRAPH_ATTRIBUTE_STRING:
            return igraph_strvector_size(attr->value.as_strvector);

        case IGRAPH_ATTRIBUTE_BOOLEAN:
            return igraph_vector_bool_size(attr->value.as_vector_bool);

        case IGRAPH_ATTRIBUTE_UNSPECIFIED:
            return 0;

        default:
            IGRAPH_ERRORF("Unsupported attribute type: %d", IGRAPH_EINVAL, (int) attr->type);
    }
}

/**
 * \function igraph_attribute_record_resize
 * \brief Resizes the value vector in an attribute record.
 *
 * </para><para>When the value vector is shorter than the desired length, it
 * will be expanded with \c IGRAPH_NAN for numeric vectors, \c false for Boolean
 * vectors and empty strings for string vectors.
 *
 * \param attr      the attribute record to update
 * \param new_size  the new size of the value vector
 * \return Error code:
 *         \c IGRAPH_ENOMEM if there is not enough memory.
 *         \c IGRAPH_EINVAL if the type of the attribute record is not specified yet.
 */
igraph_error_t igraph_attribute_record_resize(
    igraph_attribute_record_t *attr, igraph_int_t new_size
) {
    igraph_int_t i;
    igraph_vector_t *vec;
    igraph_vector_bool_t *log;
    igraph_strvector_t *str;

    IGRAPH_ASSERT(attr != NULL);

    switch (attr->type) {

        case IGRAPH_ATTRIBUTE_NUMERIC:
            vec = attr->value.as_vector;
            i = igraph_vector_size(vec);
            IGRAPH_CHECK(igraph_vector_resize(vec, new_size));
            while (i < new_size) {
                VECTOR(*vec)[i++] = attr->default_value.numeric;
            }
            break;

        case IGRAPH_ATTRIBUTE_BOOLEAN:
            log = attr->value.as_vector_bool;
            i = igraph_vector_bool_size(log);
            IGRAPH_CHECK(igraph_vector_bool_resize(log, new_size));
            while (i < new_size) {
                VECTOR(*log)[i++] = attr->default_value.boolean;
            }
            break;

        case IGRAPH_ATTRIBUTE_STRING:
            str = attr->value.as_strvector;
            if (attr->default_value.string == 0 || (*attr->default_value.string == 0)) {
                IGRAPH_CHECK(igraph_strvector_resize(str, new_size));
            } else {
                i = igraph_strvector_size(str);
                IGRAPH_CHECK(igraph_strvector_resize(str, new_size));
                while (i < new_size) {
                    IGRAPH_CHECK(igraph_strvector_set(str, i++, attr->default_value.string));
                }
            }
            break;

        case IGRAPH_ATTRIBUTE_UNSPECIFIED:
            IGRAPH_ERROR("Attribute record has no type yet.", IGRAPH_EINVAL);
            break;

        default:
            IGRAPH_ERRORF("Unsupported attribute type: %d", IGRAPH_EINVAL, (int) attr->type);
    }

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_attribute_record_set_default_numeric
 * \brief Sets the default value of the attribute to the given number.
 *
 * </para><para>
 * This function must be called for numeric attribute records only. When not
 * specified, the default value of numeric attributes is NaN.
 *
 * \param attr   the attribute record to update
 * \param value  the new default value
 * \return Error code:
 *         \c IGRAPH_EINVAL if the attribute record has a non-numeric type
 */
igraph_error_t igraph_attribute_record_set_default_numeric(
    igraph_attribute_record_t *attr, igraph_real_t value
) {
    if (attr->type != IGRAPH_ATTRIBUTE_NUMERIC) {
        return IGRAPH_EINVAL;
    }

    attr->default_value.numeric = value;
    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_attribute_record_set_default_boolean
 * \brief Sets the default value of the attribute to the given logical value.
 *
 * </para><para>
 * This function must be called for Boolean attribute records only. When not
 * specified, the default value of Boolean attributes is \c false.
 *
 * \param attr   the attribute record to update
 * \param value  the new default value
 * \return Error code:
 *         \c IGRAPH_EINVAL if the attribute record is not of Boolean type
 */
IGRAPH_EXPORT igraph_error_t igraph_attribute_record_set_default_boolean(
    igraph_attribute_record_t *attr, igraph_bool_t value
) {
    if (attr->type != IGRAPH_ATTRIBUTE_BOOLEAN) {
        return IGRAPH_EINVAL;
    }

    attr->default_value.boolean = value;
    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_attribute_record_set_default_string
 * \brief Sets the default value of the attribute to the given string.
 *
 * </para><para>
 * This function must be called for string attribute records only. When not
 * specified, the default value of string attributes is an empty string.
 *
 * \param attr   the attribute record to update
 * \param value  the new default value. \c NULL means an empty string.
 *
 * \return Error code:
 *         \c IGRAPH_ENOMEM if there is not enough memory
 *         \c IGRAPH_EINVAL if the attribute record is not of string type
 */
IGRAPH_EXPORT igraph_error_t igraph_attribute_record_set_default_string(
    igraph_attribute_record_t *attr, const char* value
) {
    char* copy;

    if (attr->type != IGRAPH_ATTRIBUTE_STRING) {
        return IGRAPH_EINVAL;
    }

    if (value && (*value != 0)) {
        copy = strdup(value);
        IGRAPH_CHECK_OOM(copy, "Insufficient memory to duplicate default value.");
    } else {
        copy = NULL;
    }

    if (attr->default_value.string) {
        igraph_free(attr->default_value.string);
    }
    attr->default_value.string = copy;

    return IGRAPH_SUCCESS;
}


/**
 * \function igraph_attribute_record_set_name
 * \brief Sets the attribute name in an attribute record.
 *
 * \param attr  the attribute record to update
 * \param name  the new name
 * \return Error code:
 *         \c IGRAPH_ENOMEM if there is not enough memory.
 */
igraph_error_t igraph_attribute_record_set_name(
    igraph_attribute_record_t *attr, const char* name
) {
    char* new_name;

    IGRAPH_ASSERT(attr != NULL);

    if (name != NULL) {
        new_name = strdup(name);
        IGRAPH_CHECK_OOM(new_name, "Insufficient memory for allocating attribute name.");
    } else {
        new_name = NULL;
    }

    if (attr->name) {
        igraph_free(attr->name);
    }

    attr->name = new_name;

    return IGRAPH_SUCCESS;
}

static void igraph_i_attribute_record_set_type(
    igraph_attribute_record_t *attr, igraph_attribute_type_t type, void *ptr
) {
    bool type_changed = attr->type != type;

    if (type_changed || attr->value.as_raw != ptr) {
        igraph_i_attribute_record_destroy_values(attr);
        attr->type = type;
        attr->value.as_raw = ptr;
    }

    if (type_changed && type == IGRAPH_ATTRIBUTE_NUMERIC) {
        IGRAPH_ASSERT(
            igraph_attribute_record_set_default_numeric(attr, IGRAPH_NAN) == IGRAPH_SUCCESS
        );
    }
}

/**
 * \function igraph_attribute_record_set_type
 * \brief Sets the type of an attribute record.
 *
 * </para><para>
 * When the new type being set is different from the old type, any values already
 * stored in the attribute record will be destroyed and a new, empty attribute
 * value vector will be allocated. When the new type is the same as the old
 * type, this function is a no-op.
 *
 * \param attr  the attribute record to update
 * \param type  the new type
 * \return Error code:
 *         \c IGRAPH_ENOMEM if there is not enough memory.
 */
igraph_error_t igraph_attribute_record_set_type(
    igraph_attribute_record_t *attr, igraph_attribute_type_t type
) {
    void* ptr;

    IGRAPH_ASSERT(attr != NULL);

    if (attr->type != type) {
        switch (type) {
            case IGRAPH_ATTRIBUTE_NUMERIC: {
                igraph_vector_t* vec = IGRAPH_CALLOC(1, igraph_vector_t);
                IGRAPH_CHECK_OOM(vec, "Insufficient memory for attribute record.");
                IGRAPH_FINALLY(igraph_free, vec);
                IGRAPH_VECTOR_INIT_FINALLY(vec, 0);
                ptr = vec;
            }
            break;

            case IGRAPH_ATTRIBUTE_STRING: {
                igraph_strvector_t* strvec = IGRAPH_CALLOC(1, igraph_strvector_t);
                IGRAPH_CHECK_OOM(strvec, "Insufficient memory for attribute record.");
                IGRAPH_FINALLY(igraph_free, strvec);
                IGRAPH_CHECK(igraph_strvector_init(strvec, 0));
                IGRAPH_FINALLY(igraph_strvector_destroy, strvec);
                ptr = strvec;
            }
            break;

            case IGRAPH_ATTRIBUTE_BOOLEAN: {
                igraph_vector_bool_t* boolvec = IGRAPH_CALLOC(1, igraph_vector_bool_t);
                IGRAPH_CHECK_OOM(boolvec, "Insufficient memory for attribute record.");
                IGRAPH_FINALLY(igraph_free, boolvec);
                IGRAPH_VECTOR_BOOL_INIT_FINALLY(boolvec, 0);
                ptr = boolvec;
            }
            break;

            default:
                IGRAPH_ERRORF("Unsupported attribute type: %d", IGRAPH_EINVAL, (int) type);
        }

        igraph_i_attribute_record_set_type(attr, type, ptr);
        IGRAPH_FINALLY_CLEAN(2);
    }

    return IGRAPH_SUCCESS;
}

#define ATTRIBUTE_RECORD_LIST
#define BASE_ATTRIBUTE_RECORD
#define CUSTOM_INIT_DESTROY
#include "igraph_pmt.h"
#include "../core/typed_list.pmt"
#include "igraph_pmt_off.h"
#undef CUSTOM_INIT_DESTROY
#undef BASE_ATTRIBUTE_RECORD
#undef ATTRIBUTE_RECORD_LIST

static igraph_error_t igraph_i_attribute_record_list_init_item(
    const igraph_attribute_record_list_t* list, igraph_attribute_record_t* item
) {
    IGRAPH_UNUSED(list);
    return igraph_attribute_record_init(item, NULL, IGRAPH_ATTRIBUTE_UNSPECIFIED);
}

static igraph_error_t igraph_i_attribute_record_list_copy_item(
    igraph_attribute_record_t* dest, const igraph_attribute_record_t* source
) {
    return igraph_attribute_record_init_copy(dest, source);
}

static void igraph_i_attribute_record_list_destroy_item(igraph_attribute_record_t* item) {
    igraph_attribute_record_destroy(item);
}

/* Should you ever want to have a thread-local attribute handler table, prepend
 * IGRAPH_THREAD_LOCAL to the following declaration and #include "config.h". */
igraph_attribute_table_t *igraph_i_attribute_table = NULL;

igraph_error_t igraph_i_attribute_init(
    igraph_t *graph, const igraph_attribute_record_list_t *attr
) {
    graph->attr = NULL;
    if (igraph_i_attribute_table) {
        IGRAPH_CHECK(igraph_i_attribute_table->init(graph, attr));
        if (graph->attr == NULL) {
            IGRAPH_ERROR("Attribute handler did not initialize attr pointer", IGRAPH_FAILURE);
        }
    }
    return IGRAPH_SUCCESS;
}

void igraph_i_attribute_destroy(igraph_t *graph) {
    if (graph->attr && igraph_i_attribute_table) {
        igraph_i_attribute_table->destroy(graph);
    }
    graph->attr = NULL;
}

igraph_error_t igraph_i_attribute_copy(igraph_t *to, const igraph_t *from, igraph_bool_t ga,
                            igraph_bool_t va, igraph_bool_t ea) {
    igraph_i_attribute_destroy(to);
    if (from->attr && igraph_i_attribute_table) {
        IGRAPH_CHECK(igraph_i_attribute_table->copy(to, from, ga, va, ea));
        if (to->attr == NULL) {
            IGRAPH_ERROR("Attribute handler did not initialize attr pointer", IGRAPH_FAILURE);
        }
    }
    return IGRAPH_SUCCESS;
}

igraph_error_t igraph_i_attribute_add_vertices(
    igraph_t *graph, igraph_int_t nv, const igraph_attribute_record_list_t *attr
) {
    if (igraph_i_attribute_table) {
        return igraph_i_attribute_table->add_vertices(graph, nv, attr);
    } else {
        return IGRAPH_SUCCESS;
    }
}

igraph_error_t igraph_i_attribute_permute_vertices(const igraph_t *graph,
                                        igraph_t *newgraph,
                                        const igraph_vector_int_t *idx) {
    /* graph and newgraph may be the same, in which case we need to support
     * in-place operations. If they are _not_ the same, it is assumed that the
     * new graph has no vertex attributes yet */
    if (igraph_i_attribute_table) {
        return igraph_i_attribute_table->permute_vertices(graph, newgraph, idx);
    } else {
        return IGRAPH_SUCCESS;
    }
}

igraph_error_t igraph_i_attribute_combine_vertices(const igraph_t *graph,
                                        igraph_t *newgraph,
                                        const igraph_vector_int_list_t *merges,
                                        const igraph_attribute_combination_t *comb) {
    /* It is assumed that the two graphs are not the same and that the new
     * graph has no vertex attributes yet. We cannot assert the latter but we
     * can assert the former */
    IGRAPH_ASSERT(graph != newgraph);
    if (igraph_i_attribute_table) {
        return igraph_i_attribute_table->combine_vertices(graph, newgraph,
                merges,
                comb);
    } else {
        return IGRAPH_SUCCESS;
    }
}

igraph_error_t igraph_i_attribute_add_edges(
    igraph_t *graph, const igraph_vector_int_t *edges,
    const igraph_attribute_record_list_t *attr
) {
    if (igraph_i_attribute_table) {
        return igraph_i_attribute_table->add_edges(graph, edges, attr);
    } else {
        return IGRAPH_SUCCESS;
    }
}

igraph_error_t igraph_i_attribute_permute_edges(const igraph_t *graph,
                                     igraph_t *newgraph,
                                     const igraph_vector_int_t *idx) {
    /* graph and newgraph may be the same, in which case we need to support
     * in-place operations. If they are _not_ the same, it is assumed that the
     * new graph has no edge attributes yet */
    if (igraph_i_attribute_table) {
        return igraph_i_attribute_table->permute_edges(graph, newgraph, idx);
    } else {
        return IGRAPH_SUCCESS;
    }
}

igraph_error_t igraph_i_attribute_combine_edges(const igraph_t *graph,
                                     igraph_t *newgraph,
                                     const igraph_vector_int_list_t *merges,
                                     const igraph_attribute_combination_t *comb) {
    /* It is assumed that the two graphs are not the same and that the new
     * graph has no eedge attributes yet. We cannot assert the latter but we
     * can assert the former */
    IGRAPH_ASSERT(graph != newgraph);
    if (igraph_i_attribute_table) {
        return igraph_i_attribute_table->combine_edges(graph, newgraph,
                merges,
                comb);
    } else {
        return IGRAPH_SUCCESS;
    }
}

igraph_error_t igraph_i_attribute_get_info(const igraph_t *graph,
                                igraph_strvector_t *gnames,
                                igraph_vector_int_t *gtypes,
                                igraph_strvector_t *vnames,
                                igraph_vector_int_t *vtypes,
                                igraph_strvector_t *enames,
                                igraph_vector_int_t *etypes) {
    if (igraph_i_attribute_table) {
        return igraph_i_attribute_table->get_info(graph, gnames, gtypes,
                vnames, vtypes,
                enames, etypes);
    } else {
        return IGRAPH_SUCCESS;
    }
}

igraph_bool_t igraph_i_attribute_has_attr(const igraph_t *graph,
        igraph_attribute_elemtype_t type,
        const char *name) {
    if (igraph_i_attribute_table) {
        return igraph_i_attribute_table->has_attr(graph, type, name);
    } else {
        return IGRAPH_SUCCESS;
    }
}

igraph_error_t igraph_i_attribute_get_type(const igraph_t *graph,
                               igraph_attribute_type_t *type,
                               igraph_attribute_elemtype_t elemtype,
                               const char *name) {
    if (igraph_i_attribute_table) {
        return igraph_i_attribute_table->get_type(graph, type, elemtype, name);
    } else {
        return IGRAPH_SUCCESS;
    }

}

igraph_error_t igraph_i_attribute_get_numeric_graph_attr(const igraph_t *graph,
        const char *name,
        igraph_vector_t *value) {
    igraph_vector_clear(value);
    if (igraph_i_attribute_table) {
        return igraph_i_attribute_table->get_numeric_graph_attr(graph, name, value);
    } else {
        return IGRAPH_SUCCESS;
    }
}

igraph_error_t igraph_i_attribute_get_numeric_vertex_attr(const igraph_t *graph,
        const char *name,
        igraph_vs_t vs,
        igraph_vector_t *value) {
    igraph_vector_clear(value);
    if (igraph_i_attribute_table) {
        return igraph_i_attribute_table->get_numeric_vertex_attr(graph, name, vs, value);
    } else {
        return IGRAPH_SUCCESS;
    }
}

igraph_error_t igraph_i_attribute_get_numeric_edge_attr(const igraph_t *graph,
        const char *name,
        igraph_es_t es,
        igraph_vector_t *value) {
    igraph_vector_clear(value);
    if (igraph_i_attribute_table) {
        return igraph_i_attribute_table->get_numeric_edge_attr(graph, name, es, value);
    } else {
        return IGRAPH_SUCCESS;
    }
}

igraph_error_t igraph_i_attribute_get_string_graph_attr(const igraph_t *graph,
        const char *name,
        igraph_strvector_t *value) {
    igraph_strvector_clear(value);
    if (igraph_i_attribute_table) {
        return igraph_i_attribute_table->get_string_graph_attr(graph, name, value);
    } else {
        return IGRAPH_SUCCESS;
    }
}

igraph_error_t igraph_i_attribute_get_string_vertex_attr(const igraph_t *graph,
        const char *name,
        igraph_vs_t vs,
        igraph_strvector_t *value) {
    igraph_strvector_clear(value);
    if (igraph_i_attribute_table) {
        return igraph_i_attribute_table->get_string_vertex_attr(graph, name, vs, value);
    } else {
        return IGRAPH_SUCCESS;
    }
}

igraph_error_t igraph_i_attribute_get_string_edge_attr(const igraph_t *graph,
        const char *name,
        igraph_es_t es,
        igraph_strvector_t *value) {
    igraph_strvector_clear(value);
    if (igraph_i_attribute_table) {
        return igraph_i_attribute_table->get_string_edge_attr(graph, name, es, value);
    } else {
        return IGRAPH_SUCCESS;
    }
}

igraph_error_t igraph_i_attribute_get_bool_graph_attr(const igraph_t *graph,
        const char *name,
        igraph_vector_bool_t *value) {
    igraph_vector_bool_clear(value);
    if (igraph_i_attribute_table) {
        return igraph_i_attribute_table->get_bool_graph_attr(graph, name, value);
    } else {
        return IGRAPH_SUCCESS;
    }
}

igraph_error_t igraph_i_attribute_get_bool_vertex_attr(const igraph_t *graph,
        const char *name,
        igraph_vs_t vs,
        igraph_vector_bool_t *value) {
    igraph_vector_bool_clear(value);
    if (igraph_i_attribute_table) {
        return igraph_i_attribute_table->get_bool_vertex_attr(graph, name, vs, value);
    } else {
        return IGRAPH_SUCCESS;
    }
}

igraph_error_t igraph_i_attribute_get_bool_edge_attr(const igraph_t *graph,
        const char *name,
        igraph_es_t es,
        igraph_vector_bool_t *value) {
    igraph_vector_bool_clear(value);
    if (igraph_i_attribute_table) {
        return igraph_i_attribute_table->get_bool_edge_attr(graph, name, es, value);
    } else {
        return IGRAPH_SUCCESS;
    }
}

/**
 * \function igraph_set_attribute_table
 * \brief Attach an attribute table.
 *
 * This function attaches attribute handling code to the igraph library.
 * Note that the attribute handler table is \em not thread-local even if
 * igraph is compiled in thread-local mode. In the vast majority of cases,
 * this is not a significant restriction.
 *
 * </para><para>
 * Attribute handlers are normally attached on program startup, and are
 * left active for the program's lifetime. This is because a graph object
 * created with a given attribute handler must not be manipulated while
 * a different attribute handler is active.
 *
 * \param table Pointer to an \ref igraph_attribute_table_t object
 *    containing the functions for attribute manipulation. Supply \c
 *    NULL here if you don't want attributes.
 * \return Pointer to the old attribute handling table.
 *
 * Time complexity: O(1).
 */

igraph_attribute_table_t *
igraph_set_attribute_table(const igraph_attribute_table_t * table) {
    igraph_attribute_table_t *old = igraph_i_attribute_table;
    igraph_i_attribute_table = (igraph_attribute_table_t*) table;
    return old;
}

igraph_bool_t igraph_has_attribute_table(void) {
    return igraph_i_attribute_table != NULL;
}


/**
 * \function igraph_attribute_combination_init
 * \brief Initialize attribute combination list.
 *
 * \param comb The uninitialized attribute combination list.
 * \return Error code.
 *
 * Time complexity: O(1)
 */
igraph_error_t igraph_attribute_combination_init(igraph_attribute_combination_t *comb) {
    IGRAPH_CHECK(igraph_vector_ptr_init(&comb->list, 0));
    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_attribute_combination_destroy
 * \brief Destroy attribute combination list.
 *
 * \param comb The attribute combination list.
 *
 * Time complexity: O(n), where n is the number of records in the
                    attribute combination list.
 */
void igraph_attribute_combination_destroy(igraph_attribute_combination_t *comb) {
    igraph_int_t i, n = igraph_vector_ptr_size(&comb->list);
    for (i = 0; i < n; i++) {
        igraph_attribute_combination_record_t *rec = VECTOR(comb->list)[i];
        if (rec->name) {
            IGRAPH_FREE(rec->name);
        }
        IGRAPH_FREE(rec);
    }
    igraph_vector_ptr_destroy(&comb->list);
}

/**
 * \function igraph_attribute_combination_add
 * \brief Add combination record to attribute combination list.
 *
 * \param comb The attribute combination list.
 * \param name The name of the attribute. If the name already exists
 *             the attribute combination record will be replaced.
 *             Use NULL to add a default combination record for all
 *             atributes not in the list.
 * \param type The type of the attribute combination. See \ref
 *             igraph_attribute_combination_type_t for the options.
 * \param func Function to be used if \p type is
 *             \c IGRAPH_ATTRIBUTE_COMBINE_FUNCTION. This function is called
 *             by the concrete attribute handler attached to igraph, and its
 *             calling signature depends completely on the attribute handler.
 *             For instance, if you are using attributes from C and you have
 *             attached the C attribute handler, you need to follow the
 *             documentation of the <link linkend="c-attribute-combination-functions">C attribute handler</link>
 *             for more details.
 * \return Error code.
 *
 * Time complexity: O(n), where n is the number of current attribute
 *                  combinations.
 */
igraph_error_t igraph_attribute_combination_add(igraph_attribute_combination_t *comb,
                                     const char *name,
                                     igraph_attribute_combination_type_t type,
                                     igraph_function_pointer_t func) {
    igraph_int_t i, n = igraph_vector_ptr_size(&comb->list);

    /* Search, in case it is already there */
    for (i = 0; i < n; i++) {
        igraph_attribute_combination_record_t *r = VECTOR(comb->list)[i];
        const char *n = r->name;
        if ( (!name && !n) ||
             (name && n && !strcmp(n, name)) ) {
            r->type = type;
            r->func = func;
            break;
        }
    }

    if (i == n) {
        /* This is a new attribute name */
        igraph_attribute_combination_record_t *rec =
            IGRAPH_CALLOC(1, igraph_attribute_combination_record_t);
        if (! rec) {
            IGRAPH_ERROR("Cannot create attribute combination data.",
                         IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */
        }
        IGRAPH_FINALLY(igraph_free, rec);
        if (! name) {
            rec->name = NULL;
        } else {
            rec->name = strdup(name);
            if (! rec->name) {
                IGRAPH_ERROR("Cannot create attribute combination data.",
                             IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */
            }
        }
        IGRAPH_FINALLY(igraph_free, (char *) rec->name); /* free() is safe on NULL */
        rec->type = type;
        rec->func = func;

        IGRAPH_CHECK(igraph_vector_ptr_push_back(&comb->list, rec));
        IGRAPH_FINALLY_CLEAN(2); /* ownership of 'rec' transferred to 'comb->list' */

    }

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_attribute_combination_remove
 * \brief Remove a record from an attribute combination list.
 *
 * \param comb The attribute combination list.
 * \param name The attribute name of the attribute combination record
 *             to remove. It will be ignored if the named attribute
 *             does not exist. It can be NULL to remove the default
 *             combination record.
 * \return Error code. This currently always returns IGRAPH_SUCCESS.
 *
 * Time complexity: O(n), where n is the number of records in the attribute
                    combination list.
 */
igraph_error_t igraph_attribute_combination_remove(igraph_attribute_combination_t *comb,
                                        const char *name) {
    igraph_int_t i, n = igraph_vector_ptr_size(&comb->list);

    /* Search, in case it is already there */
    for (i = 0; i < n; i++) {
        igraph_attribute_combination_record_t *r = VECTOR(comb->list)[i];
        const char *n = r->name;
        if ( (!name && !n) ||
             (name && n && !strcmp(n, name)) ) {
            break;
        }
    }

    if (i != n) {
        igraph_attribute_combination_record_t *r = VECTOR(comb->list)[i];
        if (r->name) {
            IGRAPH_FREE(r->name);
        }
        IGRAPH_FREE(r);
        igraph_vector_ptr_remove(&comb->list, i);
    } else {
        /* It is not there, we don't do anything */
    }

    return IGRAPH_SUCCESS;
}

igraph_error_t igraph_attribute_combination_query(const igraph_attribute_combination_t *comb,
                                       const char *name,
                                       igraph_attribute_combination_type_t *type,
                                       igraph_function_pointer_t *func) {
    igraph_int_t i, def = -1, len = igraph_vector_ptr_size(&comb->list);

    for (i = 0; i < len; i++) {
        igraph_attribute_combination_record_t *rec = VECTOR(comb->list)[i];
        const char *n = rec->name;
        if ( (!name && !n) ||
             (name && n && !strcmp(n, name)) ) {
            *type = rec->type;
            *func = rec->func;
            return IGRAPH_SUCCESS;
        }
        if (!n) {
            def = i;
        }
    }

    if (def == -1) {
        /* Did not find anything */
        *type = IGRAPH_ATTRIBUTE_COMBINE_DEFAULT;
        *func = 0;
    } else {
        igraph_attribute_combination_record_t *rec = VECTOR(comb->list)[def];
        *type = rec->type;
        *func = rec->func;
    }

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_attribute_combination
 * \brief Initialize attribute combination list and add records.
 *
 * \param comb The uninitialized attribute combination list.
 * \param ...  A list of 'name, type[, func]', where:
 * \param name The name of the attribute. If the name already exists
 *             the attribute combination record will be replaced.
 *             Use NULL to add a default combination record for all
 *             atributes not in the list.
 * \param type The type of the attribute combination. See \ref
 *             igraph_attribute_combination_type_t for the options.
 * \param func Function to be used if \p type is
 *             \c IGRAPH_ATTRIBUTE_COMBINE_FUNCTION.
 * The list is closed by setting the name to \c IGRAPH_NO_MORE_ATTRIBUTES.
 * \return Error code.
 *
 * Time complexity: O(n^2), where n is the number attribute
 *                  combinations records to add.
 *
 * \example examples/simple/igraph_attribute_combination.c
 */
igraph_error_t igraph_attribute_combination(
        igraph_attribute_combination_t *comb, ...) {

    va_list ap;

    IGRAPH_CHECK(igraph_attribute_combination_init(comb));

    va_start(ap, comb);
    while (true) {
        igraph_function_pointer_t func = NULL;
        igraph_attribute_combination_type_t type;
        const char *name;

        name = va_arg(ap, const char *);

        if (name == IGRAPH_NO_MORE_ATTRIBUTES) {
            break;
        }

        type = (igraph_attribute_combination_type_t) va_arg(ap, int);
        if (type == IGRAPH_ATTRIBUTE_COMBINE_FUNCTION) {
            func = va_arg(ap, igraph_function_pointer_t);
        }

        if (strlen(name) == 0) {
            name = 0;
        }

        igraph_error_t ret = igraph_attribute_combination_add(comb, name, type, func);
        if (ret != IGRAPH_SUCCESS) {
            va_end(ap);
            return ret;
        }
    }

    va_end(ap);

    return IGRAPH_SUCCESS;
}
