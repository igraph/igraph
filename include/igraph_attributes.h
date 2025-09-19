/*
   igraph library.
   Copyright (C) 2005-2025  The igraph development team <igraph@igraph.org>

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

#ifndef IGRAPH_ATTRIBUTES_H
#define IGRAPH_ATTRIBUTES_H

#include "igraph_decls.h"
#include "igraph_error.h"
#include "igraph_datatype.h"
#include "igraph_types.h"
#include "igraph_vector.h"
#include "igraph_strvector.h"
#include "igraph_vector_list.h"
#include "igraph_vector_ptr.h"
#include "igraph_iterators.h"

IGRAPH_BEGIN_C_DECLS

/* -------------------------------------------------- */
/* Attributes                                         */
/* -------------------------------------------------- */

/**
 * \typedef igraph_attribute_type_t
 * \brief The possible types of the attributes.
 *
 * Values of this enum are used by the attribute interface to communicate the
 * type of an attribute to igraph's C core. When igraph is integrated in a
 * high-level language, the attribute type reported by the interface may not
 * necessarily have to match the exact data type in the high-level language as
 * long as the attribute interface can provide a conversion from the native
 * high-level attribute value to one of the data types listed here. When the
 * high-level data type is complex and has no suitable conversion to one of the
 * atomic igraph attribute types (numeric, string or Boolean), the attribute
 * interface should report the attribute as having an "object" type, which is
 * ignored by the C core. See also \ref igraph_attribute_table_t.
 *
 * \enumval IGRAPH_ATTRIBUTE_UNSPECIFIED Currently used internally
 *   as a "null value" or "placeholder value" in some algorithms.
 *   Attribute records with this type must not be passed to igraph
 *   functions.
 * \enumval IGRAPH_ATTRIBUTE_NUMERIC Numeric attribute.
 * \enumval IGRAPH_ATTRIBUTE_BOOLEAN Logical values, true or false.
 * \enumval IGRAPH_ATTRIBUTE_STRING String attribute.
 * \enumval IGRAPH_ATTRIBUTE_OBJECT Custom attribute type, to be
 *   used for special data types by client applications. The R and
 *   Python interfaces use this for attributes that hold R or Python
 *   objects. Usually ignored by igraph functions.
 */
typedef enum {
    IGRAPH_ATTRIBUTE_UNSPECIFIED = 0,
    IGRAPH_ATTRIBUTE_NUMERIC = 1,
    IGRAPH_ATTRIBUTE_BOOLEAN = 2,
    IGRAPH_ATTRIBUTE_STRING = 3,
    IGRAPH_ATTRIBUTE_OBJECT = 127
} igraph_attribute_type_t;

/**
 * \typedef igraph_attribute_elemtype_t
 * \brief Types of objects to which attributes can be attached.
 *
 * \enumval IGRAPH_ATTRIBUTE_GRAPH Denotes that an attribute belongs to the
 *   entire graph.
 * \enumval IGRAPH_ATTRIBUTE_VERTEX Denotes that an attribute belongs to the
 *   vertices of a graph.
 * \enumval IGRAPH_ATTRIBUTE_EDGE Denotes that an attribute belongs to the
 *   edges of a graph.
 */
typedef enum {
    IGRAPH_ATTRIBUTE_GRAPH = 0,
    IGRAPH_ATTRIBUTE_VERTEX,
    IGRAPH_ATTRIBUTE_EDGE
} igraph_attribute_elemtype_t;

/* -------------------------------------------------- */
/* Attribute records                                  */
/* -------------------------------------------------- */

/**
 * \typedef igraph_attribute_record_t
 * \brief An attribute record holding the name, type and values of an attribute.
 *
 * This composite data type is used in the attribute interface to specify a
 * name-type-value triplet where the name is the name of a graph, vertex or
 * edge attribute, the type is the corresponding igraph type of the attribute
 * and the value is a \em vector of attribute values. Note that for graph
 * attributes we use a vector of length 1. The type of the vector depends on
 * the attribute type: it is \ref igraph_vector_t for numeric attributes,
 * \c igraph_strvector_t for string attributes and \c igraph_vector_bool_t
 * for Boolean attributes.
 *
 * </para><para>
 * The record also stores default values for the attribute. The default values
 * are used when the value vector of the record is resized with
 * \ref igraph_attribute_record_resize(). It is important that the record
 * stores \em one default value only, corresponding to the type of the
 * attribute record. The default value is \em cleared when the type of the
 * record is changed.
 */
typedef struct igraph_attribute_record_t {
    char *name;
    igraph_attribute_type_t type;
    union {
        void *as_raw;
        igraph_vector_t *as_vector;
        igraph_strvector_t *as_strvector;
        igraph_vector_bool_t *as_vector_bool;
    } value;
    union {
        igraph_real_t numeric;
        igraph_bool_t boolean;
        char *string;
    } default_value;
} igraph_attribute_record_t;

IGRAPH_EXPORT igraph_error_t igraph_attribute_record_init(
    igraph_attribute_record_t *attr, const char* name, igraph_attribute_type_t type
);
IGRAPH_EXPORT igraph_error_t igraph_attribute_record_init_copy(
    igraph_attribute_record_t *to, const igraph_attribute_record_t *from
);
IGRAPH_EXPORT igraph_error_t igraph_attribute_record_check_type(
    const igraph_attribute_record_t *attr, igraph_attribute_type_t type
);
IGRAPH_EXPORT igraph_int_t igraph_attribute_record_size(
    const igraph_attribute_record_t *attr
);
IGRAPH_EXPORT igraph_error_t igraph_attribute_record_resize(
    igraph_attribute_record_t *attr, igraph_int_t new_size
);
IGRAPH_EXPORT igraph_error_t igraph_attribute_record_set_name(
    igraph_attribute_record_t *attr, const char* name
);
IGRAPH_EXPORT igraph_error_t igraph_attribute_record_set_default_numeric(
    igraph_attribute_record_t *attr, igraph_real_t value
);
IGRAPH_EXPORT igraph_error_t igraph_attribute_record_set_default_boolean(
    igraph_attribute_record_t *attr, igraph_bool_t value
);
IGRAPH_EXPORT igraph_error_t igraph_attribute_record_set_default_string(
    igraph_attribute_record_t *attr, const char* value
);
IGRAPH_EXPORT igraph_error_t igraph_attribute_record_set_type(
    igraph_attribute_record_t *attr, igraph_attribute_type_t type
);
IGRAPH_EXPORT void igraph_attribute_record_destroy(igraph_attribute_record_t *attr);

/* -------------------------------------------------- */
/* Attribute combinations                             */
/* -------------------------------------------------- */

/**
 * \typedef igraph_attribute_combination_type_t
 * The possible types of attribute combinations.
 *
 * \enumval IGRAPH_ATTRIBUTE_COMBINE_IGNORE   Ignore old attributes, use an empty value.
 * \enumval IGRAPH_ATTRIBUTE_COMBINE_DEFAULT  Use the default way to combine attributes (decided by the attribute handler implementation).
 * \enumval IGRAPH_ATTRIBUTE_COMBINE_FUNCTION Supply your own function to combine
 *                                            attributes.
 * \enumval IGRAPH_ATTRIBUTE_COMBINE_SUM      Take the sum of the attributes.
 * \enumval IGRAPH_ATTRIBUTE_COMBINE_PROD     Take the product of the attributes.
 * \enumval IGRAPH_ATTRIBUTE_COMBINE_MIN      Take the minimum attribute.
 * \enumval IGRAPH_ATTRIBUTE_COMBINE_MAX      Take the maximum attribute.
 * \enumval IGRAPH_ATTRIBUTE_COMBINE_RANDOM   Take a random attribute.
 * \enumval IGRAPH_ATTRIBUTE_COMBINE_FIRST    Take the first attribute.
 * \enumval IGRAPH_ATTRIBUTE_COMBINE_LAST     Take the last attribute.
 * \enumval IGRAPH_ATTRIBUTE_COMBINE_MEAN     Take the mean of the attributes.
 * \enumval IGRAPH_ATTRIBUTE_COMBINE_MEDIAN   Take the median of the attributes.
 * \enumval IGRAPH_ATTRIBUTE_COMBINE_CONCAT   Concatenate the attributes.
 */
typedef enum {
    IGRAPH_ATTRIBUTE_COMBINE_IGNORE = 0,
    IGRAPH_ATTRIBUTE_COMBINE_DEFAULT = 1,
    IGRAPH_ATTRIBUTE_COMBINE_FUNCTION = 2,
    IGRAPH_ATTRIBUTE_COMBINE_SUM = 3,
    IGRAPH_ATTRIBUTE_COMBINE_PROD = 4,
    IGRAPH_ATTRIBUTE_COMBINE_MIN = 5,
    IGRAPH_ATTRIBUTE_COMBINE_MAX = 6,
    IGRAPH_ATTRIBUTE_COMBINE_RANDOM = 7,
    IGRAPH_ATTRIBUTE_COMBINE_FIRST = 8,
    IGRAPH_ATTRIBUTE_COMBINE_LAST = 9,
    IGRAPH_ATTRIBUTE_COMBINE_MEAN = 10,
    IGRAPH_ATTRIBUTE_COMBINE_MEDIAN = 11,
    IGRAPH_ATTRIBUTE_COMBINE_CONCAT = 12
} igraph_attribute_combination_type_t;

typedef void (*igraph_function_pointer_t)(void);

typedef struct igraph_attribute_combination_record_t {
    const char *name;     /* can be NULL, meaning: the rest */
    igraph_attribute_combination_type_t type;
    igraph_function_pointer_t func;
} igraph_attribute_combination_record_t;

typedef struct igraph_attribute_combination_t {
    igraph_vector_ptr_t list;
} igraph_attribute_combination_t;

#define IGRAPH_NO_MORE_ATTRIBUTES ((const char*)0)

IGRAPH_EXPORT igraph_error_t igraph_attribute_combination_init(igraph_attribute_combination_t *comb);
IGRAPH_EXPORT igraph_error_t igraph_attribute_combination(igraph_attribute_combination_t *comb, ...);
IGRAPH_EXPORT void igraph_attribute_combination_destroy(igraph_attribute_combination_t *comb);
IGRAPH_EXPORT igraph_error_t igraph_attribute_combination_add(igraph_attribute_combination_t *comb,
                                                   const char *name,
                                                   igraph_attribute_combination_type_t type,
                                                   igraph_function_pointer_t func);
IGRAPH_EXPORT igraph_error_t igraph_attribute_combination_remove(igraph_attribute_combination_t *comb,
                                                      const char *name);
IGRAPH_EXPORT igraph_error_t igraph_attribute_combination_query(const igraph_attribute_combination_t *comb,
                                                     const char *name,
                                                     igraph_attribute_combination_type_t *type,
                                                     igraph_function_pointer_t *func);

/* -------------------------------------------------- */
/* List of attribute records                          */
/* -------------------------------------------------- */

#define ATTRIBUTE_RECORD_LIST
#define BASE_ATTRIBUTE_RECORD
#include "igraph_pmt.h"
#include "igraph_typed_list_pmt.h"
#include "igraph_pmt_off.h"
#undef BASE_ATTRIBUTE_RECORD
#undef ATTRIBUTE_RECORD_LIST

/* -------------------------------------------------- */
/* Attribute handler interface                        */
/* -------------------------------------------------- */

/**
 * \struct igraph_attribute_table_t
 * \brief Table of functions to perform operations on attributes.
 *
 * This type collects the functions defining an attribute handler.
 * It has the following members:
 *
 * \member init This function is called whenever a new graph object is
 *    created, right after it is created but before any vertices or
 *    edges are added. It is supposed to set the \c attr member of the \c
 *    igraph_t object, which is guaranteed to be set to a null pointer
 *    before this function is called. It is expected to set the \c attr member
 *    to a non-null value \em or return an error code. Leaving the \c attr
 *    member at a null value while returning success is invalid and will trigger
 *    an error in the C core of igraph itself.
 * \member destroy This function is called whenever the graph object
 *    is destroyed, right before freeing the allocated memory. It is supposed
 *    to do any cleanup operations that are need to dispose of the \c attr
 *    member of the \c igraph_t object properly. The caller will set the
 *    \c attr member to a null pointer after this function returns.
 * \member copy This function is called when the C core wants to populate the
 *    attributes of a graph from another graph. The structure of the target
 *    graph is already initialized by the time this function is called, and the
 *    \c attr member of the graph is set to a null pointer. The function is
 *    supposed to populate the \c attr member of the target \c igraph_t object
 *    to a non-null value \em or return an error code. Leaving the \c attr
 *    member at a null value while returning success is invalid and will trigger
 *    an error in the C core of igraph itself.
 * \member add_vertices Called when vertices are added to a graph, after the
 *    base data structure was modified. The number of vertices that were added is
 *    supplied as an argument. The function is supposed to set up default values
 *    for each vertex attribute that is currently registered on the graph, for
 *    all the newly added vertices. Expected to return an error code.
 * \member permute_vertices Called when a new graph is created based on an
 *    existing one such that there is a mapping from the vertices of the new
 *    graph back to the vertices of the old graph (e.g. if vertices are removed
 *    from a graph). The supplied index vector defines which old vertex
 *    a new vertex corresponds to. Its length is the same as the number of
 *    vertices in the new graph, and for each new vertex it provides the ID
 *    of the corresponding vertex in the old graph. The function is supposed to
 *    set up the values of the vertex attributes of the new graph based on the
 *    attributes of the old graph and the provided index vector. Note that the
 *    old and the new graph \em may be the same, in which case it is the
 *    responsibility of the function to ensure that the operation can safely be
 *    performed in-place. If the two graph instances are \em not the same,
 *    implementors may safely assume that the new graph has no vertex attributes
 *    yet (but it may already have graph or edge attributes by the time this
 *    function is called).
 * \member combine_vertices This function is called when the creation
 *    of a new graph involves a merge (contraction, etc.) of vertices
 *    from another graph. The function is called after the new graph was created.
 *    An argument specifies how several vertices from the old graph map to a
 *    single vertex in the new graph. It is guaranteed that the old and the
 *    new graph instances are different when this callback is called.
 *    Implementors may safely assume that the new graph has no vertex attributes
 *    yet (but it may already have graph or edge attributes by the time this
 *    function is called).
 * \member add_edges Called when new edges are added to a graph, after the
 *    base data structure was modified. A vector containing the endpoints of the
 *    new edges are supplied as an argument. The function is supposed to set up
 *    default values for each edge attribute that is currently registered on the
 *    graph, for all the newly added edges. Expected to return an error code.
 * \member permute_edges Called when a new graph is created based on an
 *    existing one such that some of the edges in the new graph should copy the
 *    attributes of some edges from the old graph (this also includes the
 *    deletion of edges). The supplied index vector defines which old edge a new
 *    edge corresponds to. Its length is the same as the number of edges in the
 *    new graph, and for each edge it provides the ID of the corresponding edge
 *    in the old graph. The function is supposed to set up the values of the
 *    edge attributes of the new graph based on the attributes of the old graph
 *    and the provided index vector. Note that the old and the new graph \em may
 *    be the same, in which case it is the responsibility of the function to
 *    ensure that the operation can safely be performed in-place. If the two
 *    graph instances are \em not the same, implementors may safely assume that
 *    the new graph has no edge attributes yet (but it may already have graph or
 *    vertex attributes by the time this function is called).
 * \member combine_edges This function is called when the creation
 *    of a new graph involves a merge (contraction, etc.) of edges
 *    from another graph. The function is after the new graph was created.
 *    An argument specifies how several edges from the old graph map to a
 *    single edge in the new graph. It is guaranteed that the old and the
 *    new graph instances are different when this callback is called.
 *    Implementors may safely assume that the new graph has no edge attributes
 *    yet (but it may already have graph or vertex attributes by the time this
 *    function is called).
 * \member get_info Query the attributes of a graph, the names and
 *    types should be returned.
 * \member has_attr Check whether a graph has the named
 *    graph/vertex/edge attribute.
 * \member get_type Query the type of a graph/vertex/edge attribute.
 * \member get_numeric_graph_attr Query a numeric graph attribute. The
 *    value should be appended to the provided \p value vector. No assumptions
 *    should be made about the initial contents of the \p value vector and it is
 *    not guaranteed to be empty.
 * \member get_string_graph_attr Query a string graph attribute. The
 *    value should be appended to the provided \p value vector. No assumptions
 *    should be made about the initial contents of the \p value vector and it is
 *    not guaranteed to be empty.
 * \member get_bool_graph_attr Query a boolean graph attribute. The
 *    value should be appended to the provided \p value vector. No assumptions
 *    should be made about the initial contents of the \p value vector and it is
 *    not guaranteed to be empty.
 * \member get_numeric_vertex_attr Query a numeric vertex attribute,
 *    for the vertices included in \p vs. The attribute values should be
 *    appended to the provided \p value vector. No assumptions should be made
 *    about the initial contents of the \p value vector and it is not guaranteed
 *    to be empty.
 * \member get_string_vertex_attr Query a string vertex attribute,
 *    for the vertices included in \p vs. The attribute values should be
 *    appended to the provided \p value vector. No assumptions should be made
 *    about the initial contents of the \p value vector and it is not guaranteed
 *    to be empty.
 * \member get_bool_vertex_attr Query a boolean vertex attribute,
 *    for the vertices included in \p vs. The attribute values should be
 *    appended to the provided \p value vector. No assumptions should be made
 *    about the initial contents of the \p value vector and it is not guaranteed
 *    to be empty.
 * \member get_numeric_edge_attr Query a numeric edge attribute, for
 *    the edges included in \p es. The attribute values should be appended
 *    to the provided \p value vector. No assumptions should be made
 *    about the initial contents of the \p value vector and it is not guaranteed
 *    to be empty.
 * \member get_string_edge_attr Query a string edge attribute, for the
 *    the edges included in \p es. The attribute values should be appended
 *    to the provided \p value vector. No assumptions should be made
 *    about the initial contents of the \p value vector and it is not guaranteed
 *    to be empty.
 * \member get_bool_edge_attr Query a boolean edge attribute, for the
 *    the edges included in \p es. The attribute values should be appended
 *    to the provided \p value vector. No assumptions should be made
 *    about the initial contents of the \p value vector and it is not guaranteed
 *    to be empty.
 */

typedef struct igraph_attribute_table_t {
    igraph_error_t (*init)(igraph_t *graph, const igraph_attribute_record_list_t *attr);
    void           (*destroy)(igraph_t *graph);
    igraph_error_t (*copy)(igraph_t *to, const igraph_t *from, igraph_bool_t ga,
                           igraph_bool_t va, igraph_bool_t ea);
    igraph_error_t (*add_vertices)(
        igraph_t *graph, igraph_int_t nv,
        const igraph_attribute_record_list_t *attr
    );
    igraph_error_t (*permute_vertices)(const igraph_t *graph,
                                       igraph_t *newgraph,
                                       const igraph_vector_int_t *idx);
    igraph_error_t (*combine_vertices)(const igraph_t *graph,
                                       igraph_t *newgraph,
                                       const igraph_vector_int_list_t *merges,
                                       const igraph_attribute_combination_t *comb);
    igraph_error_t (*add_edges)(
        igraph_t *graph, const igraph_vector_int_t *edges,
        const igraph_attribute_record_list_t *attr
    );
    igraph_error_t (*permute_edges)(const igraph_t *graph,
                                    igraph_t *newgraph, const igraph_vector_int_t *idx);
    igraph_error_t (*combine_edges)(const igraph_t *graph,
                                    igraph_t *newgraph,
                                    const igraph_vector_int_list_t *merges,
                                    const igraph_attribute_combination_t *comb);
    igraph_error_t (*get_info)(const igraph_t *graph,
                               igraph_strvector_t *gnames, igraph_vector_int_t *gtypes,
                               igraph_strvector_t *vnames, igraph_vector_int_t *vtypes,
                               igraph_strvector_t *enames, igraph_vector_int_t *etypes);
    igraph_bool_t (*has_attr)(const igraph_t *graph, igraph_attribute_elemtype_t type,
                              const char *name);
    igraph_error_t (*get_type)(const igraph_t *graph, igraph_attribute_type_t *type,
                              igraph_attribute_elemtype_t elemtype, const char *name);
    igraph_error_t (*get_numeric_graph_attr)(const igraph_t *graph, const char *name,
                                             igraph_vector_t *value);
    igraph_error_t (*get_string_graph_attr)(const igraph_t *graph, const char *name,
                                            igraph_strvector_t *value);
    igraph_error_t (*get_bool_graph_attr)(const igraph_t *igraph, const char *name,
                                          igraph_vector_bool_t *value);
    igraph_error_t (*get_numeric_vertex_attr)(const igraph_t *graph, const char *name,
                                              igraph_vs_t vs,
                                              igraph_vector_t *value);
    igraph_error_t (*get_string_vertex_attr)(const igraph_t *graph, const char *name,
                                             igraph_vs_t vs,
                                             igraph_strvector_t *value);
    igraph_error_t (*get_bool_vertex_attr)(const igraph_t *graph, const char *name,
                                           igraph_vs_t vs,
                                           igraph_vector_bool_t *value);
    igraph_error_t (*get_numeric_edge_attr)(const igraph_t *graph, const char *name,
                                            igraph_es_t es,
                                            igraph_vector_t *value);
    igraph_error_t (*get_string_edge_attr)(const igraph_t *graph, const char *name,
                                           igraph_es_t es,
                                           igraph_strvector_t *value);
    igraph_error_t (*get_bool_edge_attr)(const igraph_t *graph, const char *name,
                                         igraph_es_t es,
                                         igraph_vector_bool_t *value);
} igraph_attribute_table_t;

IGRAPH_EXPORT igraph_attribute_table_t * igraph_set_attribute_table(const igraph_attribute_table_t * table);

IGRAPH_EXPORT igraph_bool_t igraph_has_attribute_table(void);

/* Experimental attribute handler in C */

IGRAPH_EXPORT extern const igraph_attribute_table_t igraph_cattribute_table;

IGRAPH_EXPORT igraph_real_t igraph_cattribute_GAN(const igraph_t *graph, const char *name);
IGRAPH_EXPORT igraph_bool_t igraph_cattribute_GAB(const igraph_t *graph, const char *name);
IGRAPH_EXPORT const char* igraph_cattribute_GAS(const igraph_t *graph, const char *name);
IGRAPH_EXPORT igraph_real_t igraph_cattribute_VAN(const igraph_t *graph, const char *name,
                                                  igraph_int_t vid);
IGRAPH_EXPORT igraph_bool_t igraph_cattribute_VAB(const igraph_t *graph, const char *name,
                                                  igraph_int_t vid);
IGRAPH_EXPORT const char* igraph_cattribute_VAS(const igraph_t *graph, const char *name,
                                                igraph_int_t vid);
IGRAPH_EXPORT igraph_real_t igraph_cattribute_EAN(const igraph_t *graph, const char *name,
                                                  igraph_int_t eid);
IGRAPH_EXPORT igraph_bool_t igraph_cattribute_EAB(const igraph_t *graph, const char *name,
                                                  igraph_int_t eid);
IGRAPH_EXPORT const char* igraph_cattribute_EAS(const igraph_t *graph, const char *name,
                                                igraph_int_t eid);

IGRAPH_EXPORT igraph_error_t igraph_cattribute_VANV(const igraph_t *graph, const char *name,
                                         igraph_vs_t vids, igraph_vector_t *result);
IGRAPH_EXPORT igraph_error_t igraph_cattribute_EANV(const igraph_t *graph, const char *name,
                                         igraph_es_t eids, igraph_vector_t *result);
IGRAPH_EXPORT igraph_error_t igraph_cattribute_VASV(const igraph_t *graph, const char *name,
                                         igraph_vs_t vids, igraph_strvector_t *result);
IGRAPH_EXPORT igraph_error_t igraph_cattribute_EASV(const igraph_t *graph, const char *name,
                                         igraph_es_t eids, igraph_strvector_t *result);
IGRAPH_EXPORT igraph_error_t igraph_cattribute_VABV(const igraph_t *graph, const char *name,
                                         igraph_vs_t vids, igraph_vector_bool_t *result);
IGRAPH_EXPORT igraph_error_t igraph_cattribute_EABV(const igraph_t *graph, const char *name,
                                         igraph_es_t eids, igraph_vector_bool_t *result);

IGRAPH_EXPORT igraph_error_t igraph_cattribute_list(const igraph_t *graph,
                                         igraph_strvector_t *gnames, igraph_vector_int_t *gtypes,
                                         igraph_strvector_t *vnames, igraph_vector_int_t *vtypes,
                                         igraph_strvector_t *enames, igraph_vector_int_t *etypes);
IGRAPH_EXPORT igraph_bool_t igraph_cattribute_has_attr(const igraph_t *graph,
                                                       igraph_attribute_elemtype_t type,
                                                       const char *name);

IGRAPH_EXPORT igraph_error_t igraph_cattribute_GAN_set(igraph_t *graph, const char *name,
                                            igraph_real_t value);
IGRAPH_EXPORT igraph_error_t igraph_cattribute_GAB_set(igraph_t *graph, const char *name,
                                            igraph_bool_t value);
IGRAPH_EXPORT igraph_error_t igraph_cattribute_GAS_set(igraph_t *graph, const char *name,
                                            const char *value);
IGRAPH_EXPORT igraph_error_t igraph_cattribute_VAN_set(igraph_t *graph, const char *name,
                                            igraph_int_t vid, igraph_real_t value);
IGRAPH_EXPORT igraph_error_t igraph_cattribute_VAB_set(igraph_t *graph, const char *name,
                                            igraph_int_t vid, igraph_bool_t value);
IGRAPH_EXPORT igraph_error_t igraph_cattribute_VAS_set(igraph_t *graph, const char *name,
                                            igraph_int_t vid, const char *value);
IGRAPH_EXPORT igraph_error_t igraph_cattribute_EAN_set(igraph_t *graph, const char *name,
                                            igraph_int_t eid, igraph_real_t value);
IGRAPH_EXPORT igraph_error_t igraph_cattribute_EAB_set(igraph_t *graph, const char *name,
                                            igraph_int_t eid, igraph_bool_t value);
IGRAPH_EXPORT igraph_error_t igraph_cattribute_EAS_set(igraph_t *graph, const char *name,
                                            igraph_int_t eid, const char *value);

IGRAPH_EXPORT igraph_error_t igraph_cattribute_VAN_setv(igraph_t *graph, const char *name,
                                             const igraph_vector_t *v);
IGRAPH_EXPORT igraph_error_t igraph_cattribute_VAB_setv(igraph_t *graph, const char *name,
                                             const igraph_vector_bool_t *v);
IGRAPH_EXPORT igraph_error_t igraph_cattribute_VAS_setv(igraph_t *graph, const char *name,
                                             const igraph_strvector_t *sv);
IGRAPH_EXPORT igraph_error_t igraph_cattribute_EAN_setv(igraph_t *graph, const char *name,
                                             const igraph_vector_t *v);
IGRAPH_EXPORT igraph_error_t igraph_cattribute_EAB_setv(igraph_t *graph, const char *name,
                                             const igraph_vector_bool_t *v);
IGRAPH_EXPORT igraph_error_t igraph_cattribute_EAS_setv(igraph_t *graph, const char *name,
                                             const igraph_strvector_t *sv);

IGRAPH_EXPORT void igraph_cattribute_remove_g(igraph_t *graph, const char *name);
IGRAPH_EXPORT void igraph_cattribute_remove_v(igraph_t *graph, const char *name);
IGRAPH_EXPORT void igraph_cattribute_remove_e(igraph_t *graph, const char *name);
IGRAPH_EXPORT void igraph_cattribute_remove_all(igraph_t *graph, igraph_bool_t g,
                                                igraph_bool_t v, igraph_bool_t e);

/**
 * \define GAN
 * Query a numeric graph attribute.
 *
 * This is shorthand for \ref igraph_cattribute_GAN().
 * \param graph The graph.
 * \param n The name of the attribute.
 * \return The value of the attribute.
 */
#define GAN(graph,n) (igraph_cattribute_GAN((graph), (n)))
/**
 * \define GAB
 * Query a boolean graph attribute.
 *
 * This is shorthand for \ref igraph_cattribute_GAB().
 * \param graph The graph.
 * \param n The name of the attribute.
 * \return The value of the attribute.
 */
#define GAB(graph,n) (igraph_cattribute_GAB((graph), (n)))
/**
 * \define GAS
 * Query a string graph attribute.
 *
 * This is shorthand for \ref igraph_cattribute_GAS().
 * \param graph The graph.
 * \param n The name of the attribute.
 * \return The value of the attribute.
 */
#define GAS(graph,n) (igraph_cattribute_GAS((graph), (n)))
/**
 * \define VAN
 * Query a numeric vertex attribute.
 *
 * This is shorthand for \ref igraph_cattribute_VAN().
 * \param graph The graph.
 * \param n The name of the attribute.
 * \param v The id of the vertex.
 * \return The value of the attribute.
 */
#define VAN(graph,n,v) (igraph_cattribute_VAN((graph), (n), (v)))
/**
 * \define VAB
 * Query a boolean vertex attribute.
 *
 * This is shorthand for \ref igraph_cattribute_VAB().
 * \param graph The graph.
 * \param n The name of the attribute.
 * \param v The id of the vertex.
 * \return The value of the attribute.
 */
#define VAB(graph,n,v) (igraph_cattribute_VAB((graph), (n), (v)))
/**
 * \define VAS
 * Query a string vertex attribute.
 *
 * This is shorthand for \ref igraph_cattribute_VAS().
 * \param graph The graph.
 * \param n The name of the attribute.
 * \param v The id of the vertex.
 * \return The value of the attribute.
 */
#define VAS(graph,n,v) (igraph_cattribute_VAS((graph), (n), (v)))
/**
 * \define VANV
 * Query a numeric vertex attribute for all vertices.
 *
 * This is a shorthand for \ref igraph_cattribute_VANV().
 * \param graph The graph.
 * \param n The name of the attribute.
 * \param vec Pointer to an initialized vector, the result is
 *        stored here. It will be resized, if needed.
 * \return Error code.
 */
#define VANV(graph,n,vec) (igraph_cattribute_VANV((graph),(n), \
                           igraph_vss_all(), (vec)))
/**
 * \define VABV
 * Query a boolean vertex attribute for all vertices.
 *
 * This is a shorthand for \ref igraph_cattribute_VABV().
 * \param graph The graph.
 * \param n The name of the attribute.
 * \param vec Pointer to an initialized boolean vector, the result is
 *        stored here. It will be resized, if needed.
 * \return Error code.
 */
#define VABV(graph,n,vec) (igraph_cattribute_VABV((graph),(n), \
                           igraph_vss_all(), (vec)))
/**
 * \define VASV
 * Query a string vertex attribute for all vertices.
 *
 * This is a shorthand for \ref igraph_cattribute_VASV().
 * \param graph The graph.
 * \param n The name of the attribute.
 * \param vec Pointer to an initialized string vector, the result is
 *        stored here. It will be resized, if needed.
 * \return Error code.
 */
#define VASV(graph,n,vec) (igraph_cattribute_VASV((graph),(n), \
                           igraph_vss_all(), (vec)))
/**
 * \define EAN
 * Query a numeric edge attribute.
 *
 * This is shorthand for \ref igraph_cattribute_EAN().
 * \param graph The graph.
 * \param n The name of the attribute.
 * \param e The id of the edge.
 * \return The value of the attribute.
 */
#define EAN(graph,n,e) (igraph_cattribute_EAN((graph), (n), (e)))
/**
 * \define EAB
 * Query a boolean edge attribute.
 *
 * This is shorthand for \ref igraph_cattribute_EAB().
 * \param graph The graph.
 * \param n The name of the attribute.
 * \param e The id of the edge.
 * \return The value of the attribute.
 */
#define EAB(graph,n,e) (igraph_cattribute_EAB((graph), (n), (e)))
/**
 * \define EAS
 * Query a string edge attribute.
 *
 * This is shorthand for \ref igraph_cattribute_EAS().
 * \param graph The graph.
 * \param n The name of the attribute.
 * \param e The id of the edge.
 * \return The value of the attribute.
 */
#define EAS(graph,n,e) (igraph_cattribute_EAS((graph), (n), (e)))
/**
 * \define EANV
 * Query a numeric edge attribute for all edges.
 *
 * This is a shorthand for \ref igraph_cattribute_EANV().
 * \param graph The graph.
 * \param n The name of the attribute.
 * \param vec Pointer to an initialized vector, the result is
 *        stored here. It will be resized, if needed.
 * \return Error code.
 */
#define EANV(graph,n,vec) (igraph_cattribute_EANV((graph),(n), \
                           igraph_ess_all(IGRAPH_EDGEORDER_ID), (vec)))
/**
 * \define EABV
 * Query a boolean edge attribute for all edges.
 *
 * This is a shorthand for \ref igraph_cattribute_EABV().
 * \param graph The graph.
 * \param n The name of the attribute.
 * \param vec Pointer to an initialized vector, the result is
 *        stored here. It will be resized, if needed.
 * \return Error code.
 */
#define EABV(graph,n,vec) (igraph_cattribute_EABV((graph),(n), \
                           igraph_ess_all(IGRAPH_EDGEORDER_ID), (vec)))

/**
 * \define EASV
 * Query a string edge attribute for all edges.
 *
 * This is a shorthand for \ref igraph_cattribute_EASV().
 * \param graph The graph.
 * \param n The name of the attribute.
 * \param vec Pointer to an initialized string vector, the result is
 *        stored here. It will be resized, if needed.
 * \return Error code.
 */
#define EASV(graph,n,vec) (igraph_cattribute_EASV((graph),(n), \
                           igraph_ess_all(IGRAPH_EDGEORDER_ID), (vec)))
/**
 * \define SETGAN
 * Set a numeric graph attribute
 *
 * This is a shorthand for \ref igraph_cattribute_GAN_set().
 * \param graph The graph.
 * \param n The name of the attribute.
 * \param value The new value of the attribute.
 * \return Error code.
 */
#define SETGAN(graph,n,value) (igraph_cattribute_GAN_set((graph),(n),(value)))
/**
 * \define SETGAB
 * Set a boolean graph attribute
 *
 * This is a shorthand for \ref igraph_cattribute_GAB_set().
 * \param graph The graph.
 * \param n The name of the attribute.
 * \param value The new value of the attribute.
 * \return Error code.
 */
#define SETGAB(graph,n,value) (igraph_cattribute_GAB_set((graph),(n),(value)))
/**
 * \define SETGAS
 * Set a string graph attribute
 *
 * This is a shorthand for \ref igraph_cattribute_GAS_set().
 * \param graph The graph.
 * \param n The name of the attribute.
 * \param value The new value of the attribute.
 * \return Error code.
 */
#define SETGAS(graph,n,value) (igraph_cattribute_GAS_set((graph),(n),(value)))
/**
 * \define SETVAN
 * Set a numeric vertex attribute
 *
 * This is a shorthand for \ref igraph_cattribute_VAN_set().
 * \param graph The graph.
 * \param n The name of the attribute.
 * \param vid Ids of the vertices to set.
 * \param value The new value of the attribute.
 * \return Error code.
 */
#define SETVAN(graph,n,vid,value) (igraph_cattribute_VAN_set((graph),(n),(vid),(value)))
/**
 * \define SETVAB
 * Set a boolean vertex attribute
 *
 * This is a shorthand for \ref igraph_cattribute_VAB_set().
 * \param graph The graph.
 * \param n The name of the attribute.
 * \param vid Ids of the vertices to set.
 * \param value The new value of the attribute.
 * \return Error code.
 */
#define SETVAB(graph,n,vid,value) (igraph_cattribute_VAB_set((graph),(n),(vid),(value)))
/**
 * \define SETVAS
 * Set a string vertex attribute
 *
 * This is a shorthand for \ref igraph_cattribute_VAS_set().
 * \param graph The graph.
 * \param n The name of the attribute.
 * \param vid Ids of the vertices to set.
 * \param value The new value of the attribute.
 * \return Error code.
 */
#define SETVAS(graph,n,vid,value) (igraph_cattribute_VAS_set((graph),(n),(vid),(value)))
/**
 * \define SETEAN
 * Set a numeric edge attribute
 *
 * This is a shorthand for \ref igraph_cattribute_EAN_set().
 * \param graph The graph.
 * \param n The name of the attribute.
 * \param eid Ids of the edges to set.
 * \param value The new value of the attribute.
 * \return Error code.
 */
#define SETEAN(graph,n,eid,value) (igraph_cattribute_EAN_set((graph),(n),(eid),(value)))
/**
 * \define SETEAB
 * Set a boolean edge attribute
 *
 * This is a shorthand for \ref igraph_cattribute_EAB_set().
 * \param graph The graph.
 * \param n The name of the attribute.
 * \param eid Ids of the edges to set.
 * \param value The new value of the attribute.
 * \return Error code.
 */
#define SETEAB(graph,n,eid,value) (igraph_cattribute_EAB_set((graph),(n),(eid),(value)))
/**
 * \define SETEAS
 * Set a string edge attribute
 *
 * This is a shorthand for \ref igraph_cattribute_EAS_set().
 * \param graph The graph.
 * \param n The name of the attribute.
 * \param eid Ids of the edges to set.
 * \param value The new value of the attribute.
 * \return Error code.
 */
#define SETEAS(graph,n,eid,value) (igraph_cattribute_EAS_set((graph),(n),(eid),(value)))

/**
 * \define SETVANV
 *  Set a numeric vertex attribute for all vertices
 *
 * This is a shorthand for \ref igraph_cattribute_VAN_setv().
 * \param graph The graph.
 * \param n The name of the attribute.
 * \param v Vector containing the new values of the attributes.
 * \return Error code.
 */
#define SETVANV(graph,n,v) (igraph_cattribute_VAN_setv((graph),(n),(v)))
/**
 * \define SETVABV
 *  Set a boolean vertex attribute for all vertices
 *
 * This is a shorthand for \ref igraph_cattribute_VAB_setv().
 * \param graph The graph.
 * \param n The name of the attribute.
 * \param v Vector containing the new values of the attributes.
 * \return Error code.
 */
#define SETVABV(graph,n,v) (igraph_cattribute_VAB_setv((graph),(n),(v)))
/**
 * \define SETVASV
 *  Set a string vertex attribute for all vertices
 *
 * This is a shorthand for \ref igraph_cattribute_VAS_setv().
 * \param graph The graph.
 * \param n The name of the attribute.
 * \param v Vector containing the new values of the attributes.
 * \return Error code.
 */
#define SETVASV(graph,n,v) (igraph_cattribute_VAS_setv((graph),(n),(v)))
/**
 * \define SETEANV
 *  Set a numeric edge attribute for all edges
 *
 * This is a shorthand for \ref igraph_cattribute_EAN_setv().
 * \param graph The graph.
 * \param n The name of the attribute.
 * \param v Vector containing the new values of the attributes.
 */
#define SETEANV(graph,n,v) (igraph_cattribute_EAN_setv((graph),(n),(v)))
/**
 * \define SETEABV
 *  Set a boolean edge attribute for all edges
 *
 * This is a shorthand for \ref igraph_cattribute_EAB_setv().
 * \param graph The graph.
 * \param n The name of the attribute.
 * \param v Vector containing the new values of the attributes.
 */
#define SETEABV(graph,n,v) (igraph_cattribute_EAB_setv((graph),(n),(v)))
/**
 * \define SETEASV
 *  Set a string edge attribute for all edges
 *
 * This is a shorthand for \ref igraph_cattribute_EAS_setv().
 * \param graph The graph.
 * \param n The name of the attribute.
 * \param v Vector containing the new values of the attributes.
 */
#define SETEASV(graph,n,v) (igraph_cattribute_EAS_setv((graph),(n),(v)))

/**
 * \define DELGA
 * Remove a graph attribute.
 *
 * A shorthand for \ref igraph_cattribute_remove_g().
 * \param graph The graph.
 * \param n The name of the attribute to remove.
 */
#define DELGA(graph,n) (igraph_cattribute_remove_g((graph),(n)))
/**
 * \define DELVA
 * Remove a vertex attribute.
 *
 * A shorthand for \ref igraph_cattribute_remove_v().
 * \param graph The graph.
 * \param n The name of the attribute to remove.
 */
#define DELVA(graph,n) (igraph_cattribute_remove_v((graph),(n)))
/**
 * \define DELEA
 * Remove an edge attribute.
 *
 * A shorthand for \ref igraph_cattribute_remove_e().
 * \param graph The graph.
 * \param n The name of the attribute to remove.
 */
#define DELEA(graph,n) (igraph_cattribute_remove_e((graph),(n)))
/**
 * \define DELGAS
 * Remove all graph attributes.
 *
 * Calls \ref igraph_cattribute_remove_all().
 * \param graph The graph.
 */
#define DELGAS(graph) (igraph_cattribute_remove_all((graph),1,0,0))
/**
 * \define DELVAS
 * Remove all vertex attributes.
 *
 * Calls \ref igraph_cattribute_remove_all().
 * \param graph The graph.
 */
#define DELVAS(graph) (igraph_cattribute_remove_all((graph),0,1,0))
/**
 * \define DELEAS
 * Remove all edge attributes.
 *
 * Calls \ref igraph_cattribute_remove_all().
 * \param graph The graph.
 */
#define DELEAS(graph) (igraph_cattribute_remove_all((graph),0,0,1))
/**
 * \define DELALL
 * Remove all attributes.
 *
 * All graph, vertex and edges attributes will be removed.
 * Calls \ref igraph_cattribute_remove_all().
 * \param graph The graph.
 */
#define DELALL(graph) (igraph_cattribute_remove_all((graph),1,1,1))

IGRAPH_END_C_DECLS

#endif
