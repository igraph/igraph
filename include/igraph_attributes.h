/* -*- mode: C -*-  */
/*
   IGraph library.
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
   Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA
   02110-1301 USA

*/

#ifndef IGRAPH_ATTRIBUTES_H
#define IGRAPH_ATTRIBUTES_H

#include "igraph_decls.h"
#include "igraph_datatype.h"
#include "igraph_types.h"
#include "igraph_vector.h"
#include "igraph_strvector.h"
#include "igraph_vector_ptr.h"
#include "igraph_iterators.h"

__BEGIN_DECLS

/* -------------------------------------------------- */
/* Attributes                                         */
/* -------------------------------------------------- */

/**
 * \section about_attributes
 *
 * <para>Attributes are numbers or strings (or basically any kind
 * of data) associated with the vertices or edges of a graph, or
 * with the graph itself. Eg. you may label vertices with symbolic names
 * or attach numeric weights to the edges of a graph. </para>
 *
 * <para>igraph attributes are designed to be flexible and extensible.
 * In igraph attributes are implemented via an interface abstraction:
 * any type implementing the functions in the interface, can be used
 * for storing vertex, edge and graph attributes. This means that
 * different attribute implementations can be used together with
 * igraph. This is reasonable: if igraph is used from Python attributes can be
 * of any Python type, from GNU R all R types are allowed. There is an
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
 * <para>By default there is no attribute interface attached to \a igraph,
 * to attach one, call \ref igraph_i_set_attribute_table with your new
 * table. </para>
 *
 */

/**
 * \typedef igraph_attribute_type_t
 * The possible types of the attributes.
 *
 * Note that this is only the
 * type communicated by the attribute interface towards igraph
 * functions. Eg. in the GNU R attribute handler, it is safe to say
 * that all complex R object attributes are strings, as long as this
 * interface is able to serialize them into strings. See also \ref
 * igraph_attribute_table_t.
 * \enumval IGRAPH_ATTRIBUTE_DEFAULT Currently not used for anything.
 * \enumval IGRAPH_ATTRIBUTE_NUMERIC Numeric attribute.
 * \enumval IGRAPH_ATTRIBUTE_BOOLEAN Logical values, true or false.
 * \enumval IGRAPH_ATTRIBUTE_STRING Attribute that can be converted to
 *   a string.
 * \enumval IGRAPH_ATTRIBUTE_R_OBJECT An R object. This is usually
 *   ignored by the igraph functions.
 * \enumval IGRAPH_ATTRIBUTE_PY_OBJECT A Python object. Usually
 *   ignored by the igraph functions.
 *
 */
typedef enum { IGRAPH_ATTRIBUTE_DEFAULT = 0,
               IGRAPH_ATTRIBUTE_NUMERIC = 1,
               IGRAPH_ATTRIBUTE_BOOLEAN = 5,
               IGRAPH_ATTRIBUTE_STRING = 2,
               IGRAPH_ATTRIBUTE_R_OBJECT = 3,
               IGRAPH_ATTRIBUTE_PY_OBJECT = 4
             } igraph_attribute_type_t;

typedef struct igraph_attribute_record_t {
    const char *name;
    igraph_attribute_type_t type;
    const void *value;
} igraph_attribute_record_t;

typedef enum { IGRAPH_ATTRIBUTE_GRAPH = 0,
               IGRAPH_ATTRIBUTE_VERTEX,
               IGRAPH_ATTRIBUTE_EDGE
             } igraph_attribute_elemtype_t;

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

DECLDIR int igraph_attribute_combination_init(igraph_attribute_combination_t *comb);
DECLDIR int igraph_attribute_combination(igraph_attribute_combination_t *comb, ...);
DECLDIR void igraph_attribute_combination_destroy(igraph_attribute_combination_t *comb);
DECLDIR int igraph_attribute_combination_add(igraph_attribute_combination_t *comb,
        const char *name,
        igraph_attribute_combination_type_t type,
        igraph_function_pointer_t func);
DECLDIR int igraph_attribute_combination_remove(igraph_attribute_combination_t *comb,
        const char *name);
DECLDIR int igraph_attribute_combination_query(const igraph_attribute_combination_t *comb,
        const char *name,
        igraph_attribute_combination_type_t *type,
        igraph_function_pointer_t *func);

/**
 * \struct igraph_attribute_table_t
 * \brief Table of functions to perform operations on attributes
 *
 * This type collects the functions defining an attribute handler.
 * It has the following members:
 * \member init This function is called whenever a new graph object is
 *    created, right after it is created but before any vertices or
 *    edges are added. It is supposed to set the \c attr member of the \c
 *    igraph_t object. It is expected to return an error code.
 * \member destroy This function is called whenever the graph object
 *    is destroyed, right before freeing the allocated memory.
 * \member copy This function is called when copying a graph with \ref
 *    igraph_copy, after the structure of the graph has been already
 *    copied. It is expected to return an error code.
 * \member add_vertices Called when vertices are added to a
 *    graph, before adding the vertices themselves.
 *    The number of vertices to add is supplied as an
 *    argument. Expected to return an error code.
 * \member permute_vertices Typically called when a new graph is
 *    created based on an existing one, e.g. if vertices are removed
 *    from a graph. The supplied index vector defines which old vertex
 *    a new vertex corresponds to. Its length must be the same as the
 *    number of vertices in the new graph.
 * \member combine_vertices This function is called when the creation
 *    of a new graph involves a merge (contraction, etc.) of vertices
 *    from another graph. The function is after the new graph was created.
 *    An argument specifies how several vertices from the old graph map to a
 *    single vertex in the new graph.
 * \member add_edges Called when new edges have been added. The number
 *    of new edges are supplied as well. It is expected to return an
 *    error code.
 * \member permute_edges Typically called when a new graph is created and
 *    some of the new edges should carry the attributes of some of the
 *    old edges. The idx vector shows the mapping between the old edges and
 *    the new ones. Its length is the same as the number of edges in the new
 *    graph, and for each edge it gives the id of the old edge (the edge in
 *    the old graph).
 * \member combine_edges This function is called when the creation
 *    of a new graph involves a merge (contraction, etc.) of edges
 *    from another graph. The function is after the new graph was created.
 *    An argument specifies how several edges from the old graph map to a
 *    single edge in the new graph.
 * \member get_info Query the attributes of a graph, the names and
 *    types should be returned.
 * \member has_attr Check whether a graph has the named
 *    graph/vertex/edge attribute.
 * \member gettype Query the type of a graph/vertex/edge attribute.
 * \member get_numeric_graph_attr Query a numeric graph attribute. The
 *    value should be placed as the first element of the \p value
 *    vector.
 * \member get_string_graph_attr Query a string graph attribute. The
 *    value should be placed as the first element of the \p value
 *    string vector.
 * \member get_bool_graph_attr Query a boolean graph attribute. The
 *    value should be placed as the first element of the \p value
 *    boolean vector.
 * \member get_numeric_vertex_attr Query a numeric vertex attribute,
 *    for the vertices included in \p vs.
 * \member get_string_vertex_attr Query a string vertex attribute,
 *    for the vertices included in \p vs.
 * \member get_bool_vertex_attr Query a boolean vertex attribute,
 *    for the vertices included in \p vs.
 * \member get_numeric_edge_attr Query a numeric edge attribute, for
 *    the edges included in \p es.
 * \member get_string_edge_attr Query a string edge attribute, for the
 *    edges included in \p es.
 * \member get_bool_edge_attr Query a boolean edge attribute, for the
 *    edges included in \p es.
 *
 * Note that the <function>get_*_*_attr</function> are allowed to
 * convert the attributes to numeric or string. E.g. if a vertex attribute
 * is a GNU R complex data type, then
 * <function>get_string_vertex_attribute</function> may serialize it
 * into a string, but this probably makes sense only if
 * <function>add_vertices</function> is able to deserialize it.
 */

typedef struct igraph_attribute_table_t {
    int (*init)(igraph_t *graph, igraph_vector_ptr_t *attr);
    void (*destroy)(igraph_t *graph);
    int (*copy)(igraph_t *to, const igraph_t *from, igraph_bool_t ga,
                igraph_bool_t va, igraph_bool_t ea);
    int (*add_vertices)(igraph_t *graph, long int nv, igraph_vector_ptr_t *attr);
    int (*permute_vertices)(const igraph_t *graph,
                            igraph_t *newgraph,
                            const igraph_vector_t *idx);
    int (*combine_vertices)(const igraph_t *graph,
                            igraph_t *newgraph,
                            const igraph_vector_ptr_t *merges,
                            const igraph_attribute_combination_t *comb);
    int (*add_edges)(igraph_t *graph, const igraph_vector_t *edges,
                     igraph_vector_ptr_t *attr);
    int (*permute_edges)(const igraph_t *graph,
                         igraph_t *newgraph, const igraph_vector_t *idx);
    int (*combine_edges)(const igraph_t *graph,
                         igraph_t *newgraph,
                         const igraph_vector_ptr_t *merges,
                         const igraph_attribute_combination_t *comb);
    int (*get_info)(const igraph_t *graph,
                    igraph_strvector_t *gnames, igraph_vector_t *gtypes,
                    igraph_strvector_t *vnames, igraph_vector_t *vtypes,
                    igraph_strvector_t *enames, igraph_vector_t *etypes);
    igraph_bool_t (*has_attr)(const igraph_t *graph, igraph_attribute_elemtype_t type,
                              const char *name);
    int (*gettype)(const igraph_t *graph, igraph_attribute_type_t *type,
                   igraph_attribute_elemtype_t elemtype, const char *name);
    int (*get_numeric_graph_attr)(const igraph_t *graph, const char *name,
                                  igraph_vector_t *value);
    int (*get_string_graph_attr)(const igraph_t *graph, const char *name,
                                 igraph_strvector_t *value);
    int (*get_bool_graph_attr)(const igraph_t *igraph, const char *name,
                               igraph_vector_bool_t *value);
    int (*get_numeric_vertex_attr)(const igraph_t *graph, const char *name,
                                   igraph_vs_t vs,
                                   igraph_vector_t *value);
    int (*get_string_vertex_attr)(const igraph_t *graph, const char *name,
                                  igraph_vs_t vs,
                                  igraph_strvector_t *value);
    int (*get_bool_vertex_attr)(const igraph_t *graph, const char *name,
                                igraph_vs_t vs,
                                igraph_vector_bool_t *value);
    int (*get_numeric_edge_attr)(const igraph_t *graph, const char *name,
                                 igraph_es_t es,
                                 igraph_vector_t *value);
    int (*get_string_edge_attr)(const igraph_t *graph, const char *name,
                                igraph_es_t es,
                                igraph_strvector_t *value);
    int (*get_bool_edge_attr)(const igraph_t *graph, const char *name,
                              igraph_es_t es,
                              igraph_vector_bool_t *value);
} igraph_attribute_table_t;

DECLDIR igraph_attribute_table_t * igraph_i_set_attribute_table(const igraph_attribute_table_t * table);

DECLDIR igraph_bool_t igraph_has_attribute_table(void);

#define IGRAPH_I_ATTRIBUTE_DESTROY(graph) \
    do {if ((graph)->attr) igraph_i_attribute_destroy(graph);} while(0)
#define IGRAPH_I_ATTRIBUTE_COPY(to,from,ga,va,ea) do { \
        int igraph_i_ret2=0; \
        if ((from)->attr) { \
            IGRAPH_CHECK(igraph_i_ret2=igraph_i_attribute_copy((to),(from),(ga),(va),(ea))); \
        } else { \
            (to)->attr = 0; \
        } \
        if (igraph_i_ret2 != 0) { \
            IGRAPH_ERROR("", igraph_i_ret2); \
        } \
    } while(0)

int igraph_i_attribute_init(igraph_t *graph, void *attr);
void igraph_i_attribute_destroy(igraph_t *graph);
int igraph_i_attribute_copy(igraph_t *to, const igraph_t *from,
                            igraph_bool_t ga, igraph_bool_t va, igraph_bool_t ea);
int igraph_i_attribute_add_vertices(igraph_t *graph, long int nv, void *attr);
int igraph_i_attribute_permute_vertices(const igraph_t *graph,
                                        igraph_t *newgraph,
                                        const igraph_vector_t *idx);
int igraph_i_attribute_combine_vertices(const igraph_t *graph,
                                        igraph_t *newgraph,
                                        const igraph_vector_ptr_t *merges,
                                        const igraph_attribute_combination_t *comb);
int igraph_i_attribute_add_edges(igraph_t *graph,
                                 const igraph_vector_t *edges, void *attr);
int igraph_i_attribute_permute_edges(const igraph_t *graph,
                                     igraph_t *newgraph,
                                     const igraph_vector_t *idx);
int igraph_i_attribute_combine_edges(const igraph_t *graph,
                                     igraph_t *newgraph,
                                     const igraph_vector_ptr_t *merges,
                                     const igraph_attribute_combination_t *comb);

int igraph_i_attribute_get_info(const igraph_t *graph,
                                igraph_strvector_t *gnames,
                                igraph_vector_t *gtypes,
                                igraph_strvector_t *vnames,
                                igraph_vector_t *vtypes,
                                igraph_strvector_t *enames,
                                igraph_vector_t *etypes);
igraph_bool_t igraph_i_attribute_has_attr(const igraph_t *graph,
        igraph_attribute_elemtype_t type,
        const char *name);
int igraph_i_attribute_gettype(const igraph_t *graph,
                               igraph_attribute_type_t *type,
                               igraph_attribute_elemtype_t elemtype,
                               const char *name);

int igraph_i_attribute_get_numeric_graph_attr(const igraph_t *graph,
        const char *name,
        igraph_vector_t *value);
int igraph_i_attribute_get_numeric_vertex_attr(const igraph_t *graph,
        const char *name,
        igraph_vs_t vs,
        igraph_vector_t *value);
int igraph_i_attribute_get_numeric_edge_attr(const igraph_t *graph,
        const char *name,
        igraph_es_t es,
        igraph_vector_t *value);
int igraph_i_attribute_get_string_graph_attr(const igraph_t *graph,
        const char *name,
        igraph_strvector_t *value);
int igraph_i_attribute_get_string_vertex_attr(const igraph_t *graph,
        const char *name,
        igraph_vs_t vs,
        igraph_strvector_t *value);
int igraph_i_attribute_get_string_edge_attr(const igraph_t *graph,
        const char *name,
        igraph_es_t es,
        igraph_strvector_t *value);
int igraph_i_attribute_get_bool_graph_attr(const igraph_t *graph,
        const char *name,
        igraph_vector_bool_t *value);
int igraph_i_attribute_get_bool_vertex_attr(const igraph_t *graph,
        const char *name,
        igraph_vs_t vs,
        igraph_vector_bool_t *value);
int igraph_i_attribute_get_bool_edge_attr(const igraph_t *graph,
        const char *name,
        igraph_es_t es,
        igraph_vector_bool_t *value);

/* Experimental attribute handler in C */

extern const igraph_attribute_table_t igraph_cattribute_table;

DECLDIR igraph_real_t igraph_cattribute_GAN(const igraph_t *graph, const char *name);
DECLDIR igraph_bool_t igraph_cattribute_GAB(const igraph_t *graph, const char *name);
DECLDIR const char* igraph_cattribute_GAS(const igraph_t *graph, const char *name);
DECLDIR igraph_real_t igraph_cattribute_VAN(const igraph_t *graph, const char *name,
        igraph_integer_t vid);
DECLDIR igraph_bool_t igraph_cattribute_VAB(const igraph_t *graph, const char *name,
        igraph_integer_t vid);
DECLDIR const char* igraph_cattribute_VAS(const igraph_t *graph, const char *name,
        igraph_integer_t vid);
DECLDIR igraph_real_t igraph_cattribute_EAN(const igraph_t *graph, const char *name,
        igraph_integer_t eid);
DECLDIR igraph_bool_t igraph_cattribute_EAB(const igraph_t *graph, const char *name,
        igraph_integer_t eid);
DECLDIR const char* igraph_cattribute_EAS(const igraph_t *graph, const char *name,
        igraph_integer_t eid);

DECLDIR int igraph_cattribute_VANV(const igraph_t *graph, const char *name,
                                   igraph_vs_t vids, igraph_vector_t *result);
DECLDIR int igraph_cattribute_EANV(const igraph_t *graph, const char *name,
                                   igraph_es_t eids, igraph_vector_t *result);
DECLDIR int igraph_cattribute_VASV(const igraph_t *graph, const char *name,
                                   igraph_vs_t vids, igraph_strvector_t *result);
DECLDIR int igraph_cattribute_EASV(const igraph_t *graph, const char *name,
                                   igraph_es_t eids, igraph_strvector_t *result);
DECLDIR int igraph_cattribute_VABV(const igraph_t *graph, const char *name,
                                   igraph_vs_t vids, igraph_vector_bool_t *result);
DECLDIR int igraph_cattribute_EABV(const igraph_t *graph, const char *name,
                                   igraph_es_t eids, igraph_vector_bool_t *result);

DECLDIR int igraph_cattribute_list(const igraph_t *graph,
                                   igraph_strvector_t *gnames, igraph_vector_t *gtypes,
                                   igraph_strvector_t *vnames, igraph_vector_t *vtypes,
                                   igraph_strvector_t *enames, igraph_vector_t *etypes);
DECLDIR igraph_bool_t igraph_cattribute_has_attr(const igraph_t *graph,
        igraph_attribute_elemtype_t type,
        const char *name);

DECLDIR int igraph_cattribute_GAN_set(igraph_t *graph, const char *name,
                                      igraph_real_t value);
DECLDIR int igraph_cattribute_GAB_set(igraph_t *graph, const char *name,
                                      igraph_bool_t value);
DECLDIR int igraph_cattribute_GAS_set(igraph_t *graph, const char *name,
                                      const char *value);
DECLDIR int igraph_cattribute_VAN_set(igraph_t *graph, const char *name,
                                      igraph_integer_t vid, igraph_real_t value);
DECLDIR int igraph_cattribute_VAB_set(igraph_t *graph, const char *name,
                                      igraph_integer_t vid, igraph_bool_t value);
DECLDIR int igraph_cattribute_VAS_set(igraph_t *graph, const char *name,
                                      igraph_integer_t vid, const char *value);
DECLDIR int igraph_cattribute_EAN_set(igraph_t *graph, const char *name,
                                      igraph_integer_t eid, igraph_real_t value);
DECLDIR int igraph_cattribute_EAB_set(igraph_t *graph, const char *name,
                                      igraph_integer_t eid, igraph_bool_t value);
DECLDIR int igraph_cattribute_EAS_set(igraph_t *graph, const char *name,
                                      igraph_integer_t eid, const char *value);

DECLDIR int igraph_cattribute_VAN_setv(igraph_t *graph, const char *name,
                                       const igraph_vector_t *v);
DECLDIR int igraph_cattribute_VAB_setv(igraph_t *graph, const char *name,
                                       const igraph_vector_bool_t *v);
DECLDIR int igraph_cattribute_VAS_setv(igraph_t *graph, const char *name,
                                       const igraph_strvector_t *sv);
DECLDIR int igraph_cattribute_EAN_setv(igraph_t *graph, const char *name,
                                       const igraph_vector_t *v);
DECLDIR int igraph_cattribute_EAB_setv(igraph_t *graph, const char *name,
                                       const igraph_vector_bool_t *v);
DECLDIR int igraph_cattribute_EAS_setv(igraph_t *graph, const char *name,
                                       const igraph_strvector_t *sv);

DECLDIR void igraph_cattribute_remove_g(igraph_t *graph, const char *name);
DECLDIR void igraph_cattribute_remove_v(igraph_t *graph, const char *name);
DECLDIR void igraph_cattribute_remove_e(igraph_t *graph, const char *name);
DECLDIR void igraph_cattribute_remove_all(igraph_t *graph, igraph_bool_t g,
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
 *  Set a numeric edge attribute for all vertices
 *
 * This is a shorthand for \ref igraph_cattribute_EAN_setv().
 * \param graph The graph.
 * \param n The name of the attribute.
 * \param v Vector containing the new values of the attributes.
 */
#define SETEANV(graph,n,v) (igraph_cattribute_EAN_setv((graph),(n),(v)))
/**
 * \define SETEABV
 *  Set a boolean edge attribute for all vertices
 *
 * This is a shorthand for \ref igraph_cattribute_EAB_setv().
 * \param graph The graph.
 * \param n The name of the attribute.
 * \param v Vector containing the new values of the attributes.
 */
#define SETEABV(graph,n,v) (igraph_cattribute_EAB_setv((graph),(n),(v)))
/**
 * \define SETEASV
 *  Set a string edge attribute for all vertices
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

__END_DECLS

#endif
