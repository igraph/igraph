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
   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
   02110-1301 USA

*/

#include "igraph_attributes.h"
#include "igraph_memory.h"

#include "graph/attributes.h"
#include "internal/hacks.h" /* strdup */

#include <string.h>
#include <stdarg.h>

/* Should you ever want to have a thread-local attribute handler table, prepend
 * IGRAPH_THREAD_LOCAL to the following declaration and #include "config.h". */
igraph_attribute_table_t *igraph_i_attribute_table = 0;

igraph_error_t igraph_i_attribute_init(igraph_t *graph, void *attr) {
    graph->attr = 0;
    if (igraph_i_attribute_table) {
        return igraph_i_attribute_table->init(graph, attr);
    } else {
        return IGRAPH_SUCCESS;
    }
}

void igraph_i_attribute_destroy(igraph_t *graph) {
    if (igraph_i_attribute_table) {
        igraph_i_attribute_table->destroy(graph);
    }
}

igraph_error_t igraph_i_attribute_copy(igraph_t *to, const igraph_t *from, igraph_bool_t ga,
                            igraph_bool_t va, igraph_bool_t ea) {
    if (igraph_i_attribute_table) {
        return igraph_i_attribute_table->copy(to, from, ga, va, ea);
    } else {
        return IGRAPH_SUCCESS;
    }
}

igraph_error_t igraph_i_attribute_add_vertices(igraph_t *graph, igraph_integer_t nv, void *attr) {
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

igraph_error_t igraph_i_attribute_add_edges(igraph_t *graph,
                                 const igraph_vector_int_t *edges, void *attr) {
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

igraph_error_t igraph_i_attribute_gettype(const igraph_t *graph,
                               igraph_attribute_type_t *type,
                               igraph_attribute_elemtype_t elemtype,
                               const char *name) {
    if (igraph_i_attribute_table) {
        return igraph_i_attribute_table->gettype(graph, type, elemtype, name);
    } else {
        return IGRAPH_SUCCESS;
    }

}

igraph_error_t igraph_i_attribute_get_numeric_graph_attr(const igraph_t *graph,
        const char *name,
        igraph_vector_t *value) {
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
    if (igraph_i_attribute_table) {
        return igraph_i_attribute_table->get_numeric_edge_attr(graph, name, es, value);
    } else {
        return IGRAPH_SUCCESS;
    }
}

igraph_error_t igraph_i_attribute_get_string_graph_attr(const igraph_t *graph,
        const char *name,
        igraph_strvector_t *value) {
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
    if (igraph_i_attribute_table) {
        return igraph_i_attribute_table->get_string_edge_attr(graph, name, es, value);
    } else {
        return IGRAPH_SUCCESS;
    }
}

igraph_error_t igraph_i_attribute_get_bool_graph_attr(const igraph_t *graph,
        const char *name,
        igraph_vector_bool_t *value) {
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

igraph_attribute_table_t *
igraph_i_set_attribute_table(const igraph_attribute_table_t * table) {
    IGRAPH_WARNING("igraph_i_set_attribute_table is deprecated, use igraph_set_attribute_table.");
    return igraph_set_attribute_table(table);
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
    igraph_integer_t i, n = igraph_vector_ptr_size(&comb->list);
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
    igraph_integer_t i, n = igraph_vector_ptr_size(&comb->list);

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
    igraph_integer_t i, n = igraph_vector_ptr_size(&comb->list);

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
    igraph_integer_t i, def = -1, len = igraph_vector_ptr_size(&comb->list);

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
