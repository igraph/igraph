/*
   IGraph library.
   Copyright (C) 2007-2022  The igraph development team

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

#include "igraph_memory.h"
#include "igraph_error.h"

#include "io/gml-tree.h"

#include <string.h>

igraph_error_t igraph_gml_tree_init_integer(igraph_gml_tree_t *t,
                                            const char *name,
                                            igraph_integer_t line,
                                            igraph_integer_t value) {

    igraph_integer_t *p;

    IGRAPH_VECTOR_PTR_INIT_FINALLY(&t->names, 1);
    IGRAPH_VECTOR_CHAR_INIT_FINALLY(&t->types, 1);
    IGRAPH_VECTOR_PTR_INIT_FINALLY(&t->children, 1);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&t->lines, 1);

    /* names */
    VECTOR(t->names)[0] = (void*) name;

    /* line number */
    VECTOR(t->lines)[0] = line;

    /* types */
    VECTOR(t->types)[0] = IGRAPH_I_GML_TREE_INTEGER;

    /* children */
    p = IGRAPH_CALLOC(1, igraph_integer_t);
    IGRAPH_CHECK_OOM(p, "Cannot create integer GML tree node.");
    *p = value;
    VECTOR(t->children)[0] = p;

    IGRAPH_FINALLY_CLEAN(4);
    return IGRAPH_SUCCESS;
}

igraph_error_t igraph_gml_tree_init_real(igraph_gml_tree_t *t,
                                         const char *name,
                                         igraph_integer_t line,
                                         igraph_real_t value) {

    igraph_real_t *p;

    IGRAPH_VECTOR_PTR_INIT_FINALLY(&t->names, 1);
    IGRAPH_VECTOR_CHAR_INIT_FINALLY(&t->types, 1);
    IGRAPH_VECTOR_PTR_INIT_FINALLY(&t->children, 1);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&t->lines, 1);

    /* names */
    VECTOR(t->names)[0] = (void*) name;

    /* line number */
    VECTOR(t->lines)[0] = line;

    /* types */
    VECTOR(t->types)[0] = IGRAPH_I_GML_TREE_REAL;

    /* children */
    p = IGRAPH_CALLOC(1, igraph_real_t);
    IGRAPH_CHECK_OOM(p, "Cannot create real GML tree node.");
    *p = value;
    VECTOR(t->children)[0] = p;

    IGRAPH_FINALLY_CLEAN(4);
    return IGRAPH_SUCCESS;
}

igraph_error_t igraph_gml_tree_init_string(igraph_gml_tree_t *t,
                                           const char *name,
                                           igraph_integer_t line,
                                           const char *value) {

    IGRAPH_VECTOR_PTR_INIT_FINALLY(&t->names, 1);
    IGRAPH_VECTOR_CHAR_INIT_FINALLY(&t->types, 1);
    IGRAPH_VECTOR_PTR_INIT_FINALLY(&t->children, 1);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&t->lines, 1);

    /* names */
    VECTOR(t->names)[0] = (void*) name;

    /* line number */
    VECTOR(t->lines)[0] = line;

    /* types */
    VECTOR(t->types)[0] = IGRAPH_I_GML_TREE_STRING;

    /* children */
    VECTOR(t->children)[0] = (void*) value;

    IGRAPH_FINALLY_CLEAN(4);
    return IGRAPH_SUCCESS;
}

igraph_error_t igraph_gml_tree_init_tree(igraph_gml_tree_t *t,
                                         const char *name,
                                         igraph_integer_t line,
                                         igraph_gml_tree_t *value) {

    IGRAPH_VECTOR_PTR_INIT_FINALLY(&t->names, 1);
    IGRAPH_VECTOR_CHAR_INIT_FINALLY(&t->types, 1);
    IGRAPH_VECTOR_PTR_INIT_FINALLY(&t->children, 1);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&t->lines, 1);

    /* names */
    VECTOR(t->names)[0] = (void*) name;

    /* line number */
    VECTOR(t->lines)[0] = line;

    /* types */
    VECTOR(t->types)[0] = IGRAPH_I_GML_TREE_TREE;

    /* children */
    VECTOR(t->children)[0] = value;

    IGRAPH_FINALLY_CLEAN(4);
    return IGRAPH_SUCCESS;

}

igraph_error_t igraph_gml_tree_init_empty(igraph_gml_tree_t *t) {
    IGRAPH_VECTOR_PTR_INIT_FINALLY(&t->names, 0);
    IGRAPH_VECTOR_CHAR_INIT_FINALLY(&t->types, 0);
    IGRAPH_VECTOR_PTR_INIT_FINALLY(&t->children, 0);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&t->lines, 0);
    IGRAPH_FINALLY_CLEAN(4);
    return IGRAPH_SUCCESS;
}

/* merge is destructive, the _second_ tree is destroyed */
igraph_error_t igraph_gml_tree_mergedest(igraph_gml_tree_t *t1, igraph_gml_tree_t *t2) {
    igraph_integer_t i, n = igraph_vector_ptr_size(&t2->children);

    for (i = 0; i < n; i++) {
        IGRAPH_CHECK(igraph_vector_ptr_push_back(&t1->names, VECTOR(t2->names)[i]));
        IGRAPH_CHECK(igraph_vector_char_push_back(&t1->types, VECTOR(t2->types)[i]));
        IGRAPH_CHECK(igraph_vector_ptr_push_back(&t1->children, VECTOR(t2->children)[i]));
        IGRAPH_CHECK(igraph_vector_int_push_back(&t1->lines, VECTOR(t2->lines)[i]));
    }

    igraph_vector_ptr_destroy(&t2->names);
    igraph_vector_char_destroy(&t2->types);
    igraph_vector_ptr_destroy(&t2->children);
    igraph_vector_int_destroy(&t2->lines);

    return IGRAPH_SUCCESS;
}

void igraph_gml_tree_destroy(igraph_gml_tree_t *t) {

    igraph_integer_t i, n = igraph_vector_ptr_size(&t->children);
    for (i = 0; i < n; i++) {
        igraph_i_gml_tree_type_t type = (igraph_i_gml_tree_type_t) VECTOR(t->types)[i];
        switch (type) {
        case IGRAPH_I_GML_TREE_TREE:
            igraph_gml_tree_destroy(VECTOR(t->children)[i]);
            IGRAPH_FREE(VECTOR(t->names)[i]);
            break;
        case IGRAPH_I_GML_TREE_INTEGER:
            IGRAPH_FREE(VECTOR(t->children)[i]);
            IGRAPH_FREE(VECTOR(t->names)[i]);
            break;
        case IGRAPH_I_GML_TREE_REAL:
            IGRAPH_FREE(VECTOR(t->children)[i]);
            IGRAPH_FREE(VECTOR(t->names)[i]);
            break;
        case IGRAPH_I_GML_TREE_STRING:
            IGRAPH_FREE(VECTOR(t->children)[i]);
            IGRAPH_FREE(VECTOR(t->names)[i]);
            break;
        case IGRAPH_I_GML_TREE_DELETED:
            break;
        }
    }
    igraph_vector_ptr_destroy(&t->names);
    igraph_vector_char_destroy(&t->types);
    igraph_vector_ptr_destroy(&t->children);
    igraph_vector_int_destroy(&t->lines);
    IGRAPH_FREE(t);
}

igraph_integer_t igraph_gml_tree_length(const igraph_gml_tree_t *t) {
    return igraph_vector_ptr_size(&t->names);
}

igraph_integer_t igraph_gml_tree_find(
    const igraph_gml_tree_t *t, const char *name, igraph_integer_t from
) {
    igraph_integer_t size = igraph_vector_ptr_size(&t->names);
    while ( from < size && (! VECTOR(t->names)[from] ||
                            strcmp(VECTOR(t->names)[from], name)) ) {
        from++;
    }

    if (from == size) {
        from = -1;
    }
    return from;
}

igraph_integer_t igraph_gml_tree_findback(
    const igraph_gml_tree_t *t, const char *name, igraph_integer_t from
) {
    while ( from >= 0 && (! VECTOR(t->names)[from] ||
                          strcmp(VECTOR(t->names)[from], name)) ) {
        from--;
    }

    return from;
}

igraph_i_gml_tree_type_t igraph_gml_tree_type(const igraph_gml_tree_t *t, igraph_integer_t pos) {
    return (igraph_i_gml_tree_type_t) VECTOR(t->types)[pos];
}

const char *igraph_gml_tree_name(const igraph_gml_tree_t *t, igraph_integer_t pos) {
    return VECTOR(t->names)[pos];
}

igraph_integer_t igraph_gml_tree_line(const igraph_gml_tree_t *t, igraph_integer_t pos) {
    return VECTOR(t->lines)[pos];
}

igraph_integer_t igraph_gml_tree_get_integer(const igraph_gml_tree_t *t,
                                             igraph_integer_t pos) {
    igraph_integer_t *i = VECTOR(t->children)[pos];
    return *i;
}

igraph_real_t igraph_gml_tree_get_real(const igraph_gml_tree_t *t,
                                       igraph_integer_t pos) {
    igraph_real_t *d = VECTOR(t->children)[pos];
    return *d;
}

const char *igraph_gml_tree_get_string(const igraph_gml_tree_t *t,
                                       igraph_integer_t pos) {
    const char *s = VECTOR(t->children)[pos];
    return s;
}

igraph_gml_tree_t *igraph_gml_tree_get_tree(const igraph_gml_tree_t *t,
                                            igraph_integer_t pos) {
    igraph_gml_tree_t *tree = VECTOR(t->children)[pos];
    return tree;
}

void igraph_gml_tree_delete(igraph_gml_tree_t *t, igraph_integer_t pos) {
    if (VECTOR(t->types)[pos] == IGRAPH_I_GML_TREE_TREE) {
        igraph_gml_tree_destroy(VECTOR(t->children)[pos]);
    }
    IGRAPH_FREE(VECTOR(t->names)[pos]);
    IGRAPH_FREE(VECTOR(t->children)[pos]);
    VECTOR(t->children)[pos] = 0;
    VECTOR(t->names)[pos] = 0;
    VECTOR(t->types)[pos] = IGRAPH_I_GML_TREE_DELETED;
}
