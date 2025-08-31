/*
   igraph library.
   Copyright (C) 2005-2012  Gabor Csardi <csardi.gabor@gmail.com>
   334 Harvard street, Cambridge, MA 02139 USA

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.

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
#include "igraph_interface.h"
#include "igraph_random.h"

#include <string.h>

/* An attribute is either a numeric vector (vector_t), a boolean vector
 * (vector_bool_t) or a string vector (strvector_t).
 * The attribute itself is stored in a struct igraph_attribute_record_t. There
 * is one such object for each attribute. The igraph_t has a pointer to an
 * igraph_i_cattribute_t, which contains three igraph_attribute_vector_list_t's,
 * each holding pointers to igraph_attribute_record_t objects. */

typedef struct igraph_i_cattributes_t {
    igraph_attribute_record_list_t gal;
    igraph_attribute_record_list_t val;
    igraph_attribute_record_list_t eal;
} igraph_i_cattributes_t;

/*
 * Finds an attribute with the given name in an attribute record list, and
 * returns the index of the record in the list, or -1 if there was no such
 * attribute.
 */
static igraph_int_t igraph_i_cattribute_find_index(
    const igraph_attribute_record_list_t *attrs, const char *name
) {
    igraph_int_t n = igraph_attribute_record_list_size(attrs);
    for (igraph_int_t i = 0; i < n; i++) {
        const igraph_attribute_record_t *rec = igraph_attribute_record_list_get_ptr(attrs, i);
        if (!strcmp(rec->name, name)) {
            return i;
        }
    }
    return -1;
}

/*
 * Finds an attribute with the given name in an attribute record list, and
 * returns the attribute record or NULL if there was no such attribute.
 * Optionally, the type of the attribute can be enforced; use
 * IGRAPH_ATTRIBUTE_UNSPECIFIED if you do not care about the type.
 */
static igraph_attribute_record_t* igraph_i_cattribute_find(
    igraph_attribute_record_list_t *attrs, const char *name,
    igraph_attribute_type_t type
) {
    igraph_int_t index = igraph_i_cattribute_find_index(attrs, name);
    igraph_attribute_record_t *rec;

    if (index >= 0) {
        rec = igraph_attribute_record_list_get_ptr(attrs, index);
        if (type == IGRAPH_ATTRIBUTE_UNSPECIFIED || type == rec->type) {
            return rec;
        }
    }

    return NULL;
}

/*
 * Same as igraph_i_cattribute_find(), but suitable for const attribute record
 * lists.
 */
static const igraph_attribute_record_t* igraph_i_cattribute_find_const(
    const igraph_attribute_record_list_t *attrs, const char *name,
    igraph_attribute_type_t type
) {
    igraph_int_t index = igraph_i_cattribute_find_index(attrs, name);
    const igraph_attribute_record_t *rec;

    if (index >= 0) {
        rec = igraph_attribute_record_list_get_ptr(attrs, index);
        if (type == IGRAPH_ATTRIBUTE_UNSPECIFIED || type == rec->type) {
            return rec;
        }
    }

    return NULL;
}

/*
 * Finds an attribute with the given name in an attribute record list, and
 * returns the attribute record in an output argument. Returns an error code
 * if there is no such attribute. Optionally, the type of the attribute can be
 * enforced; use IGRAPH_ATTRIBUTE_UNSPECIFIED if you do not care about the type.
 */
static igraph_error_t igraph_i_cattribute_find_or_return(
    igraph_attribute_record_list_t *attrs, const char *name,
    igraph_attribute_type_t type, igraph_attribute_record_t **ptr
) {
    igraph_attribute_record_t *rec;

    rec = igraph_i_cattribute_find(attrs, name, IGRAPH_ATTRIBUTE_UNSPECIFIED);
    if (!rec) {
        IGRAPH_ERRORF("Attribute '%s' does not exist.", IGRAPH_EINVAL, name);
    }

    if (type != IGRAPH_ATTRIBUTE_UNSPECIFIED) {
        IGRAPH_CHECK(igraph_attribute_record_check_type(rec, type));
    }

    if (ptr) {
        *ptr = rec;
    }

    return IGRAPH_SUCCESS;
}

/*
 * Finds an attribute with the given name in an attribute record list, and
 * returns the attribute record in an output argument. Creates a new attribute
 * with the given name if there is no such attribute. The type of the attribute
 * needs to be specified so we can create the appropriate value vector if needed.
 * You can specify a length; when the existing value vector for the attribute
 * is shorter than this length, it will be extended to ensure that it has at
 * least this many elements. You can pass 0 as the length if you do not want to
 * expand value vectors.
 */
static igraph_error_t igraph_i_cattribute_find_or_create(
    igraph_attribute_record_list_t *attrs,
    const char *name, igraph_attribute_type_t type,
    igraph_int_t length,
    igraph_attribute_record_t **ptr
) {
    igraph_attribute_record_t *rec;

    rec = igraph_i_cattribute_find(attrs, name, IGRAPH_ATTRIBUTE_UNSPECIFIED);
    if (rec) {
        if (type != IGRAPH_ATTRIBUTE_UNSPECIFIED) {
            IGRAPH_CHECK(igraph_attribute_record_check_type(rec, type));
        }
    } else {
        IGRAPH_CHECK(igraph_attribute_record_list_push_back_new(attrs, &rec));
        IGRAPH_CHECK(igraph_attribute_record_set_name(rec, name));
        IGRAPH_CHECK(igraph_attribute_record_set_type(rec, type));
    }

    if (length > 0 && igraph_attribute_record_size(rec) < length) {
        IGRAPH_CHECK(igraph_attribute_record_resize(rec, length));
    }

    if (ptr) {
        *ptr = rec;
    }

    return IGRAPH_SUCCESS;
}

/*
 * Restores attribute vector lengths to their original size after a failure.
 * This function assumes that none of the attribute vectors are shorter than origlen.
 * Some may be longer due to a partially completed size extension: these will be
 * shrunk to their original size.
 */
static void igraph_i_cattribute_revert_attribute_vector_sizes(
        igraph_attribute_record_list_t *attrlist, igraph_int_t origlen) {

    igraph_int_t no_of_attrs = igraph_attribute_record_list_size(attrlist);
    for (igraph_int_t i = 0; i < no_of_attrs; i++) {
        igraph_attribute_record_t *rec = igraph_attribute_record_list_get_ptr(attrlist, i);
        IGRAPH_ASSERT(igraph_attribute_record_size(rec) >= origlen);
        if (igraph_attribute_record_resize(rec, origlen) != IGRAPH_SUCCESS) {
            /* We should have succeeded for known attribute types because we
             * always shrink the vector so throw a fatal error if this happens */
            IGRAPH_FATAL("Unknown attribute type encountered.");
        }
    }
}

static igraph_error_t igraph_i_cattribute_init(
    igraph_t *graph, const igraph_attribute_record_list_t *attr
) {
    igraph_i_cattributes_t *nattr;

    nattr = IGRAPH_CALLOC(1, igraph_i_cattributes_t);
    IGRAPH_CHECK_OOM(nattr, "Insufficient memory to allocate attribute storage.");
    IGRAPH_FINALLY(igraph_free, nattr);

    if (attr) {
        IGRAPH_CHECK(igraph_attribute_record_list_init_copy(&nattr->gal, attr));
    } else {
        IGRAPH_CHECK(igraph_attribute_record_list_init(&nattr->gal, 0));
    }
    IGRAPH_FINALLY(igraph_attribute_record_list_destroy, &nattr->gal);

    IGRAPH_CHECK(igraph_attribute_record_list_init(&nattr->val, 0));
    IGRAPH_FINALLY(igraph_attribute_record_list_destroy, &nattr->val);

    IGRAPH_CHECK(igraph_attribute_record_list_init(&nattr->eal, 0));
    IGRAPH_FINALLY(igraph_attribute_record_list_destroy, &nattr->eal);

    graph->attr = nattr;
    IGRAPH_FINALLY_CLEAN(4);

    return IGRAPH_SUCCESS;
}

static void igraph_i_cattribute_destroy(igraph_t *graph) {
    igraph_i_cattributes_t *attr = graph->attr;
    igraph_attribute_record_list_destroy(&attr->eal);
    igraph_attribute_record_list_destroy(&attr->val);
    igraph_attribute_record_list_destroy(&attr->gal);
    IGRAPH_FREE(graph->attr); /* sets to NULL */
}

/* No reference counting here. If you use attributes in C you should
   know what you're doing. */

static igraph_error_t igraph_i_cattribute_copy(
    igraph_t *to, const igraph_t *from,
    igraph_bool_t ga, igraph_bool_t va, igraph_bool_t ea
) {
    igraph_i_cattributes_t *attrto, *attrfrom = from->attr;
    igraph_attribute_record_list_t *alto[3], *alfrom[3] = {
        &attrfrom->gal, &attrfrom->val, &attrfrom->eal
    };
    igraph_bool_t copy[3] = { ga, va, ea };

    attrto = IGRAPH_CALLOC(1, igraph_i_cattributes_t);
    IGRAPH_CHECK_OOM(attrto, "Insufficient memory to copy attributes.");
    IGRAPH_FINALLY(igraph_free, attrto);

    alto[0] = &attrto->gal;
    alto[1] = &attrto->val;
    alto[2] = &attrto->eal;

    for (igraph_int_t i = 0; i < 3; i++) {
        if (copy[i]) {
            IGRAPH_CHECK(igraph_attribute_record_list_init_copy(alto[i], alfrom[i]));
        } else {
            IGRAPH_CHECK(igraph_attribute_record_list_init(alto[i], 0));
        }
        IGRAPH_FINALLY(igraph_attribute_record_list_destroy, alto[i]);
    }

    to->attr = attrto;
    IGRAPH_FINALLY_CLEAN(4);

    return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_i_cattribute_add_vertices_or_edges_inner(
    igraph_attribute_record_list_t *val,
    igraph_int_t newlen, igraph_int_t nv,
    const igraph_attribute_record_list_t *nattr
) {
    igraph_int_t length;
    igraph_int_t nattrno = nattr == NULL ? 0 : igraph_attribute_record_list_size(nattr);
    igraph_int_t origlen = newlen - nv;

    IGRAPH_ASSERT(origlen >= 0);

    /* Find all the attributes that are newly added, and create new value vectors
     * for them in the original graph */
    for (igraph_int_t i = 0; i < nattrno; i++) {
        const igraph_attribute_record_t *nattr_entry = igraph_attribute_record_list_get_ptr(nattr, i);
        const char *nname = nattr_entry->name;
        IGRAPH_CHECK(igraph_i_cattribute_find_or_create(val, nname, nattr_entry->type, origlen, NULL));
    }

    /* Now append the new values */
    length = igraph_attribute_record_list_size(val);
    for (igraph_int_t i = 0; i < length; i++) {
        igraph_attribute_record_t *oldrec = igraph_attribute_record_list_get_ptr(val, i);
        const igraph_attribute_record_t *newrec = nattr
            ? igraph_i_cattribute_find_const(nattr, oldrec->name, oldrec->type)
            : NULL;

        IGRAPH_ASSERT(igraph_attribute_record_size(oldrec) == origlen);

        if (newrec) {
            /* This attribute is present in nattr */
            switch (oldrec->type) {
            case IGRAPH_ATTRIBUTE_NUMERIC:
                if (nv != igraph_vector_size(newrec->value.as_vector)) {
                    IGRAPH_ERROR("Invalid numeric attribute length.", IGRAPH_EINVAL);
                }
                IGRAPH_CHECK(igraph_vector_append(
                    oldrec->value.as_vector, newrec->value.as_vector
                ));
                break;
            case IGRAPH_ATTRIBUTE_STRING:
                if (nv != igraph_strvector_size(newrec->value.as_strvector)) {
                    IGRAPH_ERROR("Invalid string attribute length.", IGRAPH_EINVAL);
                }
                IGRAPH_CHECK(igraph_strvector_append(
                    oldrec->value.as_strvector, newrec->value.as_strvector
                ));
                break;
            case IGRAPH_ATTRIBUTE_BOOLEAN:
                if (nv != igraph_vector_bool_size(newrec->value.as_vector_bool)) {
                    IGRAPH_ERROR("Invalid boolean attribute length.", IGRAPH_EINVAL);
                }
                IGRAPH_CHECK(igraph_vector_bool_append(
                    oldrec->value.as_vector_bool, newrec->value.as_vector_bool
                ));
                break;
            default:
                IGRAPH_WARNINGF(
                    "Attribute '%s' with unknown type %d ignored",
                    oldrec->name, (int) oldrec->type
                );
                break;
            }
        } else {
            /* No such attribute among the new ones so just extend the length
             * of the current record */
            IGRAPH_CHECK(igraph_attribute_record_resize(oldrec, newlen));
        }

        IGRAPH_ASSERT(igraph_attribute_record_size(oldrec) == newlen);
    }

    return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_i_cattribute_add_vertices_or_edges(
    igraph_attribute_record_list_t *val,
    igraph_int_t newlen, igraph_int_t nv,
    const igraph_attribute_record_list_t *nattr
) {
    igraph_int_t origlen = newlen - nv;
    igraph_error_t err = igraph_i_cattribute_add_vertices_or_edges_inner(
        val, newlen, nv, nattr
    );

    if (err != IGRAPH_SUCCESS) {
        /* If unsuccessful, revert attribute vector sizes.
         * The following function assumes that all attributes vectors that
         * are present have a length at least as great as origlen.
         * This is true at the moment because any new attributes that are
         * added to the graph are created directly at 'origlen' instead of
         * being created at smaller sizes and resized later.
         *
         * NOTE: While this ensures that all attribute vector lengths are
         * correct, it does not ensure that no extra attributes have
         * been added to the graph. However, the presence of extra
         * attributes does not make the attribute table inconsistent
         * like the incorrect attribute vector lengths would.
         */
        igraph_i_cattribute_revert_attribute_vector_sizes(val, origlen);
    }

    return err;
}

static igraph_error_t igraph_i_cattribute_add_vertices(
    igraph_t *graph, igraph_int_t nv,
    const igraph_attribute_record_list_t *nattr
) {
    igraph_i_cattributes_t *attr = graph->attr;
    return igraph_i_cattribute_add_vertices_or_edges(&attr->val, igraph_vcount(graph), nv, nattr);
}

typedef struct {
    igraph_vector_t *numeric;
    igraph_vector_bool_t *boolean;
    igraph_vector_ptr_t *strings;
    igraph_int_t length;
} igraph_i_attribute_permutation_work_area_t;

static igraph_error_t igraph_i_attribute_permutation_work_area_init(
  igraph_i_attribute_permutation_work_area_t *work_area, igraph_int_t length
) {
    work_area->length = length;
    work_area->numeric = NULL;
    work_area->boolean = NULL;
    work_area->strings = NULL;
    return IGRAPH_SUCCESS;
}

static void igraph_i_attribute_permutation_work_area_release_stored_strvectors(
  igraph_i_attribute_permutation_work_area_t *work_area
) {
    if (work_area->strings != NULL) {
        igraph_vector_ptr_destroy_all(work_area->strings);
        IGRAPH_FREE(work_area->strings);
        work_area->strings = NULL;
    }
}

static void igraph_i_attribute_permutation_work_area_destroy(
  igraph_i_attribute_permutation_work_area_t *work_area
) {
    igraph_i_attribute_permutation_work_area_release_stored_strvectors(work_area);
    if (work_area->numeric != NULL) {
        igraph_vector_destroy(work_area->numeric);
        IGRAPH_FREE(work_area->numeric);
        work_area->numeric = NULL;
    }
    if (work_area->boolean != NULL) {
        igraph_vector_bool_destroy(work_area->boolean);
        IGRAPH_FREE(work_area->boolean);
        work_area->boolean = NULL;
    }
}

static igraph_error_t igraph_i_attribute_permutation_work_area_alloc_for_numeric(
  igraph_i_attribute_permutation_work_area_t *work_area
) {
    igraph_vector_t* vec = work_area->numeric;

    if (vec == NULL) {
        vec = IGRAPH_CALLOC(1, igraph_vector_t);
        IGRAPH_CHECK_OOM(vec, "Cannot permute attributes.");
        IGRAPH_FINALLY(igraph_free, vec);
        IGRAPH_CHECK(igraph_vector_init(vec, work_area->length));
        work_area->numeric = vec;
        IGRAPH_FINALLY_CLEAN(1);
    }

    return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_i_attribute_permutation_work_area_alloc_for_boolean(
  igraph_i_attribute_permutation_work_area_t *work_area
) {
    igraph_vector_bool_t* vec = work_area->boolean;

    if (vec == NULL) {
        vec = IGRAPH_CALLOC(1, igraph_vector_bool_t);
        IGRAPH_CHECK_OOM(vec, "Cannot permute attributes.");
        IGRAPH_FINALLY(igraph_free, vec);
        IGRAPH_CHECK(igraph_vector_bool_init(vec, work_area->length));
        work_area->boolean = vec;
        IGRAPH_FINALLY_CLEAN(1);
    }

    return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_i_attribute_permutation_work_area_alloc_for_strings(
  igraph_i_attribute_permutation_work_area_t *work_area
) {
    igraph_vector_ptr_t* vec = work_area->strings;

    if (vec == NULL) {
        vec = IGRAPH_CALLOC(1, igraph_vector_ptr_t);
        IGRAPH_CHECK_OOM(vec, "Cannot permute attributes.");
        IGRAPH_FINALLY(igraph_free, vec);
        IGRAPH_CHECK(igraph_vector_ptr_init(vec, 0));
        IGRAPH_VECTOR_PTR_SET_ITEM_DESTRUCTOR(vec, igraph_strvector_destroy);
        work_area->strings = vec;
        IGRAPH_FINALLY_CLEAN(1);
    }

    return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_i_attribute_permutation_work_area_permute_and_store_strvector(
  igraph_i_attribute_permutation_work_area_t *work_area,
  const igraph_strvector_t *vec,
  const igraph_vector_int_t *idx
) {
    igraph_strvector_t *new_vec;

    new_vec = IGRAPH_CALLOC(1, igraph_strvector_t);
    IGRAPH_CHECK_OOM(new_vec, "Cannot permute attributes.");
    IGRAPH_FINALLY(igraph_free, new_vec);
    IGRAPH_CHECK(igraph_strvector_init(new_vec, 0));
    IGRAPH_FINALLY(igraph_strvector_destroy, new_vec);
    IGRAPH_CHECK(igraph_vector_ptr_push_back(work_area->strings, new_vec));
    IGRAPH_FINALLY_CLEAN(2);

    IGRAPH_CHECK(igraph_strvector_index(vec, new_vec, idx));

    return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_i_cattribute_permute_attribute_record_list(
    igraph_attribute_record_list_t *attrs,
    igraph_attribute_record_list_t *new_attrs,
    const igraph_vector_int_t *idx
) {
    igraph_int_t no_attrs, idxlen;

    no_attrs = igraph_attribute_record_list_size(attrs);

    /* When vertices or edges  are permuted, we now assume that there are no
     * attributes in the target attribute list yet */
    IGRAPH_ASSERT(igraph_attribute_record_list_empty(new_attrs));
    IGRAPH_FINALLY(igraph_attribute_record_list_clear, new_attrs);

    idxlen = igraph_vector_int_size(idx);
    for (igraph_int_t i = 0; i < no_attrs; i++) {
        igraph_attribute_record_t *oldrec = igraph_attribute_record_list_get_ptr(attrs, i);
        igraph_attribute_type_t type = oldrec->type;

        /* Create a record for the same attribute in the new graph */
        igraph_attribute_record_t *newrec;
        IGRAPH_CHECK(igraph_i_cattribute_find_or_create(
            new_attrs, oldrec->name, oldrec->type, idxlen, &newrec
        ));

        /* The data */
        switch (type) {
        case IGRAPH_ATTRIBUTE_NUMERIC:
            IGRAPH_CHECK(igraph_vector_index(
                oldrec->value.as_vector, newrec->value.as_vector, idx
            ));
            break;
        case IGRAPH_ATTRIBUTE_BOOLEAN:
            IGRAPH_CHECK(igraph_vector_bool_index(
                oldrec->value.as_vector_bool, newrec->value.as_vector_bool, idx
            ));
            break;
        case IGRAPH_ATTRIBUTE_STRING:
            IGRAPH_CHECK(igraph_strvector_index(
                oldrec->value.as_strvector, newrec->value.as_strvector, idx
            ));
            break;
        default:
            IGRAPH_WARNINGF(
                "Attribute '%s' with unknown type %d ignored",
                oldrec->name, (int) oldrec->type
            );
        }
    }

    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_i_cattribute_permute_attribute_record_list_in_place(
    igraph_attribute_record_list_t *attrs, const igraph_vector_int_t *idx
) {
    igraph_int_t no_attrs = igraph_attribute_record_list_size(attrs);
    igraph_attribute_record_t *oldrec;
    igraph_vector_t *num, *num_work;
    igraph_strvector_t *str, str_work;
    igraph_vector_bool_t *oldbool, *bool_work;
    igraph_i_attribute_permutation_work_area_t work_area;
    igraph_int_t idx_size = igraph_vector_int_size(idx);

    /* shortcut: don't allocate anything if there are no attributes */
    if (no_attrs == 0) {
        return IGRAPH_SUCCESS;
    }

    /* do all the allocations that can potentially fail before we actually
     * start to permute the vertices to ensure that we will not ever need to
     * back out from a permutation once we've started it */
    IGRAPH_CHECK(igraph_i_attribute_permutation_work_area_init(&work_area, idx_size));
    IGRAPH_FINALLY(igraph_i_attribute_permutation_work_area_destroy, &work_area);
    for (igraph_int_t i = 0; i < no_attrs; i++) {
        oldrec = igraph_attribute_record_list_get_ptr(attrs, i);
        switch (oldrec->type) {
        case IGRAPH_ATTRIBUTE_NUMERIC:
            num = oldrec->value.as_vector;
            IGRAPH_CHECK(igraph_vector_reserve(num, idx_size));
            IGRAPH_CHECK(igraph_i_attribute_permutation_work_area_alloc_for_numeric(&work_area));
            break;

        case IGRAPH_ATTRIBUTE_BOOLEAN:
            oldbool = oldrec->value.as_vector_bool;
            IGRAPH_CHECK(igraph_vector_bool_reserve(oldbool, idx_size));
            IGRAPH_CHECK(igraph_i_attribute_permutation_work_area_alloc_for_boolean(&work_area));
            break;

        case IGRAPH_ATTRIBUTE_STRING:
            str = oldrec->value.as_strvector;
            IGRAPH_CHECK(igraph_strvector_reserve(str, idx_size));
            IGRAPH_CHECK(igraph_i_attribute_permutation_work_area_alloc_for_strings(&work_area));
            break;

        default:
            IGRAPH_WARNINGF(
                "Vertex attribute '%s' with unknown type %d ignored",
                oldrec->name, (int) oldrec->type
            );
        }
    }

    /* let's do string attributes first because these might need extra
     * allocations that can fail. The strategy is to build new igraph_strvector_t
     * instances for the permuted attributes and store them in an
     * igraph_vector_ptr_t until we are done with all of them. If any of the
     * allocations fail, we can destroy the igraph_vector_ptr_t safely */
    for (igraph_int_t i = 0; i < no_attrs; i++) {
        oldrec = igraph_attribute_record_list_get_ptr(attrs, i);
        if (oldrec->type == IGRAPH_ATTRIBUTE_STRING) {
            str = oldrec->value.as_strvector;
            IGRAPH_CHECK(
                igraph_i_attribute_permutation_work_area_permute_and_store_strvector(
                    &work_area, str, idx
                )
            );
        }
    }

    /* strings are done, and now all vectors involved in the process are
     * as large as they should be (or larger) so the operations below are not
     * supposed to fail. We can safely replace the original string attribute
     * vectors with the permuted ones, and then proceed to the remaining
     * attributes */
    for (igraph_int_t i = 0, j = 0; i < no_attrs; i++) {
        oldrec = igraph_attribute_record_list_get_ptr(attrs, i);
        if (oldrec->type != IGRAPH_ATTRIBUTE_STRING) {
            continue;
        }

        str = oldrec->value.as_strvector;
        str_work = *((igraph_strvector_t*) VECTOR(*(work_area.strings))[j]);
        *((igraph_strvector_t*) VECTOR(*(work_area.strings))[j]) = *str;
        *str = str_work;
        j++;
    }
    igraph_i_attribute_permutation_work_area_release_stored_strvectors(&work_area);

    for (igraph_int_t i = 0; i < no_attrs; i++) {
        oldrec = igraph_attribute_record_list_get_ptr(attrs, i);
        switch (oldrec->type) {
        case IGRAPH_ATTRIBUTE_NUMERIC:
            num = oldrec->value.as_vector;
            num_work = work_area.numeric;
            IGRAPH_ASSERT(num_work != NULL);
            IGRAPH_CHECK(igraph_vector_index(num, num_work, idx));
            igraph_vector_swap(num, num_work);
            break;
        case IGRAPH_ATTRIBUTE_BOOLEAN:
            oldbool = oldrec->value.as_vector_bool;
            bool_work = work_area.boolean;
            IGRAPH_ASSERT(bool_work != NULL);
            IGRAPH_CHECK(igraph_vector_bool_index(oldbool, bool_work, idx));
            igraph_vector_bool_swap(oldbool, bool_work);
            break;
        case IGRAPH_ATTRIBUTE_STRING:
            /* nothing to do */
            break;
        default:
            /* already warned */
            break;
        }
    }

    igraph_i_attribute_permutation_work_area_destroy(&work_area);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_i_cattribute_permute_vertices(
    const igraph_t *graph, igraph_t *newgraph, const igraph_vector_int_t *idx
) {
    igraph_i_cattributes_t *attr = graph->attr, *new_attr = newgraph->attr;
    igraph_attribute_record_list_t *val = &attr->val, *new_val = &new_attr->val;
    if (graph == newgraph) {
        return igraph_i_cattribute_permute_attribute_record_list_in_place(val, idx);
    } else {
        return igraph_i_cattribute_permute_attribute_record_list(val, new_val, idx);
    }
}

typedef igraph_error_t igraph_cattributes_combine_num_t(const igraph_vector_t *input,
        igraph_real_t *output);

typedef igraph_error_t igraph_cattributes_combine_str_t(const igraph_strvector_t *input,
        char **output);

typedef igraph_error_t igraph_cattributes_combine_bool_t(const igraph_vector_bool_t *input,
        igraph_bool_t *output);

static igraph_error_t igraph_i_cattributes_cn_sum(const igraph_attribute_record_t *oldrec,
                                       igraph_attribute_record_t * newrec,
                                       const igraph_vector_int_list_t *merges) {
    const igraph_vector_t *oldv = oldrec->value.as_vector;
    igraph_vector_t *newv = newrec->value.as_vector;
    igraph_int_t newlen = igraph_vector_int_list_size(merges);

    for (igraph_int_t i = 0; i < newlen; i++) {
        igraph_real_t s = 0.0;
        const igraph_vector_int_t *idx = igraph_vector_int_list_get_ptr(merges, i);
        igraph_int_t n = igraph_vector_int_size(idx);
        for (igraph_int_t j = 0; j < n; j++) {
            igraph_int_t x = VECTOR(*idx)[j];
            s += VECTOR(*oldv)[x];
        }
        VECTOR(*newv)[i] = s;
    }

    return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_i_cattributes_cn_prod(const igraph_attribute_record_t *oldrec,
                                        igraph_attribute_record_t * newrec,
                                        const igraph_vector_int_list_t *merges) {
    const igraph_vector_t *oldv = oldrec->value.as_vector;
    igraph_vector_t *newv = newrec->value.as_vector;
    igraph_int_t newlen = igraph_vector_int_list_size(merges);

    for (igraph_int_t i = 0; i < newlen; i++) {
        igraph_real_t s = 1.0;
        const igraph_vector_int_t *idx = igraph_vector_int_list_get_ptr(merges, i);
        igraph_int_t n = igraph_vector_int_size(idx);
        for (igraph_int_t j = 0; j < n; j++) {
            igraph_int_t x = VECTOR(*idx)[j];
            s *= VECTOR(*oldv)[x];
        }
        VECTOR(*newv)[i] = s;
    }

    return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_i_cattributes_cn_min(const igraph_attribute_record_t *oldrec,
                                       igraph_attribute_record_t * newrec,
                                       const igraph_vector_int_list_t *merges) {
    const igraph_vector_t *oldv = oldrec->value.as_vector;
    igraph_vector_t *newv = newrec->value.as_vector;
    igraph_int_t newlen = igraph_vector_int_list_size(merges);

    for (igraph_int_t i = 0; i < newlen; i++) {
        const igraph_vector_int_t *idx = igraph_vector_int_list_get_ptr(merges, i);
        igraph_int_t n = igraph_vector_int_size(idx);
        igraph_real_t m = n > 0 ? VECTOR(*oldv)[ VECTOR(*idx)[0] ] : IGRAPH_NAN;
        for (igraph_int_t j = 1; j < n; j++) {
            igraph_int_t x = VECTOR(*idx)[j];
            igraph_real_t val = VECTOR(*oldv)[x];
            if (val < m) {
                m = val;
            }
        }
        VECTOR(*newv)[i] = m;
    }

    return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_i_cattributes_cn_max(const igraph_attribute_record_t *oldrec,
                                       igraph_attribute_record_t * newrec,
                                       const igraph_vector_int_list_t *merges) {
    const igraph_vector_t *oldv = oldrec->value.as_vector;
    igraph_vector_t *newv = newrec->value.as_vector;
    igraph_int_t newlen = igraph_vector_int_list_size(merges);

    for (igraph_int_t i = 0; i < newlen; i++) {
        const igraph_vector_int_t *idx = igraph_vector_int_list_get_ptr(merges, i);
        igraph_int_t n = igraph_vector_int_size(idx);
        igraph_real_t m = n > 0 ? VECTOR(*oldv)[ VECTOR(*idx)[0] ] : IGRAPH_NAN;
        for (igraph_int_t j = 1; j < n; j++) {
            igraph_int_t x = VECTOR(*idx)[j];
            igraph_real_t val = VECTOR(*oldv)[x];
            if (val > m) {
                m = val;
            }
        }
        VECTOR(*newv)[i] = m;
    }

    return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_i_cattributes_cn_random(const igraph_attribute_record_t *oldrec,
                                          igraph_attribute_record_t * newrec,
                                          const igraph_vector_int_list_t *merges) {

    const igraph_vector_t *oldv = oldrec->value.as_vector;
    igraph_vector_t *newv = newrec->value.as_vector;
    igraph_int_t newlen = igraph_vector_int_list_size(merges);

    for (igraph_int_t i = 0; i < newlen; i++) {
        const igraph_vector_int_t *idx = igraph_vector_int_list_get_ptr(merges, i);
        igraph_int_t n = igraph_vector_int_size(idx);
        if (n == 0) {
            VECTOR(*newv)[i] = IGRAPH_NAN;
        } else if (n == 1) {
            VECTOR(*newv)[i] = VECTOR(*oldv)[ VECTOR(*idx)[0] ];
        } else {
            igraph_int_t r = RNG_INTEGER(0, n - 1);
            VECTOR(*newv)[i] = VECTOR(*oldv)[ VECTOR(*idx)[r] ];
        }
    }

    return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_i_cattributes_cn_first(const igraph_attribute_record_t *oldrec,
                                         igraph_attribute_record_t * newrec,
                                         const igraph_vector_int_list_t *merges) {

    const igraph_vector_t *oldv = oldrec->value.as_vector;
    igraph_vector_t *newv = newrec->value.as_vector;
    igraph_int_t newlen = igraph_vector_int_list_size(merges);

    for (igraph_int_t i = 0; i < newlen; i++) {
        const igraph_vector_int_t *idx = igraph_vector_int_list_get_ptr(merges, i);
        igraph_int_t n = igraph_vector_int_size(idx);
        if (n == 0) {
            VECTOR(*newv)[i] = IGRAPH_NAN;
        } else {
            VECTOR(*newv)[i] = VECTOR(*oldv)[ VECTOR(*idx)[0] ];
        }
    }

    return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_i_cattributes_cn_last(const igraph_attribute_record_t *oldrec,
                                        igraph_attribute_record_t * newrec,
                                        const igraph_vector_int_list_t *merges) {

    const igraph_vector_t *oldv = oldrec->value.as_vector;
    igraph_vector_t *newv = newrec->value.as_vector;
    igraph_int_t newlen = igraph_vector_int_list_size(merges);

    for (igraph_int_t i = 0; i < newlen; i++) {
        const igraph_vector_int_t *idx = igraph_vector_int_list_get_ptr(merges, i);
        igraph_int_t n = igraph_vector_int_size(idx);
        if (n == 0) {
            VECTOR(*newv)[i] = IGRAPH_NAN;
        } else {
            VECTOR(*newv)[i] = VECTOR(*oldv)[ VECTOR(*idx)[n - 1] ];
        }
    }

    return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_i_cattributes_cn_mean(const igraph_attribute_record_t *oldrec,
                                        igraph_attribute_record_t * newrec,
                                        const igraph_vector_int_list_t *merges) {
    const igraph_vector_t *oldv = oldrec->value.as_vector;
    igraph_vector_t *newv = newrec->value.as_vector;
    igraph_int_t newlen = igraph_vector_int_list_size(merges);

    for (igraph_int_t i = 0; i < newlen; i++) {
        const igraph_vector_int_t *idx = igraph_vector_int_list_get_ptr(merges, i);
        igraph_int_t n = igraph_vector_int_size(idx);
        igraph_real_t s = n > 0 ? 0.0 : IGRAPH_NAN;
        for (igraph_int_t j = 0; j < n; j++) {
            igraph_int_t x = VECTOR(*idx)[j];
            s += VECTOR(*oldv)[x];
        }
        if (n > 0) {
            s = s / n;
        }
        VECTOR(*newv)[i] = s;
    }

    return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_i_cattributes_cn_func(const igraph_attribute_record_t *oldrec,
                                        igraph_attribute_record_t *newrec,
                                        const igraph_vector_int_list_t *merges,
                                        igraph_cattributes_combine_num_t *func) {

    const igraph_vector_t *oldv = oldrec->value.as_vector;
    igraph_vector_t *newv = newrec->value.as_vector;
    igraph_int_t newlen = igraph_vector_int_list_size(merges);

    igraph_vector_t values;
    IGRAPH_VECTOR_INIT_FINALLY(&values, 0);

    for (igraph_int_t i = 0; i < newlen; i++) {
        const igraph_vector_int_t *idx = igraph_vector_int_list_get_ptr(merges, i);

        igraph_int_t n = igraph_vector_int_size(idx);
        IGRAPH_CHECK(igraph_vector_resize(&values, n));
        for (igraph_int_t j = 0; j < n; j++) {
            igraph_int_t x = VECTOR(*idx)[j];
            VECTOR(values)[j] = VECTOR(*oldv)[x];
        }

        igraph_real_t res;
        IGRAPH_CHECK(func(&values, &res));
        VECTOR(*newv)[i] = res;
    }

    igraph_vector_destroy(&values);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_i_cattributes_cb_random(const igraph_attribute_record_t *oldrec,
                                          igraph_attribute_record_t * newrec,
                                          const igraph_vector_int_list_t *merges) {

    const igraph_vector_bool_t *oldv = oldrec->value.as_vector_bool;
    igraph_vector_bool_t *newv = newrec->value.as_vector_bool;
    igraph_int_t newlen = igraph_vector_int_list_size(merges);

    for (igraph_int_t i = 0; i < newlen; i++) {
        const igraph_vector_int_t *idx = igraph_vector_int_list_get_ptr(merges, i);
        igraph_int_t n = igraph_vector_int_size(idx);
        if (n == 0) {
            VECTOR(*newv)[i] = 0;
        } else if (n == 1) {
            VECTOR(*newv)[i] = VECTOR(*oldv)[ VECTOR(*idx)[0] ];
        } else {
            igraph_int_t r = RNG_INTEGER(0, n - 1);
            VECTOR(*newv)[i] = VECTOR(*oldv)[ VECTOR(*idx)[r] ];
        }
    }

    return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_i_cattributes_cb_first(const igraph_attribute_record_t *oldrec,
                                         igraph_attribute_record_t * newrec,
                                         const igraph_vector_int_list_t *merges) {

    const igraph_vector_bool_t *oldv = oldrec->value.as_vector_bool;
    igraph_vector_bool_t *newv = newrec->value.as_vector_bool;
    igraph_int_t newlen = igraph_vector_int_list_size(merges);

    for (igraph_int_t i = 0; i < newlen; i++) {
        const igraph_vector_int_t *idx = igraph_vector_int_list_get_ptr(merges, i);
        igraph_int_t n = igraph_vector_int_size(idx);
        if (n == 0) {
            VECTOR(*newv)[i] = 0;
        } else {
            VECTOR(*newv)[i] = VECTOR(*oldv)[ VECTOR(*idx)[0] ];
        }
    }

    return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_i_cattributes_cb_last(const igraph_attribute_record_t *oldrec,
                                        igraph_attribute_record_t * newrec,
                                        const igraph_vector_int_list_t *merges) {

    const igraph_vector_bool_t *oldv = oldrec->value.as_vector_bool;
    igraph_vector_bool_t *newv = newrec->value.as_vector_bool;
    igraph_int_t newlen = igraph_vector_int_list_size(merges);

    for (igraph_int_t i = 0; i < newlen; i++) {
        const igraph_vector_int_t *idx = igraph_vector_int_list_get_ptr(merges, i);
        igraph_int_t n = igraph_vector_int_size(idx);
        if (n == 0) {
            VECTOR(*newv)[i] = 0;
        } else {
            VECTOR(*newv)[i] = VECTOR(*oldv)[ VECTOR(*idx)[n - 1] ];
        }
    }

    return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_i_cattributes_cb_all_is_true(const igraph_attribute_record_t *oldrec,
                                               igraph_attribute_record_t * newrec,
                                               const igraph_vector_int_list_t *merges) {

    const igraph_vector_bool_t *oldv = oldrec->value.as_vector_bool;
    igraph_vector_bool_t *newv = newrec->value.as_vector_bool;
    igraph_int_t newlen = igraph_vector_int_list_size(merges);

    for (igraph_int_t i = 0; i < newlen; i++) {
        const igraph_vector_int_t *idx = igraph_vector_int_list_get_ptr(merges, i);
        igraph_int_t n = igraph_vector_int_size(idx);
        VECTOR(*newv)[i] = 1;
        for (igraph_int_t j = 0; j < n; j++) {
            igraph_int_t x = VECTOR(*idx)[j];
            if (!VECTOR(*oldv)[x]) {
                VECTOR(*newv)[i] = 0;
                break;
            }
        }
    }

    return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_i_cattributes_cb_any_is_true(const igraph_attribute_record_t *oldrec,
                                               igraph_attribute_record_t * newrec,
                                               const igraph_vector_int_list_t *merges) {

    const igraph_vector_bool_t *oldv = oldrec->value.as_vector_bool;
    igraph_vector_bool_t *newv = newrec->value.as_vector_bool;
    igraph_int_t newlen = igraph_vector_int_list_size(merges);

    for (igraph_int_t i = 0; i < newlen; i++) {
        const igraph_vector_int_t *idx = igraph_vector_int_list_get_ptr(merges, i);
        igraph_int_t n = igraph_vector_int_size(idx);
        VECTOR(*newv)[i] = 0;
        for (igraph_int_t j = 0; j < n; j++) {
            igraph_int_t x = VECTOR(*idx)[j];
            if (VECTOR(*oldv)[x]) {
                VECTOR(*newv)[i] = 1;
                break;
            }
        }
    }

    return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_i_cattributes_cb_majority(const igraph_attribute_record_t *oldrec,
                                            igraph_attribute_record_t * newrec,
                                            const igraph_vector_int_list_t *merges) {

    const igraph_vector_bool_t *oldv = oldrec->value.as_vector_bool;
    igraph_vector_bool_t *newv = newrec->value.as_vector_bool;
    igraph_int_t newlen = igraph_vector_int_list_size(merges);

    for (igraph_int_t i = 0; i < newlen; i++) {
        const igraph_vector_int_t *idx = igraph_vector_int_list_get_ptr(merges, i);
        igraph_int_t n = igraph_vector_int_size(idx);

        igraph_int_t num_trues = 0;
        for (igraph_int_t j = 0; j < n; j++) {
            igraph_int_t x = VECTOR(*idx)[j];
            if (VECTOR(*oldv)[x]) {
                num_trues++;
            }
        }

        if (n % 2 != 0) {
            VECTOR(*newv)[i] = (num_trues > n / 2);
        } else {
            if (num_trues == n / 2) {
                VECTOR(*newv)[i] = RNG_BOOL();
            } else {
                VECTOR(*newv)[i] = (num_trues > n / 2);
            }
        }
    }

    return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_i_cattributes_cb_func(const igraph_attribute_record_t *oldrec,
                                        igraph_attribute_record_t *newrec,
                                        const igraph_vector_int_list_t *merges,
                                        igraph_cattributes_combine_bool_t *func) {

    const igraph_vector_bool_t *oldv = oldrec->value.as_vector_bool;
    igraph_vector_bool_t *newv = newrec->value.as_vector_bool;
    igraph_int_t newlen = igraph_vector_int_list_size(merges);

    igraph_vector_bool_t values;
    IGRAPH_VECTOR_BOOL_INIT_FINALLY(&values, 0);

    for (igraph_int_t i = 0; i < newlen; i++) {
        const igraph_vector_int_t *idx = igraph_vector_int_list_get_ptr(merges, i);

        igraph_int_t n = igraph_vector_int_size(idx);
        IGRAPH_CHECK(igraph_vector_bool_resize(&values, n));
        for (igraph_int_t j = 0; j < n; j++) {
            igraph_int_t x = VECTOR(*idx)[j];
            VECTOR(values)[j] = VECTOR(*oldv)[x];
        }

        igraph_bool_t res;
        IGRAPH_CHECK(func(&values, &res));
        VECTOR(*newv)[i] = res;
    }

    igraph_vector_bool_destroy(&values);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_i_cattributes_sn_random(const igraph_attribute_record_t *oldrec,
                                          igraph_attribute_record_t *newrec,
                                          const igraph_vector_int_list_t *merges) {

    const igraph_strvector_t *oldv = oldrec->value.as_strvector;
    igraph_strvector_t *newv = newrec->value.as_strvector;
    igraph_int_t newlen = igraph_vector_int_list_size(merges);

    for (igraph_int_t i = 0; i < newlen; i++) {
        const igraph_vector_int_t *idx = igraph_vector_int_list_get_ptr(merges, i);
        igraph_int_t n = igraph_vector_int_size(idx);
        const char *tmp;
        if (n == 0) {
            IGRAPH_CHECK(igraph_strvector_set(newv, i, ""));
        } else if (n == 1) {
            tmp = igraph_strvector_get(oldv, 0);
            IGRAPH_CHECK(igraph_strvector_set(newv, i, tmp));
        } else {
            igraph_int_t r = RNG_INTEGER(0, n - 1);
            tmp = igraph_strvector_get(oldv, r);
            IGRAPH_CHECK(igraph_strvector_set(newv, i, tmp));
        }
    }

    return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_i_cattributes_cs_first(const igraph_attribute_record_t *oldrec,
                                         igraph_attribute_record_t *newrec,
                                         const igraph_vector_int_list_t *merges) {

    const igraph_strvector_t *oldv = oldrec->value.as_strvector;
    igraph_int_t newlen = igraph_vector_int_list_size(merges);
    igraph_strvector_t *newv = newrec->value.as_strvector;

    for (igraph_int_t i = 0; i < newlen; i++) {
        const igraph_vector_int_t *idx = igraph_vector_int_list_get_ptr(merges, i);
        igraph_int_t n = igraph_vector_int_size(idx);
        if (n == 0) {
            IGRAPH_CHECK(igraph_strvector_set(newv, i, ""));
        } else {
            const char *tmp = igraph_strvector_get(oldv, VECTOR(*idx)[0]);
            IGRAPH_CHECK(igraph_strvector_set(newv, i, tmp));
        }
    }

    return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_i_cattributes_cs_last(const igraph_attribute_record_t *oldrec,
                                        igraph_attribute_record_t *newrec,
                                        const igraph_vector_int_list_t *merges) {

    const igraph_strvector_t *oldv = oldrec->value.as_strvector;
    igraph_int_t newlen = igraph_vector_int_list_size(merges);
    igraph_strvector_t *newv = newrec->value.as_strvector;

    for (igraph_int_t i = 0; i < newlen; i++) {
        const igraph_vector_int_t *idx = igraph_vector_int_list_get_ptr(merges, i);
        igraph_int_t n = igraph_vector_int_size(idx);
        if (n == 0) {
            IGRAPH_CHECK(igraph_strvector_set(newv, i, ""));
        } else {
            const char *tmp = igraph_strvector_get(oldv, VECTOR(*idx)[n - 1]);
            IGRAPH_CHECK(igraph_strvector_set(newv, i, tmp));
        }
    }

    return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_i_cattributes_cs_concat(const igraph_attribute_record_t *oldrec,
                                          igraph_attribute_record_t *newrec,
                                          const igraph_vector_int_list_t *merges) {

    const igraph_strvector_t *oldv = oldrec->value.as_strvector;
    igraph_int_t newlen = igraph_vector_int_list_size(merges);
    igraph_strvector_t *newv = newrec->value.as_strvector;

    for (igraph_int_t i = 0; i < newlen; i++) {
        const igraph_vector_int_t *idx = igraph_vector_int_list_get_ptr(merges, i);
        igraph_int_t n = igraph_vector_int_size(idx);
        size_t len = 0;
        const char *tmp;
        char *tmp2;
        for (igraph_int_t j = 0; j < n; j++) {
            tmp = igraph_strvector_get(oldv, j);
            len += strlen(tmp);
        }
        tmp2 = IGRAPH_CALLOC(len + 1, char);
        IGRAPH_CHECK_OOM(tmp2, "Cannot combine attributes.");
        IGRAPH_FINALLY(igraph_free, tmp2);
        len = 0;
        for (igraph_int_t j = 0; j < n; j++) {
            tmp = igraph_strvector_get(oldv, j);
            strcpy(tmp2 + len, tmp);
            len += strlen(tmp);
        }

        IGRAPH_CHECK(igraph_strvector_set(newv, i, tmp2));
        IGRAPH_FREE(tmp2);
        IGRAPH_FINALLY_CLEAN(1);
    }

    return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_i_cattributes_cs_func(const igraph_attribute_record_t *oldrec,
                                        igraph_attribute_record_t *newrec,
                                        const igraph_vector_int_list_t *merges,
                                        igraph_cattributes_combine_str_t *func) {

    const igraph_strvector_t *oldv = oldrec->value.as_strvector;
    igraph_strvector_t *newv = newrec->value.as_strvector;
    igraph_int_t newlen = igraph_vector_int_list_size(merges);

    igraph_strvector_t values;
    IGRAPH_STRVECTOR_INIT_FINALLY(&values, 0);

    for (igraph_int_t i = 0; i < newlen; i++) {
        const igraph_vector_int_t *idx = igraph_vector_int_list_get_ptr(merges, i);

        igraph_int_t n = igraph_vector_int_size(idx);
        IGRAPH_CHECK(igraph_strvector_resize(&values, n));
        for (igraph_int_t j = 0; j < n; j++) {
            igraph_int_t x = VECTOR(*idx)[j];
            const char *elem = igraph_strvector_get(oldv, x);
            IGRAPH_CHECK(igraph_strvector_set(newv, j, elem));
        }

        char *res;
        IGRAPH_CHECK(func(&values, &res));
        IGRAPH_FINALLY(igraph_free, res);

        IGRAPH_CHECK(igraph_strvector_set(newv, i, res));

        IGRAPH_FREE(res);
        IGRAPH_FINALLY_CLEAN(1);
    }

    igraph_strvector_destroy(&values);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}

/**
 * \section c_attribute_combination_functions
 *
 * <para>
 * The C attribute handler supports combining the attributes of multiple
 * vertices of edges into a single attribute during a vertex or edge contraction
 * operation via a user-defined function. This is achieved by setting the
 * type of the attribute combination to \c IGRAPH_ATTRIBUTE_COMBINE_FUNCTION
 * and passing in a pointer to the custom combination function when specifying
 * attribute combinations in \ref igraph_attribute_combination() or
 * \ref igraph_attribute_combination_add() . For the C attribute handler, the
 * signature of the function depends on the type of the underlying attribute.
 * For numeric attributes, use:
 * \verbatim igraph_error_t function(const igraph_vector_t *input, igraph_real_t *output); \endverbatim
 * where \p input will receive a vector containing the value of the attribute
 * for all the vertices or edges being combined, and \p output must be filled
 * by the function to the combined value. Similarly, for Boolean attributes, the
 * function takes a boolean vector in \p input and must return the combined Boolean
 * value in \p output:
 * \verbatim igraph_error_t function(const igraph_vector_bool_t *input, igraph_bool_t *output); \endverbatim
 * For string attributes, the signature is slightly different:
 * \verbatim igraph_error_t function(const igraph_strvector_t *input, char **output); \endverbatim
 * In case of strings, all strings in the input vector are \em owned by igraph
 * and must not be modified or freed in the combination handler. The string
 * returned to the caller in \p output remains owned by the caller; igraph will
 * make a copy it and store the copy in the appropriate part of the data
 * structure holding the vertex or edge attributes.
 * </para>
 */

typedef struct {
    igraph_attribute_combination_type_t type;
    union {
        igraph_function_pointer_t as_void;
        igraph_cattributes_combine_num_t *as_num;
        igraph_cattributes_combine_str_t *as_str;
        igraph_cattributes_combine_bool_t *as_bool;
    } func;
} igraph_attribute_combination_todo_item_t;

static igraph_error_t igraph_i_cattribute_combine_attribute_record_lists(
    igraph_attribute_record_list_t *attrs, igraph_attribute_record_list_t *new_attrs,
    const igraph_vector_int_list_t *merges, const igraph_attribute_combination_t *comb
) {
    igraph_int_t no_attrs = igraph_attribute_record_list_size(attrs);
    igraph_attribute_combination_todo_item_t *todo_items;

    IGRAPH_ASSERT(attrs != new_attrs);
    IGRAPH_ASSERT(igraph_attribute_record_list_empty(new_attrs));

    todo_items = IGRAPH_CALLOC(no_attrs, igraph_attribute_combination_todo_item_t);
    IGRAPH_CHECK_OOM(todo_items, "Cannot combine attributes.");
    IGRAPH_FINALLY(igraph_free, todo_items);

    for (igraph_int_t i = 0; i < no_attrs; i++) {
        const igraph_attribute_record_t *oldrec = igraph_attribute_record_list_get_ptr(attrs, i);
        const char *name = oldrec->name;
        igraph_attribute_combination_type_t type;
        igraph_function_pointer_t voidfunc;
        IGRAPH_CHECK(igraph_attribute_combination_query(comb, name, &type, &voidfunc));
        todo_items[i].type = type;
        todo_items[i].func.as_void = voidfunc;
    }

    IGRAPH_FINALLY(igraph_attribute_record_list_clear, new_attrs);

    for (igraph_int_t i = 0; i < no_attrs; i++) {
        const igraph_attribute_record_t *oldrec = igraph_attribute_record_list_get_ptr(attrs, i);
        igraph_attribute_record_t newrec;
        const char *name = oldrec->name;
        igraph_attribute_combination_todo_item_t todo_item = todo_items[i];
        igraph_attribute_type_t attr_type = oldrec->type;

        if (todo_item.type == IGRAPH_ATTRIBUTE_COMBINE_DEFAULT ||
            todo_item.type == IGRAPH_ATTRIBUTE_COMBINE_IGNORE) {
            continue;
        }

        IGRAPH_CHECK(igraph_attribute_record_init(&newrec, name, attr_type));
        IGRAPH_FINALLY(igraph_attribute_record_destroy, &newrec);

        IGRAPH_CHECK(igraph_attribute_record_resize(&newrec, igraph_vector_int_list_size(merges)));

        if (attr_type == IGRAPH_ATTRIBUTE_NUMERIC) {
            switch (todo_item.type) {
            case IGRAPH_ATTRIBUTE_COMBINE_FUNCTION:
                IGRAPH_CHECK(igraph_i_cattributes_cn_func(oldrec, &newrec, merges,
                             todo_item.func.as_num));
                break;
            case IGRAPH_ATTRIBUTE_COMBINE_SUM:
                IGRAPH_CHECK(igraph_i_cattributes_cn_sum(oldrec, &newrec, merges));
                break;
            case IGRAPH_ATTRIBUTE_COMBINE_PROD:
                IGRAPH_CHECK(igraph_i_cattributes_cn_prod(oldrec, &newrec, merges));
                break;
            case IGRAPH_ATTRIBUTE_COMBINE_MIN:
                IGRAPH_CHECK(igraph_i_cattributes_cn_min(oldrec, &newrec, merges));
                break;
            case IGRAPH_ATTRIBUTE_COMBINE_MAX:
                IGRAPH_CHECK(igraph_i_cattributes_cn_max(oldrec, &newrec, merges));
                break;
            case IGRAPH_ATTRIBUTE_COMBINE_RANDOM:
                IGRAPH_CHECK(igraph_i_cattributes_cn_random(oldrec, &newrec, merges));
                break;
            case IGRAPH_ATTRIBUTE_COMBINE_FIRST:
                IGRAPH_CHECK(igraph_i_cattributes_cn_first(oldrec, &newrec, merges));
                break;
            case IGRAPH_ATTRIBUTE_COMBINE_LAST:
                IGRAPH_CHECK(igraph_i_cattributes_cn_last(oldrec, &newrec, merges));
                break;
            case IGRAPH_ATTRIBUTE_COMBINE_MEAN:
                IGRAPH_CHECK(igraph_i_cattributes_cn_mean(oldrec, &newrec, merges));
                break;
            case IGRAPH_ATTRIBUTE_COMBINE_MEDIAN:
                IGRAPH_ERROR("Median calculation not implemented.",
                             IGRAPH_UNIMPLEMENTED);
                break;
            case IGRAPH_ATTRIBUTE_COMBINE_CONCAT:
                IGRAPH_ERROR("Cannot concatenate numeric attributes.",
                             IGRAPH_EATTRCOMBINE);
                break;
            default:
                IGRAPH_ERROR("Unknown attribute combination.",
                             IGRAPH_UNIMPLEMENTED);
                break;
            }
        } else if (attr_type == IGRAPH_ATTRIBUTE_BOOLEAN) {
            switch (todo_item.type) {
            case IGRAPH_ATTRIBUTE_COMBINE_FUNCTION:
                IGRAPH_CHECK(igraph_i_cattributes_cb_func(oldrec, &newrec, merges,
                             todo_item.func.as_bool));
                break;
            case IGRAPH_ATTRIBUTE_COMBINE_SUM:
            case IGRAPH_ATTRIBUTE_COMBINE_MAX:
                IGRAPH_CHECK(igraph_i_cattributes_cb_any_is_true(oldrec, &newrec, merges));
                break;
            case IGRAPH_ATTRIBUTE_COMBINE_PROD:
            case IGRAPH_ATTRIBUTE_COMBINE_MIN:
                IGRAPH_CHECK(igraph_i_cattributes_cb_all_is_true(oldrec, &newrec, merges));
                break;
            case IGRAPH_ATTRIBUTE_COMBINE_MEAN:
            case IGRAPH_ATTRIBUTE_COMBINE_MEDIAN:
                IGRAPH_CHECK(igraph_i_cattributes_cb_majority(oldrec, &newrec, merges));
                break;
            case IGRAPH_ATTRIBUTE_COMBINE_RANDOM:
                IGRAPH_CHECK(igraph_i_cattributes_cb_random(oldrec, &newrec, merges));
                break;
            case IGRAPH_ATTRIBUTE_COMBINE_FIRST:
                IGRAPH_CHECK(igraph_i_cattributes_cb_first(oldrec, &newrec, merges));
                break;
            case IGRAPH_ATTRIBUTE_COMBINE_LAST:
                IGRAPH_CHECK(igraph_i_cattributes_cb_last(oldrec, &newrec, merges));
                break;
            case IGRAPH_ATTRIBUTE_COMBINE_CONCAT:
                IGRAPH_ERROR("Cannot calculate concatenation of Booleans.",
                             IGRAPH_EATTRCOMBINE);
                break;
            default:
                IGRAPH_ERROR("Unknown attribute combination.",
                             IGRAPH_UNIMPLEMENTED);
                break;
            }
        } else if (attr_type == IGRAPH_ATTRIBUTE_STRING) {
            switch (todo_item.type) {
            case IGRAPH_ATTRIBUTE_COMBINE_FUNCTION:
                IGRAPH_CHECK(igraph_i_cattributes_cs_func(oldrec, &newrec, merges,
                             todo_item.func.as_str));
                break;
            case IGRAPH_ATTRIBUTE_COMBINE_SUM:
                IGRAPH_ERROR("Cannot sum strings.", IGRAPH_EATTRCOMBINE);
                break;
            case IGRAPH_ATTRIBUTE_COMBINE_PROD:
                IGRAPH_ERROR("Cannot multiply strings.", IGRAPH_EATTRCOMBINE);
                break;
            case IGRAPH_ATTRIBUTE_COMBINE_MIN:
                IGRAPH_ERROR("Cannot find minimum of strings.",
                             IGRAPH_EATTRCOMBINE);
                break;
            case IGRAPH_ATTRIBUTE_COMBINE_MAX:
                IGRAPH_ERROR("Cannot find maximum of strings.",
                             IGRAPH_EATTRCOMBINE);
                break;
            case IGRAPH_ATTRIBUTE_COMBINE_MEAN:
                IGRAPH_ERROR("Cannot calculate mean of strings.",
                             IGRAPH_EATTRCOMBINE);
                break;
            case IGRAPH_ATTRIBUTE_COMBINE_MEDIAN:
                IGRAPH_ERROR("Cannot calculate median of strings.",
                             IGRAPH_EATTRCOMBINE);
                break;
            case IGRAPH_ATTRIBUTE_COMBINE_RANDOM:
                IGRAPH_CHECK(igraph_i_cattributes_sn_random(oldrec, &newrec, merges));
                break;
            case IGRAPH_ATTRIBUTE_COMBINE_FIRST:
                IGRAPH_CHECK(igraph_i_cattributes_cs_first(oldrec, &newrec, merges));
                break;
            case IGRAPH_ATTRIBUTE_COMBINE_LAST:
                IGRAPH_CHECK(igraph_i_cattributes_cs_last(oldrec, &newrec, merges));
                break;
            case IGRAPH_ATTRIBUTE_COMBINE_CONCAT:
                IGRAPH_CHECK(igraph_i_cattributes_cs_concat(oldrec, &newrec, merges));
                break;
            default:
                IGRAPH_ERROR("Unknown attribute combination.",
                             IGRAPH_UNIMPLEMENTED);
                break;
            }
        } else {
            IGRAPH_ERROR("Unknown attribute type, this should not happen.",
                         IGRAPH_UNIMPLEMENTED);
        }

        IGRAPH_CHECK(igraph_attribute_record_list_push_back(new_attrs, &newrec));
        IGRAPH_FINALLY_CLEAN(1);  /* ownership of newrec transferred */
    }

    IGRAPH_FREE(todo_items);
    IGRAPH_FINALLY_CLEAN(2);

    return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_i_cattribute_combine_vertices(
    const igraph_t *graph, igraph_t *newgraph,
    const igraph_vector_int_list_t *merges, const igraph_attribute_combination_t *comb
) {
    igraph_i_cattributes_t *attr = graph->attr;
    igraph_i_cattributes_t *toattr = newgraph->attr;
    return igraph_i_cattribute_combine_attribute_record_lists(
        &attr->val, &toattr->val, merges, comb
    );
}

static igraph_error_t igraph_i_cattribute_add_edges(
    igraph_t *graph, const igraph_vector_int_t *edges,
    const igraph_attribute_record_list_t *nattr
) {
    igraph_int_t ne = igraph_vector_int_size(edges) / 2;
    igraph_i_cattributes_t *attr = graph->attr;
    return igraph_i_cattribute_add_vertices_or_edges(&attr->eal, igraph_ecount(graph), ne, nattr);
}

static igraph_error_t igraph_i_cattribute_permute_edges(const igraph_t *graph,
                                             igraph_t *newgraph,
                                             const igraph_vector_int_t *idx) {
    igraph_i_cattributes_t *attr = graph->attr, *new_attr = newgraph->attr;
    igraph_attribute_record_list_t *eal = &attr->eal, *new_eal = &new_attr->eal;
    if (graph == newgraph) {
        return igraph_i_cattribute_permute_attribute_record_list_in_place(eal, idx);
    } else {
        return igraph_i_cattribute_permute_attribute_record_list(eal, new_eal, idx);
    }
}

static igraph_error_t igraph_i_cattribute_combine_edges(
    const igraph_t *graph, igraph_t *newgraph,
    const igraph_vector_int_list_t *merges, const igraph_attribute_combination_t *comb
) {
    igraph_i_cattributes_t *attr = graph->attr;
    igraph_i_cattributes_t *toattr = newgraph->attr;
    return igraph_i_cattribute_combine_attribute_record_lists(
        &attr->eal, &toattr->eal, merges, comb
    );
}

static igraph_error_t igraph_i_cattribute_get_info(const igraph_t *graph,
                                        igraph_strvector_t *gnames,
                                        igraph_vector_int_t *gtypes,
                                        igraph_strvector_t *vnames,
                                        igraph_vector_int_t *vtypes,
                                        igraph_strvector_t *enames,
                                        igraph_vector_int_t *etypes) {

    igraph_strvector_t *names[3] = { gnames, vnames, enames };
    igraph_vector_int_t *types[3] = { gtypes, vtypes, etypes };
    igraph_i_cattributes_t *at = graph->attr;
    igraph_attribute_record_list_t *attr[3] = { &at->gal, &at->val, &at->eal };

    for (igraph_int_t i = 0; i < 3; i++) {
        igraph_strvector_t *n = names[i];
        igraph_vector_int_t *t = types[i];
        const igraph_attribute_record_list_t *al = attr[i];
        igraph_int_t len = igraph_attribute_record_list_size(al);

        if (n) {
            IGRAPH_CHECK(igraph_strvector_resize(n, len));
        }
        if (t) {
            IGRAPH_CHECK(igraph_vector_int_resize(t, len));
        }

        for (igraph_int_t j = 0; j < len; j++) {
            const igraph_attribute_record_t *rec = igraph_attribute_record_list_get_ptr(al, j);
            const char *name = rec->name;
            igraph_attribute_type_t type = rec->type;
            if (n) {
                IGRAPH_CHECK(igraph_strvector_set(n, j, name));
            }
            if (t) {
                VECTOR(*t)[j] = type;
            }
        }
    }

    return IGRAPH_SUCCESS;
}

static igraph_bool_t igraph_i_cattribute_has_attr(const igraph_t *graph,
                                                  igraph_attribute_elemtype_t type,
                                                  const char *name) {
    const igraph_i_cattributes_t *at = graph->attr;
    switch (type) {
    case IGRAPH_ATTRIBUTE_GRAPH:
        return igraph_i_cattribute_find_index(&at->gal, name) >= 0;
    case IGRAPH_ATTRIBUTE_VERTEX:
        return igraph_i_cattribute_find_index(&at->val, name) >= 0;
    case IGRAPH_ATTRIBUTE_EDGE:
        return igraph_i_cattribute_find_index(&at->eal, name) >= 0;
    default:
        IGRAPH_ERROR("Unknown attribute element type.", IGRAPH_EINVAL);
        break;
    }

    return false;
}

static igraph_error_t igraph_i_cattribute_get_type(const igraph_t *graph,
                                       igraph_attribute_type_t *type,
                                       igraph_attribute_elemtype_t elemtype,
                                       const char *name) {
    igraph_attribute_record_t *rec;
    igraph_i_cattributes_t *at = graph->attr;
    igraph_attribute_record_list_t *al;

    switch (elemtype) {
    case IGRAPH_ATTRIBUTE_GRAPH:
        al = &at->gal;
        break;
    case IGRAPH_ATTRIBUTE_VERTEX:
        al = &at->val;
        break;
    case IGRAPH_ATTRIBUTE_EDGE:
        al = &at->eal;
        break;
    default:
        IGRAPH_ERROR("Unknown attribute element type.", IGRAPH_EINVAL);
        break;
    }

    IGRAPH_CHECK(igraph_i_cattribute_find_or_return(al, name, IGRAPH_ATTRIBUTE_UNSPECIFIED, &rec));
    *type = rec->type;

    return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_i_cattribute_get_numeric_graph_attr(
    const igraph_t *graph, const char *name, igraph_vector_t *value
) {
    igraph_i_cattributes_t *attr = graph->attr;
    igraph_attribute_record_list_t *gal = &attr->gal;
    igraph_attribute_record_t *rec;

    IGRAPH_CHECK(igraph_i_cattribute_find_or_return(gal, name, IGRAPH_ATTRIBUTE_NUMERIC, &rec));
    IGRAPH_CHECK(igraph_vector_push_back(value, VECTOR(*rec->value.as_vector)[0]));

    return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_i_cattribute_get_bool_graph_attr(
    const igraph_t *graph, const char *name, igraph_vector_bool_t *value
) {
    igraph_i_cattributes_t *attr = graph->attr;
    igraph_attribute_record_list_t *gal = &attr->gal;
    igraph_attribute_record_t *rec;

    IGRAPH_CHECK(igraph_i_cattribute_find_or_return(gal, name, IGRAPH_ATTRIBUTE_BOOLEAN, &rec));
    IGRAPH_CHECK(igraph_vector_bool_push_back(value, VECTOR(*rec->value.as_vector_bool)[0]));

    return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_i_cattribute_get_string_graph_attr(
    const igraph_t *graph, const char *name, igraph_strvector_t *value
) {
    igraph_i_cattributes_t *attr = graph->attr;
    igraph_attribute_record_list_t *gal = &attr->gal;
    igraph_attribute_record_t *rec;

    IGRAPH_CHECK(igraph_i_cattribute_find_or_return(gal, name, IGRAPH_ATTRIBUTE_STRING, &rec));
    IGRAPH_CHECK(igraph_strvector_push_back(value, igraph_strvector_get(rec->value.as_strvector, 0)));

    return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_i_cattribute_get_numeric_vertex_attr(const igraph_t *graph,
                                                       const char *name,
                                                       igraph_vs_t vs,
                                                       igraph_vector_t *value) {
    igraph_i_cattributes_t *attr = graph->attr;
    igraph_attribute_record_list_t *val = &attr->val;
    igraph_attribute_record_t *rec;
    const igraph_vector_t *num;

    IGRAPH_CHECK(igraph_i_cattribute_find_or_return(val, name, IGRAPH_ATTRIBUTE_NUMERIC, &rec));

    num = rec->value.as_vector;
    if (igraph_vs_is_all(&vs)) {
        IGRAPH_CHECK(igraph_vector_append(value, num));
    } else {
        igraph_vit_t it;
        igraph_int_t i = igraph_vector_size(value);
        IGRAPH_CHECK(igraph_vit_create(graph, vs, &it));
        IGRAPH_FINALLY(igraph_vit_destroy, &it);
        IGRAPH_CHECK(igraph_vector_resize(value, i + IGRAPH_VIT_SIZE(it)));
        for (; !IGRAPH_VIT_END(it); IGRAPH_VIT_NEXT(it), i++) {
            igraph_int_t v = IGRAPH_VIT_GET(it);
            VECTOR(*value)[i] = VECTOR(*num)[v];
        }
        igraph_vit_destroy(&it);
        IGRAPH_FINALLY_CLEAN(1);
    }

    return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_i_cattribute_get_bool_vertex_attr(const igraph_t *graph,
                                                    const char *name,
                                                    igraph_vs_t vs,
                                                    igraph_vector_bool_t *value) {
    igraph_i_cattributes_t *attr = graph->attr;
    igraph_attribute_record_list_t *val = &attr->val;
    igraph_attribute_record_t *rec;
    const igraph_vector_bool_t *log;

    IGRAPH_CHECK(igraph_i_cattribute_find_or_return(val, name, IGRAPH_ATTRIBUTE_BOOLEAN, &rec));

    log = rec->value.as_vector_bool;
    if (igraph_vs_is_all(&vs)) {
        IGRAPH_CHECK(igraph_vector_bool_append(value, log));
    } else {
        igraph_vit_t it;
        igraph_int_t i = igraph_vector_bool_size(value);
        IGRAPH_CHECK(igraph_vit_create(graph, vs, &it));
        IGRAPH_FINALLY(igraph_vit_destroy, &it);
        IGRAPH_CHECK(igraph_vector_bool_resize(value, i + IGRAPH_VIT_SIZE(it)));
        for (; !IGRAPH_VIT_END(it); IGRAPH_VIT_NEXT(it), i++) {
            igraph_int_t v = IGRAPH_VIT_GET(it);
            VECTOR(*value)[i] = VECTOR(*log)[v];
        }
        igraph_vit_destroy(&it);
        IGRAPH_FINALLY_CLEAN(1);
    }

    return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_i_cattribute_get_string_vertex_attr(const igraph_t *graph,
                                                      const char *name,
                                                      igraph_vs_t vs,
                                                      igraph_strvector_t *value) {
    igraph_i_cattributes_t *attr = graph->attr;
    igraph_attribute_record_list_t *val = &attr->val;
    igraph_attribute_record_t *rec;
    const igraph_strvector_t *str;

    IGRAPH_CHECK(igraph_i_cattribute_find_or_return(val, name, IGRAPH_ATTRIBUTE_STRING, &rec));

    str = rec->value.as_strvector;
    if (igraph_vs_is_all(&vs)) {
        igraph_strvector_clear(value);
        IGRAPH_CHECK(igraph_strvector_append(value, str));
    } else {
        igraph_vit_t it;
        igraph_int_t i = igraph_strvector_size(value);
        IGRAPH_CHECK(igraph_vit_create(graph, vs, &it));
        IGRAPH_FINALLY(igraph_vit_destroy, &it);
        IGRAPH_CHECK(igraph_strvector_resize(value, i + IGRAPH_VIT_SIZE(it)));
        for (; !IGRAPH_VIT_END(it); IGRAPH_VIT_NEXT(it), i++) {
            igraph_int_t v = IGRAPH_VIT_GET(it);
            IGRAPH_CHECK(igraph_strvector_set(value, i, igraph_strvector_get(str, v)));
        }
        igraph_vit_destroy(&it);
        IGRAPH_FINALLY_CLEAN(1);
    }

    return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_i_cattribute_get_numeric_edge_attr(
    const igraph_t *graph, const char *name, igraph_es_t es, igraph_vector_t *value
) {
    igraph_i_cattributes_t *attr = graph->attr;
    igraph_attribute_record_list_t *eal = &attr->eal;
    igraph_attribute_record_t *rec;
    const igraph_vector_t *num;

    IGRAPH_CHECK(igraph_i_cattribute_find_or_return(eal, name, IGRAPH_ATTRIBUTE_NUMERIC, &rec));

    num = rec->value.as_vector;
    if (igraph_es_is_all(&es)) {
        IGRAPH_CHECK(igraph_vector_append(value, num));
    } else {
        igraph_eit_t it;
        igraph_int_t i = igraph_vector_size(value);
        IGRAPH_CHECK(igraph_eit_create(graph, es, &it));
        IGRAPH_FINALLY(igraph_eit_destroy, &it);
        IGRAPH_CHECK(igraph_vector_resize(value, i + IGRAPH_EIT_SIZE(it)));
        for (; !IGRAPH_EIT_END(it); IGRAPH_EIT_NEXT(it), i++) {
            igraph_int_t e = IGRAPH_EIT_GET(it);
            VECTOR(*value)[i] = VECTOR(*num)[e];
        }
        igraph_eit_destroy(&it);
        IGRAPH_FINALLY_CLEAN(1);
    }

    return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_i_cattribute_get_string_edge_attr(
    const igraph_t *graph, const char *name, igraph_es_t es,
    igraph_strvector_t *value
) {
    igraph_i_cattributes_t *attr = graph->attr;
    igraph_attribute_record_list_t *eal = &attr->eal;
    igraph_attribute_record_t *rec;
    const igraph_strvector_t *str;

    IGRAPH_CHECK(igraph_i_cattribute_find_or_return(eal, name, IGRAPH_ATTRIBUTE_STRING, &rec));

    str = rec->value.as_strvector;
    if (igraph_es_is_all(&es)) {
        IGRAPH_CHECK(igraph_strvector_append(value, str));
    } else {
        igraph_eit_t it;
        igraph_int_t i = igraph_strvector_size(value);
        IGRAPH_CHECK(igraph_eit_create(graph, es, &it));
        IGRAPH_FINALLY(igraph_eit_destroy, &it);
        IGRAPH_CHECK(igraph_strvector_resize(value, i + IGRAPH_EIT_SIZE(it)));
        for (; !IGRAPH_EIT_END(it); IGRAPH_EIT_NEXT(it), i++) {
            igraph_int_t e = IGRAPH_EIT_GET(it);
            IGRAPH_CHECK(igraph_strvector_set(value, i, igraph_strvector_get(str, e)));
        }
        igraph_eit_destroy(&it);
        IGRAPH_FINALLY_CLEAN(1);
    }

    return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_i_cattribute_get_bool_edge_attr(
    const igraph_t *graph, const char *name, igraph_es_t es,
    igraph_vector_bool_t *value
) {
    igraph_i_cattributes_t *attr = graph->attr;
    igraph_attribute_record_list_t *eal = &attr->eal;
    igraph_attribute_record_t *rec;
    const igraph_vector_bool_t *log;

    IGRAPH_CHECK(igraph_i_cattribute_find_or_return(eal, name, IGRAPH_ATTRIBUTE_BOOLEAN, &rec));

    log = rec->value.as_vector_bool;
    if (igraph_es_is_all(&es)) {
        IGRAPH_CHECK(igraph_vector_bool_append(value, log));
    } else {
        igraph_eit_t it;
        igraph_int_t i = igraph_vector_bool_size(value);
        IGRAPH_CHECK(igraph_eit_create(graph, es, &it));
        IGRAPH_FINALLY(igraph_eit_destroy, &it);
        IGRAPH_CHECK(igraph_vector_bool_resize(value, i + IGRAPH_EIT_SIZE(it)));
        for (; !IGRAPH_EIT_END(it); IGRAPH_EIT_NEXT(it), i++) {
            igraph_int_t e = IGRAPH_EIT_GET(it);
            VECTOR(*value)[i] = VECTOR(*log)[e];
        }
        igraph_eit_destroy(&it);
        IGRAPH_FINALLY_CLEAN(1);
    }

    return IGRAPH_SUCCESS;
}

/* -------------------------------------- */

const igraph_attribute_table_t igraph_cattribute_table = {
    &igraph_i_cattribute_init,
    &igraph_i_cattribute_destroy,
    &igraph_i_cattribute_copy,
    &igraph_i_cattribute_add_vertices,
    &igraph_i_cattribute_permute_vertices,
    &igraph_i_cattribute_combine_vertices,
    &igraph_i_cattribute_add_edges,
    &igraph_i_cattribute_permute_edges,
    &igraph_i_cattribute_combine_edges,
    &igraph_i_cattribute_get_info,
    &igraph_i_cattribute_has_attr,
    &igraph_i_cattribute_get_type,
    &igraph_i_cattribute_get_numeric_graph_attr,
    &igraph_i_cattribute_get_string_graph_attr,
    &igraph_i_cattribute_get_bool_graph_attr,
    &igraph_i_cattribute_get_numeric_vertex_attr,
    &igraph_i_cattribute_get_string_vertex_attr,
    &igraph_i_cattribute_get_bool_vertex_attr,
    &igraph_i_cattribute_get_numeric_edge_attr,
    &igraph_i_cattribute_get_string_edge_attr,
    &igraph_i_cattribute_get_bool_edge_attr
};

/* -------------------------------------- */

/**
 * \section cattributes
 * <para>There is an experimental attribute handler that can be used
 * from C code. In this section we show how this works. This attribute
 * handler is by default not attached (the default is no attribute
 * handler), so we first need to attach it:
 * <programlisting>
 * igraph_set_attribute_table(&amp;igraph_cattribute_table);
 * </programlisting>
 * </para>
 * <para>Now the attribute functions are available. Please note that
 * the attribute handler must be attached before you call any other
 * igraph functions, otherwise you might end up with graphs without
 * attributes and an active attribute handler, which might cause
 * unexpected program behaviour. The rule is that you attach the
 * attribute handler in the beginning of your
 * <function>main()</function> and never touch it again. Detaching
 * the attribute handler might lead to memory leaks.</para>
 *
 * <para>It is not currently possible to have attribute handlers on a
 * per-graph basis. All graphs in an application must be managed with
 * the same attribute handler. This also applies to the default case
 * when there is no attribute handler at all.</para>
 *
 * <para>The C attribute handler supports attaching real numbers, boolean
 * values and character strings as attributes. No vector values are allowed.
 * For example, vertices have a <code>name</code> attribute holding a single
 * string value for each vertex, but it is not possible to have a <code>coords</code>
 * attribute which is a vector of numbers per vertex.</para>
 *
 * <para>The functions documented in this section are specific to the C
 * attribute handler. Code using these functions will not function when
 * a different attribute handler is attached.</para>
 *
 * \example examples/simple/cattributes.c
 * \example examples/simple/cattributes2.c
 * \example examples/simple/cattributes3.c
 * \example examples/simple/cattributes4.c
 */

/**
 * \function igraph_cattribute_GAN
 * \brief Query a numeric graph attribute.
 *
 * Returns the value of the given numeric graph attribute.
 * If the attribute does not exist, a warning is issued
 * and NaN is returned.
 *
 * \param graph The input graph.
 * \param name The name of the attribute to query.
 * \return The value of the attribute.
 *
 * \sa \ref GAN for a simpler interface.
 *
 * Time complexity: O(Ag), the number of graph attributes.
 */
igraph_real_t igraph_cattribute_GAN(const igraph_t *graph, const char *name) {
    igraph_i_cattributes_t *attr = graph->attr;
    const igraph_attribute_record_t *rec;
    const igraph_vector_t *num;

    rec = igraph_i_cattribute_find(&attr->gal, name, IGRAPH_ATTRIBUTE_NUMERIC);
    if (!rec) {
        IGRAPH_WARNINGF("Graph attribute '%s' does not exist, returning default numeric attribute value.", name);
        return IGRAPH_NAN;
    }

    num = rec->value.as_vector;
    return VECTOR(*num)[0];
}

/**
 * \function igraph_cattribute_GAB
 * \brief Query a boolean graph attribute.
 *
 * Returns the value of the given boolean graph attribute.
 * If the attribute does not exist, a warning is issued
 * and false is returned.
 *
 * \param graph The input graph.
 * \param name The name of the attribute to query.
 * \return The value of the attribute.
 *
 * \sa \ref GAB for a simpler interface.
 *
 * Time complexity: O(Ag), the number of graph attributes.
 */
igraph_bool_t igraph_cattribute_GAB(const igraph_t *graph, const char *name) {
    igraph_i_cattributes_t *attr = graph->attr;
    const igraph_attribute_record_t *rec;
    const igraph_vector_bool_t *log;

    rec = igraph_i_cattribute_find(&attr->gal, name, IGRAPH_ATTRIBUTE_BOOLEAN);
    if (!rec) {
        IGRAPH_WARNINGF("Graph attribute '%s' does not exist, returning default boolean attribute value.", name);
        return false;
    }

    log = rec->value.as_vector_bool;
    return VECTOR(*log)[0];
}

/**
 * \function igraph_cattribute_GAS
 * \brief Query a string graph attribute.
 *
 * Returns a <type>const</type> pointer to the string graph attribute
 * specified in \p name. The value must not be modified.
 * If the attribute does not exist, a warning is issued and
 * an empty string is returned.
 *
 * \param graph The input graph.
 * \param name The name of the attribute to query.
 * \return The value of the attribute.
 *
 * \sa \ref GAS for a simpler interface.
 *
 * Time complexity: O(Ag), the number of graph attributes.
 */
const char *igraph_cattribute_GAS(const igraph_t *graph, const char *name) {
    igraph_i_cattributes_t *attr = graph->attr;
    const igraph_attribute_record_t *rec;
    const igraph_strvector_t *str;

    rec = igraph_i_cattribute_find(&attr->gal, name, IGRAPH_ATTRIBUTE_STRING);
    if (!rec) {
        IGRAPH_WARNINGF("Graph attribute '%s' does not exist, returning default string attribute value.", name);
        return "";
    }

    str = rec->value.as_strvector;
    return igraph_strvector_get(str, 0);
}

/**
 * \function igraph_cattribute_VAN
 * \brief Query a numeric vertex attribute.
 *
 * If the attribute does not exist, a warning is issued and
 * NaN is returned. See \ref igraph_cattribute_VANV() for
 * an error-checked version.
 *
 * \param graph The input graph.
 * \param name The name of the attribute.
 * \param vid The id of the queried vertex.
 * \return The value of the attribute.
 *
 * \sa \ref VAN macro for a simpler interface.
 *
 * Time complexity: O(Av), the number of vertex attributes.
 */
igraph_real_t igraph_cattribute_VAN(const igraph_t *graph, const char *name,
                                    igraph_int_t vid) {
    igraph_i_cattributes_t *attr = graph->attr;
    const igraph_attribute_record_t *rec;
    const igraph_vector_t *num;

    rec = igraph_i_cattribute_find(&attr->val, name, IGRAPH_ATTRIBUTE_NUMERIC);
    if (!rec) {
        IGRAPH_WARNINGF("Vertex attribute '%s' does not exist, returning default numeric attribute value.", name);
        return IGRAPH_NAN;
    }

    num = rec->value.as_vector;
    return VECTOR(*num)[vid];
}

/**
 * \function igraph_cattribute_VAB
 * \brief Query a boolean vertex attribute.
 *
 * If the vertex attribute does not exist, a warning is issued
 * and false is returned. See \ref igraph_cattribute_VABV() for
 * an error-checked version.
 *
 * \param graph The input graph.
 * \param name The name of the attribute.
 * \param vid The id of the queried vertex.
 * \return The value of the attribute.
 *
 * \sa \ref VAB macro for a simpler interface.
 *
 * Time complexity: O(Av), the number of vertex attributes.
 */
igraph_bool_t igraph_cattribute_VAB(const igraph_t *graph, const char *name,
                                    igraph_int_t vid) {
    igraph_i_cattributes_t *attr = graph->attr;
    const igraph_attribute_record_t *rec;
    const igraph_vector_bool_t *log;

    rec = igraph_i_cattribute_find(&attr->val, name, IGRAPH_ATTRIBUTE_BOOLEAN);
    if (!rec) {
        IGRAPH_WARNINGF("Vertex attribute '%s' does not exist, returning default boolean attribute value.", name);
        return false;
    }

    log = rec->value.as_vector_bool;
    return VECTOR(*log)[vid];
}

/**
 * \function igraph_cattribute_VAS
 * \brief Query a string vertex attribute.
 *
 * Returns a <type>const</type> pointer to the string vertex attribute
 * specified in \p name. The value must not be modified.
 * If the vertex attribute does not exist, a warning is issued and
 * an empty string is returned. See \ref igraph_cattribute_VASV()
 * for an error-checked version.
 *
 * \param graph The input graph.
 * \param name The name of the attribute.
 * \param vid The id of the queried vertex.
 * \return The value of the attribute.
 *
 * \sa The macro \ref VAS for a simpler interface.
 *
 * Time complexity: O(Av), the number of vertex attributes.
 */
const char *igraph_cattribute_VAS(const igraph_t *graph, const char *name,
                                  igraph_int_t vid) {
    igraph_i_cattributes_t *attr = graph->attr;
    const igraph_attribute_record_t *rec;
    const igraph_strvector_t *str;

    rec = igraph_i_cattribute_find(&attr->val, name, IGRAPH_ATTRIBUTE_STRING);
    if (!rec) {
        IGRAPH_WARNINGF("Vertex attribute '%s' does not exist, returning default string attribute value.", name);
        return "";
    }

    str = rec->value.as_strvector;
    return igraph_strvector_get(str, vid);
}

/**
 * \function igraph_cattribute_EAN
 * \brief Query a numeric edge attribute.
 *
 * If the attribute does not exist, a warning is issued and
 * NaN is returned. See \ref igraph_cattribute_EANV() for
 * an error-checked version.
 *
 * \param graph The input graph.
 * \param name The name of the attribute.
 * \param eid The id of the queried edge.
 * \return The value of the attribute.
 *
 * \sa \ref EAN for an easier interface.
 *
 * Time complexity: O(Ae), the number of edge attributes.
 */
igraph_real_t igraph_cattribute_EAN(const igraph_t *graph, const char *name,
                                    igraph_int_t eid) {
    igraph_i_cattributes_t *attr = graph->attr;
    const igraph_attribute_record_t *rec;
    const igraph_vector_t *num;

    rec = igraph_i_cattribute_find(&attr->eal, name, IGRAPH_ATTRIBUTE_NUMERIC);
    if (!rec) {
        IGRAPH_WARNINGF("Edge attribute '%s' does not exist, returning default numeric attribute value.", name);
        return IGRAPH_NAN;
    }

    num = rec->value.as_vector;
    return VECTOR(*num)[eid];
}

/**
 * \function igraph_cattribute_EAB
 * \brief Query a boolean edge attribute.
 *
 * If the edge attribute does not exist, a warning is issued and
 * false is returned. See \ref igraph_cattribute_EABV() for
 * an error-checked version.
 *
 * \param graph The input graph.
 * \param name The name of the attribute.
 * \param eid The id of the queried edge.
 * \return The value of the attribute.
 *
 * \sa \ref EAB for an easier interface.
 *
 * Time complexity: O(Ae), the number of edge attributes.
 */
igraph_bool_t igraph_cattribute_EAB(const igraph_t *graph, const char *name,
                                    igraph_int_t eid) {
    igraph_i_cattributes_t *attr = graph->attr;
    const igraph_attribute_record_t *rec;
    const igraph_vector_bool_t *log;

    rec = igraph_i_cattribute_find(&attr->eal, name, IGRAPH_ATTRIBUTE_BOOLEAN);
    if (!rec) {
        IGRAPH_WARNINGF("Edge attribute '%s' does not exist, returning default boolean attribute value.", name);
        return false;
    }

    log = rec->value.as_vector_bool;
    return VECTOR(*log)[eid];
}

/**
 * \function igraph_cattribute_EAS
 * \brief Query a string edge attribute.
 *
 * Returns a <type>const</type> pointer to the string edge attribute
 * specified in \p name. The value must not be modified.
 * If the edge attribute does not exist, a warning is issued and
 * an empty string is returned. See \ref igraph_cattribute_EASV() for
 * an error-checked version.
 *
 * \param graph The input graph.
 * \param name The name of the attribute.
 * \param eid The id of the queried edge.
 * \return The value of the attribute.
 *
 * \se \ref EAS if you want to type less.
 *
 * Time complexity: O(Ae), the number of edge attributes.
 */
const char *igraph_cattribute_EAS(const igraph_t *graph, const char *name,
                                  igraph_int_t eid) {
    igraph_i_cattributes_t *attr = graph->attr;
    const igraph_attribute_record_t *rec;
    const igraph_strvector_t *str;

    rec = igraph_i_cattribute_find(&attr->eal, name, IGRAPH_ATTRIBUTE_STRING);
    if (!rec) {
        IGRAPH_WARNINGF("Edge attribute '%s' does not exist, returning default string attribute value.", name);
        return "";
    }

    str = rec->value.as_strvector;
    return igraph_strvector_get(str, eid);
}

/**
 * \function igraph_cattribute_VANV
 * \brief Query a numeric vertex attribute for many vertices.
 *
 * \param graph The input graph.
 * \param name The name of the attribute.
 * \param vids The vertices to query.
 * \param result Pointer to an initialized vector, the result is
 *    stored here. It will be resized, if needed.
 * \return Error code.
 *
 * Time complexity: O(v), where v is the number of vertices in 'vids'.
 */

igraph_error_t igraph_cattribute_VANV(const igraph_t *graph, const char *name,
                           igraph_vs_t vids, igraph_vector_t *result) {
    igraph_vector_clear(result);
    return igraph_i_cattribute_get_numeric_vertex_attr(graph, name, vids, result);
}

/**
 * \function igraph_cattribute_VABV
 * \brief Query a boolean vertex attribute for many vertices.
 *
 * \param graph The input graph.
 * \param name The name of the attribute.
 * \param vids The vertices to query.
 * \param result Pointer to an initialized boolean vector, the result is
 *    stored here. It will be resized, if needed.
 * \return Error code.
 *
 * Time complexity: O(v), where v is the number of vertices in 'vids'.
 */

igraph_error_t igraph_cattribute_VABV(const igraph_t *graph, const char *name,
                           igraph_vs_t vids, igraph_vector_bool_t *result) {
    igraph_vector_bool_clear(result);
    return igraph_i_cattribute_get_bool_vertex_attr(graph, name, vids, result);
}

/**
 * \function igraph_cattribute_EANV
 * \brief Query a numeric edge attribute for many edges.
 *
 * \param graph The input graph.
 * \param name The name of the attribute.
 * \param eids The edges to query.
 * \param result Pointer to an initialized vector, the result is
 *    stored here. It will be resized, if needed.
 * \return Error code.
 *
 * Time complexity: O(e), where e is the number of edges in 'eids'.
 */

igraph_error_t igraph_cattribute_EANV(const igraph_t *graph, const char *name,
                           igraph_es_t eids, igraph_vector_t *result) {
    igraph_vector_clear(result);
    return igraph_i_cattribute_get_numeric_edge_attr(graph, name, eids, result);
}

/**
 * \function igraph_cattribute_EABV
 * \brief Query a boolean edge attribute for many edges.
 *
 * \param graph The input graph.
 * \param name The name of the attribute.
 * \param eids The edges to query.
 * \param result Pointer to an initialized boolean vector, the result is
 *    stored here. It will be resized, if needed.
 * \return Error code.
 *
 * Time complexity: O(e), where e is the number of edges in 'eids'.
 */

igraph_error_t igraph_cattribute_EABV(const igraph_t *graph, const char *name,
                           igraph_es_t eids, igraph_vector_bool_t *result) {
    igraph_vector_bool_clear(result);
    return igraph_i_cattribute_get_bool_edge_attr(graph, name, eids, result);
}

/**
 * \function igraph_cattribute_VASV
 * \brief Query a string vertex attribute for many vertices.
 *
 * \param graph The input graph.
 * \param name The name of the attribute.
 * \param vids The vertices to query.
 * \param result Pointer to an initialized string vector, the result
 *     is stored here. It will be resized, if needed.
 * \return Error code.
 *
 * Time complexity: O(v), where v is the number of vertices in 'vids'.
 * (We assume that the string attributes have a bounded length.)
 */

igraph_error_t igraph_cattribute_VASV(const igraph_t *graph, const char *name,
                           igraph_vs_t vids, igraph_strvector_t *result) {
    igraph_strvector_clear(result);
    return igraph_i_cattribute_get_string_vertex_attr(graph, name, vids, result);
}

/**
 * \function igraph_cattribute_EASV
 * \brief Query a string edge attribute for many edges.
 *
 * \param graph The input graph.
 * \param name The name of the attribute.
 * \param eids The edges to query.
 * \param result Pointer to an initialized string vector, the result
 *     is stored here. It will be resized, if needed.
 * \return Error code.
 *
 * Time complexity: O(e), where e is the number of edges in
 * 'eids'. (We assume that the string attributes have a bounded length.)
 */

igraph_error_t igraph_cattribute_EASV(const igraph_t *graph, const char *name,
                           igraph_es_t eids, igraph_strvector_t *result) {
    igraph_strvector_clear(result);
    return igraph_i_cattribute_get_string_edge_attr(graph, name, eids, result);
}

/**
 * \function igraph_cattribute_list
 * \brief List all attributes.
 *
 * See \ref igraph_attribute_type_t for the various attribute types.
 * \param graph The input graph.
 * \param gnames String vector, the names of the graph attributes.
 * \param gtypes Numeric vector, the types of the graph attributes.
 * \param vnames String vector, the names of the vertex attributes.
 * \param vtypes Numeric vector, the types of the vertex attributes.
 * \param enames String vector, the names of the edge attributes.
 * \param etypes Numeric vector, the types of the edge attributes.
 * \return Error code.
 *
 * Naturally, the string vector with the attribute names and the
 * numeric vector with the attribute types are in the right order,
 * i.e. the first name corresponds to the first type, etc.
 *
 * Time complexity: O(Ag+Av+Ae), the number of all attributes.
 */
igraph_error_t igraph_cattribute_list(const igraph_t *graph,
                           igraph_strvector_t *gnames, igraph_vector_int_t *gtypes,
                           igraph_strvector_t *vnames, igraph_vector_int_t *vtypes,
                           igraph_strvector_t *enames, igraph_vector_int_t *etypes) {
    return igraph_i_cattribute_get_info(graph, gnames, gtypes, vnames, vtypes,
                                        enames, etypes);
}

/**
 * \function igraph_cattribute_has_attr
 * \brief Checks whether a (graph, vertex or edge) attribute exists.
 *
 * \param graph The graph.
 * \param type The type of the attribute, \c IGRAPH_ATTRIBUTE_GRAPH,
 *        \c IGRAPH_ATTRIBUTE_VERTEX or \c IGRAPH_ATTRIBUTE_EDGE.
 * \param name Character constant, the name of the attribute.
 * \return Boolean value, \c true if the attribute exists, \c false otherwise.
 *
 * Time complexity: O(A), the number of (graph, vertex or edge)
 * attributes, assuming attribute names are not too long.
 */
igraph_bool_t igraph_cattribute_has_attr(const igraph_t *graph,
        igraph_attribute_elemtype_t type,
        const char *name) {
    return igraph_i_cattribute_has_attr(graph, type, name);
}

/**
 * \function igraph_cattribute_GAN_set
 * \brief Set a numeric graph attribute.
 *
 * \param graph The graph.
 * \param name Name of the graph attribute. If there is no such
 *   attribute yet, then it will be added.
 * \param value The (new) value of the graph attribute.
 * \return Error code.
 *
 * \se \ref SETGAN if you want to type less.
 *
 * Time complexity: O(1).
 */
igraph_error_t igraph_cattribute_GAN_set(igraph_t *graph, const char *name,
                              igraph_real_t value) {
    igraph_i_cattributes_t *attr = graph->attr;
    igraph_attribute_record_t *rec;
    igraph_vector_t *num;

    IGRAPH_CHECK(igraph_i_cattribute_find_or_create(
        &attr->gal, name, IGRAPH_ATTRIBUTE_NUMERIC, 1, &rec
    ));

    num = rec->value.as_vector;
    VECTOR(*num)[0] = value;

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_cattribute_GAB_set
 * \brief Set a boolean graph attribute.
 *
 * \param graph The graph.
 * \param name Name of the graph attribute. If there is no such
 *   attribute yet, then it will be added.
 * \param value The (new) value of the graph attribute.
 * \return Error code.
 *
 * \se \ref SETGAN if you want to type less.
 *
 * Time complexity: O(1).
 */
igraph_error_t igraph_cattribute_GAB_set(igraph_t *graph, const char *name,
                              igraph_bool_t value) {
    igraph_i_cattributes_t *attr = graph->attr;
    igraph_attribute_record_t *rec;
    igraph_vector_bool_t *log;

    IGRAPH_CHECK(igraph_i_cattribute_find_or_create(
        &attr->gal, name, IGRAPH_ATTRIBUTE_BOOLEAN, 1, &rec
    ));

    log = rec->value.as_vector_bool;
    VECTOR(*log)[0] = value;

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_cattribute_GAS_set
 * \brief Set a string graph attribute.
 *
 * \param graph The graph.
 * \param name Name of the graph attribute. If there is no such
 *   attribute yet, then it will be added.
 * \param value The (new) value of the graph attribute. It will be
 *   copied.
 * \return Error code.
 *
 * \se \ref SETGAS if you want to type less.
 *
 * Time complexity: O(1).
 */
igraph_error_t igraph_cattribute_GAS_set(igraph_t *graph, const char *name,
                              const char *value) {
    igraph_i_cattributes_t *attr = graph->attr;
    igraph_attribute_record_t *rec;
    igraph_strvector_t *str;

    IGRAPH_CHECK(igraph_i_cattribute_find_or_create(
        &attr->gal, name, IGRAPH_ATTRIBUTE_STRING, 1, &rec
    ));

    str = rec->value.as_strvector;
    IGRAPH_CHECK(igraph_strvector_set(str, 0, value));

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_cattribute_VAN_set
 * \brief Set a numeric vertex attribute.
 *
 * The attribute will be added if not present already. If present it
 * will be overwritten. The same \p value is set for all vertices
 * included in \p vid.
 * \param graph The graph.
 * \param name Name of the attribute.
 * \param vid Vertices for which to set the attribute.
 * \param value The (new) value of the attribute.
 * \return Error code.
 *
 * \sa \ref SETVAN for a simpler way.
 *
 * Time complexity: O(n), the number of vertices if the attribute is
 * new, O(|vid|) otherwise.
 */
igraph_error_t igraph_cattribute_VAN_set(igraph_t *graph, const char *name,
                              igraph_int_t vid, igraph_real_t value) {
    igraph_i_cattributes_t *attr = graph->attr;
    igraph_attribute_record_t *rec;

    IGRAPH_CHECK(igraph_i_cattribute_find_or_create(
        &attr->val, name, IGRAPH_ATTRIBUTE_NUMERIC, igraph_vcount(graph), &rec
    ));
    VECTOR(*rec->value.as_vector)[vid] = value;

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_cattribute_VAB_set
 * \brief Set a boolean vertex attribute.
 *
 * The attribute will be added if not present already. If present it
 * will be overwritten. The same \p value is set for all vertices
 * included in \p vid.
 * \param graph The graph.
 * \param name Name of the attribute.
 * \param vid Vertices for which to set the attribute.
 * \param value The (new) value of the attribute.
 * \return Error code.
 *
 * \sa \ref SETVAB for a simpler way.
 *
 * Time complexity: O(n), the number of vertices if the attribute is
 * new, O(|vid|) otherwise.
 */
igraph_error_t igraph_cattribute_VAB_set(igraph_t *graph, const char *name,
                              igraph_int_t vid, igraph_bool_t value) {
    igraph_i_cattributes_t *attr = graph->attr;
    igraph_attribute_record_t *rec;

    IGRAPH_CHECK(igraph_i_cattribute_find_or_create(
        &attr->val, name, IGRAPH_ATTRIBUTE_BOOLEAN, igraph_vcount(graph), &rec
    ));
    VECTOR(*rec->value.as_vector_bool)[vid] = value;

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_cattribute_VAS_set
 * \brief Set a string vertex attribute.
 *
 * The attribute will be added if not present already. If present it
 * will be overwritten. The same \p value is set for all vertices
 * included in \p vid.
 * \param graph The graph.
 * \param name Name of the attribute.
 * \param vid Vertices for which to set the attribute.
 * \param value The (new) value of the attribute.
 * \return Error code.
 *
 * \sa \ref SETVAS for a simpler way.
 *
 * Time complexity: O(n*l), n is the number of vertices, l is the
 * length of the string to set. If the attribute if not new then only
 * O(|vid|*l).
 */
igraph_error_t igraph_cattribute_VAS_set(igraph_t *graph, const char *name,
                              igraph_int_t vid, const char *value) {
    igraph_i_cattributes_t *attr = graph->attr;
    igraph_attribute_record_t *rec;

    IGRAPH_CHECK(igraph_i_cattribute_find_or_create(
        &attr->val, name, IGRAPH_ATTRIBUTE_STRING, igraph_vcount(graph), &rec
    ));
    IGRAPH_CHECK(igraph_strvector_set(rec->value.as_strvector, vid, value));

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_cattribute_EAN_set
 * \brief Set a numeric edge attribute.
 *
 * The attribute will be added if not present already. If present it
 * will be overwritten. The same \p value is set for all edges
 * included in \p vid.
 * \param graph The graph.
 * \param name Name of the attribute.
 * \param eid Edges for which to set the attribute.
 * \param value The (new) value of the attribute.
 * \return Error code.
 *
 * \sa \ref SETEAN for a simpler way.
 *
 * Time complexity: O(e), the number of edges if the attribute is
 * new, O(|eid|) otherwise.
 */
igraph_error_t igraph_cattribute_EAN_set(igraph_t *graph, const char *name,
                              igraph_int_t eid, igraph_real_t value) {
    igraph_i_cattributes_t *attr = graph->attr;
    igraph_attribute_record_t *rec;

    IGRAPH_CHECK(igraph_i_cattribute_find_or_create(
        &attr->eal, name, IGRAPH_ATTRIBUTE_NUMERIC, igraph_ecount(graph), &rec
    ));
    VECTOR(*rec->value.as_vector)[eid] = value;

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_cattribute_EAB_set
 * \brief Set a boolean edge attribute.
 *
 * The attribute will be added if not present already. If present it
 * will be overwritten. The same \p value is set for all edges
 * included in \p vid.
 * \param graph The graph.
 * \param name Name of the attribute.
 * \param eid Edges for which to set the attribute.
 * \param value The (new) value of the attribute.
 * \return Error code.
 *
 * \sa \ref SETEAB for a simpler way.
 *
 * Time complexity: O(e), the number of edges if the attribute is
 * new, O(|eid|) otherwise.
 */
igraph_error_t igraph_cattribute_EAB_set(igraph_t *graph, const char *name,
                              igraph_int_t eid, igraph_bool_t value) {
    igraph_i_cattributes_t *attr = graph->attr;
    igraph_attribute_record_t *rec;

    IGRAPH_CHECK(igraph_i_cattribute_find_or_create(
        &attr->eal, name, IGRAPH_ATTRIBUTE_BOOLEAN, igraph_ecount(graph), &rec
    ));
    VECTOR(*rec->value.as_vector_bool)[eid] = value;

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_cattribute_EAS_set
 * \brief Set a string edge attribute.
 *
 * The attribute will be added if not present already. If present it
 * will be overwritten. The same \p value is set for all edges
 * included in \p vid.
 * \param graph The graph.
 * \param name Name of the attribute.
 * \param eid Edges for which to set the attribute.
 * \param value The (new) value of the attribute.
 * \return Error code.
 *
 * \sa \ref SETEAS for a simpler way.
 *
 * Time complexity: O(e*l), n is the number of edges, l is the
 * length of the string to set. If the attribute if not new then only
 * O(|eid|*l).
 */
igraph_error_t igraph_cattribute_EAS_set(igraph_t *graph, const char *name,
                              igraph_int_t eid, const char *value) {
    igraph_i_cattributes_t *attr = graph->attr;
    igraph_attribute_record_t *rec;

    IGRAPH_CHECK(igraph_i_cattribute_find_or_create(
        &attr->eal, name, IGRAPH_ATTRIBUTE_STRING, igraph_ecount(graph), &rec
    ));
    IGRAPH_CHECK(igraph_strvector_set(rec->value.as_strvector, eid, value));

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_cattribute_VAN_setv
 * \brief Set a numeric vertex attribute for all vertices.
 *
 * The attribute will be added if not present yet.
 * \param graph The graph.
 * \param name Name of the attribute.
 * \param v The new attribute values. The length of this vector must
 *   match the number of vertices.
 * \return Error code.
 *
 * \sa \ref SETVANV for a simpler way.
 *
 * Time complexity: O(n), the number of vertices.
 */

igraph_error_t igraph_cattribute_VAN_setv(igraph_t *graph, const char *name,
                               const igraph_vector_t *v) {
    igraph_i_cattributes_t *attr = graph->attr;
    igraph_attribute_record_t *rec;
    igraph_int_t nv = igraph_vcount(graph);

    /* Check length first */
    if (igraph_vector_size(v) != nv) {
        IGRAPH_ERROR("Invalid vertex attribute vector length.", IGRAPH_EINVAL);
    }

    IGRAPH_CHECK(igraph_i_cattribute_find_or_create(&attr->val, name, IGRAPH_ATTRIBUTE_NUMERIC, nv, &rec));
    IGRAPH_CHECK(igraph_vector_update(rec->value.as_vector, v));

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_cattribute_VAB_setv
 * \brief Set a boolean vertex attribute for all vertices.
 *
 * The attribute will be added if not present yet.
 * \param graph The graph.
 * \param name Name of the attribute.
 * \param v The new attribute values. The length of this boolean vector must
 *   match the number of vertices.
 * \return Error code.
 *
 * \sa \ref SETVANV for a simpler way.
 *
 * Time complexity: O(n), the number of vertices.
 */

igraph_error_t igraph_cattribute_VAB_setv(igraph_t *graph, const char *name,
                               const igraph_vector_bool_t *v) {
    igraph_i_cattributes_t *attr = graph->attr;
    igraph_attribute_record_t *rec;
    igraph_int_t nv = igraph_vcount(graph);

    /* Check length first */
    if (igraph_vector_bool_size(v) != nv) {
        IGRAPH_ERROR("Invalid vertex attribute vector length.", IGRAPH_EINVAL);
    }

    IGRAPH_CHECK(igraph_i_cattribute_find_or_create(&attr->val, name, IGRAPH_ATTRIBUTE_BOOLEAN, nv, &rec));
    IGRAPH_CHECK(igraph_vector_bool_update(rec->value.as_vector_bool, v));

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_cattribute_VAS_setv
 * \brief Set a string vertex attribute for all vertices.
 *
 * The attribute will be added if not present yet.
 * \param graph The graph.
 * \param name Name of the attribute.
 * \param sv String vector, the new attribute values. The length of this vector must
 *   match the number of vertices.
 * \return Error code.
 *
 * \sa \ref SETVASV for a simpler way.
 *
 * Time complexity: O(n+l), n is the number of vertices, l is the
 * total length of the strings.
 */
igraph_error_t igraph_cattribute_VAS_setv(igraph_t *graph, const char *name,
                               const igraph_strvector_t *sv) {
    igraph_i_cattributes_t *attr = graph->attr;
    igraph_attribute_record_t *rec;
    igraph_int_t nv = igraph_vcount(graph);

    /* Check length first */
    if (igraph_strvector_size(sv) != nv) {
        IGRAPH_ERROR("Invalid vertex attribute vector length.", IGRAPH_EINVAL);
    }

    IGRAPH_CHECK(igraph_i_cattribute_find_or_create(&attr->val, name, IGRAPH_ATTRIBUTE_STRING, nv, &rec));
    IGRAPH_CHECK(igraph_strvector_update(rec->value.as_strvector, sv));

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_cattribute_EAN_setv
 * \brief Set a numeric edge attribute for all edges.
 *
 * The attribute will be added if not present yet.
 * \param graph The graph.
 * \param name Name of the attribute.
 * \param v The new attribute values. The length of this vector must
 *   match the number of edges.
 * \return Error code.
 *
 * \sa \ref SETEANV for a simpler way.
 *
 * Time complexity: O(e), the number of edges.
 */
igraph_error_t igraph_cattribute_EAN_setv(igraph_t *graph, const char *name,
                               const igraph_vector_t *v) {
    igraph_i_cattributes_t *attr = graph->attr;
    igraph_attribute_record_t *rec;
    igraph_int_t ne = igraph_ecount(graph);

    /* Check length first */
    if (igraph_vector_size(v) != ne) {
        IGRAPH_ERROR("Invalid edge attribute vector length.", IGRAPH_EINVAL);
    }

    IGRAPH_CHECK(igraph_i_cattribute_find_or_create(&attr->eal, name, IGRAPH_ATTRIBUTE_NUMERIC, ne, &rec));
    IGRAPH_CHECK(igraph_vector_update(rec->value.as_vector, v));

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_cattribute_EAB_setv
 * \brief Set a boolean edge attribute for all edges.
 *
 * The attribute will be added if not present yet.
 * \param graph The graph.
 * \param name Name of the attribute.
 * \param v The new attribute values. The length of this vector must
 *   match the number of edges.
 * \return Error code.
 *
 * \sa \ref SETEABV for a simpler way.
 *
 * Time complexity: O(e), the number of edges.
 */
igraph_error_t igraph_cattribute_EAB_setv(igraph_t *graph, const char *name,
                               const igraph_vector_bool_t *v) {
    igraph_i_cattributes_t *attr = graph->attr;
    igraph_attribute_record_t *rec;
    igraph_int_t ne = igraph_ecount(graph);

    /* Check length first */
    if (igraph_vector_bool_size(v) != ne) {
        IGRAPH_ERROR("Invalid edge attribute vector length.", IGRAPH_EINVAL);
    }

    IGRAPH_CHECK(igraph_i_cattribute_find_or_create(&attr->eal, name, IGRAPH_ATTRIBUTE_BOOLEAN, ne, &rec));
    IGRAPH_CHECK(igraph_vector_bool_update(rec->value.as_vector_bool, v));

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_cattribute_EAS_setv
 * \brief Set a string edge attribute for all edges.
 *
 * The attribute will be added if not present yet.
 * \param graph The graph.
 * \param name Name of the attribute.
 * \param sv String vector, the new attribute values. The length of this vector must
 *   match the number of edges.
 * \return Error code.
 *
 * \sa \ref SETEASV for a simpler way.
 *
 * Time complexity: O(e+l), e is the number of edges, l is the
 * total length of the strings.
 */
igraph_error_t igraph_cattribute_EAS_setv(igraph_t *graph, const char *name,
                               const igraph_strvector_t *sv) {
    igraph_i_cattributes_t *attr = graph->attr;
    igraph_attribute_record_t *rec;
    igraph_int_t ne = igraph_ecount(graph);

    /* Check length first */
    if (igraph_strvector_size(sv) != ne) {
        IGRAPH_ERROR("Invalid edge attribute vector length.", IGRAPH_EINVAL);
    }

    IGRAPH_CHECK(igraph_i_cattribute_find_or_create(&attr->eal, name, IGRAPH_ATTRIBUTE_STRING, ne, &rec));
    IGRAPH_CHECK(igraph_strvector_update(rec->value.as_strvector, sv));

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_cattribute_remove_g
 * \brief Remove a graph attribute.
 *
 * \param graph The graph object.
 * \param name Name of the graph attribute to remove.
 *
 * \sa \ref DELGA for a simpler way.
 *
 */
void igraph_cattribute_remove_g(igraph_t *graph, const char *name) {

    igraph_i_cattributes_t *attr = graph->attr;
    igraph_attribute_record_list_t *gal = &attr->gal;
    igraph_int_t j = igraph_i_cattribute_find_index(gal, name);

    if (j >= 0) {
        igraph_attribute_record_list_discard(gal, j);
    } else {
        IGRAPH_WARNING("Cannot remove non-existent graph attribute");
    }
}

/**
 * \function igraph_cattribute_remove_v
 * \brief Remove a vertex attribute.
 *
 * \param graph The graph object.
 * \param name Name of the vertex attribute to remove.
 *
 * \sa \ref DELVA for a simpler way.
 *
 */
void igraph_cattribute_remove_v(igraph_t *graph, const char *name) {

    igraph_i_cattributes_t *attr = graph->attr;
    igraph_attribute_record_list_t *val = &attr->val;
    igraph_int_t j = igraph_i_cattribute_find_index(val, name);

    if (j >= 0) {
        igraph_attribute_record_list_discard(val, j);
    } else {
        IGRAPH_WARNING("Cannot remove non-existent graph attribute");
    }
}

/**
 * \function igraph_cattribute_remove_e
 * \brief Remove an edge attribute.
 *
 * \param graph The graph object.
 * \param name Name of the edge attribute to remove.
 *
 * \sa \ref DELEA for a simpler way.
 */
void igraph_cattribute_remove_e(igraph_t *graph, const char *name) {

    igraph_i_cattributes_t *attr = graph->attr;
    igraph_attribute_record_list_t *eal = &attr->eal;
    igraph_int_t j = igraph_i_cattribute_find_index(eal, name);

    if (j >= 0) {
        igraph_attribute_record_list_discard(eal, j);
    } else {
        IGRAPH_WARNING("Cannot remove non-existent graph attribute");
    }
}

/**
 * \function igraph_cattribute_remove_all
 * \brief Remove all graph/vertex/edge attributes.
 *
 * \param graph The graph object.
 * \param g Boolean, whether to remove graph attributes.
 * \param v Boolean, whether to remove vertex attributes.
 * \param e Boolean, whether to remove edge attributes.
 *
 * \sa \ref DELGAS, \ref DELVAS, \ref DELEAS, \ref DELALL for simpler
 * ways.
 */
void igraph_cattribute_remove_all(igraph_t *graph, igraph_bool_t g,
                                  igraph_bool_t v, igraph_bool_t e) {

    igraph_i_cattributes_t *attr = graph->attr;

    if (g) {
        igraph_attribute_record_list_clear(&attr->gal);
    }
    if (v) {
        igraph_attribute_record_list_clear(&attr->val);
    }
    if (e) {
        igraph_attribute_record_list_clear(&attr->eal);
    }
}
