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
#include "igraph_interface.h"
#include "igraph_random.h"

#include "internal/hacks.h" /* strdup */

#include <string.h>

/* An attribute is either a numeric vector (vector_t), a boolean vector
 * (vector_bool_t) or a string vector (strvector_t).
 * The attribute itself is stored in a struct igraph_attribute_record_t.
 * There is one such object for each attribute. The igraph_t has a pointer
 * to an igraph_i_cattribute_t, which contains three vector_ptr_t's, each
 * holding pointers to igraph_attribute_record_t objects. */

/* This function is used for producing better error messages. */
static const char *attribute_type_name(igraph_attribute_type_t type) {
    switch (type) {
    case IGRAPH_ATTRIBUTE_UNSPECIFIED:
        return "unspecified"; /* TODO: should probably trigger a fatal error */
    case IGRAPH_ATTRIBUTE_NUMERIC:
        return "numeric";
    case IGRAPH_ATTRIBUTE_BOOLEAN:
        return "boolean";
    case IGRAPH_ATTRIBUTE_STRING:
        return "string";
    case IGRAPH_ATTRIBUTE_OBJECT:
        return "object";
    }
    /* The following line is intentionally not in a default switch label
     * so that the compiler can warn about unhandled enum values,
     * should additional attribute types ever be added in the future. */
    IGRAPH_FATALF("Invalid attribute type %d found.", (int) type);
}

static igraph_bool_t igraph_i_cattribute_find(const igraph_vector_ptr_t *ptrvec,
                                              const char *name, igraph_integer_t *idx) {
    igraph_integer_t i, n = igraph_vector_ptr_size(ptrvec);
    igraph_bool_t l = false;
    for (i = 0; !l && i < n; i++) {
        igraph_attribute_record_t *rec = VECTOR(*ptrvec)[i];
        l = !strcmp(rec->name, name);
    }
    if (idx) {
        *idx = i - 1;
    }
    return l;
}

/*
 * Restores attribute vector lengths to their original size after a failure.
 * This function assumes that none of the attribute vectors are shorter than origlen.
 * Some may be longer due to a partially completed size extension: these will be
 * shrunk to their original size.
 */
static void igraph_i_cattribute_revert_attribute_vector_sizes(
        igraph_vector_ptr_t *attrlist, igraph_integer_t origlen) {

    igraph_integer_t no_of_attrs = igraph_vector_ptr_size(attrlist);
    for (igraph_integer_t i = 0; i < no_of_attrs; i++) {
        igraph_attribute_record_t *rec = VECTOR(*attrlist)[i];
        if (rec->type == IGRAPH_ATTRIBUTE_NUMERIC) {
            igraph_vector_t *nvec = (igraph_vector_t *) rec->value;
            IGRAPH_ASSERT(igraph_vector_capacity(nvec) >= origlen);
            igraph_vector_resize(nvec, origlen); /* shrinks */
        } else if (rec->type == IGRAPH_ATTRIBUTE_BOOLEAN) {
            igraph_vector_bool_t *bvec = (igraph_vector_bool_t *) rec->value;
            IGRAPH_ASSERT(igraph_vector_bool_capacity(bvec) >= origlen);
            igraph_vector_bool_resize(bvec, origlen); /* shrinks */
        } else if (rec->type == IGRAPH_ATTRIBUTE_STRING) {
            igraph_strvector_t *svec = (igraph_strvector_t *) rec->value;
            IGRAPH_ASSERT(igraph_strvector_capacity(svec) >= origlen);
            igraph_strvector_resize(svec, origlen); /* shrinks */
        } else {
            /* Must never reach here */
            IGRAPH_FATAL("Unknown attribute type encountered.");
        }
    }
}

typedef struct igraph_i_cattributes_t {
    igraph_vector_ptr_t gal;
    igraph_vector_ptr_t val;
    igraph_vector_ptr_t eal;
} igraph_i_cattributes_t;

static igraph_error_t igraph_i_cattributes_copy_attribute_record(igraph_attribute_record_t **newrec,
                                                      const igraph_attribute_record_t *rec) {
    igraph_vector_t *num, *newnum;
    igraph_strvector_t *str, *newstr;

    *newrec = IGRAPH_CALLOC(1, igraph_attribute_record_t);
    if (!(*newrec)) {
        IGRAPH_ERROR("Cannot copy attributes", IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */
    }
    IGRAPH_FINALLY(igraph_free, *newrec);
    (*newrec)->type = rec->type;
    (*newrec)->name = strdup(rec->name);
    if (!(*newrec)->name) {
        IGRAPH_ERROR("Cannot copy attributes", IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */
    }
    IGRAPH_FINALLY(igraph_free, (void*)(*newrec)->name);
    if (rec->type == IGRAPH_ATTRIBUTE_NUMERIC) {
        num = (igraph_vector_t *)rec->value;
        newnum = IGRAPH_CALLOC(1, igraph_vector_t);
        if (!newnum) {
            IGRAPH_ERROR("Cannot copy attributes", IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */
        }
        IGRAPH_FINALLY(igraph_free, newnum);
        IGRAPH_CHECK(igraph_vector_init_copy(newnum, num));
        IGRAPH_FINALLY(igraph_vector_destroy, newnum);
        (*newrec)->value = newnum;
    } else if (rec->type == IGRAPH_ATTRIBUTE_STRING) {
        str = (igraph_strvector_t*)rec->value;
        newstr = IGRAPH_CALLOC(1, igraph_strvector_t);
        if (!newstr) {
            IGRAPH_ERROR("Cannot copy attributes", IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */
        }
        IGRAPH_FINALLY(igraph_free, newstr);
        IGRAPH_CHECK(igraph_strvector_init_copy(newstr, str));
        IGRAPH_FINALLY(igraph_strvector_destroy, newstr);
        (*newrec)->value = newstr;
    } else if (rec->type == IGRAPH_ATTRIBUTE_BOOLEAN) {
        igraph_vector_bool_t *log = (igraph_vector_bool_t*) rec->value;
        igraph_vector_bool_t *newlog = IGRAPH_CALLOC(1, igraph_vector_bool_t);
        if (!newlog) {
            IGRAPH_ERROR("Cannot copy attributes", IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */
        }
        IGRAPH_FINALLY(igraph_free, newlog);
        IGRAPH_CHECK(igraph_vector_bool_init_copy(newlog, log));
        IGRAPH_FINALLY(igraph_vector_bool_destroy, newlog);
        (*newrec)->value = newlog;
    }

    IGRAPH_FINALLY_CLEAN(4);
    return IGRAPH_SUCCESS;
}

static void igraph_i_attribute_list_destroy(igraph_vector_ptr_t *attrlist) {
    igraph_integer_t i;
    igraph_integer_t n = igraph_vector_ptr_size(attrlist);
    for (i = 0; i < n; i++) {
        igraph_attribute_record_t *rec = VECTOR(*attrlist)[i];
        if (rec) {
            if (rec->type == IGRAPH_ATTRIBUTE_NUMERIC) {
                igraph_vector_t *num = (igraph_vector_t *) rec->value;
                igraph_vector_destroy(num);
                IGRAPH_FREE(num);
            } else if (rec->type == IGRAPH_ATTRIBUTE_STRING) {
                igraph_strvector_t *str = (igraph_strvector_t *) rec->value;
                igraph_strvector_destroy(str);
                IGRAPH_FREE(str);
            } else if (rec->type == IGRAPH_ATTRIBUTE_BOOLEAN) {
                igraph_vector_bool_t *boolvec = (igraph_vector_bool_t *) rec->value;
                igraph_vector_bool_destroy(boolvec);
                IGRAPH_FREE(boolvec);
            }
            IGRAPH_FREE(rec->name);
            IGRAPH_FREE(rec);
        }
    }
    igraph_vector_ptr_destroy(attrlist);
}

static igraph_error_t igraph_i_cattribute_init(igraph_t *graph, igraph_vector_ptr_t *attr) {
    igraph_attribute_record_t *attr_rec;
    igraph_integer_t i, n;
    igraph_i_cattributes_t *nattr;

    n = attr ? igraph_vector_ptr_size(attr) : 0;

    nattr = IGRAPH_CALLOC(1, igraph_i_cattributes_t);
    if (!nattr) {
        IGRAPH_ERROR("Can't init attributes", IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */
    }
    IGRAPH_FINALLY(igraph_free, nattr);

    IGRAPH_CHECK(igraph_vector_ptr_init(&nattr->gal, n));
    IGRAPH_FINALLY(igraph_i_attribute_list_destroy, &nattr->gal);
    IGRAPH_CHECK(igraph_vector_ptr_init(&nattr->val, 0));
    IGRAPH_FINALLY(igraph_vector_ptr_destroy, &nattr->val);
    IGRAPH_CHECK(igraph_vector_ptr_init(&nattr->eal, 0));
    IGRAPH_FINALLY(igraph_vector_ptr_destroy, &nattr->eal);

    for (i = 0; i < n; i++) {
        IGRAPH_CHECK(igraph_i_cattributes_copy_attribute_record(
                         &attr_rec, VECTOR(*attr)[i]));
        VECTOR(nattr->gal)[i] = attr_rec;
    }

    graph->attr = nattr;
    IGRAPH_FINALLY_CLEAN(4);

    return IGRAPH_SUCCESS;
}

static void igraph_i_cattribute_destroy(igraph_t *graph) {
    igraph_i_cattributes_t *attr = graph->attr;
    igraph_vector_ptr_t *als[3] = { &attr->gal, &attr->val, &attr->eal };
    for (size_t a = 0; a < 3; a++) {
        igraph_i_attribute_list_destroy(als[a]);
    }
    IGRAPH_FREE(graph->attr); /* sets to NULL */
}

/* Almost the same as destroy, but we might have null pointers */

static void igraph_i_cattribute_copy_free(igraph_i_cattributes_t *attr) {
    igraph_vector_ptr_t *als[3] = { &attr->gal, &attr->val, &attr->eal };
    igraph_integer_t i, n;
    igraph_vector_t *num;
    igraph_strvector_t *str;
    igraph_vector_bool_t *boolvec;
    igraph_attribute_record_t *rec;
    for (size_t a = 0; a < 3; a++) {
        n = igraph_vector_ptr_size(als[a]);
        for (i = 0; i < n; i++) {
            rec = VECTOR(*als[a])[i];
            if (!rec) {
                continue;
            }
            if (rec->type == IGRAPH_ATTRIBUTE_NUMERIC) {
                num = (igraph_vector_t*)rec->value;
                igraph_vector_destroy(num);
                IGRAPH_FREE(num);
            } else if (rec->type == IGRAPH_ATTRIBUTE_BOOLEAN) {
                boolvec = (igraph_vector_bool_t*)rec->value;
                igraph_vector_bool_destroy(boolvec);
                IGRAPH_FREE(boolvec);
            } else if (rec->type == IGRAPH_ATTRIBUTE_STRING) {
                str = (igraph_strvector_t*)rec->value;
                igraph_strvector_destroy(str);
                IGRAPH_FREE(str);
            }
            IGRAPH_FREE(rec->name);
            IGRAPH_FREE(rec);
        }
    }
}

/* No reference counting here. If you use attributes in C you should
   know what you're doing. */

static igraph_error_t igraph_i_cattribute_copy(igraph_t *to, const igraph_t *from,
                             igraph_bool_t ga, igraph_bool_t va, igraph_bool_t ea) {
    igraph_i_cattributes_t *attrfrom = from->attr, *attrto;
    igraph_vector_ptr_t *alto[3], *alfrom[3] = { &attrfrom->gal, &attrfrom->val,
                                                 &attrfrom->eal
                                               };
    igraph_integer_t i, n;
    igraph_bool_t copy[3] = { ga, va, ea };
    to->attr = attrto = IGRAPH_CALLOC(1, igraph_i_cattributes_t);
    if (!attrto) {
        IGRAPH_ERROR("Cannot copy attributes", IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */
    }
    IGRAPH_FINALLY(igraph_free, attrto);
    IGRAPH_VECTOR_PTR_INIT_FINALLY(&attrto->gal, 0);
    IGRAPH_VECTOR_PTR_INIT_FINALLY(&attrto->val, 0);
    IGRAPH_VECTOR_PTR_INIT_FINALLY(&attrto->eal, 0);
    IGRAPH_FINALLY_CLEAN(3);
    IGRAPH_FINALLY(igraph_i_cattribute_copy_free, attrto);

    alto[0] = &attrto->gal; alto[1] = &attrto->val; alto[2] = &attrto->eal;
    for (size_t a = 0; a < 3; a++) {
        if (copy[a]) {
            n = igraph_vector_ptr_size(alfrom[a]);
            IGRAPH_CHECK(igraph_vector_ptr_resize(alto[a], n));
            igraph_vector_ptr_null(alto[a]);
            for (i = 0; i < n; i++) {
                igraph_attribute_record_t *newrec;
                IGRAPH_CHECK(igraph_i_cattributes_copy_attribute_record(&newrec,
                             VECTOR(*alfrom[a])[i]));
                VECTOR(*alto[a])[i] = newrec;
            }
        }
    }

    IGRAPH_FINALLY_CLEAN(2);
    return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_i_cattribute_add_vertices_inner(igraph_t *graph, igraph_integer_t nv,
                                                         igraph_vector_ptr_t *nattr) {

    igraph_i_cattributes_t *attr = graph->attr;
    igraph_vector_ptr_t *val = &attr->val;
    igraph_integer_t length = igraph_vector_ptr_size(val);
    igraph_integer_t nattrno = nattr == NULL ? 0 : igraph_vector_ptr_size(nattr);
    igraph_integer_t origlen = igraph_vcount(graph) - nv;
    igraph_integer_t newattrs = 0, i;
    igraph_vector_int_t news;

    /* First add the new attributes if any */
    newattrs = 0;
    IGRAPH_VECTOR_INT_INIT_FINALLY(&news, 0);
    for (i = 0; i < nattrno; i++) {
        igraph_attribute_record_t *nattr_entry = VECTOR(*nattr)[i];
        const char *nname = nattr_entry->name;
        igraph_integer_t j;
        igraph_bool_t l = igraph_i_cattribute_find(val, nname, &j);
        if (!l) {
            newattrs++;
            IGRAPH_CHECK(igraph_vector_int_push_back(&news, i));
        } else {
            /* check types */
            if (nattr_entry->type !=
                ((igraph_attribute_record_t*)VECTOR(*val)[j])->type) {
                IGRAPH_ERROR("You cannot mix attribute types", IGRAPH_EINVAL);
            }
        }
    }

    /* Add NA/empty string vectors for the existing vertices */
    if (newattrs != 0) {
        for (i = 0; i < newattrs; i++) {
            igraph_attribute_record_t *tmp = VECTOR(*nattr)[VECTOR(news)[i]];
            igraph_attribute_record_t *newrec = IGRAPH_CALLOC(1, igraph_attribute_record_t);
            igraph_attribute_type_t type = tmp->type;
            if (!newrec) {
                IGRAPH_ERROR("Cannot add attributes", IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */
            }
            IGRAPH_FINALLY(igraph_free, newrec);
            newrec->type = type;
            newrec->name = strdup(tmp->name);
            if (!newrec->name) {
                IGRAPH_ERROR("Cannot add attributes", IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */
            }
            IGRAPH_FINALLY(igraph_free, (char*)newrec->name);
            if (type == IGRAPH_ATTRIBUTE_NUMERIC) {
                igraph_vector_t *newnum = IGRAPH_CALLOC(1, igraph_vector_t);
                if (!newnum) {
                    IGRAPH_ERROR("Cannot add attributes", IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */
                }
                IGRAPH_FINALLY(igraph_free, newnum);
                IGRAPH_VECTOR_INIT_FINALLY(newnum, origlen);
                newrec->value = newnum;
                igraph_vector_fill(newnum, IGRAPH_NAN);
            } else if (type == IGRAPH_ATTRIBUTE_STRING) {
                igraph_strvector_t *newstr = IGRAPH_CALLOC(1, igraph_strvector_t);
                if (!newstr) {
                    IGRAPH_ERROR("Cannot add attributes", IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */
                }
                IGRAPH_FINALLY(igraph_free, newstr);
                IGRAPH_STRVECTOR_INIT_FINALLY(newstr, origlen);
                newrec->value = newstr;
            } else if (type == IGRAPH_ATTRIBUTE_BOOLEAN) {
                igraph_vector_bool_t *newbool = IGRAPH_CALLOC(1, igraph_vector_bool_t);
                if (!newbool) {
                    IGRAPH_ERROR("Cannot add attributes", IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */
                }
                IGRAPH_FINALLY(igraph_free, newbool);
                IGRAPH_VECTOR_BOOL_INIT_FINALLY(newbool, origlen);
                newrec->value = newbool;
                igraph_vector_bool_fill(newbool, false);
            }
            IGRAPH_CHECK(igraph_vector_ptr_push_back(val, newrec));
            IGRAPH_FINALLY_CLEAN(4);
        }
        length = igraph_vector_ptr_size(val);
    }

    /* Now append the new values */
    for (i = 0; i < length; i++) {
        igraph_attribute_record_t *oldrec = VECTOR(*val)[i];
        igraph_attribute_record_t *newrec = 0;
        const char *name = oldrec->name;
        igraph_integer_t j = -1;
        igraph_bool_t l = false;
        if (nattr) {
            l = igraph_i_cattribute_find(nattr, name, &j);
        }
        if (l) {
            /* This attribute is present in nattr */
            igraph_vector_t *oldnum, *newnum;
            igraph_strvector_t *oldstr, *newstr;
            igraph_vector_bool_t *oldbool, *newbool;
            newrec = VECTOR(*nattr)[j];
            oldnum = (igraph_vector_t*)oldrec->value;
            newnum = (igraph_vector_t*)newrec->value;
            oldstr = (igraph_strvector_t*)oldrec->value;
            newstr = (igraph_strvector_t*)newrec->value;
            oldbool = (igraph_vector_bool_t*)oldrec->value;
            newbool = (igraph_vector_bool_t*)newrec->value;
            if (oldrec->type != newrec->type) {
                IGRAPH_ERROR("Attribute types do not match.", IGRAPH_EINVAL);
            }
            switch (oldrec->type) {
            case IGRAPH_ATTRIBUTE_NUMERIC:
                if (nv != igraph_vector_size(newnum)) {
                    IGRAPH_ERROR("Invalid numeric attribute length.", IGRAPH_EINVAL);
                }
                IGRAPH_CHECK(igraph_vector_append(oldnum, newnum));
                break;
            case IGRAPH_ATTRIBUTE_STRING:
                if (nv != igraph_strvector_size(newstr)) {
                    IGRAPH_ERROR("Invalid string attribute length.", IGRAPH_EINVAL);
                }
                IGRAPH_CHECK(igraph_strvector_append(oldstr, newstr));
                break;
            case IGRAPH_ATTRIBUTE_BOOLEAN:
                if (nv != igraph_vector_bool_size(newbool)) {
                    IGRAPH_ERROR("Invalid boolean attribute length.", IGRAPH_EINVAL);
                }
                IGRAPH_CHECK(igraph_vector_bool_append(oldbool, newbool));
                break;
            default:
                IGRAPH_WARNING("Invalid attribute type.");
                break;
            }
        } else {
            /* No such attribute, append NA's */
            igraph_vector_t *oldnum = (igraph_vector_t *)oldrec->value;
            igraph_strvector_t *oldstr = (igraph_strvector_t*)oldrec->value;
            igraph_vector_bool_t *oldbool = (igraph_vector_bool_t*)oldrec->value;
            switch (oldrec->type) {
            case IGRAPH_ATTRIBUTE_NUMERIC:
                IGRAPH_CHECK(igraph_vector_resize(oldnum, origlen + nv));
                for (j = origlen; j < origlen + nv; j++) {
                    VECTOR(*oldnum)[j] = IGRAPH_NAN;
                }
                break;
            case IGRAPH_ATTRIBUTE_STRING:
                IGRAPH_CHECK(igraph_strvector_resize(oldstr, origlen + nv));
                break;
            case IGRAPH_ATTRIBUTE_BOOLEAN:
                IGRAPH_CHECK(igraph_vector_bool_resize(oldbool, origlen + nv));
                for (j = origlen; j < origlen + nv; j++) {
                    VECTOR(*oldbool)[j] = 0;
                }
                break;
            default:
                IGRAPH_WARNING("Invalid attribute type");
                break;
            }
        }
    }

    igraph_vector_int_destroy(&news);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_i_cattribute_add_vertices(igraph_t *graph, igraph_integer_t nv,
                                                       igraph_vector_ptr_t *nattr) {
    /* Record information needed to restore attribute vector sizes */
    igraph_i_cattributes_t *attr = graph->attr;
    igraph_vector_ptr_t *val = &attr->val;
    igraph_integer_t origlen = igraph_vcount(graph) - nv;

    /* Attempt adding attributes */
    igraph_error_t err = igraph_i_cattribute_add_vertices_inner(graph, nv, nattr);
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

static void igraph_i_cattribute_clear_attribute_container(igraph_vector_ptr_t *v) {
    igraph_integer_t i, n = igraph_vector_ptr_size(v);
    for (i = 0; i < n; i++) {
        igraph_attribute_record_t *rec = VECTOR(*v)[i];
        IGRAPH_FREE(rec->name);
        if (rec->type == IGRAPH_ATTRIBUTE_NUMERIC) {
            igraph_vector_t *numv = (igraph_vector_t*) rec->value;
            igraph_vector_destroy(numv);
            IGRAPH_FREE(numv);
        } else if (rec->type == IGRAPH_ATTRIBUTE_STRING) {
            igraph_strvector_t *strv = (igraph_strvector_t*) rec->value;
            igraph_strvector_destroy(strv);
            IGRAPH_FREE(strv);
        } else if (rec->type == IGRAPH_ATTRIBUTE_BOOLEAN) {
            igraph_vector_bool_t *boolv = (igraph_vector_bool_t*) rec->value;
            igraph_vector_bool_destroy(boolv);
            IGRAPH_FREE(boolv);
        }
        IGRAPH_FREE(rec);
    }
    igraph_vector_ptr_clear(v);
}

typedef struct {
    igraph_vector_t *numeric;
    igraph_vector_bool_t *boolean;
    igraph_vector_ptr_t *strings;
    igraph_integer_t length;
} igraph_i_attribute_permutation_work_area_t;

static igraph_error_t igraph_i_attribute_permutation_work_area_init(
  igraph_i_attribute_permutation_work_area_t *work_area, igraph_integer_t length
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
        IGRAPH_CHECK_OOM(vec, "Cannot permute attributes");
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
        IGRAPH_CHECK_OOM(vec, "Cannot permute attributes");
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
        IGRAPH_CHECK_OOM(vec, "Cannot permute attributes");
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
    IGRAPH_CHECK_OOM(new_vec, "Cannot permute attributes");
    IGRAPH_FINALLY(igraph_free, new_vec);
    IGRAPH_CHECK(igraph_strvector_init(new_vec, 0));
    IGRAPH_FINALLY(igraph_strvector_destroy, new_vec);
    IGRAPH_CHECK(igraph_vector_ptr_push_back(work_area->strings, new_vec));
    IGRAPH_FINALLY_CLEAN(2);

    IGRAPH_CHECK(igraph_strvector_index(vec, new_vec, idx));

    return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_i_cattribute_permute_vertices_in_place(
    igraph_t *graph, const igraph_vector_int_t *idx
) {
    igraph_i_cattributes_t *attr = graph->attr;
    igraph_vector_ptr_t *val = &attr->val;
    igraph_integer_t valno = igraph_vector_ptr_size(val);
    igraph_integer_t i, j;
    igraph_attribute_record_t *oldrec;
    igraph_vector_t *num, *num_work;
    igraph_strvector_t *str, str_work;
    igraph_vector_bool_t *oldbool, *bool_work;
    igraph_i_attribute_permutation_work_area_t work_area;
    igraph_integer_t idx_size = igraph_vector_int_size(idx);

    /* shortcut: don't allocate anything if there are no attributes */
    if (valno == 0) {
        return IGRAPH_SUCCESS;
    }

    /* do all the allocations that can potentially fail before we actually
     * start to permute the vertices to ensure that we will not ever need to
     * back out from a permutation once we've started it */
    IGRAPH_CHECK(igraph_i_attribute_permutation_work_area_init(&work_area, idx_size));
    IGRAPH_FINALLY(igraph_i_attribute_permutation_work_area_destroy, &work_area);
    for (i = 0; i < valno; i++) {
        oldrec = VECTOR(*val)[i];
        switch (oldrec->type) {
        case IGRAPH_ATTRIBUTE_NUMERIC:
            num = (igraph_vector_t*) oldrec->value;
            IGRAPH_CHECK(igraph_vector_reserve(num, idx_size));
            IGRAPH_CHECK(igraph_i_attribute_permutation_work_area_alloc_for_numeric(&work_area));
            break;

        case IGRAPH_ATTRIBUTE_BOOLEAN:
            oldbool = (igraph_vector_bool_t*) oldrec->value;
            IGRAPH_CHECK(igraph_vector_bool_reserve(oldbool, idx_size));
            IGRAPH_CHECK(igraph_i_attribute_permutation_work_area_alloc_for_boolean(&work_area));
            break;

        case IGRAPH_ATTRIBUTE_STRING:
            str = (igraph_strvector_t*) oldrec->value;
            IGRAPH_CHECK(igraph_strvector_reserve(str, idx_size));
            IGRAPH_CHECK(igraph_i_attribute_permutation_work_area_alloc_for_strings(&work_area));
            break;

        default:
            IGRAPH_WARNING("Unknown vertex attribute ignored");
        }
    }

    /* let's do string attributes first because these might need extra
     * allocations that can fail. The strategy is to build new igraph_strvector_t
     * instances for the permuted attributes and store them in an
     * igraph_vector_ptr_t until we are done with all of them. If any of the
     * allocations fail, we can destroy the igraph_vector_ptr_t safely */
    for (i = 0; i < valno; i++) {
        oldrec = VECTOR(*val)[i];
        if (oldrec->type != IGRAPH_ATTRIBUTE_STRING) {
            continue;
        }

        str = (igraph_strvector_t*) oldrec->value;
        IGRAPH_CHECK(
            igraph_i_attribute_permutation_work_area_permute_and_store_strvector(
                &work_area, str, idx
            )
        );
    }

    /* strings are done, and now all vectors involved in the process are
     * as large as they should be (or larger) so the operations below are not
     * supposed to fail. We can safely replace the original string attribute
     * vectors with the permuted ones, and then proceed to the remaining
     * attributes */
    for (i = 0, j = 0; i < valno; i++) {
        oldrec = VECTOR(*val)[i];
        if (oldrec->type != IGRAPH_ATTRIBUTE_STRING) {
            continue;
        }

        str = (igraph_strvector_t*) oldrec->value;
        str_work = *((igraph_strvector_t*) VECTOR(*(work_area.strings))[j]);
        *((igraph_strvector_t*) VECTOR(*(work_area.strings))[j]) = *str;
        *str = str_work;
        j++;
    }
    igraph_i_attribute_permutation_work_area_release_stored_strvectors(&work_area);

    for (i = 0; i < valno; i++) {
        oldrec = VECTOR(*val)[i];
        switch (oldrec->type) {
        case IGRAPH_ATTRIBUTE_NUMERIC:
            num = (igraph_vector_t*) oldrec->value;
            num_work = work_area.numeric;
            IGRAPH_ASSERT(num_work != NULL);
            IGRAPH_CHECK(igraph_vector_index(num, num_work, idx));
            work_area.numeric = num;
            oldrec->value = num_work;
            break;
        case IGRAPH_ATTRIBUTE_BOOLEAN:
            oldbool = (igraph_vector_bool_t*) oldrec->value;
            bool_work = work_area.boolean;
            IGRAPH_ASSERT(bool_work != NULL);
            IGRAPH_CHECK(igraph_vector_bool_index(oldbool, bool_work, idx));
            work_area.boolean = oldbool;
            oldrec->value = bool_work;
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
    igraph_vector_ptr_t *val = &attr->val, *new_val = &new_attr->val;
    igraph_integer_t i, valno;

    IGRAPH_ASSERT(graph == newgraph || igraph_vector_ptr_empty(new_val));

    /* Handle in-place permutation separately */
    if (graph == newgraph) {
        return igraph_i_cattribute_permute_vertices_in_place(newgraph, idx);
    }

    /* New vertex attributes */
    valno = igraph_vector_ptr_size(val);
    IGRAPH_CHECK(igraph_vector_ptr_resize(new_val, valno));
    IGRAPH_FINALLY(igraph_i_cattribute_clear_attribute_container, new_val);

    for (i = 0; i < valno; i++) {
        igraph_attribute_record_t *oldrec = VECTOR(*val)[i];
        igraph_attribute_type_t type = oldrec->type;
        igraph_vector_t *num, *newnum;
        igraph_strvector_t *str, *newstr;
        igraph_vector_bool_t *oldbool, *newbool;

        /* The record itself */
        igraph_attribute_record_t *new_rec =
            IGRAPH_CALLOC(1, igraph_attribute_record_t);
        if (! new_rec) {
            IGRAPH_ERROR("Cannot create vertex attributes", IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */
        }
        IGRAPH_FINALLY(igraph_free, new_rec);
        new_rec->name = strdup(oldrec->name);
        if (! new_rec->name) {
            IGRAPH_ERROR("Cannot create vertex attributes", IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */
        }
        IGRAPH_FINALLY(igraph_free, (char *) new_rec->name);
        new_rec->type = oldrec->type;

        /* The data */
        switch (type) {
        case IGRAPH_ATTRIBUTE_NUMERIC:
            num = (igraph_vector_t*)oldrec->value;
            newnum = IGRAPH_CALLOC(1, igraph_vector_t);
            if (!newnum) {
                IGRAPH_ERROR("Cannot permute vertex attributes", IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */
            }
            IGRAPH_FINALLY(igraph_free, newnum);
            IGRAPH_VECTOR_INIT_FINALLY(newnum, 0);
            IGRAPH_CHECK(igraph_vector_index(num, newnum, idx));
            new_rec->value = newnum;
            break;
        case IGRAPH_ATTRIBUTE_BOOLEAN:
            oldbool = (igraph_vector_bool_t*)oldrec->value;
            newbool = IGRAPH_CALLOC(1, igraph_vector_bool_t);
            if (!newbool) {
                IGRAPH_ERROR("Cannot permute vertex attributes", IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */
            }
            IGRAPH_FINALLY(igraph_free, newbool);
            IGRAPH_VECTOR_BOOL_INIT_FINALLY(newbool, 0);
            IGRAPH_CHECK(igraph_vector_bool_index(oldbool, newbool, idx));
            new_rec->value = newbool;
            break;
        case IGRAPH_ATTRIBUTE_STRING:
            str = (igraph_strvector_t*)oldrec->value;
            newstr = IGRAPH_CALLOC(1, igraph_strvector_t);
            if (!newstr) {
                IGRAPH_ERROR("Cannot permute vertex attributes", IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */
            }
            IGRAPH_FINALLY(igraph_free, newstr);
            IGRAPH_STRVECTOR_INIT_FINALLY(newstr, 0);
            IGRAPH_CHECK(igraph_strvector_index(str, newstr, idx));
            new_rec->value = newstr;
            break;
        default:
            IGRAPH_WARNING("Unknown vertex attribute ignored");
        }

        VECTOR(*new_val)[i] = new_rec;
        IGRAPH_FINALLY_CLEAN(4);
    }

    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
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
    const igraph_vector_t *oldv = oldrec->value;
    igraph_vector_t *newv = IGRAPH_CALLOC(1, igraph_vector_t);
    igraph_integer_t newlen = igraph_vector_int_list_size(merges);
    igraph_integer_t i;

    if (!newv) {
        IGRAPH_ERROR("Cannot combine attributes", IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */
    }
    IGRAPH_FINALLY(igraph_free, newv);
    IGRAPH_VECTOR_INIT_FINALLY(newv, newlen);

    for (i = 0; i < newlen; i++) {
        igraph_real_t s = 0.0;
        igraph_vector_int_t *idx = igraph_vector_int_list_get_ptr(merges, i);;
        igraph_integer_t j, n = igraph_vector_int_size(idx);
        for (j = 0; j < n; j++) {
            igraph_integer_t x = VECTOR(*idx)[j];
            s += VECTOR(*oldv)[x];
        }
        VECTOR(*newv)[i] = s;
    }

    IGRAPH_FINALLY_CLEAN(2);
    newrec->value = newv;

    return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_i_cattributes_cn_prod(const igraph_attribute_record_t *oldrec,
                                        igraph_attribute_record_t * newrec,
                                        const igraph_vector_int_list_t *merges) {
    const igraph_vector_t *oldv = oldrec->value;
    igraph_vector_t *newv = IGRAPH_CALLOC(1, igraph_vector_t);
    igraph_integer_t newlen = igraph_vector_int_list_size(merges);
    igraph_integer_t i;

    if (!newv) {
        IGRAPH_ERROR("Cannot combine attributes", IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */
    }
    IGRAPH_FINALLY(igraph_free, newv);
    IGRAPH_VECTOR_INIT_FINALLY(newv, newlen);

    for (i = 0; i < newlen; i++) {
        igraph_real_t s = 1.0;
        igraph_vector_int_t *idx = igraph_vector_int_list_get_ptr(merges, i);;
        igraph_integer_t j, n = igraph_vector_int_size(idx);
        for (j = 0; j < n; j++) {
            igraph_integer_t x = VECTOR(*idx)[j];
            s *= VECTOR(*oldv)[x];
        }
        VECTOR(*newv)[i] = s;
    }

    IGRAPH_FINALLY_CLEAN(2);
    newrec->value = newv;

    return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_i_cattributes_cn_min(const igraph_attribute_record_t *oldrec,
                                       igraph_attribute_record_t * newrec,
                                       const igraph_vector_int_list_t *merges) {
    const igraph_vector_t *oldv = oldrec->value;
    igraph_vector_t *newv = IGRAPH_CALLOC(1, igraph_vector_t);
    igraph_integer_t newlen = igraph_vector_int_list_size(merges);
    igraph_integer_t i;

    if (!newv) {
        IGRAPH_ERROR("Cannot combine attributes", IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */
    }
    IGRAPH_FINALLY(igraph_free, newv);
    IGRAPH_VECTOR_INIT_FINALLY(newv, newlen);

    for (i = 0; i < newlen; i++) {
        igraph_vector_int_t *idx = igraph_vector_int_list_get_ptr(merges, i);;
        igraph_integer_t j, n = igraph_vector_int_size(idx);
        igraph_real_t m = n > 0 ? VECTOR(*oldv)[ VECTOR(*idx)[0] ] : IGRAPH_NAN;
        for (j = 1; j < n; j++) {
            igraph_integer_t x = VECTOR(*idx)[j];
            igraph_real_t val = VECTOR(*oldv)[x];
            if (val < m) {
                m = val;
            }
        }
        VECTOR(*newv)[i] = m;
    }

    IGRAPH_FINALLY_CLEAN(2);
    newrec->value = newv;

    return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_i_cattributes_cn_max(const igraph_attribute_record_t *oldrec,
                                       igraph_attribute_record_t * newrec,
                                       const igraph_vector_int_list_t *merges) {
    const igraph_vector_t *oldv = oldrec->value;
    igraph_vector_t *newv = IGRAPH_CALLOC(1, igraph_vector_t);
    igraph_integer_t newlen = igraph_vector_int_list_size(merges);
    igraph_integer_t i;

    if (!newv) {
        IGRAPH_ERROR("Cannot combine attributes", IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */
    }
    IGRAPH_FINALLY(igraph_free, newv);
    IGRAPH_VECTOR_INIT_FINALLY(newv, newlen);

    for (i = 0; i < newlen; i++) {
        igraph_vector_int_t *idx = igraph_vector_int_list_get_ptr(merges, i);;
        igraph_integer_t j, n = igraph_vector_int_size(idx);
        igraph_real_t m = n > 0 ? VECTOR(*oldv)[ VECTOR(*idx)[0] ] : IGRAPH_NAN;
        for (j = 1; j < n; j++) {
            igraph_integer_t x = VECTOR(*idx)[j];
            igraph_real_t val = VECTOR(*oldv)[x];
            if (val > m) {
                m = val;
            }
        }
        VECTOR(*newv)[i] = m;
    }

    IGRAPH_FINALLY_CLEAN(2);
    newrec->value = newv;

    return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_i_cattributes_cn_random(const igraph_attribute_record_t *oldrec,
                                          igraph_attribute_record_t * newrec,
                                          const igraph_vector_int_list_t *merges) {

    const igraph_vector_t *oldv = oldrec->value;
    igraph_vector_t *newv = IGRAPH_CALLOC(1, igraph_vector_t);
    igraph_integer_t newlen = igraph_vector_int_list_size(merges);
    igraph_integer_t i;

    if (!newv) {
        IGRAPH_ERROR("Cannot combine attributes", IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */
    }
    IGRAPH_FINALLY(igraph_free, newv);
    IGRAPH_VECTOR_INIT_FINALLY(newv, newlen);

    RNG_BEGIN();

    for (i = 0; i < newlen; i++) {
        igraph_vector_int_t *idx = igraph_vector_int_list_get_ptr(merges, i);;
        igraph_integer_t n = igraph_vector_int_size(idx);
        if (n == 0) {
            VECTOR(*newv)[i] = IGRAPH_NAN;
        } else if (n == 1) {
            VECTOR(*newv)[i] = VECTOR(*oldv)[ VECTOR(*idx)[0] ];
        } else {
            igraph_integer_t r = RNG_INTEGER(0, n - 1);
            VECTOR(*newv)[i] = VECTOR(*oldv)[ VECTOR(*idx)[r] ];
        }
    }

    RNG_END();

    IGRAPH_FINALLY_CLEAN(2);
    newrec->value = newv;

    return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_i_cattributes_cn_first(const igraph_attribute_record_t *oldrec,
                                         igraph_attribute_record_t * newrec,
                                         const igraph_vector_int_list_t *merges) {

    const igraph_vector_t *oldv = oldrec->value;
    igraph_vector_t *newv = IGRAPH_CALLOC(1, igraph_vector_t);
    igraph_integer_t newlen = igraph_vector_int_list_size(merges);
    igraph_integer_t i;

    if (!newv) {
        IGRAPH_ERROR("Cannot combine attributes", IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */
    }
    IGRAPH_FINALLY(igraph_free, newv);
    IGRAPH_VECTOR_INIT_FINALLY(newv, newlen);

    for (i = 0; i < newlen; i++) {
        igraph_vector_int_t *idx = igraph_vector_int_list_get_ptr(merges, i);;
        igraph_integer_t n = igraph_vector_int_size(idx);
        if (n == 0) {
            VECTOR(*newv)[i] = IGRAPH_NAN;
        } else {
            VECTOR(*newv)[i] = VECTOR(*oldv)[ VECTOR(*idx)[0] ];
        }
    }

    IGRAPH_FINALLY_CLEAN(2);
    newrec->value = newv;

    return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_i_cattributes_cn_last(const igraph_attribute_record_t *oldrec,
                                        igraph_attribute_record_t * newrec,
                                        const igraph_vector_int_list_t *merges) {

    const igraph_vector_t *oldv = oldrec->value;
    igraph_vector_t *newv = IGRAPH_CALLOC(1, igraph_vector_t);
    igraph_integer_t newlen = igraph_vector_int_list_size(merges);
    igraph_integer_t i;

    if (!newv) {
        IGRAPH_ERROR("Cannot combine attributes", IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */
    }
    IGRAPH_FINALLY(igraph_free, newv);
    IGRAPH_VECTOR_INIT_FINALLY(newv, newlen);

    for (i = 0; i < newlen; i++) {
        igraph_vector_int_t *idx = igraph_vector_int_list_get_ptr(merges, i);;
        igraph_integer_t n = igraph_vector_int_size(idx);
        if (n == 0) {
            VECTOR(*newv)[i] = IGRAPH_NAN;
        } else {
            VECTOR(*newv)[i] = VECTOR(*oldv)[ VECTOR(*idx)[n - 1] ];
        }
    }

    IGRAPH_FINALLY_CLEAN(2);
    newrec->value = newv;

    return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_i_cattributes_cn_mean(const igraph_attribute_record_t *oldrec,
                                        igraph_attribute_record_t * newrec,
                                        const igraph_vector_int_list_t *merges) {
    const igraph_vector_t *oldv = oldrec->value;
    igraph_vector_t *newv = IGRAPH_CALLOC(1, igraph_vector_t);
    igraph_integer_t newlen = igraph_vector_int_list_size(merges);
    igraph_integer_t i;

    if (!newv) {
        IGRAPH_ERROR("Cannot combine attributes", IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */
    }
    IGRAPH_FINALLY(igraph_free, newv);
    IGRAPH_VECTOR_INIT_FINALLY(newv, newlen);

    for (i = 0; i < newlen; i++) {
        igraph_vector_int_t *idx = igraph_vector_int_list_get_ptr(merges, i);;
        igraph_integer_t j, n = igraph_vector_int_size(idx);
        igraph_real_t s = n > 0 ? 0.0 : IGRAPH_NAN;
        for (j = 0; j < n; j++) {
            igraph_integer_t x = VECTOR(*idx)[j];
            s += VECTOR(*oldv)[x];
        }
        if (n > 0) {
            s = s / n;
        }
        VECTOR(*newv)[i] = s;
    }

    IGRAPH_FINALLY_CLEAN(2);
    newrec->value = newv;

    return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_i_cattributes_cn_func(const igraph_attribute_record_t *oldrec,
                                        igraph_attribute_record_t *newrec,
                                        const igraph_vector_int_list_t *merges,
                                        igraph_cattributes_combine_num_t *func) {

    const igraph_vector_t *oldv = oldrec->value;
    igraph_integer_t newlen = igraph_vector_int_list_size(merges);

    igraph_vector_t *newv = IGRAPH_CALLOC(1, igraph_vector_t);
    IGRAPH_CHECK_OOM(newv, "Cannot combine attributes.");
    IGRAPH_FINALLY(igraph_free, newv);
    IGRAPH_VECTOR_INIT_FINALLY(newv, newlen);

    igraph_vector_t values;
    IGRAPH_VECTOR_INIT_FINALLY(&values, 0);

    for (igraph_integer_t i = 0; i < newlen; i++) {
        igraph_vector_int_t *idx = igraph_vector_int_list_get_ptr(merges, i);;

        igraph_integer_t n = igraph_vector_int_size(idx);
        IGRAPH_CHECK(igraph_vector_resize(&values, n));
        for (igraph_integer_t j = 0; j < n; j++) {
            igraph_integer_t x = VECTOR(*idx)[j];
            VECTOR(values)[j] = VECTOR(*oldv)[x];
        }

        igraph_real_t res;
        IGRAPH_CHECK(func(&values, &res));
        VECTOR(*newv)[i] = res;
    }

    igraph_vector_destroy(&values);
    IGRAPH_FINALLY_CLEAN(3);
    newrec->value = newv;

    return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_i_cattributes_cb_random(const igraph_attribute_record_t *oldrec,
                                          igraph_attribute_record_t * newrec,
                                          const igraph_vector_int_list_t *merges) {

    const igraph_vector_bool_t *oldv = oldrec->value;
    igraph_vector_bool_t *newv = IGRAPH_CALLOC(1, igraph_vector_bool_t);
    igraph_integer_t newlen = igraph_vector_int_list_size(merges);
    igraph_integer_t i;

    if (!newv) {
        IGRAPH_ERROR("Cannot combine attributes", IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */
    }
    IGRAPH_FINALLY(igraph_free, newv);
    IGRAPH_VECTOR_BOOL_INIT_FINALLY(newv, newlen);

    RNG_BEGIN();

    for (i = 0; i < newlen; i++) {
        igraph_vector_int_t *idx = igraph_vector_int_list_get_ptr(merges, i);;
        igraph_integer_t n = igraph_vector_int_size(idx);
        if (n == 0) {
            VECTOR(*newv)[i] = 0;
        } else if (n == 1) {
            VECTOR(*newv)[i] = VECTOR(*oldv)[ VECTOR(*idx)[0] ];
        } else {
            igraph_integer_t r = RNG_INTEGER(0, n - 1);
            VECTOR(*newv)[i] = VECTOR(*oldv)[ VECTOR(*idx)[r] ];
        }
    }

    RNG_END();

    IGRAPH_FINALLY_CLEAN(2);
    newrec->value = newv;

    return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_i_cattributes_cb_first(const igraph_attribute_record_t *oldrec,
                                         igraph_attribute_record_t * newrec,
                                         const igraph_vector_int_list_t *merges) {

    const igraph_vector_bool_t *oldv = oldrec->value;
    igraph_vector_bool_t *newv = IGRAPH_CALLOC(1, igraph_vector_bool_t);
    igraph_integer_t newlen = igraph_vector_int_list_size(merges);
    igraph_integer_t i;

    if (!newv) {
        IGRAPH_ERROR("Cannot combine attributes", IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */
    }
    IGRAPH_FINALLY(igraph_free, newv);
    IGRAPH_VECTOR_BOOL_INIT_FINALLY(newv, newlen);

    for (i = 0; i < newlen; i++) {
        igraph_vector_int_t *idx = igraph_vector_int_list_get_ptr(merges, i);;
        igraph_integer_t n = igraph_vector_int_size(idx);
        if (n == 0) {
            VECTOR(*newv)[i] = 0;
        } else {
            VECTOR(*newv)[i] = VECTOR(*oldv)[ VECTOR(*idx)[0] ];
        }
    }

    IGRAPH_FINALLY_CLEAN(2);
    newrec->value = newv;

    return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_i_cattributes_cb_last(const igraph_attribute_record_t *oldrec,
                                        igraph_attribute_record_t * newrec,
                                        const igraph_vector_int_list_t *merges) {

    const igraph_vector_bool_t *oldv = oldrec->value;
    igraph_vector_bool_t *newv = IGRAPH_CALLOC(1, igraph_vector_bool_t);
    igraph_integer_t newlen = igraph_vector_int_list_size(merges);
    igraph_integer_t i;

    if (!newv) {
        IGRAPH_ERROR("Cannot combine attributes", IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */
    }
    IGRAPH_FINALLY(igraph_free, newv);
    IGRAPH_VECTOR_BOOL_INIT_FINALLY(newv, newlen);

    for (i = 0; i < newlen; i++) {
        igraph_vector_int_t *idx = igraph_vector_int_list_get_ptr(merges, i);;
        igraph_integer_t n = igraph_vector_int_size(idx);
        if (n == 0) {
            VECTOR(*newv)[i] = 0;
        } else {
            VECTOR(*newv)[i] = VECTOR(*oldv)[ VECTOR(*idx)[n - 1] ];
        }
    }

    IGRAPH_FINALLY_CLEAN(2);
    newrec->value = newv;

    return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_i_cattributes_cb_all_is_true(const igraph_attribute_record_t *oldrec,
                                               igraph_attribute_record_t * newrec,
                                               const igraph_vector_int_list_t *merges) {

    const igraph_vector_bool_t *oldv = oldrec->value;
    igraph_vector_bool_t *newv = IGRAPH_CALLOC(1, igraph_vector_bool_t);
    igraph_integer_t newlen = igraph_vector_int_list_size(merges);
    igraph_integer_t i, j, n, x;

    if (!newv) {
        IGRAPH_ERROR("Cannot combine attributes", IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */
    }
    IGRAPH_FINALLY(igraph_free, newv);
    IGRAPH_VECTOR_BOOL_INIT_FINALLY(newv, newlen);

    for (i = 0; i < newlen; i++) {
        igraph_vector_int_t *idx = igraph_vector_int_list_get_ptr(merges, i);;
        n = igraph_vector_int_size(idx);
        VECTOR(*newv)[i] = 1;
        for (j = 0; j < n; j++) {
            x = VECTOR(*idx)[j];
            if (!VECTOR(*oldv)[x]) {
                VECTOR(*newv)[i] = 0;
                break;
            }
        }
    }

    IGRAPH_FINALLY_CLEAN(2);
    newrec->value = newv;

    return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_i_cattributes_cb_any_is_true(const igraph_attribute_record_t *oldrec,
                                               igraph_attribute_record_t * newrec,
                                               const igraph_vector_int_list_t *merges) {

    const igraph_vector_bool_t *oldv = oldrec->value;
    igraph_vector_bool_t *newv = IGRAPH_CALLOC(1, igraph_vector_bool_t);
    igraph_integer_t newlen = igraph_vector_int_list_size(merges);
    igraph_integer_t i, j, n, x;

    if (!newv) {
        IGRAPH_ERROR("Cannot combine attributes", IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */
    }
    IGRAPH_FINALLY(igraph_free, newv);
    IGRAPH_VECTOR_BOOL_INIT_FINALLY(newv, newlen);

    for (i = 0; i < newlen; i++) {
        igraph_vector_int_t *idx = igraph_vector_int_list_get_ptr(merges, i);;
        n = igraph_vector_int_size(idx);
        VECTOR(*newv)[i] = 0;
        for (j = 0; j < n; j++) {
            x = VECTOR(*idx)[j];
            if (VECTOR(*oldv)[x]) {
                VECTOR(*newv)[i] = 1;
                break;
            }
        }
    }

    IGRAPH_FINALLY_CLEAN(2);
    newrec->value = newv;

    return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_i_cattributes_cb_majority(const igraph_attribute_record_t *oldrec,
                                            igraph_attribute_record_t * newrec,
                                            const igraph_vector_int_list_t *merges) {

    const igraph_vector_bool_t *oldv = oldrec->value;
    igraph_vector_bool_t *newv = IGRAPH_CALLOC(1, igraph_vector_bool_t);
    igraph_integer_t newlen = igraph_vector_int_list_size(merges);
    igraph_integer_t i, j, n, x, num_trues;

    if (!newv) {
        IGRAPH_ERROR("Cannot combine attributes", IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */
    }
    IGRAPH_FINALLY(igraph_free, newv);
    IGRAPH_VECTOR_BOOL_INIT_FINALLY(newv, newlen);

    RNG_BEGIN();

    for (i = 0; i < newlen; i++) {
        igraph_vector_int_t *idx = igraph_vector_int_list_get_ptr(merges, i);;

        n = igraph_vector_int_size(idx);

        num_trues = 0;
        for (j = 0; j < n; j++) {
            x = VECTOR(*idx)[j];
            if (VECTOR(*oldv)[x]) {
                num_trues++;
            }
        }

        if (n % 2 != 0) {
            VECTOR(*newv)[i] = (num_trues > n / 2);
        } else {
            if (num_trues == n / 2) {
                VECTOR(*newv)[i] = (RNG_UNIF01() < 0.5);
            } else {
                VECTOR(*newv)[i] = (num_trues > n / 2);
            }
        }
    }

    RNG_END();

    IGRAPH_FINALLY_CLEAN(2);
    newrec->value = newv;

    return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_i_cattributes_cb_func(const igraph_attribute_record_t *oldrec,
                                        igraph_attribute_record_t *newrec,
                                        const igraph_vector_int_list_t *merges,
                                        igraph_cattributes_combine_bool_t *func) {

    const igraph_vector_bool_t *oldv = oldrec->value;
    igraph_integer_t newlen = igraph_vector_int_list_size(merges);

    igraph_vector_bool_t *newv = IGRAPH_CALLOC(1, igraph_vector_bool_t);
    IGRAPH_CHECK_OOM(newv, "Cannot combine attributes.");
    IGRAPH_FINALLY(igraph_free, newv);
    IGRAPH_VECTOR_BOOL_INIT_FINALLY(newv, newlen);

    igraph_vector_bool_t values;
    IGRAPH_VECTOR_BOOL_INIT_FINALLY(&values, 0);

    for (igraph_integer_t i = 0; i < newlen; i++) {
        igraph_vector_int_t *idx = igraph_vector_int_list_get_ptr(merges, i);;

        igraph_integer_t n = igraph_vector_int_size(idx);
        IGRAPH_CHECK(igraph_vector_bool_resize(&values, n));
        for (igraph_integer_t j = 0; j < n; j++) {
            igraph_integer_t x = VECTOR(*idx)[j];
            VECTOR(values)[j] = VECTOR(*oldv)[x];
        }

        igraph_bool_t res;
        IGRAPH_CHECK(func(&values, &res));
        VECTOR(*newv)[i] = res;
    }

    igraph_vector_bool_destroy(&values);
    IGRAPH_FINALLY_CLEAN(3);
    newrec->value = newv;

    return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_i_cattributes_sn_random(const igraph_attribute_record_t *oldrec,
                                          igraph_attribute_record_t *newrec,
                                          const igraph_vector_int_list_t *merges) {

    const igraph_strvector_t *oldv = oldrec->value;
    igraph_integer_t newlen = igraph_vector_int_list_size(merges);
    igraph_integer_t i;
    igraph_strvector_t *newv = IGRAPH_CALLOC(1, igraph_strvector_t);

    if (!newv) {
        IGRAPH_ERROR("Cannot combine attributes", IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */
    }
    IGRAPH_FINALLY(igraph_free, newv);
    IGRAPH_STRVECTOR_INIT_FINALLY(newv, newlen);

    RNG_BEGIN();

    for (i = 0; i < newlen; i++) {
        igraph_vector_int_t *idx = igraph_vector_int_list_get_ptr(merges, i);;
        igraph_integer_t n = igraph_vector_int_size(idx);
        const char *tmp;
        if (n == 0) {
            IGRAPH_CHECK(igraph_strvector_set(newv, i, ""));
        } else if (n == 1) {
            tmp = igraph_strvector_get(oldv, 0);
            IGRAPH_CHECK(igraph_strvector_set(newv, i, tmp));
        } else {
            igraph_integer_t r = RNG_INTEGER(0, n - 1);
            tmp = igraph_strvector_get(oldv, r);
            IGRAPH_CHECK(igraph_strvector_set(newv, i, tmp));
        }
    }

    RNG_END();

    IGRAPH_FINALLY_CLEAN(2);
    newrec->value = newv;

    return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_i_cattributes_sn_first(const igraph_attribute_record_t *oldrec,
                                         igraph_attribute_record_t *newrec,
                                         const igraph_vector_int_list_t *merges) {

    const igraph_strvector_t *oldv = oldrec->value;
    igraph_integer_t i, newlen = igraph_vector_int_list_size(merges);
    igraph_strvector_t *newv = IGRAPH_CALLOC(1, igraph_strvector_t);

    if (!newv) {
        IGRAPH_ERROR("Cannot combine attributes", IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */
    }
    IGRAPH_FINALLY(igraph_free, newv);
    IGRAPH_STRVECTOR_INIT_FINALLY(newv, newlen);

    for (i = 0; i < newlen; i++) {
        igraph_vector_int_t *idx = igraph_vector_int_list_get_ptr(merges, i);;
        igraph_integer_t n = igraph_vector_int_size(idx);
        if (n == 0) {
            IGRAPH_CHECK(igraph_strvector_set(newv, i, ""));
        } else {
            const char *tmp = igraph_strvector_get(oldv, VECTOR(*idx)[0]);
            IGRAPH_CHECK(igraph_strvector_set(newv, i, tmp));
        }
    }

    IGRAPH_FINALLY_CLEAN(2);
    newrec->value = newv;

    return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_i_cattributes_sn_last(const igraph_attribute_record_t *oldrec,
                                        igraph_attribute_record_t *newrec,
                                        const igraph_vector_int_list_t *merges) {

    const igraph_strvector_t *oldv = oldrec->value;
    igraph_integer_t i, newlen = igraph_vector_int_list_size(merges);
    igraph_strvector_t *newv = IGRAPH_CALLOC(1, igraph_strvector_t);

    if (!newv) {
        IGRAPH_ERROR("Cannot combine attributes", IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */
    }
    IGRAPH_FINALLY(igraph_free, newv);
    IGRAPH_STRVECTOR_INIT_FINALLY(newv, newlen);

    for (i = 0; i < newlen; i++) {
        igraph_vector_int_t *idx = igraph_vector_int_list_get_ptr(merges, i);;
        igraph_integer_t n = igraph_vector_int_size(idx);
        if (n == 0) {
            IGRAPH_CHECK(igraph_strvector_set(newv, i, ""));
        } else {
            const char *tmp = igraph_strvector_get(oldv, VECTOR(*idx)[n - 1]);
            IGRAPH_CHECK(igraph_strvector_set(newv, i, tmp));
        }
    }

    IGRAPH_FINALLY_CLEAN(2);
    newrec->value = newv;

    return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_i_cattributes_sn_concat(const igraph_attribute_record_t *oldrec,
                                          igraph_attribute_record_t *newrec,
                                          const igraph_vector_int_list_t *merges) {

    const igraph_strvector_t *oldv = oldrec->value;
    igraph_integer_t i, newlen = igraph_vector_int_list_size(merges);
    igraph_strvector_t *newv = IGRAPH_CALLOC(1, igraph_strvector_t);

    if (!newv) {
        IGRAPH_ERROR("Cannot combine attributes", IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */
    }
    IGRAPH_FINALLY(igraph_free, newv);
    IGRAPH_STRVECTOR_INIT_FINALLY(newv, newlen);

    for (i = 0; i < newlen; i++) {
        igraph_vector_int_t *idx = igraph_vector_int_list_get_ptr(merges, i);;
        igraph_integer_t j, n = igraph_vector_int_size(idx);
        size_t len = 0;
        const char *tmp;
        char *tmp2;
        for (j = 0; j < n; j++) {
            tmp = igraph_strvector_get(oldv, j);
            len += strlen(tmp);
        }
        tmp2 = IGRAPH_CALLOC(len + 1, char);
        if (!tmp2) {
            IGRAPH_ERROR("Cannot combine attributes", IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */
        }
        IGRAPH_FINALLY(igraph_free, tmp2);
        len = 0;
        for (j = 0; j < n; j++) {
            tmp = igraph_strvector_get(oldv, j);
            strcpy(tmp2 + len, tmp);
            len += strlen(tmp);
        }

        IGRAPH_CHECK(igraph_strvector_set(newv, i, tmp2));
        IGRAPH_FREE(tmp2);
        IGRAPH_FINALLY_CLEAN(1);
    }

    IGRAPH_FINALLY_CLEAN(2);
    newrec->value = newv;

    return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_i_cattributes_sn_func(const igraph_attribute_record_t *oldrec,
                                        igraph_attribute_record_t *newrec,
                                        const igraph_vector_int_list_t *merges,
                                        igraph_cattributes_combine_str_t *func) {

    const igraph_strvector_t *oldv = oldrec->value;
    igraph_integer_t newlen = igraph_vector_int_list_size(merges);

    igraph_strvector_t *newv = IGRAPH_CALLOC(1, igraph_strvector_t);
    IGRAPH_CHECK_OOM(newv, "Cannot combine attributes.");
    IGRAPH_FINALLY(igraph_free, newv);
    IGRAPH_STRVECTOR_INIT_FINALLY(newv, newlen);

    igraph_strvector_t values;
    IGRAPH_STRVECTOR_INIT_FINALLY(&values, 0);

    for (igraph_integer_t i = 0; i < newlen; i++) {
        igraph_vector_int_t *idx = igraph_vector_int_list_get_ptr(merges, i);;

        igraph_integer_t n = igraph_vector_int_size(idx);
        IGRAPH_CHECK(igraph_strvector_resize(&values, n));
        for (igraph_integer_t j = 0; j < n; j++) {
            igraph_integer_t x = VECTOR(*idx)[j];
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
    IGRAPH_FINALLY_CLEAN(3);
    newrec->value = newv;

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

static igraph_error_t igraph_i_cattribute_combine_vertices(const igraph_t *graph,
                                                igraph_t *newgraph,
                                                const igraph_vector_int_list_t *merges,
                                                const igraph_attribute_combination_t *comb) {

    igraph_i_cattributes_t *attr = graph->attr;
    igraph_i_cattributes_t *toattr = newgraph->attr;
    igraph_vector_ptr_t *val = &attr->val;
    igraph_vector_ptr_t *new_val = &toattr->val;
    igraph_integer_t valno = igraph_vector_ptr_size(val);
    igraph_integer_t i, j, keepno = 0;
    igraph_attribute_combination_todo_item_t *todo_items;

    IGRAPH_ASSERT(graph != newgraph);
    IGRAPH_ASSERT(igraph_vector_ptr_empty(new_val));

    todo_items = IGRAPH_CALLOC(valno, igraph_attribute_combination_todo_item_t);
    if (!todo_items) {
        IGRAPH_ERROR("Cannot combine vertex attributes",
                     IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */
    }
    IGRAPH_FINALLY(igraph_free, todo_items);

    for (i = 0; i < valno; i++) {
        igraph_attribute_record_t *oldrec = VECTOR(*val)[i];
        const char *name = oldrec->name;
        igraph_attribute_combination_type_t type;
        igraph_function_pointer_t voidfunc;
        IGRAPH_CHECK(igraph_attribute_combination_query(comb, name, &type, &voidfunc));
        todo_items[i].type = type;
        todo_items[i].func.as_void = voidfunc;
        if (type != IGRAPH_ATTRIBUTE_COMBINE_IGNORE) {
            keepno++;
        }
    }

    IGRAPH_CHECK(igraph_vector_ptr_resize(new_val, keepno));
    IGRAPH_FINALLY(igraph_i_cattribute_clear_attribute_container, new_val);

    for (i = 0, j = 0; i < valno; i++) {
        igraph_attribute_record_t *newrec, *oldrec = VECTOR(*val)[i];
        const char *name = oldrec->name;
        igraph_attribute_combination_todo_item_t todo_item = todo_items[i];
        igraph_attribute_type_t attr_type = oldrec->type;

        if (todo_item.type == IGRAPH_ATTRIBUTE_COMBINE_DEFAULT ||
            todo_item.type == IGRAPH_ATTRIBUTE_COMBINE_IGNORE) {
            continue;
        }

        newrec = IGRAPH_CALLOC(1, igraph_attribute_record_t);
        if (!newrec) {
            IGRAPH_ERROR("Cannot combine vertex attributes", IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */
        }
        IGRAPH_FINALLY(igraph_free, newrec);
        newrec->name = strdup(name);
        if (!newrec->name) {
            IGRAPH_ERROR("Cannot combine vertex attributes", IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */
        }
        IGRAPH_FINALLY(igraph_free, (char *) newrec->name);
        newrec->type = attr_type;

        if (attr_type == IGRAPH_ATTRIBUTE_NUMERIC) {
            switch (todo_item.type) {
            case IGRAPH_ATTRIBUTE_COMBINE_FUNCTION:
                IGRAPH_CHECK(igraph_i_cattributes_cn_func(oldrec, newrec, merges,
                             todo_item.func.as_num));
                break;
            case IGRAPH_ATTRIBUTE_COMBINE_SUM:
                IGRAPH_CHECK(igraph_i_cattributes_cn_sum(oldrec, newrec, merges));
                break;
            case IGRAPH_ATTRIBUTE_COMBINE_PROD:
                IGRAPH_CHECK(igraph_i_cattributes_cn_prod(oldrec, newrec, merges));
                break;
            case IGRAPH_ATTRIBUTE_COMBINE_MIN:
                IGRAPH_CHECK(igraph_i_cattributes_cn_min(oldrec, newrec, merges));
                break;
            case IGRAPH_ATTRIBUTE_COMBINE_MAX:
                IGRAPH_CHECK(igraph_i_cattributes_cn_max(oldrec, newrec, merges));
                break;
            case IGRAPH_ATTRIBUTE_COMBINE_RANDOM:
                IGRAPH_CHECK(igraph_i_cattributes_cn_random(oldrec, newrec, merges));
                break;
            case IGRAPH_ATTRIBUTE_COMBINE_FIRST:
                IGRAPH_CHECK(igraph_i_cattributes_cn_first(oldrec, newrec, merges));
                break;
            case IGRAPH_ATTRIBUTE_COMBINE_LAST:
                IGRAPH_CHECK(igraph_i_cattributes_cn_last(oldrec, newrec, merges));
                break;
            case IGRAPH_ATTRIBUTE_COMBINE_MEAN:
                IGRAPH_CHECK(igraph_i_cattributes_cn_mean(oldrec, newrec, merges));
                break;
            case IGRAPH_ATTRIBUTE_COMBINE_MEDIAN:
                IGRAPH_ERROR("Median calculation not implemented",
                             IGRAPH_UNIMPLEMENTED);
                break;
            case IGRAPH_ATTRIBUTE_COMBINE_CONCAT:
                IGRAPH_ERROR("Cannot concatenate numeric attributes",
                             IGRAPH_EATTRCOMBINE);
                break;
            default:
                IGRAPH_ERROR("Unknown attribute_combination",
                             IGRAPH_UNIMPLEMENTED);
                break;
            }
        } else if (attr_type == IGRAPH_ATTRIBUTE_BOOLEAN) {
            switch (todo_item.type) {
            case IGRAPH_ATTRIBUTE_COMBINE_FUNCTION:
                IGRAPH_CHECK(igraph_i_cattributes_cb_func(oldrec, newrec, merges,
                             todo_item.func.as_bool));
                break;
            case IGRAPH_ATTRIBUTE_COMBINE_SUM:
            case IGRAPH_ATTRIBUTE_COMBINE_MAX:
                IGRAPH_CHECK(igraph_i_cattributes_cb_any_is_true(oldrec, newrec, merges));
                break;
            case IGRAPH_ATTRIBUTE_COMBINE_PROD:
            case IGRAPH_ATTRIBUTE_COMBINE_MIN:
                IGRAPH_CHECK(igraph_i_cattributes_cb_all_is_true(oldrec, newrec, merges));
                break;
            case IGRAPH_ATTRIBUTE_COMBINE_MEAN:
            case IGRAPH_ATTRIBUTE_COMBINE_MEDIAN:
                IGRAPH_CHECK(igraph_i_cattributes_cb_majority(oldrec, newrec, merges));
                break;
            case IGRAPH_ATTRIBUTE_COMBINE_RANDOM:
                IGRAPH_CHECK(igraph_i_cattributes_cb_random(oldrec, newrec, merges));
                break;
            case IGRAPH_ATTRIBUTE_COMBINE_FIRST:
                IGRAPH_CHECK(igraph_i_cattributes_cb_first(oldrec, newrec, merges));
                break;
            case IGRAPH_ATTRIBUTE_COMBINE_LAST:
                IGRAPH_CHECK(igraph_i_cattributes_cb_last(oldrec, newrec, merges));
                break;
            case IGRAPH_ATTRIBUTE_COMBINE_CONCAT:
                IGRAPH_ERROR("Cannot calculate concatenation of Booleans",
                             IGRAPH_EATTRCOMBINE);
                break;
            default:
                IGRAPH_ERROR("Unknown attribute_combination",
                             IGRAPH_UNIMPLEMENTED);
                break;
            }
        } else if (attr_type == IGRAPH_ATTRIBUTE_STRING) {
            switch (todo_item.type) {
            case IGRAPH_ATTRIBUTE_COMBINE_FUNCTION:
                IGRAPH_CHECK(igraph_i_cattributes_sn_func(oldrec, newrec, merges,
                             todo_item.func.as_str));
                break;
            case IGRAPH_ATTRIBUTE_COMBINE_SUM:
                IGRAPH_ERROR("Cannot sum strings", IGRAPH_EATTRCOMBINE);
                break;
            case IGRAPH_ATTRIBUTE_COMBINE_PROD:
                IGRAPH_ERROR("Cannot multiply strings", IGRAPH_EATTRCOMBINE);
                break;
            case IGRAPH_ATTRIBUTE_COMBINE_MIN:
                IGRAPH_ERROR("Cannot find minimum of strings",
                             IGRAPH_EATTRCOMBINE);
                break;
            case IGRAPH_ATTRIBUTE_COMBINE_MAX:
                IGRAPH_ERROR("Cannot find maximum of strings",
                             IGRAPH_EATTRCOMBINE);
                break;
            case IGRAPH_ATTRIBUTE_COMBINE_MEAN:
                IGRAPH_ERROR("Cannot calculate mean of strings",
                             IGRAPH_EATTRCOMBINE);
                break;
            case IGRAPH_ATTRIBUTE_COMBINE_MEDIAN:
                IGRAPH_ERROR("Cannot calculate median of strings",
                             IGRAPH_EATTRCOMBINE);
                break;
            case IGRAPH_ATTRIBUTE_COMBINE_RANDOM:
                IGRAPH_CHECK(igraph_i_cattributes_sn_random(oldrec, newrec, merges));
                break;
            case IGRAPH_ATTRIBUTE_COMBINE_FIRST:
                IGRAPH_CHECK(igraph_i_cattributes_sn_first(oldrec, newrec, merges));
                break;
            case IGRAPH_ATTRIBUTE_COMBINE_LAST:
                IGRAPH_CHECK(igraph_i_cattributes_sn_last(oldrec, newrec, merges));
                break;
            case IGRAPH_ATTRIBUTE_COMBINE_CONCAT:
                IGRAPH_CHECK(igraph_i_cattributes_sn_concat(oldrec, newrec, merges));
                break;
            default:
                IGRAPH_ERROR("Unknown attribute_combination",
                             IGRAPH_UNIMPLEMENTED);
                break;
            }
        } else {
            IGRAPH_ERROR("Unknown attribute type, this should not happen",
                         IGRAPH_UNIMPLEMENTED);
        }

        VECTOR(*new_val)[j] = newrec;
        IGRAPH_FINALLY_CLEAN(2); /* newrec->name and newrec */

        j++;
    }

    IGRAPH_FREE(todo_items);
    IGRAPH_FINALLY_CLEAN(2);

    return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_i_cattribute_add_edges_inner(igraph_t *graph, const igraph_vector_int_t *edges,
                                                      igraph_vector_ptr_t *nattr) {

    igraph_i_cattributes_t *attr = graph->attr;
    igraph_vector_ptr_t *eal = &attr->eal;
    igraph_integer_t ealno = igraph_vector_ptr_size(eal);
    igraph_integer_t ne = igraph_vector_int_size(edges) / 2;
    igraph_integer_t origlen = igraph_ecount(graph) - ne;
    igraph_integer_t nattrno = nattr == 0 ? 0 : igraph_vector_ptr_size(nattr);
    igraph_vector_int_t news;
    igraph_integer_t newattrs, i;

    /* First add the new attributes if any */
    newattrs = 0;
    IGRAPH_VECTOR_INT_INIT_FINALLY(&news, 0);
    for (i = 0; i < nattrno; i++) {
        igraph_attribute_record_t *nattr_entry = VECTOR(*nattr)[i];
        const char *nname = nattr_entry->name;
        igraph_integer_t j;
        igraph_bool_t l = igraph_i_cattribute_find(eal, nname, &j);
        if (!l) {
            newattrs++;
            IGRAPH_CHECK(igraph_vector_int_push_back(&news, i));
        } else {
            /* check types */
            if (nattr_entry->type !=
                ((igraph_attribute_record_t*)VECTOR(*eal)[j])->type) {
                IGRAPH_ERROR("You cannot mix attribute types", IGRAPH_EINVAL);
            }
        }
    }

    /* Add NaN/false/"" for the existing vertices for numeric, boolean and string attributes. */
    if (newattrs != 0) {
        for (i = 0; i < newattrs; i++) {
            igraph_attribute_record_t *tmp = VECTOR(*nattr)[ VECTOR(news)[i] ];
            igraph_attribute_record_t *newrec = IGRAPH_CALLOC(1, igraph_attribute_record_t);
            igraph_attribute_type_t type = tmp->type;
            if (!newrec) {
                IGRAPH_ERROR("Cannot add attributes", IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */
            }
            IGRAPH_FINALLY(igraph_free, newrec);
            newrec->type = type;
            newrec->name = strdup(tmp->name);
            if (!newrec->name) {
                IGRAPH_ERROR("Cannot add attributes", IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */
            }
            IGRAPH_FINALLY(igraph_free, (char*)newrec->name);
            if (type == IGRAPH_ATTRIBUTE_NUMERIC) {
                igraph_vector_t *newnum = IGRAPH_CALLOC(1, igraph_vector_t);
                if (!newnum) {
                    IGRAPH_ERROR("Cannot add attributes", IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */
                }
                IGRAPH_FINALLY(igraph_free, newnum);
                IGRAPH_VECTOR_INIT_FINALLY(newnum, origlen);
                newrec->value = newnum;
                igraph_vector_fill(newnum, IGRAPH_NAN);
            } else if (type == IGRAPH_ATTRIBUTE_BOOLEAN) {
                igraph_vector_bool_t *newbool = IGRAPH_CALLOC(1, igraph_vector_bool_t);
                if (!newbool) {
                    IGRAPH_ERROR("Cannot add attributes", IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */
                }
                IGRAPH_FINALLY(igraph_free, newbool);
                IGRAPH_VECTOR_BOOL_INIT_FINALLY(newbool, origlen);
                newrec->value = newbool;
                igraph_vector_bool_fill(newbool, false);
            } else if (type == IGRAPH_ATTRIBUTE_STRING) {
                igraph_strvector_t *newstr = IGRAPH_CALLOC(1, igraph_strvector_t);
                if (!newstr) {
                    IGRAPH_ERROR("Cannot add attributes", IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */
                }
                IGRAPH_FINALLY(igraph_free, newstr);
                IGRAPH_STRVECTOR_INIT_FINALLY(newstr, origlen);
                newrec->value = newstr;
            }
            IGRAPH_CHECK(igraph_vector_ptr_push_back(eal, newrec));
            IGRAPH_FINALLY_CLEAN(4);
        }
        ealno = igraph_vector_ptr_size(eal);
    }

    /* Now append the new values */
    for (i = 0; i < ealno; i++) {
        igraph_attribute_record_t *oldrec = VECTOR(*eal)[i];
        igraph_attribute_record_t *newrec = NULL;
        const char *name = oldrec->name;
        igraph_integer_t j = -1;
        igraph_bool_t l = false;
        if (nattr) {
            l = igraph_i_cattribute_find(nattr, name, &j);
        }
        if (l) {
            /* This attribute is present in nattr */
            igraph_vector_t *oldnum, *newnum;
            igraph_strvector_t *oldstr, *newstr;
            igraph_vector_bool_t *oldbool, *newbool;
            newrec = VECTOR(*nattr)[j];
            oldnum = (igraph_vector_t*)oldrec->value;
            newnum = (igraph_vector_t*)newrec->value;
            oldstr = (igraph_strvector_t*)oldrec->value;
            newstr = (igraph_strvector_t*)newrec->value;
            oldbool = (igraph_vector_bool_t*)oldrec->value;
            newbool = (igraph_vector_bool_t*)newrec->value;
            if (oldrec->type != newrec->type) {
                IGRAPH_ERROR("Attribute types do not match.", IGRAPH_EINVAL);
            }
            switch (oldrec->type) {
            case IGRAPH_ATTRIBUTE_NUMERIC:
                if (ne != igraph_vector_size(newnum)) {
                    IGRAPH_ERROR("Invalid numeric attribute length.", IGRAPH_EINVAL);
                }
                IGRAPH_CHECK(igraph_vector_append(oldnum, newnum));
                break;
            case IGRAPH_ATTRIBUTE_STRING:
                if (ne != igraph_strvector_size(newstr)) {
                    IGRAPH_ERROR("Invalid string attribute length.", IGRAPH_EINVAL);
                }
                IGRAPH_CHECK(igraph_strvector_append(oldstr, newstr));
                break;
            case IGRAPH_ATTRIBUTE_BOOLEAN:
                if (ne != igraph_vector_bool_size(newbool)) {
                    IGRAPH_ERROR("Invalid boolean attribute length.", IGRAPH_EINVAL);
                }
                IGRAPH_CHECK(igraph_vector_bool_append(oldbool, newbool));
                break;
            default:
                IGRAPH_WARNING("Invalid attribute type.");
                break;
            }
        } else {
            /* No such attribute, append NaN/false/"". */
            igraph_vector_t *oldnum = (igraph_vector_t *)oldrec->value;
            igraph_strvector_t *oldstr = (igraph_strvector_t*)oldrec->value;
            igraph_vector_bool_t *oldbool = (igraph_vector_bool_t *)oldrec->value;
            switch (oldrec->type) {
            case IGRAPH_ATTRIBUTE_NUMERIC:
                IGRAPH_CHECK(igraph_vector_resize(oldnum, origlen + ne));
                for (j = origlen; j < origlen + ne; j++) {
                    VECTOR(*oldnum)[j] = IGRAPH_NAN;
                }
                break;
            case IGRAPH_ATTRIBUTE_STRING:
                IGRAPH_CHECK(igraph_strvector_resize(oldstr, origlen + ne));
                break;
            case IGRAPH_ATTRIBUTE_BOOLEAN:
                IGRAPH_CHECK(igraph_vector_bool_resize(oldbool, origlen + ne));
                for (j = origlen; j < origlen + ne; j++) {
                    VECTOR(*oldbool)[j] = 0;
                }
                break;
            default:
                IGRAPH_WARNING("Invalid attribute type");
                break;
            }
        }
    }

    igraph_vector_int_destroy(&news);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_i_cattribute_add_edges(igraph_t *graph, const igraph_vector_int_t *edges,
                                                    igraph_vector_ptr_t *nattr) {
    /* Record information needed to restore attribute vector sizes */
    igraph_i_cattributes_t *attr = graph->attr;
    igraph_vector_ptr_t *eal = &attr->eal;
    igraph_integer_t ne = igraph_vector_int_size(edges) / 2;
    igraph_integer_t origlen = igraph_ecount(graph) - ne;

    /* Attempt adding attributes */
    igraph_error_t err = igraph_i_cattribute_add_edges_inner(graph, edges, nattr);
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
        igraph_i_cattribute_revert_attribute_vector_sizes(eal, origlen);
    }
    return err;
}

static igraph_error_t igraph_i_cattribute_permute_edges_in_place(
    igraph_t *graph, const igraph_vector_int_t *idx
) {
    igraph_i_cattributes_t *attr = graph->attr;
    igraph_vector_ptr_t *eal = &attr->eal;
    igraph_integer_t ealno = igraph_vector_ptr_size(eal);
    igraph_integer_t i, j;
    igraph_attribute_record_t *oldrec;
    igraph_vector_t *num, *num_work;
    igraph_strvector_t *str, str_work;
    igraph_vector_bool_t *oldbool, *bool_work;
    igraph_i_attribute_permutation_work_area_t work_area;
    igraph_integer_t idx_size = igraph_vector_int_size(idx);

    /* shortcut: don't allocate anything if there are no attributes */
    if (ealno == 0) {
        return IGRAPH_SUCCESS;
    }

    IGRAPH_CHECK(igraph_i_attribute_permutation_work_area_init(&work_area, idx_size));
    IGRAPH_FINALLY(igraph_i_attribute_permutation_work_area_destroy, &work_area);
    for (i = 0; i < ealno; i++) {
        oldrec = VECTOR(*eal)[i];
        switch (oldrec->type) {
        case IGRAPH_ATTRIBUTE_NUMERIC:
            num = (igraph_vector_t*) oldrec->value;
            IGRAPH_CHECK(igraph_vector_reserve(num, idx_size));
            IGRAPH_CHECK(igraph_i_attribute_permutation_work_area_alloc_for_numeric(&work_area));
            break;

        case IGRAPH_ATTRIBUTE_BOOLEAN:
            oldbool = (igraph_vector_bool_t*) oldrec->value;
            IGRAPH_CHECK(igraph_vector_bool_reserve(oldbool, idx_size));
            IGRAPH_CHECK(igraph_i_attribute_permutation_work_area_alloc_for_boolean(&work_area));
            break;

        case IGRAPH_ATTRIBUTE_STRING:
            str = (igraph_strvector_t*) oldrec->value;
            IGRAPH_CHECK(igraph_strvector_reserve(str, idx_size));
            IGRAPH_CHECK(igraph_i_attribute_permutation_work_area_alloc_for_strings(&work_area));
            break;

        default:
            IGRAPH_WARNING("Unknown edge attribute ignored");
        }
    }

    /* let's do string attributes first because these might need extra
     * allocations that can fail. The strategy is to build new igraph_strvector_t
     * instances for the permuted attributes and store them in an
     * igraph_vector_ptr_t until we are done with all of them. If any of the
     * allocations fail, we can destroy the igraph_vector_ptr_t safely */
    for (i = 0; i < ealno; i++) {
        oldrec = VECTOR(*eal)[i];
        if (oldrec->type != IGRAPH_ATTRIBUTE_STRING) {
            continue;
        }

        str = (igraph_strvector_t*) oldrec->value;
        IGRAPH_CHECK(
            igraph_i_attribute_permutation_work_area_permute_and_store_strvector(
                &work_area, str, idx
            )
        );
    }

    /* strings are done, and now all vectors involved in the process are
     * as large as they should be (or larger) so the operations below are not
     * supposed to fail. We can safely replace the original string attribute
     * vectors with the permuted ones, and then proceed to the remaining
     * attributes */
    for (i = 0, j = 0; i < ealno; i++) {
        oldrec = VECTOR(*eal)[i];
        if (oldrec->type != IGRAPH_ATTRIBUTE_STRING) {
            continue;
        }

        str = (igraph_strvector_t*) oldrec->value;
        str_work = *((igraph_strvector_t*) VECTOR(*(work_area.strings))[j]);
        *((igraph_strvector_t*) VECTOR(*(work_area.strings))[j]) = *str;
        *str = str_work;
        j++;
    }
    igraph_i_attribute_permutation_work_area_release_stored_strvectors(&work_area);

    /* now all vectors involved in the process are as large as they should be
     * (or larger) so the operations below are not supposed to fail -- except
     * for string operations that still do some extra allocations and we are
     * not prepared for the failures of those. This must still be fixed. */
    for (i = 0; i < ealno; i++) {
        oldrec = VECTOR(*eal)[i];
        switch (oldrec->type) {
        case IGRAPH_ATTRIBUTE_NUMERIC:
            num = (igraph_vector_t*) oldrec->value;
            num_work = work_area.numeric;
            IGRAPH_ASSERT(num_work != NULL);
            IGRAPH_CHECK(igraph_vector_index(num, num_work, idx));
            work_area.numeric = num;
            oldrec->value = num_work;
            break;
        case IGRAPH_ATTRIBUTE_BOOLEAN:
            oldbool = (igraph_vector_bool_t*) oldrec->value;
            bool_work = work_area.boolean;
            IGRAPH_ASSERT(bool_work != NULL);
            IGRAPH_CHECK(igraph_vector_bool_index(oldbool, bool_work, idx));
            work_area.boolean = oldbool;
            oldrec->value = bool_work;
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

static igraph_error_t igraph_i_cattribute_permute_edges(const igraph_t *graph,
                                             igraph_t *newgraph,
                                             const igraph_vector_int_t *idx) {

    igraph_i_cattributes_t *attr = graph->attr, *new_attr = newgraph->attr;
    igraph_vector_ptr_t *eal = &attr->eal, *new_eal = &new_attr->eal;
    igraph_integer_t i, ealno;

    IGRAPH_ASSERT(graph == newgraph || igraph_vector_ptr_empty(new_eal));

    if (graph == newgraph) {
        return igraph_i_cattribute_permute_edges_in_place(newgraph, idx);
    }

    /* New edge attributes */
    ealno = igraph_vector_ptr_size(eal);
    IGRAPH_ASSERT(igraph_vector_ptr_empty(new_eal));
    IGRAPH_CHECK(igraph_vector_ptr_resize(new_eal, ealno));
    IGRAPH_FINALLY(igraph_i_cattribute_clear_attribute_container, new_eal);

    for (i = 0; i < ealno; i++) {
        igraph_attribute_record_t *oldrec = VECTOR(*eal)[i];
        igraph_attribute_type_t type = oldrec->type;
        igraph_vector_t *num, *newnum;
        igraph_strvector_t *str, *newstr;
        igraph_vector_bool_t *oldbool, *newbool;

        /* The record itself */
        igraph_attribute_record_t *new_rec =
            IGRAPH_CALLOC(1, igraph_attribute_record_t);
        if (!new_rec) {
            IGRAPH_ERROR("Cannot create edge attributes", IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */
        }
        IGRAPH_FINALLY(igraph_free, new_rec);
        new_rec->name = strdup(oldrec->name);
        if (! new_rec->name) {
            IGRAPH_ERROR("Cannot create edge attributes", IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */
        }
        IGRAPH_FINALLY(igraph_free, (char *) new_rec->name);
        new_rec->type = oldrec->type;

        switch (type) {
        case IGRAPH_ATTRIBUTE_NUMERIC:
            num = (igraph_vector_t*) oldrec->value;
            newnum = IGRAPH_CALLOC(1, igraph_vector_t);
            if (!newnum) {
                IGRAPH_ERROR("Cannot permute edge attributes", IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */
            }
            IGRAPH_FINALLY(igraph_free, newnum);
            IGRAPH_VECTOR_INIT_FINALLY(newnum, 0);
            IGRAPH_CHECK(igraph_vector_index(num, newnum, idx));
            new_rec->value = newnum;
            break;
        case IGRAPH_ATTRIBUTE_STRING:
            str = (igraph_strvector_t*)oldrec->value;
            newstr = IGRAPH_CALLOC(1, igraph_strvector_t);
            if (!newstr) {
                IGRAPH_ERROR("Cannot permute edge attributes", IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */
            }
            IGRAPH_FINALLY(igraph_free, newstr);
            IGRAPH_STRVECTOR_INIT_FINALLY(newstr, 0);
            IGRAPH_CHECK(igraph_strvector_index(str, newstr, idx));
            new_rec->value = newstr;
            break;
        case IGRAPH_ATTRIBUTE_BOOLEAN:
            oldbool = (igraph_vector_bool_t*) oldrec->value;
            newbool = IGRAPH_CALLOC(1, igraph_vector_bool_t);
            if (!newbool) {
                IGRAPH_ERROR("Cannot permute edge attributes", IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */
            }
            IGRAPH_FINALLY(igraph_free, newbool);
            IGRAPH_VECTOR_BOOL_INIT_FINALLY(newbool, 0);
            IGRAPH_CHECK(igraph_vector_bool_index(oldbool, newbool, idx));
            new_rec->value = newbool;
            break;
        default:
            IGRAPH_WARNING("Unknown edge attribute ignored");
        }
        VECTOR(*new_eal)[i] = new_rec;
        IGRAPH_FINALLY_CLEAN(4);
    }
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_i_cattribute_combine_edges(const igraph_t *graph,
                                             igraph_t *newgraph,
                                             const igraph_vector_int_list_t *merges,
                                             const igraph_attribute_combination_t *comb) {

    igraph_i_cattributes_t *attr = graph->attr;
    igraph_i_cattributes_t *toattr = newgraph->attr;
    igraph_vector_ptr_t *eal = &attr->eal;
    igraph_vector_ptr_t *new_eal = &toattr->eal;
    igraph_integer_t ealno = igraph_vector_ptr_size(eal);
    igraph_integer_t i, j, keepno = 0;
    igraph_attribute_combination_todo_item_t *todo_items;

    IGRAPH_ASSERT(graph != newgraph);
    IGRAPH_ASSERT(igraph_vector_ptr_empty(new_eal));

    todo_items = IGRAPH_CALLOC(ealno, igraph_attribute_combination_todo_item_t);
    if (!todo_items) {
        IGRAPH_ERROR("Cannot combine edge attributes", IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */
    }
    IGRAPH_FINALLY(igraph_free, todo_items);

    for (i = 0; i < ealno; i++) {
        igraph_attribute_record_t *oldrec = VECTOR(*eal)[i];
        const char *name = oldrec->name;
        igraph_attribute_combination_type_t todo;
        igraph_function_pointer_t voidfunc;
        IGRAPH_CHECK(igraph_attribute_combination_query(comb, name, &todo, &voidfunc));
        todo_items[i].type = todo;
        todo_items[i].func.as_void = voidfunc;
        if (todo != IGRAPH_ATTRIBUTE_COMBINE_IGNORE) {
            keepno++;
        }
    }

    IGRAPH_CHECK(igraph_vector_ptr_resize(new_eal, keepno));
    IGRAPH_FINALLY(igraph_i_cattribute_clear_attribute_container, new_eal);

    for (i = 0, j = 0; i < ealno; i++) {
        igraph_attribute_record_t *newrec, *oldrec = VECTOR(*eal)[i];
        const char *name = oldrec->name;
        igraph_attribute_combination_todo_item_t todo_item = todo_items[i];
        igraph_attribute_type_t attr_type = oldrec->type;

        if (todo_item.type == IGRAPH_ATTRIBUTE_COMBINE_DEFAULT ||
            todo_item.type == IGRAPH_ATTRIBUTE_COMBINE_IGNORE) {
            continue;
        }

        newrec = IGRAPH_CALLOC(1, igraph_attribute_record_t);
        if (!newrec) {
            IGRAPH_ERROR("Cannot combine edge attributes", IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */
        }
        IGRAPH_FINALLY(igraph_free, newrec);
        newrec->name = strdup(name);
        if (! newrec->name) {
            IGRAPH_ERROR("Cannot combine edge attributes", IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */
        }
        IGRAPH_FINALLY(igraph_free, (char *) newrec->name);
        newrec->type = attr_type;

        if (attr_type == IGRAPH_ATTRIBUTE_NUMERIC) {
            switch (todo_item.type) {
            case IGRAPH_ATTRIBUTE_COMBINE_FUNCTION:
                IGRAPH_CHECK(igraph_i_cattributes_cn_func(oldrec, newrec, merges,
                             todo_item.func.as_num));
                break;
            case IGRAPH_ATTRIBUTE_COMBINE_SUM:
                IGRAPH_CHECK(igraph_i_cattributes_cn_sum(oldrec, newrec, merges));
                break;
            case IGRAPH_ATTRIBUTE_COMBINE_PROD:
                IGRAPH_CHECK(igraph_i_cattributes_cn_prod(oldrec, newrec, merges));
                break;
            case IGRAPH_ATTRIBUTE_COMBINE_MIN:
                IGRAPH_CHECK(igraph_i_cattributes_cn_min(oldrec, newrec, merges));
                break;
            case IGRAPH_ATTRIBUTE_COMBINE_MAX:
                IGRAPH_CHECK(igraph_i_cattributes_cn_max(oldrec, newrec, merges));
                break;
            case IGRAPH_ATTRIBUTE_COMBINE_RANDOM:
                IGRAPH_CHECK(igraph_i_cattributes_cn_random(oldrec, newrec, merges));
                break;
            case IGRAPH_ATTRIBUTE_COMBINE_FIRST:
                IGRAPH_CHECK(igraph_i_cattributes_cn_first(oldrec, newrec, merges));
                break;
            case IGRAPH_ATTRIBUTE_COMBINE_LAST:
                IGRAPH_CHECK(igraph_i_cattributes_cn_last(oldrec, newrec, merges));
                break;
            case IGRAPH_ATTRIBUTE_COMBINE_MEAN:
                IGRAPH_CHECK(igraph_i_cattributes_cn_mean(oldrec, newrec, merges));
                break;
            case IGRAPH_ATTRIBUTE_COMBINE_MEDIAN:
                IGRAPH_ERROR("Median calculation not implemented",
                             IGRAPH_UNIMPLEMENTED);
                break;
            case IGRAPH_ATTRIBUTE_COMBINE_CONCAT:
                IGRAPH_ERROR("Cannot concatenate numeric attributes",
                             IGRAPH_EATTRCOMBINE);
                break;
            default:
                IGRAPH_ERROR("Unknown attribute_combination",
                             IGRAPH_UNIMPLEMENTED);
                break;
            }
        } else if (attr_type == IGRAPH_ATTRIBUTE_BOOLEAN) {
            switch (todo_item.type) {
            case IGRAPH_ATTRIBUTE_COMBINE_FUNCTION:
                IGRAPH_CHECK(igraph_i_cattributes_cb_func(oldrec, newrec, merges,
                             todo_item.func.as_bool));
                break;
            case IGRAPH_ATTRIBUTE_COMBINE_SUM:
            case IGRAPH_ATTRIBUTE_COMBINE_MAX:
                IGRAPH_CHECK(igraph_i_cattributes_cb_any_is_true(oldrec, newrec, merges));
                break;
            case IGRAPH_ATTRIBUTE_COMBINE_PROD:
            case IGRAPH_ATTRIBUTE_COMBINE_MIN:
                IGRAPH_CHECK(igraph_i_cattributes_cb_all_is_true(oldrec, newrec, merges));
                break;
            case IGRAPH_ATTRIBUTE_COMBINE_MEAN:
            case IGRAPH_ATTRIBUTE_COMBINE_MEDIAN:
                IGRAPH_CHECK(igraph_i_cattributes_cb_majority(oldrec, newrec, merges));
                break;
            case IGRAPH_ATTRIBUTE_COMBINE_RANDOM:
                IGRAPH_CHECK(igraph_i_cattributes_cb_random(oldrec, newrec, merges));
                break;
            case IGRAPH_ATTRIBUTE_COMBINE_FIRST:
                IGRAPH_CHECK(igraph_i_cattributes_cb_first(oldrec, newrec, merges));
                break;
            case IGRAPH_ATTRIBUTE_COMBINE_LAST:
                IGRAPH_CHECK(igraph_i_cattributes_cb_last(oldrec, newrec, merges));
                break;
            case IGRAPH_ATTRIBUTE_COMBINE_CONCAT:
                IGRAPH_ERROR("Cannot calculate concatenation of Booleans",
                             IGRAPH_EATTRCOMBINE);
                break;
            default:
                IGRAPH_ERROR("Unknown attribute_combination",
                             IGRAPH_UNIMPLEMENTED);
                break;
            }
        } else if (attr_type == IGRAPH_ATTRIBUTE_STRING) {
            switch (todo_item.type) {
            case IGRAPH_ATTRIBUTE_COMBINE_FUNCTION:
                IGRAPH_CHECK(igraph_i_cattributes_sn_func(oldrec, newrec, merges,
                             todo_item.func.as_str));
                break;
            case IGRAPH_ATTRIBUTE_COMBINE_SUM:
                IGRAPH_ERROR("Cannot sum strings", IGRAPH_EATTRCOMBINE);
                break;
            case IGRAPH_ATTRIBUTE_COMBINE_PROD:
                IGRAPH_ERROR("Cannot multiply strings", IGRAPH_EATTRCOMBINE);
                break;
            case IGRAPH_ATTRIBUTE_COMBINE_MIN:
                IGRAPH_ERROR("Cannot find minimum of strings",
                             IGRAPH_EATTRCOMBINE);
                break;
            case IGRAPH_ATTRIBUTE_COMBINE_MAX:
                IGRAPH_ERROR("Cannot find maximum of strings",
                             IGRAPH_EATTRCOMBINE);
                break;
            case IGRAPH_ATTRIBUTE_COMBINE_MEAN:
                IGRAPH_ERROR("Cannot calculate mean of strings",
                             IGRAPH_EATTRCOMBINE);
                break;
            case IGRAPH_ATTRIBUTE_COMBINE_MEDIAN:
                IGRAPH_ERROR("Cannot calculate median of strings",
                             IGRAPH_EATTRCOMBINE);
                break;
            case IGRAPH_ATTRIBUTE_COMBINE_RANDOM:
                IGRAPH_CHECK(igraph_i_cattributes_sn_random(oldrec, newrec, merges));
                break;
            case IGRAPH_ATTRIBUTE_COMBINE_FIRST:
                IGRAPH_CHECK(igraph_i_cattributes_sn_first(oldrec, newrec, merges));
                break;
            case IGRAPH_ATTRIBUTE_COMBINE_LAST:
                IGRAPH_CHECK(igraph_i_cattributes_sn_last(oldrec, newrec, merges));
                break;
            case IGRAPH_ATTRIBUTE_COMBINE_CONCAT:
                IGRAPH_CHECK(igraph_i_cattributes_sn_concat(oldrec, newrec, merges));
                break;
            default:
                IGRAPH_ERROR("Unknown attribute_combination",
                             IGRAPH_UNIMPLEMENTED);
                break;
            }
        } else {
            IGRAPH_ERROR("Unknown attribute type, this should not happen",
                         IGRAPH_UNIMPLEMENTED);
        }

        VECTOR(*new_eal)[j] = newrec;
        IGRAPH_FINALLY_CLEAN(2); /* newrec and newrc->name */

        j++;
    }

    IGRAPH_FREE(todo_items);
    IGRAPH_FINALLY_CLEAN(2);

    return IGRAPH_SUCCESS;
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
    igraph_vector_ptr_t *attr[3] = { &at->gal, &at->val, &at->eal };
    igraph_integer_t i, j;

    for (i = 0; i < 3; i++) {
        igraph_strvector_t *n = names[i];
        igraph_vector_int_t *t = types[i];
        igraph_vector_ptr_t *al = attr[i];
        igraph_integer_t len = igraph_vector_ptr_size(al);

        if (n) {
            IGRAPH_CHECK(igraph_strvector_resize(n, len));
        }
        if (t) {
            IGRAPH_CHECK(igraph_vector_int_resize(t, len));
        }

        for (j = 0; j < len; j++) {
            igraph_attribute_record_t *rec = VECTOR(*al)[j];
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
    igraph_i_cattributes_t *at = graph->attr;
    igraph_vector_ptr_t *attr[3] = { &at->gal, &at->val, &at->eal };
    igraph_integer_t attrnum;

    switch (type) {
    case IGRAPH_ATTRIBUTE_GRAPH:
        attrnum = 0;
        break;
    case IGRAPH_ATTRIBUTE_VERTEX:
        attrnum = 1;
        break;
    case IGRAPH_ATTRIBUTE_EDGE:
        attrnum = 2;
        break;
    default:
        IGRAPH_ERROR("Unknown attribute element type", IGRAPH_EINVAL);
        break;
    }

    return igraph_i_cattribute_find(attr[attrnum], name, 0);
}

static igraph_error_t igraph_i_cattribute_gettype(const igraph_t *graph,
                                       igraph_attribute_type_t *type,
                                       igraph_attribute_elemtype_t elemtype,
                                       const char *name) {
    igraph_integer_t attrnum;
    igraph_attribute_record_t *rec;
    igraph_i_cattributes_t *at = graph->attr;
    igraph_vector_ptr_t *attr[3] = { &at->gal, &at->val, &at->eal };
    igraph_vector_ptr_t *al;
    igraph_integer_t j;
    igraph_bool_t l = false;

    switch (elemtype) {
    case IGRAPH_ATTRIBUTE_GRAPH:
        attrnum = 0;
        break;
    case IGRAPH_ATTRIBUTE_VERTEX:
        attrnum = 1;
        break;
    case IGRAPH_ATTRIBUTE_EDGE:
        attrnum = 2;
        break;
    default:
        IGRAPH_ERROR("Unknown attribute element type", IGRAPH_EINVAL);
        break;
    }

    al = attr[attrnum];
    l = igraph_i_cattribute_find(al, name, &j);
    if (!l) {
        IGRAPH_ERROR("Unknown attribute", IGRAPH_EINVAL);
    }
    rec = VECTOR(*al)[j];
    *type = rec->type;

    return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_i_cattribute_get_numeric_graph_attr(const igraph_t *graph,
                                                      const char *name,
                                                      igraph_vector_t *value) {
    igraph_i_cattributes_t *attr = graph->attr;
    igraph_vector_ptr_t *gal = &attr->gal;
    igraph_integer_t j;
    igraph_attribute_record_t *rec;
    igraph_vector_t *num;
    igraph_bool_t l = igraph_i_cattribute_find(gal, name, &j);

    if (!l) {
        IGRAPH_ERRORF("The graph attribute '%s' does not exist.", IGRAPH_EINVAL, name);
    }

    rec = VECTOR(*gal)[j];
    if (rec->type != IGRAPH_ATTRIBUTE_NUMERIC) {
        IGRAPH_ERRORF("Numeric graph attribute '%s' expected, got %s.", IGRAPH_EINVAL, name, attribute_type_name(rec->type));
    }
    num = (igraph_vector_t*)rec->value;
    IGRAPH_CHECK(igraph_vector_resize(value, 1));
    VECTOR(*value)[0] = VECTOR(*num)[0];

    return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_i_cattribute_get_bool_graph_attr(const igraph_t *graph,
                                                   const char *name,
                                                   igraph_vector_bool_t *value) {
    igraph_i_cattributes_t *attr = graph->attr;
    igraph_vector_ptr_t *gal = &attr->gal;
    igraph_integer_t j;
    igraph_attribute_record_t *rec;
    igraph_vector_bool_t *log;
    igraph_bool_t l = igraph_i_cattribute_find(gal, name, &j);

    if (!l) {
        IGRAPH_ERRORF("The graph attribute '%s' does not exist.", IGRAPH_EINVAL, name);
    }

    rec = VECTOR(*gal)[j];
    if (rec->type != IGRAPH_ATTRIBUTE_BOOLEAN) {
        IGRAPH_ERRORF("Boolean graph attribute '%s' expected, got %s.", IGRAPH_EINVAL, name, attribute_type_name(rec->type));
    }
    log = (igraph_vector_bool_t*)rec->value;
    IGRAPH_CHECK(igraph_vector_bool_resize(value, 1));
    VECTOR(*value)[0] = VECTOR(*log)[0];

    return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_i_cattribute_get_string_graph_attr(const igraph_t *graph,
                                                     const char *name,
                                                     igraph_strvector_t *value) {
    igraph_i_cattributes_t *attr = graph->attr;
    igraph_vector_ptr_t *gal = &attr->gal;
    igraph_integer_t j;
    igraph_attribute_record_t *rec;
    igraph_strvector_t *str;
    igraph_bool_t l = igraph_i_cattribute_find(gal, name, &j);

    if (!l) {
        IGRAPH_ERRORF("The graph attribute '%s' does not exist.", IGRAPH_EINVAL, name);
    }

    rec = VECTOR(*gal)[j];
    if (rec->type != IGRAPH_ATTRIBUTE_STRING) {
        IGRAPH_ERRORF("String graph attribute '%s' expected, got %s.", IGRAPH_EINVAL, name, attribute_type_name(rec->type));
    }
    str = (igraph_strvector_t*)rec->value;
    IGRAPH_CHECK(igraph_strvector_resize(value, 1));
    IGRAPH_CHECK(igraph_strvector_set(value, 0, igraph_strvector_get(str, 0)));

    return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_i_cattribute_get_numeric_vertex_attr(const igraph_t *graph,
                                                       const char *name,
                                                       igraph_vs_t vs,
                                                       igraph_vector_t *value) {
    igraph_i_cattributes_t *attr = graph->attr;
    igraph_vector_ptr_t *val = &attr->val;
    igraph_integer_t j;
    igraph_attribute_record_t *rec;
    igraph_vector_t *num;
    igraph_bool_t l = igraph_i_cattribute_find(val, name, &j);

    if (!l) {
        IGRAPH_ERRORF("The vertex attribute '%s' does not exist.", IGRAPH_EINVAL, name);
    }

    rec = VECTOR(*val)[j];
    if (rec->type != IGRAPH_ATTRIBUTE_NUMERIC) {
        IGRAPH_ERRORF("Numeric vertex attribute '%s' expected, got %s.", IGRAPH_EINVAL, name, attribute_type_name(rec->type));
    }
    num = (igraph_vector_t*)rec->value;
    if (igraph_vs_is_all(&vs)) {
        igraph_vector_clear(value);
        IGRAPH_CHECK(igraph_vector_append(value, num));
    } else {
        igraph_vit_t it;
        igraph_integer_t i = 0;
        IGRAPH_CHECK(igraph_vit_create(graph, vs, &it));
        IGRAPH_FINALLY(igraph_vit_destroy, &it);
        IGRAPH_CHECK(igraph_vector_resize(value, IGRAPH_VIT_SIZE(it)));
        for (; !IGRAPH_VIT_END(it); IGRAPH_VIT_NEXT(it), i++) {
            igraph_integer_t v = IGRAPH_VIT_GET(it);
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
    igraph_vector_ptr_t *val = &attr->val;
    igraph_vit_t it;
    igraph_integer_t i, j, v;
    igraph_attribute_record_t *rec;
    igraph_vector_bool_t *log;
    igraph_bool_t l = igraph_i_cattribute_find(val, name, &j);

    if (!l) {
        IGRAPH_ERRORF("The vertex attribute '%s' does not exist.", IGRAPH_EINVAL, name);
    }

    rec = VECTOR(*val)[j];
    if (rec->type != IGRAPH_ATTRIBUTE_BOOLEAN) {
        IGRAPH_ERRORF("Boolean vertex attribute '%s' expected, got %s.", IGRAPH_EINVAL, name, attribute_type_name(rec->type));
    }
    log = (igraph_vector_bool_t*)rec->value;
    if (igraph_vs_is_all(&vs)) {
        igraph_vector_bool_clear(value);
        IGRAPH_CHECK(igraph_vector_bool_append(value, log));
    } else {
        IGRAPH_CHECK(igraph_vit_create(graph, vs, &it));
        IGRAPH_FINALLY(igraph_vit_destroy, &it);
        IGRAPH_CHECK(igraph_vector_bool_resize(value, IGRAPH_VIT_SIZE(it)));
        for (i = 0; !IGRAPH_VIT_END(it); IGRAPH_VIT_NEXT(it), i++) {
            v = IGRAPH_VIT_GET(it);
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
    igraph_vector_ptr_t *val = &attr->val;
    igraph_integer_t j;
    igraph_attribute_record_t *rec;
    igraph_strvector_t *str;
    igraph_bool_t l = igraph_i_cattribute_find(val, name, &j);

    if (!l) {
        IGRAPH_ERRORF("The vertex attribute '%s' does not exist.", IGRAPH_EINVAL, name);
    }

    rec = VECTOR(*val)[j];
    if (rec->type != IGRAPH_ATTRIBUTE_STRING) {
        IGRAPH_ERRORF("String vertex attribute '%s' expected, got %s.", IGRAPH_EINVAL, name, attribute_type_name(rec->type));
    }
    str = (igraph_strvector_t*)rec->value;
    if (igraph_vs_is_all(&vs)) {
        igraph_strvector_clear(value);
        IGRAPH_CHECK(igraph_strvector_append(value, str));
    } else {
        igraph_vit_t it;
        igraph_integer_t i = 0;
        IGRAPH_CHECK(igraph_vit_create(graph, vs, &it));
        IGRAPH_FINALLY(igraph_vit_destroy, &it);
        IGRAPH_CHECK(igraph_strvector_resize(value, IGRAPH_VIT_SIZE(it)));
        for (; !IGRAPH_VIT_END(it); IGRAPH_VIT_NEXT(it), i++) {
            igraph_integer_t v = IGRAPH_VIT_GET(it);
            const char *s = igraph_strvector_get(str, v);
            IGRAPH_CHECK(igraph_strvector_set(value, i, s));
        }
        igraph_vit_destroy(&it);
        IGRAPH_FINALLY_CLEAN(1);
    }

    return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_i_cattribute_get_numeric_edge_attr(const igraph_t *graph,
                                                     const char *name,
                                                     igraph_es_t es,
                                                     igraph_vector_t *value) {
    igraph_i_cattributes_t *attr = graph->attr;
    igraph_vector_ptr_t *eal = &attr->eal;
    igraph_integer_t j;
    igraph_attribute_record_t *rec;
    igraph_vector_t *num;
    igraph_bool_t l = igraph_i_cattribute_find(eal, name, &j);

    if (!l) {
        IGRAPH_ERRORF("The edge attribute '%s' does not exist.", IGRAPH_EINVAL, name);
    }

    rec = VECTOR(*eal)[j];
    if (rec->type != IGRAPH_ATTRIBUTE_NUMERIC) {
        IGRAPH_ERRORF("Numeric edge attribute '%s' expected, got %s.", IGRAPH_EINVAL, name, attribute_type_name(rec->type));
    }
    num = (igraph_vector_t*)rec->value;
    if (igraph_es_is_all(&es)) {
        igraph_vector_clear(value);
        IGRAPH_CHECK(igraph_vector_append(value, num));
    } else {
        igraph_eit_t it;
        igraph_integer_t i = 0;
        IGRAPH_CHECK(igraph_eit_create(graph, es, &it));
        IGRAPH_FINALLY(igraph_eit_destroy, &it);
        IGRAPH_CHECK(igraph_vector_resize(value, IGRAPH_EIT_SIZE(it)));
        for (; !IGRAPH_EIT_END(it); IGRAPH_EIT_NEXT(it), i++) {
            igraph_integer_t e = IGRAPH_EIT_GET(it);
            VECTOR(*value)[i] = VECTOR(*num)[e];
        }
        igraph_eit_destroy(&it);
        IGRAPH_FINALLY_CLEAN(1);
    }

    return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_i_cattribute_get_string_edge_attr(const igraph_t *graph,
                                                    const char *name,
                                                    igraph_es_t es,
                                                    igraph_strvector_t *value) {
    igraph_i_cattributes_t *attr = graph->attr;
    igraph_vector_ptr_t *eal = &attr->eal;
    igraph_integer_t j;
    igraph_attribute_record_t *rec;
    igraph_strvector_t *str;
    igraph_bool_t l = igraph_i_cattribute_find(eal, name, &j);

    if (!l) {
        IGRAPH_ERRORF("The edge attribute '%s' does not exist.", IGRAPH_EINVAL, name);
    }

    rec = VECTOR(*eal)[j];
    if (rec->type != IGRAPH_ATTRIBUTE_STRING) {
        IGRAPH_ERRORF("String edge attribute '%s' expected, got %s.", IGRAPH_EINVAL, name, attribute_type_name(rec->type));
    }
    str = (igraph_strvector_t*)rec->value;
    if (igraph_es_is_all(&es)) {
        igraph_strvector_clear(value);
        IGRAPH_CHECK(igraph_strvector_append(value, str));
    } else {
        igraph_eit_t it;
        igraph_integer_t i = 0;
        IGRAPH_CHECK(igraph_eit_create(graph, es, &it));
        IGRAPH_FINALLY(igraph_eit_destroy, &it);
        IGRAPH_CHECK(igraph_strvector_resize(value, IGRAPH_EIT_SIZE(it)));
        for (; !IGRAPH_EIT_END(it); IGRAPH_EIT_NEXT(it), i++) {
            igraph_integer_t e = IGRAPH_EIT_GET(it);
            const char *s = igraph_strvector_get(str, e);
            IGRAPH_CHECK(igraph_strvector_set(value, i, s));
        }
        igraph_eit_destroy(&it);
        IGRAPH_FINALLY_CLEAN(1);
    }

    return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_i_cattribute_get_bool_edge_attr(const igraph_t *graph,
                                                  const char *name,
                                                  igraph_es_t es,
                                                  igraph_vector_bool_t *value) {
    igraph_i_cattributes_t *attr = graph->attr;
    igraph_vector_ptr_t *eal = &attr->eal;
    igraph_integer_t j;
    igraph_attribute_record_t *rec;
    igraph_vector_bool_t *log;
    igraph_bool_t l = igraph_i_cattribute_find(eal, name, &j);

    if (!l) {
        IGRAPH_ERRORF("The edge attribute '%s' does not exist.", IGRAPH_EINVAL, name);
    }

    rec = VECTOR(*eal)[j];
    if (rec->type != IGRAPH_ATTRIBUTE_BOOLEAN) {
        IGRAPH_ERRORF("Boolean edge attribute '%s' expected, got %s.", IGRAPH_EINVAL, name, attribute_type_name(rec->type));
    }
    log = (igraph_vector_bool_t*)rec->value;
    if (igraph_es_is_all(&es)) {
        igraph_vector_bool_clear(value);
        IGRAPH_CHECK(igraph_vector_bool_append(value, log));
    } else {
        igraph_eit_t it;
        igraph_integer_t i = 0;
        IGRAPH_CHECK(igraph_eit_create(graph, es, &it));
        IGRAPH_FINALLY(igraph_eit_destroy, &it);
        IGRAPH_CHECK(igraph_vector_bool_resize(value, IGRAPH_EIT_SIZE(it)));
        for (; !IGRAPH_EIT_END(it); IGRAPH_EIT_NEXT(it), i++) {
            igraph_integer_t e = IGRAPH_EIT_GET(it);
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
    &igraph_i_cattribute_gettype,
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
    igraph_vector_ptr_t *gal = &attr->gal;
    igraph_integer_t j;
    igraph_attribute_record_t *rec;
    igraph_vector_t *num;
    igraph_bool_t l = igraph_i_cattribute_find(gal, name, &j);

    if (!l) {
        IGRAPH_WARNINGF("Graph attribute '%s' does not exist, returning default numeric attribute value.", name);
        return IGRAPH_NAN;
    }

    rec = VECTOR(*gal)[j];
    num = (igraph_vector_t*)rec->value;
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
    igraph_vector_ptr_t *gal = &attr->gal;
    igraph_integer_t j;
    igraph_attribute_record_t *rec;
    igraph_vector_bool_t *log;
    igraph_bool_t l = igraph_i_cattribute_find(gal, name, &j);

    if (!l) {
        IGRAPH_WARNINGF("Graph attribute '%s' does not exist, returning default boolean attribute value.", name);
        return false;
    }

    rec = VECTOR(*gal)[j];
    log = (igraph_vector_bool_t*)rec->value;
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
    igraph_vector_ptr_t *gal = &attr->gal;
    igraph_integer_t j;
    igraph_attribute_record_t *rec;
    igraph_strvector_t *str;
    igraph_bool_t l = igraph_i_cattribute_find(gal, name, &j);

    if (!l) {
        IGRAPH_WARNINGF("Graph attribute '%s' does not exist, returning default string attribute value.", name);
        return "";
    }

    rec = VECTOR(*gal)[j];
    str = (igraph_strvector_t*)rec->value;
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
                                    igraph_integer_t vid) {
    igraph_i_cattributes_t *attr = graph->attr;
    igraph_vector_ptr_t *val = &attr->val;
    igraph_integer_t j;
    igraph_attribute_record_t *rec;
    igraph_vector_t *num;
    igraph_bool_t l = igraph_i_cattribute_find(val, name, &j);

    if (!l) {
        IGRAPH_WARNINGF("Vertex attribute '%s' does not exist, returning default numeric attribute value.", name);
        return IGRAPH_NAN;
    }

    rec = VECTOR(*val)[j];
    num = (igraph_vector_t*)rec->value;
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
                                    igraph_integer_t vid) {
    igraph_i_cattributes_t *attr = graph->attr;
    igraph_vector_ptr_t *val = &attr->val;
    igraph_integer_t j;
    igraph_attribute_record_t *rec;
    igraph_vector_bool_t *log;
    igraph_bool_t l = igraph_i_cattribute_find(val, name, &j);

    if (!l) {
        IGRAPH_WARNINGF("Vertex attribute '%s' does not exist, returning default boolean attribute value.", name);
        return false;
    }

    rec = VECTOR(*val)[j];
    log = (igraph_vector_bool_t*)rec->value;
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
                                  igraph_integer_t vid) {
    igraph_i_cattributes_t *attr = graph->attr;
    igraph_vector_ptr_t *val = &attr->val;
    igraph_integer_t j;
    igraph_attribute_record_t *rec;
    igraph_strvector_t *str;
    igraph_bool_t l = igraph_i_cattribute_find(val, name, &j);

    if (!l) {
        IGRAPH_WARNINGF("Vertex attribute '%s' does not exist, returning default string attribute value.", name);
        return "";
    }

    rec = VECTOR(*val)[j];
    str = (igraph_strvector_t*)rec->value;
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
                                    igraph_integer_t eid) {
    igraph_i_cattributes_t *attr = graph->attr;
    igraph_vector_ptr_t *eal = &attr->eal;
    igraph_integer_t j;
    igraph_attribute_record_t *rec;
    igraph_vector_t *num;
    igraph_bool_t l = igraph_i_cattribute_find(eal, name, &j);

    if (!l) {
        IGRAPH_WARNINGF("Edge attribute '%s' does not exist, returning default numeric attribute value.", name);
        return IGRAPH_NAN;
    }

    rec = VECTOR(*eal)[j];
    num = (igraph_vector_t*)rec->value;
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
                                    igraph_integer_t eid) {
    igraph_i_cattributes_t *attr = graph->attr;
    igraph_vector_ptr_t *eal = &attr->eal;
    igraph_integer_t j;
    igraph_attribute_record_t *rec;
    igraph_vector_bool_t *log;
    igraph_bool_t l = igraph_i_cattribute_find(eal, name, &j);

    if (!l) {
        IGRAPH_WARNINGF("Edge attribute '%s' does not exist, returning default boolean attribute value.", name);
        return false;
    }

    rec = VECTOR(*eal)[j];
    log = (igraph_vector_bool_t*)rec->value;
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
                                  igraph_integer_t eid) {
    igraph_i_cattributes_t *attr = graph->attr;
    igraph_vector_ptr_t *eal = &attr->eal;
    igraph_integer_t j;
    igraph_attribute_record_t *rec;
    igraph_strvector_t *str;
    igraph_bool_t l = igraph_i_cattribute_find(eal, name, &j);

    if (!l) {
        IGRAPH_WARNINGF("Edge attribute '%s' does not exist, returning default string attribute value.", name);
        return "";
    }

    rec = VECTOR(*eal)[j];
    str = (igraph_strvector_t*)rec->value;
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

    return igraph_i_cattribute_get_numeric_vertex_attr(graph, name, vids,
            result);
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

    return igraph_i_cattribute_get_bool_vertex_attr(graph, name, vids,
            result);
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

    return igraph_i_cattribute_get_numeric_edge_attr(graph, name, eids,
            result);
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

    return igraph_i_cattribute_get_bool_edge_attr(graph, name, eids,
            result);
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

    return igraph_i_cattribute_get_string_vertex_attr(graph, name, vids,
            result);
}

/**
 * \function igraph_cattribute_EASV
 * \brief Query a string edge attribute for many edges.
 *
 * \param graph The input graph.
 * \param name The name of the attribute.
 * \param vids The edges to query.
 * \param result Pointer to an initialized string vector, the result
 *     is stored here. It will be resized, if needed.
 * \return Error code.
 *
 * Time complexity: O(e), where e is the number of edges in
 * 'eids'. (We assume that the string attributes have a bounded length.)
 */

igraph_error_t igraph_cattribute_EASV(const igraph_t *graph, const char *name,
                           igraph_es_t eids, igraph_strvector_t *result) {

    return igraph_i_cattribute_get_string_edge_attr(graph, name, eids,
            result);
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
 * \return Logical value, \c true if the attribute exists, \c false otherwise.
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
    igraph_vector_ptr_t *gal = &attr->gal;
    igraph_integer_t j;
    igraph_bool_t l = igraph_i_cattribute_find(gal, name, &j);

    if (l) {
        igraph_attribute_record_t *rec = VECTOR(*gal)[j];
        if (rec->type != IGRAPH_ATTRIBUTE_NUMERIC) {
            IGRAPH_ERROR("Invalid attribute type", IGRAPH_EINVAL);
        } else {
            igraph_vector_t *num = (igraph_vector_t *)rec->value;
            VECTOR(*num)[0] = value;
        }
    } else {
        igraph_attribute_record_t *rec = IGRAPH_CALLOC(1, igraph_attribute_record_t);
        igraph_vector_t *num;
        if (!rec) {
            IGRAPH_ERROR("Cannot add graph attribute", IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */
        }
        IGRAPH_FINALLY(igraph_free, rec);
        rec->name = strdup(name);
        if (!rec->name) {
            IGRAPH_ERROR("Cannot add graph attribute", IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */
        }
        IGRAPH_FINALLY(igraph_free, (char*)rec->name);
        rec->type = IGRAPH_ATTRIBUTE_NUMERIC;
        num = IGRAPH_CALLOC(1, igraph_vector_t);
        if (!num) {
            IGRAPH_ERROR("Cannot add graph attribute", IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */
        }
        IGRAPH_FINALLY(igraph_free, num);
        IGRAPH_VECTOR_INIT_FINALLY(num, 1);
        VECTOR(*num)[0] = value;
        rec->value = num;
        IGRAPH_CHECK(igraph_vector_ptr_push_back(gal, rec));
        IGRAPH_FINALLY_CLEAN(4);
    }

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
    igraph_vector_ptr_t *gal = &attr->gal;
    igraph_integer_t j;
    igraph_bool_t l = igraph_i_cattribute_find(gal, name, &j);

    if (l) {
        igraph_attribute_record_t *rec = VECTOR(*gal)[j];
        if (rec->type != IGRAPH_ATTRIBUTE_BOOLEAN) {
            IGRAPH_ERROR("Invalid attribute type", IGRAPH_EINVAL);
        } else {
            igraph_vector_bool_t *log = (igraph_vector_bool_t *)rec->value;
            VECTOR(*log)[0] = value;
        }
    } else {
        igraph_attribute_record_t *rec = IGRAPH_CALLOC(1, igraph_attribute_record_t);
        igraph_vector_bool_t *log;
        if (!rec) {
            IGRAPH_ERROR("Cannot add graph attribute", IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */
        }
        IGRAPH_FINALLY(igraph_free, rec);
        rec->name = strdup(name);
        if (!rec->name) {
            IGRAPH_ERROR("Cannot add graph attribute", IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */
        }
        IGRAPH_FINALLY(igraph_free, (char*)rec->name);
        rec->type = IGRAPH_ATTRIBUTE_BOOLEAN;
        log = IGRAPH_CALLOC(1, igraph_vector_bool_t);
        if (!log) {
            IGRAPH_ERROR("Cannot add graph attribute", IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */
        }
        IGRAPH_FINALLY(igraph_free, log);
        IGRAPH_VECTOR_BOOL_INIT_FINALLY(log, 1);
        VECTOR(*log)[0] = value;
        rec->value = log;
        IGRAPH_CHECK(igraph_vector_ptr_push_back(gal, rec));
        IGRAPH_FINALLY_CLEAN(4);
    }

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
    igraph_vector_ptr_t *gal = &attr->gal;
    igraph_integer_t j;
    igraph_bool_t l = igraph_i_cattribute_find(gal, name, &j);

    if (l) {
        igraph_attribute_record_t *rec = VECTOR(*gal)[j];
        if (rec->type != IGRAPH_ATTRIBUTE_STRING) {
            IGRAPH_ERROR("Invalid attribute type", IGRAPH_EINVAL);
        } else {
            igraph_strvector_t *str = (igraph_strvector_t*)rec->value;
            IGRAPH_CHECK(igraph_strvector_set(str, 0, value));
        }
    } else {
        igraph_attribute_record_t *rec = IGRAPH_CALLOC(1, igraph_attribute_record_t);
        igraph_strvector_t *str;
        if (!rec) {
            IGRAPH_ERROR("Cannot add graph attribute", IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */
        }
        IGRAPH_FINALLY(igraph_free, rec);
        rec->name = strdup(name);
        if (!rec->name) {
            IGRAPH_ERROR("Cannot add graph attribute", IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */
        }
        IGRAPH_FINALLY(igraph_free, (char*)rec->name);
        rec->type = IGRAPH_ATTRIBUTE_STRING;
        str = IGRAPH_CALLOC(1, igraph_strvector_t);
        if (!str) {
            IGRAPH_ERROR("Cannot add graph attribute", IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */
        }
        IGRAPH_FINALLY(igraph_free, str);
        IGRAPH_STRVECTOR_INIT_FINALLY(str, 1);
        IGRAPH_CHECK(igraph_strvector_set(str, 0, value));
        rec->value = str;
        IGRAPH_CHECK(igraph_vector_ptr_push_back(gal, rec));
        IGRAPH_FINALLY_CLEAN(4);
    }

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
                              igraph_integer_t vid, igraph_real_t value) {

    igraph_i_cattributes_t *attr = graph->attr;
    igraph_vector_ptr_t *val = &attr->val;
    igraph_integer_t j;
    igraph_bool_t l = igraph_i_cattribute_find(val, name, &j);

    if (l) {
        igraph_attribute_record_t *rec = VECTOR(*val)[j];
        if (rec->type != IGRAPH_ATTRIBUTE_NUMERIC) {
            IGRAPH_ERROR("Invalid attribute type", IGRAPH_EINVAL);
        } else {
            igraph_vector_t *num = (igraph_vector_t*)rec->value;
            VECTOR(*num)[vid] = value;
        }
    } else {
        igraph_attribute_record_t *rec = IGRAPH_CALLOC(1, igraph_attribute_record_t);
        igraph_vector_t *num;
        if (!rec) {
            IGRAPH_ERROR("Cannot add vertex attribute", IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */
        }
        IGRAPH_FINALLY(igraph_free, rec);
        rec->name = strdup(name);
        if (!rec->name) {
            IGRAPH_ERROR("Cannot add vertex attribute", IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */
        }
        IGRAPH_FINALLY(igraph_free, (char*)rec->name);
        rec->type = IGRAPH_ATTRIBUTE_NUMERIC;
        num = IGRAPH_CALLOC(1, igraph_vector_t);
        if (!num) {
            IGRAPH_ERROR("Cannot add vertex attribute", IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */
        }
        IGRAPH_FINALLY(igraph_free, num);
        IGRAPH_VECTOR_INIT_FINALLY(num, igraph_vcount(graph));
        igraph_vector_fill(num, IGRAPH_NAN);
        VECTOR(*num)[vid] = value;
        rec->value = num;
        IGRAPH_CHECK(igraph_vector_ptr_push_back(val, rec));
        IGRAPH_FINALLY_CLEAN(4);
    }

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
                              igraph_integer_t vid, igraph_bool_t value) {

    igraph_i_cattributes_t *attr = graph->attr;
    igraph_vector_ptr_t *val = &attr->val;
    igraph_integer_t j;
    igraph_bool_t l = igraph_i_cattribute_find(val, name, &j);

    if (l) {
        igraph_attribute_record_t *rec = VECTOR(*val)[j];
        if (rec->type != IGRAPH_ATTRIBUTE_BOOLEAN) {
            IGRAPH_ERROR("Invalid attribute type", IGRAPH_EINVAL);
        } else {
            igraph_vector_bool_t *log = (igraph_vector_bool_t*)rec->value;
            VECTOR(*log)[vid] = value;
        }
    } else {
        igraph_attribute_record_t *rec = IGRAPH_CALLOC(1, igraph_attribute_record_t);
        igraph_vector_bool_t *log;
        if (!rec) {
            IGRAPH_ERROR("Cannot add vertex attribute", IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */
        }
        IGRAPH_FINALLY(igraph_free, rec);
        rec->name = strdup(name);
        if (!rec->name) {
            IGRAPH_ERROR("Cannot add vertex attribute", IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */
        }
        IGRAPH_FINALLY(igraph_free, (char*)rec->name);
        rec->type = IGRAPH_ATTRIBUTE_BOOLEAN;
        log = IGRAPH_CALLOC(1, igraph_vector_bool_t);
        if (!log) {
            IGRAPH_ERROR("Cannot add vertex attribute", IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */
        }
        IGRAPH_FINALLY(igraph_free, log);
        IGRAPH_VECTOR_BOOL_INIT_FINALLY(log, igraph_vcount(graph));
        igraph_vector_bool_fill(log, false);
        VECTOR(*log)[vid] = value;
        rec->value = log;
        IGRAPH_CHECK(igraph_vector_ptr_push_back(val, rec));
        IGRAPH_FINALLY_CLEAN(4);
    }

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
                              igraph_integer_t vid, const char *value) {

    igraph_i_cattributes_t *attr = graph->attr;
    igraph_vector_ptr_t *val = &attr->val;
    igraph_integer_t j;
    igraph_bool_t l = igraph_i_cattribute_find(val, name, &j);

    if (l) {
        igraph_attribute_record_t *rec = VECTOR(*val)[j];
        if (rec->type != IGRAPH_ATTRIBUTE_STRING) {
            IGRAPH_ERROR("Invalid attribute type", IGRAPH_EINVAL);
        } else {
            igraph_strvector_t *str = (igraph_strvector_t*)rec->value;
            IGRAPH_CHECK(igraph_strvector_set(str, vid, value));
        }
    } else {
        igraph_attribute_record_t *rec = IGRAPH_CALLOC(1, igraph_attribute_record_t);
        igraph_strvector_t *str;
        if (!rec) {
            IGRAPH_ERROR("Cannot add vertex attribute", IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */
        }
        IGRAPH_FINALLY(igraph_free, rec);
        rec->name = strdup(name);
        if (!rec->name) {
            IGRAPH_ERROR("Cannot add vertex attribute", IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */
        }
        IGRAPH_FINALLY(igraph_free, (char*)rec->name);
        rec->type = IGRAPH_ATTRIBUTE_STRING;
        str = IGRAPH_CALLOC(1, igraph_strvector_t);
        if (!str) {
            IGRAPH_ERROR("Cannot add vertex attribute", IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */
        }
        IGRAPH_FINALLY(igraph_free, str);
        IGRAPH_STRVECTOR_INIT_FINALLY(str, igraph_vcount(graph));
        IGRAPH_CHECK(igraph_strvector_set(str, vid, value));
        rec->value = str;
        IGRAPH_CHECK(igraph_vector_ptr_push_back(val, rec));
        IGRAPH_FINALLY_CLEAN(4);
    }

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
                              igraph_integer_t eid, igraph_real_t value) {

    igraph_i_cattributes_t *attr = graph->attr;
    igraph_vector_ptr_t *eal = &attr->eal;
    igraph_integer_t j;
    igraph_bool_t l = igraph_i_cattribute_find(eal, name, &j);

    if (l) {
        igraph_attribute_record_t *rec = VECTOR(*eal)[j];
        if (rec->type != IGRAPH_ATTRIBUTE_NUMERIC) {
            IGRAPH_ERROR("Invalid attribute type", IGRAPH_EINVAL);
        } else {
            igraph_vector_t *num = (igraph_vector_t*)rec->value;
            VECTOR(*num)[eid] = value;
        }
    } else {
        igraph_attribute_record_t *rec = IGRAPH_CALLOC(1, igraph_attribute_record_t);
        igraph_vector_t *num;
        if (!rec) {
            IGRAPH_ERROR("Cannot add edge attribute", IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */
        }
        IGRAPH_FINALLY(igraph_free, rec);
        rec->name = strdup(name);
        if (!rec->name) {
            IGRAPH_ERROR("Cannot add edge attribute", IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */
        }
        IGRAPH_FINALLY(igraph_free, (char*)rec->name);
        rec->type = IGRAPH_ATTRIBUTE_NUMERIC;
        num = IGRAPH_CALLOC(1, igraph_vector_t);
        if (!num) {
            IGRAPH_ERROR("Cannot add edge attribute", IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */
        }
        IGRAPH_FINALLY(igraph_free, num);
        IGRAPH_VECTOR_INIT_FINALLY(num, igraph_ecount(graph));
        igraph_vector_fill(num, IGRAPH_NAN);
        VECTOR(*num)[eid] = value;
        rec->value = num;
        IGRAPH_CHECK(igraph_vector_ptr_push_back(eal, rec));
        IGRAPH_FINALLY_CLEAN(4);
    }

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
                              igraph_integer_t eid, igraph_bool_t value) {

    igraph_i_cattributes_t *attr = graph->attr;
    igraph_vector_ptr_t *eal = &attr->eal;
    igraph_integer_t j;
    igraph_bool_t l = igraph_i_cattribute_find(eal, name, &j);

    if (l) {
        igraph_attribute_record_t *rec = VECTOR(*eal)[j];
        if (rec->type != IGRAPH_ATTRIBUTE_BOOLEAN) {
            IGRAPH_ERROR("Invalid attribute type", IGRAPH_EINVAL);
        } else {
            igraph_vector_bool_t *log = (igraph_vector_bool_t*)rec->value;
            VECTOR(*log)[eid] = value;
        }
    } else {
        igraph_attribute_record_t *rec = IGRAPH_CALLOC(1, igraph_attribute_record_t);
        igraph_vector_bool_t *log;
        if (!rec) {
            IGRAPH_ERROR("Cannot add edge attribute", IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */
        }
        IGRAPH_FINALLY(igraph_free, rec);
        rec->name = strdup(name);
        if (!rec->name) {
            IGRAPH_ERROR("Cannot add edge attribute", IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */
        }
        IGRAPH_FINALLY(igraph_free, (char*)rec->name);
        rec->type = IGRAPH_ATTRIBUTE_BOOLEAN;
        log = IGRAPH_CALLOC(1, igraph_vector_bool_t);
        if (!log) {
            IGRAPH_ERROR("Cannot add edge attribute", IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */
        }
        IGRAPH_FINALLY(igraph_free, log);
        IGRAPH_VECTOR_BOOL_INIT_FINALLY(log, igraph_ecount(graph));
        igraph_vector_bool_fill(log, false);
        VECTOR(*log)[eid] = value;
        rec->value = log;
        IGRAPH_CHECK(igraph_vector_ptr_push_back(eal, rec));
        IGRAPH_FINALLY_CLEAN(4);
    }

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
                              igraph_integer_t eid, const char *value) {

    igraph_i_cattributes_t *attr = graph->attr;
    igraph_vector_ptr_t *eal = &attr->eal;
    igraph_integer_t j;
    igraph_bool_t l = igraph_i_cattribute_find(eal, name, &j);

    if (l) {
        igraph_attribute_record_t *rec = VECTOR(*eal)[j];
        if (rec->type != IGRAPH_ATTRIBUTE_STRING) {
            IGRAPH_ERROR("Invalid attribute type", IGRAPH_EINVAL);
        } else {
            igraph_strvector_t *str = (igraph_strvector_t*)rec->value;
            IGRAPH_CHECK(igraph_strvector_set(str, eid, value));
        }
    } else {
        igraph_attribute_record_t *rec = IGRAPH_CALLOC(1, igraph_attribute_record_t);
        igraph_strvector_t *str;
        if (!rec) {
            IGRAPH_ERROR("Cannot add edge attribute", IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */
        }
        IGRAPH_FINALLY(igraph_free, rec);
        rec->name = strdup(name);
        if (!rec->name) {
            IGRAPH_ERROR("Cannot add edge attribute", IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */
        }
        IGRAPH_FINALLY(igraph_free, (char*)rec->name);
        rec->type = IGRAPH_ATTRIBUTE_STRING;
        str = IGRAPH_CALLOC(1, igraph_strvector_t);
        if (!str) {
            IGRAPH_ERROR("Cannot add edge attribute", IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */
        }
        IGRAPH_FINALLY(igraph_free, str);
        IGRAPH_STRVECTOR_INIT_FINALLY(str, igraph_ecount(graph));
        IGRAPH_CHECK(igraph_strvector_set(str, eid, value));
        rec->value = str;
        IGRAPH_CHECK(igraph_vector_ptr_push_back(eal, rec));
        IGRAPH_FINALLY_CLEAN(4);
    }

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
    igraph_vector_ptr_t *val = &attr->val;
    igraph_integer_t j;
    igraph_bool_t l = igraph_i_cattribute_find(val, name, &j);

    /* Check length first */
    if (igraph_vector_size(v) != igraph_vcount(graph)) {
        IGRAPH_ERROR("Invalid vertex attribute vector length", IGRAPH_EINVAL);
    }

    if (l) {
        /* Already present, check type */
        igraph_attribute_record_t *rec = VECTOR(*val)[j];
        igraph_vector_t *num = (igraph_vector_t *)rec->value;
        if (rec->type != IGRAPH_ATTRIBUTE_NUMERIC) {
            IGRAPH_ERROR("Attribute type mismatch", IGRAPH_EINVAL);
        }
        igraph_vector_clear(num);
        IGRAPH_CHECK(igraph_vector_append(num, v));
    } else {
        /* Add it */
        igraph_attribute_record_t *rec = IGRAPH_CALLOC(1, igraph_attribute_record_t);
        igraph_vector_t *num;
        if (!rec) {
            IGRAPH_ERROR("Cannot add vertex attribute", IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */
        }
        IGRAPH_FINALLY(igraph_free, rec);
        rec->type = IGRAPH_ATTRIBUTE_NUMERIC;
        rec->name = strdup(name);
        if (!rec->name) {
            IGRAPH_ERROR("Cannot add vertex attribute", IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */
        }
        IGRAPH_FINALLY(igraph_free, (char*)rec->name);
        num = IGRAPH_CALLOC(1, igraph_vector_t);
        if (!num) {
            IGRAPH_ERROR("Cannot add vertex attribute", IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */
        }
        IGRAPH_FINALLY(igraph_free, num);
        rec->value = num;
        IGRAPH_CHECK(igraph_vector_init_copy(num, v));
        IGRAPH_FINALLY(igraph_vector_destroy, num);
        IGRAPH_CHECK(igraph_vector_ptr_push_back(val, rec));
        IGRAPH_FINALLY_CLEAN(4);
    }

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
    igraph_vector_ptr_t *val = &attr->val;
    igraph_integer_t j;
    igraph_bool_t l = igraph_i_cattribute_find(val, name, &j);

    /* Check length first */
    if (igraph_vector_bool_size(v) != igraph_vcount(graph)) {
        IGRAPH_ERROR("Invalid vertex attribute vector length", IGRAPH_EINVAL);
    }

    if (l) {
        /* Already present, check type */
        igraph_attribute_record_t *rec = VECTOR(*val)[j];
        igraph_vector_bool_t *log = (igraph_vector_bool_t *)rec->value;
        if (rec->type != IGRAPH_ATTRIBUTE_BOOLEAN) {
            IGRAPH_ERROR("Attribute type mismatch", IGRAPH_EINVAL);
        }
        igraph_vector_bool_clear(log);
        IGRAPH_CHECK(igraph_vector_bool_append(log, v));
    } else {
        /* Add it */
        igraph_attribute_record_t *rec = IGRAPH_CALLOC(1, igraph_attribute_record_t);
        igraph_vector_bool_t *log;
        if (!rec) {
            IGRAPH_ERROR("Cannot add vertex attribute", IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */
        }
        IGRAPH_FINALLY(igraph_free, rec);
        rec->type = IGRAPH_ATTRIBUTE_BOOLEAN;
        rec->name = strdup(name);
        if (!rec->name) {
            IGRAPH_ERROR("Cannot add vertex attribute", IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */
        }
        IGRAPH_FINALLY(igraph_free, (char*)rec->name);
        log = IGRAPH_CALLOC(1, igraph_vector_bool_t);
        if (!log) {
            IGRAPH_ERROR("Cannot add vertex attribute", IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */
        }
        IGRAPH_FINALLY(igraph_free, log);
        rec->value = log;
        IGRAPH_CHECK(igraph_vector_bool_init_copy(log, v));
        IGRAPH_FINALLY(igraph_vector_bool_destroy, log);
        IGRAPH_CHECK(igraph_vector_ptr_push_back(val, rec));
        IGRAPH_FINALLY_CLEAN(4);
    }

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
    igraph_vector_ptr_t *val = &attr->val;
    igraph_integer_t j;
    igraph_bool_t l = igraph_i_cattribute_find(val, name, &j);

    /* Check length first */
    if (igraph_strvector_size(sv) != igraph_vcount(graph)) {
        IGRAPH_ERROR("Invalid vertex attribute vector length", IGRAPH_EINVAL);
    }

    if (l) {
        /* Already present, check type */
        igraph_attribute_record_t *rec = VECTOR(*val)[j];
        igraph_strvector_t *str = (igraph_strvector_t *)rec->value;
        if (rec->type != IGRAPH_ATTRIBUTE_STRING) {
            IGRAPH_ERROR("Attribute type mismatch", IGRAPH_EINVAL);
        }
        igraph_strvector_clear(str);
        IGRAPH_CHECK(igraph_strvector_append(str, sv));
    } else {
        /* Add it */
        igraph_attribute_record_t *rec = IGRAPH_CALLOC(1, igraph_attribute_record_t);
        igraph_strvector_t *str;
        if (!rec) {
            IGRAPH_ERROR("Cannot add vertex attribute", IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */
        }
        IGRAPH_FINALLY(igraph_free, rec);
        rec->type = IGRAPH_ATTRIBUTE_STRING;
        rec->name = strdup(name);
        if (!rec->name) {
            IGRAPH_ERROR("Cannot add vertex attribute", IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */
        }
        IGRAPH_FINALLY(igraph_free, (char*)rec->name);
        str = IGRAPH_CALLOC(1, igraph_strvector_t);
        if (!str) {
            IGRAPH_ERROR("Cannot add vertex attribute", IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */
        }
        IGRAPH_FINALLY(igraph_free, str);
        rec->value = str;
        IGRAPH_CHECK(igraph_strvector_init_copy(str, sv));
        IGRAPH_FINALLY(igraph_strvector_destroy, str);
        IGRAPH_CHECK(igraph_vector_ptr_push_back(val, rec));
        IGRAPH_FINALLY_CLEAN(4);
    }

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
    igraph_vector_ptr_t *eal = &attr->eal;
    igraph_integer_t j;
    igraph_bool_t l = igraph_i_cattribute_find(eal, name, &j);

    /* Check length first */
    if (igraph_vector_size(v) != igraph_ecount(graph)) {
        IGRAPH_ERROR("Invalid edge attribute vector length", IGRAPH_EINVAL);
    }

    if (l) {
        /* Already present, check type */
        igraph_attribute_record_t *rec = VECTOR(*eal)[j];
        igraph_vector_t *num = (igraph_vector_t *)rec->value;
        if (rec->type != IGRAPH_ATTRIBUTE_NUMERIC) {
            IGRAPH_ERROR("Attribute type mismatch", IGRAPH_EINVAL);
        }
        igraph_vector_clear(num);
        IGRAPH_CHECK(igraph_vector_append(num, v));
    } else {
        /* Add it */
        igraph_attribute_record_t *rec = IGRAPH_CALLOC(1, igraph_attribute_record_t);
        igraph_vector_t *num;
        if (!rec) {
            IGRAPH_ERROR("Cannot add edge attribute", IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */
        }
        IGRAPH_FINALLY(igraph_free, rec);
        rec->type = IGRAPH_ATTRIBUTE_NUMERIC;
        rec->name = strdup(name);
        if (!rec->name) {
            IGRAPH_ERROR("Cannot add edge attribute", IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */
        }
        IGRAPH_FINALLY(igraph_free, (char*)rec->name);
        num = IGRAPH_CALLOC(1, igraph_vector_t);
        if (!num) {
            IGRAPH_ERROR("Cannot add edge attribute", IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */
        }
        IGRAPH_FINALLY(igraph_free, num);
        rec->value = num;
        IGRAPH_CHECK(igraph_vector_init_copy(num, v));
        IGRAPH_FINALLY(igraph_vector_destroy, num);
        IGRAPH_CHECK(igraph_vector_ptr_push_back(eal, rec));
        IGRAPH_FINALLY_CLEAN(4);
    }

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
    igraph_vector_ptr_t *eal = &attr->eal;
    igraph_integer_t j;
    igraph_bool_t l = igraph_i_cattribute_find(eal, name, &j);

    /* Check length first */
    if (igraph_vector_bool_size(v) != igraph_ecount(graph)) {
        IGRAPH_ERROR("Invalid edge attribute vector length", IGRAPH_EINVAL);
    }

    if (l) {
        /* Already present, check type */
        igraph_attribute_record_t *rec = VECTOR(*eal)[j];
        igraph_vector_bool_t *log = (igraph_vector_bool_t *)rec->value;
        if (rec->type != IGRAPH_ATTRIBUTE_BOOLEAN) {
            IGRAPH_ERROR("Attribute type mismatch", IGRAPH_EINVAL);
        }
        igraph_vector_bool_clear(log);
        IGRAPH_CHECK(igraph_vector_bool_append(log, v));
    } else {
        /* Add it */
        igraph_attribute_record_t *rec = IGRAPH_CALLOC(1, igraph_attribute_record_t);
        igraph_vector_bool_t *log;
        if (!rec) {
            IGRAPH_ERROR("Cannot add edge attribute", IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */
        }
        IGRAPH_FINALLY(igraph_free, rec);
        rec->type = IGRAPH_ATTRIBUTE_BOOLEAN;
        rec->name = strdup(name);
        if (!rec->name) {
            IGRAPH_ERROR("Cannot add edge attribute", IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */
        }
        IGRAPH_FINALLY(igraph_free, (char*)rec->name);
        log = IGRAPH_CALLOC(1, igraph_vector_bool_t);
        if (!log) {
            IGRAPH_ERROR("Cannot add edge attribute", IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */
        }
        IGRAPH_FINALLY(igraph_free, log);
        rec->value = log;
        IGRAPH_CHECK(igraph_vector_bool_init_copy(log, v));
        IGRAPH_FINALLY(igraph_vector_bool_destroy, log);
        IGRAPH_CHECK(igraph_vector_ptr_push_back(eal, rec));
        IGRAPH_FINALLY_CLEAN(4);
    }

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
    igraph_vector_ptr_t *eal = &attr->eal;
    igraph_integer_t j;
    igraph_bool_t l = igraph_i_cattribute_find(eal, name, &j);

    /* Check length first */
    if (igraph_strvector_size(sv) != igraph_ecount(graph)) {
        IGRAPH_ERROR("Invalid edge attribute vector length", IGRAPH_EINVAL);
    }

    if (l) {
        /* Already present, check type */
        igraph_attribute_record_t *rec = VECTOR(*eal)[j];
        igraph_strvector_t *str = (igraph_strvector_t *)rec->value;
        if (rec->type != IGRAPH_ATTRIBUTE_STRING) {
            IGRAPH_ERROR("Attribute type mismatch", IGRAPH_EINVAL);
        }
        igraph_strvector_clear(str);
        IGRAPH_CHECK(igraph_strvector_append(str, sv));
    } else {
        /* Add it */
        igraph_attribute_record_t *rec = IGRAPH_CALLOC(1, igraph_attribute_record_t);
        igraph_strvector_t *str;
        if (!rec) {
            IGRAPH_ERROR("Cannot add vertex attribute", IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */
        }
        IGRAPH_FINALLY(igraph_free, rec);
        rec->type = IGRAPH_ATTRIBUTE_STRING;
        rec->name = strdup(name);
        if (!rec->name) {
            IGRAPH_ERROR("Cannot add vertex attribute", IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */
        }
        IGRAPH_FINALLY(igraph_free, (char*)rec->name);
        str = IGRAPH_CALLOC(1, igraph_strvector_t);
        if (!str) {
            IGRAPH_ERROR("Cannot add vertex attribute", IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */
        }
        IGRAPH_FINALLY(igraph_free, str);
        rec->value = str;
        IGRAPH_CHECK(igraph_strvector_init_copy(str, sv));
        IGRAPH_FINALLY(igraph_strvector_destroy, str);
        IGRAPH_CHECK(igraph_vector_ptr_push_back(eal, rec));
        IGRAPH_FINALLY_CLEAN(4);
    }

    return IGRAPH_SUCCESS;
}

static void igraph_i_cattribute_free_rec(igraph_attribute_record_t *rec) {

    if (rec->type == IGRAPH_ATTRIBUTE_NUMERIC) {
        igraph_vector_t *num = (igraph_vector_t*)rec->value;
        igraph_vector_destroy(num);
    } else if (rec->type == IGRAPH_ATTRIBUTE_STRING) {
        igraph_strvector_t *str = (igraph_strvector_t*)rec->value;
        igraph_strvector_destroy(str);
    } else if (rec->type == IGRAPH_ATTRIBUTE_BOOLEAN) {
        igraph_vector_bool_t *boolvec = (igraph_vector_bool_t*)rec->value;
        igraph_vector_bool_destroy(boolvec);
    }
    IGRAPH_FREE(rec->name);
    IGRAPH_FREE(rec->value);
    IGRAPH_FREE(rec);
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
    igraph_vector_ptr_t *gal = &attr->gal;
    igraph_integer_t j;
    igraph_bool_t l = igraph_i_cattribute_find(gal, name, &j);

    if (l) {
        igraph_i_cattribute_free_rec(VECTOR(*gal)[j]);
        igraph_vector_ptr_remove(gal, j);
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
    igraph_vector_ptr_t *val = &attr->val;
    igraph_integer_t j;
    igraph_bool_t l = igraph_i_cattribute_find(val, name, &j);

    if (l) {
        igraph_i_cattribute_free_rec(VECTOR(*val)[j]);
        igraph_vector_ptr_remove(val, j);
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
 *
 */
void igraph_cattribute_remove_e(igraph_t *graph, const char *name) {

    igraph_i_cattributes_t *attr = graph->attr;
    igraph_vector_ptr_t *eal = &attr->eal;
    igraph_integer_t j;
    igraph_bool_t l = igraph_i_cattribute_find(eal, name, &j);

    if (l) {
        igraph_i_cattribute_free_rec(VECTOR(*eal)[j]);
        igraph_vector_ptr_remove(eal, j);
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
        igraph_vector_ptr_t *gal = &attr->gal;
        igraph_integer_t i, n = igraph_vector_ptr_size(gal);
        for (i = 0; i < n; i++) {
            igraph_i_cattribute_free_rec(VECTOR(*gal)[i]);
        }
        igraph_vector_ptr_clear(gal);
    }
    if (v) {
        igraph_vector_ptr_t *val = &attr->val;
        igraph_integer_t i, n = igraph_vector_ptr_size(val);
        for (i = 0; i < n; i++) {
            igraph_i_cattribute_free_rec(VECTOR(*val)[i]);
        }
        igraph_vector_ptr_clear(val);
    }
    if (e) {
        igraph_vector_ptr_t *eal = &attr->eal;
        igraph_integer_t i, n = igraph_vector_ptr_size(eal);
        for (i = 0; i < n; i++) {
            igraph_i_cattribute_free_rec(VECTOR(*eal)[i]);
        }
        igraph_vector_ptr_clear(eal);
    }
}
