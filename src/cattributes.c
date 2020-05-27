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
#include "igraph_math.h"
#include "igraph_interface.h"
#include "igraph_random.h"
#include "config.h"

#include <string.h>

/* An attribute is either a numeric vector (vector_t) or a string
   vector (strvector_t). The attribute itself is stored in a
   struct igraph_attribute_record_t, there is one such object for each
   attribute. The igraph_t has a pointer to an array of three
   vector_ptr_t's which contains pointers to
   igraph_i_cattribute_t's. Graph attributes are first, then vertex
   and edge attributes. */

igraph_bool_t igraph_i_cattribute_find(const igraph_vector_ptr_t *ptrvec,
                                       const char *name, long int *idx) {
    long int i, n = igraph_vector_ptr_size(ptrvec);
    igraph_bool_t l = 0;
    for (i = 0; !l && i < n; i++) {
        igraph_attribute_record_t *rec = VECTOR(*ptrvec)[i];
        l = !strcmp(rec->name, name);
    }
    if (idx) {
        *idx = i - 1;
    }
    return l;
}

typedef struct igraph_i_cattributes_t {
    igraph_vector_ptr_t gal;
    igraph_vector_ptr_t val;
    igraph_vector_ptr_t eal;
} igraph_i_cattributes_t;

int igraph_i_cattributes_copy_attribute_record(igraph_attribute_record_t **newrec,
        const igraph_attribute_record_t *rec) {
    igraph_vector_t *num, *newnum;
    igraph_strvector_t *str, *newstr;

    *newrec = igraph_Calloc(1, igraph_attribute_record_t);
    if (!(*newrec)) {
        IGRAPH_ERROR("Cannot copy attributes", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, *newrec);
    (*newrec)->type = rec->type;
    (*newrec)->name = strdup(rec->name);
    if (!(*newrec)->name) {
        IGRAPH_ERROR("Cannot copy attributes", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, (void*)(*newrec)->name);
    if (rec->type == IGRAPH_ATTRIBUTE_NUMERIC) {
        num = (igraph_vector_t *)rec->value;
        newnum = igraph_Calloc(1, igraph_vector_t);
        if (!newnum) {
            IGRAPH_ERROR("Cannot copy attributes", IGRAPH_ENOMEM);
        }
        IGRAPH_FINALLY(igraph_free, newnum);
        IGRAPH_CHECK(igraph_vector_copy(newnum, num));
        IGRAPH_FINALLY(igraph_vector_destroy, newnum);
        (*newrec)->value = newnum;
    } else if (rec->type == IGRAPH_ATTRIBUTE_STRING) {
        str = (igraph_strvector_t*)rec->value;
        newstr = igraph_Calloc(1, igraph_strvector_t);
        if (!newstr) {
            IGRAPH_ERROR("Cannot copy attributes", IGRAPH_ENOMEM);
        }
        IGRAPH_FINALLY(igraph_free, newstr);
        IGRAPH_CHECK(igraph_strvector_copy(newstr, str));
        IGRAPH_FINALLY(igraph_strvector_destroy, newstr);
        (*newrec)->value = newstr;
    } else if (rec->type == IGRAPH_ATTRIBUTE_BOOLEAN) {
        igraph_vector_bool_t *log = (igraph_vector_bool_t*) rec->value;
        igraph_vector_bool_t *newlog = igraph_Calloc(1, igraph_vector_bool_t);
        if (!newlog) {
            IGRAPH_ERROR("Cannot copy attributes", IGRAPH_ENOMEM);
        }
        IGRAPH_FINALLY(igraph_free, newlog);
        IGRAPH_CHECK(igraph_vector_bool_copy(newlog, log));
        IGRAPH_FINALLY(igraph_vector_bool_destroy, newlog);
        (*newrec)->value = newlog;
    }

    IGRAPH_FINALLY_CLEAN(4);
    return 0;
}


int igraph_i_cattribute_init(igraph_t *graph, igraph_vector_ptr_t *attr) {
    igraph_attribute_record_t *attr_rec;
    long int i, n;
    igraph_i_cattributes_t *nattr;

    n = attr ? igraph_vector_ptr_size(attr) : 0;

    nattr = igraph_Calloc(1, igraph_i_cattributes_t);
    if (!nattr) {
        IGRAPH_ERROR("Can't init attributes", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, nattr);

    IGRAPH_CHECK(igraph_vector_ptr_init(&nattr->gal, n));
    IGRAPH_FINALLY(igraph_vector_ptr_destroy, &nattr->gal);
    IGRAPH_CHECK(igraph_vector_ptr_init(&nattr->val, 0));
    IGRAPH_FINALLY(igraph_vector_ptr_destroy, &nattr->val);
    IGRAPH_CHECK(igraph_vector_ptr_init(&nattr->eal, 0));
    IGRAPH_FINALLY_CLEAN(3);

    for (i = 0; i < n; i++) {
        IGRAPH_CHECK(igraph_i_cattributes_copy_attribute_record(
                         &attr_rec, VECTOR(*attr)[i]));
        VECTOR(nattr->gal)[i] = attr_rec;
    }

    graph->attr = nattr;

    return 0;
}

void igraph_i_cattribute_destroy(igraph_t *graph) {
    igraph_i_cattributes_t *attr = graph->attr;
    igraph_vector_ptr_t *als[3] = { &attr->gal, &attr->val, &attr->eal };
    long int i, n, a;
    igraph_vector_t *num;
    igraph_strvector_t *str;
    igraph_vector_bool_t *boolvec;
    igraph_attribute_record_t *rec;
    for (a = 0; a < 3; a++) {
        n = igraph_vector_ptr_size(als[a]);
        for (i = 0; i < n; i++) {
            rec = VECTOR(*als[a])[i];
            if (rec) {
                if (rec->type == IGRAPH_ATTRIBUTE_NUMERIC) {
                    num = (igraph_vector_t*)rec->value;
                    igraph_vector_destroy(num);
                    igraph_free(num);
                } else if (rec->type == IGRAPH_ATTRIBUTE_STRING) {
                    str = (igraph_strvector_t*)rec->value;
                    igraph_strvector_destroy(str);
                    igraph_free(str);
                } else if (rec->type == IGRAPH_ATTRIBUTE_BOOLEAN) {
                    boolvec = (igraph_vector_bool_t*)rec->value;
                    igraph_vector_bool_destroy(boolvec);
                    igraph_free(boolvec);
                }
                igraph_free((char*)rec->name);
                igraph_free(rec);
            }
        }
    }
    igraph_vector_ptr_destroy(&attr->gal);
    igraph_vector_ptr_destroy(&attr->val);
    igraph_vector_ptr_destroy(&attr->eal);
    igraph_free(graph->attr);
    graph->attr = 0;
}

/* Almost the same as destroy, but we might have null pointers */

void igraph_i_cattribute_copy_free(igraph_i_cattributes_t *attr) {
    igraph_vector_ptr_t *als[3] = { &attr->gal, &attr->val, &attr->eal };
    long int i, n, a;
    igraph_vector_t *num;
    igraph_strvector_t *str;
    igraph_vector_bool_t *boolvec;
    igraph_attribute_record_t *rec;
    for (a = 0; a < 3; a++) {
        n = igraph_vector_ptr_size(als[a]);
        for (i = 0; i < n; i++) {
            rec = VECTOR(*als[a])[i];
            if (!rec) {
                continue;
            }
            if (rec->type == IGRAPH_ATTRIBUTE_NUMERIC) {
                num = (igraph_vector_t*)rec->value;
                igraph_vector_destroy(num);
                igraph_free(num);
            } else if (rec->type == IGRAPH_ATTRIBUTE_BOOLEAN) {
                boolvec = (igraph_vector_bool_t*)rec->value;
                igraph_vector_bool_destroy(boolvec);
                igraph_free(boolvec);
            } else if (rec->type == IGRAPH_ATTRIBUTE_STRING) {
                str = (igraph_strvector_t*)rec->value;
                igraph_strvector_destroy(str);
                igraph_free(str);
            }
            igraph_free((char*)rec->name);
            igraph_free(rec);
        }
    }
}

/* No reference counting here. If you use attributes in C you should
   know what you're doing. */

int igraph_i_cattribute_copy(igraph_t *to, const igraph_t *from,
                             igraph_bool_t ga, igraph_bool_t va, igraph_bool_t ea) {
    igraph_i_cattributes_t *attrfrom = from->attr, *attrto;
    igraph_vector_ptr_t *alto[3], *alfrom[3] = { &attrfrom->gal, &attrfrom->val,
                                                 &attrfrom->eal
                                               };
    long int i, n, a;
    igraph_bool_t copy[3] = { ga, va, ea };
    to->attr = attrto = igraph_Calloc(1, igraph_i_cattributes_t);
    if (!attrto) {
        IGRAPH_ERROR("Cannot copy attributes", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, attrto);
    IGRAPH_VECTOR_PTR_INIT_FINALLY(&attrto->gal, 0);
    IGRAPH_VECTOR_PTR_INIT_FINALLY(&attrto->val, 0);
    IGRAPH_VECTOR_PTR_INIT_FINALLY(&attrto->eal, 0);
    IGRAPH_FINALLY_CLEAN(3);
    IGRAPH_FINALLY(igraph_i_cattribute_copy_free, attrto);

    alto[0] = &attrto->gal; alto[1] = &attrto->val; alto[2] = &attrto->eal;
    for (a = 0; a < 3; a++) {
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
    return 0;
}

int igraph_i_cattribute_add_vertices(igraph_t *graph, long int nv,
                                     igraph_vector_ptr_t *nattr) {

    igraph_i_cattributes_t *attr = graph->attr;
    igraph_vector_ptr_t *val = &attr->val;
    long int length = igraph_vector_ptr_size(val);
    long int nattrno = nattr == NULL ? 0 : igraph_vector_ptr_size(nattr);
    long int origlen = igraph_vcount(graph) - nv;
    long int newattrs = 0, i;
    igraph_vector_t news;

    /* First add the new attributes if any */
    newattrs = 0;
    IGRAPH_VECTOR_INIT_FINALLY(&news, 0);
    for (i = 0; i < nattrno; i++) {
        igraph_attribute_record_t *nattr_entry = VECTOR(*nattr)[i];
        const char *nname = nattr_entry->name;
        long int j;
        igraph_bool_t l = igraph_i_cattribute_find(val, nname, &j);
        if (!l) {
            newattrs++;
            IGRAPH_CHECK(igraph_vector_push_back(&news, i));
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
            igraph_attribute_record_t *tmp = VECTOR(*nattr)[(long int)VECTOR(news)[i]];
            igraph_attribute_record_t *newrec = igraph_Calloc(1, igraph_attribute_record_t);
            igraph_attribute_type_t type = tmp->type;
            if (!newrec) {
                IGRAPH_ERROR("Cannot add attributes", IGRAPH_ENOMEM);
            }
            IGRAPH_FINALLY(igraph_free, newrec);
            newrec->type = type;
            newrec->name = strdup(tmp->name);
            if (!newrec->name) {
                IGRAPH_ERROR("Cannot add attributes", IGRAPH_ENOMEM);
            }
            IGRAPH_FINALLY(igraph_free, (char*)newrec->name);
            if (type == IGRAPH_ATTRIBUTE_NUMERIC) {
                igraph_vector_t *newnum = igraph_Calloc(1, igraph_vector_t);
                if (!newnum) {
                    IGRAPH_ERROR("Cannot add attributes", IGRAPH_ENOMEM);
                }
                IGRAPH_FINALLY(igraph_free, newnum);
                IGRAPH_VECTOR_INIT_FINALLY(newnum, origlen);
                newrec->value = newnum;
                igraph_vector_fill(newnum, IGRAPH_NAN);
            } else if (type == IGRAPH_ATTRIBUTE_STRING) {
                igraph_strvector_t *newstr = igraph_Calloc(1, igraph_strvector_t);
                if (!newstr) {
                    IGRAPH_ERROR("Cannot add attributes", IGRAPH_ENOMEM);
                }
                IGRAPH_FINALLY(igraph_free, newstr);
                IGRAPH_STRVECTOR_INIT_FINALLY(newstr, origlen);
                newrec->value = newstr;
            } else if (type == IGRAPH_ATTRIBUTE_BOOLEAN) {
                igraph_vector_bool_t *newbool = igraph_Calloc(1, igraph_vector_bool_t);
                if (!newbool) {
                    IGRAPH_ERROR("Cannot add attributes", IGRAPH_ENOMEM);
                }
                IGRAPH_FINALLY(igraph_free, newbool);
                IGRAPH_CHECK(igraph_vector_bool_init(newbool, origlen));
                IGRAPH_FINALLY(igraph_vector_bool_destroy, newbool);
                newrec->value = newbool;
                igraph_vector_bool_fill(newbool, 0);
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
        long int j;
        igraph_bool_t l = 0;
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
                IGRAPH_ERROR("Attribute types do not match", IGRAPH_EINVAL);
            }
            switch (oldrec->type) {
            case IGRAPH_ATTRIBUTE_NUMERIC:
                if (nv != igraph_vector_size(newnum)) {
                    IGRAPH_ERROR("Invalid numeric attribute length", IGRAPH_EINVAL);
                }
                IGRAPH_CHECK(igraph_vector_append(oldnum, newnum));
                break;
            case IGRAPH_ATTRIBUTE_STRING:
                if (nv != igraph_strvector_size(newstr)) {
                    IGRAPH_ERROR("Invalid string attribute length", IGRAPH_EINVAL);
                }
                IGRAPH_CHECK(igraph_strvector_append(oldstr, newstr));
                break;
            case IGRAPH_ATTRIBUTE_BOOLEAN:
                if (nv != igraph_vector_bool_size(newbool)) {
                    IGRAPH_ERROR("Invalid Boolean attribute length", IGRAPH_EINVAL);
                }
                IGRAPH_CHECK(igraph_vector_bool_append(oldbool, newbool));
                break;
            default:
                IGRAPH_WARNING("Invalid attribute type");
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

    igraph_vector_destroy(&news);
    IGRAPH_FINALLY_CLEAN(1);

    return 0;
}

void igraph_i_cattribute_permute_free(igraph_vector_ptr_t *v) {
    long int i, n = igraph_vector_ptr_size(v);
    for (i = 0; i < n; i++) {
        igraph_attribute_record_t *rec = VECTOR(*v)[i];
        igraph_Free(rec->name);
        if (rec->type == IGRAPH_ATTRIBUTE_NUMERIC) {
            igraph_vector_t *numv = (igraph_vector_t*) rec->value;
            igraph_vector_destroy(numv);
            igraph_Free(numv);
        } else if (rec->type == IGRAPH_ATTRIBUTE_STRING) {
            igraph_strvector_t *strv = (igraph_strvector_t*) rec->value;
            igraph_strvector_destroy(strv);
            igraph_Free(strv);
        } else if (rec->type == IGRAPH_ATTRIBUTE_BOOLEAN) {
            igraph_vector_bool_t *boolv = (igraph_vector_bool_t*) rec->value;
            igraph_vector_bool_destroy(boolv);
            igraph_Free(boolv);
        }
        igraph_Free(rec);
    }
    igraph_vector_ptr_clear(v);
}

int igraph_i_cattribute_permute_vertices(const igraph_t *graph,
        igraph_t *newgraph,
        const igraph_vector_t *idx) {

    if (graph == newgraph) {

        igraph_i_cattributes_t *attr = graph->attr;
        igraph_vector_ptr_t *val = &attr->val;
        long int valno = igraph_vector_ptr_size(val);
        long int i;

        for (i = 0; i < valno; i++) {
            igraph_attribute_record_t *oldrec = VECTOR(*val)[i];
            igraph_attribute_type_t type = oldrec->type;
            igraph_vector_t *num, *newnum;
            igraph_strvector_t *str, *newstr;
            igraph_vector_bool_t *oldbool, *newbool;
            switch (type) {
            case IGRAPH_ATTRIBUTE_NUMERIC:
                num = (igraph_vector_t*) oldrec->value;
                newnum = igraph_Calloc(1, igraph_vector_t);
                if (!newnum) {
                    IGRAPH_ERROR("Cannot permute vertex attributes", IGRAPH_ENOMEM);
                }
                IGRAPH_VECTOR_INIT_FINALLY(newnum, 0);
                igraph_vector_index(num, newnum, idx);
                oldrec->value = newnum;
                igraph_vector_destroy(num);
                igraph_Free(num);
                IGRAPH_FINALLY_CLEAN(1);
                break;
            case IGRAPH_ATTRIBUTE_BOOLEAN:
                oldbool = (igraph_vector_bool_t*) oldrec->value;
                newbool = igraph_Calloc(1, igraph_vector_bool_t);
                if (!newbool) {
                    IGRAPH_ERROR("Cannot permute vertex attributes", IGRAPH_ENOMEM);
                }
                IGRAPH_CHECK(igraph_vector_bool_init(newbool, 0));
                IGRAPH_FINALLY(igraph_vector_bool_destroy, newbool);
                igraph_vector_bool_index(oldbool, newbool, idx);
                oldrec->value = newbool;
                igraph_vector_bool_destroy(oldbool);
                igraph_Free(oldbool);
                IGRAPH_FINALLY_CLEAN(1);
                break;
            case IGRAPH_ATTRIBUTE_STRING:
                str = (igraph_strvector_t*)oldrec->value;
                newstr = igraph_Calloc(1, igraph_strvector_t);
                if (!newstr) {
                    IGRAPH_ERROR("Cannot permute vertex attributes", IGRAPH_ENOMEM);
                }
                IGRAPH_CHECK(igraph_strvector_init(newstr, 0));
                IGRAPH_FINALLY(igraph_strvector_destroy, newstr);
                igraph_strvector_index(str, newstr, idx);
                oldrec->value = newstr;
                igraph_strvector_destroy(str);
                igraph_Free(str);
                IGRAPH_FINALLY_CLEAN(1);
                break;
            default:
                IGRAPH_WARNING("Unknown edge attribute ignored");
            }
        }

    } else {
        igraph_i_cattributes_t *attr = graph->attr;
        igraph_vector_ptr_t *val = &attr->val;
        long int valno = igraph_vector_ptr_size(val);
        long int i;

        /* New vertex attributes */
        igraph_i_cattributes_t *new_attr = newgraph->attr;
        igraph_vector_ptr_t *new_val = &new_attr->val;
        if (igraph_vector_ptr_size(new_val) != 0) {
            IGRAPH_ERROR("Vertex attributes were already copied",
                         IGRAPH_EATTRIBUTES);
        }
        IGRAPH_CHECK(igraph_vector_ptr_resize(new_val, valno));

        IGRAPH_FINALLY(igraph_i_cattribute_permute_free, new_val);

        for (i = 0; i < valno; i++) {
            igraph_attribute_record_t *oldrec = VECTOR(*val)[i];
            igraph_attribute_type_t type = oldrec->type;
            igraph_vector_t *num, *newnum;
            igraph_strvector_t *str, *newstr;
            igraph_vector_bool_t *oldbool, *newbool;

            /* The record itself */
            igraph_attribute_record_t *new_rec =
                igraph_Calloc(1, igraph_attribute_record_t);
            if (!new_rec) {
                IGRAPH_ERROR("Cannot create vertex attributes", IGRAPH_ENOMEM);
            }
            new_rec->name = strdup(oldrec->name);
            new_rec->type = oldrec->type;
            VECTOR(*new_val)[i] = new_rec;

            /* The data */
            switch (type) {
            case IGRAPH_ATTRIBUTE_NUMERIC:
                num = (igraph_vector_t*)oldrec->value;
                newnum = igraph_Calloc(1, igraph_vector_t);
                if (!newnum) {
                    IGRAPH_ERROR("Cannot permute vertex attributes", IGRAPH_ENOMEM);
                }
                IGRAPH_VECTOR_INIT_FINALLY(newnum, 0);
                igraph_vector_index(num, newnum, idx);
                new_rec->value = newnum;
                IGRAPH_FINALLY_CLEAN(1);
                break;
            case IGRAPH_ATTRIBUTE_BOOLEAN:
                oldbool = (igraph_vector_bool_t*)oldrec->value;
                newbool = igraph_Calloc(1, igraph_vector_bool_t);
                if (!newbool) {
                    IGRAPH_ERROR("Cannot permute vertex attributes", IGRAPH_ENOMEM);
                }
                IGRAPH_CHECK(igraph_vector_bool_init(newbool, 0));
                IGRAPH_FINALLY(igraph_vector_bool_destroy, newbool);
                igraph_vector_bool_index(oldbool, newbool, idx);
                new_rec->value = newbool;
                IGRAPH_FINALLY_CLEAN(1);
                break;
            case IGRAPH_ATTRIBUTE_STRING:
                str = (igraph_strvector_t*)oldrec->value;
                newstr = igraph_Calloc(1, igraph_strvector_t);
                if (!newstr) {
                    IGRAPH_ERROR("Cannot permute vertex attributes", IGRAPH_ENOMEM);
                }
                IGRAPH_CHECK(igraph_strvector_init(newstr, 0));
                IGRAPH_FINALLY(igraph_strvector_destroy, newstr);
                igraph_strvector_index(str, newstr, idx);
                new_rec->value = newstr;
                IGRAPH_FINALLY_CLEAN(1);
                break;
            default:
                IGRAPH_WARNING("Unknown vertex attribute ignored");
            }
        }
    }

    IGRAPH_FINALLY_CLEAN(1);
    return 0;
}

typedef int igraph_cattributes_combine_num_t(const igraph_vector_t *input,
        igraph_real_t *output);

typedef int igraph_cattributes_combine_str_t(const igraph_strvector_t *input,
        char **output);

typedef int igraph_cattributes_combine_bool_t(const igraph_vector_bool_t *input,
        igraph_bool_t *output);

int igraph_i_cattributes_cn_sum(const igraph_attribute_record_t *oldrec,
                                igraph_attribute_record_t * newrec,
                                const igraph_vector_ptr_t *merges) {
    const igraph_vector_t *oldv = oldrec->value;
    igraph_vector_t *newv = igraph_Calloc(1, igraph_vector_t);
    long int newlen = igraph_vector_ptr_size(merges);
    long int i;

    if (!newv) {
        IGRAPH_ERROR("Cannot combine attributes", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, newv);
    IGRAPH_VECTOR_INIT_FINALLY(newv, newlen);

    for (i = 0; i < newlen; i++) {
        igraph_real_t s = 0.0;
        igraph_vector_t *idx = VECTOR(*merges)[i];
        long int j, n = igraph_vector_size(idx);
        for (j = 0; j < n; j++) {
            long int x = (long int) VECTOR(*idx)[j];
            s += VECTOR(*oldv)[x];
        }
        VECTOR(*newv)[i] = s;
    }

    IGRAPH_FINALLY_CLEAN(2);
    newrec->value = newv;

    return 0;
}

int igraph_i_cattributes_cn_prod(const igraph_attribute_record_t *oldrec,
                                 igraph_attribute_record_t * newrec,
                                 const igraph_vector_ptr_t *merges) {
    const igraph_vector_t *oldv = oldrec->value;
    igraph_vector_t *newv = igraph_Calloc(1, igraph_vector_t);
    long int newlen = igraph_vector_ptr_size(merges);
    long int i;

    if (!newv) {
        IGRAPH_ERROR("Cannot combine attributes", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, newv);
    IGRAPH_VECTOR_INIT_FINALLY(newv, newlen);

    for (i = 0; i < newlen; i++) {
        igraph_real_t s = 1.0;
        igraph_vector_t *idx = VECTOR(*merges)[i];
        long int j, n = igraph_vector_size(idx);
        for (j = 0; j < n; j++) {
            long int x = (long int) VECTOR(*idx)[j];
            s *= VECTOR(*oldv)[x];
        }
        VECTOR(*newv)[i] = s;
    }

    IGRAPH_FINALLY_CLEAN(2);
    newrec->value = newv;

    return 0;
}

int igraph_i_cattributes_cn_min(const igraph_attribute_record_t *oldrec,
                                igraph_attribute_record_t * newrec,
                                const igraph_vector_ptr_t *merges) {
    const igraph_vector_t *oldv = oldrec->value;
    igraph_vector_t *newv = igraph_Calloc(1, igraph_vector_t);
    long int newlen = igraph_vector_ptr_size(merges);
    long int i;
    igraph_real_t nan = IGRAPH_NAN;

    if (!newv) {
        IGRAPH_ERROR("Cannot combine attributes", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, newv);
    IGRAPH_VECTOR_INIT_FINALLY(newv, newlen);

    for (i = 0; i < newlen; i++) {
        igraph_vector_t *idx = VECTOR(*merges)[i];
        long int j, n = igraph_vector_size(idx);
        igraph_real_t m = n > 0 ? VECTOR(*oldv)[ (long int) VECTOR(*idx)[0] ] : nan;
        for (j = 1; j < n; j++) {
            long int x = (long int) VECTOR(*idx)[j];
            igraph_real_t val = VECTOR(*oldv)[x];
            if (val < m) {
                m = val;
            }
        }
        VECTOR(*newv)[i] = m;
    }

    IGRAPH_FINALLY_CLEAN(2);
    newrec->value = newv;

    return 0;
}

int igraph_i_cattributes_cn_max(const igraph_attribute_record_t *oldrec,
                                igraph_attribute_record_t * newrec,
                                const igraph_vector_ptr_t *merges) {
    const igraph_vector_t *oldv = oldrec->value;
    igraph_vector_t *newv = igraph_Calloc(1, igraph_vector_t);
    long int newlen = igraph_vector_ptr_size(merges);
    long int i;
    igraph_real_t nan = IGRAPH_NAN;

    if (!newv) {
        IGRAPH_ERROR("Cannot combine attributes", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, newv);
    IGRAPH_VECTOR_INIT_FINALLY(newv, newlen);

    for (i = 0; i < newlen; i++) {
        igraph_vector_t *idx = VECTOR(*merges)[i];
        long int j, n = igraph_vector_size(idx);
        igraph_real_t m = n > 0 ? VECTOR(*oldv)[ (long int) VECTOR(*idx)[0] ] : nan;
        for (j = 1; j < n; j++) {
            long int x = (long int) VECTOR(*idx)[j];
            igraph_real_t val = VECTOR(*oldv)[x];
            if (val > m) {
                m = val;
            }
        }
        VECTOR(*newv)[i] = m;
    }

    IGRAPH_FINALLY_CLEAN(2);
    newrec->value = newv;

    return 0;
}

int igraph_i_cattributes_cn_random(const igraph_attribute_record_t *oldrec,
                                   igraph_attribute_record_t * newrec,
                                   const igraph_vector_ptr_t *merges) {

    const igraph_vector_t *oldv = oldrec->value;
    igraph_vector_t *newv = igraph_Calloc(1, igraph_vector_t);
    long int newlen = igraph_vector_ptr_size(merges);
    long int i;
    igraph_real_t nan = IGRAPH_NAN;

    if (!newv) {
        IGRAPH_ERROR("Cannot combine attributes", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, newv);
    IGRAPH_VECTOR_INIT_FINALLY(newv, newlen);

    RNG_BEGIN();

    for (i = 0; i < newlen; i++) {
        igraph_vector_t *idx = VECTOR(*merges)[i];
        long int n = igraph_vector_size(idx);
        if (n == 0) {
            VECTOR(*newv)[i] = nan;
        } else if (n == 1) {
            VECTOR(*newv)[i] = VECTOR(*oldv)[ (long int) VECTOR(*idx)[0] ];
        } else {
            long int r = RNG_INTEGER(0, n - 1);
            VECTOR(*newv)[i] = VECTOR(*oldv)[ (long int) VECTOR(*idx)[r] ];
        }
    }

    RNG_END();

    IGRAPH_FINALLY_CLEAN(2);
    newrec->value = newv;

    return 0;
}

int igraph_i_cattributes_cn_first(const igraph_attribute_record_t *oldrec,
                                  igraph_attribute_record_t * newrec,
                                  const igraph_vector_ptr_t *merges) {

    const igraph_vector_t *oldv = oldrec->value;
    igraph_vector_t *newv = igraph_Calloc(1, igraph_vector_t);
    long int newlen = igraph_vector_ptr_size(merges);
    long int i;
    igraph_real_t nan = IGRAPH_NAN;

    if (!newv) {
        IGRAPH_ERROR("Cannot combine attributes", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, newv);
    IGRAPH_VECTOR_INIT_FINALLY(newv, newlen);

    for (i = 0; i < newlen; i++) {
        igraph_vector_t *idx = VECTOR(*merges)[i];
        long int n = igraph_vector_size(idx);
        if (n == 0) {
            VECTOR(*newv)[i] = nan;
        } else {
            VECTOR(*newv)[i] = VECTOR(*oldv)[ (long int) VECTOR(*idx)[0] ];
        }
    }

    IGRAPH_FINALLY_CLEAN(2);
    newrec->value = newv;

    return 0;
}

int igraph_i_cattributes_cn_last(const igraph_attribute_record_t *oldrec,
                                 igraph_attribute_record_t * newrec,
                                 const igraph_vector_ptr_t *merges) {

    const igraph_vector_t *oldv = oldrec->value;
    igraph_vector_t *newv = igraph_Calloc(1, igraph_vector_t);
    long int newlen = igraph_vector_ptr_size(merges);
    long int i;
    igraph_real_t nan = IGRAPH_NAN;

    if (!newv) {
        IGRAPH_ERROR("Cannot combine attributes", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, newv);
    IGRAPH_VECTOR_INIT_FINALLY(newv, newlen);

    for (i = 0; i < newlen; i++) {
        igraph_vector_t *idx = VECTOR(*merges)[i];
        long int n = igraph_vector_size(idx);
        if (n == 0) {
            VECTOR(*newv)[i] = nan;
        } else {
            VECTOR(*newv)[i] = VECTOR(*oldv)[ (long int) VECTOR(*idx)[n - 1] ];
        }
    }

    IGRAPH_FINALLY_CLEAN(2);
    newrec->value = newv;

    return 0;
}

int igraph_i_cattributes_cn_mean(const igraph_attribute_record_t *oldrec,
                                 igraph_attribute_record_t * newrec,
                                 const igraph_vector_ptr_t *merges) {
    const igraph_vector_t *oldv = oldrec->value;
    igraph_vector_t *newv = igraph_Calloc(1, igraph_vector_t);
    long int newlen = igraph_vector_ptr_size(merges);
    long int i;
    igraph_real_t nan = IGRAPH_NAN;

    if (!newv) {
        IGRAPH_ERROR("Cannot combine attributes", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, newv);
    IGRAPH_VECTOR_INIT_FINALLY(newv, newlen);

    for (i = 0; i < newlen; i++) {
        igraph_vector_t *idx = VECTOR(*merges)[i];
        long int j, n = igraph_vector_size(idx);
        igraph_real_t s = n > 0 ? 0.0 : nan;
        for (j = 0; j < n; j++) {
            long int x = (long int) VECTOR(*idx)[j];
            s += VECTOR(*oldv)[x];
        }
        if (n > 0) {
            s = s / n;
        }
        VECTOR(*newv)[i] = s;
    }

    IGRAPH_FINALLY_CLEAN(2);
    newrec->value = newv;

    return 0;
}

int igraph_i_cattributes_cn_func(const igraph_attribute_record_t *oldrec,
                                 igraph_attribute_record_t *newrec,
                                 const igraph_vector_ptr_t *merges,
                                 igraph_cattributes_combine_num_t *func) {

    const igraph_vector_t *oldv = oldrec->value;
    long int newlen = igraph_vector_ptr_size(merges);
    long int i;
    igraph_vector_t *newv = igraph_Calloc(1, igraph_vector_t);
    igraph_vector_t values;

    if (!newv) {
        IGRAPH_ERROR("Cannot combine attributes", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, newv);
    IGRAPH_VECTOR_INIT_FINALLY(newv, newlen);

    IGRAPH_VECTOR_INIT_FINALLY(&values, 0);

    for (i = 0; i < newlen; i++) {
        igraph_vector_t *idx = VECTOR(*merges)[i];
        long int j, n = igraph_vector_size(idx);
        igraph_real_t res;
        IGRAPH_CHECK(igraph_vector_resize(&values, n));
        for (j = 0; j < n; j++) {
            long int x = (long int) VECTOR(*idx)[j];
            VECTOR(values)[j] = VECTOR(*oldv)[x];
        }
        IGRAPH_CHECK(func(&values, &res));
        VECTOR(*newv)[i] = res;
    }

    igraph_vector_destroy(&values);
    IGRAPH_FINALLY_CLEAN(3);
    newrec->value = newv;

    return 0;
}

int igraph_i_cattributes_cb_random(const igraph_attribute_record_t *oldrec,
                                   igraph_attribute_record_t * newrec,
                                   const igraph_vector_ptr_t *merges) {

    const igraph_vector_bool_t *oldv = oldrec->value;
    igraph_vector_bool_t *newv = igraph_Calloc(1, igraph_vector_bool_t);
    long int newlen = igraph_vector_ptr_size(merges);
    long int i;

    if (!newv) {
        IGRAPH_ERROR("Cannot combine attributes", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, newv);
    IGRAPH_CHECK(igraph_vector_bool_init(newv, newlen));
    IGRAPH_FINALLY(igraph_vector_bool_destroy, newv);

    RNG_BEGIN();

    for (i = 0; i < newlen; i++) {
        igraph_vector_t *idx = VECTOR(*merges)[i];
        long int n = igraph_vector_size(idx);
        if (n == 0) {
            VECTOR(*newv)[i] = 0;
        } else if (n == 1) {
            VECTOR(*newv)[i] = VECTOR(*oldv)[ (long int) VECTOR(*idx)[0] ];
        } else {
            long int r = RNG_INTEGER(0, n - 1);
            VECTOR(*newv)[i] = VECTOR(*oldv)[ (long int) VECTOR(*idx)[r] ];
        }
    }

    RNG_END();

    IGRAPH_FINALLY_CLEAN(2);
    newrec->value = newv;

    return 0;
}

int igraph_i_cattributes_cb_first(const igraph_attribute_record_t *oldrec,
                                  igraph_attribute_record_t * newrec,
                                  const igraph_vector_ptr_t *merges) {

    const igraph_vector_bool_t *oldv = oldrec->value;
    igraph_vector_bool_t *newv = igraph_Calloc(1, igraph_vector_bool_t);
    long int newlen = igraph_vector_ptr_size(merges);
    long int i;

    if (!newv) {
        IGRAPH_ERROR("Cannot combine attributes", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, newv);
    IGRAPH_CHECK(igraph_vector_bool_init(newv, newlen));
    IGRAPH_FINALLY(igraph_vector_bool_destroy, newv);

    for (i = 0; i < newlen; i++) {
        igraph_vector_t *idx = VECTOR(*merges)[i];
        long int n = igraph_vector_size(idx);
        if (n == 0) {
            VECTOR(*newv)[i] = 0;
        } else {
            VECTOR(*newv)[i] = VECTOR(*oldv)[ (long int) VECTOR(*idx)[0] ];
        }
    }

    IGRAPH_FINALLY_CLEAN(2);
    newrec->value = newv;

    return 0;
}

int igraph_i_cattributes_cb_last(const igraph_attribute_record_t *oldrec,
                                 igraph_attribute_record_t * newrec,
                                 const igraph_vector_ptr_t *merges) {

    const igraph_vector_bool_t *oldv = oldrec->value;
    igraph_vector_bool_t *newv = igraph_Calloc(1, igraph_vector_bool_t);
    long int newlen = igraph_vector_ptr_size(merges);
    long int i;

    if (!newv) {
        IGRAPH_ERROR("Cannot combine attributes", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, newv);
    IGRAPH_CHECK(igraph_vector_bool_init(newv, newlen));
    IGRAPH_FINALLY(igraph_vector_bool_destroy, newv);

    for (i = 0; i < newlen; i++) {
        igraph_vector_t *idx = VECTOR(*merges)[i];
        long int n = igraph_vector_size(idx);
        if (n == 0) {
            VECTOR(*newv)[i] = 0;
        } else {
            VECTOR(*newv)[i] = VECTOR(*oldv)[ (long int) VECTOR(*idx)[n - 1] ];
        }
    }

    IGRAPH_FINALLY_CLEAN(2);
    newrec->value = newv;

    return 0;
}

int igraph_i_cattributes_cb_all_is_true(const igraph_attribute_record_t *oldrec,
                                        igraph_attribute_record_t * newrec,
                                        const igraph_vector_ptr_t *merges) {

    const igraph_vector_bool_t *oldv = oldrec->value;
    igraph_vector_bool_t *newv = igraph_Calloc(1, igraph_vector_bool_t);
    long int newlen = igraph_vector_ptr_size(merges);
    long int i, j, n, x;

    if (!newv) {
        IGRAPH_ERROR("Cannot combine attributes", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, newv);
    IGRAPH_CHECK(igraph_vector_bool_init(newv, newlen));
    IGRAPH_FINALLY(igraph_vector_bool_destroy, newv);

    for (i = 0; i < newlen; i++) {
        igraph_vector_t *idx = VECTOR(*merges)[i];
        n = igraph_vector_size(idx);
        VECTOR(*newv)[i] = 1;
        for (j = 0; j < n; j++) {
            x = (long int) VECTOR(*idx)[j];
            if (!VECTOR(*oldv)[x]) {
                VECTOR(*newv)[i] = 0;
                break;
            }
        }
    }

    IGRAPH_FINALLY_CLEAN(2);
    newrec->value = newv;

    return 0;
}

int igraph_i_cattributes_cb_any_is_true(const igraph_attribute_record_t *oldrec,
                                        igraph_attribute_record_t * newrec,
                                        const igraph_vector_ptr_t *merges) {

    const igraph_vector_bool_t *oldv = oldrec->value;
    igraph_vector_bool_t *newv = igraph_Calloc(1, igraph_vector_bool_t);
    long int newlen = igraph_vector_ptr_size(merges);
    long int i, j, n, x;

    if (!newv) {
        IGRAPH_ERROR("Cannot combine attributes", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, newv);
    IGRAPH_CHECK(igraph_vector_bool_init(newv, newlen));
    IGRAPH_FINALLY(igraph_vector_bool_destroy, newv);

    for (i = 0; i < newlen; i++) {
        igraph_vector_t *idx = VECTOR(*merges)[i];
        n = igraph_vector_size(idx);
        VECTOR(*newv)[i] = 0;
        for (j = 0; j < n; j++) {
            x = (long int) VECTOR(*idx)[j];
            if (VECTOR(*oldv)[x]) {
                VECTOR(*newv)[i] = 1;
                break;
            }
        }
    }

    IGRAPH_FINALLY_CLEAN(2);
    newrec->value = newv;

    return 0;
}

int igraph_i_cattributes_cb_majority(const igraph_attribute_record_t *oldrec,
                                     igraph_attribute_record_t * newrec,
                                     const igraph_vector_ptr_t *merges) {

    const igraph_vector_bool_t *oldv = oldrec->value;
    igraph_vector_bool_t *newv = igraph_Calloc(1, igraph_vector_bool_t);
    long int newlen = igraph_vector_ptr_size(merges);
    long int i, j, n, x, num_trues;

    if (!newv) {
        IGRAPH_ERROR("Cannot combine attributes", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, newv);
    IGRAPH_CHECK(igraph_vector_bool_init(newv, newlen));
    IGRAPH_FINALLY(igraph_vector_bool_destroy, newv);

    RNG_BEGIN();

    for (i = 0; i < newlen; i++) {
        igraph_vector_t *idx = VECTOR(*merges)[i];

        n = igraph_vector_size(idx);

        num_trues = 0;
        for (j = 0; j < n; j++) {
            x = (long int) VECTOR(*idx)[j];
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

    return 0;
}

int igraph_i_cattributes_cb_func(const igraph_attribute_record_t *oldrec,
                                 igraph_attribute_record_t *newrec,
                                 const igraph_vector_ptr_t *merges,
                                 igraph_cattributes_combine_bool_t *func) {

    const igraph_vector_bool_t *oldv = oldrec->value;
    long int newlen = igraph_vector_ptr_size(merges);
    long int i;
    igraph_vector_bool_t *newv = igraph_Calloc(1, igraph_vector_bool_t);
    igraph_vector_bool_t values;

    if (!newv) {
        IGRAPH_ERROR("Cannot combine attributes", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, newv);
    IGRAPH_CHECK(igraph_vector_bool_init(newv, newlen));
    IGRAPH_FINALLY(igraph_vector_bool_destroy, newv);

    IGRAPH_CHECK(igraph_vector_bool_init(&values, 0));
    IGRAPH_FINALLY(igraph_vector_bool_destroy, newv);

    for (i = 0; i < newlen; i++) {
        igraph_vector_t *idx = VECTOR(*merges)[i];
        long int j, n = igraph_vector_size(idx);
        igraph_bool_t res;

        IGRAPH_CHECK(igraph_vector_bool_resize(&values, n));
        for (j = 0; j < n; j++) {
            long int x = (long int) VECTOR(*idx)[j];
            VECTOR(values)[j] = VECTOR(*oldv)[x];
        }

        IGRAPH_CHECK(func(&values, &res));
        VECTOR(*newv)[i] = res;
    }

    igraph_vector_bool_destroy(&values);
    IGRAPH_FINALLY_CLEAN(3);
    newrec->value = newv;

    return 0;
}

int igraph_i_cattributes_sn_random(const igraph_attribute_record_t *oldrec,
                                   igraph_attribute_record_t *newrec,
                                   const igraph_vector_ptr_t *merges) {

    const igraph_strvector_t *oldv = oldrec->value;
    long int newlen = igraph_vector_ptr_size(merges);
    long int i;
    igraph_strvector_t *newv = igraph_Calloc(1, igraph_strvector_t);

    if (!newv) {
        IGRAPH_ERROR("Cannot combine attributes", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, newv);
    IGRAPH_CHECK(igraph_strvector_init(newv, newlen));
    IGRAPH_FINALLY(igraph_strvector_destroy, newv);

    RNG_BEGIN();

    for (i = 0; i < newlen; i++) {
        igraph_vector_t *idx = VECTOR(*merges)[i];
        long int n = igraph_vector_size(idx);
        char *tmp;
        if (n == 0) {
            IGRAPH_CHECK(igraph_strvector_set(newv, i, ""));
        } else if (n == 1) {
            igraph_strvector_get(oldv, 0, &tmp);
            IGRAPH_CHECK(igraph_strvector_set(newv, i, tmp));
        } else {
            long int r = RNG_INTEGER(0, n - 1);
            igraph_strvector_get(oldv, r, &tmp);
            IGRAPH_CHECK(igraph_strvector_set(newv, i, tmp));
        }
    }

    RNG_END();

    IGRAPH_FINALLY_CLEAN(2);
    newrec->value = newv;

    return 0;
}

int igraph_i_cattributes_sn_first(const igraph_attribute_record_t *oldrec,
                                  igraph_attribute_record_t *newrec,
                                  const igraph_vector_ptr_t *merges) {

    const igraph_strvector_t *oldv = oldrec->value;
    long int i, newlen = igraph_vector_ptr_size(merges);
    igraph_strvector_t *newv = igraph_Calloc(1, igraph_strvector_t);

    if (!newv) {
        IGRAPH_ERROR("Cannot combine attributes", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, newv);
    IGRAPH_CHECK(igraph_strvector_init(newv, newlen));
    IGRAPH_FINALLY(igraph_strvector_destroy, newv);

    for (i = 0; i < newlen; i++) {
        igraph_vector_t *idx = VECTOR(*merges)[i];
        long int n = igraph_vector_size(idx);
        if (n == 0) {
            IGRAPH_CHECK(igraph_strvector_set(newv, i, ""));
        } else {
            char *tmp;
            igraph_strvector_get(oldv, (long int) VECTOR(*idx)[0], &tmp);
            IGRAPH_CHECK(igraph_strvector_set(newv, i, tmp));
        }
    }

    IGRAPH_FINALLY_CLEAN(2);
    newrec->value = newv;

    return 0;
}

int igraph_i_cattributes_sn_last(const igraph_attribute_record_t *oldrec,
                                 igraph_attribute_record_t *newrec,
                                 const igraph_vector_ptr_t *merges) {

    const igraph_strvector_t *oldv = oldrec->value;
    long int i, newlen = igraph_vector_ptr_size(merges);
    igraph_strvector_t *newv = igraph_Calloc(1, igraph_strvector_t);

    if (!newv) {
        IGRAPH_ERROR("Cannot combine attributes", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, newv);
    IGRAPH_CHECK(igraph_strvector_init(newv, newlen));
    IGRAPH_FINALLY(igraph_strvector_destroy, newv);

    for (i = 0; i < newlen; i++) {
        igraph_vector_t *idx = VECTOR(*merges)[i];
        long int n = igraph_vector_size(idx);
        if (n == 0) {
            IGRAPH_CHECK(igraph_strvector_set(newv, i, ""));
        } else {
            char *tmp;
            igraph_strvector_get(oldv, (long int) VECTOR(*idx)[n - 1], &tmp);
            IGRAPH_CHECK(igraph_strvector_set(newv, i, tmp));
        }
    }

    IGRAPH_FINALLY_CLEAN(2);
    newrec->value = newv;

    return 0;
}

int igraph_i_cattributes_sn_concat(const igraph_attribute_record_t *oldrec,
                                   igraph_attribute_record_t *newrec,
                                   const igraph_vector_ptr_t *merges) {

    const igraph_strvector_t *oldv = oldrec->value;
    long int i, newlen = igraph_vector_ptr_size(merges);
    igraph_strvector_t *newv = igraph_Calloc(1, igraph_strvector_t);

    if (!newv) {
        IGRAPH_ERROR("Cannot combine attributes", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, newv);
    IGRAPH_CHECK(igraph_strvector_init(newv, newlen));
    IGRAPH_FINALLY(igraph_strvector_destroy, newv);

    for (i = 0; i < newlen; i++) {
        igraph_vector_t *idx = VECTOR(*merges)[i];
        long int j, n = igraph_vector_size(idx);
        size_t len = 0;
        char *tmp, *tmp2;
        for (j = 0; j < n; j++) {
            igraph_strvector_get(oldv, j, &tmp);
            len += strlen(tmp);
        }
        tmp2 = igraph_Calloc(len + 1, char);
        if (!tmp2) {
            IGRAPH_ERROR("Cannot combine attributes", IGRAPH_ENOMEM);
        }
        IGRAPH_FINALLY(igraph_free, tmp2);
        len = 0;
        for (j = 0; j < n; j++) {
            igraph_strvector_get(oldv, j, &tmp);
            strcpy(tmp2 + len, tmp);
            len += strlen(tmp);
        }

        IGRAPH_CHECK(igraph_strvector_set(newv, i, tmp2));
        igraph_Free(tmp2);
        IGRAPH_FINALLY_CLEAN(1);
    }

    IGRAPH_FINALLY_CLEAN(2);
    newrec->value = newv;

    return 0;
}

int igraph_i_cattributes_sn_func(const igraph_attribute_record_t *oldrec,
                                 igraph_attribute_record_t *newrec,
                                 const igraph_vector_ptr_t *merges,
                                 igraph_cattributes_combine_str_t *func) {

    const igraph_strvector_t *oldv = oldrec->value;
    long int newlen = igraph_vector_ptr_size(merges);
    long int i;
    igraph_strvector_t *newv = igraph_Calloc(1, igraph_strvector_t);
    igraph_strvector_t values;

    if (!newv) {
        IGRAPH_ERROR("Cannot combine attributes", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, newv);
    IGRAPH_CHECK(igraph_strvector_init(newv, newlen));
    IGRAPH_FINALLY(igraph_strvector_destroy, newv);

    IGRAPH_CHECK(igraph_strvector_init(newv, 0));
    IGRAPH_FINALLY(igraph_strvector_destroy, &values);

    for (i = 0; i < newlen; i++) {
        igraph_vector_t *idx = VECTOR(*merges)[i];
        long int j, n = igraph_vector_size(idx);
        char *res;
        IGRAPH_CHECK(igraph_strvector_resize(&values, n));
        for (j = 0; j < n; j++) {
            long int x = (long int) VECTOR(*idx)[j];
            char *elem;
            igraph_strvector_get(oldv, x, &elem);
            IGRAPH_CHECK(igraph_strvector_set(newv, j, elem));
        }
        IGRAPH_CHECK(func(&values, &res));
        IGRAPH_FINALLY(igraph_free, res);
        IGRAPH_CHECK(igraph_strvector_set(newv, i, res));
        IGRAPH_FINALLY_CLEAN(1);
        igraph_Free(res);
    }

    igraph_strvector_destroy(&values);
    IGRAPH_FINALLY_CLEAN(3);
    newrec->value = newv;

    return 0;
}


int igraph_i_cattribute_combine_vertices(const igraph_t *graph,
        igraph_t *newgraph,
        const igraph_vector_ptr_t *merges,
        const igraph_attribute_combination_t *comb) {

    igraph_i_cattributes_t *attr = graph->attr;
    igraph_i_cattributes_t *toattr = newgraph->attr;
    igraph_vector_ptr_t *val = &attr->val;
    igraph_vector_ptr_t *new_val = &toattr->val;
    long int valno = igraph_vector_ptr_size(val);
    long int i, j, keepno = 0;
    int *TODO;
    igraph_function_pointer_t *funcs;

    TODO = igraph_Calloc(valno, int);
    if (!TODO) {
        IGRAPH_ERROR("Cannot combine vertex attributes",
                     IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, TODO);
    funcs = igraph_Calloc(valno, igraph_function_pointer_t);
    if (!funcs) {
        IGRAPH_ERROR("Cannot combine vertex attributes",
                     IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, funcs);

    for (i = 0; i < valno; i++) {
        igraph_attribute_record_t *oldrec = VECTOR(*val)[i];
        const char *name = oldrec->name;
        igraph_attribute_combination_type_t todo;
        igraph_function_pointer_t voidfunc;
        igraph_attribute_combination_query(comb, name, &todo, &voidfunc);
        TODO[i] = todo;
        funcs[i] = voidfunc;
        if (todo != IGRAPH_ATTRIBUTE_COMBINE_IGNORE) {
            keepno++;
        }
    }

    IGRAPH_CHECK(igraph_vector_ptr_resize(new_val, keepno));
    IGRAPH_FINALLY(igraph_i_cattribute_permute_free, new_val);

    for (i = 0, j = 0; i < valno; i++) {
        igraph_attribute_record_t *newrec, *oldrec = VECTOR(*val)[i];
        const char *name = oldrec->name;
        igraph_attribute_combination_type_t todo =
            (igraph_attribute_combination_type_t) (TODO[i]);
        igraph_attribute_type_t type = oldrec->type;
        igraph_cattributes_combine_num_t *numfunc =
            (igraph_cattributes_combine_num_t*) funcs[i];
        igraph_cattributes_combine_str_t *strfunc =
            (igraph_cattributes_combine_str_t*) funcs[i];
        igraph_cattributes_combine_bool_t *boolfunc =
            (igraph_cattributes_combine_bool_t*) funcs[i];

        if (todo == IGRAPH_ATTRIBUTE_COMBINE_DEFAULT ||
            todo == IGRAPH_ATTRIBUTE_COMBINE_IGNORE) {
            continue;
        }

        newrec = igraph_Calloc(1, igraph_attribute_record_t);
        if (!newrec) {
            IGRAPH_ERROR("Cannot combine vertex attributes",
                         IGRAPH_ENOMEM);
        }
        newrec->name = strdup(name);
        newrec->type = type;
        VECTOR(*new_val)[j] = newrec;

        if (type == IGRAPH_ATTRIBUTE_NUMERIC) {
            switch (todo) {
            case IGRAPH_ATTRIBUTE_COMBINE_FUNCTION:
                IGRAPH_CHECK(igraph_i_cattributes_cn_func(oldrec, newrec, merges,
                             numfunc));
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
        } else if (type == IGRAPH_ATTRIBUTE_BOOLEAN) {
            switch (todo) {
            case IGRAPH_ATTRIBUTE_COMBINE_FUNCTION:
                IGRAPH_CHECK(igraph_i_cattributes_cb_func(oldrec, newrec, merges,
                             boolfunc));
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
        } else if (type == IGRAPH_ATTRIBUTE_STRING) {
            switch (todo) {
            case IGRAPH_ATTRIBUTE_COMBINE_FUNCTION:
                IGRAPH_CHECK(igraph_i_cattributes_sn_func(oldrec, newrec, merges,
                             strfunc));
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

        j++;
    }

    igraph_free(funcs);
    igraph_free(TODO);
    IGRAPH_FINALLY_CLEAN(2);

    return 0;
}

/* void igraph_i_cattribute_delete_vertices(igraph_t *graph, */
/*                     const igraph_vector_t *eidx, */
/*                     const igraph_vector_t *vidx) { */

/*   igraph_i_cattributes_t *attr=graph->attr; */
/*   igraph_vector_ptr_t *val=&attr->val; */
/*   igraph_vector_ptr_t *eal=&attr->eal; */
/*   long int valno=igraph_vector_ptr_size(val); */
/*   long int ealno=igraph_vector_ptr_size(eal); */
/*   long int i; */
/*   long int origlen, newlen; */

/*   /\* Vertices *\/ */
/*   origlen=igraph_vector_size(vidx); */
/*   newlen=0; */
/*   for (i=0; i<origlen; i++) { */
/*     if (VECTOR(*vidx)[i]>0) { */
/*       newlen++; */
/*     } */
/*   } */
/*   for (i=0; i<valno; i++) { */
/*     igraph_attribute_record_t *oldrec=VECTOR(*val)[i]; */
/*     igraph_attribute_type_t type=oldrec->type; */
/*     igraph_vector_t *num=(igraph_vector_t*)oldrec->value; */
/*     igraph_strvector_t *str=(igraph_strvector_t*)oldrec->value; */
/*     switch (type) { */
/*     case IGRAPH_ATTRIBUTE_NUMERIC: */
/*       igraph_vector_permdelete(num, vidx, origlen-newlen); */
/*       break; */
/*     case IGRAPH_ATTRIBUTE_STRING: */
/*       igraph_strvector_permdelete(str, vidx, origlen-newlen); */
/*       break; */
/*     default: */
/*       IGRAPH_WARNING("Unknown vertex attribute ignored"); */
/*     } */
/*   } */

/*   /\* Edges *\/ */
/*   origlen=igraph_vector_size(eidx); */
/*   newlen=0; */
/*   for (i=0; i<origlen; i++) { */
/*     if (VECTOR(*eidx)[i]>0) { */
/*       newlen++; */
/*     } */
/*   } */
/*   for (i=0; i<ealno; i++) { */
/*     igraph_attribute_record_t *oldrec=VECTOR(*eal)[i]; */
/*     igraph_attribute_type_t type=oldrec->type; */
/*     igraph_vector_t *num=(igraph_vector_t*)oldrec->value; */
/*     igraph_strvector_t *str=(igraph_strvector_t*)oldrec->value; */
/*     switch (type) { */
/*     case IGRAPH_ATTRIBUTE_NUMERIC: */
/*       igraph_vector_permdelete(num, eidx, origlen-newlen); */
/*       break; */
/*     case IGRAPH_ATTRIBUTE_STRING: */
/*       igraph_strvector_permdelete(str, eidx, origlen-newlen); */
/*       break; */
/*     default: */
/*       IGRAPH_WARNING("Unknown edge attribute ignored"); */
/*     } */
/*   } */
/* } */

int igraph_i_cattribute_add_edges(igraph_t *graph, const igraph_vector_t *edges,
                                  igraph_vector_ptr_t *nattr) {

    igraph_i_cattributes_t *attr = graph->attr;
    igraph_vector_ptr_t *eal = &attr->eal;
    long int ealno = igraph_vector_ptr_size(eal);
    long int ne = igraph_vector_size(edges) / 2;
    long int origlen = igraph_ecount(graph) - ne;
    long int nattrno = nattr == 0 ? 0 : igraph_vector_ptr_size(nattr);
    igraph_vector_t news;
    long int newattrs, i;

    /* First add the new attributes if any */
    newattrs = 0;
    IGRAPH_VECTOR_INIT_FINALLY(&news, 0);
    for (i = 0; i < nattrno; i++) {
        igraph_attribute_record_t *nattr_entry = VECTOR(*nattr)[i];
        const char *nname = nattr_entry->name;
        long int j;
        igraph_bool_t l = igraph_i_cattribute_find(eal, nname, &j);
        if (!l) {
            newattrs++;
            IGRAPH_CHECK(igraph_vector_push_back(&news, i));
        } else {
            /* check types */
            if (nattr_entry->type !=
                ((igraph_attribute_record_t*)VECTOR(*eal)[j])->type) {
                IGRAPH_ERROR("You cannot mix attribute types", IGRAPH_EINVAL);
            }
        }
    }

    /* Add NA/empty string vectors for the existing vertices */
    if (newattrs != 0) {
        for (i = 0; i < newattrs; i++) {
            igraph_attribute_record_t *tmp = VECTOR(*nattr)[(long int)VECTOR(news)[i]];
            igraph_attribute_record_t *newrec = igraph_Calloc(1, igraph_attribute_record_t);
            igraph_attribute_type_t type = tmp->type;
            if (!newrec) {
                IGRAPH_ERROR("Cannot add attributes", IGRAPH_ENOMEM);
            }
            IGRAPH_FINALLY(igraph_free, newrec);
            newrec->type = type;
            newrec->name = strdup(tmp->name);
            if (!newrec->name) {
                IGRAPH_ERROR("Cannot add attributes", IGRAPH_ENOMEM);
            }
            IGRAPH_FINALLY(igraph_free, (char*)newrec->name);
            if (type == IGRAPH_ATTRIBUTE_NUMERIC) {
                igraph_vector_t *newnum = igraph_Calloc(1, igraph_vector_t);
                if (!newnum) {
                    IGRAPH_ERROR("Cannot add attributes", IGRAPH_ENOMEM);
                }
                IGRAPH_FINALLY(igraph_free, newnum);
                IGRAPH_VECTOR_INIT_FINALLY(newnum, origlen);
                newrec->value = newnum;
                igraph_vector_fill(newnum, IGRAPH_NAN);
            } else if (type == IGRAPH_ATTRIBUTE_BOOLEAN) {
                igraph_vector_bool_t *newbool = igraph_Calloc(1, igraph_vector_bool_t);
                if (!newbool) {
                    IGRAPH_ERROR("Cannot add attributes", IGRAPH_ENOMEM);
                }
                IGRAPH_FINALLY(igraph_free, newbool);
                IGRAPH_CHECK(igraph_vector_bool_init(newbool, origlen));
                IGRAPH_FINALLY(igraph_vector_bool_destroy, newbool);
                newrec->value = newbool;
                igraph_vector_bool_fill(newbool, 0);
            } else if (type == IGRAPH_ATTRIBUTE_STRING) {
                igraph_strvector_t *newstr = igraph_Calloc(1, igraph_strvector_t);
                if (!newstr) {
                    IGRAPH_ERROR("Cannot add attributes", IGRAPH_ENOMEM);
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
        igraph_attribute_record_t *newrec = 0;
        const char *name = oldrec->name;
        long int j;
        igraph_bool_t l = 0;
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
                IGRAPH_ERROR("Attribute types do not match", IGRAPH_EINVAL);
            }
            switch (oldrec->type) {
            case IGRAPH_ATTRIBUTE_NUMERIC:
                if (ne != igraph_vector_size(newnum)) {
                    IGRAPH_ERROR("Invalid numeric attribute length", IGRAPH_EINVAL);
                }
                IGRAPH_CHECK(igraph_vector_append(oldnum, newnum));
                break;
            case IGRAPH_ATTRIBUTE_STRING:
                if (ne != igraph_strvector_size(newstr)) {
                    IGRAPH_ERROR("Invalid string attribute length", IGRAPH_EINVAL);
                }
                IGRAPH_CHECK(igraph_strvector_append(oldstr, newstr));
                break;
            case IGRAPH_ATTRIBUTE_BOOLEAN:
                if (ne != igraph_vector_bool_size(newbool)) {
                    IGRAPH_ERROR("Invalid Boolean attribute length", IGRAPH_EINVAL);
                }
                IGRAPH_CHECK(igraph_vector_bool_append(oldbool, newbool));
                break;
            default:
                IGRAPH_WARNING("Invalid attribute type");
                break;
            }
        } else {
            /* No such attribute, append NA's */
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

    igraph_vector_destroy(&news);
    IGRAPH_FINALLY_CLEAN(1);

    return 0;
}

/* void igraph_i_cattribute_delete_edges(igraph_t *graph, const igraph_vector_t *idx) { */

/*   igraph_i_cattributes_t *attr=graph->attr; */
/*   igraph_vector_ptr_t *eal=&attr->eal; */
/*   long int ealno=igraph_vector_ptr_size(eal); */
/*   long int i; */
/*   long int origlen=igraph_vector_size(idx), newlen; */

/*   newlen=0; */
/*   for (i=0; i<origlen; i++) { */
/*     if (VECTOR(*idx)[i]>0) { */
/*       newlen++; */
/*     } */
/*   } */
/*   for (i=0; i<ealno; i++) { */
/*     igraph_attribute_record_t *oldrec=VECTOR(*eal)[i]; */
/*     igraph_attribute_type_t type=oldrec->type; */
/*     igraph_vector_t *num=(igraph_vector_t*)oldrec->value; */
/*     igraph_strvector_t *str=(igraph_strvector_t*)oldrec->value; */
/*     switch (type) { */
/*     case IGRAPH_ATTRIBUTE_NUMERIC: */
/*       igraph_vector_permdelete(num, idx, origlen-newlen); */
/*       break; */
/*     case IGRAPH_ATTRIBUTE_STRING: */
/*       igraph_strvector_permdelete(str, idx, origlen-newlen); */
/*       break; */
/*     default: */
/*       IGRAPH_WARNING("Unknown edge attribute ignored"); */
/*     } */
/*   } */

/* } */

int igraph_i_cattribute_permute_edges(const igraph_t *graph,
                                      igraph_t *newgraph,
                                      const igraph_vector_t *idx) {

    if (graph == newgraph) {

        igraph_i_cattributes_t *attr = graph->attr;
        igraph_vector_ptr_t *eal = &attr->eal;
        long int ealno = igraph_vector_ptr_size(eal);
        long int i;

        for (i = 0; i < ealno; i++) {
            igraph_attribute_record_t *oldrec = VECTOR(*eal)[i];
            igraph_attribute_type_t type = oldrec->type;
            igraph_vector_t *num, *newnum;
            igraph_strvector_t *str, *newstr;
            igraph_vector_bool_t *oldbool, *newbool;
            switch (type) {
            case IGRAPH_ATTRIBUTE_NUMERIC:
                num = (igraph_vector_t*) oldrec->value;
                newnum = igraph_Calloc(1, igraph_vector_t);
                if (!newnum) {
                    IGRAPH_ERROR("Cannot permute edge attributes", IGRAPH_ENOMEM);
                }
                IGRAPH_VECTOR_INIT_FINALLY(newnum, 0);
                igraph_vector_index(num, newnum, idx);
                oldrec->value = newnum;
                igraph_vector_destroy(num);
                igraph_Free(num);
                IGRAPH_FINALLY_CLEAN(1);
                break;
            case IGRAPH_ATTRIBUTE_BOOLEAN:
                oldbool = (igraph_vector_bool_t*) oldrec->value;
                newbool = igraph_Calloc(1, igraph_vector_bool_t);
                if (!newbool) {
                    IGRAPH_ERROR("Cannot permute edge attributes", IGRAPH_ENOMEM);
                }
                IGRAPH_CHECK(igraph_vector_bool_init(newbool, 0));
                IGRAPH_FINALLY(igraph_vector_bool_destroy, newbool);
                igraph_vector_bool_index(oldbool, newbool, idx);
                oldrec->value = newbool;
                igraph_vector_bool_destroy(oldbool);
                igraph_Free(oldbool);
                IGRAPH_FINALLY_CLEAN(1);
                break;
            case IGRAPH_ATTRIBUTE_STRING:
                str = (igraph_strvector_t*)oldrec->value;
                newstr = igraph_Calloc(1, igraph_strvector_t);
                if (!newstr) {
                    IGRAPH_ERROR("Cannot permute edge attributes", IGRAPH_ENOMEM);
                }
                IGRAPH_CHECK(igraph_strvector_init(newstr, 0));
                IGRAPH_FINALLY(igraph_strvector_destroy, newstr);
                igraph_strvector_index(str, newstr, idx);
                oldrec->value = newstr;
                igraph_strvector_destroy(str);
                igraph_Free(str);
                IGRAPH_FINALLY_CLEAN(1);
                break;
            default:
                IGRAPH_WARNING("Unknown edge attribute ignored");
            }
        }

    } else {

        igraph_i_cattributes_t *attr = graph->attr;
        igraph_vector_ptr_t *eal = &attr->eal;
        long int ealno = igraph_vector_ptr_size(eal);
        long int i;

        /* New edge attributes */
        igraph_i_cattributes_t *new_attr = newgraph->attr;
        igraph_vector_ptr_t *new_eal = &new_attr->eal;
        IGRAPH_CHECK(igraph_vector_ptr_resize(new_eal, ealno));

        IGRAPH_FINALLY(igraph_i_cattribute_permute_free, new_eal);

        for (i = 0; i < ealno; i++) {
            igraph_attribute_record_t *oldrec = VECTOR(*eal)[i];
            igraph_attribute_type_t type = oldrec->type;
            igraph_vector_t *num, *newnum;
            igraph_strvector_t *str, *newstr;
            igraph_vector_bool_t *oldbool, *newbool;

            /* The record itself */
            igraph_attribute_record_t *new_rec =
                igraph_Calloc(1, igraph_attribute_record_t);
            if (!new_rec) {
                IGRAPH_ERROR("Cannot create edge attributes", IGRAPH_ENOMEM);
            }
            new_rec->name = strdup(oldrec->name);
            new_rec->type = oldrec->type;
            VECTOR(*new_eal)[i] = new_rec;

            switch (type) {
            case IGRAPH_ATTRIBUTE_NUMERIC:
                num = (igraph_vector_t*) oldrec->value;
                newnum = igraph_Calloc(1, igraph_vector_t);
                if (!newnum) {
                    IGRAPH_ERROR("Cannot permute edge attributes", IGRAPH_ENOMEM);
                }
                IGRAPH_VECTOR_INIT_FINALLY(newnum, 0);
                igraph_vector_index(num, newnum, idx);
                new_rec->value = newnum;
                IGRAPH_FINALLY_CLEAN(1);
                break;
            case IGRAPH_ATTRIBUTE_STRING:
                str = (igraph_strvector_t*)oldrec->value;
                newstr = igraph_Calloc(1, igraph_strvector_t);
                if (!newstr) {
                    IGRAPH_ERROR("Cannot permute edge attributes", IGRAPH_ENOMEM);
                }
                IGRAPH_CHECK(igraph_strvector_init(newstr, 0));
                IGRAPH_FINALLY(igraph_strvector_destroy, newstr);
                igraph_strvector_index(str, newstr, idx);
                new_rec->value = newstr;
                IGRAPH_FINALLY_CLEAN(1);
                break;
            case IGRAPH_ATTRIBUTE_BOOLEAN:
                oldbool = (igraph_vector_bool_t*) oldrec->value;
                newbool = igraph_Calloc(1, igraph_vector_bool_t);
                if (!newbool) {
                    IGRAPH_ERROR("Cannot permute edge attributes", IGRAPH_ENOMEM);
                }
                IGRAPH_CHECK(igraph_vector_bool_init(newbool, 0));
                IGRAPH_FINALLY(igraph_vector_bool_destroy, newbool);
                igraph_vector_bool_index(oldbool, newbool, idx);
                new_rec->value = newbool;
                IGRAPH_FINALLY_CLEAN(1);
                break;
            default:
                IGRAPH_WARNING("Unknown edge attribute ignored");
            }
        }
        IGRAPH_FINALLY_CLEAN(1);
    }

    return 0;
}

int igraph_i_cattribute_combine_edges(const igraph_t *graph,
                                      igraph_t *newgraph,
                                      const igraph_vector_ptr_t *merges,
                                      const igraph_attribute_combination_t *comb) {

    igraph_i_cattributes_t *attr = graph->attr;
    igraph_i_cattributes_t *toattr = newgraph->attr;
    igraph_vector_ptr_t *eal = &attr->eal;
    igraph_vector_ptr_t *new_eal = &toattr->eal;
    long int ealno = igraph_vector_ptr_size(eal);
    long int i, j, keepno = 0;
    int *TODO;
    igraph_function_pointer_t *funcs;

    TODO = igraph_Calloc(ealno, int);
    if (!TODO) {
        IGRAPH_ERROR("Cannot combine edge attributes",
                     IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, TODO);
    funcs = igraph_Calloc(ealno, igraph_function_pointer_t);
    if (!funcs) {
        IGRAPH_ERROR("Cannot combine edge attributes",
                     IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, funcs);

    for (i = 0; i < ealno; i++) {
        igraph_attribute_record_t *oldrec = VECTOR(*eal)[i];
        const char *name = oldrec->name;
        igraph_attribute_combination_type_t todo;
        igraph_function_pointer_t voidfunc;
        igraph_attribute_combination_query(comb, name, &todo, &voidfunc);
        TODO[i] = todo;
        funcs[i] = voidfunc;
        if (todo != IGRAPH_ATTRIBUTE_COMBINE_IGNORE) {
            keepno++;
        }
    }

    IGRAPH_CHECK(igraph_vector_ptr_resize(new_eal, keepno));
    IGRAPH_FINALLY(igraph_i_cattribute_permute_free, new_eal);

    for (i = 0, j = 0; i < ealno; i++) {
        igraph_attribute_record_t *newrec, *oldrec = VECTOR(*eal)[i];
        const char *name = oldrec->name;
        igraph_attribute_combination_type_t todo =
            (igraph_attribute_combination_type_t) (TODO[i]);
        igraph_attribute_type_t type = oldrec->type;
        igraph_cattributes_combine_num_t *numfunc =
            (igraph_cattributes_combine_num_t*) funcs[i];
        igraph_cattributes_combine_str_t *strfunc =
            (igraph_cattributes_combine_str_t*) funcs[i];
        igraph_cattributes_combine_bool_t *boolfunc =
            (igraph_cattributes_combine_bool_t*) funcs[i];

        if (todo == IGRAPH_ATTRIBUTE_COMBINE_DEFAULT ||
            todo == IGRAPH_ATTRIBUTE_COMBINE_IGNORE) {
            continue;
        }

        newrec = igraph_Calloc(1, igraph_attribute_record_t);
        if (!newrec) {
            IGRAPH_ERROR("Cannot combine edge attributes",
                         IGRAPH_ENOMEM);
        }
        newrec->name = strdup(name);
        newrec->type = type;
        VECTOR(*new_eal)[j] = newrec;

        if (type == IGRAPH_ATTRIBUTE_NUMERIC) {
            switch (todo) {
            case IGRAPH_ATTRIBUTE_COMBINE_FUNCTION:
                IGRAPH_CHECK(igraph_i_cattributes_cn_func(oldrec, newrec, merges,
                             numfunc));
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
        } else if (type == IGRAPH_ATTRIBUTE_BOOLEAN) {
            switch (todo) {
            case IGRAPH_ATTRIBUTE_COMBINE_FUNCTION:
                IGRAPH_CHECK(igraph_i_cattributes_cb_func(oldrec, newrec, merges,
                             boolfunc));
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
        } else if (type == IGRAPH_ATTRIBUTE_STRING) {
            switch (todo) {
            case IGRAPH_ATTRIBUTE_COMBINE_FUNCTION:
                IGRAPH_CHECK(igraph_i_cattributes_sn_func(oldrec, newrec, merges,
                             strfunc));
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

        j++;
    }

    igraph_free(funcs);
    igraph_free(TODO);
    IGRAPH_FINALLY_CLEAN(2);

    return 0;
}

int igraph_i_cattribute_get_info(const igraph_t *graph,
                                 igraph_strvector_t *gnames,
                                 igraph_vector_t *gtypes,
                                 igraph_strvector_t *vnames,
                                 igraph_vector_t *vtypes,
                                 igraph_strvector_t *enames,
                                 igraph_vector_t *etypes) {

    igraph_strvector_t *names[3] = { gnames, vnames, enames };
    igraph_vector_t *types[3] = { gtypes, vtypes, etypes };
    igraph_i_cattributes_t *at = graph->attr;
    igraph_vector_ptr_t *attr[3] = { &at->gal, &at->val, &at->eal };
    long int i, j;

    for (i = 0; i < 3; i++) {
        igraph_strvector_t *n = names[i];
        igraph_vector_t *t = types[i];
        igraph_vector_ptr_t *al = attr[i];
        long int len = igraph_vector_ptr_size(al);

        if (n) {
            IGRAPH_CHECK(igraph_strvector_resize(n, len));
        }
        if (t) {
            IGRAPH_CHECK(igraph_vector_resize(t, len));
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

    return 0;
}

igraph_bool_t igraph_i_cattribute_has_attr(const igraph_t *graph,
        igraph_attribute_elemtype_t type,
        const char *name) {
    igraph_i_cattributes_t *at = graph->attr;
    igraph_vector_ptr_t *attr[3] = { &at->gal, &at->val, &at->eal };
    long int attrnum;

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

int igraph_i_cattribute_gettype(const igraph_t *graph,
                                igraph_attribute_type_t *type,
                                igraph_attribute_elemtype_t elemtype,
                                const char *name) {
    long int attrnum;
    igraph_attribute_record_t *rec;
    igraph_i_cattributes_t *at = graph->attr;
    igraph_vector_ptr_t *attr[3] = { &at->gal, &at->val, &at->eal };
    igraph_vector_ptr_t *al;
    long int j;
    igraph_bool_t l = 0;

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

    return 0;
}

int igraph_i_cattribute_get_numeric_graph_attr(const igraph_t *graph,
        const char *name,
        igraph_vector_t *value) {
    igraph_i_cattributes_t *attr = graph->attr;
    igraph_vector_ptr_t *gal = &attr->gal;
    long int j;
    igraph_attribute_record_t *rec;
    igraph_vector_t *num;
    igraph_bool_t l = igraph_i_cattribute_find(gal, name, &j);

    if (!l) {
        IGRAPH_ERROR("Unknown attribute", IGRAPH_EINVAL);
    }

    rec = VECTOR(*gal)[j];
    num = (igraph_vector_t*)rec->value;
    IGRAPH_CHECK(igraph_vector_resize(value, 1));
    VECTOR(*value)[0] = VECTOR(*num)[0];

    return 0;
}

int igraph_i_cattribute_get_bool_graph_attr(const igraph_t *graph,
        const char *name,
        igraph_vector_bool_t *value) {
    igraph_i_cattributes_t *attr = graph->attr;
    igraph_vector_ptr_t *gal = &attr->gal;
    long int j;
    igraph_attribute_record_t *rec;
    igraph_vector_bool_t *log;
    igraph_bool_t l = igraph_i_cattribute_find(gal, name, &j);

    if (!l) {
        IGRAPH_ERROR("Unknown attribute", IGRAPH_EINVAL);
    }

    rec = VECTOR(*gal)[j];
    log = (igraph_vector_bool_t*)rec->value;
    IGRAPH_CHECK(igraph_vector_bool_resize(value, 1));
    VECTOR(*value)[0] = VECTOR(*log)[0];

    return 0;
}

int igraph_i_cattribute_get_string_graph_attr(const igraph_t *graph,
        const char *name,
        igraph_strvector_t *value) {
    igraph_i_cattributes_t *attr = graph->attr;
    igraph_vector_ptr_t *gal = &attr->gal;
    long int j;
    igraph_attribute_record_t *rec;
    igraph_strvector_t *str;
    igraph_bool_t l = igraph_i_cattribute_find(gal, name, &j);

    if (!l) {
        IGRAPH_ERROR("Unknown attribute", IGRAPH_EINVAL);
    }

    rec = VECTOR(*gal)[j];
    str = (igraph_strvector_t*)rec->value;
    IGRAPH_CHECK(igraph_strvector_resize(value, 1));
    IGRAPH_CHECK(igraph_strvector_set(value, 0, STR(*str, 0)));

    return 0;
}

int igraph_i_cattribute_get_numeric_vertex_attr(const igraph_t *graph,
        const char *name,
        igraph_vs_t vs,
        igraph_vector_t *value) {
    igraph_i_cattributes_t *attr = graph->attr;
    igraph_vector_ptr_t *val = &attr->val;
    long int j;
    igraph_attribute_record_t *rec;
    igraph_vector_t *num;
    igraph_bool_t l = igraph_i_cattribute_find(val, name, &j);

    if (!l) {
        IGRAPH_ERROR("Unknown attribute", IGRAPH_EINVAL);
    }

    rec = VECTOR(*val)[j];
    num = (igraph_vector_t*)rec->value;
    if (igraph_vs_is_all(&vs)) {
        igraph_vector_clear(value);
        IGRAPH_CHECK(igraph_vector_append(value, num));
    } else {
        igraph_vit_t it;
        long int i = 0;
        IGRAPH_CHECK(igraph_vit_create(graph, vs, &it));
        IGRAPH_FINALLY(igraph_vit_destroy, &it);
        IGRAPH_CHECK(igraph_vector_resize(value, IGRAPH_VIT_SIZE(it)));
        for (; !IGRAPH_VIT_END(it); IGRAPH_VIT_NEXT(it), i++) {
            long int v = IGRAPH_VIT_GET(it);
            VECTOR(*value)[i] = VECTOR(*num)[v];
        }
        igraph_vit_destroy(&it);
        IGRAPH_FINALLY_CLEAN(1);
    }

    return 0;
}

int igraph_i_cattribute_get_bool_vertex_attr(const igraph_t *graph,
        const char *name,
        igraph_vs_t vs,
        igraph_vector_bool_t *value) {
    igraph_i_cattributes_t *attr = graph->attr;
    igraph_vector_ptr_t *val = &attr->val;
    igraph_vit_t it;
    long int i, j, v;
    igraph_attribute_record_t *rec;
    igraph_vector_bool_t *log;
    igraph_bool_t l = igraph_i_cattribute_find(val, name, &j);

    if (!l) {
        IGRAPH_ERROR("Unknown attribute", IGRAPH_EINVAL);
    }

    rec = VECTOR(*val)[j];
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

    return 0;
}

int igraph_i_cattribute_get_string_vertex_attr(const igraph_t *graph,
        const char *name,
        igraph_vs_t vs,
        igraph_strvector_t *value) {
    igraph_i_cattributes_t *attr = graph->attr;
    igraph_vector_ptr_t *val = &attr->val;
    long int j;
    igraph_attribute_record_t *rec;
    igraph_strvector_t *str;
    igraph_bool_t l = igraph_i_cattribute_find(val, name, &j);

    if (!l) {
        IGRAPH_ERROR("Unknown attribute", IGRAPH_EINVAL);
    }

    rec = VECTOR(*val)[j];
    str = (igraph_strvector_t*)rec->value;
    if (igraph_vs_is_all(&vs)) {
        igraph_strvector_resize(value, 0);
        IGRAPH_CHECK(igraph_strvector_append(value, str));
    } else {
        igraph_vit_t it;
        long int i = 0;
        IGRAPH_CHECK(igraph_vit_create(graph, vs, &it));
        IGRAPH_FINALLY(igraph_vit_destroy, &it);
        IGRAPH_CHECK(igraph_strvector_resize(value, IGRAPH_VIT_SIZE(it)));
        for (; !IGRAPH_VIT_END(it); IGRAPH_VIT_NEXT(it), i++) {
            long int v = IGRAPH_VIT_GET(it);
            char *s;
            igraph_strvector_get(str, v, &s);
            IGRAPH_CHECK(igraph_strvector_set(value, i, s));
        }
        igraph_vit_destroy(&it);
        IGRAPH_FINALLY_CLEAN(1);
    }

    return 0;
}

int igraph_i_cattribute_get_numeric_edge_attr(const igraph_t *graph,
        const char *name,
        igraph_es_t es,
        igraph_vector_t *value) {
    igraph_i_cattributes_t *attr = graph->attr;
    igraph_vector_ptr_t *eal = &attr->eal;
    long int j;
    igraph_attribute_record_t *rec;
    igraph_vector_t *num;
    igraph_bool_t l = igraph_i_cattribute_find(eal, name, &j);

    if (!l) {
        IGRAPH_ERROR("Unknown attribute", IGRAPH_EINVAL);
    }

    rec = VECTOR(*eal)[j];
    num = (igraph_vector_t*)rec->value;
    if (igraph_es_is_all(&es)) {
        igraph_vector_clear(value);
        IGRAPH_CHECK(igraph_vector_append(value, num));
    } else {
        igraph_eit_t it;
        long int i = 0;
        IGRAPH_CHECK(igraph_eit_create(graph, es, &it));
        IGRAPH_FINALLY(igraph_eit_destroy, &it);
        IGRAPH_CHECK(igraph_vector_resize(value, IGRAPH_EIT_SIZE(it)));
        for (; !IGRAPH_EIT_END(it); IGRAPH_EIT_NEXT(it), i++) {
            long int e = IGRAPH_EIT_GET(it);
            VECTOR(*value)[i] = VECTOR(*num)[e];
        }
        igraph_eit_destroy(&it);
        IGRAPH_FINALLY_CLEAN(1);
    }

    return 0;
}

int igraph_i_cattribute_get_string_edge_attr(const igraph_t *graph,
        const char *name,
        igraph_es_t es,
        igraph_strvector_t *value) {
    igraph_i_cattributes_t *attr = graph->attr;
    igraph_vector_ptr_t *eal = &attr->eal;
    long int j;
    igraph_attribute_record_t *rec;
    igraph_strvector_t *str;
    igraph_bool_t l = igraph_i_cattribute_find(eal, name, &j);

    if (!l) {
        IGRAPH_ERROR("Unknown attribute", IGRAPH_EINVAL);
    }

    rec = VECTOR(*eal)[j];
    str = (igraph_strvector_t*)rec->value;
    if (igraph_es_is_all(&es)) {
        igraph_strvector_resize(value, 0);
        IGRAPH_CHECK(igraph_strvector_append(value, str));
    } else {
        igraph_eit_t it;
        long int i = 0;
        IGRAPH_CHECK(igraph_eit_create(graph, es, &it));
        IGRAPH_FINALLY(igraph_eit_destroy, &it);
        IGRAPH_CHECK(igraph_strvector_resize(value, IGRAPH_EIT_SIZE(it)));
        for (; !IGRAPH_EIT_END(it); IGRAPH_EIT_NEXT(it), i++) {
            long int e = IGRAPH_EIT_GET(it);
            char *s;
            igraph_strvector_get(str, e, &s);
            IGRAPH_CHECK(igraph_strvector_set(value, i, s));
        }
        igraph_eit_destroy(&it);
        IGRAPH_FINALLY_CLEAN(1);
    }

    return 0;
}

int igraph_i_cattribute_get_bool_edge_attr(const igraph_t *graph,
        const char *name,
        igraph_es_t es,
        igraph_vector_bool_t *value) {
    igraph_i_cattributes_t *attr = graph->attr;
    igraph_vector_ptr_t *eal = &attr->eal;
    long int j;
    igraph_attribute_record_t *rec;
    igraph_vector_bool_t *log;
    igraph_bool_t l = igraph_i_cattribute_find(eal, name, &j);

    if (!l) {
        IGRAPH_ERROR("Unknown attribute", IGRAPH_EINVAL);
    }

    rec = VECTOR(*eal)[j];
    log = (igraph_vector_bool_t*)rec->value;
    if (igraph_es_is_all(&es)) {
        igraph_vector_bool_clear(value);
        IGRAPH_CHECK(igraph_vector_bool_append(value, log));
    } else {
        igraph_eit_t it;
        long int i = 0;
        IGRAPH_CHECK(igraph_eit_create(graph, es, &it));
        IGRAPH_FINALLY(igraph_eit_destroy, &it);
        IGRAPH_CHECK(igraph_vector_bool_resize(value, IGRAPH_EIT_SIZE(it)));
        for (; !IGRAPH_EIT_END(it); IGRAPH_EIT_NEXT(it), i++) {
            long int e = IGRAPH_EIT_GET(it);
            VECTOR(*value)[i] = VECTOR(*log)[e];
        }
        igraph_eit_destroy(&it);
        IGRAPH_FINALLY_CLEAN(1);
    }

    return 0;
}

/* -------------------------------------- */

const igraph_attribute_table_t igraph_cattribute_table = {
    &igraph_i_cattribute_init, &igraph_i_cattribute_destroy,
    &igraph_i_cattribute_copy, &igraph_i_cattribute_add_vertices,
    &igraph_i_cattribute_permute_vertices,
    &igraph_i_cattribute_combine_vertices, &igraph_i_cattribute_add_edges,
    &igraph_i_cattribute_permute_edges,
    &igraph_i_cattribute_combine_edges,
    &igraph_i_cattribute_get_info,
    &igraph_i_cattribute_has_attr, &igraph_i_cattribute_gettype,
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
 * igraph_i_set_attribute_table(&amp;igraph_cattribute_table);
 * </programlisting>
 * </para>
 * <para>Now the attribute functions are available. Please note that
 * the attribute handler must be attached before you call any other
 * igraph functions, otherwise you might end up with graphs without
 * attributes and an active attribute handler, which might cause
 * unexpected program behaviour. The rule is that you attach the
 * attribute handler in the beginning of your
 * <function>main()</function> and never touch it again. (Detaching
 * the attribute handler might lead to memory leaks.)</para>
 *
 * <para>It is not currently possible to have attribute handlers on a
 * per-graph basis. All graphs in an application must be managed with
 * the same attribute handler. (Including the default case when there
 * is no attribute handler at all.</para>
 *
 * <para>The C attribute handler supports attaching real numbers and
 * character strings as attributes. No vectors are allowed, ie. every
 * vertex might have an attribute called <code>name</code>, but it is
 * not possible to have a <code>coords</code> graph (or other)
 * attribute which is a vector of numbers.</para>
 *
 * \example examples/simple/cattributes.c
 * \example examples/simple/cattributes2.c
 * \example examples/simple/cattributes3.c
 * \example examples/simple/cattributes4.c
 */

/**
 * \function igraph_cattribute_GAN
 * Query a numeric graph attribute.
 *
 * Returns the value of the given numeric graph attribute.
 * The attribute must exist, otherwise an error is triggered.
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
    long int j;
    igraph_attribute_record_t *rec;
    igraph_vector_t *num;
    igraph_bool_t l = igraph_i_cattribute_find(gal, name, &j);

    if (!l) {
        igraph_error("Unknown attribute", __FILE__, __LINE__, IGRAPH_EINVAL);
        return 0;
    }

    rec = VECTOR(*gal)[j];
    num = (igraph_vector_t*)rec->value;
    return VECTOR(*num)[0];
}

/**
 * \function igraph_cattribute_GAB
 * Query a boolean graph attribute.
 *
 * Returns the value of the given numeric graph attribute.
 * The attribute must exist, otherwise an error is triggered.
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
    long int j;
    igraph_attribute_record_t *rec;
    igraph_vector_bool_t *log;
    igraph_bool_t l = igraph_i_cattribute_find(gal, name, &j);

    if (!l) {
        igraph_error("Unknown attribute", __FILE__, __LINE__, IGRAPH_EINVAL);
        return 0;
    }

    rec = VECTOR(*gal)[j];
    log = (igraph_vector_bool_t*)rec->value;
    return VECTOR(*log)[0];
}

/**
 * \function igraph_cattribute_GAS
 * Query a string graph attribute.
 *
 * Returns a <type>const</type> pointer to the string graph attribute
 * specified in \p name.
 * The attribute must exist, otherwise an error is triggered.
 * \param graph The input graph.
 * \param name The name of the attribute to query.
 * \return The value of the attribute.
 *
 * \sa \ref GAS for a simpler interface.
 *
 * Time complexity: O(Ag), the number of graph attributes.
 */
const char* igraph_cattribute_GAS(const igraph_t *graph, const char *name) {

    igraph_i_cattributes_t *attr = graph->attr;
    igraph_vector_ptr_t *gal = &attr->gal;
    long int j;
    igraph_attribute_record_t *rec;
    igraph_strvector_t *str;
    igraph_bool_t l = igraph_i_cattribute_find(gal, name, &j);

    if (!l) {
        igraph_error("Unknown attribute", __FILE__, __LINE__, IGRAPH_EINVAL);
        return 0;
    }

    rec = VECTOR(*gal)[j];
    str = (igraph_strvector_t*)rec->value;
    return STR(*str, 0);
}

/**
 * \function igraph_cattribute_VAN
 * Query a numeric vertex attribute.
 *
 * The attribute must exist, otherwise an error is triggered.
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
    long int j;
    igraph_attribute_record_t *rec;
    igraph_vector_t *num;
    igraph_bool_t l = igraph_i_cattribute_find(val, name, &j);

    if (!l) {
        igraph_error("Unknown attribute", __FILE__, __LINE__, IGRAPH_EINVAL);
        return 0;
    }

    rec = VECTOR(*val)[j];
    num = (igraph_vector_t*)rec->value;
    return VECTOR(*num)[(long int)vid];
}

/**
 * \function igraph_cattribute_VAB
 * Query a boolean vertex attribute.
 *
 * The attribute must exist, otherwise an error is triggered.
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
    long int j;
    igraph_attribute_record_t *rec;
    igraph_vector_bool_t *log;
    igraph_bool_t l = igraph_i_cattribute_find(val, name, &j);

    if (!l) {
        igraph_error("Unknown attribute", __FILE__, __LINE__, IGRAPH_EINVAL);
        return 0;
    }

    rec = VECTOR(*val)[j];
    log = (igraph_vector_bool_t*)rec->value;
    return VECTOR(*log)[(long int)vid];
}

/**
 * \function igraph_cattribute_VAS
 * Query a string vertex attribute.
 *
 * The attribute must exist, otherwise an error is triggered.
 * \param graph The input graph.
 * \param name The name of the attribute.
 * \param vid The id of the queried vertex.
 * \return The value of the attribute.
 *
 * \sa The macro \ref VAS for a simpler interface.
 *
 * Time complexity: O(Av), the number of vertex attributes.
 */
const char* igraph_cattribute_VAS(const igraph_t *graph, const char *name,
                                  igraph_integer_t vid) {
    igraph_i_cattributes_t *attr = graph->attr;
    igraph_vector_ptr_t *val = &attr->val;
    long int j;
    igraph_attribute_record_t *rec;
    igraph_strvector_t *str;
    igraph_bool_t l = igraph_i_cattribute_find(val, name, &j);

    if (!l) {
        igraph_error("Unknown attribute", __FILE__, __LINE__, IGRAPH_EINVAL);
        return 0;
    }

    rec = VECTOR(*val)[j];
    str = (igraph_strvector_t*)rec->value;
    return STR(*str, (long int)vid);
}

/**
 * \function igraph_cattribute_EAN
 * Query a numeric edge attribute.
 *
 * The attribute must exist, otherwise an error is triggered.
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
    long int j;
    igraph_attribute_record_t *rec;
    igraph_vector_t *num;
    igraph_bool_t l = igraph_i_cattribute_find(eal, name, &j);

    if (!l) {
        igraph_error("Unknown attribute", __FILE__, __LINE__, IGRAPH_EINVAL);
        return 0;
    }

    rec = VECTOR(*eal)[j];
    num = (igraph_vector_t*)rec->value;
    return VECTOR(*num)[(long int)eid];
}

/**
 * \function igraph_cattribute_EAB
 * Query a boolean edge attribute.
 *
 * The attribute must exist, otherwise an error is triggered.
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
    long int j;
    igraph_attribute_record_t *rec;
    igraph_vector_bool_t *log;
    igraph_bool_t l = igraph_i_cattribute_find(eal, name, &j);

    if (!l) {
        igraph_error("Unknown attribute", __FILE__, __LINE__, IGRAPH_EINVAL);
        return 0;
    }

    rec = VECTOR(*eal)[j];
    log = (igraph_vector_bool_t*)rec->value;
    return VECTOR(*log)[(long int)eid];
}

/**
 * \function igraph_cattribute_EAS
 * Query a string edge attribute.
 *
 * The attribute must exist, otherwise an error is triggered.
 * \param graph The input graph.
 * \param name The name of the attribute.
 * \param eid The id of the queried edge.
 * \return The value of the attribute.
 *
 * \se \ref EAS if you want to type less.
 *
 * Time complexity: O(Ae), the number of edge attributes.
 */
const char* igraph_cattribute_EAS(const igraph_t *graph, const char *name,
                                  igraph_integer_t eid) {
    igraph_i_cattributes_t *attr = graph->attr;
    igraph_vector_ptr_t *eal = &attr->eal;
    long int j;
    igraph_attribute_record_t *rec;
    igraph_strvector_t *str;
    igraph_bool_t l = igraph_i_cattribute_find(eal, name, &j);

    if (!l) {
        igraph_error("Unknown attribute", __FILE__, __LINE__, IGRAPH_EINVAL);
        return 0;
    }

    rec = VECTOR(*eal)[j];
    str = (igraph_strvector_t*)rec->value;
    return STR(*str, (long int)eid);
}

/**
 * \function igraph_cattribute_VANV
 * Query a numeric vertex attribute for many vertices
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

int igraph_cattribute_VANV(const igraph_t *graph, const char *name,
                           igraph_vs_t vids, igraph_vector_t *result) {

    return igraph_i_cattribute_get_numeric_vertex_attr(graph, name, vids,
            result);
}

/**
 * \function igraph_cattribute_VABV
 * Query a boolean vertex attribute for many vertices
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

int igraph_cattribute_VABV(const igraph_t *graph, const char *name,
                           igraph_vs_t vids, igraph_vector_bool_t *result) {

    return igraph_i_cattribute_get_bool_vertex_attr(graph, name, vids,
            result);
}

/**
 * \function igraph_cattribute_EANV
 * Query a numeric edge attribute for many edges
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

int igraph_cattribute_EANV(const igraph_t *graph, const char *name,
                           igraph_es_t eids, igraph_vector_t *result) {

    return igraph_i_cattribute_get_numeric_edge_attr(graph, name, eids,
            result);
}

/**
 * \function igraph_cattribute_EABV
 * Query a boolean edge attribute for many edges
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

int igraph_cattribute_EABV(const igraph_t *graph, const char *name,
                           igraph_es_t eids, igraph_vector_bool_t *result) {

    return igraph_i_cattribute_get_bool_edge_attr(graph, name, eids,
            result);
}

/**
 * \function igraph_cattribute_VASV
 * Query a string vertex attribute for many vertices
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

int igraph_cattribute_VASV(const igraph_t *graph, const char *name,
                           igraph_vs_t vids, igraph_strvector_t *result) {

    return igraph_i_cattribute_get_string_vertex_attr(graph, name, vids,
            result);
}

/**
 * \function igraph_cattribute_EASV
 * Query a string edge attribute for many edges
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

int igraph_cattribute_EASV(const igraph_t *graph, const char *name,
                           igraph_es_t eids, igraph_strvector_t *result) {

    return igraph_i_cattribute_get_string_edge_attr(graph, name, eids,
            result);
}

/**
 * \function igraph_cattribute_list
 * List all attributes
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
int igraph_cattribute_list(const igraph_t *graph,
                           igraph_strvector_t *gnames, igraph_vector_t *gtypes,
                           igraph_strvector_t *vnames, igraph_vector_t *vtypes,
                           igraph_strvector_t *enames, igraph_vector_t *etypes) {
    return igraph_i_cattribute_get_info(graph, gnames, gtypes, vnames, vtypes,
                                        enames, etypes);
}

/**
 * \function igraph_cattribute_has_attr
 * Checks whether a (graph, vertex or edge) attribute exists
 *
 * \param graph The graph.
 * \param type The type of the attribute, \c IGRAPH_ATTRIBUTE_GRAPH,
 *        \c IGRAPH_ATTRIBUTE_VERTEX or \c IGRAPH_ATTRIBUTE_EDGE.
 * \param name Character constant, the name of the attribute.
 * \return Logical value, TRUE if the attribute exists, FALSE otherwise.
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
 * Set a numeric graph attribute
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
int igraph_cattribute_GAN_set(igraph_t *graph, const char *name,
                              igraph_real_t value) {

    igraph_i_cattributes_t *attr = graph->attr;
    igraph_vector_ptr_t *gal = &attr->gal;
    long int j;
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
        igraph_attribute_record_t *rec = igraph_Calloc(1, igraph_attribute_record_t);
        igraph_vector_t *num;
        if (!rec) {
            IGRAPH_ERROR("Cannot add graph attribute", IGRAPH_ENOMEM);
        }
        IGRAPH_FINALLY(igraph_free, rec);
        rec->name = strdup(name);
        if (!rec->name) {
            IGRAPH_ERROR("Cannot add graph attribute", IGRAPH_ENOMEM);
        }
        IGRAPH_FINALLY(igraph_free, (char*)rec->name);
        rec->type = IGRAPH_ATTRIBUTE_NUMERIC;
        num = igraph_Calloc(1, igraph_vector_t);
        if (!num) {
            IGRAPH_ERROR("Cannot add graph attribute", IGRAPH_ENOMEM);
        }
        IGRAPH_FINALLY(igraph_free, num);
        IGRAPH_VECTOR_INIT_FINALLY(num, 1);
        VECTOR(*num)[0] = value;
        rec->value = num;
        IGRAPH_CHECK(igraph_vector_ptr_push_back(gal, rec));
        IGRAPH_FINALLY_CLEAN(4);
    }

    return 0;
}

/**
 * \function igraph_cattribute_GAB_set
 * Set a boolean graph attribute
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
int igraph_cattribute_GAB_set(igraph_t *graph, const char *name,
                              igraph_bool_t value) {

    igraph_i_cattributes_t *attr = graph->attr;
    igraph_vector_ptr_t *gal = &attr->gal;
    long int j;
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
        igraph_attribute_record_t *rec = igraph_Calloc(1, igraph_attribute_record_t);
        igraph_vector_bool_t *log;
        if (!rec) {
            IGRAPH_ERROR("Cannot add graph attribute", IGRAPH_ENOMEM);
        }
        IGRAPH_FINALLY(igraph_free, rec);
        rec->name = strdup(name);
        if (!rec->name) {
            IGRAPH_ERROR("Cannot add graph attribute", IGRAPH_ENOMEM);
        }
        IGRAPH_FINALLY(igraph_free, (char*)rec->name);
        rec->type = IGRAPH_ATTRIBUTE_BOOLEAN;
        log = igraph_Calloc(1, igraph_vector_bool_t);
        if (!log) {
            IGRAPH_ERROR("Cannot add graph attribute", IGRAPH_ENOMEM);
        }
        IGRAPH_FINALLY(igraph_free, log);
        IGRAPH_CHECK(igraph_vector_bool_init(log, 1));
        IGRAPH_FINALLY(igraph_vector_bool_destroy, log);
        VECTOR(*log)[0] = value;
        rec->value = log;
        IGRAPH_CHECK(igraph_vector_ptr_push_back(gal, rec));
        IGRAPH_FINALLY_CLEAN(4);
    }

    return 0;
}

/**
 * \function igraph_cattribute_GAS_set
 * Set a string graph attribute.
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
int igraph_cattribute_GAS_set(igraph_t *graph, const char *name,
                              const char *value) {

    igraph_i_cattributes_t *attr = graph->attr;
    igraph_vector_ptr_t *gal = &attr->gal;
    long int j;
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
        igraph_attribute_record_t *rec = igraph_Calloc(1, igraph_attribute_record_t);
        igraph_strvector_t *str;
        if (!rec) {
            IGRAPH_ERROR("Cannot add graph attribute", IGRAPH_ENOMEM);
        }
        IGRAPH_FINALLY(igraph_free, rec);
        rec->name = strdup(name);
        if (!rec->name) {
            IGRAPH_ERROR("Cannot add graph attribute", IGRAPH_ENOMEM);
        }
        IGRAPH_FINALLY(igraph_free, (char*)rec->name);
        rec->type = IGRAPH_ATTRIBUTE_STRING;
        str = igraph_Calloc(1, igraph_strvector_t);
        if (!str) {
            IGRAPH_ERROR("Cannot add graph attribute", IGRAPH_ENOMEM);
        }
        IGRAPH_FINALLY(igraph_free, str);
        IGRAPH_STRVECTOR_INIT_FINALLY(str, 1);
        IGRAPH_CHECK(igraph_strvector_set(str, 0, value));
        rec->value = str;
        IGRAPH_CHECK(igraph_vector_ptr_push_back(gal, rec));
        IGRAPH_FINALLY_CLEAN(4);
    }

    return 0;
}

/**
 * \function igraph_cattribute_VAN_set
 * Set a numeric vertex attribute
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
int igraph_cattribute_VAN_set(igraph_t *graph, const char *name,
                              igraph_integer_t vid, igraph_real_t value) {

    igraph_i_cattributes_t *attr = graph->attr;
    igraph_vector_ptr_t *val = &attr->val;
    long int j;
    igraph_bool_t l = igraph_i_cattribute_find(val, name, &j);

    if (l) {
        igraph_attribute_record_t *rec = VECTOR(*val)[j];
        if (rec->type != IGRAPH_ATTRIBUTE_NUMERIC) {
            IGRAPH_ERROR("Invalid attribute type", IGRAPH_EINVAL);
        } else {
            igraph_vector_t *num = (igraph_vector_t*)rec->value;
            VECTOR(*num)[(long int)vid] = value;
        }
    } else {
        igraph_attribute_record_t *rec = igraph_Calloc(1, igraph_attribute_record_t);
        igraph_vector_t *num;
        if (!rec) {
            IGRAPH_ERROR("Cannot add vertex attribute", IGRAPH_ENOMEM);
        }
        IGRAPH_FINALLY(igraph_free, rec);
        rec->name = strdup(name);
        if (!rec->name) {
            IGRAPH_ERROR("Cannot add vertex attribute", IGRAPH_ENOMEM);
        }
        IGRAPH_FINALLY(igraph_free, (char*)rec->name);
        rec->type = IGRAPH_ATTRIBUTE_NUMERIC;
        num = igraph_Calloc(1, igraph_vector_t);
        if (!num) {
            IGRAPH_ERROR("Cannot add vertex attribute", IGRAPH_ENOMEM);
        }
        IGRAPH_FINALLY(igraph_free, num);
        IGRAPH_VECTOR_INIT_FINALLY(num, igraph_vcount(graph));
        igraph_vector_fill(num, IGRAPH_NAN);
        VECTOR(*num)[(long int)vid] = value;
        rec->value = num;
        IGRAPH_CHECK(igraph_vector_ptr_push_back(val, rec));
        IGRAPH_FINALLY_CLEAN(4);
    }

    return 0;
}

/**
 * \function igraph_cattribute_VAB_set
 * Set a boolean vertex attribute
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
int igraph_cattribute_VAB_set(igraph_t *graph, const char *name,
                              igraph_integer_t vid, igraph_bool_t value) {

    igraph_i_cattributes_t *attr = graph->attr;
    igraph_vector_ptr_t *val = &attr->val;
    long int j;
    igraph_bool_t l = igraph_i_cattribute_find(val, name, &j);

    if (l) {
        igraph_attribute_record_t *rec = VECTOR(*val)[j];
        if (rec->type != IGRAPH_ATTRIBUTE_BOOLEAN) {
            IGRAPH_ERROR("Invalid attribute type", IGRAPH_EINVAL);
        } else {
            igraph_vector_bool_t *log = (igraph_vector_bool_t*)rec->value;
            VECTOR(*log)[(long int)vid] = value;
        }
    } else {
        igraph_attribute_record_t *rec = igraph_Calloc(1, igraph_attribute_record_t);
        igraph_vector_bool_t *log;
        if (!rec) {
            IGRAPH_ERROR("Cannot add vertex attribute", IGRAPH_ENOMEM);
        }
        IGRAPH_FINALLY(igraph_free, rec);
        rec->name = strdup(name);
        if (!rec->name) {
            IGRAPH_ERROR("Cannot add vertex attribute", IGRAPH_ENOMEM);
        }
        IGRAPH_FINALLY(igraph_free, (char*)rec->name);
        rec->type = IGRAPH_ATTRIBUTE_BOOLEAN;
        log = igraph_Calloc(1, igraph_vector_bool_t);
        if (!log) {
            IGRAPH_ERROR("Cannot add vertex attribute", IGRAPH_ENOMEM);
        }
        IGRAPH_FINALLY(igraph_free, log);
        IGRAPH_CHECK(igraph_vector_bool_init(log, igraph_vcount(graph)));
        IGRAPH_FINALLY(igraph_vector_bool_destroy, log);
        igraph_vector_bool_fill(log, 0);
        VECTOR(*log)[(long int)vid] = value;
        rec->value = log;
        IGRAPH_CHECK(igraph_vector_ptr_push_back(val, rec));
        IGRAPH_FINALLY_CLEAN(4);
    }

    return 0;
}

/**
 * \function igraph_cattribute_VAS_set
 * Set a string vertex attribute
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
int igraph_cattribute_VAS_set(igraph_t *graph, const char *name,
                              igraph_integer_t vid, const char *value) {

    igraph_i_cattributes_t *attr = graph->attr;
    igraph_vector_ptr_t *val = &attr->val;
    long int j;
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
        igraph_attribute_record_t *rec = igraph_Calloc(1, igraph_attribute_record_t);
        igraph_strvector_t *str;
        if (!rec) {
            IGRAPH_ERROR("Cannot add vertex attribute", IGRAPH_ENOMEM);
        }
        IGRAPH_FINALLY(igraph_free, rec);
        rec->name = strdup(name);
        if (!rec->name) {
            IGRAPH_ERROR("Cannot add vertex attribute", IGRAPH_ENOMEM);
        }
        IGRAPH_FINALLY(igraph_free, (char*)rec->name);
        rec->type = IGRAPH_ATTRIBUTE_STRING;
        str = igraph_Calloc(1, igraph_strvector_t);
        if (!str) {
            IGRAPH_ERROR("Cannot add vertex attribute", IGRAPH_ENOMEM);
        }
        IGRAPH_FINALLY(igraph_free, str);
        IGRAPH_STRVECTOR_INIT_FINALLY(str, igraph_vcount(graph));
        IGRAPH_CHECK(igraph_strvector_set(str, vid, value));
        rec->value = str;
        IGRAPH_CHECK(igraph_vector_ptr_push_back(val, rec));
        IGRAPH_FINALLY_CLEAN(4);
    }

    return 0;
}

/**
 * \function igraph_cattribute_EAN_set
 * Set a numeric edge attribute
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
int igraph_cattribute_EAN_set(igraph_t *graph, const char *name,
                              igraph_integer_t eid, igraph_real_t value) {

    igraph_i_cattributes_t *attr = graph->attr;
    igraph_vector_ptr_t *eal = &attr->eal;
    long int j;
    igraph_bool_t l = igraph_i_cattribute_find(eal, name, &j);

    if (l) {
        igraph_attribute_record_t *rec = VECTOR(*eal)[j];
        if (rec->type != IGRAPH_ATTRIBUTE_NUMERIC) {
            IGRAPH_ERROR("Invalid attribute type", IGRAPH_EINVAL);
        } else {
            igraph_vector_t *num = (igraph_vector_t*)rec->value;
            VECTOR(*num)[(long int)eid] = value;
        }
    } else {
        igraph_attribute_record_t *rec = igraph_Calloc(1, igraph_attribute_record_t);
        igraph_vector_t *num;
        if (!rec) {
            IGRAPH_ERROR("Cannot add edge attribute", IGRAPH_ENOMEM);
        }
        IGRAPH_FINALLY(igraph_free, rec);
        rec->name = strdup(name);
        if (!rec->name) {
            IGRAPH_ERROR("Cannot add edge attribute", IGRAPH_ENOMEM);
        }
        IGRAPH_FINALLY(igraph_free, (char*)rec->name);
        rec->type = IGRAPH_ATTRIBUTE_NUMERIC;
        num = igraph_Calloc(1, igraph_vector_t);
        if (!num) {
            IGRAPH_ERROR("Cannot add edge attribute", IGRAPH_ENOMEM);
        }
        IGRAPH_FINALLY(igraph_free, num);
        IGRAPH_VECTOR_INIT_FINALLY(num, igraph_ecount(graph));
        igraph_vector_fill(num, IGRAPH_NAN);
        VECTOR(*num)[(long int)eid] = value;
        rec->value = num;
        IGRAPH_CHECK(igraph_vector_ptr_push_back(eal, rec));
        IGRAPH_FINALLY_CLEAN(4);
    }

    return 0;
}

/**
 * \function igraph_cattribute_EAB_set
 * Set a boolean edge attribute
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
int igraph_cattribute_EAB_set(igraph_t *graph, const char *name,
                              igraph_integer_t eid, igraph_bool_t value) {

    igraph_i_cattributes_t *attr = graph->attr;
    igraph_vector_ptr_t *eal = &attr->eal;
    long int j;
    igraph_bool_t l = igraph_i_cattribute_find(eal, name, &j);

    if (l) {
        igraph_attribute_record_t *rec = VECTOR(*eal)[j];
        if (rec->type != IGRAPH_ATTRIBUTE_BOOLEAN) {
            IGRAPH_ERROR("Invalid attribute type", IGRAPH_EINVAL);
        } else {
            igraph_vector_bool_t *log = (igraph_vector_bool_t*)rec->value;
            VECTOR(*log)[(long int)eid] = value;
        }
    } else {
        igraph_attribute_record_t *rec = igraph_Calloc(1, igraph_attribute_record_t);
        igraph_vector_bool_t *log;
        if (!rec) {
            IGRAPH_ERROR("Cannot add edge attribute", IGRAPH_ENOMEM);
        }
        IGRAPH_FINALLY(igraph_free, rec);
        rec->name = strdup(name);
        if (!rec->name) {
            IGRAPH_ERROR("Cannot add edge attribute", IGRAPH_ENOMEM);
        }
        IGRAPH_FINALLY(igraph_free, (char*)rec->name);
        rec->type = IGRAPH_ATTRIBUTE_BOOLEAN;
        log = igraph_Calloc(1, igraph_vector_bool_t);
        if (!log) {
            IGRAPH_ERROR("Cannot add edge attribute", IGRAPH_ENOMEM);
        }
        IGRAPH_FINALLY(igraph_free, log);
        IGRAPH_CHECK(igraph_vector_bool_init(log, igraph_ecount(graph)));
        IGRAPH_FINALLY(igraph_vector_bool_destroy, log);
        igraph_vector_bool_fill(log, 0);
        VECTOR(*log)[(long int)eid] = value;
        rec->value = log;
        IGRAPH_CHECK(igraph_vector_ptr_push_back(eal, rec));
        IGRAPH_FINALLY_CLEAN(4);
    }

    return 0;
}

/**
 * \function igraph_cattribute_EAS_set
 * Set a string edge attribute
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
int igraph_cattribute_EAS_set(igraph_t *graph, const char *name,
                              igraph_integer_t eid, const char *value) {

    igraph_i_cattributes_t *attr = graph->attr;
    igraph_vector_ptr_t *eal = &attr->eal;
    long int j;
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
        igraph_attribute_record_t *rec = igraph_Calloc(1, igraph_attribute_record_t);
        igraph_strvector_t *str;
        if (!rec) {
            IGRAPH_ERROR("Cannot add edge attribute", IGRAPH_ENOMEM);
        }
        IGRAPH_FINALLY(igraph_free, rec);
        rec->name = strdup(name);
        if (!rec->name) {
            IGRAPH_ERROR("Cannot add edge attribute", IGRAPH_ENOMEM);
        }
        IGRAPH_FINALLY(igraph_free, (char*)rec->name);
        rec->type = IGRAPH_ATTRIBUTE_STRING;
        str = igraph_Calloc(1, igraph_strvector_t);
        if (!str) {
            IGRAPH_ERROR("Cannot add edge attribute", IGRAPH_ENOMEM);
        }
        IGRAPH_FINALLY(igraph_free, str);
        IGRAPH_STRVECTOR_INIT_FINALLY(str, igraph_ecount(graph));
        IGRAPH_CHECK(igraph_strvector_set(str, eid, value));
        rec->value = str;
        IGRAPH_CHECK(igraph_vector_ptr_push_back(eal, rec));
        IGRAPH_FINALLY_CLEAN(4);
    }

    return 0;
}

/**
 * \function igraph_cattribute_VAN_setv
 * Set a numeric vertex attribute for all vertices.
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

int igraph_cattribute_VAN_setv(igraph_t *graph, const char *name,
                               const igraph_vector_t *v) {
    igraph_i_cattributes_t *attr = graph->attr;
    igraph_vector_ptr_t *val = &attr->val;
    long int j;
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
        igraph_attribute_record_t *rec = igraph_Calloc(1, igraph_attribute_record_t);
        igraph_vector_t *num;
        if (!rec) {
            IGRAPH_ERROR("Cannot add vertex attribute", IGRAPH_ENOMEM);
        }
        IGRAPH_FINALLY(igraph_free, rec);
        rec->type = IGRAPH_ATTRIBUTE_NUMERIC;
        rec->name = strdup(name);
        if (!rec->name) {
            IGRAPH_ERROR("Cannot add vertex attribute", IGRAPH_ENOMEM);
        }
        IGRAPH_FINALLY(igraph_free, (char*)rec->name);
        num = igraph_Calloc(1, igraph_vector_t);
        if (!num) {
            IGRAPH_ERROR("Cannot add vertex attribute", IGRAPH_ENOMEM);
        }
        IGRAPH_FINALLY(igraph_free, num);
        rec->value = num;
        IGRAPH_CHECK(igraph_vector_copy(num, v));
        IGRAPH_FINALLY(igraph_vector_destroy, num);
        IGRAPH_CHECK(igraph_vector_ptr_push_back(val, rec));
        IGRAPH_FINALLY_CLEAN(4);
    }

    return 0;
}
/**
 * \function igraph_cattribute_VAB_setv
 * Set a boolean vertex attribute for all vertices.
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

int igraph_cattribute_VAB_setv(igraph_t *graph, const char *name,
                               const igraph_vector_bool_t *v) {
    igraph_i_cattributes_t *attr = graph->attr;
    igraph_vector_ptr_t *val = &attr->val;
    long int j;
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
        igraph_attribute_record_t *rec = igraph_Calloc(1, igraph_attribute_record_t);
        igraph_vector_bool_t *log;
        if (!rec) {
            IGRAPH_ERROR("Cannot add vertex attribute", IGRAPH_ENOMEM);
        }
        IGRAPH_FINALLY(igraph_free, rec);
        rec->type = IGRAPH_ATTRIBUTE_BOOLEAN;
        rec->name = strdup(name);
        if (!rec->name) {
            IGRAPH_ERROR("Cannot add vertex attribute", IGRAPH_ENOMEM);
        }
        IGRAPH_FINALLY(igraph_free, (char*)rec->name);
        log = igraph_Calloc(1, igraph_vector_bool_t);
        if (!log) {
            IGRAPH_ERROR("Cannot add vertex attribute", IGRAPH_ENOMEM);
        }
        IGRAPH_FINALLY(igraph_free, log);
        rec->value = log;
        IGRAPH_CHECK(igraph_vector_bool_copy(log, v));
        IGRAPH_FINALLY(igraph_vector_destroy, log);
        IGRAPH_CHECK(igraph_vector_ptr_push_back(val, rec));
        IGRAPH_FINALLY_CLEAN(4);
    }

    return 0;
}

/**
 * \function igraph_cattribute_VAS_setv
 * Set a string vertex attribute for all vertices.
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
int igraph_cattribute_VAS_setv(igraph_t *graph, const char *name,
                               const igraph_strvector_t *sv) {

    igraph_i_cattributes_t *attr = graph->attr;
    igraph_vector_ptr_t *val = &attr->val;
    long int j;
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
        igraph_attribute_record_t *rec = igraph_Calloc(1, igraph_attribute_record_t);
        igraph_strvector_t *str;
        if (!rec) {
            IGRAPH_ERROR("Cannot add vertex attribute", IGRAPH_ENOMEM);
        }
        IGRAPH_FINALLY(igraph_free, rec);
        rec->type = IGRAPH_ATTRIBUTE_STRING;
        rec->name = strdup(name);
        if (!rec->name) {
            IGRAPH_ERROR("Cannot add vertex attribute", IGRAPH_ENOMEM);
        }
        IGRAPH_FINALLY(igraph_free, (char*)rec->name);
        str = igraph_Calloc(1, igraph_strvector_t);
        if (!str) {
            IGRAPH_ERROR("Cannot add vertex attribute", IGRAPH_ENOMEM);
        }
        IGRAPH_FINALLY(igraph_free, str);
        rec->value = str;
        IGRAPH_CHECK(igraph_strvector_copy(str, sv));
        IGRAPH_FINALLY(igraph_strvector_destroy, str);
        IGRAPH_CHECK(igraph_vector_ptr_push_back(val, rec));
        IGRAPH_FINALLY_CLEAN(4);
    }

    return 0;
}

/**
 * \function igraph_cattribute_EAN_setv
 * Set a numeric edge attribute for all vertices.
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
int igraph_cattribute_EAN_setv(igraph_t *graph, const char *name,
                               const igraph_vector_t *v) {

    igraph_i_cattributes_t *attr = graph->attr;
    igraph_vector_ptr_t *eal = &attr->eal;
    long int j;
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
        igraph_attribute_record_t *rec = igraph_Calloc(1, igraph_attribute_record_t);
        igraph_vector_t *num;
        if (!rec) {
            IGRAPH_ERROR("Cannot add edge attribute", IGRAPH_ENOMEM);
        }
        IGRAPH_FINALLY(igraph_free, rec);
        rec->type = IGRAPH_ATTRIBUTE_NUMERIC;
        rec->name = strdup(name);
        if (!rec->name) {
            IGRAPH_ERROR("Cannot add edge attribute", IGRAPH_ENOMEM);
        }
        IGRAPH_FINALLY(igraph_free, (char*)rec->name);
        num = igraph_Calloc(1, igraph_vector_t);
        if (!num) {
            IGRAPH_ERROR("Cannot add edge attribute", IGRAPH_ENOMEM);
        }
        IGRAPH_FINALLY(igraph_free, num);
        rec->value = num;
        IGRAPH_CHECK(igraph_vector_copy(num, v));
        IGRAPH_FINALLY(igraph_vector_destroy, num);
        IGRAPH_CHECK(igraph_vector_ptr_push_back(eal, rec));
        IGRAPH_FINALLY_CLEAN(4);
    }

    return 0;
}

/**
 * \function igraph_cattribute_EAB_setv
 * Set a boolean edge attribute for all vertices.
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
int igraph_cattribute_EAB_setv(igraph_t *graph, const char *name,
                               const igraph_vector_bool_t *v) {

    igraph_i_cattributes_t *attr = graph->attr;
    igraph_vector_ptr_t *eal = &attr->eal;
    long int j;
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
        igraph_attribute_record_t *rec = igraph_Calloc(1, igraph_attribute_record_t);
        igraph_vector_bool_t *log;
        if (!rec) {
            IGRAPH_ERROR("Cannot add edge attribute", IGRAPH_ENOMEM);
        }
        IGRAPH_FINALLY(igraph_free, rec);
        rec->type = IGRAPH_ATTRIBUTE_BOOLEAN;
        rec->name = strdup(name);
        if (!rec->name) {
            IGRAPH_ERROR("Cannot add edge attribute", IGRAPH_ENOMEM);
        }
        IGRAPH_FINALLY(igraph_free, (char*)rec->name);
        log = igraph_Calloc(1, igraph_vector_bool_t);
        if (!log) {
            IGRAPH_ERROR("Cannot add edge attribute", IGRAPH_ENOMEM);
        }
        IGRAPH_FINALLY(igraph_free, log);
        rec->value = log;
        IGRAPH_CHECK(igraph_vector_bool_copy(log, v));
        IGRAPH_FINALLY(igraph_vector_bool_destroy, log);
        IGRAPH_CHECK(igraph_vector_ptr_push_back(eal, rec));
        IGRAPH_FINALLY_CLEAN(4);
    }

    return 0;
}

/**
 * \function igraph_cattribute_EAS_setv
 * Set a string edge attribute for all vertices.
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
int igraph_cattribute_EAS_setv(igraph_t *graph, const char *name,
                               const igraph_strvector_t *sv) {

    igraph_i_cattributes_t *attr = graph->attr;
    igraph_vector_ptr_t *eal = &attr->eal;
    long int j;
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
        igraph_attribute_record_t *rec = igraph_Calloc(1, igraph_attribute_record_t);
        igraph_strvector_t *str;
        if (!rec) {
            IGRAPH_ERROR("Cannot add vertex attribute", IGRAPH_ENOMEM);
        }
        IGRAPH_FINALLY(igraph_free, rec);
        rec->type = IGRAPH_ATTRIBUTE_STRING;
        rec->name = strdup(name);
        if (!rec->name) {
            IGRAPH_ERROR("Cannot add vertex attribute", IGRAPH_ENOMEM);
        }
        IGRAPH_FINALLY(igraph_free, (char*)rec->name);
        str = igraph_Calloc(1, igraph_strvector_t);
        if (!str) {
            IGRAPH_ERROR("Cannot add vertex attribute", IGRAPH_ENOMEM);
        }
        IGRAPH_FINALLY(igraph_free, str);
        rec->value = str;
        IGRAPH_CHECK(igraph_strvector_copy(str, sv));
        IGRAPH_FINALLY(igraph_strvector_destroy, str);
        IGRAPH_CHECK(igraph_vector_ptr_push_back(eal, rec));
        IGRAPH_FINALLY_CLEAN(4);
    }

    return 0;
}

void igraph_i_cattribute_free_rec(igraph_attribute_record_t *rec) {

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
    igraph_Free(rec->name);
    igraph_Free(rec->value);
    igraph_Free(rec);
}

/**
 * \function igraph_cattribute_remove_g
 * Remove a graph attribute
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
    long int j;
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
 * Remove a vertex attribute
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
    long int j;
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
 * Remove an edge attribute
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
    long int j;
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
 * Remove all graph/vertex/edge attributes
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
        long int i, n = igraph_vector_ptr_size(gal);
        for (i = 0; i < n; i++) {
            igraph_i_cattribute_free_rec(VECTOR(*gal)[i]);
        }
        igraph_vector_ptr_clear(gal);
    }
    if (v) {
        igraph_vector_ptr_t *val = &attr->val;
        long int i, n = igraph_vector_ptr_size(val);
        for (i = 0; i < n; i++) {
            igraph_i_cattribute_free_rec(VECTOR(*val)[i]);
        }
        igraph_vector_ptr_clear(val);
    }
    if (e) {
        igraph_vector_ptr_t *eal = &attr->eal;
        long int i, n = igraph_vector_ptr_size(eal);
        for (i = 0; i < n; i++) {
            igraph_i_cattribute_free_rec(VECTOR(*eal)[i]);
        }
        igraph_vector_ptr_clear(eal);
    }
}
