/* -*- mode: C -*-  */
/*
   IGraph R package.
   Copyright (C) 2006-2012  Gabor Csardi <csardi.gabor@gmail.com>
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

#include "igraph_foreign.h"
#include "igraph_attributes.h"
#include "igraph_interface.h"
#include "igraph_memory.h"

#include "core/math.h"
#include "core/trie.h"
#include "graph/attributes.h"
#include "internal/hacks.h" /* strcasecmp */

#include "config.h"

#include <locale.h>
#include <math.h>    /* isnan */
#include <string.h>
#include <stdarg.h>  /* va_start & co */

#define GRAPHML_NAMESPACE_URI "http://graphml.graphdrawing.org/xmlns"

#if HAVE_LIBXML == 1
#include <libxml/encoding.h>
#include <libxml/parser.h>

xmlEntity blankEntityStruct = {
#ifndef XML_WITHOUT_CORBA
    0,
#endif
    XML_ENTITY_DECL,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    XML_EXTERNAL_GENERAL_PARSED_ENTITY,
    0,
    0,
    0,
    0,
    0,
    1
};

xmlEntityPtr blankEntity = &blankEntityStruct;

#define GRAPHML_PARSE_ERROR_WITH_CODE(state, msg, code) do {  \
        if (state->successful) {                                    \
            igraph_error(msg, IGRAPH_FILE_BASENAME, __LINE__, code);              \
            igraph_i_graphml_sax_handler_error(state, msg);           \
        }                                                           \
    } while (0)
#define GRAPHML_PARSE_ERROR(state, msg) \
    GRAPHML_PARSE_ERROR_WITH_CODE(state, msg, IGRAPH_PARSEERROR)
#define RETURN_GRAPHML_PARSE_ERROR_WITH_CODE(state, msg, code) do {  \
        GRAPHML_PARSE_ERROR_WITH_CODE(state, msg, code);            \
        return;                                                     \
    } while (1)
#define RETURN_GRAPHML_PARSE_ERROR(state, msg) do {           \
        GRAPHML_PARSE_ERROR(state, msg);                            \
        return;                                                     \
    } while (1)

/* TODO: proper error handling */

typedef struct igraph_i_graphml_attribute_record_t {
    const char *id;           /* GraphML id */
    enum { I_GRAPHML_BOOLEAN, I_GRAPHML_INTEGER, I_GRAPHML_LONG,
           I_GRAPHML_FLOAT, I_GRAPHML_DOUBLE, I_GRAPHML_STRING,
           I_GRAPHML_UNKNOWN_TYPE
         } type; /* GraphML type */
    union {
        igraph_real_t as_numeric;
        igraph_bool_t as_boolean;
        char* as_string;
    } default_value;   /* Default value of the attribute, if any */
    igraph_attribute_record_t record;
} igraph_i_graphml_attribute_record_t;

struct igraph_i_graphml_parser_state {
    enum { START, INSIDE_GRAPHML, INSIDE_GRAPH, INSIDE_NODE, INSIDE_EDGE,
           INSIDE_KEY, INSIDE_DEFAULT, INSIDE_DATA, FINISH, UNKNOWN, ERROR
         } st;
    igraph_t *g;
    igraph_trie_t node_trie;
    igraph_strvector_t edgeids;
    igraph_vector_t edgelist;
    igraph_vector_int_t prev_state_stack;
    unsigned int unknown_depth;
    int index;
    igraph_bool_t successful, edges_directed, destroyed;
    igraph_trie_t v_names;
    igraph_vector_ptr_t v_attrs;
    igraph_trie_t e_names;
    igraph_vector_ptr_t e_attrs;
    igraph_trie_t g_names;
    igraph_vector_ptr_t g_attrs;
    igraph_i_graphml_attribute_record_t* current_attr_record;
    xmlChar *data_key;
    igraph_attribute_elemtype_t data_type;
    char *error_message;
    char *data_char;
    long int act_node;
    igraph_bool_t ignore_namespaces;
};

static void igraph_i_report_unhandled_attribute_target(const char* target,
        const char* file, int line) {
    igraph_warningf("Attribute target '%s' is not handled; ignoring corresponding "
                    "attribute specifications", file, line, 0, target);
}

static igraph_real_t igraph_i_graphml_parse_numeric(const char* char_data,
                                                    igraph_real_t default_value) {
    double result;

    if (char_data == 0) {
        return default_value;
    }

    if (sscanf(char_data, "%lf", &result) == 0) {
        return default_value;
    }

    return result;
}

static igraph_bool_t igraph_i_graphml_parse_boolean(const char* char_data,
                                                    igraph_bool_t default_value) {
    int value;
    if (char_data == 0) {
        return default_value;
    }
    if (!strcasecmp("true", char_data)) {
        return 1;
    }
    if (!strcasecmp("yes", char_data)) {
        return 1;
    }
    if (!strcasecmp("false", char_data)) {
        return 0;
    }
    if (!strcasecmp("no", char_data)) {
        return 0;
    }
    if (sscanf(char_data, "%d", &value) == 0) {
        return default_value;
    }
    return value != 0;
}

static void igraph_i_graphml_attribute_record_destroy(igraph_i_graphml_attribute_record_t* rec) {
    if (rec->record.type == IGRAPH_ATTRIBUTE_NUMERIC) {
        if (rec->record.value != 0) {
            igraph_vector_destroy((igraph_vector_t*)rec->record.value);
            IGRAPH_FREE(rec->record.value);
        }
    } else if (rec->record.type == IGRAPH_ATTRIBUTE_STRING) {
        if (rec->record.value != 0) {
            igraph_strvector_destroy((igraph_strvector_t*)rec->record.value);
            if (rec->default_value.as_string != 0) {
                IGRAPH_FREE(rec->default_value.as_string);
            }
            IGRAPH_FREE(rec->record.value);
        }
    } else if (rec->record.type == IGRAPH_ATTRIBUTE_BOOLEAN) {
        if (rec->record.value != 0) {
            igraph_vector_bool_destroy((igraph_vector_bool_t*)rec->record.value);
            IGRAPH_FREE(rec->record.value);
        }
    }
    if (rec->id != 0) {
        IGRAPH_FREE(rec->id);
    }
    if (rec->record.name != 0) {
        IGRAPH_FREE(rec->record.name);
    }
}

static void igraph_i_graphml_destroy_state(struct igraph_i_graphml_parser_state* state) {
    if (state->destroyed) {
        return;
    }
    state->destroyed = 1;

    igraph_trie_destroy(&state->node_trie);
    igraph_strvector_destroy(&state->edgeids);
    igraph_trie_destroy(&state->v_names);
    igraph_trie_destroy(&state->e_names);
    igraph_trie_destroy(&state->g_names);
    igraph_vector_destroy(&state->edgelist);
    igraph_vector_int_destroy(&state->prev_state_stack);

    if (state->error_message) {
        free(state->error_message);
    }
    if (state->data_key) {
        free(state->data_key);
    }
    if (state->data_char) {
        free(state->data_char);
    }

    igraph_vector_ptr_destroy_all(&state->v_attrs);
    igraph_vector_ptr_destroy_all(&state->e_attrs);
    igraph_vector_ptr_destroy_all(&state->g_attrs);

    IGRAPH_FINALLY_CLEAN(1);
}

static void igraph_i_graphml_sax_handler_error(void *state0, const char* msg, ...) {
    struct igraph_i_graphml_parser_state *state =
        (struct igraph_i_graphml_parser_state*)state0;
    va_list ap;

    va_start(ap, msg);

    if (state->error_message == 0) {
        state->error_message = IGRAPH_CALLOC(4096, char);
    }

    state->successful = 0;
    state->st = ERROR;
    vsnprintf(state->error_message, 4096, msg, ap);

    va_end(ap);
}

static xmlEntityPtr igraph_i_graphml_sax_handler_get_entity(void *state0,
        const xmlChar* name) {
    xmlEntityPtr predef = xmlGetPredefinedEntity(name);
    IGRAPH_UNUSED(state0);
    if (predef != NULL) {
        return predef;
    }
    IGRAPH_WARNING("unknown XML entity found\n");
    return blankEntity;
}

static void igraph_i_graphml_handle_unknown_start_tag(struct igraph_i_graphml_parser_state *state) {
    if (state->st != UNKNOWN) {
        igraph_vector_int_push_back(&state->prev_state_stack, state->st);
        state->st = UNKNOWN;
        state->unknown_depth = 1;
    } else {
        state->unknown_depth++;
    }
}

static void igraph_i_graphml_sax_handler_start_document(void *state0) {
    struct igraph_i_graphml_parser_state *state =
        (struct igraph_i_graphml_parser_state*)state0;
    int ret;

    state->st = START;
    state->successful = 1;
    state->edges_directed = 0;
    state->destroyed = 0;
    state->data_key = 0;
    state->error_message = 0;
    state->data_char = 0;
    state->unknown_depth = 0;
    state->ignore_namespaces = 0;

    ret = igraph_vector_int_init(&state->prev_state_stack, 0);
    if (ret) {
        RETURN_GRAPHML_PARSE_ERROR_WITH_CODE(state, "Cannot parse GraphML file", ret);
    }
    ret = igraph_vector_int_reserve(&state->prev_state_stack, 32);
    if (ret) {
        RETURN_GRAPHML_PARSE_ERROR_WITH_CODE(state, "Cannot parse GraphML file", ret);
    }
    IGRAPH_FINALLY(igraph_vector_int_destroy, &state->prev_state_stack);

    ret = igraph_vector_ptr_init(&state->v_attrs, 0);
    if (ret) {
        RETURN_GRAPHML_PARSE_ERROR_WITH_CODE(state, "Cannot parse GraphML file", ret);
    }
    IGRAPH_VECTOR_PTR_SET_ITEM_DESTRUCTOR(&state->v_attrs,
                                          igraph_i_graphml_attribute_record_destroy);
    IGRAPH_FINALLY(igraph_vector_ptr_destroy, &state->v_attrs);

    ret = igraph_vector_ptr_init(&state->e_attrs, 0);
    if (ret) {
        RETURN_GRAPHML_PARSE_ERROR_WITH_CODE(state, "Cannot parse GraphML file", ret);
    }
    IGRAPH_VECTOR_PTR_SET_ITEM_DESTRUCTOR(&state->e_attrs,
                                          igraph_i_graphml_attribute_record_destroy);
    IGRAPH_FINALLY(igraph_vector_ptr_destroy, &state->e_attrs);

    ret = igraph_vector_ptr_init(&state->g_attrs, 0);
    if (ret) {
        RETURN_GRAPHML_PARSE_ERROR_WITH_CODE(state, "Cannot parse GraphML file", ret);
    }
    IGRAPH_VECTOR_PTR_SET_ITEM_DESTRUCTOR(&state->g_attrs,
                                          igraph_i_graphml_attribute_record_destroy);
    IGRAPH_FINALLY(igraph_vector_ptr_destroy, &state->g_attrs);

    ret = igraph_vector_init(&state->edgelist, 0);
    if (ret) {
        RETURN_GRAPHML_PARSE_ERROR_WITH_CODE(state, "Cannot parse GraphML file", ret);
    }
    IGRAPH_FINALLY(igraph_vector_destroy, &state->edgelist);

    ret = igraph_trie_init(&state->node_trie, 1);
    if (ret) {
        RETURN_GRAPHML_PARSE_ERROR_WITH_CODE(state, "Cannot parse GraphML file", ret);
    }
    IGRAPH_FINALLY(igraph_trie_destroy, &state->node_trie);

    ret = igraph_strvector_init(&state->edgeids, 0);
    if (ret) {
        RETURN_GRAPHML_PARSE_ERROR_WITH_CODE(state, "Cannot parse GraphML file", ret);
    }
    IGRAPH_FINALLY(igraph_strvector_destroy, &state->edgeids);

    ret = igraph_trie_init(&state->v_names, 0);
    if (ret) {
        RETURN_GRAPHML_PARSE_ERROR_WITH_CODE(state, "Cannot parse GraphML file", ret);
    }
    IGRAPH_FINALLY(igraph_trie_destroy, &state->v_names);

    ret = igraph_trie_init(&state->e_names, 0);
    if (ret) {
        RETURN_GRAPHML_PARSE_ERROR_WITH_CODE(state, "Cannot parse GraphML file", ret);
    }
    IGRAPH_FINALLY(igraph_trie_destroy, &state->e_names);

    ret = igraph_trie_init(&state->g_names, 0);
    if (ret) {
        RETURN_GRAPHML_PARSE_ERROR_WITH_CODE(state, "Cannot parse GraphML file", ret);
    }
    IGRAPH_FINALLY(igraph_trie_destroy, &state->g_names);

    IGRAPH_FINALLY_CLEAN(10);
    IGRAPH_FINALLY(igraph_i_graphml_destroy_state, state);
}

static void igraph_i_graphml_sax_handler_end_document(void *state0) {
    struct igraph_i_graphml_parser_state *state =
        (struct igraph_i_graphml_parser_state*)state0;
    long i, l;
    int r;
    igraph_attribute_record_t idrec, eidrec;
    const char *idstr = "id";
    igraph_bool_t already_has_vertex_id = 0, already_has_edge_id = 0;

    if (!state->successful) {
        return;
    }

    if (state->index < 0) {

        igraph_vector_ptr_t vattr, eattr, gattr;
        long int esize = igraph_vector_ptr_size(&state->e_attrs);
        const void **tmp;
        r = igraph_vector_ptr_init(&vattr,
                                   igraph_vector_ptr_size(&state->v_attrs) + 1);
        if (r) {
            igraph_error("Cannot parse GraphML file", IGRAPH_FILE_BASENAME, __LINE__, r);
            igraph_i_graphml_sax_handler_error(state, "Cannot parse GraphML file");
            return;
        }
        IGRAPH_FINALLY(igraph_vector_ptr_destroy, &vattr);
        if (igraph_strvector_size(&state->edgeids) != 0) {
            esize++;
        }
        r = igraph_vector_ptr_init(&eattr, esize);
        if (r) {
            igraph_error("Cannot parse GraphML file", IGRAPH_FILE_BASENAME, __LINE__, r);
            igraph_i_graphml_sax_handler_error(state, "Cannot parse GraphML file");
            return;
        }
        IGRAPH_FINALLY(igraph_vector_ptr_destroy, &eattr);
        r = igraph_vector_ptr_init(&gattr, igraph_vector_ptr_size(&state->g_attrs));
        if (r) {
            igraph_error("Cannot parse GraphML file", IGRAPH_FILE_BASENAME, __LINE__, r);
            igraph_i_graphml_sax_handler_error(state, "Cannot parse GraphML file");
            return;
        }
        IGRAPH_FINALLY(igraph_vector_ptr_destroy, &gattr);

        for (i = 0; i < igraph_vector_ptr_size(&state->v_attrs); i++) {
            igraph_i_graphml_attribute_record_t *graphmlrec =
                VECTOR(state->v_attrs)[i];
            igraph_attribute_record_t *rec = &graphmlrec->record;

            /* Check that the name of the vertex attribute is not 'id'.
            If it is then we cannot the complimentary 'id' attribute. */
            if (! strcmp(rec->name, idstr)) {
                already_has_vertex_id = 1;
            }

            if (rec->type == IGRAPH_ATTRIBUTE_NUMERIC) {
                igraph_vector_t *vec = (igraph_vector_t*)rec->value;
                long int origsize = igraph_vector_size(vec);
                long int nodes = igraph_trie_size(&state->node_trie);
                igraph_vector_resize(vec, nodes);
                for (l = origsize; l < nodes; l++) {
                    VECTOR(*vec)[l] = graphmlrec->default_value.as_numeric;
                }
            } else if (rec->type == IGRAPH_ATTRIBUTE_STRING) {
                igraph_strvector_t *strvec = (igraph_strvector_t*)rec->value;
                long int origsize = igraph_strvector_size(strvec);
                long int nodes = igraph_trie_size(&state->node_trie);
                igraph_strvector_resize(strvec, nodes);
                for (l = origsize; l < nodes; l++) {
                    igraph_strvector_set(strvec, l, graphmlrec->default_value.as_string);
                }
            } else if (rec->type == IGRAPH_ATTRIBUTE_BOOLEAN) {
                igraph_vector_bool_t *boolvec = (igraph_vector_bool_t*)rec->value;
                long int origsize = igraph_vector_bool_size(boolvec);
                long int nodes = igraph_trie_size(&state->node_trie);
                igraph_vector_bool_resize(boolvec, nodes);
                for (l = origsize; l < nodes; l++) {
                    VECTOR(*boolvec)[l] = graphmlrec->default_value.as_boolean;
                }
            }
            VECTOR(vattr)[i] = rec;
        }
        if (!already_has_vertex_id) {
            idrec.name = idstr;
            idrec.type = IGRAPH_ATTRIBUTE_STRING;
            tmp = &idrec.value;
            igraph_trie_getkeys(&state->node_trie, (const igraph_strvector_t **)tmp);
            VECTOR(vattr)[i] = &idrec;
        } else {
            igraph_vector_ptr_pop_back(&vattr);
        }

        for (i = 0; i < igraph_vector_ptr_size(&state->e_attrs); i++) {
            igraph_i_graphml_attribute_record_t *graphmlrec =
                VECTOR(state->e_attrs)[i];
            igraph_attribute_record_t *rec = &graphmlrec->record;

            if (! strcmp(rec->name, idstr)) {
                already_has_edge_id = 1;
            }

            if (rec->type == IGRAPH_ATTRIBUTE_NUMERIC) {
                igraph_vector_t *vec = (igraph_vector_t*)rec->value;
                long int origsize = igraph_vector_size(vec);
                long int edges = igraph_vector_size(&state->edgelist) / 2;
                igraph_vector_resize(vec, edges);
                for (l = origsize; l < edges; l++) {
                    VECTOR(*vec)[l] = graphmlrec->default_value.as_numeric;
                }
            } else if (rec->type == IGRAPH_ATTRIBUTE_STRING) {
                igraph_strvector_t *strvec = (igraph_strvector_t*)rec->value;
                long int origsize = igraph_strvector_size(strvec);
                long int edges = igraph_vector_size(&state->edgelist) / 2;
                igraph_strvector_resize(strvec, edges);
                for (l = origsize; l < edges; l++) {
                    igraph_strvector_set(strvec, l, graphmlrec->default_value.as_string);
                }
            } else if (rec->type == IGRAPH_ATTRIBUTE_BOOLEAN) {
                igraph_vector_bool_t *boolvec = (igraph_vector_bool_t*)rec->value;
                long int origsize = igraph_vector_bool_size(boolvec);
                long int edges = igraph_vector_size(&state->edgelist) / 2;
                igraph_vector_bool_resize(boolvec, edges);
                for (l = origsize; l < edges; l++) {
                    VECTOR(*boolvec)[l] = graphmlrec->default_value.as_boolean;
                }
            }
            VECTOR(eattr)[i] = rec;
        }
        if (igraph_strvector_size(&state->edgeids) != 0) {
            if (!already_has_edge_id) {
                long int origsize = igraph_strvector_size(&state->edgeids);
                eidrec.name = idstr;
                eidrec.type = IGRAPH_ATTRIBUTE_STRING;
                igraph_strvector_resize(&state->edgeids,
                                        igraph_vector_size(&state->edgelist) / 2);
                for (; origsize < igraph_strvector_size(&state->edgeids); origsize++) {
                    igraph_strvector_set(&state->edgeids, origsize, "");
                }
                eidrec.value = &state->edgeids;
                VECTOR(eattr)[(long int)igraph_vector_ptr_size(&eattr) - 1] = &eidrec;
            } else {
                igraph_vector_ptr_pop_back(&eattr);
                IGRAPH_WARNING("Could not add edge ids, "
                               "there is already an 'id' edge attribute");
            }
        }

        for (i = 0; i < igraph_vector_ptr_size(&state->g_attrs); i++) {
            igraph_i_graphml_attribute_record_t *graphmlrec =
                VECTOR(state->g_attrs)[i];
            igraph_attribute_record_t *rec = &graphmlrec->record;
            if (rec->type == IGRAPH_ATTRIBUTE_NUMERIC) {
                igraph_vector_t *vec = (igraph_vector_t*)rec->value;
                long int origsize = igraph_vector_size(vec);
                igraph_vector_resize(vec, 1);
                for (l = origsize; l < 1; l++) {
                    VECTOR(*vec)[l] = graphmlrec->default_value.as_numeric;
                }
            } else if (rec->type == IGRAPH_ATTRIBUTE_STRING) {
                igraph_strvector_t *strvec = (igraph_strvector_t*)rec->value;
                long int origsize = igraph_strvector_size(strvec);
                igraph_strvector_resize(strvec, 1);
                for (l = origsize; l < 1; l++) {
                    igraph_strvector_set(strvec, l, graphmlrec->default_value.as_string);
                }
            } else if (rec->type == IGRAPH_ATTRIBUTE_BOOLEAN) {
                igraph_vector_bool_t *boolvec = (igraph_vector_bool_t*)rec->value;
                long int origsize = igraph_vector_bool_size(boolvec);
                igraph_vector_bool_resize(boolvec, 1);
                for (l = origsize; l < 1; l++) {
                    VECTOR(*boolvec)[l] = graphmlrec->default_value.as_boolean;
                }
            }
            VECTOR(gattr)[i] = rec;
        }

        igraph_empty_attrs(state->g, 0, state->edges_directed, &gattr);
        igraph_add_vertices(state->g, (igraph_integer_t)
                            igraph_trie_size(&state->node_trie), &vattr);
        igraph_add_edges(state->g, &state->edgelist, &eattr);

        igraph_vector_ptr_destroy(&vattr);
        igraph_vector_ptr_destroy(&eattr);
        igraph_vector_ptr_destroy(&gattr);
        IGRAPH_FINALLY_CLEAN(3);
    }

    igraph_i_graphml_destroy_state(state);
}

#define toXmlChar(a)   (BAD_CAST(a))
#define fromXmlChar(a) ((char *)(a)) /* not the most elegant way... */

#define XML_ATTR_LOCALNAME(it) (*(it))
#define XML_ATTR_PREFIX(it) (*(it+1))
#define XML_ATTR_URI(it) (*(it+2))
#define XML_ATTR_VALUE_START(it) (*(it+3))
#define XML_ATTR_VALUE_END(it) (*(it+4))
#define XML_ATTR_VALUE(it) *(it+3), (*(it+4))-(*(it+3))

static igraph_i_graphml_attribute_record_t* igraph_i_graphml_add_attribute_key(
        const xmlChar** attrs, int nb_attrs,
        struct igraph_i_graphml_parser_state *state) {
    xmlChar **it;
    xmlChar *localname;
    igraph_trie_t *trie = 0;
    igraph_vector_ptr_t *ptrvector = 0;
    long int id;
    unsigned short int skip = 0;
    int i, ret;
    igraph_i_graphml_attribute_record_t *rec;

    if (!state->successful) {
        return 0;
    }

    rec = IGRAPH_CALLOC(1, igraph_i_graphml_attribute_record_t);
    if (rec == 0) {
        GRAPHML_PARSE_ERROR_WITH_CODE(state, "Cannot parse GraphML file", IGRAPH_ENOMEM);
        return 0;
    }
    IGRAPH_FINALLY(igraph_free, rec);

    rec->type = I_GRAPHML_UNKNOWN_TYPE;

    for (i = 0, it = (xmlChar**)attrs; i < nb_attrs; i++, it += 5) {
        if (XML_ATTR_URI(it) != 0 &&
            !xmlStrEqual(toXmlChar(GRAPHML_NAMESPACE_URI), XML_ATTR_URI(it))) {
            continue;
        }

        localname = XML_ATTR_LOCALNAME(it);

        if (xmlStrEqual(localname, toXmlChar("id"))) {
            rec->id = fromXmlChar(xmlStrndup(XML_ATTR_VALUE(it)));
        } else if (xmlStrEqual(localname, toXmlChar("attr.name"))) {
            rec->record.name = fromXmlChar(xmlStrndup(XML_ATTR_VALUE(it)));
        } else if (xmlStrEqual(localname, toXmlChar("attr.type"))) {
            if (!xmlStrncmp(toXmlChar("boolean"), XML_ATTR_VALUE(it))) {
                rec->type = I_GRAPHML_BOOLEAN;
                rec->record.type = IGRAPH_ATTRIBUTE_BOOLEAN;
                rec->default_value.as_boolean = 0;
            } else if (!xmlStrncmp(toXmlChar("string"), XML_ATTR_VALUE(it))) {
                rec->type = I_GRAPHML_STRING;
                rec->record.type = IGRAPH_ATTRIBUTE_STRING;
                rec->default_value.as_string = strdup("");
            } else if (!xmlStrncmp(toXmlChar("float"), XML_ATTR_VALUE(it))) {
                rec->type = I_GRAPHML_FLOAT;
                rec->record.type = IGRAPH_ATTRIBUTE_NUMERIC;
                rec->default_value.as_numeric = IGRAPH_NAN;
            } else if (!xmlStrncmp(toXmlChar("double"), XML_ATTR_VALUE(it))) {
                rec->type = I_GRAPHML_DOUBLE;
                rec->record.type = IGRAPH_ATTRIBUTE_NUMERIC;
                rec->default_value.as_numeric = IGRAPH_NAN;
            } else if (!xmlStrncmp(toXmlChar("int"), XML_ATTR_VALUE(it))) {
                rec->type = I_GRAPHML_INTEGER;
                rec->record.type = IGRAPH_ATTRIBUTE_NUMERIC;
                rec->default_value.as_numeric = IGRAPH_NAN;
            } else if (!xmlStrncmp(toXmlChar("long"), XML_ATTR_VALUE(it))) {
                rec->type = I_GRAPHML_LONG;
                rec->record.type = IGRAPH_ATTRIBUTE_NUMERIC;
                rec->default_value.as_numeric = IGRAPH_NAN;
            } else {
                GRAPHML_PARSE_ERROR(state,
                                    "Cannot parse GraphML file, unknown attribute type");
                return 0;
            }
        } else if (xmlStrEqual(*it, toXmlChar("for"))) {
            /* graph, vertex or edge attribute? */
            if (!xmlStrncmp(toXmlChar("graph"), XML_ATTR_VALUE(it))) {
                trie = &state->g_names;
                ptrvector = &state->g_attrs;
            } else if (!xmlStrncmp(toXmlChar("node"), XML_ATTR_VALUE(it))) {
                trie = &state->v_names;
                ptrvector = &state->v_attrs;
            } else if (!xmlStrncmp(toXmlChar("edge"), XML_ATTR_VALUE(it))) {
                trie = &state->e_names;
                ptrvector = &state->e_attrs;
            } else if (!xmlStrncmp(toXmlChar("graphml"), XML_ATTR_VALUE(it))) {
                igraph_i_report_unhandled_attribute_target("graphml", IGRAPH_FILE_BASENAME, __LINE__);
                skip = 1;
            } else if (!xmlStrncmp(toXmlChar("hyperedge"), XML_ATTR_VALUE(it))) {
                igraph_i_report_unhandled_attribute_target("hyperedge", IGRAPH_FILE_BASENAME, __LINE__);
                skip = 1;
            } else if (!xmlStrncmp(toXmlChar("port"), XML_ATTR_VALUE(it))) {
                igraph_i_report_unhandled_attribute_target("port", IGRAPH_FILE_BASENAME, __LINE__);
                skip = 1;
            } else if (!xmlStrncmp(toXmlChar("endpoint"), XML_ATTR_VALUE(it))) {
                igraph_i_report_unhandled_attribute_target("endpoint", IGRAPH_FILE_BASENAME, __LINE__);
                skip = 1;
            } else if (!xmlStrncmp(toXmlChar("all"), XML_ATTR_VALUE(it))) {
                /* TODO: we should handle this */
                igraph_i_report_unhandled_attribute_target("all", IGRAPH_FILE_BASENAME, __LINE__);
                skip = 1;
            } else {
                GRAPHML_PARSE_ERROR(state,
                                    "Cannot parse GraphML file, unknown value in the 'for' attribute of a <key> tag");
                return 0;
            }
        }
    }

    /* throw an error if there is no ID; this is a clear violation of the GraphML
     * DTD */
    if (rec->id == 0) {
        GRAPHML_PARSE_ERROR(state, "Found <key> tag with no 'id' attribute");
        return 0;
    }

    /* in case of a missing attr.name attribute, use the id as the attribute name */
    if (rec->record.name == 0) {
        rec->record.name = strdup(rec->id);
    }

    /* if the attribute type is missing, throw an error */
    if (!skip && rec->type == I_GRAPHML_UNKNOWN_TYPE) {
        igraph_warningf("Ignoring <key id=\"%s\"> because of a missing or unknown 'attr.type' attribute", IGRAPH_FILE_BASENAME, __LINE__, 0, rec->id);
        skip = 1;
    }

    /* if the value of the 'for' attribute was unknown, throw an error */
    if (!skip && trie == 0) {
        GRAPHML_PARSE_ERROR(state,
                            "Cannot parse GraphML file, missing 'for' attribute in a <key> tag");
        return 0;
    }

    /* if the code above requested skipping the attribute, free everything and
     * return */
    if (skip) {
        igraph_free(rec);
        IGRAPH_FINALLY_CLEAN(1);
        return 0;
    }

    /* add to trie, attribues */
    igraph_trie_get(trie, rec->id, &id);
    if (id != igraph_trie_size(trie) - 1) {
        GRAPHML_PARSE_ERROR(state, "Cannot parse GraphML file, duplicate attribute");
        return 0;
    }

    ret = igraph_vector_ptr_push_back(ptrvector, rec);
    if (ret) {
        GRAPHML_PARSE_ERROR_WITH_CODE(state, "Cannot read GraphML file", ret);
        return 0;
    }

    /* Ownership of 'rec' is now taken by ptrvector so we can clean the
     * finally stack */
    IGRAPH_FINALLY_CLEAN(1);  /* rec */

    /* create the attribute values */
    switch (rec->record.type) {
        igraph_vector_t *vec;
        igraph_vector_bool_t *boolvec;
        igraph_strvector_t *strvec;
    case IGRAPH_ATTRIBUTE_BOOLEAN:
        boolvec = IGRAPH_CALLOC(1, igraph_vector_bool_t);
        if (boolvec == 0) {
            GRAPHML_PARSE_ERROR_WITH_CODE(state, "Cannot parse GraphML file", IGRAPH_ENOMEM);
            return 0;
        }
        rec->record.value = boolvec;
        igraph_vector_bool_init(boolvec, 0);
        break;
    case IGRAPH_ATTRIBUTE_NUMERIC:
        vec = IGRAPH_CALLOC(1, igraph_vector_t);
        if (vec == 0) {
            GRAPHML_PARSE_ERROR_WITH_CODE(state, "Cannot parse GraphML file", IGRAPH_ENOMEM);
            return 0;
        }
        rec->record.value = vec;
        igraph_vector_init(vec, 0);
        break;
    case IGRAPH_ATTRIBUTE_STRING:
        strvec = IGRAPH_CALLOC(1, igraph_strvector_t);
        if (strvec == 0) {
            GRAPHML_PARSE_ERROR_WITH_CODE(state, "Cannot parse GraphML file", IGRAPH_ENOMEM);
            return 0;
        }
        rec->record.value = strvec;
        igraph_strvector_init(strvec, 0);
        break;
    default: break;
    }

    return rec;
}

static void igraph_i_graphml_attribute_data_setup(struct igraph_i_graphml_parser_state *state,
                                                  const xmlChar **attrs,
                                                  int nb_attrs,
                                                  igraph_attribute_elemtype_t type) {
    xmlChar **it;
    int i;

    if (!state->successful) {
        return;
    }

    for (i = 0, it = (xmlChar**)attrs; i < nb_attrs; i++, it += 5) {
        if (XML_ATTR_URI(it) != 0 &&
            !xmlStrEqual(toXmlChar(GRAPHML_NAMESPACE_URI), XML_ATTR_URI(it))) {
            continue;
        }

        if (xmlStrEqual(*it, toXmlChar("key"))) {
            if (state->data_key) {
                free(state->data_key);
            }
            state->data_key = xmlStrndup(XML_ATTR_VALUE(it));
            if (state->data_char) {
                free(state->data_char);
            }
            state->data_char = 0;
            state->data_type = type;
        } else {
            /* ignore */
        }
    }
}

static void igraph_i_graphml_append_to_data_char(struct igraph_i_graphml_parser_state *state,
                                                 const xmlChar *data, int len) {
    long int data_char_new_start = 0;

    if (!state->successful) {
        return;
    }

    if (state->data_char) {
        data_char_new_start = (long int) strlen(state->data_char);
        state->data_char = IGRAPH_REALLOC(state->data_char,
                                          (size_t)(data_char_new_start + len + 1), char);
    } else {
        state->data_char = IGRAPH_CALLOC((size_t) len + 1, char);
    }
    if (state->data_char == 0) {
        RETURN_GRAPHML_PARSE_ERROR_WITH_CODE(state, "Cannot parse GraphML file", IGRAPH_ENOMEM);
    }
    memcpy(state->data_char + data_char_new_start, data,
           (size_t) len * sizeof(xmlChar));
    state->data_char[data_char_new_start + len] = '\0';
}

static void igraph_i_graphml_attribute_data_finish(struct igraph_i_graphml_parser_state *state) {
    const char *key = fromXmlChar(state->data_key);
    igraph_attribute_elemtype_t type = state->data_type;
    igraph_trie_t *trie = 0;
    igraph_vector_ptr_t *ptrvector = 0;
    igraph_i_graphml_attribute_record_t *graphmlrec;
    igraph_attribute_record_t *rec;
    long int recid, id = 0;
    int ret;

    switch (type) {
    case IGRAPH_ATTRIBUTE_GRAPH:
        trie = &state->g_names;
        ptrvector = &state->g_attrs;
        id = 0;
        break;
    case IGRAPH_ATTRIBUTE_VERTEX:
        trie = &state->v_names;
        ptrvector = &state->v_attrs;
        id = state->act_node;
        break;
    case IGRAPH_ATTRIBUTE_EDGE:
        trie = &state->e_names;
        ptrvector = &state->e_attrs;
        id = igraph_vector_size(&state->edgelist) / 2 - 1; /* hack */
        break;
    default:
        /* impossible */
        break;
    }

    if (key == 0) {
        /* no key specified, issue a warning */
        igraph_warningf(
            "missing attribute key in a <data> tag, ignoring attribute",
            IGRAPH_FILE_BASENAME, __LINE__, 0,
            key
        );
        IGRAPH_FREE(state->data_char);
        return;
    }

    igraph_trie_check(trie, key, &recid);
    if (recid < 0) {
        /* no such attribute key, issue a warning */
        igraph_warningf(
            "unknown attribute key '%s' in a <data> tag, ignoring attribute",
            IGRAPH_FILE_BASENAME, __LINE__, 0,
            key
        );
        IGRAPH_FREE(state->data_char);
        return;
    }

    graphmlrec = VECTOR(*ptrvector)[recid];
    rec = &graphmlrec->record;

    switch (rec->type) {
        igraph_vector_bool_t *boolvec;
        igraph_vector_t *vec;
        igraph_strvector_t *strvec;
        long int s, i;
        const char* strvalue;
    case IGRAPH_ATTRIBUTE_BOOLEAN:
        boolvec = (igraph_vector_bool_t *)rec->value;
        s = igraph_vector_bool_size(boolvec);
        if (id >= s) {
            ret = igraph_vector_bool_resize(boolvec, id + 1);
            if (ret) {
                RETURN_GRAPHML_PARSE_ERROR_WITH_CODE(state, "Cannot parse GraphML file", ret);
            }
            for (i = s; i < id; i++) {
                VECTOR(*boolvec)[i] = graphmlrec->default_value.as_boolean;
            }
        }
        VECTOR(*boolvec)[id] = igraph_i_graphml_parse_boolean(state->data_char,
                               graphmlrec->default_value.as_boolean);
        break;
    case IGRAPH_ATTRIBUTE_NUMERIC:
        vec = (igraph_vector_t *)rec->value;
        s = igraph_vector_size(vec);
        if (id >= s) {
            ret = igraph_vector_resize(vec, id + 1);
            if (ret) {
                RETURN_GRAPHML_PARSE_ERROR_WITH_CODE(state, "Cannot parse GraphML file", ret);
            }
            for (i = s; i < id; i++) {
                VECTOR(*vec)[i] = graphmlrec->default_value.as_numeric;
            }
        }
        VECTOR(*vec)[id] = igraph_i_graphml_parse_numeric(state->data_char,
                           graphmlrec->default_value.as_numeric);
        break;
    case IGRAPH_ATTRIBUTE_STRING:
        strvec = (igraph_strvector_t *)rec->value;
        s = igraph_strvector_size(strvec);
        if (id >= s) {
            ret = igraph_strvector_resize(strvec, id + 1);
            if (ret) {
                RETURN_GRAPHML_PARSE_ERROR_WITH_CODE(state, "Cannot parse GraphML file", ret);
            }
            strvalue = graphmlrec->default_value.as_string;
            for (i = s; i < id; i++) {
                igraph_strvector_set(strvec, i, strvalue);
            }
        }
        if (state->data_char) {
            strvalue = state->data_char;
        } else {
            strvalue = graphmlrec->default_value.as_string;
        }
        ret = igraph_strvector_set(strvec, id, strvalue);
        if (ret) {
            RETURN_GRAPHML_PARSE_ERROR_WITH_CODE(state, "Cannot parse GraphML file", ret);
        }
        break;
    default:
        break;
    }

    if (state->data_char) {
        IGRAPH_FREE(state->data_char);
    }
}

static void igraph_i_graphml_attribute_default_value_finish(
        struct igraph_i_graphml_parser_state *state) {
    igraph_i_graphml_attribute_record_t *graphmlrec = state->current_attr_record;

    if (graphmlrec == 0) {
        igraph_warning("state->current_attr_record was null where it should have been "
                       "non-null; this is probably a bug. Please notify the developers!",
                       IGRAPH_FILE_BASENAME, __LINE__, 0);
        return;
    }

    if (state->data_char == 0) {
        return;
    }

    switch (graphmlrec->record.type) {
    case IGRAPH_ATTRIBUTE_BOOLEAN:
        graphmlrec->default_value.as_boolean = igraph_i_graphml_parse_boolean(
                state->data_char, 0);
        break;
    case IGRAPH_ATTRIBUTE_NUMERIC:
        graphmlrec->default_value.as_numeric = igraph_i_graphml_parse_numeric(
                state->data_char, IGRAPH_NAN);
        break;
    case IGRAPH_ATTRIBUTE_STRING:
        if (state->data_char) {
            if (graphmlrec->default_value.as_string != 0) {
                free(graphmlrec->default_value.as_string);
            }
            graphmlrec->default_value.as_string = strdup(state->data_char);
        }
        break;
    default:
        break;
    }

    if (state->data_char) {
        IGRAPH_FREE(state->data_char);
    }
}

static void igraph_i_graphml_sax_handler_start_element_ns(
        void *state0, const xmlChar* localname, const xmlChar* prefix,
        const xmlChar* uri, int nb_namespaces, const xmlChar** namespaces,
        int nb_attributes, int nb_defaulted, const xmlChar** attributes) {
    struct igraph_i_graphml_parser_state *state =
        (struct igraph_i_graphml_parser_state*)state0;
    xmlChar** it;
    char* attr_value;
    long int id1, id2;
    int i;
    igraph_bool_t tag_is_unknown = 0;

    IGRAPH_UNUSED(prefix);
    IGRAPH_UNUSED(nb_namespaces);
    IGRAPH_UNUSED(namespaces);
    IGRAPH_UNUSED(nb_defaulted);

    if (!state->successful) {
        return;
    }

    if (uri) {
        if (!xmlStrEqual(toXmlChar(GRAPHML_NAMESPACE_URI), uri)) {
            /* Tag is in a different namespace, so treat it as an unknown start
             * tag irrespectively of our state */
            tag_is_unknown = 1;
        }
    } else {
        /* No namespace URI. If we are in lenient mode, accept it and proceed
         * as if we are in the GraphML namespace to handle lots of naive
         * non-namespace-aware GraphML files floating out there. If we are not
         * in lenient mode _but_ we are in the START state, accept it as well
         * and see whether the root tag is <graphml> (in which case we will
         * enter lenient mode). Otherwise, reject the tag */
        if (!state->ignore_namespaces && state->st != START) {
            tag_is_unknown = 1;
        }
    }

    if (tag_is_unknown) {
        igraph_i_graphml_handle_unknown_start_tag(state);
        return;
    }

    switch (state->st) {
    case START:
        /* If we are in the START state and received a graphml tag,
         * change to INSIDE_GRAPHML state. Otherwise, change to UNKNOWN. */
        if (xmlStrEqual(localname, toXmlChar("graphml"))) {
            if (uri == 0) {
                state->ignore_namespaces = 1;
            }
            state->st = INSIDE_GRAPHML;
        } else {
            igraph_i_graphml_handle_unknown_start_tag(state);
        }
        break;

    case INSIDE_GRAPHML:
        /* If we are in the INSIDE_GRAPHML state and received a graph tag,
         * change to INSIDE_GRAPH state if the state->index counter reached
         * zero (this is to handle multiple graphs in the same file).
         * Otherwise, change to UNKNOWN. */
        if (xmlStrEqual(localname, toXmlChar("graph"))) {
            if (state->index == 0) {
                state->st = INSIDE_GRAPH;
                for (i = 0, it = (xmlChar**)attributes; i < nb_attributes; i++, it += 5) {
                    if (XML_ATTR_URI(it) != 0 &&
                        !xmlStrEqual(toXmlChar(GRAPHML_NAMESPACE_URI), XML_ATTR_URI(it))) {
                        /* Attribute is from a different namespace, so skip it */
                        continue;
                    }
                    if (xmlStrEqual(*it, toXmlChar("edgedefault"))) {
                        if (!xmlStrncmp(toXmlChar("directed"), XML_ATTR_VALUE(it))) {
                            state->edges_directed = 1;
                        } else if (!xmlStrncmp(toXmlChar("undirected"), XML_ATTR_VALUE(it))) {
                            state->edges_directed = 0;
                        }
                    }
                }
            }
            state->index--;
        } else if (xmlStrEqual(localname, toXmlChar("key"))) {
            state->current_attr_record =
                igraph_i_graphml_add_attribute_key(attributes, nb_attributes, state);
            state->st = INSIDE_KEY;
        } else {
            igraph_i_graphml_handle_unknown_start_tag(state);
        }
        break;

    case INSIDE_KEY:
        /* If we are in the INSIDE_KEY state, check for default tag */
        if (xmlStrEqual(localname, toXmlChar("default"))) {
            state->st = INSIDE_DEFAULT;
        } else {
            igraph_i_graphml_handle_unknown_start_tag(state);
        }
        break;

    case INSIDE_DEFAULT:
        /* If we are in the INSIDE_DEFAULT state, every further tag will be unknown */
        igraph_i_graphml_handle_unknown_start_tag(state);
        break;

    case INSIDE_GRAPH:
        /* If we are in the INSIDE_GRAPH state, check for node and edge tags */
        if (xmlStrEqual(localname, toXmlChar("edge"))) {
            id1 = -1; id2 = -1;
            for (i = 0, it = (xmlChar**)attributes; i < nb_attributes; i++, it += 5) {
                if (XML_ATTR_URI(it) != 0 &&
                    !xmlStrEqual(toXmlChar(GRAPHML_NAMESPACE_URI), XML_ATTR_URI(it))) {
                    /* Attribute is from a different namespace, so skip it */
                    continue;
                }
                if (xmlStrEqual(*it, toXmlChar("source"))) {
                    attr_value = fromXmlChar(xmlStrndup(XML_ATTR_VALUE(it)));
                    igraph_trie_get(&state->node_trie, attr_value, &id1);
                    free(attr_value);
                } else if (xmlStrEqual(*it, toXmlChar("target"))) {
                    attr_value = fromXmlChar(xmlStrndup(XML_ATTR_VALUE(it)));
                    igraph_trie_get(&state->node_trie, attr_value, &id2);
                    free(attr_value);
                } else if (xmlStrEqual(*it, toXmlChar("id"))) {
                    long int edges = igraph_vector_size(&state->edgelist) / 2 + 1;
                    long int origsize = igraph_strvector_size(&state->edgeids);
                    attr_value = fromXmlChar(xmlStrndup(XML_ATTR_VALUE(it)));
                    igraph_strvector_resize(&state->edgeids, edges);
                    for (; origsize < edges - 1; origsize++) {
                        igraph_strvector_set(&state->edgeids, origsize, "");
                    }
                    igraph_strvector_set(&state->edgeids, edges - 1, attr_value);
                    free(attr_value);
                }
            }
            if (id1 >= 0 && id2 >= 0) {
                igraph_vector_push_back(&state->edgelist, id1);
                igraph_vector_push_back(&state->edgelist, id2);
            } else {
                igraph_i_graphml_sax_handler_error(state, "Edge with missing source or target encountered");
                return;
            }
            state->st = INSIDE_EDGE;
        } else if (xmlStrEqual(localname, toXmlChar("node"))) {
            id1 = -1;
            for (i = 0, it = (xmlChar**)attributes; i < nb_attributes; i++, it += 5) {
                if (XML_ATTR_URI(it) != 0 &&
                    !xmlStrEqual(toXmlChar(GRAPHML_NAMESPACE_URI), XML_ATTR_URI(it))) {
                    /* Attribute is from a different namespace, so skip it */
                    continue;
                }
                if (xmlStrEqual(XML_ATTR_LOCALNAME(it), toXmlChar("id"))) {
                    attr_value = fromXmlChar(xmlStrndup(XML_ATTR_VALUE(it)));
                    igraph_trie_get(&state->node_trie, attr_value, &id1);
                    free(attr_value);
                    break;
                }
            }
            if (id1 >= 0) {
                state->act_node = id1;
            } else {
                state->act_node = -1;
                igraph_i_graphml_sax_handler_error(state, "Node with missing id encountered");
                return;
            }
            state->st = INSIDE_NODE;
        } else if (xmlStrEqual(localname, toXmlChar("data"))) {
            igraph_i_graphml_attribute_data_setup(state, attributes, nb_attributes,
                                                  IGRAPH_ATTRIBUTE_GRAPH);
            igraph_vector_int_push_back(&state->prev_state_stack, state->st);
            state->st = INSIDE_DATA;
        } else {
            igraph_i_graphml_handle_unknown_start_tag(state);
        }
        break;

    case INSIDE_NODE:
        if (xmlStrEqual(localname, toXmlChar("data"))) {
            igraph_i_graphml_attribute_data_setup(state, attributes, nb_attributes,
                                                  IGRAPH_ATTRIBUTE_VERTEX);
            igraph_vector_int_push_back(&state->prev_state_stack, state->st);
            state->st = INSIDE_DATA;
        }
        break;

    case INSIDE_EDGE:
        if (xmlStrEqual(localname, toXmlChar("data"))) {
            igraph_i_graphml_attribute_data_setup(state, attributes, nb_attributes,
                                                  IGRAPH_ATTRIBUTE_EDGE);
            igraph_vector_int_push_back(&state->prev_state_stack, state->st);
            state->st = INSIDE_DATA;
        }
        break;

    case INSIDE_DATA:
        /* We do not expect any new tags within a <data> tag */
        igraph_i_graphml_handle_unknown_start_tag(state);
        break;

    case UNKNOWN:
        igraph_i_graphml_handle_unknown_start_tag(state);
        break;

    default:
        break;
    }
}

static void igraph_i_graphml_sax_handler_end_element_ns(
        void *state0,
        const xmlChar* localname, const xmlChar* prefix,
        const xmlChar* uri) {
    struct igraph_i_graphml_parser_state *state =
        (struct igraph_i_graphml_parser_state*)state0;

    if (!state->successful) {
        return;
    }

    IGRAPH_UNUSED(localname);
    IGRAPH_UNUSED(prefix);
    IGRAPH_UNUSED(uri);

    switch (state->st) {
    case INSIDE_GRAPHML:
        state->st = FINISH;
        break;

    case INSIDE_GRAPH:
        state->st = INSIDE_GRAPHML;
        break;

    case INSIDE_KEY:
        state->current_attr_record = 0;
        state->st = INSIDE_GRAPHML;
        break;

    case INSIDE_DEFAULT:
        igraph_i_graphml_attribute_default_value_finish(state);
        state->st = INSIDE_KEY;
        break;

    case INSIDE_NODE:
        state->st = INSIDE_GRAPH;
        break;

    case INSIDE_EDGE:
        state->st = INSIDE_GRAPH;
        break;

    case INSIDE_DATA:
        igraph_i_graphml_attribute_data_finish(state);
        state->st = igraph_vector_int_pop_back(&state->prev_state_stack);
        break;

    case UNKNOWN:
        state->unknown_depth--;
        if (!state->unknown_depth) {
            state->st = igraph_vector_int_pop_back(&state->prev_state_stack);
        }
        break;

    default:
        break;
    }
}

static void igraph_i_graphml_sax_handler_chars(void* state0, const xmlChar* ch, int len) {
    struct igraph_i_graphml_parser_state *state =
        (struct igraph_i_graphml_parser_state*)state0;

    if (!state->successful) {
        return;
    }

    switch (state->st) {
    case INSIDE_KEY:
        break;

    case INSIDE_DATA:
    case INSIDE_DEFAULT:
        igraph_i_graphml_append_to_data_char(state, ch, len);
        break;

    default:
        /* just ignore it */
        break;
    }
}

static xmlSAXHandler igraph_i_graphml_sax_handler = {
    /* internalSubset = */ 0,
    /* isStandalone = */ 0,
    /* hasInternalSubset = */ 0,
    /* hasExternalSubset = */ 0,
    /* resolveEntity = */ 0,
    /* getEntity = */ igraph_i_graphml_sax_handler_get_entity,
    /* entityDecl = */ 0,
    /* notationDecl = */ 0,
    /* attributeDecl = */ 0,
    /* elementDecl = */ 0,
    /* unparsedEntityDecl = */ 0,
    /* setDocumentLocator = */ 0,
    /* startDocument = */ igraph_i_graphml_sax_handler_start_document,
    /* endDocument = */ igraph_i_graphml_sax_handler_end_document,
    /* startElement = */ 0,
    /* endElement = */ 0,
    /* reference = */ 0,
    /* characters = */ igraph_i_graphml_sax_handler_chars,
    /* ignorableWhitespaceFunc = */ 0,
    /* processingInstruction = */ 0,
    /* comment = */ 0,
    /* warning = */ igraph_i_graphml_sax_handler_error,
    /* error = */ igraph_i_graphml_sax_handler_error,
    /* fatalError = */ igraph_i_graphml_sax_handler_error,
    /* getParameterEntity = */ 0,
    /* cdataBlock = */ 0,
    /* externalSubset = */ 0,
    /* initialized = */ XML_SAX2_MAGIC,
    /* _private = */ 0,
    /* startElementNs = */ igraph_i_graphml_sax_handler_start_element_ns,
    /* endElementNs = */ igraph_i_graphml_sax_handler_end_element_ns,
    /* serror = */ 0
};

#endif

#define IS_FORBIDDEN_CONTROL_CHAR(x) ((x) < ' ' && (x) != '\t' && (x) != '\r' && (x) != '\n')

static int igraph_i_xml_escape(char* src, char** dest) {
    long int destlen = 0;
    char *s, *d;
    unsigned char ch;

    for (s = src; *s; s++, destlen++) {
        ch = (unsigned char)(*s);
        if (ch == '&') {
            destlen += 4;
        } else if (ch == '<') {
            destlen += 3;
        } else if (ch == '>') {
            destlen += 3;
        } else if (ch == '"') {
            destlen += 5;
        } else if (ch == '\'') {
            destlen += 5;
        } else if (IS_FORBIDDEN_CONTROL_CHAR(ch)) {
            char msg[4096];
            snprintf(msg, 4096, "Forbidden control character 0x%02X found in igraph_i_xml_escape",
                     ch);
            IGRAPH_ERROR(msg, IGRAPH_EINVAL);
        }
    }
    *dest = IGRAPH_CALLOC(destlen + 1, char);
    if (!*dest) {
        IGRAPH_ERROR("Not enough memory", IGRAPH_ENOMEM);
    }
    for (s = src, d = *dest; *s; s++, d++) {
        ch = (unsigned char)(*s);
        switch (ch) {
        case '&':
            strcpy(d, "&amp;"); d += 4; break;
        case '<':
            strcpy(d, "&lt;"); d += 3; break;
        case '>':
            strcpy(d, "&gt;"); d += 3; break;
        case '"':
            strcpy(d, "&quot;"); d += 5; break;
        case '\'':
            strcpy(d, "&apos;"); d += 5; break;
        default:
            *d = ch;
        }
    }
    *d = 0;
    return 0;
}

/**
 * \ingroup loadsave
 * \function igraph_read_graph_graphml
 * \brief Reads a graph from a GraphML file.
 *
 * </para><para>
 * GraphML is an XML-based file format for representing various types of
 * graphs. Currently only the most basic import functionality is implemented
 * in igraph: it can read GraphML files without nested graphs and hyperedges.
 * Attributes of the graph are loaded only if an attribute interface
 * is attached, i.e. if you use igraph from R or Python.
 *
 * </para><para>
 * Graph attribute names are taken from the <code>attr.name</code> attributes of the
 * \c key tags in the GraphML file. Since <code>attr.name</code> is not mandatory,
 * igraph will fall back to the \c id attribute of the \c key tag if
 * <code>attr.name</code> is missing.
 *
 * \param graph Pointer to an uninitialized graph object.
 * \param instream A stream, it should be readable.
 * \param index If the GraphML file contains more than one graph, the one
 *              specified by this index will be loaded. Indices start from
 *              zero, so supply zero here if your GraphML file contains only
 *              a single graph.
 *
 * \return Error code:
 *         \c IGRAPH_PARSEERROR: if there is a
 *         problem reading the file, or the file is syntactically
 *         incorrect.
 *         \c IGRAPH_UNIMPLEMENTED: the GraphML functionality was disabled
 *         at compile-time
 *
 * \example examples/simple/graphml.c
 */
int igraph_read_graph_graphml(igraph_t *graph, FILE *instream,
                              int index) {

#if HAVE_LIBXML == 1
    xmlParserCtxtPtr ctxt;
    struct igraph_i_graphml_parser_state state;
    int res;
    char buffer[4096];

    if (index < 0) {
        IGRAPH_ERROR("Graph index must be non-negative", IGRAPH_EINVAL);
    }

    xmlInitParser();

    /* Create a progressive parser context */
    state.g = graph;
    state.index = index < 0 ? 0 : index;
    res = (int) fread(buffer, 1, 4096, instream);
    ctxt = xmlCreatePushParserCtxt(&igraph_i_graphml_sax_handler,
                                   &state,
                                   buffer,
                                   res,
                                   NULL);
    /*   ctxt=xmlCreateIOParserCtxt(&igraph_i_graphml_sax_handler, &state, */
    /*               igraph_i_libxml2_read_callback, */
    /*               igraph_i_libxml2_close_callback, */
    /*               instream, XML_CHAR_ENCODING_NONE); */
    if (ctxt == NULL) {
        IGRAPH_ERROR("Can't create progressive parser context", IGRAPH_PARSEERROR);
    }

    /* Set parsing options */
    if (xmlCtxtUseOptions(ctxt,
                          XML_PARSE_NOENT | XML_PARSE_NOBLANKS |
                          XML_PARSE_NONET | XML_PARSE_NSCLEAN |
                          XML_PARSE_NOCDATA | XML_PARSE_HUGE
                         )) {
        IGRAPH_ERROR("Cannot set options for the parser context", IGRAPH_EINVAL);
    }

    /* Parse the file */
    while ((res = (int) fread(buffer, 1, 4096, instream)) > 0) {
        xmlParseChunk(ctxt, buffer, res, 0);
        if (!state.successful) {
            break;
        }
    }
    xmlParseChunk(ctxt, buffer, res, 1);

    /* Free the context */
    xmlFreeParserCtxt(ctxt);
    if (!state.successful) {
        if (state.error_message != 0) {
            IGRAPH_ERROR(state.error_message, IGRAPH_PARSEERROR);
        } else {
            IGRAPH_ERROR("Malformed GraphML file", IGRAPH_PARSEERROR);
        }
    }
    if (state.index >= 0) {
        IGRAPH_ERROR("Graph index was too large", IGRAPH_EINVAL);
    }

    return 0;
#else
    IGRAPH_UNUSED(graph);
    IGRAPH_UNUSED(instream);
    IGRAPH_UNUSED(index);

    IGRAPH_ERROR("GraphML support is disabled", IGRAPH_UNIMPLEMENTED);
#endif
}

/**
 * \ingroup loadsave
 * \function igraph_write_graph_graphml
 * \brief Writes the graph to a file in GraphML format
 *
 * </para><para>
 * GraphML is an XML-based file format for representing various types of
 * graphs. See the GraphML Primer (http://graphml.graphdrawing.org/primer/graphml-primer.html)
 * for detailed format description.
 *
 * \param graph The graph to write.
 * \param outstream The stream object to write to, it should be
 *        writable.
 * \param prefixattr Logical value, whether to put a prefix in front of the
 *        attribute names to ensure uniqueness if the graph has vertex and
 *        edge (or graph) attributes with the same name.
 * \return Error code:
 *         \c IGRAPH_EFILE if there is an error
 *         writing the file.
 *
 * Time complexity: O(|V|+|E|) otherwise. All
 * file operations are expected to have time complexity
 * O(1).
 *
 * \example examples/simple/graphml.c
 */
int igraph_write_graph_graphml(const igraph_t *graph, FILE *outstream,
                               igraph_bool_t prefixattr) {
    int ret;
    igraph_integer_t l, vc;
    igraph_eit_t it;
    igraph_strvector_t gnames, vnames, enames;
    igraph_vector_t gtypes, vtypes, etypes;
    long int i;
    igraph_vector_t numv;
    igraph_strvector_t strv;
    igraph_vector_bool_t boolv;
    const char *gprefix = prefixattr ? "g_" : "";
    const char *vprefix = prefixattr ? "v_" : "";
    const char *eprefix = prefixattr ? "e_" : "";

    /* set standard C locale lest we sometimes get commas instead of dots */
    char *saved_locale = strdup(setlocale(LC_NUMERIC, NULL));
    if (saved_locale == NULL) {
        IGRAPH_ERROR("Not enough memory", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, saved_locale);
    setlocale(LC_NUMERIC, "C");

    ret = fprintf(outstream, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
    if (ret < 0) {
        IGRAPH_ERROR("Write failed", IGRAPH_EFILE);
    }
    ret = fprintf(outstream, "<graphml xmlns=\"%s\"\n", GRAPHML_NAMESPACE_URI);
    if (ret < 0) {
        IGRAPH_ERROR("Write failed", IGRAPH_EFILE);
    }
    ret = fprintf(outstream, "         xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"\n");
    if (ret < 0) {
        IGRAPH_ERROR("Write failed", IGRAPH_EFILE);
    }
    ret = fprintf(outstream, "         xsi:schemaLocation=\"%s\n", GRAPHML_NAMESPACE_URI);
    if (ret < 0) {
        IGRAPH_ERROR("Write failed", IGRAPH_EFILE);
    }
    ret = fprintf(outstream, "         %s/1.0/graphml.xsd\">\n", GRAPHML_NAMESPACE_URI);
    if (ret < 0) {
        IGRAPH_ERROR("Write failed", IGRAPH_EFILE);
    }
    ret = fprintf(outstream, "<!-- Created by igraph -->\n");
    if (ret < 0) {
        IGRAPH_ERROR("Write failed", IGRAPH_EFILE);
    }

    /* dump the <key> elements if any */

    IGRAPH_VECTOR_INIT_FINALLY(&numv, 1);
    IGRAPH_STRVECTOR_INIT_FINALLY(&strv, 1);
    IGRAPH_VECTOR_BOOL_INIT_FINALLY(&boolv, 1);

    IGRAPH_STRVECTOR_INIT_FINALLY(&gnames, 0);
    IGRAPH_STRVECTOR_INIT_FINALLY(&vnames, 0);
    IGRAPH_STRVECTOR_INIT_FINALLY(&enames, 0);
    IGRAPH_VECTOR_INIT_FINALLY(&gtypes, 0);
    IGRAPH_VECTOR_INIT_FINALLY(&vtypes, 0);
    IGRAPH_VECTOR_INIT_FINALLY(&etypes, 0);
    igraph_i_attribute_get_info(graph,
                                &gnames, &gtypes,
                                &vnames, &vtypes,
                                &enames, &etypes);

    /* graph attributes */
    for (i = 0; i < igraph_vector_size(&gtypes); i++) {
        char *name, *name_escaped;
        igraph_strvector_get(&gnames, i, &name);
        IGRAPH_CHECK(igraph_i_xml_escape(name, &name_escaped));
        if (VECTOR(gtypes)[i] == IGRAPH_ATTRIBUTE_STRING) {
            ret = fprintf(outstream, "  <key id=\"%s%s\" for=\"graph\" attr.name=\"%s\" attr.type=\"string\"/>\n", gprefix, name_escaped, name_escaped);
            if (ret < 0) {
                IGRAPH_ERROR("Write failed", IGRAPH_EFILE);
            }
        } else if (VECTOR(gtypes)[i] == IGRAPH_ATTRIBUTE_NUMERIC) {
            ret = fprintf(outstream, "  <key id=\"%s%s\" for=\"graph\" attr.name=\"%s\" attr.type=\"double\"/>\n", gprefix, name_escaped, name_escaped);
            if (ret < 0) {
                IGRAPH_ERROR("Write failed", IGRAPH_EFILE);
            }
        } else if (VECTOR(gtypes)[i] == IGRAPH_ATTRIBUTE_BOOLEAN) {
            ret = fprintf(outstream, "  <key id=\"%s%s\" for=\"graph\" attr.name=\"%s\" attr.type=\"boolean\"/>\n", gprefix, name_escaped, name_escaped);
            if (ret < 0) {
                IGRAPH_ERROR("Write failed", IGRAPH_EFILE);
            }
        }
        IGRAPH_FREE(name_escaped);
    }

    /* vertex attributes */
    for (i = 0; i < igraph_vector_size(&vtypes); i++) {
        char *name, *name_escaped;
        igraph_strvector_get(&vnames, i, &name);
        IGRAPH_CHECK(igraph_i_xml_escape(name, &name_escaped));
        if (VECTOR(vtypes)[i] == IGRAPH_ATTRIBUTE_STRING) {
            ret = fprintf(outstream, "  <key id=\"%s%s\" for=\"node\" attr.name=\"%s\" attr.type=\"string\"/>\n", vprefix, name_escaped, name_escaped);
            if (ret < 0) {
                IGRAPH_ERROR("Write failed", IGRAPH_EFILE);
            }
        } else if (VECTOR(vtypes)[i] == IGRAPH_ATTRIBUTE_NUMERIC) {
            ret = fprintf(outstream, "  <key id=\"%s%s\" for=\"node\" attr.name=\"%s\" attr.type=\"double\"/>\n", vprefix, name_escaped, name_escaped);
            if (ret < 0) {
                IGRAPH_ERROR("Write failed", IGRAPH_EFILE);
            }
        } else if (VECTOR(vtypes)[i] == IGRAPH_ATTRIBUTE_BOOLEAN) {
            ret = fprintf(outstream, "  <key id=\"%s%s\" for=\"node\" attr.name=\"%s\" attr.type=\"boolean\"/>\n", vprefix, name_escaped, name_escaped);
            if (ret < 0) {
                IGRAPH_ERROR("Write failed", IGRAPH_EFILE);
            }
        }
        IGRAPH_FREE(name_escaped);
    }

    /* edge attributes */
    for (i = 0; i < igraph_vector_size(&etypes); i++) {
        char *name, *name_escaped;
        igraph_strvector_get(&enames, i, &name);
        IGRAPH_CHECK(igraph_i_xml_escape(name, &name_escaped));
        if (VECTOR(etypes)[i] == IGRAPH_ATTRIBUTE_STRING) {
            ret = fprintf(outstream, "  <key id=\"%s%s\" for=\"edge\" attr.name=\"%s\" attr.type=\"string\"/>\n", eprefix, name_escaped, name_escaped);
            if (ret < 0) {
                IGRAPH_ERROR("Write failed", IGRAPH_EFILE);
            }
        } else if (VECTOR(etypes)[i] == IGRAPH_ATTRIBUTE_NUMERIC) {
            ret = fprintf(outstream, "  <key id=\"%s%s\" for=\"edge\" attr.name=\"%s\" attr.type=\"double\"/>\n", eprefix, name_escaped, name_escaped);
            if (ret < 0) {
                IGRAPH_ERROR("Write failed", IGRAPH_EFILE);
            }
        } else if (VECTOR(etypes)[i] == IGRAPH_ATTRIBUTE_BOOLEAN) {
            ret = fprintf(outstream, "  <key id=\"%s%s\" for=\"edge\" attr.name=\"%s\" attr.type=\"boolean\"/>\n", eprefix, name_escaped, name_escaped);
            if (ret < 0) {
                IGRAPH_ERROR("Write failed", IGRAPH_EFILE);
            }
        }
        IGRAPH_FREE(name_escaped);
    }

    ret = fprintf(outstream, "  <graph id=\"G\" edgedefault=\"%s\">\n", (igraph_is_directed(graph) ? "directed" : "undirected"));
    if (ret < 0) {
        IGRAPH_ERROR("Write failed", IGRAPH_EFILE);
    }

    /* Write the graph atributes before anything else */

    for (i = 0; i < igraph_vector_size(&gtypes); i++) {
        char *name, *name_escaped;
        if (VECTOR(gtypes)[i] == IGRAPH_ATTRIBUTE_NUMERIC) {
            igraph_strvector_get(&gnames, i, &name);
            IGRAPH_CHECK(igraph_i_attribute_get_numeric_graph_attr(graph, name, &numv));
            if (!isnan(VECTOR(numv)[0])) {
                IGRAPH_CHECK(igraph_i_xml_escape(name, &name_escaped));
                ret = fprintf(outstream, "    <data key=\"%s%s\">", gprefix, name_escaped);
                IGRAPH_FREE(name_escaped);
                if (ret < 0) {
                    IGRAPH_ERROR("Write failed", IGRAPH_EFILE);
                }
                ret = igraph_real_fprintf_precise(outstream, VECTOR(numv)[0]);
                if (ret < 0) {
                    IGRAPH_ERROR("Write failed", IGRAPH_EFILE);
                }
                ret = fprintf(outstream, "</data>\n");
                if (ret < 0) {
                    IGRAPH_ERROR("Write failed", IGRAPH_EFILE);
                }
            }
        } else if (VECTOR(gtypes)[i] == IGRAPH_ATTRIBUTE_STRING) {
            char *s, *s_escaped;
            igraph_strvector_get(&gnames, i, &name);
            IGRAPH_CHECK(igraph_i_xml_escape(name, &name_escaped));
            ret = fprintf(outstream, "    <data key=\"%s%s\">", gprefix,
                          name_escaped);
            IGRAPH_FREE(name_escaped);
            IGRAPH_CHECK(igraph_i_attribute_get_string_graph_attr(graph, name, &strv));
            igraph_strvector_get(&strv, 0, &s);
            IGRAPH_CHECK(igraph_i_xml_escape(s, &s_escaped));
            ret = fprintf(outstream, "%s", s_escaped);
            IGRAPH_FREE(s_escaped);
            if (ret < 0) {
                IGRAPH_ERROR("Write failed", IGRAPH_EFILE);
            }
            ret = fprintf(outstream, "</data>\n");
            if (ret < 0) {
                IGRAPH_ERROR("Write failed", IGRAPH_EFILE);
            }
        } else if (VECTOR(gtypes)[i] == IGRAPH_ATTRIBUTE_BOOLEAN) {
            igraph_strvector_get(&gnames, i, &name);
            IGRAPH_CHECK(igraph_i_attribute_get_bool_graph_attr(graph, name, &boolv));
            IGRAPH_CHECK(igraph_i_xml_escape(name, &name_escaped));
            ret = fprintf(outstream, "    <data key=\"%s%s\">%s</data>\n",
                          gprefix, name_escaped, VECTOR(boolv)[0] ? "true" : "false");
            IGRAPH_FREE(name_escaped);
            if (ret < 0) {
                IGRAPH_ERROR("Write failed", IGRAPH_EFILE);
            }
        }
    }

    /* Let's dump the nodes first */
    vc = igraph_vcount(graph);
    for (l = 0; l < vc; l++) {
        char *name, *name_escaped;
        ret = fprintf(outstream, "    <node id=\"n%ld\">\n", (long)l);

        if (ret < 0) {
            IGRAPH_ERROR("Write failed", IGRAPH_EFILE);
        }

        for (i = 0; i < igraph_vector_size(&vtypes); i++) {
            if (VECTOR(vtypes)[i] == IGRAPH_ATTRIBUTE_NUMERIC) {
                igraph_strvector_get(&vnames, i, &name);
                IGRAPH_CHECK(igraph_i_attribute_get_numeric_vertex_attr(graph, name,
                             igraph_vss_1(l), &numv));
                if (!isnan(VECTOR(numv)[0])) {
                    IGRAPH_CHECK(igraph_i_xml_escape(name, &name_escaped));
                    ret = fprintf(outstream, "      <data key=\"%s%s\">", vprefix, name_escaped);
                    IGRAPH_FREE(name_escaped);
                    if (ret < 0) {
                        IGRAPH_ERROR("Write failed", IGRAPH_EFILE);
                    }
                    ret = igraph_real_fprintf_precise(outstream, VECTOR(numv)[0]);
                    if (ret < 0) {
                        IGRAPH_ERROR("Write failed", IGRAPH_EFILE);
                    }
                    ret = fprintf(outstream, "</data>\n");
                    if (ret < 0) {
                        IGRAPH_ERROR("Write failed", IGRAPH_EFILE);
                    }
                }
            } else if (VECTOR(vtypes)[i] == IGRAPH_ATTRIBUTE_STRING) {
                char *s, *s_escaped;
                igraph_strvector_get(&vnames, i, &name);
                IGRAPH_CHECK(igraph_i_xml_escape(name, &name_escaped));
                ret = fprintf(outstream, "      <data key=\"%s%s\">", vprefix,
                              name_escaped);
                IGRAPH_FREE(name_escaped);
                IGRAPH_CHECK(igraph_i_attribute_get_string_vertex_attr(graph, name,
                             igraph_vss_1(l), &strv));
                igraph_strvector_get(&strv, 0, &s);
                IGRAPH_CHECK(igraph_i_xml_escape(s, &s_escaped));
                ret = fprintf(outstream, "%s", s_escaped);
                IGRAPH_FREE(s_escaped);
                if (ret < 0) {
                    IGRAPH_ERROR("Write failed", IGRAPH_EFILE);
                }
                ret = fprintf(outstream, "</data>\n");
                if (ret < 0) {
                    IGRAPH_ERROR("Write failed", IGRAPH_EFILE);
                }
            } else if (VECTOR(vtypes)[i] == IGRAPH_ATTRIBUTE_BOOLEAN) {
                igraph_strvector_get(&vnames, i, &name);
                IGRAPH_CHECK(igraph_i_attribute_get_bool_vertex_attr(graph, name,
                             igraph_vss_1(l), &boolv));
                IGRAPH_CHECK(igraph_i_xml_escape(name, &name_escaped));
                ret = fprintf(outstream, "      <data key=\"%s%s\">%s</data>\n",
                              vprefix, name_escaped, VECTOR(boolv)[0] ? "true" : "false");
                IGRAPH_FREE(name_escaped);
                if (ret < 0) {
                    IGRAPH_ERROR("Write failed", IGRAPH_EFILE);
                }
            }
        }

        ret = fprintf(outstream, "    </node>\n");
        if (ret < 0) {
            IGRAPH_ERROR("Write failed", IGRAPH_EFILE);
        }
    }

    /* Now the edges */
    IGRAPH_CHECK(igraph_eit_create(graph, igraph_ess_all(0), &it));
    IGRAPH_FINALLY(igraph_eit_destroy, &it);
    while (!IGRAPH_EIT_END(it)) {
        igraph_integer_t from, to;
        char *name, *name_escaped;
        long int edge = IGRAPH_EIT_GET(it);
        igraph_edge(graph, (igraph_integer_t) edge, &from, &to);
        ret = fprintf(outstream, "    <edge source=\"n%ld\" target=\"n%ld\">\n",
                      (long int)from, (long int)to);
        if (ret < 0) {
            IGRAPH_ERROR("Write failed", IGRAPH_EFILE);
        }

        for (i = 0; i < igraph_vector_size(&etypes); i++) {
            if (VECTOR(etypes)[i] == IGRAPH_ATTRIBUTE_NUMERIC) {
                igraph_strvector_get(&enames, i, &name);
                IGRAPH_CHECK(igraph_i_attribute_get_numeric_edge_attr(graph, name,
                             igraph_ess_1((igraph_integer_t) edge), &numv));
                if (!isnan(VECTOR(numv)[0])) {
                    IGRAPH_CHECK(igraph_i_xml_escape(name, &name_escaped));
                    ret = fprintf(outstream, "      <data key=\"%s%s\">", eprefix, name_escaped);
                    IGRAPH_FREE(name_escaped);
                    if (ret < 0) {
                        IGRAPH_ERROR("Write failed", IGRAPH_EFILE);
                    }
                    ret = igraph_real_fprintf_precise(outstream, VECTOR(numv)[0]);
                    if (ret < 0) {
                        IGRAPH_ERROR("Write failed", IGRAPH_EFILE);
                    }
                    ret = fprintf(outstream, "</data>\n");
                    if (ret < 0) {
                        IGRAPH_ERROR("Write failed", IGRAPH_EFILE);
                    }
                }
            } else if (VECTOR(etypes)[i] == IGRAPH_ATTRIBUTE_STRING) {
                char *s, *s_escaped;
                igraph_strvector_get(&enames, i, &name);
                IGRAPH_CHECK(igraph_i_xml_escape(name, &name_escaped));
                ret = fprintf(outstream, "      <data key=\"%s%s\">", eprefix,
                              name_escaped);
                IGRAPH_FREE(name_escaped);
                IGRAPH_CHECK(igraph_i_attribute_get_string_edge_attr(graph, name,
                             igraph_ess_1((igraph_integer_t) edge), &strv));
                igraph_strvector_get(&strv, 0, &s);
                IGRAPH_CHECK(igraph_i_xml_escape(s, &s_escaped));
                ret = fprintf(outstream, "%s", s_escaped);
                IGRAPH_FREE(s_escaped);
                if (ret < 0) {
                    IGRAPH_ERROR("Write failed", IGRAPH_EFILE);
                }
                ret = fprintf(outstream, "</data>\n");
                if (ret < 0) {
                    IGRAPH_ERROR("Write failed", IGRAPH_EFILE);
                }
            } else if (VECTOR(etypes)[i] == IGRAPH_ATTRIBUTE_BOOLEAN) {
                igraph_strvector_get(&enames, i, &name);
                IGRAPH_CHECK(igraph_i_attribute_get_bool_edge_attr(graph, name,
                             igraph_ess_1((igraph_integer_t) edge), &boolv));
                IGRAPH_CHECK(igraph_i_xml_escape(name, &name_escaped));
                ret = fprintf(outstream, "      <data key=\"%s%s\">%s</data>\n",
                              eprefix, name_escaped, VECTOR(boolv)[0] ? "true" : "false");
                IGRAPH_FREE(name_escaped);
                if (ret < 0) {
                    IGRAPH_ERROR("Write failed", IGRAPH_EFILE);
                }
            }
        }

        ret = fprintf(outstream, "    </edge>\n");
        if (ret < 0) {
            IGRAPH_ERROR("Write failed", IGRAPH_EFILE);
        }
        IGRAPH_EIT_NEXT(it);
    }
    igraph_eit_destroy(&it);
    IGRAPH_FINALLY_CLEAN(1);

    ret = fprintf(outstream, "  </graph>\n");
    if (ret < 0) {
        IGRAPH_ERROR("Write failed", IGRAPH_EFILE);
    }
    fprintf(outstream, "</graphml>\n");
    if (ret < 0) {
        IGRAPH_ERROR("Write failed", IGRAPH_EFILE);
    }

    /* reset locale to whatever was before this function */
    setlocale(LC_NUMERIC, saved_locale);

    igraph_free(saved_locale);
    igraph_strvector_destroy(&gnames);
    igraph_strvector_destroy(&vnames);
    igraph_strvector_destroy(&enames);
    igraph_vector_destroy(&gtypes);
    igraph_vector_destroy(&vtypes);
    igraph_vector_destroy(&etypes);
    igraph_vector_destroy(&numv);
    igraph_strvector_destroy(&strv);
    igraph_vector_bool_destroy(&boolv);
    IGRAPH_FINALLY_CLEAN(10);

    return 0;
}
