/*
   IGraph library.
   Copyright (C) 2006-2023  The igraph development team <igraph@igraph.org>

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

#include "igraph_foreign.h"

#include "igraph_attributes.h"
#include "igraph_interface.h"
#include "igraph_memory.h"

#include "core/interruption.h"
#include "core/trie.h"
#include "graph/attributes.h"
#include "internal/hacks.h" /* strcasecmp & strdup */
#include "io/parse_utils.h"

#include "config.h"

#include <ctype.h>   /* isdigit */
#include <math.h>    /* isnan */
#include <string.h>
#include <stdarg.h>  /* va_start & co */

#define GRAPHML_NAMESPACE_URI "http://graphml.graphdrawing.org/xmlns"

#if HAVE_LIBXML == 1
#include <libxml/globals.h>
#include <libxml/parser.h>

xmlEntity blankEntityStruct = {
#ifndef XML_WITHOUT_CORBA
    NULL, /* _private */
#endif
    XML_ENTITY_DECL, /* type */
    NULL, /* name */
    NULL, /* children */
    NULL, /* last */
    NULL, /* parent */
    NULL, /* next */
    NULL, /* prev */
    NULL, /* doc */

    NULL, /* orig */
    NULL, /* content */
    0,    /* length */
    XML_EXTERNAL_GENERAL_PARSED_ENTITY, /* etype */
    NULL, /* ExternalID */
    NULL, /* SystemID */

    NULL, /* nexte */
    NULL, /* URI */
    0,    /* owner */
#if LIBXML_VERSION < 21100   /* Versions < 2.11.0: */
    1     /* checked */
#else                        /* Starting with verson 2.11.0: */
    1,    /* flags */
    0     /* expandedSize */
#endif
};

xmlEntityPtr blankEntity = &blankEntityStruct;

#define toXmlChar(a)   (BAD_CAST(a))
#define fromXmlChar(a) ((char *)(a)) /* not the most elegant way... */

#define GRAPHML_PARSE_ERROR_WITH_CODE(state, msg, code) \
    do {  \
        if (state->successful) { \
            igraph_i_graphml_sax_handler_error(state, msg); \
        } \
    } while (0)
#define RETURN_GRAPHML_PARSE_ERROR_WITH_CODE(state, msg, code) \
    do { \
        GRAPHML_PARSE_ERROR_WITH_CODE(state, msg, code); \
        return; \
    } while (0)

typedef struct igraph_i_graphml_attribute_record_t {
    const char *id;           /* GraphML id */
    enum { I_GRAPHML_BOOLEAN, I_GRAPHML_INTEGER, I_GRAPHML_LONG,
           I_GRAPHML_FLOAT, I_GRAPHML_DOUBLE, I_GRAPHML_STRING,
           I_GRAPHML_UNKNOWN_TYPE
         } type; /* GraphML type */
    union {
        igraph_real_t as_numeric;
        igraph_bool_t as_boolean;
        char *as_string;
    } default_value;   /* Default value of the attribute, if any */
    igraph_attribute_record_t record;
} igraph_i_graphml_attribute_record_t;

typedef enum {
    START, INSIDE_GRAPHML, INSIDE_GRAPH, INSIDE_NODE, INSIDE_EDGE,
    INSIDE_KEY, INSIDE_DEFAULT, INSIDE_DATA, FINISH, UNKNOWN, ERROR
} igraph_i_graphml_parser_state_index_t;

struct igraph_i_graphml_parser_state {
    igraph_i_graphml_parser_state_index_t st;
    igraph_t *g;
    igraph_trie_t node_trie;
    igraph_strvector_t edgeids;
    igraph_vector_int_t edgelist;
    igraph_vector_int_t prev_state_stack;
    unsigned int unknown_depth;
    igraph_integer_t index;
    igraph_bool_t successful;
    igraph_bool_t edges_directed;
    igraph_trie_t v_attr_ids;
    igraph_vector_ptr_t v_attrs;
    igraph_trie_t e_attr_ids;
    igraph_vector_ptr_t e_attrs;
    igraph_trie_t g_attr_ids;
    igraph_vector_ptr_t g_attrs;
    igraph_i_graphml_attribute_record_t* current_attr_record;
    xmlChar *data_key;
    igraph_attribute_elemtype_t data_type;
    char *error_message;
    igraph_vector_char_t data_char;
    igraph_integer_t act_node;
    igraph_bool_t ignore_namespaces;
};

static void igraph_i_report_unhandled_attribute_target(const char* target,
        const char* file, int line) {
    igraph_warningf("Attribute target '%s' is not handled; ignoring corresponding "
                    "attribute specifications.", file, line, target);
}

static igraph_error_t igraph_i_graphml_parse_numeric(
    const char* char_data, igraph_real_t* result, igraph_real_t default_value
) {
    const char* trimmed;
    size_t trimmed_length;

    if (char_data == 0) {
        *result = default_value;
        return IGRAPH_SUCCESS;
    }

    igraph_i_trim_whitespace(char_data, strlen(char_data), &trimmed, &trimmed_length);
    if (trimmed_length > 0) {
        IGRAPH_CHECK(igraph_i_parse_real(trimmed, trimmed_length, result));
    } else {
        *result = default_value;
    }

    return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_i_graphml_parse_boolean(
    const char* char_data, igraph_bool_t* result, igraph_bool_t default_value
) {
    igraph_integer_t value;
    const char* trimmed;
    size_t trimmed_length;

    if (char_data == 0) {
        *result = default_value;
        return IGRAPH_SUCCESS;
    }

    igraph_i_trim_whitespace(char_data, strlen(char_data), &trimmed, &trimmed_length);

    if (trimmed_length == 4 && !strncasecmp(trimmed, "true", trimmed_length)) {
        *result = 1;
        return IGRAPH_SUCCESS;
    }

    if (trimmed_length == 3 && !strncasecmp(trimmed, "yes", trimmed_length)) {
        *result = 1;
        return IGRAPH_SUCCESS;
    }

    if (trimmed_length == 5 && !strncasecmp(trimmed, "false", trimmed_length)) {
        *result = 0;
        return IGRAPH_SUCCESS;
    }

    if (trimmed_length == 2 && !strncasecmp(trimmed, "no", trimmed_length)) {
        *result = 0;
        return IGRAPH_SUCCESS;
    }

    if (trimmed_length > 0) {
        if (isdigit(trimmed[0])) {
            IGRAPH_CHECK(igraph_i_parse_integer(trimmed, trimmed_length, &value));
        } else {
            IGRAPH_ERRORF("Cannot parse '%.*s' as Boolean value.", IGRAPH_PARSEERROR,
                          (int) trimmed_length, trimmed);
        }
    } else {
        *result = default_value;
        return IGRAPH_SUCCESS;
    }

    *result = value != 0;
    return IGRAPH_SUCCESS;
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
            IGRAPH_FREE(rec->record.value);
        }
        if (rec->default_value.as_string != 0) {
            IGRAPH_FREE(rec->default_value.as_string);
        }
    } else if (rec->record.type == IGRAPH_ATTRIBUTE_BOOLEAN) {
        if (rec->record.value != 0) {
            igraph_vector_bool_destroy((igraph_vector_bool_t*)rec->record.value);
            IGRAPH_FREE(rec->record.value);
        }
    } else if (rec->record.type == IGRAPH_ATTRIBUTE_UNSPECIFIED) {
        /* no value was set */
    }
    if (rec->id != NULL) {
        xmlFree((void *) rec->id);
        rec->id = NULL;
    }
    if (rec->record.name != 0) {
        IGRAPH_FREE(rec->record.name);
    }
}

static igraph_error_t igraph_i_graphml_parser_state_init(struct igraph_i_graphml_parser_state* state, igraph_t* graph, igraph_integer_t index) {
    memset(state, 0, sizeof(struct igraph_i_graphml_parser_state));

    state->g = graph;
    state->index = index < 0 ? 0 : index;
    state->successful = 1;
    state->error_message = NULL;

    IGRAPH_CHECK(igraph_vector_int_init(&state->prev_state_stack, 0));
    IGRAPH_CHECK(igraph_vector_int_reserve(&state->prev_state_stack, 32));
    IGRAPH_FINALLY(igraph_vector_int_destroy, &state->prev_state_stack);

    IGRAPH_CHECK(igraph_vector_ptr_init(&state->v_attrs, 0));
    IGRAPH_VECTOR_PTR_SET_ITEM_DESTRUCTOR(&state->v_attrs,
                                          igraph_i_graphml_attribute_record_destroy);
    IGRAPH_FINALLY(igraph_vector_ptr_destroy_all, &state->v_attrs);

    IGRAPH_CHECK(igraph_vector_ptr_init(&state->e_attrs, 0));
    IGRAPH_VECTOR_PTR_SET_ITEM_DESTRUCTOR(&state->e_attrs,
                                          igraph_i_graphml_attribute_record_destroy);
    IGRAPH_FINALLY(igraph_vector_ptr_destroy_all, &state->e_attrs);

    IGRAPH_CHECK(igraph_vector_ptr_init(&state->g_attrs, 0));
    IGRAPH_VECTOR_PTR_SET_ITEM_DESTRUCTOR(&state->g_attrs,
                                          igraph_i_graphml_attribute_record_destroy);
    IGRAPH_FINALLY(igraph_vector_ptr_destroy_all, &state->g_attrs);

    IGRAPH_CHECK(igraph_vector_int_init(&state->edgelist, 0));
    IGRAPH_FINALLY(igraph_vector_int_destroy, &state->edgelist);

    IGRAPH_CHECK(igraph_trie_init(&state->node_trie, 1));
    IGRAPH_FINALLY(igraph_trie_destroy, &state->node_trie);

    IGRAPH_CHECK(igraph_strvector_init(&state->edgeids, 0));
    IGRAPH_FINALLY(igraph_strvector_destroy, &state->edgeids);

    IGRAPH_CHECK(igraph_trie_init(&state->v_attr_ids, 0));
    IGRAPH_FINALLY(igraph_trie_destroy, &state->v_attr_ids);

    IGRAPH_CHECK(igraph_trie_init(&state->e_attr_ids, 0));
    IGRAPH_FINALLY(igraph_trie_destroy, &state->e_attr_ids);

    IGRAPH_CHECK(igraph_trie_init(&state->g_attr_ids, 0));
    IGRAPH_FINALLY(igraph_trie_destroy, &state->g_attr_ids);

    IGRAPH_VECTOR_CHAR_INIT_FINALLY(&state->data_char, 0);

    IGRAPH_FINALLY_CLEAN(11);

    return IGRAPH_SUCCESS;
}

static void igraph_i_graphml_parser_state_destroy(struct igraph_i_graphml_parser_state* state) {
    igraph_trie_destroy(&state->node_trie);
    igraph_strvector_destroy(&state->edgeids);
    igraph_trie_destroy(&state->v_attr_ids);
    igraph_trie_destroy(&state->e_attr_ids);
    igraph_trie_destroy(&state->g_attr_ids);
    igraph_vector_int_destroy(&state->edgelist);
    igraph_vector_int_destroy(&state->prev_state_stack);

    igraph_vector_ptr_destroy_all(&state->v_attrs);
    igraph_vector_ptr_destroy_all(&state->e_attrs);
    igraph_vector_ptr_destroy_all(&state->g_attrs);

    igraph_vector_char_destroy(&state->data_char);

    if (state->data_key) {
        xmlFree((void *) state->data_key);
        state->data_key = NULL;
    }

    if (state->error_message) {
        IGRAPH_FREE(state->error_message);
    }
}

static void igraph_i_graphml_parser_state_set_error_from_varargs(
    struct igraph_i_graphml_parser_state *state, const char* msg, va_list args
) {
    const size_t max_error_message_length = 4096;

    state->successful = 0;
    state->st = ERROR;

    if (state->error_message == 0) {
        /* ownership of state->error_message passed on immediately to
         * state so the state destructor is responsible for freeing it */
        state->error_message = IGRAPH_CALLOC(max_error_message_length, char);
    }

    /* we need to guard against state->error_message == 0, which may happen
     * if the memory allocation for the error message itself failed */
    if (state->error_message != 0) {
        vsnprintf(state->error_message, max_error_message_length, msg, args);
    }
}

static void igraph_i_graphml_parser_state_set_error_from_xmlerror(
    struct igraph_i_graphml_parser_state *state, const xmlError *error
) {
    const size_t max_error_message_length = 4096;

    state->successful = 0;
    state->st = ERROR;

    if (state->error_message == 0) {
        /* ownership of state->error_message passed on immediately to
         * state so the state destructor is responsible for freeing it */
        state->error_message = IGRAPH_CALLOC(max_error_message_length, char);
    }

    /* we need to guard against state->error_message == 0, which may happen
     * if the memory allocation for the error message itself failed */
    if (state->error_message != 0) {
        snprintf(state->error_message, max_error_message_length, "Line %d: %s",
            error->line, error->message);
    }
}

static void igraph_i_graphml_sax_handler_error(void *state0, const char* msg, ...) {
    struct igraph_i_graphml_parser_state *state =
        (struct igraph_i_graphml_parser_state*)state0;
    va_list args;
    va_start(args, msg);
    igraph_i_graphml_parser_state_set_error_from_varargs(state, msg, args);
    va_end(args);
}

static xmlEntityPtr igraph_i_graphml_sax_handler_get_entity(void *state0,
        const xmlChar* name) {
    xmlEntityPtr predef = xmlGetPredefinedEntity(name);
    const char* entityName;

    IGRAPH_UNUSED(state0);
    if (predef != NULL) {
        return predef;
    }

    entityName = fromXmlChar(name);
    IGRAPH_WARNINGF("Unknown XML entity found: '%s'.", entityName);

    return blankEntity;
}

static igraph_error_t igraph_i_graphml_handle_unknown_start_tag(struct igraph_i_graphml_parser_state *state) {
    if (state->st != UNKNOWN) {
        IGRAPH_CHECK(igraph_vector_int_push_back(&state->prev_state_stack, state->st));
        state->st = UNKNOWN;
        state->unknown_depth = 1;
    } else {
        state->unknown_depth++;
    }
    return IGRAPH_SUCCESS;
}

static void igraph_i_graphml_sax_handler_start_document(void *state0) {
    struct igraph_i_graphml_parser_state *state =
        (struct igraph_i_graphml_parser_state*)state0;

    state->st = START;
    state->successful = 1;
    state->edges_directed = 0;
    state->data_key = NULL;
    state->unknown_depth = 0;
    state->ignore_namespaces = 0;
}

static igraph_error_t igraph_i_graphml_parser_state_finish_parsing(struct igraph_i_graphml_parser_state *state) {
    igraph_integer_t i, l;
    igraph_attribute_record_t idrec, eidrec;
    const char *idstr = "id";
    igraph_bool_t already_has_vertex_id = false, already_has_edge_id = false;
    igraph_vector_ptr_t vattr, eattr, gattr;
    igraph_integer_t esize;

    IGRAPH_ASSERT(state->successful);

    /* check that we have found and parsed the graph the user is interested in */
    IGRAPH_ASSERT(state->index < 0);

    IGRAPH_CHECK(igraph_vector_ptr_init(&vattr, igraph_vector_ptr_size(&state->v_attrs) + 1)); /* +1 for 'id' */
    IGRAPH_FINALLY(igraph_vector_ptr_destroy, &vattr);
    igraph_vector_ptr_resize(&vattr, 0); /* will be filled with push_back() */

    esize = igraph_vector_ptr_size(&state->e_attrs);
    if (igraph_strvector_size(&state->edgeids) != 0) {
        esize++;
    }
    IGRAPH_CHECK(igraph_vector_ptr_init(&eattr, esize));
    IGRAPH_FINALLY(igraph_vector_ptr_destroy, &eattr);
    igraph_vector_ptr_resize(&eattr, 0); /* will be filled with push_back() */

    IGRAPH_CHECK(igraph_vector_ptr_init(&gattr, igraph_vector_ptr_size(&state->g_attrs)));
    IGRAPH_FINALLY(igraph_vector_ptr_destroy, &gattr);
    igraph_vector_ptr_resize(&gattr, 0); /* will be filled with push_back() */

    for (i = 0; i < igraph_vector_ptr_size(&state->v_attrs); i++) {
        igraph_i_graphml_attribute_record_t *graphmlrec =
            VECTOR(state->v_attrs)[i];
        igraph_attribute_record_t *rec = &graphmlrec->record;

        /* Check that the name of the vertex attribute is not 'id'.
         * If it is then we cannot add the complementary 'id' attribute. */
        if (! strcmp(rec->name, idstr)) {
            already_has_vertex_id = 1;
        }

        if (rec->type == IGRAPH_ATTRIBUTE_NUMERIC) {
            igraph_vector_t *vec = (igraph_vector_t*)rec->value;
            igraph_integer_t origsize = igraph_vector_size(vec);
            igraph_integer_t nodes = igraph_trie_size(&state->node_trie);
            IGRAPH_CHECK(igraph_vector_resize(vec, nodes));
            for (l = origsize; l < nodes; l++) {
                VECTOR(*vec)[l] = graphmlrec->default_value.as_numeric;
            }
        } else if (rec->type == IGRAPH_ATTRIBUTE_STRING) {
            igraph_strvector_t *strvec = (igraph_strvector_t*)rec->value;
            igraph_integer_t origsize = igraph_strvector_size(strvec);
            igraph_integer_t nodes = igraph_trie_size(&state->node_trie);
            IGRAPH_CHECK(igraph_strvector_resize(strvec, nodes));
            for (l = origsize; l < nodes; l++) {
                IGRAPH_CHECK(igraph_strvector_set(strvec, l, graphmlrec->default_value.as_string));
            }
        } else if (rec->type == IGRAPH_ATTRIBUTE_BOOLEAN) {
            igraph_vector_bool_t *boolvec = (igraph_vector_bool_t*)rec->value;
            igraph_integer_t origsize = igraph_vector_bool_size(boolvec);
            igraph_integer_t nodes = igraph_trie_size(&state->node_trie);
            IGRAPH_CHECK(igraph_vector_bool_resize(boolvec, nodes));
            for (l = origsize; l < nodes; l++) {
                VECTOR(*boolvec)[l] = graphmlrec->default_value.as_boolean;
            }
        } else if (rec->type == IGRAPH_ATTRIBUTE_UNSPECIFIED) {
            continue; /* skipped attribute */
        }
        igraph_vector_ptr_push_back(&vattr, rec); /* reserved */
    }
    if (!already_has_vertex_id) {
        idrec.name = idstr;
        idrec.type = IGRAPH_ATTRIBUTE_STRING;
        idrec.value = igraph_i_trie_borrow_keys(&state->node_trie);
        igraph_vector_ptr_push_back(&vattr, &idrec); /* reserved */
    } else {
        IGRAPH_WARNING("Could not add vertex ids, there is already an 'id' vertex attribute.");
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
            igraph_integer_t origsize = igraph_vector_size(vec);
            igraph_integer_t edges = igraph_vector_int_size(&state->edgelist) / 2;
            IGRAPH_CHECK(igraph_vector_resize(vec, edges));
            for (l = origsize; l < edges; l++) {
                VECTOR(*vec)[l] = graphmlrec->default_value.as_numeric;
            }
        } else if (rec->type == IGRAPH_ATTRIBUTE_STRING) {
            igraph_strvector_t *strvec = (igraph_strvector_t*)rec->value;
            igraph_integer_t origsize = igraph_strvector_size(strvec);
            igraph_integer_t edges = igraph_vector_int_size(&state->edgelist) / 2;
            IGRAPH_CHECK(igraph_strvector_resize(strvec, edges));
            for (l = origsize; l < edges; l++) {
                IGRAPH_CHECK(igraph_strvector_set(strvec, l, graphmlrec->default_value.as_string));
            }
        } else if (rec->type == IGRAPH_ATTRIBUTE_BOOLEAN) {
            igraph_vector_bool_t *boolvec = (igraph_vector_bool_t*)rec->value;
            igraph_integer_t origsize = igraph_vector_bool_size(boolvec);
            igraph_integer_t edges = igraph_vector_int_size(&state->edgelist) / 2;
            IGRAPH_CHECK(igraph_vector_bool_resize(boolvec, edges));
            for (l = origsize; l < edges; l++) {
                VECTOR(*boolvec)[l] = graphmlrec->default_value.as_boolean;
            }
        } else if (rec->type == IGRAPH_ATTRIBUTE_UNSPECIFIED) {
            continue; /* skipped attribute */
        }
        igraph_vector_ptr_push_back(&eattr, rec); /* reserved */
    }
    if (igraph_strvector_size(&state->edgeids) != 0) {
        if (!already_has_edge_id) {
            igraph_integer_t origsize = igraph_strvector_size(&state->edgeids);
            eidrec.name = idstr;
            eidrec.type = IGRAPH_ATTRIBUTE_STRING;
            IGRAPH_CHECK(igraph_strvector_resize(&state->edgeids, igraph_vector_int_size(&state->edgelist) / 2));
            for (; origsize < igraph_strvector_size(&state->edgeids); origsize++) {
                IGRAPH_CHECK(igraph_strvector_set(&state->edgeids, origsize, ""));
            }
            eidrec.value = &state->edgeids;
            igraph_vector_ptr_push_back(&eattr, &eidrec); /* reserved */
        } else {
            IGRAPH_WARNING("Could not add edge ids, there is already an 'id' edge attribute.");
        }
    }

    for (i = 0; i < igraph_vector_ptr_size(&state->g_attrs); i++) {
        igraph_i_graphml_attribute_record_t *graphmlrec =
            VECTOR(state->g_attrs)[i];
        igraph_attribute_record_t *rec = &graphmlrec->record;
        if (rec->type == IGRAPH_ATTRIBUTE_NUMERIC) {
            igraph_vector_t *vec = (igraph_vector_t*)rec->value;
            igraph_integer_t origsize = igraph_vector_size(vec);
            IGRAPH_CHECK(igraph_vector_resize(vec, 1));
            for (l = origsize; l < 1; l++) {
                VECTOR(*vec)[l] = graphmlrec->default_value.as_numeric;
            }
        } else if (rec->type == IGRAPH_ATTRIBUTE_STRING) {
            igraph_strvector_t *strvec = (igraph_strvector_t*)rec->value;
            igraph_integer_t origsize = igraph_strvector_size(strvec);
            IGRAPH_CHECK(igraph_strvector_resize(strvec, 1));
            for (l = origsize; l < 1; l++) {
                IGRAPH_CHECK(igraph_strvector_set(strvec, l, graphmlrec->default_value.as_string));
            }
        } else if (rec->type == IGRAPH_ATTRIBUTE_BOOLEAN) {
            igraph_vector_bool_t *boolvec = (igraph_vector_bool_t*)rec->value;
            igraph_integer_t origsize = igraph_vector_bool_size(boolvec);
            IGRAPH_CHECK(igraph_vector_bool_resize(boolvec, 1));
            for (l = origsize; l < 1; l++) {
                VECTOR(*boolvec)[l] = graphmlrec->default_value.as_boolean;
            }
        } else if (rec->type == IGRAPH_ATTRIBUTE_UNSPECIFIED) {
            continue; /* skipped attribute */
        }
        igraph_vector_ptr_push_back(&gattr, rec); /* reserved */
    }

    IGRAPH_CHECK(igraph_empty_attrs(state->g, 0, state->edges_directed, &gattr));
    IGRAPH_FINALLY(igraph_destroy, state->g); /* because the next two lines may fail as well */
    IGRAPH_CHECK(igraph_add_vertices(state->g, igraph_trie_size(&state->node_trie), &vattr));
    IGRAPH_CHECK(igraph_add_edges(state->g, &state->edgelist, &eattr));
    IGRAPH_FINALLY_CLEAN(1); /* graph construction completed successfully */

    igraph_vector_ptr_destroy(&vattr);
    igraph_vector_ptr_destroy(&eattr);
    igraph_vector_ptr_destroy(&gattr);
    IGRAPH_FINALLY_CLEAN(3);

    return IGRAPH_SUCCESS;
}

/* See https://gnome.pages.gitlab.gnome.org/libxml2/devhelp/libxml2-parser.html#startElementNsSAX2Func */
#define XML_ATTR_LOCALNAME(it) it[0]
#define XML_ATTR_PREFIX(it) it[1]
#define XML_ATTR_URI(it) it[2]
#define XML_ATTR_VALUE_START(it) it[3]
#define XML_ATTR_VALUE_END(it) it[4]
#define XML_ATTR_VALUE_LENGTH(it) (size_t)(it[4] - it[3])
#define XML_ATTR_VALUE(it) it[3], (int)(it[4] - it[3]) /* for use in strnxxx()-style functions that take a char * and a length */
#define XML_ATTR_VALUE_PF(it) (int)(it[4] - it[3]), it[3] /* for use in printf-style function with "%.*s" */

static igraph_bool_t xmlAttrValueEqual(xmlChar** attr, const char* expected) {
    size_t expected_length = strlen(expected);
    return (
        expected_length == XML_ATTR_VALUE_LENGTH(attr) &&
        !xmlStrncmp(toXmlChar(expected), XML_ATTR_VALUE(attr))
    );
}

static igraph_error_t igraph_i_graphml_add_attribute_key(
    igraph_i_graphml_attribute_record_t** record,
    const xmlChar** attrs, int nb_attrs,
    struct igraph_i_graphml_parser_state *state
) {

    /* This function must return in three possible ways:
     *
     * - a proper newly allocated attribute record is returned in 'record' and
     *   the function returns IGRAPH_SUCCESS; the parser will process the attribute
     * - NULL is returned in 'record' and the function returns an igraph error
     *   code; the parser will handle the error
     * - NULL is returned in 'record', but the function itself returns
     *   IGRAPH_SUCCESS; the parser will skip the attribute
     *
     * The caller should be prepared to handle all three cases.
     */

    xmlChar **it;
    xmlChar *localname;
    xmlChar *xmlStr;
    igraph_trie_t *trie = NULL;
    igraph_vector_ptr_t *ptrvector = NULL;
    igraph_integer_t i, n;
    igraph_integer_t id;
    igraph_i_graphml_attribute_record_t *rec = NULL;
    igraph_bool_t skip = false;

    if (!state->successful) {
        /* Parser is already in an error state */
        goto exit;
    }

    rec = IGRAPH_CALLOC(1, igraph_i_graphml_attribute_record_t);
    if (rec == NULL) {
        IGRAPH_ERROR("Cannot allocate attribute record.", IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */
    }
    IGRAPH_FINALLY(igraph_free, rec);
    IGRAPH_FINALLY(igraph_i_graphml_attribute_record_destroy, rec);

    rec->type = I_GRAPHML_UNKNOWN_TYPE;

    for (i = 0, it = (xmlChar**)attrs; i < nb_attrs; i++, it += 5) {
        if (XML_ATTR_URI(it) != 0 &&
            !xmlStrEqual(toXmlChar(GRAPHML_NAMESPACE_URI), XML_ATTR_URI(it))) {
            continue;
        }

        localname = XML_ATTR_LOCALNAME(it);

        if (xmlStrEqual(localname, toXmlChar("id"))) {
            xmlStr = xmlStrndup(XML_ATTR_VALUE(it));
            IGRAPH_CHECK_OOM(xmlStr, "Cannot duplicate value of 'id' attribute.");
            rec->id = fromXmlChar(xmlStr);
            xmlStr = NULL;
        } else if (xmlStrEqual(localname, toXmlChar("attr.name"))) {
            xmlStr = xmlStrndup(XML_ATTR_VALUE(it));
            IGRAPH_CHECK_OOM(xmlStr, "Cannot duplicate value of 'attr.name' attribute.");
            rec->record.name = fromXmlChar(xmlStr);
            xmlStr = NULL;
        } else if (xmlStrEqual(localname, toXmlChar("attr.type"))) {
            if (xmlAttrValueEqual(it, "boolean")) {
                rec->type = I_GRAPHML_BOOLEAN;
                rec->record.type = IGRAPH_ATTRIBUTE_BOOLEAN;
                rec->default_value.as_boolean = 0;
            } else if (xmlAttrValueEqual(it, "string")) {
                char *str = strdup("");
                IGRAPH_CHECK_OOM(str, "Cannot allocate new empty string.");
                rec->type = I_GRAPHML_STRING;
                rec->record.type = IGRAPH_ATTRIBUTE_STRING;
                rec->default_value.as_string = str;
            } else if (xmlAttrValueEqual(it, "float")) {
                rec->type = I_GRAPHML_FLOAT;
                rec->record.type = IGRAPH_ATTRIBUTE_NUMERIC;
                rec->default_value.as_numeric = IGRAPH_NAN;
            } else if (xmlAttrValueEqual(it, "double")) {
                rec->type = I_GRAPHML_DOUBLE;
                rec->record.type = IGRAPH_ATTRIBUTE_NUMERIC;
                rec->default_value.as_numeric = IGRAPH_NAN;
            } else if (xmlAttrValueEqual(it, "int")) {
                rec->type = I_GRAPHML_INTEGER;
                rec->record.type = IGRAPH_ATTRIBUTE_NUMERIC;
                rec->default_value.as_numeric = IGRAPH_NAN;
            } else if (xmlAttrValueEqual(it, "long")) {
                rec->type = I_GRAPHML_LONG;
                rec->record.type = IGRAPH_ATTRIBUTE_NUMERIC;
                rec->default_value.as_numeric = IGRAPH_NAN;
            } else {
                IGRAPH_ERRORF("Unknown attribute type '%.*s'.", IGRAPH_PARSEERROR,
                              XML_ATTR_VALUE_PF(it));
            }
        } else if (xmlStrEqual(*it, toXmlChar("for"))) {
            /* graph, vertex or edge attribute? */
            if (xmlAttrValueEqual(it, "graph")) {
                trie = &state->g_attr_ids;
                ptrvector = &state->g_attrs;
            } else if (xmlAttrValueEqual(it, "node")) {
                trie = &state->v_attr_ids;
                ptrvector = &state->v_attrs;
            } else if (xmlAttrValueEqual(it, "edge")) {
                trie = &state->e_attr_ids;
                ptrvector = &state->e_attrs;
            } else if (xmlAttrValueEqual(it, "graphml")) {
                igraph_i_report_unhandled_attribute_target("graphml", IGRAPH_FILE_BASENAME, __LINE__);
                skip = 1;
            } else if (xmlAttrValueEqual(it, "hyperedge")) {
                igraph_i_report_unhandled_attribute_target("hyperedge", IGRAPH_FILE_BASENAME, __LINE__);
                skip = 1;
            } else if (xmlAttrValueEqual(it, "port")) {
                igraph_i_report_unhandled_attribute_target("port", IGRAPH_FILE_BASENAME, __LINE__);
                skip = 1;
            } else if (xmlAttrValueEqual(it, "endpoint")) {
                igraph_i_report_unhandled_attribute_target("endpoint", IGRAPH_FILE_BASENAME, __LINE__);
                skip = 1;
            } else if (xmlAttrValueEqual(it, "all")) {
                /* TODO: we should handle this */
                igraph_i_report_unhandled_attribute_target("all", IGRAPH_FILE_BASENAME, __LINE__);
                skip = 1;
            } else {
                IGRAPH_ERRORF("Unknown value '%.*s' in the 'for' attribute of a <key> tag.", IGRAPH_PARSEERROR,
                              XML_ATTR_VALUE_PF(it));
            }
        }
    }

    /* throw an error if there is no ID; this is a clear violation of the GraphML DTD */
    if (rec->id == NULL) {
        IGRAPH_ERROR("Found <key> tag with no 'id' attribute.", IGRAPH_PARSEERROR);
    }

    /* throw an error if the ID is an empty string; this is also a clear violation of the GraphML DTD */
    if (*(rec->id) == '\0') {
        IGRAPH_ERROR("Found <key> tag with an empty 'id' attribute.", IGRAPH_PARSEERROR);
    }

    /* in case of a missing attr.name attribute, use the id as the attribute name */
    if (rec->record.name == NULL) {
        rec->record.name = strdup(rec->id);
        IGRAPH_CHECK_OOM(rec->record.name, "Cannot duplicate attribute ID as name.");
    }

    /* if the attribute type is missing, ignore the attribute with a warning */
    if (!skip && rec->type == I_GRAPHML_UNKNOWN_TYPE) {
        IGRAPH_WARNINGF("Ignoring <key id=\"%s\"> because of a missing 'attr.type' attribute.", rec->id);
        skip = 1;
    }

    /* if the value of the 'for' attribute was unknown, throw an error */
    if (!skip && trie == 0) {
        IGRAPH_ERROR("Missing 'for' attribute in a <key> tag.", IGRAPH_PARSEERROR);
    }

    /* If attribute is skipped, proceed according to the type of the associated graph element. */
    if (skip) {
        if (trie == 0) {
            /* Attribute was skipped because it is not for a node, edge or the graph.
             * Free everything and return. */
            if (rec) {
                igraph_i_graphml_attribute_record_destroy(rec);
                IGRAPH_FREE(rec);
            }
            IGRAPH_FINALLY_CLEAN(2);
            goto exit;
        } else {
            /* If the skipped attribute was for a supported graph element, we add it
             * as "UNSPECIFIED" so that we can avoid reporting "unknown attribute" warnings
             * later. */
            rec->record.type = IGRAPH_ATTRIBUTE_UNSPECIFIED;
        }
    }

    /* check if we have already seen this ID */
    IGRAPH_CHECK(igraph_trie_check(trie, rec->id, &id));
    if (id >= 0) {
        IGRAPH_ERRORF("Duplicate attribute ID found: '%s'.", IGRAPH_PARSEERROR, rec->id);
    }

    /* check if we have already seen this attribute name */
    n = igraph_vector_ptr_size(ptrvector);
    for (i = 0; i < n; i++) {
        if (!strcmp(
            rec->record.name,
            ((igraph_i_graphml_attribute_record_t*) igraph_vector_ptr_get(ptrvector, i))->record.name
        )) {
            IGRAPH_ERRORF(
                "Duplicate attribute name found: '%s' (for <key id='%s'>).",
                IGRAPH_PARSEERROR, rec->record.name, rec->id
            );
        }
    }

    /* add to trie, attributes */
    IGRAPH_CHECK(igraph_trie_get(trie, rec->id, &id));
    IGRAPH_CHECK(igraph_vector_ptr_push_back(ptrvector, rec));

    /* Ownership of 'rec' is now taken by ptrvector so we are not responsible
     * for destroying and freeing it any more */
    IGRAPH_FINALLY_CLEAN(2);

    /* create the attribute values */
    switch (rec->record.type) {
        igraph_vector_t *vec;
        igraph_vector_bool_t *boolvec;
        igraph_strvector_t *strvec;
    case IGRAPH_ATTRIBUTE_BOOLEAN:
        boolvec = IGRAPH_CALLOC(1, igraph_vector_bool_t);
        IGRAPH_CHECK_OOM(boolvec, "Cannot allocate value vector for Boolean attribute.");
        IGRAPH_FINALLY(igraph_free, boolvec);
        IGRAPH_CHECK(igraph_vector_bool_init(boolvec, 0));
        rec->record.value = boolvec;
        IGRAPH_FINALLY_CLEAN(1);
        break;
    case IGRAPH_ATTRIBUTE_NUMERIC:
        vec = IGRAPH_CALLOC(1, igraph_vector_t);
        IGRAPH_CHECK_OOM(vec, "Cannot allocate value vector for numeric attribute.");
        IGRAPH_FINALLY(igraph_free, vec);
        IGRAPH_CHECK(igraph_vector_init(vec, 0));
        rec->record.value = vec;
        IGRAPH_FINALLY_CLEAN(1);
        break;
    case IGRAPH_ATTRIBUTE_STRING:
        strvec = IGRAPH_CALLOC(1, igraph_strvector_t);
        IGRAPH_CHECK_OOM(strvec, "Cannot allocate value vector for string attribute.");
        IGRAPH_FINALLY(igraph_free, strvec);
        IGRAPH_CHECK(igraph_strvector_init(strvec, 0));
        rec->record.value = strvec;
        IGRAPH_FINALLY_CLEAN(1);
        break;
    case IGRAPH_ATTRIBUTE_UNSPECIFIED:
        rec->record.value = NULL;
        break;
    default:
        IGRAPH_FATAL("Unexpected attribute type.");
    }

exit:
    *record = rec;
    return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_i_graphml_attribute_data_setup(
    struct igraph_i_graphml_parser_state *state, const xmlChar **attrs,
    int nb_attrs, igraph_attribute_elemtype_t type
) {
    xmlChar **it;
    int i;

    if (!state->successful) {
        return IGRAPH_SUCCESS;
    }

    for (i = 0, it = (xmlChar**)attrs; i < nb_attrs; i++, it += 5) {
        if (XML_ATTR_URI(it) != 0 &&
            !xmlStrEqual(toXmlChar(GRAPHML_NAMESPACE_URI), XML_ATTR_URI(it))) {
            continue;
        }

        if (xmlStrEqual(*it, toXmlChar("key"))) {
            if (state->data_key) {
                xmlFree((void *) state->data_key);
                state->data_key = NULL;
            }
            state->data_key = xmlStrndup(XML_ATTR_VALUE(it));
            if (state->data_key == 0) {
                return IGRAPH_ENOMEM; /* LCOV_EXCL_LINE */
            }
            igraph_vector_char_clear(&state->data_char);
            state->data_type = type;
        } else {
            /* ignore */
        }
    }

    return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_i_graphml_append_to_data_char(
    struct igraph_i_graphml_parser_state *state, const xmlChar *data, int len
) {
    if (!state->successful) {
        return IGRAPH_SUCCESS;
    }

    /* vector_push_back() minimizes reallocations by doubling the size of the buffer,
     * while vector_append() would only allocate as much additional memory as needed. */
    IGRAPH_STATIC_ASSERT(sizeof(char) == sizeof(xmlChar));
    for (int i=0; i < len; i++) {
        IGRAPH_CHECK(igraph_vector_char_push_back(&state->data_char, data[i]));
    }

    return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_i_graphml_attribute_data_finish(struct igraph_i_graphml_parser_state *state) {
    const char *key = fromXmlChar(state->data_key);
    igraph_attribute_elemtype_t type = state->data_type;
    igraph_trie_t *trie = NULL;
    igraph_vector_ptr_t *ptrvector = NULL;
    igraph_i_graphml_attribute_record_t *graphmlrec;
    igraph_attribute_record_t *rec;
    igraph_integer_t recid, id = 0;
    igraph_error_t result = IGRAPH_SUCCESS;

    switch (type) {
    case IGRAPH_ATTRIBUTE_GRAPH:
        trie = &state->g_attr_ids;
        ptrvector = &state->g_attrs;
        id = 0;
        break;
    case IGRAPH_ATTRIBUTE_VERTEX:
        trie = &state->v_attr_ids;
        ptrvector = &state->v_attrs;
        id = state->act_node;
        break;
    case IGRAPH_ATTRIBUTE_EDGE:
        trie = &state->e_attr_ids;
        ptrvector = &state->e_attrs;
        id = igraph_vector_int_size(&state->edgelist) / 2 - 1; /* hack */
        break;
    default:
        IGRAPH_FATAL("Unexpected attribute element type.");
    }

    if (key == 0) {
        /* no key specified, issue a warning */
        IGRAPH_WARNING("Missing attribute key in a <data> tag, ignoring attribute.");
        goto exit;
    }

    IGRAPH_CHECK(igraph_trie_check(trie, key, &recid));
    if (recid < 0) {
        /* no such attribute key, issue a warning */
        IGRAPH_WARNINGF(
            "Unknown attribute key '%s' in a <data> tag, ignoring attribute.",
            key
        );
        goto exit;
    }

    graphmlrec = VECTOR(*ptrvector)[recid];
    rec = &graphmlrec->record;

    switch (rec->type) {
        igraph_vector_bool_t *boolvec;
        igraph_vector_t *vec;
        igraph_strvector_t *strvec;
        igraph_integer_t s, i;
        const char* strvalue;

    case IGRAPH_ATTRIBUTE_BOOLEAN:
        boolvec = (igraph_vector_bool_t *)rec->value;
        s = igraph_vector_bool_size(boolvec);
        if (id >= s) {
            IGRAPH_CHECK(igraph_vector_bool_resize(boolvec, id + 1));
            for (i = s; i < id; i++) {
                VECTOR(*boolvec)[i] = graphmlrec->default_value.as_boolean;
            }
        }

        /* Add null terminator */
        IGRAPH_CHECK(igraph_vector_char_push_back(&state->data_char, '\x00'));
        IGRAPH_CHECK(igraph_i_graphml_parse_boolean(
            VECTOR(state->data_char), VECTOR(*boolvec) + id,  graphmlrec->default_value.as_boolean
        ));
        break;

    case IGRAPH_ATTRIBUTE_NUMERIC:
        vec = (igraph_vector_t *)rec->value;
        s = igraph_vector_size(vec);
        if (id >= s) {
            IGRAPH_CHECK(igraph_vector_resize(vec, id + 1));
            for (i = s; i < id; i++) {
                VECTOR(*vec)[i] = graphmlrec->default_value.as_numeric;
            }
        }

        /* Add null terminator */
        IGRAPH_CHECK(igraph_vector_char_push_back(&state->data_char, '\x00'));
        IGRAPH_CHECK(igraph_i_graphml_parse_numeric(
            VECTOR(state->data_char), VECTOR(*vec) + id,
            graphmlrec->default_value.as_numeric
        ));
        break;

    case IGRAPH_ATTRIBUTE_STRING:
        strvec = (igraph_strvector_t *)rec->value;
        s = igraph_strvector_size(strvec);
        if (id >= s) {
            IGRAPH_CHECK(igraph_strvector_resize(strvec, id + 1));
            strvalue = graphmlrec->default_value.as_string;
            for (i = s; i < id; i++) {
                IGRAPH_CHECK(igraph_strvector_set(strvec, i, strvalue));
            }
        }
        if (igraph_vector_char_size(&state->data_char) > 0) {
            /* Ensure that the vector ends with a null terminator */
            IGRAPH_CHECK(igraph_vector_char_push_back(&state->data_char, '\x00'));
            strvalue = VECTOR(state->data_char);
        } else {
            strvalue = graphmlrec->default_value.as_string;
        }
        IGRAPH_CHECK(igraph_strvector_set(strvec, id, strvalue));
        break;

    case IGRAPH_ATTRIBUTE_UNSPECIFIED:
        break;

    default:
        IGRAPH_FATAL("Unexpected attribute type.");
    }

exit:
    igraph_vector_char_clear(&state->data_char);

    return result;
}

static igraph_error_t igraph_i_graphml_attribute_default_value_finish(struct igraph_i_graphml_parser_state *state) {
    igraph_i_graphml_attribute_record_t *graphmlrec = state->current_attr_record;
    igraph_error_t result = IGRAPH_SUCCESS;
    char* str = 0;

    IGRAPH_ASSERT(state->current_attr_record != NULL);

    if (igraph_vector_char_size(&state->data_char) == 0) {
        return IGRAPH_SUCCESS;
    }

    switch (graphmlrec->record.type) {
    case IGRAPH_ATTRIBUTE_BOOLEAN:
        /* Add null terminator */
        IGRAPH_CHECK(igraph_vector_char_push_back(&state->data_char, '\x00'));
        IGRAPH_CHECK(igraph_i_graphml_parse_boolean(
            VECTOR(state->data_char), &graphmlrec->default_value.as_boolean, 0
        ));
        break;
    case IGRAPH_ATTRIBUTE_NUMERIC:
        /* Add null terminator */
        IGRAPH_CHECK(igraph_vector_char_push_back(&state->data_char, '\x00'));
        IGRAPH_CHECK(igraph_i_graphml_parse_numeric(
            VECTOR(state->data_char), &graphmlrec->default_value.as_numeric, IGRAPH_NAN
        ));
        break;
    case IGRAPH_ATTRIBUTE_STRING:
        /* Add null terminator */
        IGRAPH_CHECK(igraph_vector_char_push_back(&state->data_char, '\x00'));
        str = strdup(VECTOR(state->data_char));
        IGRAPH_CHECK_OOM(str, "Cannot allocate memory for string attribute.");

        if (graphmlrec->default_value.as_string != 0) {
            IGRAPH_FREE(graphmlrec->default_value.as_string);
        }
        graphmlrec->default_value.as_string = str;
        str = NULL;
        break;
    case IGRAPH_ATTRIBUTE_UNSPECIFIED:
        break;
    default:
        IGRAPH_FATAL("Unexpected attribute type.");
    }

    igraph_vector_char_clear(&state->data_char);

    return result;
}

static igraph_error_t igraph_i_graphml_sax_handler_start_element_ns_inner(
        struct igraph_i_graphml_parser_state* state, const xmlChar* localname, const xmlChar* prefix,
        const xmlChar* uri, int nb_namespaces, const xmlChar** namespaces,
        int nb_attributes, int nb_defaulted, const xmlChar** attributes) {
    xmlChar** it;
    xmlChar* attr_value = 0;
    igraph_integer_t id1, id2;
    int i;
    igraph_bool_t tag_is_unknown = false;

    IGRAPH_UNUSED(prefix);
    IGRAPH_UNUSED(nb_namespaces);
    IGRAPH_UNUSED(namespaces);
    IGRAPH_UNUSED(nb_defaulted);

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
        IGRAPH_CHECK(igraph_i_graphml_handle_unknown_start_tag(state));
        goto exit;
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
            IGRAPH_CHECK(igraph_i_graphml_handle_unknown_start_tag(state));
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
                        if (xmlAttrValueEqual(it, "directed")) {
                            state->edges_directed = 1;
                        } else if (xmlAttrValueEqual(it, "undirected")) {
                            state->edges_directed = 0;
                        }
                    }
                }
            }
            state->index--;
        } else if (xmlStrEqual(localname, toXmlChar("key"))) {
            IGRAPH_CHECK(
                igraph_i_graphml_add_attribute_key(
                    &state->current_attr_record,
                    attributes, nb_attributes, state
                )
            );
            /* NULL is okay here for state->current_attr_record -- we should have
             * triggered an error in the parser already if we returned NULL, and
             * the rest of the code is prepared to handle NULLs */
            state->st = INSIDE_KEY;
        } else {
            IGRAPH_CHECK(igraph_i_graphml_handle_unknown_start_tag(state));
        }
        break;

    case INSIDE_KEY:
        /* If we are in the INSIDE_KEY state and we are not skipping the current
         * attribute, check for default tag */
        if (state->current_attr_record != NULL && xmlStrEqual(localname, toXmlChar("default"))) {
            state->st = INSIDE_DEFAULT;
        } else {
            IGRAPH_CHECK(igraph_i_graphml_handle_unknown_start_tag(state));
        }
        break;

    case INSIDE_DEFAULT:
        /* If we are in the INSIDE_DEFAULT state, every further tag will be unknown */
        IGRAPH_CHECK(igraph_i_graphml_handle_unknown_start_tag(state));
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
                    attr_value = xmlStrndup(XML_ATTR_VALUE(it));
                    if (attr_value == 0) {
                        IGRAPH_ERROR("Cannot copy value of edge source attribute.", IGRAPH_ENOMEM);
                    }
                    IGRAPH_FINALLY(xmlFree, attr_value);

                    IGRAPH_CHECK(igraph_trie_get(&state->node_trie, fromXmlChar(attr_value), &id1));

                    xmlFree(attr_value); attr_value = NULL;
                    IGRAPH_FINALLY_CLEAN(1);
                } else if (xmlStrEqual(*it, toXmlChar("target"))) {
                    attr_value = xmlStrndup(XML_ATTR_VALUE(it));
                    if (attr_value == 0) {
                        IGRAPH_ERROR("Cannot copy value of edge target attribute.", IGRAPH_ENOMEM);
                    }
                    IGRAPH_FINALLY(xmlFree, attr_value);

                    IGRAPH_CHECK(igraph_trie_get(&state->node_trie, fromXmlChar(attr_value), &id2));

                    xmlFree(attr_value); attr_value = NULL;
                    IGRAPH_FINALLY_CLEAN(1);
                } else if (xmlStrEqual(*it, toXmlChar("id"))) {
                    igraph_integer_t edges = igraph_vector_int_size(&state->edgelist) / 2 + 1;
                    igraph_integer_t origsize = igraph_strvector_size(&state->edgeids);

                    attr_value = xmlStrndup(XML_ATTR_VALUE(it));
                    if (attr_value == 0) {
                        IGRAPH_ERROR("Cannot copy value of edge ID attribute.", IGRAPH_ENOMEM);
                    }
                    IGRAPH_FINALLY(xmlFree, attr_value);

                    IGRAPH_CHECK(igraph_strvector_resize(&state->edgeids, edges));

                    for (; origsize < edges - 1; origsize++) {
                        IGRAPH_CHECK(igraph_strvector_set(&state->edgeids, origsize, ""));
                    }

                    IGRAPH_CHECK(igraph_strvector_set(&state->edgeids, edges - 1, fromXmlChar(attr_value)));

                    xmlFree(attr_value); attr_value = NULL;
                    IGRAPH_FINALLY_CLEAN(1);
                }
            }
            if (id1 >= 0 && id2 >= 0) {
                IGRAPH_CHECK(igraph_vector_int_push_back(&state->edgelist, id1));
                IGRAPH_CHECK(igraph_vector_int_push_back(&state->edgelist, id2));
            } else {
                IGRAPH_ERROR("Edge with missing source or target encountered.", IGRAPH_PARSEERROR);
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
                    attr_value = xmlStrndup(XML_ATTR_VALUE(it));
                    if (attr_value == 0) {
                        IGRAPH_ERROR("Cannot copy value of node ID attribute.", IGRAPH_ENOMEM);
                    }
                    IGRAPH_FINALLY(xmlFree, attr_value);

                    IGRAPH_CHECK(igraph_trie_get(&state->node_trie, fromXmlChar(attr_value), &id1));

                    xmlFree(attr_value); attr_value = NULL;
                    IGRAPH_FINALLY_CLEAN(1);
                    break;
                }
            }
            if (id1 >= 0) {
                state->act_node = id1;
            } else {
                state->act_node = -1;
                IGRAPH_ERROR("Node with missing ID encountered.", IGRAPH_PARSEERROR);
            }
            state->st = INSIDE_NODE;
        } else if (xmlStrEqual(localname, toXmlChar("data"))) {
            IGRAPH_CHECK(igraph_i_graphml_attribute_data_setup(
                state, attributes, nb_attributes, IGRAPH_ATTRIBUTE_GRAPH
            ));
            IGRAPH_CHECK(igraph_vector_int_push_back(&state->prev_state_stack, state->st));
            state->st = INSIDE_DATA;
        } else {
            IGRAPH_CHECK(igraph_i_graphml_handle_unknown_start_tag(state));
        }
        break;

    case INSIDE_NODE:
        if (xmlStrEqual(localname, toXmlChar("data"))) {
            IGRAPH_CHECK(igraph_i_graphml_attribute_data_setup(
                state, attributes, nb_attributes, IGRAPH_ATTRIBUTE_VERTEX
            ));
            IGRAPH_CHECK(igraph_vector_int_push_back(&state->prev_state_stack, state->st));
            state->st = INSIDE_DATA;
        } else {
            IGRAPH_CHECK(igraph_i_graphml_handle_unknown_start_tag(state));
        }
        break;

    case INSIDE_EDGE:
        if (xmlStrEqual(localname, toXmlChar("data"))) {
            IGRAPH_CHECK(igraph_i_graphml_attribute_data_setup(
                state, attributes, nb_attributes, IGRAPH_ATTRIBUTE_EDGE
            ));
            IGRAPH_CHECK(igraph_vector_int_push_back(&state->prev_state_stack, state->st));
            state->st = INSIDE_DATA;
        } else {
            IGRAPH_CHECK(igraph_i_graphml_handle_unknown_start_tag(state));
        }
        break;

    case INSIDE_DATA:
        /* We do not expect any new tags within a <data> tag */
        IGRAPH_CHECK(igraph_i_graphml_handle_unknown_start_tag(state));
        break;

    case UNKNOWN:
        IGRAPH_CHECK(igraph_i_graphml_handle_unknown_start_tag(state));
        break;

    case FINISH:
        break;

    default:
        IGRAPH_FATALF("Unexpected GraphML reader state %d.", (int) state->st);
    }

exit:
    return IGRAPH_SUCCESS;
}

static void igraph_i_graphml_sax_handler_start_element_ns(
        void *state0, const xmlChar* localname, const xmlChar* prefix,
        const xmlChar* uri, int nb_namespaces, const xmlChar** namespaces,
        int nb_attributes, int nb_defaulted, const xmlChar** attributes) {
    struct igraph_i_graphml_parser_state *state =
        (struct igraph_i_graphml_parser_state*)state0;
    igraph_error_t result;

    if (!state->successful) {
        return;
    }

    result = igraph_i_graphml_sax_handler_start_element_ns_inner(
        state, localname, prefix, uri, nb_namespaces, namespaces,
        nb_attributes, nb_defaulted, attributes
    );

    if (result != IGRAPH_SUCCESS) {
        RETURN_GRAPHML_PARSE_ERROR_WITH_CODE(state, "Cannot parse GraphML file.", result);
    }
}

static igraph_error_t igraph_i_graphml_sax_handler_end_element_ns_inner(
        struct igraph_i_graphml_parser_state* state,
        const xmlChar* localname, const xmlChar* prefix,
        const xmlChar* uri
) {
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
        state->current_attr_record = NULL;
        state->st = INSIDE_GRAPHML;
        break;

    case INSIDE_DEFAULT:
        IGRAPH_CHECK(igraph_i_graphml_attribute_default_value_finish(state));
        state->st = INSIDE_KEY;
        break;

    case INSIDE_NODE:
        state->st = INSIDE_GRAPH;
        break;

    case INSIDE_EDGE:
        state->st = INSIDE_GRAPH;
        break;

    case INSIDE_DATA:
        IGRAPH_CHECK(igraph_i_graphml_attribute_data_finish(state));
        IGRAPH_ASSERT(!igraph_vector_int_empty(&state->prev_state_stack));
        state->st = (igraph_i_graphml_parser_state_index_t) igraph_vector_int_pop_back(&state->prev_state_stack);
        break;

    case UNKNOWN:
        state->unknown_depth--;
        if (!state->unknown_depth) {
            IGRAPH_ASSERT(!igraph_vector_int_empty(&state->prev_state_stack));
            state->st = (igraph_i_graphml_parser_state_index_t) igraph_vector_int_pop_back(&state->prev_state_stack);
        }
        break;

    case FINISH:
        break;

    default:
        IGRAPH_FATALF("Unexpected GraphML reader state %d.", (int) state->st);
    }

    return IGRAPH_SUCCESS;
}

static void igraph_i_graphml_sax_handler_end_element_ns(
        void *state0,
        const xmlChar* localname, const xmlChar* prefix,
        const xmlChar* uri) {
    struct igraph_i_graphml_parser_state *state =
        (struct igraph_i_graphml_parser_state*)state0;
    igraph_error_t result;

    if (!state->successful) {
        return;
    }

    result = igraph_i_graphml_sax_handler_end_element_ns_inner(
        state, localname, prefix, uri
    );

    if (result != IGRAPH_SUCCESS) {
        RETURN_GRAPHML_PARSE_ERROR_WITH_CODE(state, "Cannot parse GraphML file.", result);
    }
}

static void igraph_i_graphml_sax_handler_chars(void* state0, const xmlChar* ch, int len) {
    struct igraph_i_graphml_parser_state *state =
        (struct igraph_i_graphml_parser_state*)state0;
    igraph_error_t result = IGRAPH_SUCCESS;

    if (!state->successful) {
        return;
    }

    switch (state->st) {
    case INSIDE_KEY:
        break;

    case INSIDE_DATA:
    case INSIDE_DEFAULT:
        result = igraph_i_graphml_append_to_data_char(state, ch, len);
        break;

    default:
        /* just ignore it */
        break;
    }

    if (result != IGRAPH_SUCCESS) {
        RETURN_GRAPHML_PARSE_ERROR_WITH_CODE(state, "Cannot parse GraphML file.", result);
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
    /* endDocument = */ 0,
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

#endif // HAVE_LIBXML == 1

#define IS_FORBIDDEN_CONTROL_CHAR(x) ((x) < ' ' && (x) != '\t' && (x) != '\r' && (x) != '\n')

static igraph_error_t igraph_i_xml_escape(const char* src, char** dest) {
    igraph_integer_t destlen = 0;
    const char *s;
    char *d;
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
            IGRAPH_ERRORF("Forbidden control character 0x%02X found in igraph_i_xml_escape.", IGRAPH_EINVAL, ch);
        }
    }
    *dest = IGRAPH_CALLOC(destlen + 1, char);
    if (!*dest) {
        IGRAPH_ERROR("Not enough memory.", IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */
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
    return IGRAPH_SUCCESS;
}

#if HAVE_LIBXML == 1
static void igraph_i_libxml_generic_error_handler(void* ctx, const char* msg, ...) {
    struct igraph_i_graphml_parser_state* state = (struct igraph_i_graphml_parser_state*) ctx;
    va_list args;
    va_start(args, msg);
    igraph_i_graphml_parser_state_set_error_from_varargs(state, msg, args);
    va_end(args);
}

#if LIBXML_VERSION < 21200
static void igraph_i_libxml_structured_error_handler(void* ctx, xmlError *error) {
#else
static void igraph_i_libxml_structured_error_handler(void* ctx, const xmlError *error) {
#endif
    struct igraph_i_graphml_parser_state* state = (struct igraph_i_graphml_parser_state*) ctx;
    igraph_i_graphml_parser_state_set_error_from_xmlerror(state, error);
}
#endif // HAVE_LIBXML == 1

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
 * is attached, see \ref igraph_set_attribute_table(). String attrribute values
 * are returned in UTF-8 encoding.
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
igraph_error_t igraph_read_graph_graphml(igraph_t *graph, FILE *instream, igraph_integer_t index) {

#if HAVE_LIBXML == 1
    xmlParserCtxtPtr ctxt;
    xmlGenericErrorFunc libxml_old_generic_error_handler;
    void* libxml_old_generic_error_context;
    xmlStructuredErrorFunc libxml_old_structured_error_handler;
    void* libxml_old_structured_error_context;
    xmlDocPtr doc;

    struct igraph_i_graphml_parser_state state;
    int res;
    char buffer[4096];
    igraph_bool_t parsing_successful;
    char* error_message;

    if (index < 0) {
        IGRAPH_ERROR("Graph index must be non-negative.", IGRAPH_EINVAL);
    }

    xmlInitParser();

    IGRAPH_CHECK(igraph_i_graphml_parser_state_init(&state, graph, index));
    IGRAPH_FINALLY(igraph_i_graphml_parser_state_destroy, &state);

    /* Create a progressive parser context and use the first 4K to detect the
     * encoding */
    res = (int) fread(buffer, 1, sizeof(buffer), instream);
    if (res < (int) sizeof buffer && !feof(instream)) {
        IGRAPH_ERROR("IO error while reading GraphML data.", IGRAPH_EFILE);
    }

    /* Retrieve the current libxml2 error handlers and temporarily replace them
     * with ones that do not print anything to stdout/stderr */
    libxml_old_generic_error_handler = xmlGenericError;
    libxml_old_generic_error_context = xmlGenericErrorContext;
    libxml_old_structured_error_handler = xmlStructuredError;
    libxml_old_structured_error_context = xmlStructuredErrorContext;
    xmlSetGenericErrorFunc(&state, &igraph_i_libxml_generic_error_handler);
    xmlSetStructuredErrorFunc(&state, &igraph_i_libxml_structured_error_handler);

    /* Okay, parsing will start now. The parser might do things that eventually
     * trigger the igraph error handler, but we want the parser state to
     * survive whatever happens here. So, we put a barrier on the FINALLY stack
     * that prevents IGRAPH_ERROR() from freeing the parser state, and then we
     * do this ourselves when needed */
    IGRAPH_FINALLY_ENTER();
    {
        ctxt = xmlCreatePushParserCtxt(&igraph_i_graphml_sax_handler,
                                       &state,
                                       buffer,
                                       res,
                                       NULL);
        if (ctxt) {
            if (xmlCtxtUseOptions(ctxt,
                                  XML_PARSE_NOBLANKS |
                                  XML_PARSE_NONET | XML_PARSE_NSCLEAN |
                                  XML_PARSE_NOCDATA | XML_PARSE_HUGE
                                  )) {
                xmlFreeParserCtxt(ctxt);
                ctxt = NULL;
            }
        }

        /* Do the parsing */
        if (ctxt) {
            while ((res = (int) fread(buffer, 1, sizeof buffer, instream)) > 0) {
                xmlParseChunk(ctxt, buffer, res, 0);
                if (!state.successful) {
                    break;
                }
                IGRAPH_ALLOW_INTERRUPTION();
            }
            xmlParseChunk(ctxt, buffer, res, 1);
        }
    }
    IGRAPH_FINALLY_EXIT();

    /* Restore error handlers */
    xmlSetGenericErrorFunc(libxml_old_generic_error_context, libxml_old_generic_error_handler);
    xmlSetStructuredErrorFunc(libxml_old_structured_error_context, libxml_old_structured_error_handler);

    /* Free the context */
    if (ctxt) {
        doc = ctxt->myDoc;
        xmlFreeParserCtxt(ctxt);
        if (doc) {
            /* In theory this should not be necessary, but it looks like certain malformed
             * GraphML files leave a partially-parsed doc in memory */
            xmlFreeDoc(doc);
        }
    } else {
        /* We could not create the context earlier so no parsing was done */
        IGRAPH_ERROR("Cannot create XML parser context.", IGRAPH_FAILURE);
    }

    /* Extract the error message from the parser state (if any), and make a
     * copy so we can safely destroy the parser state before triggering the
     * error */
    parsing_successful = state.successful;
    error_message = parsing_successful || state.error_message == NULL ? NULL : strdup(state.error_message);

    /* ...and we can also put the error message pointer on the FINALLY stack */
    if (error_message != NULL) {
        IGRAPH_FINALLY(igraph_free, error_message);
    }

    /* Trigger the stored error if needed */
    if (!parsing_successful) {
        if (error_message != NULL) {
            size_t len = strlen(error_message);
            if (error_message[len-1] == '\n') {
                error_message[len-1] = '\0';
            }
            IGRAPH_ERROR(error_message, IGRAPH_PARSEERROR);
        } else {
            IGRAPH_ERROR("Malformed GraphML file.", IGRAPH_PARSEERROR);
        }
    }

    /* Did we actually manage to reach the graph to be parsed, given its index?
     * If not, that's an error as well. */
    if (state.index >= 0) {
        IGRAPH_ERROR("Graph index was too large.", IGRAPH_EINVAL);
    }

    /* Okay, everything seems good. We can now take the parser state and
     * construct our graph from the data gathered during the parsing */
    IGRAPH_CHECK(igraph_i_graphml_parser_state_finish_parsing(&state));

    igraph_i_graphml_parser_state_destroy(&state);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
#else // HAVE_LIBXML == 1
    IGRAPH_UNUSED(graph);
    IGRAPH_UNUSED(instream);
    IGRAPH_UNUSED(index);

    IGRAPH_ERROR("GraphML support is disabled.", IGRAPH_UNIMPLEMENTED);
#endif // HAVE_LIBXML == 1
}

/**
 * \ingroup loadsave
 * \function igraph_write_graph_graphml
 * \brief Writes the graph to a file in GraphML format.
 *
 * GraphML is an XML-based file format for representing various types of
 * graphs. See the GraphML Primer (http://graphml.graphdrawing.org/primer/graphml-primer.html)
 * for detailed format description.
 *
 * </para><para>
 * When a numerical attribute value is NaN, it will be omitted from the file.
 *
 * </para><para>
 * This function assumes that non-ASCII characters in attribute names and string
 * attribute values are UTF-8 encoded. If this is not the case, the resulting
 * XML file will be invalid.
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
igraph_error_t igraph_write_graph_graphml(const igraph_t *graph, FILE *outstream,
                               igraph_bool_t prefixattr) {
    int ret;
    igraph_integer_t l, vc;
    igraph_eit_t it;
    igraph_strvector_t gnames, vnames, enames;
    igraph_vector_int_t gtypes, vtypes, etypes;
    igraph_integer_t i;
    igraph_vector_t numv;
    igraph_strvector_t strv;
    igraph_vector_bool_t boolv;
    const char *gprefix = prefixattr ? "g_" : "";
    const char *vprefix = prefixattr ? "v_" : "";
    const char *eprefix = prefixattr ? "e_" : "";

    ret = fprintf(outstream, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
    if (ret < 0) {
        IGRAPH_ERROR("Write failed.", IGRAPH_EFILE);
    }
    ret = fprintf(outstream, "<graphml xmlns=\"%s\"\n", GRAPHML_NAMESPACE_URI);
    if (ret < 0) {
        IGRAPH_ERROR("Write failed.", IGRAPH_EFILE);
    }
    ret = fprintf(outstream, "         xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"\n");
    if (ret < 0) {
        IGRAPH_ERROR("Write failed.", IGRAPH_EFILE);
    }
    ret = fprintf(outstream, "         xsi:schemaLocation=\"%s\n", GRAPHML_NAMESPACE_URI);
    if (ret < 0) {
        IGRAPH_ERROR("Write failed.", IGRAPH_EFILE);
    }
    ret = fprintf(outstream, "         %s/1.0/graphml.xsd\">\n", GRAPHML_NAMESPACE_URI);
    if (ret < 0) {
        IGRAPH_ERROR("Write failed.", IGRAPH_EFILE);
    }
    ret = fprintf(outstream, "<!-- Created by igraph -->\n");
    if (ret < 0) {
        IGRAPH_ERROR("Write failed.", IGRAPH_EFILE);
    }

    /* dump the <key> elements if any */

    IGRAPH_VECTOR_INIT_FINALLY(&numv, 1);
    IGRAPH_STRVECTOR_INIT_FINALLY(&strv, 1);
    IGRAPH_VECTOR_BOOL_INIT_FINALLY(&boolv, 1);

    IGRAPH_STRVECTOR_INIT_FINALLY(&gnames, 0);
    IGRAPH_STRVECTOR_INIT_FINALLY(&vnames, 0);
    IGRAPH_STRVECTOR_INIT_FINALLY(&enames, 0);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&gtypes, 0);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&vtypes, 0);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&etypes, 0);
    igraph_i_attribute_get_info(graph,
                                &gnames, &gtypes,
                                &vnames, &vtypes,
                                &enames, &etypes);

    /* graph attributes */
    for (i = 0; i < igraph_vector_int_size(&gtypes); i++) {
        const char *name; char *name_escaped;
        name = igraph_strvector_get(&gnames, i);
        IGRAPH_CHECK(igraph_i_xml_escape(name, &name_escaped));
        IGRAPH_FINALLY(igraph_free, name_escaped);
        if (VECTOR(gtypes)[i] == IGRAPH_ATTRIBUTE_STRING) {
            ret = fprintf(outstream, "  <key id=\"%s%s\" for=\"graph\" attr.name=\"%s\" attr.type=\"string\"/>\n", gprefix, name_escaped, name_escaped);
            if (ret < 0) {
                IGRAPH_ERROR("Write failed.", IGRAPH_EFILE);
            }
        } else if (VECTOR(gtypes)[i] == IGRAPH_ATTRIBUTE_NUMERIC) {
            ret = fprintf(outstream, "  <key id=\"%s%s\" for=\"graph\" attr.name=\"%s\" attr.type=\"double\"/>\n", gprefix, name_escaped, name_escaped);
            if (ret < 0) {
                IGRAPH_ERROR("Write failed.", IGRAPH_EFILE);
            }
        } else if (VECTOR(gtypes)[i] == IGRAPH_ATTRIBUTE_BOOLEAN) {
            ret = fprintf(outstream, "  <key id=\"%s%s\" for=\"graph\" attr.name=\"%s\" attr.type=\"boolean\"/>\n", gprefix, name_escaped, name_escaped);
            if (ret < 0) {
                IGRAPH_ERROR("Write failed.", IGRAPH_EFILE);
            }
        }
        IGRAPH_FREE(name_escaped);
        IGRAPH_FINALLY_CLEAN(1);
    }

    /* vertex attributes */
    for (i = 0; i < igraph_vector_int_size(&vtypes); i++) {
        const char *name; char *name_escaped;
        name = igraph_strvector_get(&vnames, i);
        IGRAPH_CHECK(igraph_i_xml_escape(name, &name_escaped));
        IGRAPH_FINALLY(igraph_free, name_escaped);
        if (VECTOR(vtypes)[i] == IGRAPH_ATTRIBUTE_STRING) {
            ret = fprintf(outstream, "  <key id=\"%s%s\" for=\"node\" attr.name=\"%s\" attr.type=\"string\"/>\n", vprefix, name_escaped, name_escaped);
            if (ret < 0) {
                IGRAPH_ERROR("Write failed.", IGRAPH_EFILE);
            }
        } else if (VECTOR(vtypes)[i] == IGRAPH_ATTRIBUTE_NUMERIC) {
            ret = fprintf(outstream, "  <key id=\"%s%s\" for=\"node\" attr.name=\"%s\" attr.type=\"double\"/>\n", vprefix, name_escaped, name_escaped);
            if (ret < 0) {
                IGRAPH_ERROR("Write failed.", IGRAPH_EFILE);
            }
        } else if (VECTOR(vtypes)[i] == IGRAPH_ATTRIBUTE_BOOLEAN) {
            ret = fprintf(outstream, "  <key id=\"%s%s\" for=\"node\" attr.name=\"%s\" attr.type=\"boolean\"/>\n", vprefix, name_escaped, name_escaped);
            if (ret < 0) {
                IGRAPH_ERROR("Write failed.", IGRAPH_EFILE);
            }
        }
        IGRAPH_FREE(name_escaped);
        IGRAPH_FINALLY_CLEAN(1);
    }

    /* edge attributes */
    for (i = 0; i < igraph_vector_int_size(&etypes); i++) {
        const char *name; char *name_escaped;
        name = igraph_strvector_get(&enames, i);
        IGRAPH_CHECK(igraph_i_xml_escape(name, &name_escaped));
        IGRAPH_FINALLY(igraph_free, name_escaped);
        if (VECTOR(etypes)[i] == IGRAPH_ATTRIBUTE_STRING) {
            ret = fprintf(outstream, "  <key id=\"%s%s\" for=\"edge\" attr.name=\"%s\" attr.type=\"string\"/>\n", eprefix, name_escaped, name_escaped);
            if (ret < 0) {
                IGRAPH_ERROR("Write failed.", IGRAPH_EFILE);
            }
        } else if (VECTOR(etypes)[i] == IGRAPH_ATTRIBUTE_NUMERIC) {
            ret = fprintf(outstream, "  <key id=\"%s%s\" for=\"edge\" attr.name=\"%s\" attr.type=\"double\"/>\n", eprefix, name_escaped, name_escaped);
            if (ret < 0) {
                IGRAPH_ERROR("Write failed.", IGRAPH_EFILE);
            }
        } else if (VECTOR(etypes)[i] == IGRAPH_ATTRIBUTE_BOOLEAN) {
            ret = fprintf(outstream, "  <key id=\"%s%s\" for=\"edge\" attr.name=\"%s\" attr.type=\"boolean\"/>\n", eprefix, name_escaped, name_escaped);
            if (ret < 0) {
                IGRAPH_ERROR("Write failed.", IGRAPH_EFILE);
            }
        }
        IGRAPH_FREE(name_escaped);
        IGRAPH_FINALLY_CLEAN(1);
    }

    ret = fprintf(outstream, "  <graph id=\"G\" edgedefault=\"%s\">\n", (igraph_is_directed(graph) ? "directed" : "undirected"));
    if (ret < 0) {
        IGRAPH_ERROR("Write failed.", IGRAPH_EFILE);
    }

    /* Write the graph atributes before anything else */

    for (i = 0; i < igraph_vector_int_size(&gtypes); i++) {
        const char *name; char *name_escaped;
        if (VECTOR(gtypes)[i] == IGRAPH_ATTRIBUTE_NUMERIC) {
            name = igraph_strvector_get(&gnames, i);
            IGRAPH_CHECK(igraph_i_attribute_get_numeric_graph_attr(graph, name, &numv));
            if (!isnan(VECTOR(numv)[0])) {
                IGRAPH_CHECK(igraph_i_xml_escape(name, &name_escaped));
                ret = fprintf(outstream, "    <data key=\"%s%s\">", gprefix, name_escaped);
                IGRAPH_FREE(name_escaped);
                if (ret < 0) {
                    IGRAPH_ERROR("Write failed.", IGRAPH_EFILE);
                }
                ret = igraph_real_fprintf_precise(outstream, VECTOR(numv)[0]);
                if (ret < 0) {
                    IGRAPH_ERROR("Write failed.", IGRAPH_EFILE);
                }
                ret = fprintf(outstream, "</data>\n");
                if (ret < 0) {
                    IGRAPH_ERROR("Write failed.", IGRAPH_EFILE);
                }
            }
        } else if (VECTOR(gtypes)[i] == IGRAPH_ATTRIBUTE_STRING) {
            const char *s;
            char *s_escaped;
            name = igraph_strvector_get(&gnames, i);
            IGRAPH_CHECK(igraph_i_xml_escape(name, &name_escaped));
            ret = fprintf(outstream, "    <data key=\"%s%s\">", gprefix,
                          name_escaped);
            IGRAPH_FREE(name_escaped);
            if (ret < 0) {
                IGRAPH_ERROR("Write failed.", IGRAPH_EFILE);
            }
            IGRAPH_CHECK(igraph_i_attribute_get_string_graph_attr(graph, name, &strv));
            s = igraph_strvector_get(&strv, 0);
            IGRAPH_CHECK(igraph_i_xml_escape(s, &s_escaped));
            ret = fprintf(outstream, "%s", s_escaped);
            IGRAPH_FREE(s_escaped);
            if (ret < 0) {
                IGRAPH_ERROR("Write failed.", IGRAPH_EFILE);
            }
            ret = fprintf(outstream, "</data>\n");
            if (ret < 0) {
                IGRAPH_ERROR("Write failed.", IGRAPH_EFILE);
            }
        } else if (VECTOR(gtypes)[i] == IGRAPH_ATTRIBUTE_BOOLEAN) {
            name = igraph_strvector_get(&gnames, i);
            IGRAPH_CHECK(igraph_i_attribute_get_bool_graph_attr(graph, name, &boolv));
            IGRAPH_CHECK(igraph_i_xml_escape(name, &name_escaped));
            ret = fprintf(outstream, "    <data key=\"%s%s\">%s</data>\n",
                          gprefix, name_escaped, VECTOR(boolv)[0] ? "true" : "false");
            IGRAPH_FREE(name_escaped);
            if (ret < 0) {
                IGRAPH_ERROR("Write failed.", IGRAPH_EFILE);
            }
        }
    }

    /* Let's dump the nodes first */
    vc = igraph_vcount(graph);
    for (l = 0; l < vc; l++) {
        const char *name; char *name_escaped;
        ret = fprintf(outstream, "    <node id=\"n%" IGRAPH_PRId "\">\n", l);

        if (ret < 0) {
            IGRAPH_ERROR("Write failed.", IGRAPH_EFILE);
        }

        for (i = 0; i < igraph_vector_int_size(&vtypes); i++) {
            if (VECTOR(vtypes)[i] == IGRAPH_ATTRIBUTE_NUMERIC) {
                name = igraph_strvector_get(&vnames, i);
                IGRAPH_CHECK(igraph_i_attribute_get_numeric_vertex_attr(graph, name,
                             igraph_vss_1(l), &numv));
                if (!isnan(VECTOR(numv)[0])) {
                    IGRAPH_CHECK(igraph_i_xml_escape(name, &name_escaped));
                    ret = fprintf(outstream, "      <data key=\"%s%s\">", vprefix, name_escaped);
                    IGRAPH_FREE(name_escaped);
                    if (ret < 0) {
                        IGRAPH_ERROR("Write failed.", IGRAPH_EFILE);
                    }
                    ret = igraph_real_fprintf_precise(outstream, VECTOR(numv)[0]);
                    if (ret < 0) {
                        IGRAPH_ERROR("Write failed.", IGRAPH_EFILE);
                    }
                    ret = fprintf(outstream, "</data>\n");
                    if (ret < 0) {
                        IGRAPH_ERROR("Write failed.", IGRAPH_EFILE);
                    }
                }
            } else if (VECTOR(vtypes)[i] == IGRAPH_ATTRIBUTE_STRING) {
                const char *s;
                char *s_escaped;
                name = igraph_strvector_get(&vnames, i);
                IGRAPH_CHECK(igraph_i_xml_escape(name, &name_escaped));
                ret = fprintf(outstream, "      <data key=\"%s%s\">", vprefix,
                              name_escaped);
                IGRAPH_FREE(name_escaped);
                if (ret < 0) {
                    IGRAPH_ERROR("Write failed.", IGRAPH_EFILE);
                }
                IGRAPH_CHECK(igraph_i_attribute_get_string_vertex_attr(graph, name,
                             igraph_vss_1(l), &strv));
                s = igraph_strvector_get(&strv, 0);
                IGRAPH_CHECK(igraph_i_xml_escape(s, &s_escaped));
                ret = fprintf(outstream, "%s", s_escaped);
                IGRAPH_FREE(s_escaped);
                if (ret < 0) {
                    IGRAPH_ERROR("Write failed.", IGRAPH_EFILE);
                }
                ret = fprintf(outstream, "</data>\n");
                if (ret < 0) {
                    IGRAPH_ERROR("Write failed.", IGRAPH_EFILE);
                }
            } else if (VECTOR(vtypes)[i] == IGRAPH_ATTRIBUTE_BOOLEAN) {
                name = igraph_strvector_get(&vnames, i);
                IGRAPH_CHECK(igraph_i_attribute_get_bool_vertex_attr(graph, name,
                             igraph_vss_1(l), &boolv));
                IGRAPH_CHECK(igraph_i_xml_escape(name, &name_escaped));
                ret = fprintf(outstream, "      <data key=\"%s%s\">%s</data>\n",
                              vprefix, name_escaped, VECTOR(boolv)[0] ? "true" : "false");
                IGRAPH_FREE(name_escaped);
                if (ret < 0) {
                    IGRAPH_ERROR("Write failed.", IGRAPH_EFILE);
                }
            }
        }

        ret = fprintf(outstream, "    </node>\n");
        if (ret < 0) {
            IGRAPH_ERROR("Write failed.", IGRAPH_EFILE);
        }
    }

    /* Now the edges */
    IGRAPH_CHECK(igraph_eit_create(graph, igraph_ess_all(IGRAPH_EDGEORDER_ID), &it));
    IGRAPH_FINALLY(igraph_eit_destroy, &it);
    while (!IGRAPH_EIT_END(it)) {
        igraph_integer_t from, to;
        const char *name; char *name_escaped;
        igraph_integer_t edge = IGRAPH_EIT_GET(it);
        igraph_edge(graph, edge, &from, &to);
        ret = fprintf(outstream, "    <edge source=\"n%" IGRAPH_PRId "\" target=\"n%" IGRAPH_PRId "\">\n",
                      from, to);
        if (ret < 0) {
            IGRAPH_ERROR("Write failed.", IGRAPH_EFILE);
        }

        for (i = 0; i < igraph_vector_int_size(&etypes); i++) {
            if (VECTOR(etypes)[i] == IGRAPH_ATTRIBUTE_NUMERIC) {
                name = igraph_strvector_get(&enames, i);
                IGRAPH_CHECK(igraph_i_attribute_get_numeric_edge_attr(graph, name,
                             igraph_ess_1(edge), &numv));
                if (!isnan(VECTOR(numv)[0])) {
                    IGRAPH_CHECK(igraph_i_xml_escape(name, &name_escaped));
                    ret = fprintf(outstream, "      <data key=\"%s%s\">", eprefix, name_escaped);
                    IGRAPH_FREE(name_escaped);
                    if (ret < 0) {
                        IGRAPH_ERROR("Write failed.", IGRAPH_EFILE);
                    }
                    ret = igraph_real_fprintf_precise(outstream, VECTOR(numv)[0]);
                    if (ret < 0) {
                        IGRAPH_ERROR("Write failed.", IGRAPH_EFILE);
                    }
                    ret = fprintf(outstream, "</data>\n");
                    if (ret < 0) {
                        IGRAPH_ERROR("Write failed.", IGRAPH_EFILE);
                    }
                }
            } else if (VECTOR(etypes)[i] == IGRAPH_ATTRIBUTE_STRING) {
                const char *s;
                char *s_escaped;
                name = igraph_strvector_get(&enames, i);
                IGRAPH_CHECK(igraph_i_xml_escape(name, &name_escaped));
                ret = fprintf(outstream, "      <data key=\"%s%s\">", eprefix,
                              name_escaped);
                IGRAPH_FREE(name_escaped);
                if (ret < 0) {
                    IGRAPH_ERROR("Write failed.", IGRAPH_EFILE);
                }
                IGRAPH_CHECK(igraph_i_attribute_get_string_edge_attr(graph, name,
                             igraph_ess_1(edge), &strv));
                s = igraph_strvector_get(&strv, 0);
                IGRAPH_CHECK(igraph_i_xml_escape(s, &s_escaped));
                ret = fprintf(outstream, "%s", s_escaped);
                IGRAPH_FREE(s_escaped);
                if (ret < 0) {
                    IGRAPH_ERROR("Write failed.", IGRAPH_EFILE);
                }
                ret = fprintf(outstream, "</data>\n");
                if (ret < 0) {
                    IGRAPH_ERROR("Write failed.", IGRAPH_EFILE);
                }
            } else if (VECTOR(etypes)[i] == IGRAPH_ATTRIBUTE_BOOLEAN) {
                name = igraph_strvector_get(&enames, i);
                IGRAPH_CHECK(igraph_i_attribute_get_bool_edge_attr(graph, name,
                             igraph_ess_1(edge), &boolv));
                IGRAPH_CHECK(igraph_i_xml_escape(name, &name_escaped));
                ret = fprintf(outstream, "      <data key=\"%s%s\">%s</data>\n",
                              eprefix, name_escaped, VECTOR(boolv)[0] ? "true" : "false");
                IGRAPH_FREE(name_escaped);
                if (ret < 0) {
                    IGRAPH_ERROR("Write failed.", IGRAPH_EFILE);
                }
            }
        }

        ret = fprintf(outstream, "    </edge>\n");
        if (ret < 0) {
            IGRAPH_ERROR("Write failed.", IGRAPH_EFILE);
        }
        IGRAPH_EIT_NEXT(it);
    }
    igraph_eit_destroy(&it);
    IGRAPH_FINALLY_CLEAN(1);

    ret = fprintf(outstream, "  </graph>\n");
    if (ret < 0) {
        IGRAPH_ERROR("Write failed.", IGRAPH_EFILE);
    }
    fprintf(outstream, "</graphml>\n");
    if (ret < 0) {
        IGRAPH_ERROR("Write failed.", IGRAPH_EFILE);
    }

    igraph_strvector_destroy(&gnames);
    igraph_strvector_destroy(&vnames);
    igraph_strvector_destroy(&enames);
    igraph_vector_int_destroy(&gtypes);
    igraph_vector_int_destroy(&vtypes);
    igraph_vector_int_destroy(&etypes);
    igraph_vector_destroy(&numv);
    igraph_strvector_destroy(&strv);
    igraph_vector_bool_destroy(&boolv);
    IGRAPH_FINALLY_CLEAN(9);

    return IGRAPH_SUCCESS;
}
