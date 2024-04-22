/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2005-2022  The igraph development team

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
#include "igraph_version.h"

#include "core/trie.h"
#include "graph/attributes.h"
#include "internal/hacks.h" /* strdup, strncasecmp */
#include "math/safe_intop.h"

#include "io/gml-header.h"
#include "io/parsers/gml-parser.h"

#include <ctype.h>
#include <time.h>
#include <string.h>

int igraph_gml_yylex_init_extra(igraph_i_gml_parsedata_t *user_defined, void *scanner);
int igraph_gml_yylex_destroy(void *scanner);
int igraph_gml_yyparse(igraph_i_gml_parsedata_t *context);
void igraph_gml_yyset_in(FILE *in_str, void *yyscanner);

/* Checks if a null-terminated string needs encoding or decoding.
 *
 * Encoding is needed when an " or & character is present.
 *
 * Decoding is needed when an &xyz; style entity is present, so it's sufficient to look
 * for & characters. " characters are never present in the raw strings returned by the
 * GML parser, so we can use the same function to detect the need for either encoding
 * or decoding.
 */
static igraph_bool_t needs_coding(const char *str) {
    while (*str) {
        if (*str == '&' || *str == '"') {
            return true;
        }
        str++;
    }
    return false;
}

/* Encode & and " character in 'src' to &amp; and &quot;
 * '*dest' must be deallocated by the caller.
 */
static igraph_error_t entity_encode(const char *src, char **dest, igraph_bool_t only_quot) {
    igraph_integer_t destlen = 0;
    const char *s;
    char *d;

    for (s = src; *s != '\0'; s++, destlen++) {
        switch (*s) {
        case '&': /* &amp; */
            if (! only_quot) {
                destlen += 4;
            }
            break;
        case '"': /* &quot; */
            destlen += 5; break;
        }
    }
    *dest = IGRAPH_CALLOC(destlen + 1, char);
    IGRAPH_CHECK_OOM(dest, "Not enough memory to encode string for GML export.");
    for (s = src, d = *dest; *s != '\0'; s++, d++) {
        switch (*s) {
        case '&':
            if (! only_quot) {
                strcpy(d, "&amp;");
                d += 4;
            } else {
                *d = *s;
            }
            break;
        case '"':
            strcpy(d, "&quot;"); d += 5; break;
        default:
            *d = *s;
        }
    }
    *d = '\0';
    return IGRAPH_SUCCESS;
}

/* Decode the five standard predefined XML entities. Unknown entities or stray & characters
 * will be passed through unchanged. '*dest' must be deallocated by the caller.
 * If '*warned' is false, warnings will be issued for unsupported entities and
 * '*warned' will be set to true. This is to prevent a flood of warnings in some files.
 */
static igraph_error_t entity_decode(const char *src, char **dest, igraph_bool_t *warned) {
    const char *entity_names[] = {
        "&quot;", "&amp;", "&apos;", "&lt;", "&gt;"
    };

    const char entity_values[] = {
        '"', '&', '\'', '<', '>'
    };

    const int entity_count = sizeof entity_values / sizeof entity_values[0];

    const char *s;
    char *d;
    size_t len = strlen(src);
    *dest = IGRAPH_CALLOC(len+1, char); /* at most as much storage needed as for 'src' */
    IGRAPH_CHECK_OOM(dest, "Not enough memory to decode string during GML import.");

    for (s = src, d = *dest; *s != '\0';) {
        if (*s == '&') {
            int i;
            for (i=0; i < entity_count; i++) {
                size_t entity_len = strlen(entity_names[i]);
                if (!strncasecmp(s, entity_names[i], entity_len)) {
                    *d++ = entity_values[i];
                    s += entity_len;
                    break;
                }
            }
            /* None of the known entities match, report warning and pass through unchanged. */
            if (i == entity_count) {
                if (! *warned) {
                    const int max_entity_name_length = 34;
                    int j = 0;
                    while (s[j] != '\0' && s[j] != ';' && j < max_entity_name_length) {
                        j++;
                    }
                    if (s[j] == '\0' || j == max_entity_name_length) {
                        IGRAPH_WARNING("Unterminated entity or stray & character found, will be returned verbatim.");
                    } else {
                        IGRAPH_WARNINGF("One or more unknown entities will be returned verbatim (%.*s).", j+1, s);
                    }
                    *warned = true; /* warn only once */
                }
                *d++ = *s++;
            }
        } else {
            *d++ = *s++;
        }
    }
    *d = '\0';

    return IGRAPH_SUCCESS;
}

static void igraph_i_gml_destroy_attrs(igraph_vector_ptr_t **ptr) {

    igraph_vector_ptr_t *vec;
    for (igraph_integer_t i = 0; i < 3; i++) {
        vec = ptr[i];
        for (igraph_integer_t j = 0; j < igraph_vector_ptr_size(vec); j++) {
            igraph_attribute_record_t *atrec = VECTOR(*vec)[j];
            if (atrec->type == IGRAPH_ATTRIBUTE_NUMERIC) {
                igraph_vector_t *value = (igraph_vector_t*)atrec->value;
                if (value != 0) {
                    igraph_vector_destroy(value);
                    IGRAPH_FREE(value);
                }
            } else if (atrec->type == IGRAPH_ATTRIBUTE_STRING) {
                igraph_strvector_t *value = (igraph_strvector_t*)atrec->value;
                if (value != 0) {
                    igraph_strvector_destroy(value);
                    IGRAPH_FREE(value);
                }
            } else {
                /* Some empty attribute records may have been created for composite attributes */
            }
            IGRAPH_FREE(atrec->name);
            IGRAPH_FREE(atrec);
        }
        igraph_vector_ptr_destroy(vec);
    }
}

static igraph_real_t igraph_i_gml_toreal(igraph_gml_tree_t *node, igraph_integer_t pos) {
    igraph_i_gml_tree_type_t type = igraph_gml_tree_type(node, pos);

    switch (type) {
    case IGRAPH_I_GML_TREE_INTEGER:
        return igraph_gml_tree_get_integer(node, pos);
    case IGRAPH_I_GML_TREE_REAL:
        return igraph_gml_tree_get_real(node, pos);
    case IGRAPH_I_GML_TREE_TREE:
        return IGRAPH_NAN; /* default value of NaN when composite is ignored */
    default:
        /* Must never reach here, regardless of the contents of the GML file. */
        IGRAPH_FATALF("Unexpected node type in GML tree, line %" IGRAPH_PRId ".",
                      igraph_gml_tree_line(node, pos)); /* LCOV_EXCL_LINE */
    }
}

static const char *igraph_i_gml_tostring(igraph_gml_tree_t *node, igraph_integer_t pos) {
    igraph_i_gml_tree_type_t type = igraph_gml_tree_type(node, pos);
    static char tmp[100];
    const char *p = tmp;
    igraph_integer_t i;
    igraph_real_t d;

    switch (type) {
    case IGRAPH_I_GML_TREE_INTEGER:
        i = igraph_gml_tree_get_integer(node, pos);
        snprintf(tmp, sizeof(tmp) / sizeof(char), "%" IGRAPH_PRId, i);
        break;
    case IGRAPH_I_GML_TREE_REAL:
        d = igraph_gml_tree_get_real(node, pos);
        igraph_real_snprintf_precise(tmp, sizeof(tmp) / sizeof(char), d);
        break;
    case IGRAPH_I_GML_TREE_STRING:
        p = igraph_gml_tree_get_string(node, pos);
        break;
    case IGRAPH_I_GML_TREE_TREE:
        tmp[0] = '\0'; /* default value of "" when composite is ignored */
        break;
    default:
        /* Must never reach here, regardless of the contents of the GML file. */
        IGRAPH_FATALF("Unexpected node type in GML tree, line %" IGRAPH_PRId ".",
                      igraph_gml_tree_line(node, pos)); /* LCOV_EXCL_LINE */
    }

    return p;
}

igraph_error_t igraph_i_gml_parsedata_init(igraph_i_gml_parsedata_t *context) {
    context->depth = 0;
    context->scanner = NULL;
    context->tree = NULL;
    context->errmsg[0] = '\0';
    context->igraph_errno = IGRAPH_SUCCESS;

    return IGRAPH_SUCCESS;
}

void igraph_i_gml_parsedata_destroy(igraph_i_gml_parsedata_t *context) {
    if (context->tree != NULL) {
        igraph_gml_tree_destroy(context->tree);
        context->tree = NULL;
    }

    if (context->scanner != NULL) {
        (void) igraph_gml_yylex_destroy(context->scanner);
        context->scanner = NULL;
    }
}

/* Takes a vector of attribute records and removes those elements
 * whose type is unspecified, i.e. IGRAPH_ATTRIBUTE_UNSPECIFIED. */
static void prune_unknown_attributes(igraph_vector_ptr_t *attrs) {
    igraph_integer_t i, j;
    for (i = 0, j = 0; i < igraph_vector_ptr_size(attrs); i++) {
        igraph_attribute_record_t *atrec = VECTOR(*attrs)[i];
        if (atrec->type == IGRAPH_ATTRIBUTE_UNSPECIFIED) {
            IGRAPH_FREE(atrec->name);
            IGRAPH_FREE(atrec);
        } else {
            VECTOR(*attrs)[j++] = VECTOR(*attrs)[i];
        }
    }
    igraph_vector_ptr_resize(attrs, j); /* shrinks */
}

/* Converts an integer id to an optionally prefixed string id. */
static const char *strid(igraph_integer_t id, const char *prefix) {
    static char name[100];
    snprintf(name, sizeof(name) / sizeof(char) - 1, "%s%" IGRAPH_PRId, prefix, id);
    return name;
}

/* Creates an empty attribute record or if it exists, updates its type as needed.
 * 'name' is the attribute name. 'type' is the current type in the GML tree,
 * which will determine the igraph attribute type to use. */
static igraph_error_t create_or_update_attribute(const char *name,
                                                 igraph_i_gml_tree_type_t type,
                                                 igraph_trie_t *attrnames,
                                                 igraph_vector_ptr_t *attrs) {

    igraph_integer_t trieid, triesize = igraph_trie_size(attrnames);
    IGRAPH_CHECK(igraph_trie_get(attrnames, name, &trieid));
    if (trieid == triesize) {
        /* new attribute */
        igraph_attribute_record_t *atrec = IGRAPH_CALLOC(1, igraph_attribute_record_t);
        IGRAPH_CHECK_OOM(atrec, "Cannot read GML file.");
        IGRAPH_FINALLY(igraph_free, atrec);

        atrec->name = strdup(name);
        IGRAPH_CHECK_OOM(atrec->name, "Cannot read GML file.");
        IGRAPH_FINALLY(igraph_free, (char *) atrec->name);

        if (type == IGRAPH_I_GML_TREE_INTEGER || type == IGRAPH_I_GML_TREE_REAL) {
            atrec->type = IGRAPH_ATTRIBUTE_NUMERIC;
        } else if (type == IGRAPH_I_GML_TREE_STRING) {
            atrec->type = IGRAPH_ATTRIBUTE_STRING;
        } else {
            atrec->type = IGRAPH_ATTRIBUTE_UNSPECIFIED;
        }
        IGRAPH_CHECK(igraph_vector_ptr_push_back(attrs, atrec));
        IGRAPH_FINALLY_CLEAN(2);
    } else {
        /* already seen, should we update type? */
        igraph_attribute_record_t *atrec = VECTOR(*attrs)[trieid];
        igraph_attribute_type_t type1 = atrec->type;
        if (type == IGRAPH_I_GML_TREE_STRING) {
            atrec->type = IGRAPH_ATTRIBUTE_STRING;
        } else if (type1 == IGRAPH_ATTRIBUTE_UNSPECIFIED) {
            if (type == IGRAPH_I_GML_TREE_INTEGER || type == IGRAPH_I_GML_TREE_REAL) {
                atrec->type = IGRAPH_ATTRIBUTE_NUMERIC;
            }
        }
    }

    return IGRAPH_SUCCESS;
}

/* Allocates the contents of attribute records stored in 'attrs'.
 * 'no_of_items' is the length of attribute vectors, i.e. no_of_nodes,
 * no_of_edges, or 1 for vertex, edge and graph attributes.
 * The 'kind' parameter can be "vertex", "edge" or "graph", and
 * is used solely for showing better warning messages. */
static igraph_error_t allocate_attributes(igraph_vector_ptr_t *attrs,
                                          igraph_integer_t no_of_items,
                                          const char *kind) {

    igraph_integer_t i, n = igraph_vector_ptr_size(attrs);
    for (i = 0; i < n; i++) {
        igraph_attribute_record_t *atrec = VECTOR(*attrs)[i];
        igraph_attribute_type_t type = atrec->type;
        if (type == IGRAPH_ATTRIBUTE_NUMERIC) {
            igraph_vector_t *p = IGRAPH_CALLOC(1, igraph_vector_t);
            IGRAPH_CHECK_OOM(p, "Cannot read GML file.");
            IGRAPH_FINALLY(igraph_free, p);
            IGRAPH_CHECK(igraph_vector_init(p, no_of_items));
            igraph_vector_fill(p, IGRAPH_NAN); /* use NaN as default */
            atrec->value = p;
            IGRAPH_FINALLY_CLEAN(1);
        } else if (type == IGRAPH_ATTRIBUTE_STRING) {
            igraph_strvector_t *p = IGRAPH_CALLOC(1, igraph_strvector_t);
            IGRAPH_CHECK_OOM(p, "Cannot read GML file.");
            IGRAPH_FINALLY(igraph_free, p);
            IGRAPH_CHECK(igraph_strvector_init(p, no_of_items));
            atrec->value = p;
            IGRAPH_FINALLY_CLEAN(1);
        } else if (type == IGRAPH_ATTRIBUTE_UNSPECIFIED) {
            IGRAPH_WARNINGF("Composite %s attribute '%s' ignored in GML file.", kind, atrec->name);
        } else {
            /* Must never reach here. */
            IGRAPH_FATAL("Unexpected attribute type.");
        }
    }
    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_read_graph_gml
 * \brief Read a graph in GML format.
 *
 * GML is a simple textual format, see
 * https://web.archive.org/web/20190207140002/http://www.fim.uni-passau.de/index.php?id=17297%26L=1
 * for details.
 *
 * </para><para>
 * Although all syntactically correct GML can be parsed,
 * we implement only a subset of this format. Some attributes might be
 * ignored. Here is a list of all the differences:
 * \olist
 * \oli Only attributes with a simple type are used: integer, real or
 *      string. If an attribute is composite, i.e. an array or a record,
 *      then it is ignored. When some values of the attribute are simple and
 *      some compound, the composite ones are replaced with a default value
 *      (NaN for numeric, <code>""</code> for string).
 * \oli <code>comment</code> fields are not ignored. They are treated as any
 *      other field and converted to attributes.
 * \oli Top level attributes except for <code>Version</code> and the
 *      first <code>graph</code> attribute are completely ignored.
 * \oli There is no maximum line length or maximum keyword length.
 * \oli Only the \c quot, \c amp, \c apos, \c lt and \c gt character entities
 *      are supported. Any other entity is passed through unchanged by the reader
 *      after issuing a warning, and is expected to be decoded by the user.
 * \oli We allow <code>inf</code>, <code>-inf</code> and <code>nan</code>
 *      (not a number) as a real number. This is case insensitive, so
 *      <code>nan</code>, <code>NaN</code> and <code>NAN</code> are equivalent.
 * \endolist
 *
 * </para><para> Please contact us if you cannot live with these
 * limitations of the GML parser.
 *
 * \param graph Pointer to an uninitialized graph object.
 * \param instream The stream to read the GML file from.
 * \return Error code.
 *
 * Time complexity: should be proportional to the length of the file.
 *
 * \sa \ref igraph_read_graph_graphml() for a more modern format,
 * \ref igraph_write_graph_gml() for writing GML files.
 *
 * \example examples/simple/gml.c
 */
igraph_error_t igraph_read_graph_gml(igraph_t *graph, FILE *instream) {

    igraph_integer_t i;
    igraph_integer_t no_of_nodes = 0, no_of_edges = 0;
    igraph_integer_t node_no;
    igraph_trie_t trie;
    igraph_vector_int_t edges;
    igraph_bool_t directed = IGRAPH_UNDIRECTED;
    igraph_bool_t has_directed = false;
    igraph_gml_tree_t *gtree;
    igraph_integer_t gidx;
    igraph_trie_t vattrnames;
    igraph_trie_t eattrnames;
    igraph_trie_t gattrnames;
    igraph_vector_ptr_t gattrs = IGRAPH_VECTOR_PTR_NULL,
                        vattrs = IGRAPH_VECTOR_PTR_NULL,
                        eattrs = IGRAPH_VECTOR_PTR_NULL;
    igraph_vector_ptr_t *attrs[3];
    igraph_integer_t edgeptr = 0;
    igraph_i_gml_parsedata_t context;
    igraph_bool_t entity_warned = false; /* used to warn at most once about unsupported entities */

    attrs[0] = &gattrs; attrs[1] = &vattrs; attrs[2] = &eattrs;

    IGRAPH_CHECK(igraph_i_gml_parsedata_init(&context));
    IGRAPH_FINALLY(igraph_i_gml_parsedata_destroy, &context);

    igraph_gml_yylex_init_extra(&context, &context.scanner);

    igraph_gml_yyset_in(instream, context.scanner);

    /* Protect 'context' from being destroyed before returning from yyparse() */
    IGRAPH_FINALLY_ENTER();
    int err = igraph_gml_yyparse(&context);
    IGRAPH_FINALLY_EXIT();
    switch (err) {
    case 0: /* success */
        break;
    case 1: /* parse error */
        if (context.errmsg[0] != '\0') {
            IGRAPH_ERROR(context.errmsg, IGRAPH_PARSEERROR);
        } else if (context.igraph_errno != IGRAPH_SUCCESS) {
            IGRAPH_ERROR("", context.igraph_errno);
        } else {
            IGRAPH_ERROR("Cannot read GML file.", IGRAPH_PARSEERROR);
        }
        break;
    case 2: /* out of memory */
        IGRAPH_ERROR("Cannot read GML file.", IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */
        break;
    default: /* must never reach here */
        /* Hint: This will usually be triggered if an IGRAPH_CHECK() is used in a Bison
         * action instead of an IGRAPH_YY_CHECK(), resulting in an igraph errno being
         * returned in place of a Bison error code.
         * TODO: What if future Bison versions introduce error codes other than 0, 1 and 2?
         */
        IGRAPH_FATALF("Parser returned unexpected error code (%d) when reading GML file.", err);  /* LCOV_EXCL_LINE */
    }

    /* Check version, if present, integer and not '1' then ignored */
    i = igraph_gml_tree_find(context.tree, "Version", 0);
    if (i >= 0 &&
        igraph_gml_tree_type(context.tree, i) == IGRAPH_I_GML_TREE_INTEGER &&
        igraph_gml_tree_get_integer(context.tree, i) != 1) {
        IGRAPH_WARNINGF("Unknown GML version: %" IGRAPH_PRId ". "
                        "Parsing will continue assuming GML version 1, but may fail.",
                        igraph_gml_tree_get_integer(context.tree, i));
    }

    /* Get the graph */
    gidx = igraph_gml_tree_find(context.tree, "graph", 0);
    if (gidx == -1) {
        IGRAPH_ERROR("No 'graph' object in GML file.", IGRAPH_PARSEERROR);
    }
    if (igraph_gml_tree_type(context.tree, gidx) !=
        IGRAPH_I_GML_TREE_TREE) {
        IGRAPH_ERRORF("Invalid type for 'graph' object in GML file, line %" IGRAPH_PRId ".", IGRAPH_PARSEERROR,
                      igraph_gml_tree_line(context.tree, gidx));
    }
    gtree = igraph_gml_tree_get_tree(context.tree, gidx);

    IGRAPH_FINALLY(igraph_i_gml_destroy_attrs, attrs);
    IGRAPH_CHECK(igraph_vector_ptr_init(&gattrs, 0));
    IGRAPH_CHECK(igraph_vector_ptr_init(&vattrs, 0));
    IGRAPH_CHECK(igraph_vector_ptr_init(&eattrs, 0));

    IGRAPH_TRIE_INIT_FINALLY(&trie, 0);
    IGRAPH_TRIE_INIT_FINALLY(&vattrnames, 0);
    IGRAPH_TRIE_INIT_FINALLY(&eattrnames, 0);
    IGRAPH_TRIE_INIT_FINALLY(&gattrnames, 0);

    /* Now we go over all objects in the graph to
     *  - collect the attribute names and types
     *  - collect node IDs
     *  - set directedness
     *  - do some checks which the following code relies on
     *
     * The 'id' fields of 'node' objects are converted into strings, so that they
     * can be inserted into a trie and re-encoded as consecutive integers starting
     * at 0. The GML spec allows isolated nodes with no 'id' field. These get a
     * generated string id of the form "n123" consisting of "n" and their count
     * (i.e. ordinal position) within the GML file.
     *
     * We use an attribute type value of IGRAPH_ATTRIBUTE_UNSPECIFIED to mark attribute
     * records which correspond to composite GML values and must therefore be removed
     * before creating the graph.
     */
    node_no = 0;
    for (i = 0; i < igraph_gml_tree_length(gtree); i++) {
        const char *name = igraph_gml_tree_name(gtree, i);
        if (!strcmp(name, "node")) {
            igraph_gml_tree_t *node;
            igraph_bool_t hasid;
            node_no++;
            no_of_nodes++;
            if (igraph_gml_tree_type(gtree, i) != IGRAPH_I_GML_TREE_TREE) {
                IGRAPH_ERRORF("'node' is not a list in GML file, line %" IGRAPH_PRId ".", IGRAPH_PARSEERROR,
                              igraph_gml_tree_line(gtree, i));
            }
            node = igraph_gml_tree_get_tree(gtree, i);
            hasid = false;
            for (igraph_integer_t j = 0; j < igraph_gml_tree_length(node); j++) {
                const char *name = igraph_gml_tree_name(node, j);
                igraph_i_gml_tree_type_t type = igraph_gml_tree_type(node, j);
                IGRAPH_CHECK(create_or_update_attribute(name, type, &vattrnames, &vattrs));
                /* check id */
                if (!strcmp(name, "id")) {
                    igraph_integer_t id, trie_id;
                    igraph_integer_t trie_size = igraph_trie_size(&trie);
                    if (hasid) {
                        /* A 'node' must not have more than one 'id' field.
                         * This error cannot be relaxed into a warning because all ids we find are
                         * added to the trie, and eventually converted to igraph vertex ids. */
                        IGRAPH_ERRORF("Node has multiple 'id' fields in GML file, line %" IGRAPH_PRId ".",
                                      IGRAPH_PARSEERROR,
                                      igraph_gml_tree_line(node, j));
                    }
                    if (type != IGRAPH_I_GML_TREE_INTEGER) {
                        IGRAPH_ERRORF("Non-integer node id in GML file, line %" IGRAPH_PRId ".", IGRAPH_PARSEERROR,
                                      igraph_gml_tree_line(node, j));
                    }
                    id = igraph_gml_tree_get_integer(node, j);
                    IGRAPH_CHECK(igraph_trie_get(&trie, strid(id, ""), &trie_id));
                    if (trie_id != trie_size) {
                        /* This id has already been seen in a previous node. */
                        IGRAPH_ERRORF("Duplicate node id in GML file, line %" IGRAPH_PRId ".", IGRAPH_PARSEERROR,
                                      igraph_gml_tree_line(node, j));
                    }
                    hasid = true;
                }
            }
            if (!hasid) {
                /* Isolated nodes are allowed not to have an id.
                 * We generate an "n"-prefixed string id to be used in the trie. */
                igraph_integer_t trie_id;
                IGRAPH_CHECK(igraph_trie_get(&trie, strid(node_no, "n"), &trie_id));
            }
        } else if (!strcmp(name, "edge")) {
            igraph_gml_tree_t *edge;
            igraph_bool_t has_source = false, has_target = false;
            no_of_edges++;
            if (igraph_gml_tree_type(gtree, i) != IGRAPH_I_GML_TREE_TREE) {
                IGRAPH_ERRORF("'edge' is not a list in GML file, line %" IGRAPH_PRId ".", IGRAPH_PARSEERROR,
                              igraph_gml_tree_line(gtree, i));
            }
            edge = igraph_gml_tree_get_tree(gtree, i);
            for (igraph_integer_t j = 0; j < igraph_gml_tree_length(edge); j++) {
                const char *name = igraph_gml_tree_name(edge, j);
                igraph_i_gml_tree_type_t type = igraph_gml_tree_type(edge, j);
                if (!strcmp(name, "source")) {
                    if (has_source) {
                        /* An edge must not have more than one 'source' field.
                         * This could be relaxed to a warning, but we keep it as an error
                         * for consistency with the handling of duplicate node 'id' field,
                         * and because it indicates a serious corruption in the GML file. */
                        IGRAPH_ERRORF("Duplicate 'source' in an edge in GML file, line %" IGRAPH_PRId ".",
                                      IGRAPH_PARSEERROR,
                                      igraph_gml_tree_line(edge, j));
                    }
                    has_source = true;
                    if (type != IGRAPH_I_GML_TREE_INTEGER) {
                        IGRAPH_ERRORF("Non-integer 'source' for an edge in GML file, line %" IGRAPH_PRId ".",
                                      IGRAPH_PARSEERROR,
                                      igraph_gml_tree_line(edge, j));
                    }
                } else if (!strcmp(name, "target")) {
                    if (has_target) {
                        /* An edge must not have more than one 'target' field. */
                        IGRAPH_ERRORF("Duplicate 'target' in an edge in GML file, line %" IGRAPH_PRId ".",
                                      IGRAPH_PARSEERROR,
                                      igraph_gml_tree_line(edge, j));
                    }
                    has_target = true;
                    if (type != IGRAPH_I_GML_TREE_INTEGER) {
                        IGRAPH_ERRORF("Non-integer 'target' for an edge in GML file, line %" IGRAPH_PRId ".",
                                      IGRAPH_PARSEERROR,
                                      igraph_gml_tree_line(edge, j));
                    }
                } else {
                    IGRAPH_CHECK(create_or_update_attribute(name, type, &eattrnames, &eattrs));
                }
            } /* for */
            if (!has_source) {
                IGRAPH_ERRORF("No 'source' for edge in GML file, line %" IGRAPH_PRId ".", IGRAPH_PARSEERROR,
                              igraph_gml_tree_line(gtree, i));
            }
            if (!has_target) {
                IGRAPH_ERRORF("No 'target' for edge in GML file, line %" IGRAPH_PRId ".", IGRAPH_PARSEERROR,
                              igraph_gml_tree_line(gtree, i));
            }
        } else if (! strcmp(name, "directed")) {
            /* Set directedness of graph. */
            if (has_directed) {
                /* Be tolerant of duplicate entries, but do show a warning. */
                IGRAPH_WARNINGF("Duplicate 'directed' field in 'graph', line %" IGRAPH_PRId ". "
                                "Ignoring previous 'directed' fields.",
                                igraph_gml_tree_line(gtree, i));
            }
            if (igraph_gml_tree_type(gtree, i) == IGRAPH_I_GML_TREE_INTEGER) {
                igraph_integer_t dir = igraph_gml_tree_get_integer(gtree, i);
                if (dir != 0 && dir != 1) {
                    IGRAPH_WARNINGF(
                        "Invalid value %" IGRAPH_PRId " for 'directed' attribute on line %" IGRAPH_PRId ", should be 0 or 1.",
                        dir, igraph_gml_tree_line(gtree, i));
                }
                if (dir) {
                    directed = IGRAPH_DIRECTED;
                }
                has_directed = true;
            } else {
                IGRAPH_WARNINGF("Invalid type for 'directed' attribute on line %" IGRAPH_PRId ", assuming undirected.",
                                igraph_gml_tree_line(gtree, i));
            }
        } else {
            /* Add the rest of items as graph attributes. */
            igraph_i_gml_tree_type_t type = igraph_gml_tree_type(gtree, i);
            IGRAPH_CHECK(create_or_update_attribute(name, type, &gattrnames, &gattrs));
        }
    }

    /* At this point, all nodes must have an id (from the file or generated) stored
     * in the trie. Any condition that violates this should have been caught during
     * the preceding checks. */
    IGRAPH_ASSERT(igraph_trie_size(&trie) == no_of_nodes);

    /* Now we allocate the vectors and strvectors for the attributes */
    IGRAPH_CHECK(allocate_attributes(&vattrs, no_of_nodes, "vertex"));
    IGRAPH_CHECK(allocate_attributes(&eattrs, no_of_edges, "edge"));
    IGRAPH_CHECK(allocate_attributes(&gattrs, 1, "graph"));

    /* Add edges, edge attributes and vertex attributes */
    IGRAPH_VECTOR_INT_INIT_FINALLY(&edges, no_of_edges * 2);
    node_no = 0;
    for (i = 0; i < igraph_gml_tree_length(gtree); i++) {
        const char *name;
        name = igraph_gml_tree_name(gtree, i);
        if (!strcmp(name, "node")) {
            igraph_gml_tree_t *node = igraph_gml_tree_get_tree(gtree, i);
            igraph_integer_t iidx = igraph_gml_tree_find(node, "id", 0);
            igraph_integer_t trie_id;
            const char *sid;
            node_no++;
            if (iidx < 0) {
                /* Isolated node with no id field, use n-prefixed generated id */
                sid = strid(node_no, "n");
            } else {
                sid = strid(igraph_gml_tree_get_integer(node, iidx), "");
            }
            IGRAPH_CHECK(igraph_trie_get(&trie, sid, &trie_id));
            for (igraph_integer_t j = 0; j < igraph_gml_tree_length(node); j++) {
                const char *aname = igraph_gml_tree_name(node, j);
                igraph_attribute_record_t *atrec;
                igraph_attribute_type_t type;
                igraph_integer_t ai;
                IGRAPH_CHECK(igraph_trie_get(&vattrnames, aname, &ai));
                atrec = VECTOR(vattrs)[ai];
                type = atrec->type;
                if (type == IGRAPH_ATTRIBUTE_NUMERIC) {
                    igraph_vector_t *v = (igraph_vector_t *) atrec->value;
                    VECTOR(*v)[trie_id] = igraph_i_gml_toreal(node, j);
                } else if (type == IGRAPH_ATTRIBUTE_STRING) {
                    igraph_strvector_t *v = (igraph_strvector_t *) atrec->value;
                    const char *value = igraph_i_gml_tostring(node, j);
                    if (needs_coding(value)) {
                        char *value_decoded;
                        IGRAPH_CHECK(entity_decode(value, &value_decoded, &entity_warned));
                        IGRAPH_FINALLY(igraph_free, value_decoded);
                        IGRAPH_CHECK(igraph_strvector_set(v, trie_id, value_decoded));
                        IGRAPH_FREE(value_decoded);
                        IGRAPH_FINALLY_CLEAN(1);
                    } else {
                        IGRAPH_CHECK(igraph_strvector_set(v, trie_id, value));
                    }
                } else {
                    /* Ignored composite attribute */
                }
            }
        } else if (!strcmp(name, "edge")) {
            igraph_gml_tree_t *edge;
            igraph_integer_t from, to, fromidx = 0, toidx = 0;
            edge = igraph_gml_tree_get_tree(gtree, i);
            for (igraph_integer_t j = 0; j < igraph_gml_tree_length(edge); j++) {
                const char *aname = igraph_gml_tree_name(edge, j);
                if (!strcmp(aname, "source")) {
                    fromidx = igraph_gml_tree_find(edge, "source", 0);
                } else if (!strcmp(aname, "target")) {
                    toidx = igraph_gml_tree_find(edge, "target", 0);
                } else {
                    igraph_integer_t edgeid = edgeptr / 2;
                    igraph_integer_t ai;
                    igraph_attribute_record_t *atrec;
                    igraph_attribute_type_t type;
                    IGRAPH_CHECK(igraph_trie_get(&eattrnames, aname, &ai));
                    atrec = VECTOR(eattrs)[ai];
                    type = atrec->type;
                    if (type == IGRAPH_ATTRIBUTE_NUMERIC) {
                        igraph_vector_t *v = (igraph_vector_t *) atrec->value;
                        VECTOR(*v)[edgeid] = igraph_i_gml_toreal(edge, j);
                    } else if (type == IGRAPH_ATTRIBUTE_STRING) {
                        igraph_strvector_t *v = (igraph_strvector_t *) atrec->value;
                        const char *value = igraph_i_gml_tostring(edge, j);
                        if (needs_coding(value)) {
                            char *value_decoded;
                            IGRAPH_CHECK(entity_decode(value, &value_decoded, &entity_warned));
                            IGRAPH_FINALLY(igraph_free, value_decoded);
                            IGRAPH_CHECK(igraph_strvector_set(v, edgeid, value_decoded));
                            IGRAPH_FREE(value_decoded);
                            IGRAPH_FINALLY_CLEAN(1);
                        } else {
                            IGRAPH_CHECK(igraph_strvector_set(v, edgeid, value));
                        }
                    } else {
                        /* Ignored composite attribute */
                    }
                }
            }
            from = igraph_gml_tree_get_integer(edge, fromidx);
            to = igraph_gml_tree_get_integer(edge, toidx);
            IGRAPH_CHECK(igraph_trie_check(&trie, strid(from, ""), &from));
            if (from < 0) {
                IGRAPH_ERRORF("Unknown source node id found in an edge in GML file, line %" IGRAPH_PRId ".",
                             IGRAPH_PARSEERROR, igraph_gml_tree_line(edge, fromidx));
            }
            IGRAPH_CHECK(igraph_trie_check(&trie, strid(to, ""), &to));
            if (to < 0) {
                IGRAPH_ERRORF("Unknown target node id found in an edge in GML file, line %" IGRAPH_PRId ".",
                             IGRAPH_PARSEERROR, igraph_gml_tree_line(edge, toidx));
            }
            VECTOR(edges)[edgeptr++] = from;
            VECTOR(edges)[edgeptr++] = to;
        } else if (! strcmp(name, "directed")) {
            /* Nothing to do for 'directed' field, already handled earlier. */
        } else {
            /* Set the rest as graph attributes */
            igraph_integer_t ai;
            igraph_attribute_record_t *atrec;
            igraph_attribute_type_t type;
            IGRAPH_CHECK(igraph_trie_get(&gattrnames, name, &ai));
            atrec = VECTOR(gattrs)[ai];
            type = atrec->type;
            if (type == IGRAPH_ATTRIBUTE_NUMERIC) {
                igraph_vector_t *v = (igraph_vector_t *) atrec->value;
                VECTOR(*v)[0] = igraph_i_gml_toreal(gtree, i);
            } else if (type == IGRAPH_ATTRIBUTE_STRING) {
                igraph_strvector_t *v = (igraph_strvector_t *) atrec->value;
                const char *value = igraph_i_gml_tostring(gtree, i);
                if (needs_coding(value)) {
                    char *value_decoded;
                    IGRAPH_CHECK(entity_decode(value, &value_decoded, &entity_warned));
                    IGRAPH_FINALLY(igraph_free, value_decoded);
                    IGRAPH_CHECK(igraph_strvector_set(v, 0, value_decoded));
                    IGRAPH_FREE(value_decoded);
                    IGRAPH_FINALLY_CLEAN(1);
                } else {
                    IGRAPH_CHECK(igraph_strvector_set(v, 0, value));
                }
            } else {
                /* Ignored composite attribute */
            }
        }
    }

    /* Remove composite attributes */
    prune_unknown_attributes(&vattrs);
    prune_unknown_attributes(&eattrs);
    prune_unknown_attributes(&gattrs);

    igraph_trie_destroy(&trie);
    igraph_trie_destroy(&gattrnames);
    igraph_trie_destroy(&vattrnames);
    igraph_trie_destroy(&eattrnames);
    IGRAPH_FINALLY_CLEAN(4);

    IGRAPH_CHECK(igraph_empty_attrs(graph, 0, directed, &gattrs));
    IGRAPH_FINALLY(igraph_destroy, graph);
    IGRAPH_CHECK(igraph_add_vertices(graph, no_of_nodes, &vattrs));
    IGRAPH_CHECK(igraph_add_edges(graph, &edges, &eattrs));
    IGRAPH_FINALLY_CLEAN(1); /* do not destroy 'graph', just pop it from the stack */

    igraph_vector_int_destroy(&edges);
    igraph_i_gml_destroy_attrs(attrs);
    igraph_i_gml_parsedata_destroy(&context);
    IGRAPH_FINALLY_CLEAN(3);

    return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_i_gml_convert_to_key(const char *orig, char **key) {
    char strno[50];
    size_t i, len = strlen(orig), newlen = 0, plen = 0;

    /* do we need a prefix? */
    if (len == 0 || !isalpha(orig[0])) {
        snprintf(strno, sizeof(strno) - 1, "igraph");
        plen = newlen = strlen(strno);
    }
    for (i = 0; i < len; i++) {
        if (isalnum(orig[i])) {
            newlen++;
        }
    }
    *key = IGRAPH_CALLOC(newlen + 1, char);
    IGRAPH_CHECK_OOM(*key, "Writing GML format failed.");
    memcpy(*key, strno, plen * sizeof(char));
    for (i = 0; i < len; i++) {
        if (isalnum(orig[i])) {
            (*key)[plen++] = orig[i];
        }
    }
    (*key)[newlen] = '\0';

    return IGRAPH_SUCCESS;
}

/* Checks if a vector is free of duplicates. Since NaN == NaN is false, duplicate NaN values
 * will not be detected. */
static igraph_error_t igraph_i_vector_is_duplicate_free(const igraph_vector_t *v, igraph_bool_t *res) {
    igraph_vector_t u;
    igraph_integer_t n = igraph_vector_size(v);

    if (n < 2) {
        *res = true;
        return IGRAPH_SUCCESS;
    }

    IGRAPH_CHECK(igraph_vector_init_copy(&u, v));
    IGRAPH_FINALLY(igraph_vector_destroy, &u);
    igraph_vector_sort(&u);

    *res = true;
    for (igraph_integer_t i=1; i < n; i++) {
        if (VECTOR(u)[i-1] == VECTOR(u)[i]) {
            *res = false;
            break;
        }
    }

    igraph_vector_destroy(&u);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}

#define CHECK(cmd) do { int ret=cmd; if (ret<0) IGRAPH_ERROR("Writing GML format failed.", IGRAPH_EFILE); } while (0)

/**
 * \function igraph_write_graph_gml
 * \brief Write the graph to a stream in GML format.
 *
 * GML is a quite general textual format, see
 * https://web.archive.org/web/20190207140002/http://www.fim.uni-passau.de/index.php?id=17297%26L=1
 * for details.
 *
 * </para><para>
 * The graph, vertex and edges attributes are written to the
 * file as well, if they are numeric or string. Boolean attributes are converted
 * to numeric, with 0 and 1 used for false and true, respectively.
 * NaN values of numeric attributes are skipped, as NaN is not part of the GML
 * specification and other software may not be able to read files containing them.
 * This is consistent with \ref igraph_read_graph_gml(), which produces NaN
 * when an attribute value is missing. In contrast with NaN, infinite values
 * are retained. Ensure that none of the numeric attributes values are infinite
 * to produce a conformant GML file that can be read by other software.
 *
 * </para><para>
 * As igraph is more forgiving about attribute names, it might
 * be necessary to simplify the them before writing to the GML file.
 * This way we'll have a syntactically correct GML file. The following
 * simple procedure is performed on each attribute name: first the alphanumeric
 * characters are extracted, the others are ignored. Then if the first character
 * is not a letter then the attribute name is prefixed with <quote>igraph</quote>.
 * Note that this might result identical names for two attributes, igraph
 * does not check this.
 *
 * </para><para>
 * The <quote>id</quote> vertex attribute is treated specially.
 * If the <parameter>id</parameter> argument is not \c NULL then it should be a numeric
 * vector with the vertex IDs and the <quote>id</quote> vertex attribute is
 * ignored (if there is one). If <parameter>id</parameter> is \c NULL and there is a
 * numeric <quote>id</quote> vertex attribute, it will be used instead. If ids
 * are not specified in either way then the regular igraph vertex IDs are used.
 * If some of the supplied id values are invalid (non-integer or NaN), all supplied
 * id are ignored and igraph vertex IDs are used instead.
 *
 * </para><para>
 * Note that whichever way vertex IDs are specified, their uniqueness is not checked.
 *
 * </para><para>
 * If the graph has edge attributes that become <quote>source</quote>
 * or <quote>target</quote> after encoding, or the graph has an attribute that becomes
 * <quote>directed</quote>, they will be ignored with a warning. GML uses these attributes
 * to specify the edge endpoints, and the graph directedness, so we cannot write them
 * to the file. Rename them before calling this function if you want to preserve them.
 *
 * \param graph The graph to write to the stream.
 * \param outstream The stream to write the file to.
 * \param options Set of <code>|</code>-combinable boolean flags for writing the GML file.
 *        \clist
 *          \cli 0
 *               All options turned off.
 *          \cli IGRAPH_WRITE_GML_DEFAULT_SW
 *               Default options, currently equivalent to 0. May change in future versions.
 *          \cli IGRAPH_WRITE_GML_ENCODE_ONLY_QUOT_SW
 *               Do not encode any other characters than " as entities. Specifically, this
 *               option prevents the encoding of &amp;. Useful when re-exporting a graph
 *               that was read from a GML file in which igraph could not interpret all entities,
 *               and thus passed them through without decoding.
 *        \endclist
 * \param id Either <code>NULL</code> or a numeric vector with the vertex IDs.
 *        See details above.
 * \param creator An optional string to write to the stream in the creator line.
 *        If \c NULL, the igraph version with the current date and time is added.
 *        If <code>""</code>, the creator line is omitted. Otherwise, the
 *        supplied string is used verbatim.
 * \return Error code.
 *
 * Time complexity: should be proportional to the number of characters written
 * to the file.
 *
 * \sa \ref igraph_read_graph_gml() for reading GML files,
 * \ref igraph_read_graph_graphml() for a more modern format.
 *
 * \example examples/simple/gml.c
 */

igraph_error_t igraph_write_graph_gml(const igraph_t *graph, FILE *outstream,
                                      igraph_write_gml_sw_t options,
                                      const igraph_vector_t *id, const char *creator) {
    igraph_strvector_t gnames, vnames, enames; /* attribute names */
    igraph_vector_int_t gtypes, vtypes, etypes; /* attribute types */
    igraph_integer_t gattr_no, vattr_no, eattr_no; /* attribute counts */
    igraph_vector_t numv;
    igraph_strvector_t strv;
    igraph_vector_bool_t boolv;
    igraph_integer_t i;
    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_integer_t no_of_edges = igraph_ecount(graph);

    /* Each element is a bit field used to prevent showing more
     * than one warning for each vertex or edge attribute. */
    igraph_vector_int_t warning_shown;

    igraph_vector_t v_myid;
    const igraph_vector_t *myid = id;

    /* Creator line */
    if (creator == NULL) {
        time_t curtime = time(0);
        char *timestr = ctime(&curtime);
        timestr[strlen(timestr) - 1] = '\0'; /* nicely remove \n */

        CHECK(fprintf(outstream,
                      "Creator \"igraph version %s %s\"\n",
                      IGRAPH_VERSION, timestr));
    } else if (creator[0] == '\0') {
        /* creator == "", omit Creator line */
    } else {
        if (needs_coding(creator)) {
            char *d;
            IGRAPH_CHECK(entity_encode(creator, &d, IGRAPH_WRITE_GML_ENCODE_ONLY_QUOT_SW & options));
            IGRAPH_FINALLY(igraph_free, d);
            CHECK(fprintf(outstream,
                          "Creator \"%s\"\n",
                          creator));
            IGRAPH_FREE(d);
            IGRAPH_FINALLY_CLEAN(1);
        } else {
            CHECK(fprintf(outstream,
                          "Creator \"%s\"\n",
                          creator));
        }
    }

    /* Version line */
    CHECK(fprintf(outstream, "Version 1\n"));

    /* The graph */
    CHECK(fprintf(outstream, "graph\n[\n"));

    IGRAPH_STRVECTOR_INIT_FINALLY(&gnames, 0);
    IGRAPH_STRVECTOR_INIT_FINALLY(&vnames, 0);
    IGRAPH_STRVECTOR_INIT_FINALLY(&enames, 0);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&gtypes, 0);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&vtypes, 0);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&etypes, 0);
    IGRAPH_CHECK(igraph_i_attribute_get_info(graph,
                 &gnames, &gtypes,
                 &vnames, &vtypes,
                 &enames, &etypes));
    gattr_no = igraph_vector_int_size(&gtypes);
    vattr_no = igraph_vector_int_size(&vtypes);
    eattr_no = igraph_vector_int_size(&etypes);

    IGRAPH_VECTOR_INIT_FINALLY(&numv, 1);
    IGRAPH_STRVECTOR_INIT_FINALLY(&strv, 1);
    IGRAPH_VECTOR_BOOL_INIT_FINALLY(&boolv, 1);

    /* Check whether there is an 'id' node attribute if the supplied is 0 */
    if (!id) {
        igraph_bool_t found = false;
        for (i = 0; i < igraph_vector_int_size(&vtypes); i++) {
            const char *n = igraph_strvector_get(&vnames, i);
            if (!strcmp(n, "id") && VECTOR(vtypes)[i] == IGRAPH_ATTRIBUTE_NUMERIC) {
                found = true; break;
            }
        }
        if (found) {
            IGRAPH_VECTOR_INIT_FINALLY(&v_myid, no_of_nodes);
            IGRAPH_CHECK(igraph_i_attribute_get_numeric_vertex_attr(graph, "id",
                         igraph_vss_all(),
                         &v_myid));
            myid = &v_myid;
        }
    }

    /* Scan id vector for invalid values. If any are found, all ids are ignored.
     * Invalid values may occur as a result of reading a GML file in which some
     * nodes did not have an id, or by adding new vertices to a graph with an "id"
     * attribute. In this case, the "id" attribute will contain NaN values.
     */
    if (myid) {
        if (igraph_vector_size(myid) != no_of_nodes) {
            IGRAPH_ERROR("Size of id vector must match vertex count.", IGRAPH_EINVAL);
        }
        for (i = 0; i < no_of_nodes; ++i) {
            igraph_real_t val = VECTOR(*myid)[i];
            igraph_real_t trunc_val = trunc(val);
            if (! (val == trunc_val && igraph_i_is_real_representable_as_integer(trunc_val))) {
                IGRAPH_WARNINGF("%g is not a valid integer id for GML files, ignoring all supplied ids.", val);
                if (myid == &v_myid) {
                    igraph_vector_destroy(&v_myid);
                    IGRAPH_FINALLY_CLEAN(1);
                }
                myid = NULL;
                break;
            }
        }
    }

    if (myid) {
        igraph_bool_t duplicate_free;
        IGRAPH_CHECK(igraph_i_vector_is_duplicate_free(myid, &duplicate_free));
        if (! duplicate_free) {
            IGRAPH_WARNING("Duplicate id values found, ignoring supplies ids.");
            if (myid == &v_myid) {
                igraph_vector_destroy(&v_myid);
                IGRAPH_FINALLY_CLEAN(1);
            }
            myid = NULL;
        }
    }

    /* directedness */
    CHECK(fprintf(outstream, "  directed %i\n", igraph_is_directed(graph) ? 1 : 0));

    /* Graph attributes first */
    for (i = 0; i < gattr_no; i++) {
        const char *name;
        char *newname;
        name = igraph_strvector_get(&gnames, i);
        IGRAPH_CHECK(igraph_i_gml_convert_to_key(name, &newname));
        IGRAPH_FINALLY(igraph_free, newname);
        if (!strcmp(newname, "directed")|| !strcmp(newname, "edge") || !strcmp(newname, "node")) {
            IGRAPH_WARNINGF("The graph attribute '%s' was ignored while writing GML format.", name);
        } else {
            if (VECTOR(gtypes)[i] == IGRAPH_ATTRIBUTE_NUMERIC) {
                IGRAPH_CHECK(igraph_i_attribute_get_numeric_graph_attr(graph, name, &numv));
                /* Treat NaN as missing, skip writing it. GML does not officially support NaN. */
                if (! isnan(VECTOR(numv)[0])) {
                    if (! isfinite(VECTOR(numv)[0])) {
                        IGRAPH_WARNINGF("Infinite value in numeric graph attribute '%s'. "
                                        "Produced GML file will not be conformant.", name);
                    }
                    CHECK(fprintf(outstream, "  %s ", newname));
                    CHECK(igraph_real_fprintf_precise(outstream, VECTOR(numv)[0]));
                    CHECK(fputc('\n', outstream));
                }
            } else if (VECTOR(gtypes)[i] == IGRAPH_ATTRIBUTE_STRING) {
                const char *s;
                IGRAPH_CHECK(igraph_i_attribute_get_string_graph_attr(graph, name, &strv));
                s = igraph_strvector_get(&strv, 0);
                if (needs_coding(s)) {
                    char *d;
                    IGRAPH_CHECK(entity_encode(s, &d, IGRAPH_WRITE_GML_ENCODE_ONLY_QUOT_SW & options));
                    IGRAPH_FINALLY(igraph_free, d);
                    CHECK(fprintf(outstream, "  %s \"%s\"\n", newname, d));
                    IGRAPH_FREE(d);
                    IGRAPH_FINALLY_CLEAN(1);
                } else {
                    CHECK(fprintf(outstream, "  %s \"%s\"\n", newname, s));
                }
            } else if (VECTOR(gtypes)[i] == IGRAPH_ATTRIBUTE_BOOLEAN) {
                IGRAPH_CHECK(igraph_i_attribute_get_bool_graph_attr(graph, name, &boolv));
                CHECK(fprintf(outstream, "  %s %d\n", newname, VECTOR(boolv)[0] ? 1 : 0));
                IGRAPH_WARNING("A boolean graph attribute was converted to numeric.");
            } else {
                IGRAPH_WARNING("A non-numeric, non-string, non-boolean graph attribute ignored.");
            }
        }
        IGRAPH_FREE(newname);
        IGRAPH_FINALLY_CLEAN(1);
    }

    /* Macros used to work with the bit fiels in 'warning_shown',
     * and avoid showing warnings more than once for each attribute. */
#define GETBIT(k, i) ((k) & (1 << i))
#define SETBIT(k, i) ((k) |= (1 << i))
#define WARN_ONCE(attrno, bit, warn) \
    do { \
        igraph_integer_t *p = &VECTOR(warning_shown)[attrno]; \
        if (! GETBIT(*p, bit)) { \
            warn; \
            SETBIT(*p, bit); \
        } \
    } while (0)


    /* Now come the vertices */
    IGRAPH_VECTOR_INT_INIT_FINALLY(&warning_shown, vattr_no);
    for (i = 0; i < no_of_nodes; i++) {
        igraph_integer_t j;
        CHECK(fprintf(outstream, "  node\n  [\n"));
        /* id */
        CHECK(fprintf(outstream, "    id %" IGRAPH_PRId "\n", myid ? (igraph_integer_t)VECTOR(*myid)[i] : i));
        /* other attributes */
        for (j = 0; j < vattr_no; j++) {
            igraph_attribute_type_t type = (igraph_attribute_type_t) VECTOR(vtypes)[j];
            const char *name;
            char *newname;
            name = igraph_strvector_get(&vnames, j);
            if (!strcmp(name, "id")) {
                /* No warning, the presence of this attribute is expected, and is handled specially. */
                continue;
            }
            IGRAPH_CHECK(igraph_i_gml_convert_to_key(name, &newname));
            IGRAPH_FINALLY(igraph_free, newname);
            if (!strcmp(newname, "id")) {
                /* In case an attribute name would conflict with 'id' only after encoding. */
                WARN_ONCE(j, 0,
                          IGRAPH_WARNINGF("The vertex attribute '%s' was ignored while writing GML format.", name));
            } else {
                if (type == IGRAPH_ATTRIBUTE_NUMERIC) {
                    IGRAPH_CHECK(igraph_i_attribute_get_numeric_vertex_attr(graph, name,
                                 igraph_vss_1(i), &numv));
                    /* Treat NaN as missing, skip writing it. GML does not officially support NaN. */
                    if (! isnan(VECTOR(numv)[0])) {
                        if (! isfinite(VECTOR(numv)[0])) {
                            WARN_ONCE(j, 3,
                                      IGRAPH_WARNINGF("Infinite value in numeric vertex attribute '%s'. "
                                                      "Produced GML file will not be conformant.", name));
                        }
                        CHECK(fprintf(outstream, "    %s ", newname));
                        CHECK(igraph_real_fprintf_precise(outstream, VECTOR(numv)[0]));
                        CHECK(fputc('\n', outstream));
                    }
                } else if (type == IGRAPH_ATTRIBUTE_STRING) {
                    const char *s;
                    IGRAPH_CHECK(igraph_i_attribute_get_string_vertex_attr(graph, name,
                                 igraph_vss_1(i), &strv));
                    s = igraph_strvector_get(&strv, 0);
                    if (needs_coding(s)) {
                        char *d;
                        IGRAPH_CHECK(entity_encode(s, &d, IGRAPH_WRITE_GML_ENCODE_ONLY_QUOT_SW & options));
                        IGRAPH_FINALLY(igraph_free, d);
                        CHECK(fprintf(outstream, "    %s \"%s\"\n", newname, d));
                        IGRAPH_FREE(d);
                        IGRAPH_FINALLY_CLEAN(1);
                    } else {
                        CHECK(fprintf(outstream, "    %s \"%s\"\n", newname, s));
                    }
                } else if (type == IGRAPH_ATTRIBUTE_BOOLEAN) {
                    IGRAPH_CHECK(igraph_i_attribute_get_bool_vertex_attr(graph, name,
                                 igraph_vss_1(i), &boolv));
                    CHECK(fprintf(outstream, "    %s %d\n", newname, VECTOR(boolv)[0] ? 1 : 0));
                    WARN_ONCE(j, 1,
                              IGRAPH_WARNINGF("The boolean vertex attribute '%s' was converted to numeric.", name));
                } else {
                    WARN_ONCE(j, 2,
                              IGRAPH_WARNINGF("The non-numeric, non-string, non-boolean vertex attribute '%s' was ignored.", name));
                }
            }
            IGRAPH_FREE(newname);
            IGRAPH_FINALLY_CLEAN(1);
        }
        CHECK(fprintf(outstream, "  ]\n"));
    }

    /* The edges too */
    IGRAPH_CHECK(igraph_vector_int_resize(&warning_shown, eattr_no));
    igraph_vector_int_fill(&warning_shown, 0);
    for (i = 0; i < no_of_edges; i++) {
        igraph_integer_t from = IGRAPH_FROM(graph, i);
        igraph_integer_t to = IGRAPH_TO(graph, i);
        igraph_integer_t j;
        CHECK(fprintf(outstream, "  edge\n  [\n"));
        /* source and target */
        CHECK(fprintf(outstream, "    source %" IGRAPH_PRId "\n",
                      myid ? (igraph_integer_t)VECTOR(*myid)[from] : from));
        CHECK(fprintf(outstream, "    target %" IGRAPH_PRId "\n",
                      myid ? (igraph_integer_t)VECTOR(*myid)[to] : to));

        /* other attributes */
        for (j = 0; j < eattr_no; j++) {
            igraph_attribute_type_t type = (igraph_attribute_type_t) VECTOR(etypes)[j];
            const char *name;
            char *newname;
            name = igraph_strvector_get(&enames, j);
            IGRAPH_CHECK(igraph_i_gml_convert_to_key(name, &newname));
            IGRAPH_FINALLY(igraph_free, newname);
            if (!strcmp(newname, "source") || !strcmp(newname, "target")) {
                WARN_ONCE(j, 0,
                          IGRAPH_WARNINGF("The edge attribute '%s' was ignored while writing GML format.", name));
            } else {
                if (type == IGRAPH_ATTRIBUTE_NUMERIC) {
                    IGRAPH_CHECK(igraph_i_attribute_get_numeric_edge_attr(graph, name,
                                 igraph_ess_1(i), &numv));
                    /* Treat NaN as missing, skip writing it. GML does not officially support NaN. */
                    if (! isnan(VECTOR(numv)[0])) {
                        if (! isfinite(VECTOR(numv)[0])) {
                            WARN_ONCE(j, 3,
                                      IGRAPH_WARNINGF("Infinite value in numeric edge attribute '%s'. "
                                                      "Produced GML file will not be conformant.", name));
                        }
                        CHECK(fprintf(outstream, "    %s ", newname));
                        CHECK(igraph_real_fprintf_precise(outstream, VECTOR(numv)[0]));
                        CHECK(fputc('\n', outstream));
                    }
                } else if (type == IGRAPH_ATTRIBUTE_STRING) {
                    const char *s;
                    IGRAPH_CHECK(igraph_i_attribute_get_string_edge_attr(graph, name,
                                 igraph_ess_1(i), &strv));
                    s = igraph_strvector_get(&strv, 0);
                    if (needs_coding(s)) {
                        char *d;
                        IGRAPH_CHECK(entity_encode(s, &d, IGRAPH_WRITE_GML_ENCODE_ONLY_QUOT_SW & options));
                        IGRAPH_FINALLY(igraph_free, d);
                        CHECK(fprintf(outstream, "    %s \"%s\"\n", newname, d));
                        IGRAPH_FREE(d);
                        IGRAPH_FINALLY_CLEAN(1);
                    } else {
                        CHECK(fprintf(outstream, "    %s \"%s\"\n", newname, s));
                    }
                } else if (type == IGRAPH_ATTRIBUTE_BOOLEAN) {
                    IGRAPH_CHECK(igraph_i_attribute_get_bool_edge_attr(graph, name,
                                 igraph_ess_1(i), &boolv));
                    CHECK(fprintf(outstream, "    %s %d\n", newname, VECTOR(boolv)[0] ? 1 : 0));
                    WARN_ONCE(j, 1,
                              IGRAPH_WARNINGF("The boolean edge attribute '%s' was converted to numeric.", name));
                } else {
                    WARN_ONCE(j, 2,
                              IGRAPH_WARNINGF("The non-numeric, non-string, non-boolean edge attribute '%s' was ignored.", name));
                }
            }
            IGRAPH_FREE(newname);
            IGRAPH_FINALLY_CLEAN(1);
        }
        CHECK(fprintf(outstream, "  ]\n"));
    }

    CHECK(fprintf(outstream, "]\n"));

#undef GETBIT
#undef SETBIT
#undef WARN_ONCE

    if (&v_myid == myid) {
        igraph_vector_destroy(&v_myid);
        IGRAPH_FINALLY_CLEAN(1);
    }

    igraph_vector_int_destroy(&warning_shown);
    igraph_vector_bool_destroy(&boolv);
    igraph_strvector_destroy(&strv);
    igraph_vector_destroy(&numv);
    igraph_vector_int_destroy(&etypes);
    igraph_vector_int_destroy(&vtypes);
    igraph_vector_int_destroy(&gtypes);
    igraph_strvector_destroy(&enames);
    igraph_strvector_destroy(&vnames);
    igraph_strvector_destroy(&gnames);
    IGRAPH_FINALLY_CLEAN(10);

    return IGRAPH_SUCCESS;
}

#undef CHECK
