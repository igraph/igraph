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
#include "internal/hacks.h" /* strdup */
#include "io/gml-header.h"

#include <ctype.h>
#include <time.h>
#include <string.h>

int igraph_gml_yylex_init_extra (igraph_i_gml_parsedata_t* user_defined,
                                 void* scanner);
void igraph_gml_yylex_destroy (void *scanner );
int igraph_gml_yyparse (igraph_i_gml_parsedata_t* context);
void igraph_gml_yyset_in  (FILE * in_str, void* yyscanner );

static void igraph_i_gml_destroy_attrs(igraph_vector_ptr_t **ptr) {
    igraph_integer_t i;
    igraph_vector_ptr_t *vec;
    for (i = 0; i < 3; i++) {
        igraph_integer_t j;
        vec = ptr[i];
        for (j = 0; j < igraph_vector_ptr_size(vec); j++) {
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
                      igraph_gml_tree_line(node, pos));
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
                      igraph_gml_tree_line(node, pos));
    }

    return p;
}

igraph_error_t igraph_i_gml_parsedata_init(igraph_i_gml_parsedata_t* context) {
    context->eof = 0;
    context->depth = 0;
    context->scanner = 0;
    context->tree = 0;
    context->errmsg[0] = '\0';
    context->igraph_errno = IGRAPH_SUCCESS;

    return IGRAPH_SUCCESS;
}

void igraph_i_gml_parsedata_destroy(igraph_i_gml_parsedata_t* context) {
    if (context->tree != 0) {
        igraph_gml_tree_destroy(context->tree);
        context->tree = 0;
    }

    if (context->scanner != 0) {
        igraph_gml_yylex_destroy(context->scanner);
        context->scanner = 0;
    }
}

/**
 * \function igraph_read_graph_gml
 * \brief Read a graph in GML format.
 *
 * GML is a simple textual format, see
 * http://www.fim.uni-passau.de/en/fim/faculty/chairs/theoretische-informatik/projects.html for details.
 *
 * </para><para>
 * Although all syntactically correct GML can be parsed,
 * we implement only a subset of this format, some attributes might be
 * ignored. Here is a list of all the differences:
 * \olist
 * \oli Only <code>node</code> and <code>edge</code> attributes are
 *      used, and only if they have a simple type: integer, real or
 *      string. So if an attribute is an array or a record, then it is
 *      ignored. This is also true if only some values of the
 *      attribute are complex.
 * \oli Top level attributes except for <code>Version</code> and the
 *      first <code>graph</code> attribute are completely ignored.
 * \oli Graph attributes except for <code>node</code> and
 *      <code>edge</code> are completely ignored.
 * \oli There is no maximum line length.
 * \oli There is no maximum keyword length.
 * \oli Character entities in strings are not interpreted.
 * \oli We allow <code>inf</code> (infinity) and <code>nan</code>
 *      (not a number) as a real number. This is case insensitive, so
 *      <code>nan</code>, <code>NaN</code> and <code>NAN</code> are equal.
 * \endolist
 *
 * </para><para> Please contact us if you cannot live with these
 * limitations of the GML parser.
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

    igraph_integer_t i, p;
    igraph_integer_t no_of_nodes = 0, no_of_edges = 0;
    igraph_integer_t node_no;
    igraph_trie_t trie;
    igraph_vector_int_t edges;
    igraph_bool_t directed = IGRAPH_UNDIRECTED;
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
        if (context.errmsg[0] != 0) {
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

    /* Is it directed? */
    i = igraph_gml_tree_find(gtree, "directed", 0);
    if (i >= 0) {
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
        } else {
            IGRAPH_WARNINGF("Invalid type for 'directed' attribute on line %" IGRAPH_PRId ", assuming undirected.",
                            igraph_gml_tree_line(gtree, i));
        }
    }

    /* Now we go over all objects in the graph and collect the attribute names and
     * types. Plus we collect node IDs. We also do some checks.
     *
     * The 'id' fields of 'node' objects are converted into strings, so that they
     * can be inserted into a trie and re-encoded as consecutive integers starting
     * at 0. The GML spec allows isolated nodes with no 'id' field. These get a
     * generated string id of the form "n123" consisting of "n" and their count
     * (i.e. ordinal position) within the GML file.
     *
     * We use an attribute type value of -1 to mark attribute records which
     * correspond to composite GML values and must therefore be removed before
     * creating the graph.
     */
    node_no = 0;
    for (i = 0; i < igraph_gml_tree_length(gtree); i++) {
        char cname[100];
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
            hasid = 0;
            for (igraph_integer_t j = 0; j < igraph_gml_tree_length(node); j++) {
                const char *name = igraph_gml_tree_name(node, j);
                igraph_i_gml_tree_type_t type = igraph_gml_tree_type(node, j);
                igraph_integer_t trieid, triesize = igraph_trie_size(&vattrnames);
                IGRAPH_CHECK(igraph_trie_get(&vattrnames, name, &trieid));
                if (trieid == triesize) {
                    /* new attribute */
                    igraph_attribute_record_t *atrec = IGRAPH_CALLOC(1, igraph_attribute_record_t);
                    if (!atrec) {
                        IGRAPH_ERROR("Cannot read GML file.", IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */
                    }
                    IGRAPH_FINALLY(igraph_free, atrec);
                    atrec->name = strdup(name);
                    if (! atrec->name) {
                        IGRAPH_ERROR("Cannot read GML file.", IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */
                    }
                    IGRAPH_FINALLY(igraph_free, (char *) atrec->name);
                    if (type == IGRAPH_I_GML_TREE_INTEGER || type == IGRAPH_I_GML_TREE_REAL) {
                        atrec->type = IGRAPH_ATTRIBUTE_NUMERIC;
                    } else if (type == IGRAPH_I_GML_TREE_STRING) {
                        atrec->type = IGRAPH_ATTRIBUTE_STRING;
                    } else {
                        atrec->type = -1;
                    }
                    IGRAPH_CHECK(igraph_vector_ptr_push_back(&vattrs, atrec));
                    IGRAPH_FINALLY_CLEAN(2);
                } else {
                    /* already seen, should we update type? */
                    igraph_attribute_record_t *atrec = VECTOR(vattrs)[trieid];
                    igraph_attribute_type_t type1 = atrec->type;
                    if (type == IGRAPH_I_GML_TREE_STRING) {
                        atrec->type = IGRAPH_ATTRIBUTE_STRING;
                    } else if (type1 == -1) {
                        if (type == IGRAPH_I_GML_TREE_INTEGER || type == IGRAPH_I_GML_TREE_REAL) {
                            atrec->type = IGRAPH_ATTRIBUTE_NUMERIC;
                        }
                    }
                }
                /* check id */
                if (!strcmp(name, "id")) {
                    igraph_integer_t id, trie_id;
                    igraph_integer_t trie_size =igraph_trie_size(&trie);
                    if (hasid) {
                        /* A 'node' must not have more than one 'id' field. */
                        IGRAPH_ERRORF("Node has multiple id fields in GML file, line %" IGRAPH_PRId ".", IGRAPH_PARSEERROR,
                                      igraph_gml_tree_line(node, j));
                    }
                    if (igraph_gml_tree_type(node, j) != IGRAPH_I_GML_TREE_INTEGER) {
                        IGRAPH_ERRORF("Non-integer node id in GML file, line %" IGRAPH_PRId ".", IGRAPH_PARSEERROR,
                                      igraph_gml_tree_line(node, j));
                    }
                    id = igraph_gml_tree_get_integer(node, j);
                    snprintf(cname, sizeof(cname) / sizeof(char) - 1, "%" IGRAPH_PRId, id);
                    IGRAPH_CHECK(igraph_trie_get(&trie, cname, &trie_id));
                    if (trie_id != trie_size) {
                        /* This id has already been seen in a previous node. */
                        IGRAPH_ERRORF("Duplicate node id in GML file, line %" IGRAPH_PRId ".", IGRAPH_PARSEERROR,
                                      igraph_gml_tree_line(node, j));
                    }
                    hasid = 1;
                }
            }
            if (!hasid) {
                /* Isolated nodes are allowed not to have an id.
                 * We generate an "n"-prefixed string id to be used in the trie. */
                igraph_integer_t trie_id;
                snprintf(cname, sizeof(cname) / sizeof(char) - 1, "n%" IGRAPH_PRId, node_no);
                IGRAPH_CHECK(igraph_trie_get(&trie, cname, &trie_id));
            }
        } else if (!strcmp(name, "edge")) {
            igraph_gml_tree_t *edge;
            igraph_bool_t has_source = 0, has_target = 0;
            no_of_edges++;
            if (igraph_gml_tree_type(gtree, i) != IGRAPH_I_GML_TREE_TREE) {
                IGRAPH_ERRORF("'edge' is not a list in GML file, line %" IGRAPH_PRId ".", IGRAPH_PARSEERROR,
                              igraph_gml_tree_line(gtree, i));
            }
            edge = igraph_gml_tree_get_tree(gtree, i);
            has_source = has_target = 0;
            for (igraph_integer_t j = 0; j < igraph_gml_tree_length(edge); j++) {
                const char *name = igraph_gml_tree_name(edge, j);
                if (!strcmp(name, "source")) {
                    has_source = 1;
                    if (igraph_gml_tree_type(edge, j) != IGRAPH_I_GML_TREE_INTEGER) {
                        IGRAPH_ERRORF("Non-integer 'source' for an edge in GML file, line %" IGRAPH_PRId ".",
                                      IGRAPH_PARSEERROR,
                                      igraph_gml_tree_line(edge, j));
                    }
                } else if (!strcmp(name, "target")) {
                    has_target = 1;
                    if (igraph_gml_tree_type(edge, j) != IGRAPH_I_GML_TREE_INTEGER) {
                        IGRAPH_ERRORF("Non-integer 'target' for an edge in GML file, line %" IGRAPH_PRId ".",
                                      IGRAPH_PARSEERROR,
                                      igraph_gml_tree_line(edge, j));
                    }
                } else {
                    igraph_integer_t trieid, triesize = igraph_trie_size(&eattrnames);
                    igraph_i_gml_tree_type_t type = igraph_gml_tree_type(edge, j);
                    IGRAPH_CHECK(igraph_trie_get(&eattrnames, name, &trieid));
                    if (trieid == triesize) {
                        /* new attribute */
                        igraph_attribute_record_t *atrec = IGRAPH_CALLOC(1, igraph_attribute_record_t);
                        if (!atrec) {
                            IGRAPH_ERROR("Cannot read GML file.", IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */
                        }
                        IGRAPH_FINALLY(igraph_free, atrec);
                        atrec->name = strdup(name);
                        if (! atrec->name) {
                            IGRAPH_ERROR("Cannot read GML file.", IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */
                        }
                        IGRAPH_FINALLY(igraph_free, (char *) atrec->name);
                        if (type == IGRAPH_I_GML_TREE_INTEGER || type == IGRAPH_I_GML_TREE_REAL) {
                            atrec->type = IGRAPH_ATTRIBUTE_NUMERIC;
                        } else if (type == IGRAPH_I_GML_TREE_STRING) {
                            atrec->type = IGRAPH_ATTRIBUTE_STRING;
                        } else {
                            atrec->type = -1;
                        }
                        IGRAPH_CHECK(igraph_vector_ptr_push_back(&eattrs, atrec));
                        IGRAPH_FINALLY_CLEAN(2);
                    } else {
                        /* already seen, should we update type? */
                        igraph_attribute_record_t *atrec = VECTOR(eattrs)[trieid];
                        igraph_attribute_type_t type1 = atrec->type;
                        if (type == IGRAPH_I_GML_TREE_STRING) {
                            atrec->type = IGRAPH_ATTRIBUTE_STRING;
                        } else if (type1 == -1) {
                            if (type == IGRAPH_I_GML_TREE_INTEGER || type == IGRAPH_I_GML_TREE_REAL) {
                                atrec->type = IGRAPH_ATTRIBUTE_NUMERIC;
                            }
                        }
                    }
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
        } else {
            /* anything to do? Maybe add as graph attribute.... */
        }
    }

    /* At this point, all nodes must have an id (from the file or generated) stored
     * in the trie. Any condition that violates this should have been caught during
     * the preceding checks. */
    IGRAPH_ASSERT(igraph_trie_size(&trie) == no_of_nodes);

    /* Now we allocate the vectors and strvectors for the attributes */
    for (i = 0; i < igraph_vector_ptr_size(&vattrs); i++) {
        igraph_attribute_record_t *atrec = VECTOR(vattrs)[i];
        igraph_attribute_type_t type = atrec->type;
        if (type == IGRAPH_ATTRIBUTE_NUMERIC) {
            igraph_vector_t *p = IGRAPH_CALLOC(1, igraph_vector_t);
            if (! p) {
                IGRAPH_ERROR("Cannot read GML file.", IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */
            }
            IGRAPH_FINALLY(igraph_free, p);
            IGRAPH_CHECK(igraph_vector_init(p, no_of_nodes));
            igraph_vector_fill(p, IGRAPH_NAN); /* use NaN as default */
            atrec->value = p;
            IGRAPH_FINALLY_CLEAN(1);
        } else if (type == IGRAPH_ATTRIBUTE_STRING) {
            igraph_strvector_t *p = IGRAPH_CALLOC(1, igraph_strvector_t);
            if (! p) {
                IGRAPH_ERROR("Cannot read GML file.", IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */
            }
            IGRAPH_FINALLY(igraph_free, p);
            IGRAPH_CHECK(igraph_strvector_init(p, no_of_nodes));
            atrec->value = p;
            IGRAPH_FINALLY_CLEAN(1);
        } else {
            IGRAPH_WARNINGF("Composite vertex attribute '%s' ignored in GML file.", atrec->name);
        }
    }

    for (i = 0; i < igraph_vector_ptr_size(&eattrs); i++) {
        igraph_attribute_record_t *atrec = VECTOR(eattrs)[i];
        igraph_attribute_type_t  type = atrec->type;
        if (type == IGRAPH_ATTRIBUTE_NUMERIC) {
            igraph_vector_t *p = IGRAPH_CALLOC(1, igraph_vector_t);
            if (! p) {
                IGRAPH_ERROR("Cannot read GML file.", IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */
            }
            IGRAPH_FINALLY(igraph_free, p);
            IGRAPH_CHECK(igraph_vector_init(p, no_of_edges));
            igraph_vector_fill(p, IGRAPH_NAN); /* use NaN as default */
            atrec->value = p;
            IGRAPH_FINALLY_CLEAN(1);
        } else if (type == IGRAPH_ATTRIBUTE_STRING) {
            igraph_strvector_t *p = IGRAPH_CALLOC(1, igraph_strvector_t);
            if (! p) {
                IGRAPH_ERROR("Cannot read GML file.", IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */
            }
            IGRAPH_FINALLY(igraph_free, p);
            IGRAPH_CHECK(igraph_strvector_init(p, no_of_edges));
            atrec->value = p;
            IGRAPH_FINALLY_CLEAN(1);
        } else {
            IGRAPH_WARNINGF("Composite edge attribute '%s' ignored in GML file.", atrec->name);
        }
    }

    /* Ok, now the edges, attributes too */
    IGRAPH_VECTOR_INT_INIT_FINALLY(&edges, no_of_edges * 2);
    p = -1;
    while ( (p = igraph_gml_tree_find(gtree, "edge", p + 1)) != -1) {
        igraph_gml_tree_t *edge;
        igraph_integer_t from, to, fromidx = 0, toidx = 0;
        char name[100];
        igraph_integer_t j;
        edge = igraph_gml_tree_get_tree(gtree, p);
        for (j = 0; j < igraph_gml_tree_length(edge); j++) {
            const char *n = igraph_gml_tree_name(edge, j);
            if (!strcmp(n, "source")) {
                fromidx = igraph_gml_tree_find(edge, "source", 0);
            } else if (!strcmp(n, "target")) {
                toidx = igraph_gml_tree_find(edge, "target", 0);
            } else {
                igraph_integer_t edgeid = edgeptr / 2;
                igraph_integer_t trieidx;
                igraph_attribute_record_t *atrec;
                igraph_attribute_type_t type;
                IGRAPH_CHECK(igraph_trie_get(&eattrnames, n, &trieidx));
                atrec = VECTOR(eattrs)[trieidx];
                type = atrec->type;
                if (type == IGRAPH_ATTRIBUTE_NUMERIC) {
                    igraph_vector_t *v = (igraph_vector_t *) atrec->value;
                    VECTOR(*v)[edgeid] = igraph_i_gml_toreal(edge, j);
                } else if (type == IGRAPH_ATTRIBUTE_STRING) {
                    igraph_strvector_t *v = (igraph_strvector_t *) atrec->value;
                    const char *value = igraph_i_gml_tostring(edge, j);
                    IGRAPH_CHECK(igraph_strvector_set(v, edgeid, value));
                } else {
                    /* Ignored composite attribute */
                }
            }
        }
        from = igraph_gml_tree_get_integer(edge, fromidx);
        to = igraph_gml_tree_get_integer(edge, toidx);
        snprintf(name, sizeof(name) / sizeof(char) - 1, "%" IGRAPH_PRId, from);
        IGRAPH_CHECK(igraph_trie_check(&trie, name, &from));
        if (from < 0) {
            IGRAPH_ERRORF("Unknown source node id found in an edge in GML file, line %" IGRAPH_PRId ".",
                         IGRAPH_PARSEERROR, igraph_gml_tree_line(edge, fromidx));
        }
        snprintf(name, sizeof(name) / sizeof(char) - 1, "%" IGRAPH_PRId, to);
        IGRAPH_CHECK(igraph_trie_check(&trie, name, &to));
        if (to < 0) {
            IGRAPH_ERRORF("Unknown target node id found in an edge in GML file, line %" IGRAPH_PRId ".",
                         IGRAPH_PARSEERROR, igraph_gml_tree_line(edge, toidx));
        }
        VECTOR(edges)[edgeptr++] = from;
        VECTOR(edges)[edgeptr++] = to;
    }

    /* and add vertex attributes */
    node_no = 0;
    for (i = 0; i < igraph_gml_tree_length(gtree); i++) {
        const char *n;
        char name[100];
        igraph_integer_t j, k;
        n = igraph_gml_tree_name(gtree, i);
        if (!strcmp(n, "node")) {
            igraph_gml_tree_t *node = igraph_gml_tree_get_tree(gtree, i);
            igraph_integer_t iidx = igraph_gml_tree_find(node, "id", 0);
            igraph_integer_t trie_id;
            node_no++;
            if (iidx < 0) {
                /* Isolated node with no id field, use generated id */
                snprintf(name, sizeof(name) / sizeof(char) - 1, "n%" IGRAPH_PRId, node_no);
            } else {
                igraph_integer_t id;
                id = igraph_gml_tree_get_integer(node, iidx);
                snprintf(name, sizeof(name) / sizeof(char) - 1, "%" IGRAPH_PRId, id);
            }
            IGRAPH_CHECK(igraph_trie_get(&trie, name, &trie_id));
            for (j = 0; j < igraph_gml_tree_length(node); j++) {
                const char *aname = igraph_gml_tree_name(node, j);
                igraph_attribute_record_t *atrec;
                igraph_attribute_type_t type;
                IGRAPH_CHECK(igraph_trie_get(&vattrnames, aname, &k));
                atrec = VECTOR(vattrs)[k];
                type = atrec->type;
                if (type == IGRAPH_ATTRIBUTE_NUMERIC) {
                    igraph_vector_t *v = (igraph_vector_t *)atrec->value;
                    VECTOR(*v)[trie_id] = igraph_i_gml_toreal(node, j);
                } else if (type == IGRAPH_ATTRIBUTE_STRING) {
                    igraph_strvector_t *v = (igraph_strvector_t *)atrec->value;
                    const char *value = igraph_i_gml_tostring(node, j);
                    IGRAPH_CHECK(igraph_strvector_set(v, trie_id, value));
                } else {
                    /* Ignored composite attribute */
                }
            }
        }
    }

    /* Remove composite attributes */
    {
        igraph_integer_t j;

        for (i = 0, j = 0; i < igraph_vector_ptr_size(&vattrs); i++) {
            igraph_attribute_record_t *atrec = VECTOR(vattrs)[i];
            if (atrec->type == -1) {
                IGRAPH_FREE(atrec->name);
                IGRAPH_FREE(atrec);
            } else {
                VECTOR(vattrs)[j++] = VECTOR(vattrs)[i];
            }
        }
        igraph_vector_ptr_resize(&vattrs, j);

        for (i = 0, j = 0; i < igraph_vector_ptr_size(&eattrs); i++) {
            igraph_attribute_record_t *atrec = VECTOR(eattrs)[i];
            if (atrec->type == -1) {
                IGRAPH_FREE(atrec->name);
                IGRAPH_FREE(atrec);
            } else {
                VECTOR(eattrs)[j++] = VECTOR(eattrs)[i];
            }
        }
        igraph_vector_ptr_resize(&eattrs, j);
    }

    igraph_trie_destroy(&trie);
    igraph_trie_destroy(&gattrnames);
    igraph_trie_destroy(&vattrnames);
    igraph_trie_destroy(&eattrnames);
    IGRAPH_FINALLY_CLEAN(4);

    IGRAPH_CHECK(igraph_empty_attrs(graph, 0, directed, 0)); /* TODO https://github.com/igraph/igraph/issues/174 */
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
    int no = 1;
    char strno[50];
    size_t i, len = strlen(orig), newlen = 0, plen = 0;

    /* do we need a prefix? */
    if (len == 0 || !isalpha(orig[0])) {
        no++;
        snprintf(strno, sizeof(strno) - 1, "igraph");
        plen = newlen = strlen(strno);
    }
    for (i = 0; i < len; i++) {
        if (isalnum(orig[i])) {
            newlen++;
        }
    }
    *key = IGRAPH_CALLOC(newlen + 1, char);
    if (! *key) {
        IGRAPH_ERROR("Writing GML format failed.", IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */
    }
    memcpy(*key, strno, plen * sizeof(char));
    for (i = 0; i < len; i++) {
        if (isalnum(orig[i])) {
            (*key)[plen++] = orig[i];
        }
    }
    (*key)[newlen] = '\0';

    return IGRAPH_SUCCESS;
}

#define CHECK(cmd) do { ret=cmd; if (ret<0) IGRAPH_ERROR("Writing GML format failed.", IGRAPH_EFILE); } while (0)

/**
 * \function igraph_write_graph_gml
 * \brief Write the graph to a stream in GML format.
 *
 * GML is a quite general textual format, see
 * http://www.fim.uni-passau.de/en/fim/faculty/chairs/theoretische-informatik/projects.html for details.
 *
 * </para><para> The graph, vertex and edges attributes are written to the
 * file as well, if they are numeric or string.
 *
 * </para><para> As igraph is more forgiving about attribute names, it might
 * be necessary to simplify the them before writing to the GML file.
 * This way we'll have a syntactically correct GML file. The following
 * simple procedure is performed on each attribute name: first the alphanumeric
 * characters are extracted, the others are ignored. Then if the first character
 * is not a letter then the attribute name is prefixed with <quote>igraph</quote>.
 * Note that this might result identical names for two attributes, igraph
 * does not check this.
 *
 * </para><para> The <quote>id</quote> vertex attribute is treated specially.
 * If the <parameter>id</parameter> argument is not \c NULL then it should be a numeric
 * vector with the vertex IDs and the <quote>id</quote> vertex attribute is
 * ignored (if there is one). If <parameter>id</parameter> is \c NULL and there is a
 * numeric <quote>id</quote> vertex attribute, it will be used instead. If ids
 * are not specified in either way then the regular igraph vertex IDs are used.
 * If some of the supplied id values are invalid (non-integer or NaN), all supplied
 * id are ignored and igraph vertex IDs are used instead.
 *
 * </para><para> Note that whichever way vertex IDs are specified, their
 * uniqueness is not checked.
 *
 * </para><para> If the graph has edge attributes named <quote>source</quote>
 * or <quote>target</quote> they're silently ignored. GML uses these attributes
 * to specify the edges, so we cannot write them to the file. Rename them
 * before calling this function if you want to preserve them.
 *
 * \param graph The graph to write to the stream.
 * \param outstream The stream to write the file to.
 * \param id Either <code>NULL</code> or a numeric vector with the vertex IDs.
 *        See details above.
 * \param creator An optional string to append to the creator line.
 *        If this is \c NULL then the current date and time is added.
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
                           const igraph_vector_t *id, const char *creator) {
    igraph_error_t ret;
    igraph_strvector_t gnames, vnames, enames;
    igraph_vector_int_t gtypes, vtypes, etypes;
    igraph_vector_t numv;
    igraph_strvector_t strv;
    igraph_vector_bool_t boolv;
    igraph_integer_t i;
    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_integer_t no_of_edges = igraph_ecount(graph);

    igraph_vector_t v_myid;
    const igraph_vector_t *myid = id;

    time_t curtime = time(0);
    char *timestr = ctime(&curtime);
    timestr[strlen(timestr) - 1] = '\0'; /* nicely remove \n */

    CHECK(fprintf(outstream,
                  "Creator \"igraph version %s %s\"\nVersion 1\ngraph\n[\n",
                  IGRAPH_VERSION, creator ? creator : timestr));

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

    IGRAPH_VECTOR_INIT_FINALLY(&numv, 1);
    IGRAPH_STRVECTOR_INIT_FINALLY(&strv, 1);
    IGRAPH_VECTOR_BOOL_INIT_FINALLY(&boolv, 1);

    /* Check whether there is an 'id' node attribute if the supplied is 0 */
    if (!id) {
        igraph_bool_t found = 0;
        for (i = 0; i < igraph_vector_int_size(&vtypes); i++) {
            const char *n = igraph_strvector_get(&vnames, i);
            if (!strcmp(n, "id") && VECTOR(vtypes)[i] == IGRAPH_ATTRIBUTE_NUMERIC) {
                found = 1; break;
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
     * nodes did not have an id. In this case, the "id" attribute created by
     * igraph will contain a NaN value.
     *
     * TODO: Check that ids are unique?
     */
    if (myid) {
        if (igraph_vector_size(myid) != no_of_nodes) {
            IGRAPH_ERROR("Size of id vector must match vertex count.", IGRAPH_EINVAL);
        }
        for (i = 0; i < no_of_nodes; ++i) {
            igraph_real_t val = VECTOR(*myid)[i];
            if (val != (igraph_integer_t) val) {
                IGRAPH_WARNINGF("%g is not a valid integer id for GML files, ignoring all supplied ids.", val);
                myid = NULL;
                break;
            }
        }
    }

    /* directedness */
    CHECK(fprintf(outstream, "  directed %i\n", igraph_is_directed(graph) ? 1 : 0));

    /* Graph attributes first */
    for (i = 0; i < igraph_vector_int_size(&gtypes); i++) {
        const char *name;
        char *newname;
        name = igraph_strvector_get(&gnames, i);
        IGRAPH_CHECK(igraph_i_gml_convert_to_key(name, &newname));
        IGRAPH_FINALLY(igraph_free, newname);
        if (VECTOR(gtypes)[i] == IGRAPH_ATTRIBUTE_NUMERIC) {
            IGRAPH_CHECK(igraph_i_attribute_get_numeric_graph_attr(graph, name, &numv));
            CHECK(fprintf(outstream, "  %s ", newname));
            CHECK(igraph_real_fprintf_precise(outstream, VECTOR(numv)[0]));
            CHECK(fputc('\n', outstream));
        } else if (VECTOR(gtypes)[i] == IGRAPH_ATTRIBUTE_STRING) {
            const char *s;
            IGRAPH_CHECK(igraph_i_attribute_get_string_graph_attr(graph, name, &strv));
            s = igraph_strvector_get(&strv, 0);
            CHECK(fprintf(outstream, "  %s \"%s\"\n", newname, s));
        } else if (VECTOR(gtypes)[i] == IGRAPH_ATTRIBUTE_BOOLEAN) {
            IGRAPH_CHECK(igraph_i_attribute_get_bool_graph_attr(graph, name, &boolv));
            CHECK(fprintf(outstream, "  %s %d\n", newname, VECTOR(boolv)[0] ? 1 : 0));
            IGRAPH_WARNING("A boolean graph attribute was converted to numeric");
        } else {
            IGRAPH_WARNING("A non-numeric, non-string, non-boolean graph attribute ignored");
        }
        IGRAPH_FREE(newname);
        IGRAPH_FINALLY_CLEAN(1);
    }

    /* Now come the vertices */
    for (i = 0; i < no_of_nodes; i++) {
        igraph_integer_t j;
        CHECK(fprintf(outstream, "  node\n  [\n"));
        /* id */
        CHECK(fprintf(outstream, "    id %" IGRAPH_PRId "\n", myid ? (igraph_integer_t)VECTOR(*myid)[i] : i));
        /* other attributes */
        for (j = 0; j < igraph_vector_int_size(&vtypes); j++) {
            igraph_attribute_type_t type = (igraph_attribute_type_t) VECTOR(vtypes)[j];
            const char *name;
            char *newname;
            name = igraph_strvector_get(&vnames, j);
            if (!strcmp(name, "id")) {
                continue;
            }
            IGRAPH_CHECK(igraph_i_gml_convert_to_key(name, &newname));
            IGRAPH_FINALLY(igraph_free, newname);
            if (type == IGRAPH_ATTRIBUTE_NUMERIC) {
                IGRAPH_CHECK(igraph_i_attribute_get_numeric_vertex_attr(graph, name,
                             igraph_vss_1(i), &numv));
                CHECK(fprintf(outstream, "    %s ", newname));
                CHECK(igraph_real_fprintf_precise(outstream, VECTOR(numv)[0]));
                CHECK(fputc('\n', outstream));
            } else if (type == IGRAPH_ATTRIBUTE_STRING) {
                const char *s;
                IGRAPH_CHECK(igraph_i_attribute_get_string_vertex_attr(graph, name,
                             igraph_vss_1(i), &strv));
                s = igraph_strvector_get(&strv, 0);
                CHECK(fprintf(outstream, "    %s \"%s\"\n", newname, s));
            } else if (type == IGRAPH_ATTRIBUTE_BOOLEAN) {
                IGRAPH_CHECK(igraph_i_attribute_get_bool_vertex_attr(graph, name,
                             igraph_vss_1(i), &boolv));
                CHECK(fprintf(outstream, "    %s %d\n", newname, VECTOR(boolv)[0] ? 1 : 0));
                IGRAPH_WARNING("A boolean vertex attribute was converted to numeric");
            } else {
                IGRAPH_WARNING("A non-numeric, non-string, non-boolean edge attribute was ignored");
            }
            IGRAPH_FREE(newname);
            IGRAPH_FINALLY_CLEAN(1);
        }
        CHECK(fprintf(outstream, "  ]\n"));
    }

    /* The edges too */
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
        for (j = 0; j < igraph_vector_int_size(&etypes); j++) {
            igraph_attribute_type_t type = (igraph_attribute_type_t) VECTOR(etypes)[j];
            const char *name;
            char *newname;
            name = igraph_strvector_get(&enames, j);
            if (!strcmp(name, "source") || !strcmp(name, "target")) {
                continue;
            }
            IGRAPH_CHECK(igraph_i_gml_convert_to_key(name, &newname));
            IGRAPH_FINALLY(igraph_free, newname);
            if (type == IGRAPH_ATTRIBUTE_NUMERIC) {
                IGRAPH_CHECK(igraph_i_attribute_get_numeric_edge_attr(graph, name,
                             igraph_ess_1(i), &numv));
                CHECK(fprintf(outstream, "    %s ", newname));
                CHECK(igraph_real_fprintf_precise(outstream, VECTOR(numv)[0]));
                CHECK(fputc('\n', outstream));
            } else if (type == IGRAPH_ATTRIBUTE_STRING) {
                const char *s;
                IGRAPH_CHECK(igraph_i_attribute_get_string_edge_attr(graph, name,
                             igraph_ess_1(i), &strv));
                s = igraph_strvector_get(&strv, 0);
                CHECK(fprintf(outstream, "    %s \"%s\"\n", newname, s));
            } else if (type == IGRAPH_ATTRIBUTE_BOOLEAN) {
                IGRAPH_CHECK(igraph_i_attribute_get_bool_edge_attr(graph, name,
                             igraph_ess_1(i), &boolv));
                CHECK(fprintf(outstream, "    %s %d\n", newname, VECTOR(boolv)[0] ? 1 : 0));
                IGRAPH_WARNING("A boolean edge attribute was converted to numeric");
            } else {
                IGRAPH_WARNING("A non-numeric, non-string, non-boolean edge attribute was ignored");
            }
            IGRAPH_FREE(newname);
            IGRAPH_FINALLY_CLEAN(1);
        }
        CHECK(fprintf(outstream, "  ]\n"));
    }

    CHECK(fprintf(outstream, "]\n"));

    if (&v_myid == myid) {
        igraph_vector_destroy(&v_myid);
        IGRAPH_FINALLY_CLEAN(1);
    }

    igraph_vector_bool_destroy(&boolv);
    igraph_strvector_destroy(&strv);
    igraph_vector_destroy(&numv);
    igraph_vector_int_destroy(&etypes);
    igraph_vector_int_destroy(&vtypes);
    igraph_vector_int_destroy(&gtypes);
    igraph_strvector_destroy(&enames);
    igraph_strvector_destroy(&vnames);
    igraph_strvector_destroy(&gnames);
    IGRAPH_FINALLY_CLEAN(9);

    return IGRAPH_SUCCESS;
}

#undef CHECK
