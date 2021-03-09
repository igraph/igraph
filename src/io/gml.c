/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2005-2020  The igraph development team

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
    long int i;
    igraph_vector_ptr_t *vec;
    for (i = 0; i < 3; i++) {
        long int j;
        vec = ptr[i];
        for (j = 0; j < igraph_vector_ptr_size(vec); j++) {
            igraph_attribute_record_t *atrec = VECTOR(*vec)[j];
            if (atrec->type == IGRAPH_ATTRIBUTE_NUMERIC) {
                igraph_vector_t *value = (igraph_vector_t*)atrec->value;
                if (value != 0) {
                    igraph_vector_destroy(value);
                    IGRAPH_FREE(value);
                }
            } else {
                igraph_strvector_t *value = (igraph_strvector_t*)atrec->value;
                if (value != 0) {
                    igraph_strvector_destroy(value);
                    IGRAPH_FREE(value);
                }
            }
            IGRAPH_FREE(atrec->name);
            IGRAPH_FREE(atrec);
        }
        igraph_vector_ptr_destroy(vec);
    }
}

static int igraph_i_gml_toreal(igraph_gml_tree_t *node, long int pos, igraph_real_t *result) {

    igraph_real_t value = 0.0;
    int type = igraph_gml_tree_type(node, pos);

    switch (type) {
    case IGRAPH_I_GML_TREE_INTEGER:
        value = igraph_gml_tree_get_integer(node, pos);
        break;
    case IGRAPH_I_GML_TREE_REAL:
        value = igraph_gml_tree_get_real(node, pos);
        break;
    default:
        IGRAPH_ERROR("Internal error while parsing GML file.", IGRAPH_FAILURE);
        break;
    }

    *result = value;
    return IGRAPH_SUCCESS;
}

static const char *igraph_i_gml_tostring(igraph_gml_tree_t *node, long int pos) {

    int type = igraph_gml_tree_type(node, pos);
    static char tmp[256];
    const char *p = tmp;
    long int i;
    igraph_real_t d;

    switch (type) {
    case IGRAPH_I_GML_TREE_INTEGER:
        i = igraph_gml_tree_get_integer(node, pos);
        snprintf(tmp, sizeof(tmp) / sizeof(char), "%li", i);
        break;
    case IGRAPH_I_GML_TREE_REAL:
        d = igraph_gml_tree_get_real(node, pos);
        igraph_real_snprintf_precise(tmp, sizeof(tmp) / sizeof(char), d);
        break;
    case IGRAPH_I_GML_TREE_STRING:
        p = igraph_gml_tree_get_string(node, pos);
        break;
    default:
        break;
    }

    return p;
}

int igraph_i_gml_parsedata_init(igraph_i_gml_parsedata_t* context) {
    context->eof = 0;
    context->scanner = 0;
    context->tree = 0;
    context->errmsg[0] = 0;

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
int igraph_read_graph_gml(igraph_t *graph, FILE *instream) {

    long int i, p;
    long int no_of_nodes = 0, no_of_edges = 0;
    igraph_trie_t trie;
    igraph_vector_t edges;
    igraph_bool_t directed = IGRAPH_UNDIRECTED;
    igraph_gml_tree_t *gtree;
    long int gidx;
    igraph_trie_t vattrnames;
    igraph_trie_t eattrnames;
    igraph_trie_t gattrnames;
    igraph_vector_ptr_t gattrs = IGRAPH_VECTOR_PTR_NULL,
                        vattrs = IGRAPH_VECTOR_PTR_NULL, eattrs = IGRAPH_VECTOR_PTR_NULL;
    igraph_vector_ptr_t *attrs[3];
    long int edgeptr = 0;
    igraph_i_gml_parsedata_t context;

    attrs[0] = &gattrs; attrs[1] = &vattrs; attrs[2] = &eattrs;

    IGRAPH_CHECK(igraph_i_gml_parsedata_init(&context));
    IGRAPH_FINALLY(igraph_i_gml_parsedata_destroy, &context);

    igraph_gml_yylex_init_extra(&context, &context.scanner);

    igraph_gml_yyset_in(instream, context.scanner);

    i = igraph_gml_yyparse(&context);
    if (i != 0) {
        if (context.errmsg[0] != 0) {
            IGRAPH_ERROR(context.errmsg, IGRAPH_PARSEERROR);
        } else {
            IGRAPH_ERROR("Cannot read GML file.", IGRAPH_PARSEERROR);
        }
    }

    IGRAPH_VECTOR_INIT_FINALLY(&edges, 0);

    /* Check version, if present, integer and not '1' then ignored */
    i = igraph_gml_tree_find(context.tree, "Version", 0);
    if (i >= 0 &&
        igraph_gml_tree_type(context.tree, i) == IGRAPH_I_GML_TREE_INTEGER &&
        igraph_gml_tree_get_integer(context.tree, i) != 1) {
        IGRAPH_ERROR("Unknown GML version.", IGRAPH_UNIMPLEMENTED);
        /* RETURN HERE!!!! */
    }

    /* get the graph */
    gidx = igraph_gml_tree_find(context.tree, "graph", 0);
    if (gidx == -1) {
        IGRAPH_ERROR("No 'graph' object in GML file.", IGRAPH_PARSEERROR);
    }
    if (igraph_gml_tree_type(context.tree, gidx) !=
        IGRAPH_I_GML_TREE_TREE) {
        IGRAPH_ERROR("Invalid type for 'graph' object in GML file.", IGRAPH_PARSEERROR);
    }
    gtree = igraph_gml_tree_get_tree(context.tree, gidx);

    IGRAPH_FINALLY(igraph_i_gml_destroy_attrs, attrs);
    igraph_vector_ptr_init(&gattrs, 0);
    igraph_vector_ptr_init(&vattrs, 0);
    igraph_vector_ptr_init(&eattrs, 0);

    IGRAPH_TRIE_INIT_FINALLY(&trie, 0);
    IGRAPH_TRIE_INIT_FINALLY(&vattrnames, 0);
    IGRAPH_TRIE_INIT_FINALLY(&eattrnames, 0);
    IGRAPH_TRIE_INIT_FINALLY(&gattrnames, 0);

    /* Is is directed? */
    i = igraph_gml_tree_find(gtree, "directed", 0);
    if (i >= 0 && igraph_gml_tree_type(gtree, i) == IGRAPH_I_GML_TREE_INTEGER) {
        if (igraph_gml_tree_get_integer(gtree, i) == 1) {
            directed = IGRAPH_DIRECTED;
        }
    }

    /* Now we go over all objects in the graph and collect the attribute names and
       types. Plus we collect node ids. We also do some checks. */
    for (i = 0; i < igraph_gml_tree_length(gtree); i++) {
        long int j;
        char cname[100];
        const char *name = igraph_gml_tree_name(gtree, i);
        if (!strcmp(name, "node")) {
            igraph_gml_tree_t *node;
            igraph_bool_t hasid;
            no_of_nodes++;
            if (igraph_gml_tree_type(gtree, i) != IGRAPH_I_GML_TREE_TREE) {
                IGRAPH_ERROR("'node' is not a list in GML file.", IGRAPH_PARSEERROR);
            }
            node = igraph_gml_tree_get_tree(gtree, i);
            hasid = 0;
            for (j = 0; j < igraph_gml_tree_length(node); j++) {
                const char *name = igraph_gml_tree_name(node, j);
                long int trieid, triesize = igraph_trie_size(&vattrnames);
                IGRAPH_CHECK(igraph_trie_get(&vattrnames, name, &trieid));
                if (trieid == triesize) {
                    /* new attribute */
                    igraph_attribute_record_t *atrec = IGRAPH_CALLOC(1, igraph_attribute_record_t);
                    int type = igraph_gml_tree_type(node, j);
                    if (!atrec) {
                        IGRAPH_ERROR("Cannot read GML file.", IGRAPH_ENOMEM);
                    }
                    IGRAPH_CHECK(igraph_vector_ptr_push_back(&vattrs, atrec));
                    atrec->name = strdup(name);
                    if (type == IGRAPH_I_GML_TREE_INTEGER || type == IGRAPH_I_GML_TREE_REAL) {
                        atrec->type = IGRAPH_ATTRIBUTE_NUMERIC;
                    } else {
                        atrec->type = IGRAPH_ATTRIBUTE_STRING;
                    }
                } else {
                    /* already seen, should we update type? */
                    igraph_attribute_record_t *atrec = VECTOR(vattrs)[trieid];
                    int type1 = atrec->type;
                    int type2 = igraph_gml_tree_type(node, j);
                    if (type1 == IGRAPH_ATTRIBUTE_NUMERIC && type2 == IGRAPH_I_GML_TREE_STRING) {
                        atrec->type = IGRAPH_ATTRIBUTE_STRING;
                    }
                }
                /* check id */
                if (!hasid && !strcmp(name, "id")) {
                    long int id;
                    if (igraph_gml_tree_type(node, j) != IGRAPH_I_GML_TREE_INTEGER) {
                        IGRAPH_ERROR("Non-integer node id in GML file.", IGRAPH_PARSEERROR);
                    }
                    id = igraph_gml_tree_get_integer(node, j);
                    snprintf(cname, sizeof(cname) / sizeof(char) -1, "%li", id);
                    IGRAPH_CHECK(igraph_trie_get(&trie, cname, &id));
                    hasid = 1;
                }
            }
            if (!hasid) {
                IGRAPH_ERROR("Node without 'id' while parsing GML file.", IGRAPH_PARSEERROR);
            }
        } else if (!strcmp(name, "edge")) {
            igraph_gml_tree_t *edge;
            igraph_bool_t has_source = 0, has_target = 0;
            no_of_edges++;
            if (igraph_gml_tree_type(gtree, i) != IGRAPH_I_GML_TREE_TREE) {
                IGRAPH_ERROR("'edge' is not a list in GML file.", IGRAPH_PARSEERROR);
            }
            edge = igraph_gml_tree_get_tree(gtree, i);
            has_source = has_target = 0;
            for (j = 0; j < igraph_gml_tree_length(edge); j++) {
                const char *name = igraph_gml_tree_name(edge, j);
                if (!strcmp(name, "source")) {
                    has_source = 1;
                    if (igraph_gml_tree_type(edge, j) != IGRAPH_I_GML_TREE_INTEGER) {
                        IGRAPH_ERROR("Non-integer 'source' for an edge in GML file.",
                                     IGRAPH_PARSEERROR);
                    }
                } else if (!strcmp(name, "target")) {
                    has_target = 1;
                    if (igraph_gml_tree_type(edge, j) != IGRAPH_I_GML_TREE_INTEGER) {
                        IGRAPH_ERROR("Non-integer 'source' for an edge in GML file.",
                                     IGRAPH_PARSEERROR);
                    }
                } else {
                    long int trieid, triesize = igraph_trie_size(&eattrnames);
                    IGRAPH_CHECK(igraph_trie_get(&eattrnames, name, &trieid));
                    if (trieid == triesize) {
                        /* new attribute */
                        igraph_attribute_record_t *atrec = IGRAPH_CALLOC(1, igraph_attribute_record_t);
                        int type = igraph_gml_tree_type(edge, j);
                        if (!atrec) {
                            IGRAPH_ERROR("Cannot read GML file.", IGRAPH_ENOMEM);
                        }
                        IGRAPH_CHECK(igraph_vector_ptr_push_back(&eattrs, atrec));
                        atrec->name = strdup(name);
                        if (type == IGRAPH_I_GML_TREE_INTEGER || type == IGRAPH_I_GML_TREE_REAL) {
                            atrec->type = IGRAPH_ATTRIBUTE_NUMERIC;
                        } else {
                            atrec->type = IGRAPH_ATTRIBUTE_STRING;
                        }
                    } else {
                        /* already seen, should we update type? */
                        igraph_attribute_record_t *atrec = VECTOR(eattrs)[trieid];
                        int type1 = atrec->type;
                        int type2 = igraph_gml_tree_type(edge, j);
                        if (type1 == IGRAPH_ATTRIBUTE_NUMERIC && type2 == IGRAPH_I_GML_TREE_STRING) {
                            atrec->type = IGRAPH_ATTRIBUTE_STRING;
                        }
                    }
                }
            } /* for */
            if (!has_source) {
                IGRAPH_ERROR("No 'source' for edge in GML file.", IGRAPH_PARSEERROR);
            }
            if (!has_target) {
                IGRAPH_ERROR("No 'target' for edge in GML file.", IGRAPH_PARSEERROR);
            }
        } else {
            /* anything to do? Maybe add as graph attribute.... */
        }
    }

    /* check vertex id uniqueness */
    if (igraph_trie_size(&trie) != no_of_nodes) {
        IGRAPH_ERROR("Node 'id' not unique in GML file.", IGRAPH_PARSEERROR);
    }

    /* now we allocate the vectors and strvectors for the attributes */
    for (i = 0; i < igraph_vector_ptr_size(&vattrs); i++) {
        igraph_attribute_record_t *atrec = VECTOR(vattrs)[i];
        int type = atrec->type;
        if (type == IGRAPH_ATTRIBUTE_NUMERIC) {
            igraph_vector_t *p = IGRAPH_CALLOC(1, igraph_vector_t);
            atrec->value = p;
            IGRAPH_CHECK(igraph_vector_init(p, no_of_nodes));
        } else if (type == IGRAPH_ATTRIBUTE_STRING) {
            igraph_strvector_t *p = IGRAPH_CALLOC(1, igraph_strvector_t);
            atrec->value = p;
            IGRAPH_CHECK(igraph_strvector_init(p, no_of_nodes));
        } else {
            IGRAPH_WARNING("A composite attribute was ignored in the GML file.");
        }
    }

    for (i = 0; i < igraph_vector_ptr_size(&eattrs); i++) {
        igraph_attribute_record_t *atrec = VECTOR(eattrs)[i];
        int type = atrec->type;
        if (type == IGRAPH_ATTRIBUTE_NUMERIC) {
            igraph_vector_t *p = IGRAPH_CALLOC(1, igraph_vector_t);
            atrec->value = p;
            IGRAPH_CHECK(igraph_vector_init(p, no_of_edges));
        } else if (type == IGRAPH_ATTRIBUTE_STRING) {
            igraph_strvector_t *p = IGRAPH_CALLOC(1, igraph_strvector_t);
            atrec->value = p;
            IGRAPH_CHECK(igraph_strvector_init(p, no_of_edges));
        } else {
            IGRAPH_WARNING("A composite attribute was ignored in the GML file.");
        }
    }

    /* Ok, now the edges, attributes too */
    IGRAPH_CHECK(igraph_vector_resize(&edges, no_of_edges * 2));
    p = -1;
    while ( (p = igraph_gml_tree_find(gtree, "edge", p + 1)) != -1) {
        igraph_gml_tree_t *edge;
        long int from, to, fromidx = 0, toidx = 0;
        char name[100];
        long int j;
        edge = igraph_gml_tree_get_tree(gtree, p);
        for (j = 0; j < igraph_gml_tree_length(edge); j++) {
            const char *n = igraph_gml_tree_name(edge, j);
            if (!strcmp(n, "source")) {
                fromidx = igraph_gml_tree_find(edge, "source", 0);
            } else if (!strcmp(n, "target")) {
                toidx = igraph_gml_tree_find(edge, "target", 0);
            } else {
                long int edgeid = edgeptr / 2;
                long int trieidx;
                igraph_attribute_record_t *atrec;
                int type;
                igraph_trie_get(&eattrnames, n, &trieidx);
                atrec = VECTOR(eattrs)[trieidx];
                type = atrec->type;
                if (type == IGRAPH_ATTRIBUTE_NUMERIC) {
                    igraph_vector_t *v = (igraph_vector_t *)atrec->value;
                    IGRAPH_CHECK(igraph_i_gml_toreal(edge, j, VECTOR(*v) + edgeid));
                } else if (type == IGRAPH_ATTRIBUTE_STRING) {
                    igraph_strvector_t *v = (igraph_strvector_t *)atrec->value;
                    const char *value = igraph_i_gml_tostring(edge, j);
                    IGRAPH_CHECK(igraph_strvector_set(v, edgeid, value));
                }
            }
        }
        from = igraph_gml_tree_get_integer(edge, fromidx);
        to = igraph_gml_tree_get_integer(edge, toidx);
        snprintf(name, sizeof(name) / sizeof(char) -1, "%li", from);
        IGRAPH_CHECK(igraph_trie_get(&trie, name, &from));
        snprintf(name, sizeof(name) / sizeof(char) -1, "%li", to);
        IGRAPH_CHECK(igraph_trie_get(&trie, name, &to));
        if (igraph_trie_size(&trie) != no_of_nodes) {
            IGRAPH_ERROR("Unknown node id found in an edge in GML file.", IGRAPH_PARSEERROR);
        }
        VECTOR(edges)[edgeptr++] = from;
        VECTOR(edges)[edgeptr++] = to;
    }

    /* and add vertex attributes */
    for (i = 0; i < igraph_gml_tree_length(gtree); i++) {
        const char *n;
        char name[100];
        long int j, k;
        n = igraph_gml_tree_name(gtree, i);
        if (!strcmp(n, "node")) {
            igraph_gml_tree_t *node = igraph_gml_tree_get_tree(gtree, i);
            long int iidx = igraph_gml_tree_find(node, "id", 0);
            long int id = igraph_gml_tree_get_integer(node, iidx);
            snprintf(name, sizeof(name) / sizeof(char) -1, "%li", id);
            igraph_trie_get(&trie, name, &id);
            for (j = 0; j < igraph_gml_tree_length(node); j++) {
                const char *aname = igraph_gml_tree_name(node, j);
                igraph_attribute_record_t *atrec;
                int type;
                igraph_trie_get(&vattrnames, aname, &k);
                atrec = VECTOR(vattrs)[k];
                type = atrec->type;
                if (type == IGRAPH_ATTRIBUTE_NUMERIC) {
                    igraph_vector_t *v = (igraph_vector_t *)atrec->value;
                    IGRAPH_CHECK(igraph_i_gml_toreal(node, j, VECTOR(*v) + id));
                } else if (type == IGRAPH_ATTRIBUTE_STRING) {
                    igraph_strvector_t *v = (igraph_strvector_t *)atrec->value;
                    const char *value = igraph_i_gml_tostring(node, j);
                    IGRAPH_CHECK(igraph_strvector_set(v, id, value));
                }
            }
        }
    }

    igraph_trie_destroy(&trie);
    igraph_trie_destroy(&gattrnames);
    igraph_trie_destroy(&vattrnames);
    igraph_trie_destroy(&eattrnames);
    IGRAPH_FINALLY_CLEAN(4);

    IGRAPH_CHECK(igraph_empty_attrs(graph, 0, directed, 0)); /* TODO */
    IGRAPH_CHECK(igraph_add_vertices(graph, (igraph_integer_t) no_of_nodes,
                                     &vattrs));
    IGRAPH_CHECK(igraph_add_edges(graph, &edges, &eattrs));

    igraph_i_gml_destroy_attrs(attrs);
    igraph_vector_destroy(&edges);
    igraph_i_gml_parsedata_destroy(&context);
    IGRAPH_FINALLY_CLEAN(3);

    return IGRAPH_SUCCESS;
}

static int igraph_i_gml_convert_to_key(const char *orig, char **key) {
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
        IGRAPH_ERROR("Writing GML format failed.", IGRAPH_ENOMEM);
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
 * \brief Write the graph to a stream in GML format
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
 * If the <parameter>id</parameter> argument is not 0 then it should be a numeric
 * vector with the vertex ids and the <quote>id</quote> vertex attribute is
 * ignored (if there is one). If <parameter>id</parameter> is 0 and there is a
 * numeric <quote>id</quote> vertex attribute that is used instead. If ids
 * are not specified in either way then the regular igraph vertex ids are used.
 *
 * </para><para> Note that whichever way vertex ids are specified, their
 * uniqueness is not checked.
 *
 * </para><para> If the graph has edge attributes named <quote>source</quote>
 * or <quote>target</quote> they're silently ignored. GML uses these attributes
 * to specify the edges, so we cannot write them to the file. Rename them
 * before calling this function if you want to preserve them.
 * \param graph The graph to write to the stream.
 * \param outstream The stream to write the file to.
 * \param id Either <code>NULL</code> or a numeric vector with the vertex ids.
 *        See details above.
 * \param creator An optional string to write to the stream in the creator line.
 *        If this is 0 then the current date and time is added.
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

int igraph_write_graph_gml(const igraph_t *graph, FILE *outstream,
                           const igraph_vector_t *id, const char *creator) {
    int ret;
    igraph_strvector_t gnames, vnames, enames;
    igraph_vector_t gtypes, vtypes, etypes;
    igraph_vector_t numv;
    igraph_strvector_t strv;
    igraph_vector_bool_t boolv;
    long int i;
    long int no_of_nodes = igraph_vcount(graph);
    long int no_of_edges = igraph_ecount(graph);

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
    IGRAPH_VECTOR_INIT_FINALLY(&gtypes, 0);
    IGRAPH_VECTOR_INIT_FINALLY(&vtypes, 0);
    IGRAPH_VECTOR_INIT_FINALLY(&etypes, 0);
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
        for (i = 0; i < igraph_vector_size(&vtypes); i++) {
            char *n;
            igraph_strvector_get(&vnames, i, &n);
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

    /* directedness */
    CHECK(fprintf(outstream, "  directed %i\n", igraph_is_directed(graph) ? 1 : 0));

    /* Graph attributes first */
    for (i = 0; i < igraph_vector_size(&gtypes); i++) {
        char *name, *newname;
        igraph_strvector_get(&gnames, i, &name);
        IGRAPH_CHECK(igraph_i_gml_convert_to_key(name, &newname));
        IGRAPH_FINALLY(igraph_free, newname);
        if (VECTOR(gtypes)[i] == IGRAPH_ATTRIBUTE_NUMERIC) {
            IGRAPH_CHECK(igraph_i_attribute_get_numeric_graph_attr(graph, name, &numv));
            CHECK(fprintf(outstream, "  %s ", newname));
            CHECK(igraph_real_fprintf_precise(outstream, VECTOR(numv)[0]));
            CHECK(fputc('\n', outstream));
        } else if (VECTOR(gtypes)[i] == IGRAPH_ATTRIBUTE_STRING) {
            char *s;
            IGRAPH_CHECK(igraph_i_attribute_get_string_graph_attr(graph, name, &strv));
            igraph_strvector_get(&strv, 0, &s);
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
        long int j;
        CHECK(fprintf(outstream, "  node\n  [\n"));
        /* id */
        CHECK(fprintf(outstream, "    id %li\n", myid ? (long int)VECTOR(*myid)[i] : i));
        /* other attributes */
        for (j = 0; j < igraph_vector_size(&vtypes); j++) {
            int type = (int) VECTOR(vtypes)[j];
            char *name, *newname;
            igraph_strvector_get(&vnames, j, &name);
            if (!strcmp(name, "id")) {
                continue;
            }
            IGRAPH_CHECK(igraph_i_gml_convert_to_key(name, &newname));
            IGRAPH_FINALLY(igraph_free, newname);
            if (type == IGRAPH_ATTRIBUTE_NUMERIC) {
                IGRAPH_CHECK(igraph_i_attribute_get_numeric_vertex_attr(graph, name,
                             igraph_vss_1((igraph_integer_t) i), &numv));
                CHECK(fprintf(outstream, "    %s ", newname));
                CHECK(igraph_real_fprintf_precise(outstream, VECTOR(numv)[0]));
                CHECK(fputc('\n', outstream));
            } else if (type == IGRAPH_ATTRIBUTE_STRING) {
                char *s;
                IGRAPH_CHECK(igraph_i_attribute_get_string_vertex_attr(graph, name,
                             igraph_vss_1((igraph_integer_t) i), &strv));
                igraph_strvector_get(&strv, 0, &s);
                CHECK(fprintf(outstream, "    %s \"%s\"\n", newname, s));
            } else if (type == IGRAPH_ATTRIBUTE_BOOLEAN) {
                IGRAPH_CHECK(igraph_i_attribute_get_bool_vertex_attr(graph, name,
                             igraph_vss_1((igraph_integer_t) i), &boolv));
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
        long int from = IGRAPH_FROM(graph, i);
        long int to = IGRAPH_TO(graph, i);
        long int j;
        CHECK(fprintf(outstream, "  edge\n  [\n"));
        /* source and target */
        CHECK(fprintf(outstream, "    source %li\n",
                      myid ? (long int)VECTOR(*myid)[from] : from));
        CHECK(fprintf(outstream, "    target %li\n",
                      myid ? (long int)VECTOR(*myid)[to] : to));

        /* other attributes */
        for (j = 0; j < igraph_vector_size(&etypes); j++) {
            int type = (int) VECTOR(etypes)[j];
            char *name, *newname;
            igraph_strvector_get(&enames, j, &name);
            if (!strcmp(name, "source") || !strcmp(name, "target")) {
                continue;
            }
            IGRAPH_CHECK(igraph_i_gml_convert_to_key(name, &newname));
            IGRAPH_FINALLY(igraph_free, newname);
            if (type == IGRAPH_ATTRIBUTE_NUMERIC) {
                IGRAPH_CHECK(igraph_i_attribute_get_numeric_edge_attr(graph, name,
                             igraph_ess_1((igraph_integer_t) i), &numv));
                CHECK(fprintf(outstream, "    %s ", newname));
                CHECK(igraph_real_fprintf_precise(outstream, VECTOR(numv)[0]));
                CHECK(fputc('\n', outstream));
            } else if (type == IGRAPH_ATTRIBUTE_STRING) {
                char *s;
                IGRAPH_CHECK(igraph_i_attribute_get_string_edge_attr(graph, name,
                             igraph_ess_1((igraph_integer_t) i), &strv));
                igraph_strvector_get(&strv, 0, &s);
                CHECK(fprintf(outstream, "    %s \"%s\"\n", newname, s));
            } else if (type == IGRAPH_ATTRIBUTE_BOOLEAN) {
                IGRAPH_CHECK(igraph_i_attribute_get_bool_edge_attr(graph, name,
                             igraph_ess_1((igraph_integer_t) i), &boolv));
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
    igraph_vector_destroy(&etypes);
    igraph_vector_destroy(&vtypes);
    igraph_vector_destroy(&gtypes);
    igraph_strvector_destroy(&enames);
    igraph_strvector_destroy(&vnames);
    igraph_strvector_destroy(&gnames);
    IGRAPH_FINALLY_CLEAN(9);

    return IGRAPH_SUCCESS;
}

#undef CHECK
