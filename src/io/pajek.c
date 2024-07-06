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
#include "igraph_error.h"
#include "igraph_interface.h"
#include "igraph_memory.h"

#include "graph/attributes.h"

#include "internal/hacks.h" /* IGRAPH_STATIC_ASSERT */

#include "io/pajek-header.h"
#include "io/parsers/pajek-parser.h"

#include <ctype.h>
#include <string.h>

int igraph_pajek_yylex_init_extra(igraph_i_pajek_parsedata_t *user_defined,
                                  void *scanner);
int igraph_pajek_yylex_destroy(void *scanner);
int igraph_pajek_yyparse(igraph_i_pajek_parsedata_t *context);
void igraph_pajek_yyset_in(FILE *in_str, void *yyscanner);

/* for IGRAPH_FINALLY, which assumes that destructor functions return void */
void igraph_pajek_yylex_destroy_wrapper (void *scanner ) {
    (void) igraph_pajek_yylex_destroy(scanner);
}

void igraph_i_pajek_destroy_attr_vector(igraph_vector_ptr_t *attrs) {
    const igraph_integer_t attr_count = igraph_vector_ptr_size(attrs);
    for (igraph_integer_t i = 0; i < attr_count; i++) {
        igraph_attribute_record_t *rec = VECTOR(*attrs)[i];
        if (rec->type == IGRAPH_ATTRIBUTE_NUMERIC) {
            igraph_vector_t *vec = (igraph_vector_t*) rec->value;
            igraph_vector_destroy(vec);
            IGRAPH_FREE(vec);
        } else if (rec->type == IGRAPH_ATTRIBUTE_BOOLEAN) {
            igraph_vector_bool_t *vec = (igraph_vector_bool_t*) rec->value;
            igraph_vector_bool_destroy(vec);
            IGRAPH_FREE(vec);
        } else if (rec->type == IGRAPH_ATTRIBUTE_STRING) {
            igraph_strvector_t *strvec = (igraph_strvector_t *)rec->value;
            igraph_strvector_destroy(strvec);
            IGRAPH_FREE(strvec);
        } else {
            /* Must never reach here */
            IGRAPH_FATAL("Unknown attribute type encountered.");
        }
        IGRAPH_FREE(rec->name);
        IGRAPH_FREE(rec);
    }
    igraph_vector_ptr_destroy(attrs);
}

/**
 * \function igraph_read_graph_pajek
 * \brief Reads a file in Pajek format.
 *
 * Only a subset of the Pajek format is implemented. This is partially
 * because there is no formal specification for this format, but also because
 * <command>igraph</command> does not support some Pajek features, like
 * mixed graphs.
 *
 * </para><para>
 * Starting from version 0.6.1 igraph reads bipartite (two-mode)
 * graphs from Pajek files and adds the \c type Boolean vertex attribute for
 * them. Warnings are given for invalid edges, i.e. edges connecting
 * vertices of the same type.
 *
 * </para><para>
 * The list of the current limitations:
 * \olist
 * \oli Only <filename>.net</filename> files are supported, Pajek
 * project files (<filename>.paj</filename>) are not.
 * \oli Temporal networks (i.e. with time events) are not supported.
 * \oli Graphs with both directed and non-directed edges are not
 * supported, as they cannot be represented in <command>igraph</command>.
 * \oli Only Pajek networks are supported; permutations, hierarchies,
 * clusters and vectors are not.
 * \oli Multi-relational networks (i.e. networks with multiple edge
 * types) are not supported.
 * \oli Unicode characters encoded as <code>&amp;#dddd;</code>, or newlines
 * encoded as <code>\n</code> will not be decoded.
 * \endolist
 *
 * </para><para>
 * If an attribute handler is installed,
 * <command>igraph</command> also reads the vertex and edge attributes
 * from the file. Most attributes are renamed to be more informative:
 * \c color instead of \c c, \c xfact instead of \c x_fact,
 * \c yfact instead of y_fact, \c labeldist instead of \c lr,
 * \c labeldegree2 instead of \c lphi, \c framewidth instead of \c bw,
 * \c fontsize instead of \c fos, \c rotation instead of \c phi,
 * \c radius instead of \c r, \c diamondratio instead of \c q,
 * \c labeldegree instead of \c la,
 * \c color instead of \c ic, \c framecolor instead of \c bc,
 * \c labelcolor instead of \c lc; these belong to vertices.
 *
 * </para><para>
 * Edge attributes are also renamed, \c s to \c arrowsize,
 * \c w to \c edgewidth, \c h1 to \c hook1, \c h2 to \c hook2,
 * \c a1 to \c angle1, \c a2 to \c angle2, \c k1 to
 * \c velocity1, \c k2 to \c velocity2, \c ap to \c arrowpos,
 * \c lp to \c labelpos, \c lr to \c labelangle,
 * \c lphi to \c labelangle2, \c la to \c labeldegree,
 * \c fos to \c fontsize, \c a to \c arrowtype, \c p to \c linepattern,
 * \c l to \c label, \c lc to \c labelcolor, \c c to \c color.
 *
 * </para><para>
 * Unknown vertex or edge parameters are read as string vertex
 * or edge attributes. If the parameter name conflicts with one
 * the standard attribute names mentioned above, a <code>_</code>
 * character is appended to it to avoid conflict.
 *
 * </para><para>
 * In addition the following vertex attributes might be added: \c id
 * and \c name are added (with the same value) if there are vertex IDs in the
 * file. \c id is deprecated in favour of \c name and will no longer be used
 * by future versions of igraph. \c x and \c y, and potentially \c z are also
 * added if there are vertex coordinates in the file.
 *
 * </para><para>
 * The \c weight edge attribute will be added if there are edge weights present.
 *
 * </para><para>
 * See the Pajek homepage:
 * http://vlado.fmf.uni-lj.si/pub/networks/pajek/ for more info on
 * Pajek. The Pajek manual,
 * http://vlado.fmf.uni-lj.si/pub/networks/pajek/doc/pajekman.pdf,
 * and http://mrvar.fdv.uni-lj.si/pajek/DrawEPS.htm
 * have information on the Pajek file format. There is additional
 * useful information and sample files at
 * http://mrvar.fdv.uni-lj.si/pajek/history.htm
 *
 * \param graph Pointer to an uninitialized graph object.
 * \param file An already opened file handler.
 * \return Error code.
 *
 * Time complexity: O(|V|+|E|+|A|), |V| is the number of vertices, |E|
 * the number of edges, |A| the number of attributes (vertex + edge)
 * in the graph if there are attribute handlers installed.
 *
 * \sa \ref igraph_write_graph_pajek() for writing Pajek files, \ref
 * igraph_read_graph_graphml() for reading GraphML files.
 *
 * \example examples/simple/foreign.c
 */

igraph_error_t igraph_read_graph_pajek(igraph_t *graph, FILE *instream) {

    igraph_vector_int_t edges;
    igraph_trie_t vattrnames;
    igraph_vector_ptr_t vattrs;
    igraph_trie_t eattrnames;
    igraph_vector_ptr_t eattrs;
    igraph_integer_t i, j;
    igraph_i_pajek_parsedata_t context;

    IGRAPH_VECTOR_INT_INIT_FINALLY(&edges, 0);

    IGRAPH_TRIE_INIT_FINALLY(&vattrnames, 1);
    IGRAPH_CHECK(igraph_vector_ptr_init(&vattrs, 0));
    IGRAPH_FINALLY(igraph_i_pajek_destroy_attr_vector, &vattrs);

    IGRAPH_TRIE_INIT_FINALLY(&eattrnames, 1);
    IGRAPH_CHECK(igraph_vector_ptr_init(&eattrs, 0));
    IGRAPH_FINALLY(igraph_i_pajek_destroy_attr_vector, &eattrs);

    context.directed = false; /* assume undirected until an element implying directedness is encountered */
    context.vector = &edges;
    context.vcount = -1;
    context.vertexid = 0;
    context.vertex_attribute_names = &vattrnames;
    context.vertex_attributes = &vattrs;
    context.edge_attribute_names = &eattrnames;
    context.edge_attributes = &eattrs;
    context.actedge = 0;
    context.eof = false;
    context.errmsg[0] = '\0';
    context.igraph_errno = IGRAPH_SUCCESS;

    igraph_pajek_yylex_init_extra(&context, &context.scanner);
    IGRAPH_FINALLY(igraph_pajek_yylex_destroy_wrapper, context.scanner);

    igraph_pajek_yyset_in(instream, context.scanner);

    /* Use ENTER/EXIT to avoid destroying context.scanner before this function returns */
    IGRAPH_FINALLY_ENTER();
    int err = igraph_pajek_yyparse(&context);
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
            IGRAPH_ERROR("Cannot read Pajek file.", IGRAPH_PARSEERROR);
        }
        break;
    case 2: /* out of memory */
        IGRAPH_ERROR("Cannot read Pajek file.", IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */
        break;
    default: /* must never reach here */
        /* Hint: This will usually be triggered if an IGRAPH_CHECK() is used in a Bison
         * action instead of an IGRAPH_YY_CHECK(), resulting in an igraph errno being
         * returned in place of a Bison error code.
         * TODO: What if future Bison versions introduce error codes other than 0, 1 and 2?
         */
        IGRAPH_FATALF("Parser returned unexpected error code (%d) when reading Pajek file.", err);
    }

    /* Prepare attributes */
    const igraph_integer_t eattr_count = igraph_vector_ptr_size(&eattrs);
    for (i = 0; i < eattr_count; i++) {
        igraph_attribute_record_t *rec = VECTOR(eattrs)[i];
        if (rec->type == IGRAPH_ATTRIBUTE_NUMERIC) {
            igraph_vector_t *vec = (igraph_vector_t*)rec->value;
            igraph_integer_t origsize = igraph_vector_size(vec);
            IGRAPH_CHECK(igraph_vector_resize(vec, context.actedge));
            for (j = origsize; j < context.actedge; j++) {
                VECTOR(*vec)[j] = IGRAPH_NAN;
            }
        } else if (rec->type == IGRAPH_ATTRIBUTE_BOOLEAN) {
            /* Boolean attributes are not currently added by the parser.
             * This section is here for future-proofing. */
            igraph_vector_bool_t *vec = (igraph_vector_bool_t*)rec->value;
            igraph_integer_t origsize = igraph_vector_bool_size(vec);
            IGRAPH_CHECK(igraph_vector_bool_resize(vec, context.actedge));
            for (j = origsize; j < context.actedge; j++) {
                VECTOR(*vec)[j] = 0;
            }
        } else if (rec->type == IGRAPH_ATTRIBUTE_STRING) {
            igraph_strvector_t *strvec = (igraph_strvector_t*)rec->value;
            /* strvector_resize() adds empty strings */
            IGRAPH_CHECK(igraph_strvector_resize(strvec, context.actedge));
        } else {
            /* Must never reach here */
            IGRAPH_FATAL("Unknown attribute type encountered.");
        }
    }

    /* Create graph */
    IGRAPH_CHECK(igraph_empty(graph, 0, context.directed));
    IGRAPH_FINALLY(igraph_destroy, graph);
    IGRAPH_CHECK(igraph_add_vertices(graph, context.vcount, &vattrs));
    IGRAPH_CHECK(igraph_add_edges(graph, &edges, &eattrs));

    igraph_vector_int_destroy(&edges);
    igraph_i_pajek_destroy_attr_vector(&eattrs);
    igraph_trie_destroy(&eattrnames);
    igraph_i_pajek_destroy_attr_vector(&vattrs);
    igraph_trie_destroy(&vattrnames);
    igraph_pajek_yylex_destroy(context.scanner);
    IGRAPH_FINALLY_CLEAN(7); /* +1 for 'graph' */

    return IGRAPH_SUCCESS;
}

/***** Writing Pajek files *****/

/* Order matters here! */
#define V_ID                0
#define V_X                 1
#define V_Y                 2
#define V_Z                 3
#define V_SHAPE             4
#define V_XFACT             5
#define V_YFACT             6
#define V_LABELDIST         7
#define V_LABELDEGREE2      8
#define V_FRAMEWIDTH        9
#define V_FONTSIZE         10
#define V_ROTATION         11
#define V_RADIUS           12
#define V_DIAMONDRATIO     13
#define V_LABELDEGREE      14
#define V_FONT             15
#define V_URL              16
#define V_COLOR            17
#define V_FRAMECOLOR       18
#define V_LABELCOLOR       19
#define V_LAST             20

#define E_WEIGHT            0
#define E_LAST              1

/* Pajek encodes newlines as \n, and any unicode character can be encoded
 * in the form &#hhhh;. Therefore we encode quotation marks as &#34; */
static igraph_error_t igraph_i_pajek_escape(const char* src, char** dest) {
    igraph_integer_t destlen = 0;
    igraph_bool_t need_escape = false;

    /* Determine whether the string contains characters to be escaped */
    const char *s;
    char *d;
    for (s = src; *s; s++, destlen++) {
        if (*s == '\n' || *s == '\r') {
            need_escape = true;
            destlen++;
        } else if (*s == '"') {
            need_escape = true;
            destlen += 4;
        } else if (!isalnum(*s)) {
            need_escape = true;
        }
    }

    if (!need_escape) {
        /* At this point, we know that the string does not contain any chars
         * that would warrant encoding. Therefore, we simply quote it and
         * return the quoted string. This is necessary because Pajek uses some
         * reserved words in its format (like 'c' standing for color) and they
         * have to be quoted as well.
         */
        *dest = IGRAPH_CALLOC(destlen + 3, char);
        CHECK_OOM_WP(*dest);

        d = *dest;
        strcpy(d + 1, src);
        d[0] = d[destlen + 1] = '"';
        d[destlen + 2] = 0;
        return IGRAPH_SUCCESS;
    }

    *dest = IGRAPH_CALLOC(destlen + 3, char);
    CHECK_OOM_WP(*dest);

    d = *dest;
    *d = '"'; d++;

    for (s = src; *s; s++, d++) {
        switch (*s) {
        /* Encode quotation marks as &#34;, as they would otherwise signify
           the end/beginning of a string. */
        case '"':
            strcpy(d, "&#34;"); d += 4; break;
            break;
        /* Encode both CR and LF as \n, as neither should apear in a quoted string.
           \n is the _only_ escape sequence Pajek understands. */
        case '\n':
        case '\r':
            *d = '\\'; d++;
            *d = 'n';
            break;
        default:
            *d = *s;
        }
    }
    *d = '"'; d++; *d = 0;

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_write_graph_pajek
 * \brief Writes a graph to a file in Pajek format.
 *
 * Writes files in the native format of the Pajek software. This format
 * is not recommended for data exchange or archival. It is meant solely
 * for interoperability with Pajek.
 *
 * </para><para>
 * The Pajek vertex and edge parameters (like color) are determined by
 * the attributes of the vertices and edges. Of course this requires
 * an attribute handler to be installed. The names of the
 * corresponding vertex and edge attributes are listed at \ref
 * igraph_read_graph_pajek(), e.g. the \c color vertex attributes
 * determines the color (\c c in Pajek) parameter.
 *
 * </para><para>
 * Vertex and edge attributes that do not correspond to any documented
 * Pajek parameter are discarded.
 *
 * </para><para>
 * As of version 0.6.1 igraph writes bipartite graphs into Pajek files
 * correctly, i.e. they will be also bipartite when read into Pajek.
 * As Pajek is less flexible for bipartite graphs (the numeric IDs of
 * the vertices must be sorted according to vertex type), igraph might
 * need to reorder the vertices when writing a bipartite Pajek file.
 * This effectively means that numeric vertex IDs usually change when
 * a bipartite graph is written to a Pajek file, and then read back
 * into igraph.
 *
 * </para><para>
 * Early versions of Pajek supported only Windows-style line endings
 * in Pajek files, but recent versions support both Windows and Unix
 * line endings. igraph therefore uses the platform-native line endings
 * when the input file is opened in text mode, and uses Unix-style
 * line endings when the input file is opened in binary mode. If you
 * are using an old version of Pajek, you are on Unix and you are having
 * problems reading files written by igraph on a Windows machine, convert the
 * line endings manually with a text editor or with \c unix2dos or \c iconv
 * from the command line).
 *
 * </para><para>
 * Pajek will only interpret UTF-8 encoded files if they contain a byte-order
 * mark (BOM) at the beginning. igraph is agnostic of string attribute encodings
 * and therefore it will never write a BOM. You need to add this manually
 * if/when necessary.
 *
 * \param graph The graph object to write.
 * \param outstream The file to write to. It should be opened and writable.
 * \return Error code.
 *
 * Time complexity: O(|V|+|E|+|A|), |V| is the number of vertices, |E|
 * is the number of edges, |A| the number of attributes (vertex +
 * edge) in the graph if there are attribute handlers installed.
 *
 * \sa \ref igraph_read_graph_pajek() for reading Pajek graphs, \ref
 * igraph_write_graph_graphml() for writing a graph in GraphML format,
 * this suites <command>igraph</command> graphs better.
 *
 * \example examples/simple/igraph_write_graph_pajek.c
 */

igraph_error_t igraph_write_graph_pajek(const igraph_t *graph, FILE *outstream) {
    igraph_integer_t no_of_nodes = igraph_vcount(graph);

    igraph_attribute_type_t vtypes[V_LAST], etypes[E_LAST];
    igraph_bool_t write_vertex_attrs = false;

    /* Same order as the #define's */
    const char *vnames[] = { "id", "x", "y", "z", "shape", "xfact", "yfact",
        "labeldist", "labeldegree2", "framewidth",
        "fontsize", "rotation", "radius",
        "diamondratio", "labeldegree",
        "font", "url", "color", "framecolor",
        "labelcolor"
    };
    IGRAPH_STATIC_ASSERT(sizeof(vnames) / sizeof(vnames[0]) == V_LAST);

    /* Arrays called xxx[]  are igraph attribute names,
     *               xxx2[] are the corresponding Pajek names. */
    const char *vnumnames[] = { "xfact", "yfact", "labeldist",
                                "labeldegree2", "framewidth", "fontsize",
                                "rotation", "radius", "diamondratio",
                                "labeldegree"
                              };
    const char *vnumnames2[] = { "x_fact", "y_fact", "lr", "lphi", "bw",
                                 "fos", "phi", "r", "q", "la"
                               };
    IGRAPH_STATIC_ASSERT(sizeof(vnumnames) == sizeof(vnumnames2));

    const char *vstrnames[] = { "font", "url", "color", "framecolor",
                                "labelcolor"
                              };
    const char *vstrnames2[] = { "font", "url", "ic", "bc", "lc" };
    IGRAPH_STATIC_ASSERT(sizeof(vstrnames) == sizeof(vstrnames2));

    /* Same order as the #define's */
    const char *enames[] = { "weight" };
    IGRAPH_STATIC_ASSERT(sizeof(enames) / sizeof(enames[0]) == E_LAST);

    const char *enumnames[] = { "arrowsize", "edgewidth", "hook1", "hook2",
                                "angle1", "angle2", "velocity1", "velocity2",
                                "arrowpos", "labelpos", "labelangle",
                                "labelangle2", "labeldegree", "fontsize"
                              };
    const char *enumnames2[] = { "s", "w", "h1", "h2", "a1", "a2", "k1", "k2",
                                 "ap", "lp", "lr", "lphi", "la", "fos"
                               };
    IGRAPH_STATIC_ASSERT(sizeof(enumnames) == sizeof(enumnames2));

    const char *estrnames[] = { "arrowtype", "linepattern", "label",
                                "labelcolor", "color", "font"
                              };
    const char *estrnames2[] = { "a", "p", "l", "lc", "c", "font" };
    IGRAPH_STATIC_ASSERT(sizeof(estrnames) == sizeof(estrnames2));

    /* Newer Pajek versions support both Unix and Windows-style line endings,
     * so we just use Unix style. This will get converted to CRLF on Windows
     * when the file is opened in text mode */
    const char *newline = "\n";

    igraph_es_t es;
    igraph_eit_t eit;

    igraph_vector_t numv;
    igraph_strvector_t strv;

    igraph_vector_int_t ex_numa;
    igraph_vector_int_t ex_stra;
    igraph_vector_int_t vx_numa;
    igraph_vector_int_t vx_stra;

    const char *s;
    char *escaped;

    igraph_bool_t bipartite = false;
    igraph_vector_int_t bip_index, bip_index2;
    igraph_vector_bool_t bvec;
    igraph_integer_t notop = 0, nobottom = 0;

    IGRAPH_VECTOR_INIT_FINALLY(&numv, 1);
    IGRAPH_STRVECTOR_INIT_FINALLY(&strv, 1);

    IGRAPH_VECTOR_INT_INIT_FINALLY(&ex_numa, 0);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&ex_stra, 0);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&vx_numa, 0);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&vx_stra, 0);

    /* Check if graph is bipartite, i.e. whether it has a Boolean 'type' vertex attribute. */
    if (igraph_i_attribute_has_attr(graph, IGRAPH_ATTRIBUTE_VERTEX, "type")) {
        igraph_attribute_type_t type_type;
        IGRAPH_CHECK(igraph_i_attribute_gettype(graph, &type_type, IGRAPH_ATTRIBUTE_VERTEX, "type"));
        if (type_type == IGRAPH_ATTRIBUTE_BOOLEAN) {
            bipartite = true; write_vertex_attrs = true;
            /* Count top and bottom vertices, we go over them twice,
            because we want to keep their original order */
            IGRAPH_VECTOR_INT_INIT_FINALLY(&bip_index, no_of_nodes);
            IGRAPH_VECTOR_INT_INIT_FINALLY(&bip_index2, no_of_nodes);
            IGRAPH_VECTOR_BOOL_INIT_FINALLY(&bvec, 1);
            for (igraph_integer_t i = 0; i < no_of_nodes; i++) {
                IGRAPH_CHECK(igraph_i_attribute_get_bool_vertex_attr(graph,
                             "type", igraph_vss_1(i), &bvec));
                if (VECTOR(bvec)[0]) {
                    notop++;
                } else {
                    nobottom++;
                }
            }
            for (igraph_integer_t i = 0, bptr = 0, tptr = nobottom; i < no_of_nodes; i++) {
                IGRAPH_CHECK(igraph_i_attribute_get_bool_vertex_attr(graph,
                             "type", igraph_vss_1(i), &bvec));
                if (VECTOR(bvec)[0]) {
                    VECTOR(bip_index)[tptr] = i;
                    VECTOR(bip_index2)[i] = tptr;
                    tptr++;
                } else {
                    VECTOR(bip_index)[bptr] = i;
                    VECTOR(bip_index2)[i] = bptr;
                    bptr++;
                }
            }
            igraph_vector_bool_destroy(&bvec);
            IGRAPH_FINALLY_CLEAN(1);
        }
    }

    /* Write header */
    if (bipartite) {
        if (fprintf(outstream, "*Vertices %" IGRAPH_PRId " %" IGRAPH_PRId "%s", no_of_nodes, nobottom,
                    newline) < 0) {
            IGRAPH_ERROR("Cannot write pajek file.", IGRAPH_EFILE);
        }
    } else {
        if (fprintf(outstream, "*Vertices %" IGRAPH_PRId "%s", no_of_nodes, newline) < 0) {
            IGRAPH_ERROR("Cannot write pajek file.", IGRAPH_EFILE);
        }
    }

    /* Check the vertex attributes, and determine if we need to write them. */
    memset(vtypes, 0, sizeof(vtypes[0])*V_LAST);
    for (igraph_integer_t i = 0; i < V_LAST; i++) {
        if (igraph_i_attribute_has_attr(graph, IGRAPH_ATTRIBUTE_VERTEX, vnames[i])) {
            IGRAPH_CHECK(igraph_i_attribute_gettype(
                             graph, &vtypes[i], IGRAPH_ATTRIBUTE_VERTEX, vnames[i]));
            write_vertex_attrs = true;
        } else {
            vtypes[i] = (igraph_attribute_type_t) -1;
        }
    }
    for (igraph_integer_t i = 0; i < (igraph_integer_t) (sizeof(vnumnames) / sizeof(vnumnames[0])); i++) {
        igraph_attribute_type_t type;
        if (igraph_i_attribute_has_attr(graph, IGRAPH_ATTRIBUTE_VERTEX, vnumnames[i])) {
            IGRAPH_CHECK(igraph_i_attribute_gettype(
                             graph, &type, IGRAPH_ATTRIBUTE_VERTEX, vnumnames[i]));
            if (type == IGRAPH_ATTRIBUTE_NUMERIC) {
                IGRAPH_CHECK(igraph_vector_int_push_back(&vx_numa, i));
            }
        }
    }
    for (igraph_integer_t i = 0; i < (igraph_integer_t) (sizeof(vstrnames) / sizeof(vstrnames[0])); i++) {
        igraph_attribute_type_t type;
        if (igraph_i_attribute_has_attr(graph, IGRAPH_ATTRIBUTE_VERTEX, vstrnames[i])) {
            IGRAPH_CHECK(igraph_i_attribute_gettype(
                             graph, &type, IGRAPH_ATTRIBUTE_VERTEX, vstrnames[i]));
            if (type == IGRAPH_ATTRIBUTE_STRING) {
                IGRAPH_CHECK(igraph_vector_int_push_back(&vx_stra, i));
            }
        }
    }

    /* Write vertices */
    if (write_vertex_attrs) {
        for (igraph_integer_t i = 0; i < no_of_nodes; i++) {
            igraph_integer_t id = bipartite ? VECTOR(bip_index)[i] : i;

            /* vertex ID */
            fprintf(outstream, "%" IGRAPH_PRId, i + 1);
            if (vtypes[V_ID] == IGRAPH_ATTRIBUTE_NUMERIC) {
                IGRAPH_CHECK(igraph_i_attribute_get_numeric_vertex_attr(
                                 graph, vnames[V_ID], igraph_vss_1(id), &numv));
                fputs(" \"", outstream);
                igraph_real_fprintf_precise(outstream, VECTOR(numv)[0]);
                fputc('"', outstream);
            } else if (vtypes[V_ID] == IGRAPH_ATTRIBUTE_STRING) {
                IGRAPH_CHECK(igraph_i_attribute_get_string_vertex_attr(
                                 graph, vnames[V_ID], igraph_vss_1(id), &strv));
                s = igraph_strvector_get(&strv, 0);
                IGRAPH_CHECK(igraph_i_pajek_escape(s, &escaped));
                fprintf(outstream, " %s", escaped);
                IGRAPH_FREE(escaped);
            } else {
                fprintf(outstream, " \"%" IGRAPH_PRId "\"", id + 1);
            }

            /* coordinates */
            if (vtypes[V_X] == IGRAPH_ATTRIBUTE_NUMERIC &&
                vtypes[V_Y] == IGRAPH_ATTRIBUTE_NUMERIC) {
                IGRAPH_CHECK(igraph_i_attribute_get_numeric_vertex_attr(
                                 graph, vnames[V_X], igraph_vss_1(id), &numv));
                fputc(' ', outstream);
                igraph_real_fprintf_precise(outstream, VECTOR(numv)[0]);
                IGRAPH_CHECK(igraph_i_attribute_get_numeric_vertex_attr(
                                 graph, vnames[V_Y], igraph_vss_1(id), &numv));
                fputc(' ', outstream);
                igraph_real_fprintf_precise(outstream, VECTOR(numv)[0]);
                if (vtypes[V_Z] == IGRAPH_ATTRIBUTE_NUMERIC) {
                    IGRAPH_CHECK(igraph_i_attribute_get_numeric_vertex_attr(graph, vnames[V_Z],
                            igraph_vss_1(id), &numv));
                    fputc(' ', outstream);
                    igraph_real_fprintf_precise(outstream, VECTOR(numv)[0]);
                }
            }

            /* shape */
            if (vtypes[V_SHAPE] == IGRAPH_ATTRIBUTE_STRING) {
                IGRAPH_CHECK(igraph_i_attribute_get_string_vertex_attr(
                                 graph, vnames[V_SHAPE], igraph_vss_1(id), &strv));
                s = igraph_strvector_get(&strv, 0);
                IGRAPH_CHECK(igraph_i_pajek_escape(s, &escaped));
                fprintf(outstream, " %s", escaped);
                IGRAPH_FREE(escaped);
            }

            /* numeric parameters */
            for (igraph_integer_t j = 0; j < igraph_vector_int_size(&vx_numa); j++) {
                igraph_integer_t idx = VECTOR(vx_numa)[j];
                IGRAPH_CHECK(igraph_i_attribute_get_numeric_vertex_attr(
                                 graph, vnumnames[idx], igraph_vss_1(id), &numv));
                fprintf(outstream, " %s ", vnumnames2[idx]);
                igraph_real_fprintf_precise(outstream, VECTOR(numv)[0]);
            }

            /* string parameters */
            for (igraph_integer_t j = 0; j < igraph_vector_int_size(&vx_stra); j++) {
                igraph_integer_t idx = VECTOR(vx_stra)[j];
                IGRAPH_CHECK(igraph_i_attribute_get_string_vertex_attr(
                                 graph, vstrnames[idx], igraph_vss_1(id), &strv));
                s = igraph_strvector_get(&strv, 0);
                IGRAPH_CHECK(igraph_i_pajek_escape(s, &escaped));
                fprintf(outstream, " %s %s", vstrnames2[idx], escaped);
                IGRAPH_FREE(escaped);
            }

            /* trailing newline */
            fprintf(outstream, "%s", newline);
        }
    }

    /* edges header */
    if (igraph_is_directed(graph)) {
        fprintf(outstream, "*Arcs%s", newline);
    } else {
        fprintf(outstream, "*Edges%s", newline);
    }

    IGRAPH_CHECK(igraph_es_all(&es, IGRAPH_EDGEORDER_ID));
    IGRAPH_FINALLY(igraph_es_destroy, &es);
    IGRAPH_CHECK(igraph_eit_create(graph, es, &eit));
    IGRAPH_FINALLY(igraph_eit_destroy, &eit);

    /* Check edge attributes */
    /* TODO: refactor and simplify since only "weight" is relevant */
    for (igraph_integer_t i = 0; i < E_LAST; i++) {
        if (igraph_i_attribute_has_attr(graph, IGRAPH_ATTRIBUTE_EDGE, enames[i])) {
            IGRAPH_CHECK(igraph_i_attribute_gettype(
                             graph, &etypes[i], IGRAPH_ATTRIBUTE_EDGE, enames[i]));
        } else {
            etypes[i] = (igraph_attribute_type_t) -1;
        }
    }
    for (igraph_integer_t i = 0; i < (igraph_integer_t) (sizeof(enumnames) / sizeof(enumnames[0])); i++) {
        igraph_attribute_type_t type;
        if (igraph_i_attribute_has_attr(graph, IGRAPH_ATTRIBUTE_EDGE, enumnames[i])) {
            IGRAPH_CHECK(igraph_i_attribute_gettype(
                             graph, &type, IGRAPH_ATTRIBUTE_EDGE, enumnames[i]));
            if (type == IGRAPH_ATTRIBUTE_NUMERIC) {
                IGRAPH_CHECK(igraph_vector_int_push_back(&ex_numa, i));
            }
        }
    }
    for (igraph_integer_t i = 0; i < (igraph_integer_t) (sizeof(estrnames) / sizeof(estrnames[0])); i++) {
        igraph_attribute_type_t type;
        if (igraph_i_attribute_has_attr(graph, IGRAPH_ATTRIBUTE_EDGE, estrnames[i])) {
            IGRAPH_CHECK(igraph_i_attribute_gettype(
                             graph, &type, IGRAPH_ATTRIBUTE_EDGE, estrnames[i]));
            if (type == IGRAPH_ATTRIBUTE_STRING) {
                IGRAPH_CHECK(igraph_vector_int_push_back(&ex_stra, i));
            }
        }
    }

    for (igraph_integer_t i = 0; !IGRAPH_EIT_END(eit); IGRAPH_EIT_NEXT(eit), i++) {
        igraph_integer_t edge = IGRAPH_EIT_GET(eit);
        igraph_integer_t from, to;
        igraph_edge(graph, edge, &from,  &to);
        if (bipartite) {
            from = VECTOR(bip_index2)[from];
            to  = VECTOR(bip_index2)[to];
        }
        fprintf(outstream, "%" IGRAPH_PRId " %" IGRAPH_PRId , from + 1, to + 1);

        /* Weights */
        if (etypes[E_WEIGHT] == IGRAPH_ATTRIBUTE_NUMERIC) {
            IGRAPH_CHECK(igraph_i_attribute_get_numeric_edge_attr(
                             graph, enames[E_WEIGHT], igraph_ess_1(edge), &numv));
            fputc(' ', outstream);
            igraph_real_fprintf_precise(outstream, VECTOR(numv)[0]);
        }

        /* numeric parameters */
        for (igraph_integer_t j = 0; j < igraph_vector_int_size(&ex_numa); j++) {
            igraph_integer_t idx = VECTOR(ex_numa)[j];
            IGRAPH_CHECK(igraph_i_attribute_get_numeric_edge_attr(
                             graph, enumnames[idx], igraph_ess_1(edge), &numv));
            fprintf(outstream, " %s ", enumnames2[idx]);
            igraph_real_fprintf_precise(outstream, VECTOR(numv)[0]);
        }

        /* string parameters */
        for (igraph_integer_t j = 0; j < igraph_vector_int_size(&ex_stra); j++) {
            igraph_integer_t idx = VECTOR(ex_stra)[j];
            IGRAPH_CHECK(igraph_i_attribute_get_string_edge_attr(
                             graph, estrnames[idx], igraph_ess_1(edge), &strv));
            s = igraph_strvector_get(&strv, 0);
            IGRAPH_CHECK(igraph_i_pajek_escape(s, &escaped));
            fprintf(outstream, " %s %s", estrnames2[idx], escaped);
            IGRAPH_FREE(escaped);
        }

        /* trailing newline */
        fprintf(outstream, "%s", newline);
    }

    igraph_eit_destroy(&eit);
    igraph_es_destroy(&es);
    IGRAPH_FINALLY_CLEAN(2);

    if (bipartite) {
        igraph_vector_int_destroy(&bip_index2);
        igraph_vector_int_destroy(&bip_index);
        IGRAPH_FINALLY_CLEAN(2);
    }

    igraph_vector_int_destroy(&ex_numa);
    igraph_vector_int_destroy(&ex_stra);
    igraph_vector_int_destroy(&vx_numa);
    igraph_vector_int_destroy(&vx_stra);
    igraph_strvector_destroy(&strv);
    igraph_vector_destroy(&numv);
    IGRAPH_FINALLY_CLEAN(6);
    return IGRAPH_SUCCESS;
}
