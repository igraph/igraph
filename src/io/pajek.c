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

#include "pajek-header.h"

#include <ctype.h>
#include <string.h>

int igraph_pajek_yylex_init_extra(igraph_i_pajek_parsedata_t* user_defined,
                                  void* scanner);
void igraph_pajek_yylex_destroy (void *scanner );
int igraph_pajek_yyparse (igraph_i_pajek_parsedata_t* context);
void igraph_pajek_yyset_in  (FILE * in_str, void* yyscanner );

/**
 * \function igraph_read_graph_pajek
 * \brief Reads a file in Pajek format
 *
 * \param graph Pointer to an uninitialized graph object.
 * \param file An already opened file handler.
 * \return Error code.
 *
 * </para><para>
 * Only a subset of the Pajek format is implemented. This is partially
 * because this format is not very well documented, but also because
 * <command>igraph</command> does not support some Pajek features, like
 * multigraphs.
 *
 * </para><para>
 * Starting from version 0.6.1 igraph reads bipartite (two-mode)
 * graphs from Pajek files and add the \c type vertex attribute for them.
 * Warnings are given for invalid edges, i.e. edges connecting
 * vertices of the same type.
 *
 * </para><para>
 * The list of the current limitations:
 * \olist
 * \oli Only <filename>.net</filename> files are supported, Pajek
 * project files (<filename>.paj</filename>) are not. These might be
 * supported in the future if there is need for it.
 * \oli Time events networks are not supported.
 * \oli Hypergraphs (i.e. graphs with non-binary edges) are not
 * supported.
 * \oli Graphs with both directed and non-directed edges are not
 * supported, are they cannot be represented in
 * <command>igraph</command>.
 * \oli Only Pajek networks are supported, permutations, hierarchies,
 * clusters and vectors are not.
 * \oli Graphs with multiple edge sets are not supported.
 * \endolist
 *
 * </para><para>
 * If there are attribute handlers installed,
 * <command>igraph</command> also reads the vertex and edge attributes
 * from the file. Most attributes are renamed to be more informative:
 * \c color instead of \c c, \c xfact instead of \c x_fact,
 * \c yfact instead of y_fact, \c labeldist instead of \c lr,
 * \c labeldegree2 instead of \c lphi, \c framewidth instead of \c bw,
 * \c fontsize
 * instead of \c fos, \c rotation instead of \c phi, \c radius instead
 * of \c r,
 * \c diamondratio instead of \c q, \c labeldegree instead of \c la,
 * \c vertexsize
 * instead of \c size, \c color instead of \c ic, \c framecolor instead of
 * \c bc, \c labelcolor instead of \c lc, these belong to vertices.
 *
 * </para><para>
 * Edge attributes are also renamed, \c s to \c arrowsize, \c w
 * to \c edgewidth, \c h1 to \c hook1, \c h2 to \c hook2,
 * \c a1 to \c angle1, \c a2 to \c angle2, \c k1 to
 * \c velocity1, \c k2 to \c velocity2, \c ap to \c
 * arrowpos, \c lp to \c labelpos, \c lr to
 * \c labelangle, \c lphi to \c labelangle2, \c la to \c
 * labeldegree, \c fos to
 * \c fontsize, \c a to \c arrowtype, \c p to \c
 * linepattern, \c l to \c label, \c lc to
 * \c labelcolor, \c c to \c color.
 *
 * </para><para>
 * In addition the following vertex attributes might be added: \c id
 * if there are vertex ids in the file, \c x and \c y or \c x
 * and \c y and \c z if there are vertex coordinates in the file.
 *
 * </para><para>The \c weight edge attribute might be
 * added if there are edge weights present.
 *
 * </para><para>
 * See the pajek homepage:
 * http://vlado.fmf.uni-lj.si/pub/networks/pajek/ for more info on
 * Pajek and the Pajek manual:
 * http://vlado.fmf.uni-lj.si/pub/networks/pajek/doc/pajekman.pdf for
 * information on the Pajek file format.
 *
 * </para><para>
 * Time complexity: O(|V|+|E|+|A|), |V| is the number of vertices, |E|
 * the number of edges, |A| the number of attributes (vertex + edge)
 * in the graph if there are attribute handlers installed.
 *
 * \sa \ref igraph_write_graph_pajek() for writing Pajek files, \ref
 * igraph_read_graph_graphml() for reading GraphML files.
 *
 * \example examples/simple/foreign.c
 */

int igraph_read_graph_pajek(igraph_t *graph, FILE *instream) {

    igraph_vector_t edges;
    igraph_trie_t vattrnames;
    igraph_vector_ptr_t vattrs;
    igraph_trie_t eattrnames;
    igraph_vector_ptr_t eattrs;
    long int i, j;
    igraph_i_pajek_parsedata_t context;

    IGRAPH_VECTOR_INIT_FINALLY(&edges, 0);

    IGRAPH_TRIE_INIT_FINALLY(&vattrnames, 1);
    IGRAPH_VECTOR_PTR_INIT_FINALLY(&vattrs, 0);
    IGRAPH_TRIE_INIT_FINALLY(&eattrnames, 1);
    IGRAPH_VECTOR_PTR_INIT_FINALLY(&eattrs, 0);

    context.vector = &edges;
    context.mode = 0;
    context.vcount = -1;
    context.vertexid = 0;
    context.vertex_attribute_names = &vattrnames;
    context.vertex_attributes = &vattrs;
    context.edge_attribute_names = &eattrnames;
    context.edge_attributes = &eattrs;
    context.actedge = 0;
    context.eof = 0;

    igraph_pajek_yylex_init_extra(&context, &context.scanner);
    IGRAPH_FINALLY(igraph_pajek_yylex_destroy, context.scanner);

    igraph_pajek_yyset_in(instream, context.scanner);

    if (igraph_pajek_yyparse(&context)) {
        if (context.errmsg[0] != 0) {
            IGRAPH_ERROR(context.errmsg, IGRAPH_PARSEERROR);
        } else {
            IGRAPH_ERROR("Cannot read Pajek file", IGRAPH_PARSEERROR);
        }
    }

    if (context.vcount < 0) {
        IGRAPH_ERROR("invalid vertex count in Pajek file", IGRAPH_EINVAL);
    }
    if (context.vcount2 < 0) {
        IGRAPH_ERROR("invalid 2-mode vertex count in Pajek file", IGRAPH_EINVAL);
    }

    for (i = 0; i < igraph_vector_ptr_size(&eattrs); i++) {
        igraph_attribute_record_t *rec = VECTOR(eattrs)[i];
        if (rec->type == IGRAPH_ATTRIBUTE_NUMERIC) {
            igraph_vector_t *vec = (igraph_vector_t*)rec->value;
            long int origsize = igraph_vector_size(vec);
            igraph_vector_resize(vec, context.actedge);
            for (j = origsize; j < context.actedge; j++) {
                VECTOR(*vec)[j] = IGRAPH_NAN;
            }
        } else if (rec->type == IGRAPH_ATTRIBUTE_STRING) {
            igraph_strvector_t *strvec = (igraph_strvector_t*)rec->value;
            long int origsize = igraph_strvector_size(strvec);
            igraph_strvector_resize(strvec, context.actedge);
            for (j = origsize; j < context.actedge; j++) {
                igraph_strvector_set(strvec, j, "");
            }
        }
    }

    IGRAPH_CHECK(igraph_empty(graph, 0, context.directed));
    IGRAPH_FINALLY(igraph_destroy, graph);
    IGRAPH_CHECK(igraph_add_vertices(graph, context.vcount, &vattrs));
    IGRAPH_CHECK(igraph_add_edges(graph, &edges, &eattrs));

    for (i = 0; i < igraph_vector_ptr_size(&vattrs); i++) {
        igraph_attribute_record_t *rec = VECTOR(vattrs)[i];
        if (rec->type == IGRAPH_ATTRIBUTE_NUMERIC) {
            igraph_vector_t *vec = (igraph_vector_t*) rec->value;
            igraph_vector_destroy(vec);
            IGRAPH_FREE(vec);
        } else if (rec->type == IGRAPH_ATTRIBUTE_STRING) {
            igraph_strvector_t *strvec = (igraph_strvector_t *)rec->value;
            igraph_strvector_destroy(strvec);
            IGRAPH_FREE(strvec);
        }
        igraph_free( (char*)(rec->name));
        IGRAPH_FREE(rec);
    }

    for (i = 0; i < igraph_vector_ptr_size(&eattrs); i++) {
        igraph_attribute_record_t *rec = VECTOR(eattrs)[i];
        if (rec->type == IGRAPH_ATTRIBUTE_NUMERIC) {
            igraph_vector_t *vec = (igraph_vector_t*) rec->value;
            igraph_vector_destroy(vec);
            IGRAPH_FREE(vec);
        } else if (rec->type == IGRAPH_ATTRIBUTE_STRING) {
            igraph_strvector_t *strvec = (igraph_strvector_t *)rec->value;
            igraph_strvector_destroy(strvec);
            IGRAPH_FREE(strvec);
        }
        igraph_free( (char*)(rec->name));
        IGRAPH_FREE(rec);
    }

    igraph_vector_destroy(&edges);
    igraph_vector_ptr_destroy(&eattrs);
    igraph_trie_destroy(&eattrnames);
    igraph_vector_ptr_destroy(&vattrs);
    igraph_trie_destroy(&vattrnames);
    igraph_pajek_yylex_destroy(context.scanner);

    IGRAPH_FINALLY_CLEAN(7);
    return 0;
}

/* Order matters here! */
#define V_ID                0
#define V_X                 1
#define V_Y                 2
#define V_Z                 3
#define V_SHAPE             4
#define V_XFACT             5
#define V_YFACT             6
#define V_COLOR_RED         7
#define V_COLOR_GREEN       8
#define V_COLOR_BLUE        9
#define V_FRAMECOLOR_RED   10
#define V_FRAMECOLOR_GREEN 11
#define V_FRAMECOLOR_BLUE  12
#define V_LABELCOLOR_RED   13
#define V_LABELCOLOR_GREEN 14
#define V_LABELCOLOR_BLUE  15
#define V_LABELDIST        16
#define V_LABELDEGREE2     17
#define V_FRAMEWIDTH       18
#define V_FONTSIZE         19
#define V_ROTATION         20
#define V_RADIUS           21
#define V_DIAMONDRATIO     22
#define V_LABELDEGREE      23
#define V_VERTEXSIZE       24
#define V_FONT             25
#define V_URL              26
#define V_COLOR            27
#define V_FRAMECOLOR       28
#define V_LABELCOLOR       29
#define V_LAST             30

#define E_WEIGHT            0
#define E_COLOR_RED         1
#define E_COLOR_GREEN       2
#define E_COLOR_BLUE        3
#define E_ARROWSIZE         4
#define E_EDGEWIDTH         5
#define E_HOOK1             6
#define E_HOOK2             7
#define E_ANGLE1            8
#define E_ANGLE2            9
#define E_VELOCITY1        10
#define E_VELOCITY2        11
#define E_ARROWPOS         12
#define E_LABELPOS         13
#define E_LABELANGLE       14
#define E_LABELANGLE2      15
#define E_LABELDEGREE      16
#define E_FONTSIZE         17
#define E_ARROWTYPE        18
#define E_LINEPATTERN      19
#define E_LABEL            20
#define E_LABELCOLOR       21
#define E_COLOR            22
#define E_LAST             23

static int igraph_i_pajek_escape(char* src, char** dest) {
    long int destlen = 0;
    igraph_bool_t need_escape = 0;

    /* Determine whether the string contains characters to be escaped */
    char *s, *d;
    for (s = src; *s; s++, destlen++) {
        if (*s == '\\') {
            need_escape = 1;
            destlen++;
        } else if (*s == '"') {
            need_escape = 1;
            destlen++;
        } else if (!isalnum(*s)) {
            need_escape = 1;
        }
    }

    if (!need_escape) {
        /* At this point, we know that the string does not contain any chars
         * that would warrant escaping. Therefore, we simply quote it and
         * return the quoted string. This is necessary because Pajek uses some
         * reserved words in its format (like 'c' standing for color) and they
         * have to be quoted as well.
         */
        *dest = IGRAPH_CALLOC(destlen + 3, char);
        if (!*dest) {
            IGRAPH_ERROR("Not enough memory", IGRAPH_ENOMEM);
        }

        d = *dest;
        strcpy(d + 1, src);
        d[0] = d[destlen + 1] = '"';
        d[destlen + 2] = 0;
        return IGRAPH_SUCCESS;
    }

    *dest = IGRAPH_CALLOC(destlen + 3, char);
    if (!*dest) {
        IGRAPH_ERROR("Not enough memory", IGRAPH_ENOMEM);
    }

    d = *dest;
    *d = '"'; d++;

    for (s = src; *s; s++, d++) {
        switch (*s) {
        case '\\':
        case '"':
            *d = '\\'; d++;
            *d = *s;
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
 * </para><para>
 * The Pajek vertex and edge parameters (like color) are determined by
 * the attributes of the vertices and edges, of course this requires
 * an attribute handler to be installed. The names of the
 * corresponding vertex and edge attributes are listed at \ref
 * igraph_read_graph_pajek(), e.g. the \c color vertex attributes
 * determines the color (\c c in Pajek) parameter.
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
 * \param graph The graph object to write.
 * \param outstream The file to write to. It should be opened and
 * writable. Make sure that you open the file in binary format if you use MS Windows,
 * otherwise end of line characters will be messed up. (igraph will be able
 * to read back these messed up files, but Pajek won't.)
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

int igraph_write_graph_pajek(const igraph_t *graph, FILE *outstream) {
    long int no_of_nodes = igraph_vcount(graph);
    long int i, j;

    igraph_attribute_type_t vtypes[V_LAST], etypes[E_LAST];
    igraph_bool_t write_vertex_attrs = 0;

    /* Same order as the #define's */
    const char *vnames[] = { "id", "x", "y", "z", "shape", "xfact", "yfact",
                             "", "", "", "", "", "", "", "", "",
                             "labeldist", "labeldegree2", "framewidth",
                             "fontsize", "rotation", "radius",
                             "diamondratio", "labeldegree", "vertexsize",
                             "font", "url", "color", "framecolor",
                             "labelcolor"
                           };

    const char *vnumnames[] = { "xfact", "yfact", "labeldist",
                                "labeldegree2", "framewidth", "fontsize",
                                "rotation", "radius", "diamondratio",
                                "labeldegree", "vertexsize"
                              };
    const char *vnumnames2[] = { "x_fact", "y_fact", "lr", "lphi", "bw",
                                 "fos", "phi", "r", "q", "la", "size"
                               };
    const char *vstrnames[] = { "font", "url", "color", "framecolor",
                                "labelcolor"
                              };
    const char *vstrnames2[] = { "font", "url", "ic", "bc", "lc" };

    const char *enames[] = { "weight", "", "", "",
                             "arrowsize", "edgewidth", "hook1", "hook2",
                             "angle1", "angle2", "velocity1", "velocity2",
                             "arrowpos", "labelpos", "labelangle",
                             "labelangle2", "labeldegree", "fontsize",
                             "arrowtype", "linepattern", "label", "labelcolor",
                             "color"
                           };
    const char *enumnames[] = { "arrowsize", "edgewidth", "hook1", "hook2",
                                "angle1", "angle2", "velocity1", "velocity2",
                                "arrowpos", "labelpos", "labelangle",
                                "labelangle2", "labeldegree", "fontsize"
                              };
    const char *enumnames2[] = { "s", "w", "h1", "h2", "a1", "a2", "k1", "k2",
                                 "ap", "lp", "lr", "lphi", "la", "fos"
                               };
    const char *estrnames[] = { "arrowtype", "linepattern", "label",
                                "labelcolor", "color"
                              };
    const char *estrnames2[] = { "a", "p", "l", "lc", "c" };

    /* Newer Pajek versions support both Unix and Windows-style line endings,
     * so we just use Unix style. This will get converted to CRLF on Windows
     * when the file is opened in text mode */
    const char *newline = "\n";

    igraph_es_t es;
    igraph_eit_t eit;

    igraph_vector_t numv;
    igraph_strvector_t strv;

    igraph_vector_t ex_numa;
    igraph_vector_t ex_stra;
    igraph_vector_t vx_numa;
    igraph_vector_t vx_stra;

    char *s, *escaped;

    igraph_bool_t bipartite = 0;
    igraph_vector_int_t bip_index, bip_index2;
    igraph_vector_bool_t bvec;
    long int notop = 0, nobottom = 0;

    IGRAPH_VECTOR_INIT_FINALLY(&numv, 1);
    IGRAPH_STRVECTOR_INIT_FINALLY(&strv, 1);

    IGRAPH_VECTOR_INIT_FINALLY(&ex_numa, 0);
    IGRAPH_VECTOR_INIT_FINALLY(&ex_stra, 0);
    IGRAPH_VECTOR_INIT_FINALLY(&vx_numa, 0);
    IGRAPH_VECTOR_INIT_FINALLY(&vx_stra, 0);

    /* Check if graph is bipartite */
    if (igraph_i_attribute_has_attr(graph, IGRAPH_ATTRIBUTE_VERTEX, "type")) {
        igraph_attribute_type_t type_type;
        igraph_i_attribute_gettype(graph, &type_type, IGRAPH_ATTRIBUTE_VERTEX,
                                   "type");
        if (type_type == IGRAPH_ATTRIBUTE_BOOLEAN) {
            int bptr = 0, tptr = 0;
            bipartite = 1; write_vertex_attrs = 1;
            /* Count top and bottom vertices, we go over them twice,
            because we want to keep their original order */
            IGRAPH_CHECK(igraph_vector_int_init(&bip_index, no_of_nodes));
            IGRAPH_FINALLY(igraph_vector_int_destroy, &bip_index);
            IGRAPH_CHECK(igraph_vector_int_init(&bip_index2, no_of_nodes));
            IGRAPH_FINALLY(igraph_vector_int_destroy, &bip_index2);
            IGRAPH_CHECK(igraph_vector_bool_init(&bvec, 1));
            IGRAPH_FINALLY(igraph_vector_bool_destroy, &bvec);
            for (i = 0; i < no_of_nodes; i++) {
                IGRAPH_CHECK(igraph_i_attribute_get_bool_vertex_attr(graph,
                             "type", igraph_vss_1((igraph_integer_t) i), &bvec));
                if (VECTOR(bvec)[0]) {
                    notop++;
                } else {
                    nobottom++;
                }
            }
            for (i = 0, bptr = 0, tptr = (int) nobottom; i < no_of_nodes; i++) {
                IGRAPH_CHECK(igraph_i_attribute_get_bool_vertex_attr(graph,
                             "type", igraph_vss_1((igraph_integer_t) i), &bvec));
                if (VECTOR(bvec)[0]) {
                    VECTOR(bip_index)[tptr] = (int) i;
                    VECTOR(bip_index2)[i] = tptr;
                    tptr++;
                } else {
                    VECTOR(bip_index)[bptr] = (int) i;
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
        if (fprintf(outstream, "*Vertices %li %li%s", no_of_nodes, nobottom,
                    newline) < 0) {
            IGRAPH_ERROR("Cannot write pajek file", IGRAPH_EFILE);
        }
    } else {
        if (fprintf(outstream, "*Vertices %li%s", no_of_nodes, newline) < 0) {
            IGRAPH_ERROR("Cannot write pajek file", IGRAPH_EFILE);
        }
    }

    /* Check the vertex attributes */
    memset(vtypes, 0, sizeof(vtypes[0])*V_LAST);
    for (i = 0; i < V_LAST; i++) {
        if (igraph_i_attribute_has_attr(graph, IGRAPH_ATTRIBUTE_VERTEX,
                                        vnames[i])) {
            igraph_i_attribute_gettype(graph, &vtypes[i], IGRAPH_ATTRIBUTE_VERTEX,
                                       vnames[i]);
            write_vertex_attrs = 1;
        } else {
            vtypes[i] = (igraph_attribute_type_t) -1;
        }
    }
    for (i = 0; i < (long int) (sizeof(vnumnames) / sizeof(const char*)); i++) {
        igraph_attribute_type_t type;
        if (igraph_i_attribute_has_attr(graph, IGRAPH_ATTRIBUTE_VERTEX,
                                        vnumnames[i])) {
            igraph_i_attribute_gettype(graph, &type, IGRAPH_ATTRIBUTE_VERTEX,
                                       vnumnames[i]);
            if (type == IGRAPH_ATTRIBUTE_NUMERIC) {
                IGRAPH_CHECK(igraph_vector_push_back(&vx_numa, i));
            }
        }
    }
    for (i = 0; i < (long int) (sizeof(vstrnames) / sizeof(const char*)); i++) {
        igraph_attribute_type_t type;
        if (igraph_i_attribute_has_attr(graph, IGRAPH_ATTRIBUTE_VERTEX,
                                        vstrnames[i])) {
            igraph_i_attribute_gettype(graph, &type, IGRAPH_ATTRIBUTE_VERTEX,
                                       vstrnames[i]);
            if (type == IGRAPH_ATTRIBUTE_STRING) {
                IGRAPH_CHECK(igraph_vector_push_back(&vx_stra, i));
            }
        }
    }

    /* Write vertices */
    if (write_vertex_attrs) {
        for (i = 0; i < no_of_nodes; i++) {
            long int id = bipartite ? VECTOR(bip_index)[i] : i;

            /* vertex id */
            fprintf(outstream, "%li", i + 1);
            if (vtypes[V_ID] == IGRAPH_ATTRIBUTE_NUMERIC) {
                igraph_i_attribute_get_numeric_vertex_attr(graph, vnames[V_ID],
                        igraph_vss_1((igraph_integer_t) id), &numv);
                fputs(" \"", outstream);
                igraph_real_fprintf_precise(outstream, VECTOR(numv)[0]);
                fputc('"', outstream);
            } else if (vtypes[V_ID] == IGRAPH_ATTRIBUTE_STRING) {
                igraph_i_attribute_get_string_vertex_attr(graph, vnames[V_ID],
                        igraph_vss_1((igraph_integer_t) id), &strv);
                igraph_strvector_get(&strv, 0, &s);
                IGRAPH_CHECK(igraph_i_pajek_escape(s, &escaped));
                fprintf(outstream, " %s", escaped);
                IGRAPH_FREE(escaped);
            } else {
                fprintf(outstream, " \"%li\"", id + 1);
            }

            /* coordinates */
            if (vtypes[V_X] == IGRAPH_ATTRIBUTE_NUMERIC &&
                vtypes[V_Y] == IGRAPH_ATTRIBUTE_NUMERIC) {
                igraph_i_attribute_get_numeric_vertex_attr(graph, vnames[V_X],
                        igraph_vss_1((igraph_integer_t) id), &numv);
                fputc(' ', outstream);
                igraph_real_fprintf_precise(outstream, VECTOR(numv)[0]);
                igraph_i_attribute_get_numeric_vertex_attr(graph, vnames[V_Y],
                        igraph_vss_1((igraph_integer_t) id), &numv);
                fputc(' ', outstream);
                igraph_real_fprintf_precise(outstream, VECTOR(numv)[0]);
                if (vtypes[V_Z] == IGRAPH_ATTRIBUTE_NUMERIC) {
                    igraph_i_attribute_get_numeric_vertex_attr(graph, vnames[V_Z],
                            igraph_vss_1((igraph_integer_t) id), &numv);
                    fputc(' ', outstream);
                    igraph_real_fprintf_precise(outstream, VECTOR(numv)[0]);
                }
            }

            /* shape */
            if (vtypes[V_SHAPE] == IGRAPH_ATTRIBUTE_STRING) {
                igraph_i_attribute_get_string_vertex_attr(graph, vnames[V_SHAPE],
                        igraph_vss_1((igraph_integer_t) id), &strv);
                igraph_strvector_get(&strv, 0, &s);
                IGRAPH_CHECK(igraph_i_pajek_escape(s, &escaped));
                fprintf(outstream, " %s", escaped);
                IGRAPH_FREE(escaped);
            }

            /* numeric parameters */
            for (j = 0; j < igraph_vector_size(&vx_numa); j++) {
                int idx = (int) VECTOR(vx_numa)[j];
                igraph_i_attribute_get_numeric_vertex_attr(graph, vnumnames[idx],
                        igraph_vss_1((igraph_integer_t) id), &numv);
                fprintf(outstream, " %s ", vnumnames2[idx]);
                igraph_real_fprintf_precise(outstream, VECTOR(numv)[0]);
            }

            /* string parameters */
            for (j = 0; j < igraph_vector_size(&vx_stra); j++) {
                int idx = (int) VECTOR(vx_stra)[j];
                igraph_i_attribute_get_string_vertex_attr(graph, vstrnames[idx],
                        igraph_vss_1((igraph_integer_t) id), &strv);
                igraph_strvector_get(&strv, 0, &s);
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
    for (i = 0; i < E_LAST; i++) {
        if (igraph_i_attribute_has_attr(graph, IGRAPH_ATTRIBUTE_EDGE,
                                        enames[i])) {
            igraph_i_attribute_gettype(graph, &etypes[i], IGRAPH_ATTRIBUTE_EDGE,
                                       enames[i]);
        } else {
            etypes[i] = (igraph_attribute_type_t) -1;
        }
    }
    for (i = 0; i < (long int) (sizeof(enumnames) / sizeof(const char*)); i++) {
        igraph_attribute_type_t type;
        if (igraph_i_attribute_has_attr(graph, IGRAPH_ATTRIBUTE_EDGE,
                                        enumnames[i])) {
            igraph_i_attribute_gettype(graph, &type, IGRAPH_ATTRIBUTE_EDGE,
                                       enumnames[i]);
            if (type == IGRAPH_ATTRIBUTE_NUMERIC) {
                IGRAPH_CHECK(igraph_vector_push_back(&ex_numa, i));
            }
        }
    }
    for (i = 0; i < (long int) (sizeof(estrnames) / sizeof(const char*)); i++) {
        igraph_attribute_type_t type;
        if (igraph_i_attribute_has_attr(graph, IGRAPH_ATTRIBUTE_EDGE,
                                        estrnames[i])) {
            igraph_i_attribute_gettype(graph, &type, IGRAPH_ATTRIBUTE_EDGE,
                                       estrnames[i]);
            if (type == IGRAPH_ATTRIBUTE_STRING) {
                IGRAPH_CHECK(igraph_vector_push_back(&ex_stra, i));
            }
        }
    }

    for (i = 0; !IGRAPH_EIT_END(eit); IGRAPH_EIT_NEXT(eit), i++) {
        long int edge = IGRAPH_EIT_GET(eit);
        igraph_integer_t from, to;
        igraph_edge(graph, (igraph_integer_t) edge, &from,  &to);
        if (bipartite) {
            from = VECTOR(bip_index2)[from];
            to  = VECTOR(bip_index2)[to];
        }
        fprintf(outstream, "%li %li", (long int) from + 1, (long int) to + 1);

        /* Weights */
        if (etypes[E_WEIGHT] == IGRAPH_ATTRIBUTE_NUMERIC) {
            igraph_i_attribute_get_numeric_edge_attr(graph, enames[E_WEIGHT],
                    igraph_ess_1((igraph_integer_t) edge), &numv);
            fputc(' ', outstream);
            igraph_real_fprintf_precise(outstream, VECTOR(numv)[0]);
        }

        /* numeric parameters */
        for (j = 0; j < igraph_vector_size(&ex_numa); j++) {
            int idx = (int) VECTOR(ex_numa)[j];
            igraph_i_attribute_get_numeric_edge_attr(graph, enumnames[idx],
                    igraph_ess_1((igraph_integer_t) edge), &numv);
            fprintf(outstream, " %s ", enumnames2[idx]);
            igraph_real_fprintf_precise(outstream, VECTOR(numv)[0]);
        }

        /* string parameters */
        for (j = 0; j < igraph_vector_size(&ex_stra); j++) {
            int idx = (int) VECTOR(ex_stra)[j];
            igraph_i_attribute_get_string_edge_attr(graph, estrnames[idx],
                                                    igraph_ess_1((igraph_integer_t) edge), &strv);
            igraph_strvector_get(&strv, 0, &s);
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

    igraph_vector_destroy(&ex_numa);
    igraph_vector_destroy(&ex_stra);
    igraph_vector_destroy(&vx_numa);
    igraph_vector_destroy(&vx_stra);
    igraph_strvector_destroy(&strv);
    igraph_vector_destroy(&numv);
    IGRAPH_FINALLY_CLEAN(6);
    return 0;
}
