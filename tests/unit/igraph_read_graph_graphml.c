/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2006-2022  The igraph development team

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
   Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA
   02110-1301 USA

*/
#include <igraph.h>
#include <stdio.h>
#include <unistd.h>     /* unlink */

void custom_warning_handler (const char *reason, const char *file,
                             int line) {
    IGRAPH_UNUSED(file);
    IGRAPH_UNUSED(line);
    printf("Warning: %s\n", reason);
}

void dump_graph(const char* header, const igraph_t* g) {
    fputs(header, stdout);
    printf("Vertices: %" IGRAPH_PRId "\n", igraph_vcount(g));
    printf("Edges: %" IGRAPH_PRId "\n", igraph_ecount(g));
    printf("Directed: %i\n", igraph_is_directed(g) ? 1 : 0);
    igraph_write_graph_edgelist(g, stdout);
}

void dump_vertex_attribute_bool(const char* name, const igraph_t* g) {
    igraph_integer_t i, n = igraph_vcount(g);

    printf("Vertex attribute '%s':", name);
    for (i = 0; i < n; i++) {
        printf(" %s", VAB(g, name, i) ? "true" : "false");
    }
    printf("\n");
}

void dump_vertex_attribute_numeric(const char* name, const igraph_t* g) {
    igraph_integer_t i, n = igraph_vcount(g);

    printf("Vertex attribute '%s':", name);
    for (i = 0; i < n; i++) {
        printf(" ");
        igraph_real_printf_precise(VAN(g, name, i));
    }
    printf("\n");
}

void dump_vertex_attribute_string(const char* name, const igraph_t* g) {
    igraph_integer_t i, n = igraph_vcount(g);

    printf("Vertex attribute '%s':", name);
    for (i = 0; i < n; i++) {
        printf(" '%s'", VAS(g, name, i));
    }
    printf("\n");
}

int main() {
    igraph_t g;
    igraph_error_handler_t* oldhandler;
    igraph_warning_handler_t* oldwarnhandler;
    int result;
    FILE *ifile, *ofile;

    igraph_set_attribute_table(&igraph_cattribute_table);

    /* GraphML */
    ifile = fopen("test.graphml", "r");
    if (ifile == 0) {
        return 10;
    }

    oldhandler = igraph_set_error_handler(igraph_error_handler_ignore);
    oldwarnhandler = igraph_set_warning_handler(custom_warning_handler);
    if ((result = igraph_read_graph_graphml(&g, ifile, 0))) {
        /* maybe it is simply disabled at compile-time */
        if (result == IGRAPH_UNIMPLEMENTED) {
            return 77;
        }
        return 1;
    }
    igraph_set_error_handler(oldhandler);

    fclose(ifile);

    /* Write it back */
    ofile = fopen("test2.graphml", "w");
    /* If we can't create the test file, just skip the test */
    if (ofile) {
        if ((result = igraph_write_graph_graphml(&g, ofile, /*prefixattr=*/ 1))) {
            printf("Received unexpected return code: %d\n", result);
            return 1;
        }
        fclose(ofile);
        unlink("test2.graphml");
    }
    dump_graph("Directed graph:\n", &g);
    igraph_destroy(&g);

    /* The same with undirected graph */
    ifile = fopen("test.graphml", "r");
    if ((result = igraph_read_graph_graphml(&g, ifile, 0))) {
        printf("Received unexpected return code: %d\n", result);
        return 1;
    }
    fclose(ifile);
    dump_graph("Undirected graph:\n", &g);
    igraph_destroy(&g);

    /* Test a GraphML file with default attributes */
    ifile = fopen("graphml-default-attrs.xml", "r");
    if ((result = igraph_read_graph_graphml(&g, ifile, 0))) {
        printf("Received unexpected return code: %d\n", result);
        return 1;
    }
    fclose(ifile);
    dump_graph("Graph with default attributes:\n", &g);
    dump_vertex_attribute_bool("type", &g);
    dump_vertex_attribute_string("gender", &g);
    dump_vertex_attribute_numeric("age", &g);
    dump_vertex_attribute_bool("retired", &g);
    igraph_destroy(&g);

    /* Test a GraphML file with namespaces */
    ifile = fopen("graphml-namespace.xml", "r");
    if ((result = igraph_read_graph_graphml(&g, ifile, 0))) {
        printf("Received unexpected return code: %d\n", result);
        return 1;
    }
    fclose(ifile);
    dump_graph("Graph with namespace:\n", &g);
    igraph_destroy(&g);

    /* Test a not-really-valid GraphML file as it has no namespace information */
    ifile = fopen("graphml-lenient.xml", "r");
    if ((result = igraph_read_graph_graphml(&g, ifile, 0))) {
        printf("Received unexpected return code: %d\n", result);
        return 1;
    }
    fclose(ifile);
    dump_graph("Graph without namespace information:\n", &g);
    igraph_destroy(&g);

    /* Test a GraphML file with excess whitespace around attribute values
     * (which we attempt to handle gracefully) */
    ifile = fopen("graphml-whitespace.xml", "r");
    if ((result = igraph_read_graph_graphml(&g, ifile, 0))) {
        printf("Received unexpected return code: %d\n", result);
        return 1;
    }
    fclose(ifile);
    dump_graph("Graph with whitespace in attributes:\n", &g);
    dump_vertex_attribute_bool("type", &g);
    dump_vertex_attribute_string("name", &g);
    dump_vertex_attribute_numeric("weight", &g);
    igraph_destroy(&g);

    /* Restore the old error handler */
    igraph_set_error_handler(igraph_error_handler_abort);

    /* Restore the old warning handler */
    igraph_set_warning_handler(oldwarnhandler);

    /* There were sometimes problems with this file */
    /* Only if called from R though, and only on random occasions, once in every
       ten reads. Do testing here doesn't make much sense, but if we have the file
       then let's do it anyway. */
    ifile = fopen("graphml-hsa05010.xml", "r");
    igraph_read_graph_graphml(&g, ifile, 0);
    fclose(ifile);
    igraph_destroy(&g);

    return 0;
}
