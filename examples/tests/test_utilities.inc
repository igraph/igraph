#ifndef TEST_UTILITIES_INC
#define TEST_UTILITIES_INC

/*
 * This file contains functions that are useful when writing tests.
 * Include it in the test program using #include "test_utilities.inc"
 */

#include <igraph.h>
#include <stdio.h>

/* Print elements of a vector. Use parentheses to make it clear when a vector has size zero. */
void print_vector(const igraph_vector_t *v, FILE *f) {
  long i;
  fprintf(f, "(");
  for (i=0; i < igraph_vector_size(v); i++) {
    fprintf(f, " %f", VECTOR(*v)[i]);
  }
  fprintf(f, " )\n");
}


/* Round elements of a vector to integers and print them. */
/* This is meant to be used when the elements of a vector are integer values. */
void print_vector_round(const igraph_vector_t *v, FILE *f) {
  long i;
  fprintf(f, "(");
  for (i=0; i < igraph_vector_size(v); i++) {
    fprintf(f, " %li", (long int) VECTOR(*v)[i]);
  }
  fprintf(f, " )\n");
}


/* Print elements of an integer vector */
void print_vector_int(const igraph_vector_int_t *v, FILE *f) {
  long i;
  fprintf(f, "(");
  for (i=0; i < igraph_vector_int_size(v); i++) {
    fprintf(f, " %d", VECTOR(*v)[i]);
  }
  fprintf(f, " )\n");
}


/* Print a graph. Use brackets to make it obvious when the edge list is empty. */
void print_graph(const igraph_t *graph, FILE *f) {
    long ecount = igraph_ecount(graph);
    long vcount = igraph_vcount(graph);
    long i;

    fprintf(f, "directed: %s\n", igraph_is_directed(graph) ? "true" : "false");
    fprintf(f, "vcount: %ld\n", vcount);
    fprintf(f, "edges: {\n");
    for (i=0; i < ecount; ++i)
        fprintf(f, "%d %d\n", IGRAPH_FROM(graph, i), IGRAPH_TO(graph, i));
    fprintf(f, "}\n");
}

#endif /* TEST_UTILITIES_INC */
