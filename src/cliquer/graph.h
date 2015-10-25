
#ifndef CLIQUER_GRAPH_H
#define CLIQUER_GRAPH_H

#include "set.h"

typedef struct _graph_t graph_t;
struct _graph_t {
	int n;             /* Vertices numbered 0...n-1 */
	set_t *edges;      /* A list of n sets (the edges). */
	int *weights;      /* A list of n vertex weights. */
};


#define GRAPH_IS_EDGE_FAST(g,i,j)  (SET_CONTAINS_FAST((g)->edges[(i)],(j)))
#define GRAPH_IS_EDGE(g,i,j) (((i)<((g)->n))?SET_CONTAINS((g)->edges[(i)], \
							  (j)):FALSE)
#define GRAPH_ADD_EDGE(g,i,j) do {            \
	SET_ADD_ELEMENT((g)->edges[(i)],(j)); \
	SET_ADD_ELEMENT((g)->edges[(j)],(i)); \
} while (FALSE)
#define GRAPH_DEL_EDGE(g,i,j) do {            \
	SET_DEL_ELEMENT((g)->edges[(i)],(j)); \
	SET_DEL_ELEMENT((g)->edges[(j)],(i)); \
} while (FALSE)


extern graph_t *graph_new(int n);
extern void graph_free(graph_t *g);
extern void graph_resize(graph_t *g, int size);
extern void graph_crop(graph_t *g);

extern boolean graph_weighted(graph_t *g);
extern int graph_edge_count(graph_t *g);

/*
extern graph_t *graph_read_dimacs(FILE *fp);
extern graph_t *graph_read_dimacs_file(char *file);
extern boolean graph_write_dimacs_ascii(graph_t *g, char *comment,FILE *fp);
extern boolean graph_write_dimacs_ascii_file(graph_t *g,char *comment,
					     char *file);
extern boolean graph_write_dimacs_binary(graph_t *g, char *comment,FILE *fp);
extern boolean graph_write_dimacs_binary_file(graph_t *g, char *comment,
					      char *file);
*/

extern void graph_print(graph_t *g);
extern boolean graph_test(graph_t *g, FILE *output);
extern int graph_test_regular(graph_t *g);

UNUSED_FUNCTION INLINE
static int graph_subgraph_weight(graph_t *g,set_t s) {
	int i,j;
	int count=0;
	setelement e;

	for (i=0; i<SET_ARRAY_LENGTH(s); i++) {
		if (s[i]) {
			e=s[i];
			for (j=0; j<ELEMENTSIZE; j++) {
				if (e&1)
					count+=g->weights[i*ELEMENTSIZE+j];
				e = e>>1;
			}
		}
	}
	return count;
}

UNUSED_FUNCTION INLINE
static int graph_vertex_degree(graph_t *g, int v) {
	return set_size(g->edges[v]);
}

#endif /* !CLIQUER_GRAPH_H */
