
#ifndef CLIQUER_REORDER_H
#define CLIQUER_REORDER_H

#include "set.h"
#include "graph.h"

extern void reorder_set(set_t s,int *order);
extern void reorder_graph(graph_t *g, int *order);
extern int *reorder_duplicate(int *order,int n);
extern void reorder_invert(int *order,int n);
extern void reorder_reverse(int *order,int n);
extern int *reorder_ident(int n);
extern boolean reorder_is_bijection(int *order,int n);


#define reorder_by_default reorder_by_greedy_coloring
extern int *reorder_by_greedy_coloring(graph_t *g, boolean weighted);
extern int *reorder_by_weighted_greedy_coloring(graph_t *g, boolean weighted);
extern int *reorder_by_unweighted_greedy_coloring(graph_t *g,boolean weighted);
extern int *reorder_by_degree(graph_t *g, boolean weighted);
extern int *reorder_by_random(graph_t *g, boolean weighted);
extern int *reorder_by_ident(graph_t *g, boolean weighted);
extern int *reorder_by_reverse(graph_t *g, boolean weighted);

#endif /* !CLIQUER_REORDER_H */
