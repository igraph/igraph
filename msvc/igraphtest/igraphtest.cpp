// igraphtest.cpp : Defines the entry point for the console application.

#define IGRAPH_STATIC 1

#include <stdio.h>
#include "igraph.h"

int main(int argc, char* argv[])
{
	igraph_t g;
	FILE *outfile;

	igraph_barabasi_game(&g, 1000, /*power=*/ 1,/*m=*/2, /*outseq=*/0, 
			     /*outpref=*/0, /*A=*/ 1, /*directed=*/1, 
			     IGRAPH_BARABASI_PSUMTREE, 
			     /*start_from=*/ 0);
	igraph_simplify(&g, /*multiple=*/1, /*loops=*/0, /*edge_comb=*/ 0);
	outfile=fopen("out.txt", "w");
	igraph_write_graph_edgelist(&g, outfile);
	fclose(outfile);
	igraph_destroy(&g);
	return 0;
}

