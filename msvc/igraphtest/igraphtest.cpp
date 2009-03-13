// igraphtest.cpp : Defines the entry point for the console application.

#include <stdio.h>
#include "igraph.h"

int main(int argc, char* argv[])
{
	igraph_t g;
	FILE *outfile;

	igraph_barabasi_game(&g, 1000, /*m=*/2, /*outseq=*/0, /*outpref=*/0, /*directed=*/1);
	igraph_simplify(&g, /*multiple=*/1, /*loops=*/0);
	outfile=fopen("out.txt", "w");
	igraph_write_graph_edgelist(&g, outfile);
	fclose(outfile);
	igraph_destroy(&g);
	return 0;
}

