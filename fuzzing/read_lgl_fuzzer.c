/*
  Heuristic graph coloring algorithms.
  Copyright (C) 2020 Szabolcs Horvat <szhorvat@gmail.com>

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

#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include "igraph.h"
#include <stdio.h>

int LLVMFuzzerTestOneInput(const uint8_t *data, size_t size){
	if(size<5) return 0;

	igraph_set_error_handler(igraph_error_handler_ignore);

	// Create input file	
	char filename[256];
	sprintf(filename, "/tmp/libfuzzer.gml");
	FILE *fp = fopen(filename, "wb");
	if (!fp) return 0;
	fwrite(data, size, 1, fp);
	fclose(fp);
	
	// Read input file
	FILE *ifile;
	ifile = fopen("/tmp/libfuzzer.gml", "r");
	if(ifile == 0){
		unlink(filename);
		return 0;
	}
	
	// Do the fuzzing	
	igraph_t g;
	igraph_read_graph_gml(&g, ifile);
	fclose(ifile);
	
	// Clean up
	igraph_destroy(&g);
	unlink(filename);
	return 0;
}
