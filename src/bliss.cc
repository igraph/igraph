/*
 Copyright (C) 2003-2006 Tommi Junttila

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License version 2
 as published by the Free Software Foundation.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA.
*/

/* FSF address fixed in the above notice on 1 Oct 2009 by Tamas Nepusz */

#include "bliss_timer.hh"
#include "bliss_graph.hh"
#include "bliss_kqueue.hh"
#include "bliss_utils.hh"

#include "igraph_types.h"
#include "igraph_topology.h"

#include <string.h>

using namespace igraph;
using namespace std;

/**
 * \function igraph_canonical_permutation
 * Canonical permutation using BLISS
 * 
 * This function computes the canonical permutation which transforms
 * the graph into a canonical form by using the BLISS algorithm.
 * 
 * \param graph The input graph, it is treated as undirected and the
 *    multiple edges are ignored.
 * \param labeling Pointer to a vector, the result is stored here. The
 *    permutation takes vertex 0 to the first element of the vector,
 *    vertex 1 to the second, etc. The vector will be resized as
 *    needed.
 * \param sh The split heuristics to be used in BLISS. See \ref
 *    igraph_bliss_sh_t.
 * \param info If not \c NULL then information on BLISS internals is
 *    stored here. See \ref igraph_bliss_info_t.
 * \return Error code.
 * 
 * Time complexity: exponential, in practice it is fast for many graphs.
 */
int igraph_canonical_permutation(const igraph_t *graph, igraph_vector_t *labeling,
				 igraph_bliss_sh_t sh, igraph_bliss_info_t *info) {
  Graph *g = Graph::from_igraph(graph);
  Stats stats;
  const unsigned int N=g->get_nof_vertices();
  unsigned int gsh=Graph::sh_flm;

  switch (sh) { 
  case IGRAPH_BLISS_F:    gsh= Graph::sh_f;   break;
  case IGRAPH_BLISS_FL:   gsh= Graph::sh_fl;  break;
  case IGRAPH_BLISS_FS:   gsh= Graph::sh_fs;  break;
  case IGRAPH_BLISS_FM:   gsh= Graph::sh_fm;  break;
  case IGRAPH_BLISS_FLM:  gsh= Graph::sh_flm; break;
  case IGRAPH_BLISS_FSM:  gsh= Graph::sh_fsm; break;
  }

  g->set_splitting_heuristics(gsh);
  const unsigned int *cl = g->canonical_form(stats);
  IGRAPH_CHECK(igraph_vector_resize(labeling, N));
  for (unsigned int i=0; i<N; i++) {
    VECTOR(*labeling)[i] = cl[i];
  }
  delete g;

  if (info) { 
    info->nof_nodes      = stats.nof_nodes;
    info->nof_leaf_nodes = stats.nof_leaf_nodes;
    info->nof_bad_nodes  = stats.nof_bad_nodes;
    info->nof_canupdates = stats.nof_canupdates;
    info->max_level      = stats.max_level;
    stats.group_size.tostring(&info->group_size);
  }
  
  return 0;
}

/**
 * \function igraph_automorphisms
 * Number of automorphisms using BLISS
 * 
 * The number of automorphisms of a graph is computed using BLISS. The
 * result is returned as part of the \p info structure, in tag \c
 * group_size. It is returned as a string, as it can be very high even
 * for relatively small graphs. If the GNU MP library is used then
 * this number is exact, otherwise a <type>long double</type> is used
 * and it is only approximate. See also \ref igraph_bliss_info_t.
 * 
 * \param graph The input graph, it is treated as undirected and the
 *    multiple edges are ignored.
 * \param sh The split heuristics to be used in BLISS. See \ref
 *    igraph_bliss_sh_t.
 * \param info The result is stored here, in particular in the \c
 *    group_size tag of \p info.
 * \return Error code.
 * 
 * Time complexity: exponential, in practice it is fast for many graphs.
 */
int igraph_automorphisms(const igraph_t *graph,
			 igraph_bliss_sh_t sh, igraph_bliss_info_t *info) {
  
  Graph *g = Graph::from_igraph(graph);
  Stats stats;
  unsigned int gsh=Graph::sh_flm;

  switch (sh) { 
  case IGRAPH_BLISS_F:    gsh= Graph::sh_f;   break;
  case IGRAPH_BLISS_FL:   gsh= Graph::sh_fl;  break;
  case IGRAPH_BLISS_FS:   gsh= Graph::sh_fs;  break;
  case IGRAPH_BLISS_FM:   gsh= Graph::sh_fm;  break;
  case IGRAPH_BLISS_FLM:  gsh= Graph::sh_flm; break;
  case IGRAPH_BLISS_FSM:  gsh= Graph::sh_fsm; break;
  }

  g->set_splitting_heuristics(gsh);
  g->find_automorphisms(stats);

  if (info) { 
    info->nof_nodes      = stats.nof_nodes;
    info->nof_leaf_nodes = stats.nof_leaf_nodes;
    info->nof_bad_nodes  = stats.nof_bad_nodes;
    info->nof_canupdates = stats.nof_canupdates;
    info->max_level      = stats.max_level;
    stats.group_size.tostring(&info->group_size);
  }
  delete g;
  
  return 0;
}

bool bliss_verbose = false;
// FILE *bliss_verbstr = stdout;

namespace igraph {

typedef enum {FORMAT_BIN = 0, FORMAT_ADJ} Format;
// static Format input_format;

// static char *infilename = 0;

// static bool opt_canonize = false;
// static char *opt_output_can_file = 0;
// static unsigned int sh = Graph::sh_fm;

// static void usage(FILE *fp, char *argv0)
// {
//   char *program_name;
  
//   program_name = strrchr(argv0, '/');
  
//   if(program_name) program_name++;
//   else program_name = argv0;
  
//   if(!*program_name) program_name = "bliss";
//   fprintf(fp, "bliss, version 0.35, compiled " __DATE__ "\n");
//   fprintf(fp, "Copyright 2003-2006 Tommi Junttila\n");
//   fprintf(fp,
// "%s [<graph file>]\n"
// "\n"
// "  -can        compute canonical form\n"
// "  -ocan=f     compute canonical form and output it in file f\n"
// //"  -v        switch verbose mode on\n"
// "  -sh=x       select splitting heuristics, where x is\n"
// "                f    first non-singleton cell\n"
// "                fl   first largest non-singleton cell\n"
// "                fs   first smallest non-singleton cell\n"
// "                fm   first maximally non-trivially connected non-singleton cell [default]\n"
// "                flm  first largest maximally non-trivially connected non-singleton cell\n"
// "                fsm  first smallest maximally non-trivially connected non-singleton cell\n"
//           ,program_name);
// }


// static void parse_options(int argc, char ** argv)
// {
//   for(int i = 1; i < argc; i++)
//     {
//       //if(strcmp(argv[i], "-v") == 0 || strcmp(argv[i], "-verbose") == 0)
//       //bliss_verbose = true;
//       /*
//       if(strcmp(argv[i], "-bin") == 0)
// 	input_format = FORMAT_BIN;
//       else if(strcmp(argv[i], "-adj") == 0)
// 	input_format = FORMAT_ADJ;
//       */
//       if(strcmp(argv[i], "-can") == 0)
// 	opt_canonize = true;
//       else if((strncmp(argv[i], "-ocan=", 6) == 0) && (strlen(argv[i]) > 6))
// 	{
// 	  opt_canonize = true;
// 	  opt_output_can_file = argv[i]+6;
// 	}
//       else if(strcmp(argv[i], "-sh=f") == 0)
// 	sh = Graph::sh_f;
//       else if(strcmp(argv[i], "-sh=fs") == 0)
// 	sh = Graph::sh_fs;
//       else if(strcmp(argv[i], "-sh=fl") == 0)
// 	sh = Graph::sh_fl;
//       else if(strcmp(argv[i], "-sh=fm") == 0)
// 	sh = Graph::sh_fm;
//       else if(strcmp(argv[i], "-sh=fsm") == 0)
// 	sh = Graph::sh_fsm;
//       else if(strcmp(argv[i], "-sh=flm") == 0)
// 	sh = Graph::sh_flm;
//       else if(argv[i][0] == '-') {
// 	fprintf(stderr, "unknown command line argument `%s'\n", argv[i]);
// 	usage(stderr, argv[0]);
// 	exit(1);
//       }
//       else {
// 	if(infilename) {
// 	  fprintf(stderr, "too many file arguments\n");
// 	  usage(stderr, argv[0]);
// 	  exit(1);
// 	}
// 	else {
// 	  infilename = argv[i];
// 	}
//       }
//     }
// }			    

}

// using namespace igraph;

// int main(int argc, char **argv)
// {
//   Timer t;
//   t.start();

//   parse_options(argc, argv);

//   Graph *g = 0;
  
//   FILE *infile = stdin;
//   if(infilename) {
//     if(input_format == FORMAT_BIN)
//       infile = fopen(infilename, "rb");
//     else
//       infile = fopen(infilename, "r");
//     if(!infile) {
//       fprintf(stderr, "cannot not open `%s' for input\n", infilename);
//       exit(1); }
//   }

//   g = Graph::read_dimacs(infile);

//   if(infile != stdin)
//     fclose(infile);

//   if(!g)
//     return 0;
  
//   fprintf(stdout, "Graph read in %.2fs\n", t.get_intermediate());
  
// #ifdef DEBUG_PRINT_DOT
//   g->to_dot("debug_graph.dot");
// #endif

//   Stats stats;

//   g->set_splitting_heuristics(sh);

//   if(opt_canonize)
//     {
//       const unsigned int *cl = g->canonical_form(stats);
//       //fprintf(stdout, "Canonical labeling: ");
//       //print_permutation(stdout, g->get_nof_vertices(), cl);
//       //fprintf(stdout, "\n");
//       if(opt_output_can_file)
// 	{
// 	  Graph *cf = g->permute(cl);
// 	  FILE *fp = fopen(opt_output_can_file, "w");
// 	  if(!fp)
// 	    {
// 	      fprintf(stderr, "Can not open '%s' for outputting the canonical form", opt_output_can_file);
// 	      exit(1);
// 	    }
// 	  cf->print_dimacs(fp);
// 	  fclose(fp);
// 	  delete cf;
// 	}
//     }
//   else
//     {
//       g->find_automorphisms(stats);
//     }

//   printf("Nodes:\t\t%lu\n", stats.nof_nodes);
//   printf("Leaf nodes:\t%lu\n", stats.nof_leaf_nodes);
//   printf("Bad nodes:\t%lu\n", stats.nof_bad_nodes);
//   printf("Canrep updates:\t%lu\n", stats.nof_canupdates);
//   printf("Generators:\t%lu\n", stats.nof_generators);
//   printf("Max level:\t%lu\n", stats.max_level);
//   printf("|Aut|:\t\t"); stats.group_size.print(stdout); printf("\n");

//   t.stop();
//   printf("Total time:\t%.2fs\n", t.get_duration());

//   delete g;

//   return 0;
// }

