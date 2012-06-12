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

#include <cstdio>
#include <cassert>
#include <cctype>
#include <set>
#include <list>
#include <algorithm>
#include "bliss_defs.hh"
#include "bliss_timer.hh"
#include "bliss_graph.hh"
#include "bliss_partition.hh"
#include <climits>		// INT_MAX, etc

#include "igraph_datatype.h"
#include "igraph_interface.h"
#include "igraph_topology.h"
#include "igraph_statusbar.h"

using namespace std;

extern bool bliss_verbose;
// extern FILE *bliss_verbstr;

namespace igraph {

static const bool should_not_happen = false;

/*-------------------------------------------------------------------------
 *
 * Constructor and destructor routines for the abstract graph class
 *
 *-------------------------------------------------------------------------*/


AbstractGraph::AbstractGraph()
{
  /* Initialize stuff */
  first_path_labeling = 0;
  first_path_labeling_inv = 0;
  best_path_labeling = 0;
  best_path_labeling_inv = 0;
  first_path_automorphism = 0;
  best_path_automorphism = 0;
  //certificate = 0;
  in_search = false;
}


AbstractGraph::~AbstractGraph()
{
  if(first_path_labeling) {
    free(first_path_labeling); first_path_labeling = 0; }
  if(first_path_labeling_inv) {
    free(first_path_labeling_inv); first_path_labeling_inv = 0; }
  if(best_path_labeling) {
    free(best_path_labeling); best_path_labeling = 0; }
  if(best_path_labeling_inv) {
    free(best_path_labeling_inv); best_path_labeling_inv = 0; }
  if(first_path_automorphism) {
    free(first_path_automorphism); first_path_automorphism = 0; }
  if(best_path_automorphism) {
    free(best_path_automorphism); best_path_automorphism = 0; }
  //if(certificate) {
  //  free(certificate); certificate = 0; }
  while(!long_prune_fixed.empty())
    {
      delete long_prune_fixed.back();
      long_prune_fixed.pop_back();
    }
  while(!long_prune_mcrs.empty())
    {
      delete long_prune_mcrs.back();
      long_prune_mcrs.pop_back();
    }
}





/*-------------------------------------------------------------------------
 *
 * Routines for refinement to equitable partition
 *
 *-------------------------------------------------------------------------*/


void AbstractGraph::refine_to_equitable()
{
  assert(p.splitting_queue.is_empty());

  for(Cell *cell = p.first_cell; cell; cell = cell->next)
    {
      p.add_in_splitting_queue(cell);
    }

  return do_refine_to_equitable();
}


void AbstractGraph::refine_to_equitable(Cell *cell1)
{
  DEBUG_ASSERT(cell1->length == 1);

#ifdef EXPENSIVE_CONSISTENCY_CHECKS
  for(Cell *cell = p.first_cell; cell; cell = cell->next) {
    assert(cell->in_splitting_queue == false);
    assert(cell->in_neighbour_heap == false);
  }
#endif
  
  assert(p.splitting_queue.is_empty());

  p.add_in_splitting_queue(cell1);

  return do_refine_to_equitable();
}


void AbstractGraph::refine_to_equitable(Cell *cell1, Cell *cell2)
{
  DEBUG_ASSERT(cell1->length == 1);
  DEBUG_ASSERT(cell2->length == 1);

#ifdef EXPENSIVE_CONSISTENCY_CHECKS
  for(Cell *cell = p.first_cell; cell; cell = cell->next) {
    assert(cell->in_splitting_queue == false);
    assert(cell->in_neighbour_heap == false);
  }
#endif
  
  assert(p.splitting_queue.is_empty());

  p.add_in_splitting_queue(cell1);
  p.add_in_splitting_queue(cell2);

  return do_refine_to_equitable();
}


void AbstractGraph::do_refine_to_equitable()
{
  assert(!p.splitting_queue.is_empty());
  assert(neighbour_heap.is_empty());

  eqref_hash.reset();

  while(!p.splitting_queue.is_empty())
    {
      Cell *cell = p.splitting_queue.pop_front();
      DEBUG_ASSERT(cell->in_splitting_queue);
      cell->in_splitting_queue = false;

      if(cell->length == 1)
	{
	  if(in_search) {
	    if(first_path_automorphism) {
	      /* Build the (potential) automorphism on-the-fly */
	      assert(first_path_labeling_inv);
	      first_path_automorphism[first_path_labeling_inv[cell->first]] =
		p.elements[cell->first];
	    }
	    if(best_path_automorphism)
	      {
		/* Build the (potential) automorphism on-the-fly */
		assert(best_path_labeling_inv);
		best_path_automorphism[best_path_labeling_inv[cell->first]] =
		  p.elements[cell->first];
	      }
	  }
	  
	  bool worse = split_neighbourhood_of_unit_cell(cell);
	  if(in_search && worse)
	    goto worse_exit;
	}
      else
	{
	  split_neighbourhood_of_cell(cell);
	}
    }

  eqref_worse_than_certificate = false;
  return;

 worse_exit:
  /* Clear splitting_queue */
  p.clear_splitting_queue();
  eqref_worse_than_certificate = true;
  return;
}





/*-------------------------------------------------------------------------
 *
 * Routines for handling the canonical labeling
 *
 *-------------------------------------------------------------------------*/


void AbstractGraph::update_labeling(unsigned int * const labeling)
{
  const unsigned int N = get_nof_vertices();
  unsigned int *ep = p.elements;
  for(unsigned int i = 0; i < N; i++, ep++)
    labeling[*ep] = i;
}


void AbstractGraph::update_labeling_and_its_inverse(unsigned int * const labeling,
						    unsigned int * const labeling_inv)
{
  const unsigned int N = get_nof_vertices();
  unsigned int *ep = p.elements;
  unsigned int *clip = labeling_inv;

  for(unsigned int i = 0; i < N; ) {
    labeling[*ep] = i;
    i++;
    *clip = *ep;
    ep++;
    clip++;
  }
}



/*-------------------------------------------------------------------------
 *
 * Routines for handling automorphisms
 *
 *-------------------------------------------------------------------------*/


void AbstractGraph::reset_permutation(unsigned int *perm)
{
  const unsigned int N = get_nof_vertices();
  for(unsigned int i = 0; i < N; i++, perm++)
    *perm = i;
}


bool AbstractGraph::is_automorphism(unsigned int * const perm)
{
  assert(should_not_happen);
  return false;
}





/*-------------------------------------------------------------------------
 *
 * Long prune code
 *
 *-------------------------------------------------------------------------*/

void AbstractGraph::long_prune_init()
{
  const unsigned int N = get_nof_vertices();
  long_prune_temp.clear();
  long_prune_temp.resize(N);
#ifdef DEBUG
  for(unsigned int i = 0; i < N; i++)
    assert(long_prune_temp[i] == false);
#endif
  const unsigned int nof_fitting_in_max_mem =
    (long_prune_options_max_mem * 1024 * 1024) / (((N * 2) / 8)+1);
  long_prune_max_stored_autss = long_prune_options_max_stored_auts;
  /* Had some problems with g++ in using (a<b)?a:b when constants involved,
     so had to make this in a stupid way...*/
  if(nof_fitting_in_max_mem < long_prune_options_max_stored_auts)
    long_prune_max_stored_autss = nof_fitting_in_max_mem;

  while(!long_prune_fixed.empty())
    {
      delete long_prune_fixed.back();
      long_prune_fixed.pop_back();
    }
  while(!long_prune_mcrs.empty())
    {
      delete long_prune_mcrs.back();
      long_prune_mcrs.pop_back();
    }
  for(unsigned int i = 0; i < long_prune_max_stored_autss; i++)
    {
      long_prune_fixed.push_back(new std::vector<bool>(N));
      long_prune_mcrs.push_back(new std::vector<bool>(N));
    }
  long_prune_begin = 0;
  long_prune_end = 0;
}

void AbstractGraph::long_prune_swap(const unsigned int i, const unsigned int j)
{
  assert(long_prune_begin <= long_prune_end);
  assert(long_prune_end - long_prune_end <= long_prune_max_stored_autss);
  assert(i >= long_prune_begin);
  assert(i < long_prune_end);
  assert(j >= long_prune_begin);
  assert(j < long_prune_end);
  const unsigned int real_i = i % long_prune_max_stored_autss;
  const unsigned int real_j = j % long_prune_max_stored_autss;
  std::vector<bool> * tmp = long_prune_fixed[real_i];
  long_prune_fixed[real_i] = long_prune_fixed[real_j];
  long_prune_fixed[real_j] = tmp;
  tmp = long_prune_mcrs[real_i];
  long_prune_mcrs[real_i] = long_prune_mcrs[real_j];
  long_prune_mcrs[real_j] = tmp;
}

std::vector<bool> &AbstractGraph::long_prune_get_fixed(const unsigned int index)
{
  assert(long_prune_begin <= long_prune_end);
  assert(long_prune_end - long_prune_end <= long_prune_max_stored_autss);
  assert(index >= long_prune_begin);
  assert(index < long_prune_end);
  return *long_prune_fixed[index % long_prune_max_stored_autss];
}

std::vector<bool> &AbstractGraph::long_prune_get_mcrs(const unsigned int index)
{
  assert(long_prune_begin <= long_prune_end);
  assert(long_prune_end - long_prune_end <= long_prune_max_stored_autss);
  assert(index >= long_prune_begin);
  assert(index < long_prune_end);
  return *long_prune_mcrs[index % long_prune_max_stored_autss];
}


void AbstractGraph::long_prune_add_automorphism(const unsigned int *aut)
{
  if(long_prune_max_stored_autss == 0)
    return;

  const unsigned int N = get_nof_vertices();

#ifdef DEBUG
  assert(long_prune_temp.size() == N);
  for(unsigned int i = 0; i < N; i++)
    assert(long_prune_temp[i] == false);
#endif

  DEBUG_ASSERT(long_prune_fixed.size() == long_prune_mcrs.size());
  assert(long_prune_begin <= long_prune_end);
  assert(long_prune_end - long_prune_end <= long_prune_max_stored_autss);
  if(long_prune_end - long_prune_begin == long_prune_max_stored_autss)
    {
      long_prune_begin++;
    }
  long_prune_end++;
  std::vector<bool> &fixed = long_prune_get_fixed(long_prune_end-1);
  std::vector<bool> &mcrs = long_prune_get_mcrs(long_prune_end-1);

  for(unsigned int i = 0; i < N; i++)
    {
      fixed[i] = (aut[i] == i);
      if(!long_prune_temp[i])
	{
	  mcrs[i] = true;
	  unsigned int j = aut[i];
	  while(j != i)
	    {
	      assert(i <= j);
	      long_prune_temp[j] = true;
	      j = aut[j];
	    }
	}
      else
	{
	  mcrs[i] = false;
	}
      long_prune_temp[i] = false;
    }


#ifdef DEBUG
  for(unsigned int i = 0; i < N; i++)
    assert(long_prune_temp[i] == false);
#endif
}


/*-------------------------------------------------------------------------
 *
 * Routines for handling orbit information
 *
 *-------------------------------------------------------------------------*/


void AbstractGraph::update_orbit_information(Orbit &o, const unsigned int *p)
{
  const unsigned int N = get_nof_vertices();
  for(unsigned int i = 0; i < N; i++)
    if(p[i] != i)
      o.merge_orbits(i, p[i]);
}





/*-------------------------------------------------------------------------
 *
 * Print a permutation in cycle notation
 *
 *-------------------------------------------------------------------------*/


void AbstractGraph::print_permutation(FILE *fp, const unsigned int *perm)
{
  const unsigned int N = get_nof_vertices();
  for(unsigned int i = 0; i < N; i++) {
    unsigned int j = perm[i];
    if(j == i)
      continue;
    bool is_first = true;
    while(j != i) {
      if(j < i) {
	is_first = false;
	break;
      }
      j = perm[j];
    }
    if(!is_first)
      continue;
    fprintf(fp, "(%u,", i);
    j = perm[i];
    while(j != i) {
      fprintf(fp, "%u", j);
      j = perm[j];
      if(j != i)
	fprintf(fp, ",");
    }
    fprintf(fp, ")");
  }
}





/*-------------------------------------------------------------------------
 *
 * The actual backtracking search
 *
 *-------------------------------------------------------------------------*/


typedef struct {
  int split_element;
  unsigned int split_cell_first;
  unsigned int refinement_stack_size;
  unsigned int certificate_index;

  bool in_first_path;
  bool in_best_path;
  bool equal_to_first_path;
  int cmp_to_best_path;

  bool needs_long_prune;
  unsigned int long_prune_begin;
  std::set<unsigned int, std::less<unsigned int> > long_prune_redundant;
  
  EqrefHash eqref_hash;
  unsigned int subcertificate_length;
} LevelInfo;



typedef struct t_path_info {
  unsigned int splitting_element;
  unsigned int certificate_index;
  unsigned int subcertificate_length;
  EqrefHash eqref_hash;
} PathInfo;


void AbstractGraph::search(const bool canonical, Stats &stats)
{
  const unsigned int N = get_nof_vertices();

  // const bool write_automorphisms = 0;

  unsigned int all_same_level = UINT_MAX;

  p.graph = this;

  /*
   * Must be done!
   */
  remove_duplicate_edges();

  /*
   * Reset search statistics
   */
  stats.group_size.assign(1);
  stats.nof_nodes = 1;
  stats.nof_leaf_nodes = 1;
  stats.nof_bad_nodes = 0;
  stats.nof_canupdates = 0;
  stats.nof_generators = 0;
  stats.max_level = 0;

  if(first_path_labeling)
    {
      free(first_path_labeling);
      first_path_labeling = 0;
    }
  if(first_path_labeling_inv)
    {
      free(first_path_labeling_inv);
      first_path_labeling_inv = 0;
   }
  if(first_path_automorphism)
    {
      free(first_path_automorphism);
      first_path_automorphism = 0;
   }

  if(best_path_labeling)
    {
      free(best_path_labeling);
      best_path_labeling = 0;
    }
  if(best_path_labeling_inv)
    {
      free(best_path_labeling_inv);
      best_path_labeling_inv = 0;
   }
  if(best_path_automorphism)
    {
      free(best_path_automorphism);
      best_path_automorphism = 0;
   }

  if(N == 0)
    return;

  p.init(N);
  neighbour_heap.init(N);

  in_search = false;

  p.level = 0;

  Timer t1;
  t1.start();

  make_initial_equitable_partition();

#if defined(VERIFY_EQUITABLEDNESS)
  assert(is_equitable());
#endif

  t1.stop();
  
  igraph_statusf("Initial partition computed in %.2fs", 0,
		 t1.get_duration());
  
  /*
   * Allocate space for the labelings
   */
  if(first_path_labeling)
    free(first_path_labeling);
  first_path_labeling = (unsigned int*)calloc(N, sizeof(unsigned int));
  if(best_path_labeling)
    free(best_path_labeling);
  best_path_labeling = (unsigned int*)calloc(N, sizeof(unsigned int));

  /*
   * Are there any non-singleton cells?
   */
  if(p.is_discrete())
    {
      update_labeling(best_path_labeling);
      return;
    }

  //p.print_signature(stderr); fprintf(stderr, "\n");

  /*
   * Allocate space for the inverses of the labelings
   */
  if(first_path_labeling_inv)
    free(first_path_labeling_inv);
  first_path_labeling_inv = (unsigned int*)calloc(N, sizeof(unsigned int));
  if(best_path_labeling_inv)
    free(best_path_labeling_inv);
  best_path_labeling_inv = (unsigned int*)calloc(N, sizeof(unsigned int));


  /*
   * Allocate space for the automorphisms
   */
  if(first_path_automorphism) free(first_path_automorphism);
  first_path_automorphism = (unsigned int*)malloc(N * sizeof(unsigned int));
  if(best_path_automorphism) free(best_path_automorphism);
  best_path_automorphism = (unsigned int*)malloc(N * sizeof(unsigned int));


  /*
   * Initialize orbit information
   */
  first_path_orbits.init(N);
  best_path_orbits.init(N);

  /*
   * Initialize certificate memory
   */
  initialize_certificate();
  //assert(certificate);
  assert(certificate_index == 0);

  LevelInfo info;
  std::vector<LevelInfo> search_stack;
  std::vector<PathInfo> first_path_info;
  std::vector<PathInfo> best_path_info;

  search_stack.clear();
  p.refinement_stack.clean();
  assert(neighbour_heap.is_empty());

  /*
   * Initialize long prune
   */
  long_prune_init();

  /*
   * Build the first level info
   */
  info.split_cell_first = find_next_cell_to_be_splitted(p.first_cell)->first;
  info.split_element = -1;
  info.refinement_stack_size = p.refinement_stack.size();
  info.certificate_index = 0;
  info.in_first_path = false;
  info.in_best_path = false;
  info.long_prune_begin = 0;
  search_stack.push_back(info);

  /*
   * Set status and global flags for search related procedures
   */
  in_search = true;
  refine_compare_certificate = false;
  stats.nof_leaf_nodes = 0;


#ifdef PRINT_SEARCH_TREE_DOT
  dotty_output = fopen("debug_stree.dot", "w");
  fprintf(dotty_output, "digraph stree {\n");
  fprintf(dotty_output, "\"n\" [label=\"");
  fprintf(dotty_output, "M"); //p.print(dotty_output);
  fprintf(dotty_output, "\"];\n");
#endif

  p.consistency_check();

  /*
   * The actual backtracking search
   */
  while(!search_stack.empty()) 
    {
      info = search_stack.back();
      search_stack.pop_back();

      p.consistency_check();

      /*
       * Restore partition, certificate index, and split cell
       */
      p.unrefine(p.level, info.refinement_stack_size);
      assert(info.certificate_index <= certificate_size);
      certificate_index = info.certificate_index;
      certificate_current_path.resize(certificate_index);
      Cell * const cell = p.element_to_cell_map[p.elements[info.split_cell_first]];
      assert(cell->length > 1);

      p.consistency_check();

      if(p.level > 0 && !info.in_first_path)
	{
	  if(info.split_element == -1)
	    {
	      info.needs_long_prune = true;
	    }
	  else if(info.needs_long_prune)
	    {
	      info.needs_long_prune = false;
	      /* THIS IS A QUITE HORRIBLE HACK! */
	      unsigned int begin = (info.long_prune_begin>long_prune_begin)?info.long_prune_begin:long_prune_begin;
	      for(unsigned int i = begin; i < long_prune_end; i++)
		{
		  const std::vector<bool> &fixed = long_prune_get_fixed(i);
		  bool fixes_all = true;
		  for(unsigned int l = 0; l < p.level; l++)
		    {
		      if(fixed[search_stack[l].split_element] == false)
			{
			  fixes_all = false;
			  break;
			}
		    }
		  if(!fixes_all)
		    {
		      long_prune_swap(begin, i);
		      begin++;
		      info.long_prune_begin = begin;
		      continue;
		    }
		  const std::vector<bool> &mcrs = long_prune_get_mcrs(i);
		  unsigned int *ep = p.elements + cell->first;
		  for(unsigned int j = cell->length; j > 0; j--, ep++) {
		    if(mcrs[*ep] == false)
		      {
			info.long_prune_redundant.insert(*ep);
		      }
		  }
		}
	    }
	}

      /*
       * Find the next smallest element in cell
       */
      unsigned int next_split_element = UINT_MAX;
      unsigned int *next_split_element_pos = 0;
      unsigned int *ep = p.elements + cell->first;
      if(info.in_first_path)
	{
	  /* Find the next larger splitting element that is a mor */
	  for(unsigned int i = cell->length; i > 0; i--, ep++) {
	    if((int)(*ep) > info.split_element &&
	       *ep < next_split_element &&
	       first_path_orbits.is_minimal_representative(*ep)) {
	      next_split_element = *ep;
	      next_split_element_pos = ep;
	    }
	  }
	}
      else if(info.in_best_path)
	{
	  /* Find the next larger splitting element that is a mor */
	  for(unsigned int i = cell->length; i > 0; i--, ep++) {
	    if((int)(*ep) > info.split_element &&
	       *ep < next_split_element &&
	       best_path_orbits.is_minimal_representative(*ep) &&
	       (info.long_prune_redundant.find(*ep) ==
		info.long_prune_redundant.end())) {
	      next_split_element = *ep;
	      next_split_element_pos = ep;
	    }
	  }
	}
      else
	{
	  /* Find the next larger splitting element */
	  for(unsigned int i = cell->length; i > 0; i--, ep++) {
	    if((int)(*ep) > info.split_element &&
	       *ep < next_split_element &&
	       (info.long_prune_redundant.find(*ep) ==
		info.long_prune_redundant.end())) {
	      next_split_element = *ep;
	      next_split_element_pos = ep;
	    }
	  }
	}
      if(next_split_element == UINT_MAX)
	{
	  /*
	   * No more splitting elements (unexplored children) in the cell
	   */
	  /* Update group size if required */
	  if(info.in_first_path == true) {
	    const unsigned int index =
	      first_path_orbits.orbit_size(first_path_info[p.level].splitting_element);
	    stats.group_size.multiply(index);
	    /*
	     * Update all_same_level
	     */
	    if(index == cell->length && all_same_level == p.level+1)
	      all_same_level = p.level;
	    igraph_statusf("Level %u: orbits=%u, index=%u/%u, "
			   "all_same_level=%u", 0, 
			   p.level,
			   first_path_orbits.nof_orbits(),
			   index, cell->length,
			   all_same_level);
	  }
	  /* Backtrack to the previous level */
	  p.level--;
	  continue;
	}

      /* Split on smallest */
      info.split_element = next_split_element;
      
      /*
       * Save the current search situation
       */
      search_stack.push_back(info);

      /*
       * No more in the first path
       */
      info.in_first_path = false;
      /*
       * No more in the best path
       */
      info.in_best_path = false;

      p.level++;
      stats.nof_nodes++;
      if(p.level > stats.max_level)
	stats.max_level = p.level;

      p.consistency_check();

      /*
       * Move the split element to be the last in the cell
       */
      *next_split_element_pos = p.elements[cell->first + cell->length - 1];
      p.in_pos[*next_split_element_pos] = next_split_element_pos;
      p.elements[cell->first + cell->length - 1] = next_split_element;
      p.in_pos[next_split_element] = p.elements+ cell->first + cell->length -1;
      /*
       * Split the cell in two:
       * the last element in the cell (split element) forms a singleton cell
       */
      Cell * const new_cell = p.aux_split_in_two(cell, cell->length - 1);
      p.element_to_cell_map[p.elements[new_cell->first]] = new_cell;
      p.consistency_check();

      /*      
      const bool prev_equal_to_first_path = info.equal_to_first_path;
      const int prev_cmp_to_best_path = info.cmp_to_best_path;
      */
      //assert(!(!info.equal_to_first_path && info.cmp_to_best_path < 0));

      if(!first_path_info.empty())
	{
	  refine_equal_to_first = info.equal_to_first_path;
	  if(refine_equal_to_first)
	    refine_first_path_subcertificate_end =
	      first_path_info[p.level-1].certificate_index +
	      first_path_info[p.level-1].subcertificate_length;
	  if(canonical)
	    {
	      refine_cmp_to_best = info.cmp_to_best_path;
	      if(refine_cmp_to_best == 0)
		refine_best_path_subcertificate_end =
		  best_path_info[p.level-1].certificate_index +
		  best_path_info[p.level-1].subcertificate_length;
	    }
	  else
	    refine_cmp_to_best = -1;
	}
      /*
       * Refine the new partition to equitable
       */
      if(cell->length == 1)
	refine_to_equitable(cell, new_cell);
      else 
	refine_to_equitable(new_cell);

      p.consistency_check();


#ifdef PRINT_SEARCH_TREE_DOT
      fprintf(dotty_output, "\"n");
      for(unsigned int i = 0; i < search_stack.size(); i++) {
	fprintf(dotty_output, "%u", search_stack[i].split_element);
	if(i < search_stack.size() - 1) fprintf(dotty_output, ".");
      }
      fprintf(dotty_output, "\"");
      fprintf(dotty_output, " [label=\"");
      fprintf(dotty_output, "%u",cell->first); /*p.print(dotty_output);*/
      fprintf(dotty_output, "\"]");
      if(!first_path_info.empty() && canonical && refine_cmp_to_best > 0) {
	fprintf(dotty_output, "[color=green]");
      }
      fprintf(dotty_output, ";\n");
      
      fprintf(dotty_output, "\"n");
      for(unsigned int i = 0; i < search_stack.size() - 1; i++) {
	fprintf(dotty_output, "%u", search_stack[i].split_element);
	if(i < search_stack.size() - 2) fprintf(dotty_output, ".");
      }
      fprintf(dotty_output, "\" -> \"n");
      for(unsigned int i = 0; i < search_stack.size(); i++) {
	fprintf(dotty_output, "%u", search_stack[i].split_element);
	if(i < search_stack.size() - 1) fprintf(dotty_output, ".");
      }
      fprintf(dotty_output, "\" [label=\"%d\"];\n", next_split_element);
#endif

      /*
      if(prev_cmp_to_best_path == 0 && refine_cmp_to_best < 0)
	fprintf(stderr, "BP- ");
      if(prev_cmp_to_best_path == 0 && refine_cmp_to_best > 0)
	fprintf(stderr, "BP+ ");
      */

      if(p.is_discrete())
	{
	  /* Update statistics */
	  stats.nof_leaf_nodes++;
	  /*
	    if(stats.nof_leaf_nodes % 100 == 0) {
	    fprintf(stdout, "Nodes: %lu, Leafs: %lu, Bad: %lu\n",
	    stats.nof_nodes, stats.nof_leaf_nodes,
	    stats.nof_bad_nodes);
	    fflush(stdout);
	    }
	  */
	}

      if(!first_path_info.empty())
	{
	  /* We are no longer on the first path */
	  assert(best_path_info.size() > 0);
	  assert(certificate_current_path.size() >= certificate_index);
	  const unsigned int subcertificate_length = 
	    certificate_current_path.size() - certificate_index;
	  if(refine_equal_to_first)
	    {
	      /* Was equal to the first path so far */
	      assert(first_path_info.size() >= p.level);
	      PathInfo &first_pinfo = first_path_info[p.level-1];
	      assert(first_pinfo.certificate_index == certificate_index);
	      if(subcertificate_length != first_pinfo.subcertificate_length)
		{
		  refine_equal_to_first = false;
		}
	      else if(first_pinfo.eqref_hash.cmp(eqref_hash) != 0)
		{
		  refine_equal_to_first = false;
		}
	    }
	  if(canonical && (refine_cmp_to_best == 0))
	    {
	      /* Was equal to the best path so far */
	      assert(best_path_info.size() >= p.level);
	      PathInfo &best_pinfo = best_path_info[p.level-1];
	      assert(best_pinfo.certificate_index == certificate_index);
	      if(subcertificate_length < best_pinfo.subcertificate_length)
		{
		  refine_cmp_to_best = -1;
		  //fprintf(stderr, "BSCL- ");
		}
	      else if(subcertificate_length > best_pinfo.subcertificate_length)
		{
		  refine_cmp_to_best = 1;
		  //fprintf(stderr, "BSCL+ ");
		}
	      else if(best_pinfo.eqref_hash.cmp(eqref_hash) > 0)
		{
		  refine_cmp_to_best = -1;
		  //fprintf(stderr, "BHL- ");
		}
	      else if(best_pinfo.eqref_hash.cmp(eqref_hash) < 0)
		{
		  refine_cmp_to_best = 1;
		  //fprintf(stderr, "BHL+ ");
		}
	    }
	  if(refine_equal_to_first == false &&
	     (!canonical || (refine_cmp_to_best < 0)))
	    {
	      /* Backtrack */
#ifdef PRINT_SEARCH_TREE_DOT
	      fprintf(dotty_output, "\"n");
	      for(unsigned int i = 0; i < search_stack.size(); i++) {
		fprintf(dotty_output, "%u", search_stack[i].split_element);
		if(i < search_stack.size() - 1) fprintf(dotty_output, ".");
	      }
	      fprintf(dotty_output, "\" [color=red];\n");
#endif
	      stats.nof_bad_nodes++;
	      if(search_stack.back().equal_to_first_path == true &&
		 p.level > all_same_level)
		{
		  assert(all_same_level >= 1);
		  for(unsigned int i = all_same_level;
		      i < search_stack.size();
		      i++)
		    {
		      search_stack[i].equal_to_first_path = false;
		    }
		}
	      while(!search_stack.empty())
		{
		  p.level--;
		  LevelInfo &info2 = search_stack.back();
		  if(!(info2.equal_to_first_path == false &&
		       (!canonical || (info2.cmp_to_best_path < 0))))
		    break;
		  search_stack.pop_back();
		}
	      continue;
	    }
	}

#if defined(VERIFY_EQUITABLEDNESS)
      /* The new partition should be equitable */
      assert(is_equitable());
#endif

      info.equal_to_first_path = refine_equal_to_first;
      info.cmp_to_best_path = refine_cmp_to_best;

      certificate_index = certificate_current_path.size();

      search_stack.back().eqref_hash = eqref_hash;
      search_stack.back().subcertificate_length =
	certificate_index - info.certificate_index;


      if(!p.is_discrete())
	{
	  /*
	   * An internal, non-leaf node
	   */
	  /* Build the next node info */
	  /* Find the next cell to be splitted */
	  assert(cell == p.element_to_cell_map[p.elements[info.split_cell_first]]);
	  Cell * const next_split_cell = find_next_cell_to_be_splitted(cell);
	  assert(next_split_cell);
	  /* Copy current info to the search stack */
	  search_stack.push_back(info);
	  LevelInfo &new_info = search_stack.back();
	  new_info.split_cell_first = next_split_cell->first;
	  new_info.split_element = -1;
	  new_info.certificate_index = certificate_index;
	  new_info.refinement_stack_size = p.refinement_stack.size();
	  new_info.long_prune_redundant.clear();
	  new_info.long_prune_begin = info.long_prune_begin;
	  continue;
	}

      /*
       * A leaf node
       */
      assert(certificate_index == certificate_size);

      if(first_path_info.empty())
	{
	  /* The first path, update first_path and best_path */
	  //fprintf(stdout, "Level %u: FIRST\n", p.level); fflush(stdout);
	  stats.nof_canupdates++;
	  /*
	   * Update labelings and their inverses
	   */
	  update_labeling_and_its_inverse(first_path_labeling,
					  first_path_labeling_inv);
	  update_labeling_and_its_inverse(best_path_labeling,
					  best_path_labeling_inv);
	  /*
	   * Reset automorphism array
	   */
	  reset_permutation(first_path_automorphism);
	  reset_permutation(best_path_automorphism);
	  /*
	   * Reset orbit information
	   */
	  first_path_orbits.reset();
	  best_path_orbits.reset();
	  /*
	   * Reset group size
	   */
	  stats.group_size.assign(1);
	  /*
	   * Reset all_same_level
	   */
	  all_same_level = p.level;
	  /*
	   * Mark the current path to be the first and best one and save it
	   */
	  const unsigned int base_size = search_stack.size();
	  assert(p.level == base_size);
	  best_path_info.clear();
	  //fprintf(stdout, " New base is: ");
	  for(unsigned int i = 0; i < base_size; i++) {
	    search_stack[i].in_first_path = true;
	    search_stack[i].in_best_path = true;
	    search_stack[i].equal_to_first_path = true;
	    search_stack[i].cmp_to_best_path = 0;
	    PathInfo path_info;
	    path_info.splitting_element = search_stack[i].split_element;
	    path_info.certificate_index = search_stack[i].certificate_index;
	    path_info.eqref_hash = search_stack[i].eqref_hash;
	    path_info.subcertificate_length = search_stack[i].subcertificate_length;
	    first_path_info.push_back(path_info);
	    best_path_info.push_back(path_info);
	    //fprintf(stdout, "%u ", search_stack[i].split_element);
	  }
	  //fprintf(stdout, "\n"); fflush(stdout);
	  certificate_first_path = certificate_current_path;
	  certificate_best_path = certificate_current_path;

	  refine_compare_certificate = true;
	  /*
	   * Backtrack to the previous level
	   */
	  p.level--;
	  continue;
	}

      DEBUG_ASSERT(first_path_info.size() > 0);

      //fprintf(stdout, "Level %u: LEAF %d %d\n", p.level, info.equal_to_first_path, info.cmp_to_best_path); fflush(stdout);

      if(info.equal_to_first_path)
	{
	  /*
	   * An automorphism found: aut[i] = elements[first_path_labeling[i]]
	   */
	  assert(!info.in_first_path);
	  //fprintf(stdout, "A"); fflush(stdout);
	  
#ifdef PRINT_SEARCH_TREE_DOT
	  fprintf(dotty_output, "\"n");
	  for(unsigned int i = 0; i < search_stack.size(); i++) {
	    fprintf(dotty_output, "%u", search_stack[i].split_element);
	    if(i < search_stack.size() - 1) fprintf(dotty_output, ".");
	  }
	  fprintf(dotty_output, "\" [color=blue];\n");
#endif
	  
#if defined(DEBUG)
	  /* Verify that the automorphism is correctly built */
	  for(unsigned int i = 0; i < N; i++)
	    assert(first_path_automorphism[i] ==
		   p.elements[first_path_labeling[i]]);
#endif
	  
#if defined(VERIFY_AUTOMORPHISMS)
	  /* Verify that it really is an automorphism */
	  assert(is_automorphism(first_path_automorphism));
#endif

	  long_prune_add_automorphism(first_path_automorphism);

	  /*
	   * Update orbit information
	   */
	  update_orbit_information(first_path_orbits, first_path_automorphism);
	  
	  /*
	   * Compute backjumping level
	   */
	  unsigned int backjumping_level = 0;
	  for(unsigned int i = search_stack.size(); i > 0; i--) {
	    const unsigned int split_element =
	      search_stack[backjumping_level].split_element;
	    if(first_path_automorphism[split_element] != split_element)
	      break;
	    backjumping_level++;
	  }
	  assert(backjumping_level < p.level);
	  /*
	   * Go back to backjumping_level
	   */
	  p.level = backjumping_level;
	  search_stack.resize(p.level + 1);
	  
	  // if(write_automorphisms)
	  //   {
	  //     print_permutation(stdout, first_path_automorphism);
	  //     fprintf(stdout, "\n");
	  //   }
	  stats.nof_generators++;
	  continue;
	}

      assert(canonical);
      assert(info.cmp_to_best_path >= 0);
      if(info.cmp_to_best_path > 0)
	{
	  /*
	   * A new, better representative found
	   */
	  //fprintf(stdout, "Level %u: NEW BEST\n", p.level); fflush(stdout);
	  stats.nof_canupdates++;
	  /*
	   * Update canonical labeling and its inverse
	   */
	  update_labeling_and_its_inverse(best_path_labeling,
					  best_path_labeling_inv);
	  /* Reset best path automorphism */
	  reset_permutation(best_path_automorphism);
	  /* Reset best path orbit structure */
	  best_path_orbits.reset();
	  /*
	   * Mark the current path to be the best one and save it
	   */
	  const unsigned int base_size = search_stack.size();
	  assert(p.level == base_size);
	  best_path_info.clear();
	  //fprintf(stdout, " New base is: ");
	  for(unsigned int i = 0; i < base_size; i++) {
	    search_stack[i].cmp_to_best_path = 0;
	    search_stack[i].in_best_path = true;
	    PathInfo path_info;
	    path_info.splitting_element = search_stack[i].split_element;
	    path_info.certificate_index = search_stack[i].certificate_index;
	    path_info.eqref_hash = search_stack[i].eqref_hash;
	    path_info.subcertificate_length = search_stack[i].subcertificate_length;
	    best_path_info.push_back(path_info);
	    //fprintf(stdout, "%u ", search_stack[i].split_element);
	  }
	  certificate_best_path = certificate_current_path;
	  //fprintf(stdout, "\n"); fflush(stdout);
	  /*
	   * Backtrack to the previous level
	   */
	  p.level--;
	  continue;
	}

      {
	//fprintf(stderr, "BAUT ");
	/*
	 * Equal to the previous best path
	 */
#if defined(DEBUG)
	/* Verify that the automorphism is correctly built */
	for(unsigned int i = 0; i < N; i++)
	  assert(best_path_automorphism[i] ==
		 p.elements[best_path_labeling[i]]);
#endif
	
#if defined(VERIFY_AUTOMORPHISMS)
	/* Verify that it really is an automorphism */
	assert(is_automorphism(best_path_automorphism));
#endif
      
	unsigned int gca_level_with_first = 0;
	for(unsigned int i = search_stack.size(); i > 0; i--) {
	  if((int)first_path_info[gca_level_with_first].splitting_element !=
	     search_stack[gca_level_with_first].split_element)
	    break;
	  gca_level_with_first++;
	}
	assert(gca_level_with_first < p.level);

	unsigned int gca_level_with_best = 0;
	for(unsigned int i = search_stack.size(); i > 0; i--) {
	  if((int)best_path_info[gca_level_with_best].splitting_element !=
	     search_stack[gca_level_with_best].split_element)
	    break;
	  gca_level_with_best++;
	}
	assert(gca_level_with_best < p.level);

	long_prune_add_automorphism(best_path_automorphism);
	    
	/*
	 * Update orbit information
	 */
	update_orbit_information(best_path_orbits, best_path_automorphism);

	/*
	 * Update orbit information
	 */
	const unsigned int nof_old_orbits = first_path_orbits.nof_orbits();
	update_orbit_information(first_path_orbits, best_path_automorphism);
	if(nof_old_orbits != first_path_orbits.nof_orbits())
	  {
	    // if(write_automorphisms)
	    //   {
	    // 	print_permutation(stdout, best_path_automorphism);
	    // 	fprintf(stdout, "\n");
	    //   }
	    stats.nof_generators++;
	  }
	  
	/*
	 * Compute backjumping level
	 */
	unsigned int backjumping_level = p.level - 1;
	if(!first_path_orbits.is_minimal_representative(search_stack[gca_level_with_first].split_element))
	  {
	    backjumping_level = gca_level_with_first;
	    /*fprintf(stderr, "bj1: %u %u\n", p.level, backjumping_level);*/
	  }
	else
	  {
	    assert(!best_path_orbits.is_minimal_representative(search_stack[gca_level_with_best].split_element));
	    backjumping_level = gca_level_with_best;
	    /*fprintf(stderr, "bj2: %u %u\n", p.level, backjumping_level);*/
	  }
	/* Backtrack */
	search_stack.resize(backjumping_level + 1);
	p.level = backjumping_level;
	continue;
      }
    }

#ifdef PRINT_SEARCH_TREE_DOT
  fprintf(dotty_output, "}\n");
  fclose(dotty_output);
#endif
}




void AbstractGraph::find_automorphisms(Stats &stats)
{
  search(false, stats);

  if(first_path_labeling)
    {
      free(first_path_labeling);
      first_path_labeling = 0;
    }
  if(best_path_labeling)
    {
      free(best_path_labeling);
      best_path_labeling = 0;
    }
}


const unsigned int *AbstractGraph::canonical_form(Stats &stats)
{
  search(true, stats);

  return best_path_labeling;
}




/*-------------------------------------------------------------------------
 *
 * Routines for undirected graphs
 *
 *-------------------------------------------------------------------------*/

Graph::Vertex::Vertex()
{
  label = 1;
  nof_edges = 0;
}


Graph::Vertex::~Vertex()
{
  ;
}


void Graph::Vertex::add_edge(const unsigned int other_vertex)
{
  edges.push_back(other_vertex);
  nof_edges++;
  DEBUG_ASSERT(nof_edges == edges.size());
}


void Graph::Vertex::remove_duplicate_edges(bool * const duplicate_array)
{
  for(std::vector<unsigned int>::iterator iter = edges.begin();
      iter != edges.end(); )
    {
      const unsigned int dest_vertex = *iter;
      if(duplicate_array[dest_vertex] == true)
	{
	  /* A duplicate edge found! */
	  iter = edges.erase(iter);
	  nof_edges--;
	  DEBUG_ASSERT(nof_edges == edges.size());
	}
      else
	{
	  /* Not seen earlier, mark as seen */
	  duplicate_array[dest_vertex] = true;
	  iter++;
	}
    }

  /* Clear duplicate_array */
  for(std::vector<unsigned int>::iterator iter = edges.begin();
      iter != edges.end();
      iter++)
    {
      duplicate_array[*iter] = false;
    }
}





/*-------------------------------------------------------------------------
 *
 * Constructor and destructor for undirected graphs
 *
 *-------------------------------------------------------------------------*/


Graph::Graph(const unsigned int nof_vertices)
{
  vertices.resize(nof_vertices);
  sh = sh_flm;
}


Graph::~Graph()
{
  ;
}


unsigned int Graph::add_vertex(const unsigned int new_label)
{
  const unsigned int new_vertex_num = vertices.size();
  vertices.resize(new_vertex_num + 1);
  vertices.back().label = new_label;
  return new_vertex_num;
}


void Graph::add_edge(const unsigned int vertex1, const unsigned int vertex2)
{
  //fprintf(stderr, "(%u,%u) ", vertex1, vertex2);
  assert(vertex1 < vertices.size());
  assert(vertex2 < vertices.size());
  vertices[vertex1].add_edge(vertex2);
  vertices[vertex2].add_edge(vertex1);
}


void Graph::change_label(const unsigned int vertex,
			   const unsigned int new_label)
{
  assert(vertex < vertices.size());
  vertices[vertex].label = new_label;
}





/*-------------------------------------------------------------------------
 *
 * Read graph in the DIMACS format
 *
 *-------------------------------------------------------------------------*/

// Graph *Graph::read_dimacs(FILE *fp)
// {
//   Graph *g = 0;
//   unsigned int nof_vertices, nof_edges;
//   unsigned int line_num = 1;
//   int c;
  
//   /* read comments and problem line*/
//   while(1) {
//     c = getc(fp);
//     if(c == 'c') {
//       while((c = getc(fp)) != '\n') {
//         if(c == EOF) {
//           fprintf(stderr, "error in line %u: not in DIMACS format\n",
//                   line_num);
//           goto error_exit;
//         }
//       }
//       line_num++;
//       continue;
//     }
//     if(c == 'p') {
//       if(fscanf(fp, " edge %u %u\n", &nof_vertices, &nof_edges) != 2) {
//         fprintf(stderr, "error in line %u: not in DIMACS format\n",
//                 line_num);
//         goto error_exit; }
//       line_num++;
//       break;
//     }
//     fprintf(stderr, "error in line %u: not in DIMACS format\n", line_num);
//     goto error_exit;
//   }
  
//   if(nof_vertices <= 0) {
//     fprintf(stderr, "error: no vertices\n");
//     goto error_exit;
//   }
// #if 0
//   if(nof_edges <= 0) {
//     fprintf(stderr, "error: no edges\n");
//     goto error_exit;
//   }
// #endif
//   if(bliss_verbose) {
//     fprintf(bliss_verbstr, "Instance has %d vertices and %d edges\n",
//             nof_vertices, nof_edges);
//     fflush(bliss_verbstr);
//   }

//   g = new Graph(nof_vertices);

//   //
//   // Read vertex labels
//   //
//   if(bliss_verbose) {
//     fprintf(bliss_verbstr, "Reading vertex labels...\n");
//     fflush(bliss_verbstr); }
//   while(1) {
//     c = getc(fp);
//     if(c != 'n') {
//       ungetc(c, fp);
//       break;
//     }
//     ungetc(c, fp);
//     unsigned int vertex, label;
//     if(fscanf(fp, "n %u %u\n", &vertex, &label) != 2) {
//       fprintf(stderr, "error in line %u: not in DIMACS format\n",
// 	      line_num);
//       goto error_exit;
//     }
//     if(vertex > nof_vertices) {
//       fprintf(stderr, "error in line %u: not in DIMACS format\n",
// 	      line_num);
//       goto error_exit;
//     }
//     line_num++;
//     g->change_label(vertex - 1, label);
//   }
//   if(bliss_verbose) {
//     fprintf(bliss_verbstr, "Done\n");
//     fflush(bliss_verbstr); }

//   //
//   // Read edges
//   //
//   if(bliss_verbose) {
//     fprintf(bliss_verbstr, "Reading edges...\n");
//     fflush(bliss_verbstr); }
//   for(unsigned i = 0; i < nof_edges; i++) {
//     unsigned int from, to;
//     if(fscanf(fp, "e %u %u\n", &from, &to) != 2) {
//       fprintf(stderr, "error in line %u: not in DIMACS format\n",
// 	      line_num);
//       goto error_exit;
//     }
//     if(from > nof_vertices || to > nof_vertices) {
//       fprintf(stderr, "error in line %u: not in DIMACS format\n",
// 	      line_num);
//       goto error_exit;
//     }
//     line_num++;
//     g->add_edge(from - 1, to - 1);
//   }
//   if(bliss_verbose) {
//     fprintf(bliss_verbstr, "Done\n");
//     fflush(bliss_verbstr);
//   }

//   return g;

//  error_exit:
//   if(g)
//     delete g;
//   return 0;

// }

Graph *Graph::from_igraph(const igraph_t *graph) {
  
  unsigned int nof_vertices= (unsigned int)igraph_vcount(graph);
  unsigned int nof_edges= (unsigned int)igraph_ecount(graph);
  Graph *g=new Graph(nof_vertices);
//   for (unsigned int i=0; i<nof_vertices; i++) {
//     g->change_label(i, i);
//   }
  for (unsigned int i=0; i<nof_edges; i++) {
    g->add_edge((unsigned int)IGRAPH_FROM(graph, i), 
		(unsigned int)IGRAPH_TO(graph, i));
  }  
  return g;
}

void Graph::print_dimacs(FILE *fp)
{
  unsigned int nof_edges = 0;
  for(unsigned int i = 0; i < get_nof_vertices(); i++)
    {
      Vertex &v = vertices[i];
      for(std::vector<unsigned int>::const_iterator ei = v.edges.begin();
	  ei != v.edges.end();
	  ei++)
	{
	  const unsigned int dest_i = *ei;
	  if(dest_i < i)
	    continue;
	  nof_edges++;
	}
    }

  fprintf(fp, "p edge %u %u\n", get_nof_vertices(), nof_edges);
  for(unsigned int i = 0; i < get_nof_vertices(); i++)
    {
      Vertex &v = vertices[i];
      if(v.label != 1)
	{
	  fprintf(fp, "n %u %u\n", i+1, v.label);
	}
    }
  for(unsigned int i = 0; i < get_nof_vertices(); i++)
    {
      Vertex &v = vertices[i];
      for(std::vector<unsigned int>::const_iterator ei = v.edges.begin();
	  ei != v.edges.end();
	  ei++)
	{
	  const unsigned int dest_i = *ei;
	  if(dest_i < i)
	    continue;
	  fprintf(fp, "e %u %u\n", i+1, dest_i+1);
	}
    }
}




Graph *Graph::permute(const unsigned int *perm)
{
  Graph *g = new Graph(get_nof_vertices());
  for(unsigned int i = 0; i < get_nof_vertices(); i++)
    {
      Vertex &v = vertices[i];
      Vertex &permuted_v = g->vertices[perm[i]];
      permuted_v.label = v.label;
      for(std::vector<unsigned int>::const_iterator ei = v.edges.begin();
	  ei != v.edges.end();
	  ei++)
	{
	  const unsigned int dest_v = *ei;
	  permuted_v.add_edge(perm[dest_v]);
	}
      std::sort(permuted_v.edges.begin(), permuted_v.edges.end());
    }
  return g;
}





/*-------------------------------------------------------------------------
 *
 * Print graph in graphviz format
 *
 *-------------------------------------------------------------------------*/


void Graph::to_dot(const char *file_name)
{
  FILE *fp = fopen(file_name, "w");
  if(fp)
    to_dot(fp);
  fclose(fp);
}


void Graph::to_dot(FILE *fp)
{
  remove_duplicate_edges();

  fprintf(fp, "graph g {\n");

  unsigned int vnum = 0;
  for(std::vector<Vertex>::iterator vi = vertices.begin();
      vi != vertices.end();
      vi++, vnum++)
    {
      Vertex &v = *vi;
      fprintf(fp, "v%u [label=\"%u:%u\"];\n", vnum, vnum, v.label);
      for(std::vector<unsigned int>::const_iterator ei = v.edges.begin();
	  ei != v.edges.end();
	  ei++)
	{
	  const unsigned int vnum2 = *ei;
	  if(vnum2 > vnum)
	    fprintf(fp, "v%u -- v%u\n", vnum, vnum2);
	}
    }

  fprintf(fp, "}\n");
}





void Graph::remove_duplicate_edges()
{
  bool *duplicate_array = (bool*)calloc(vertices.size(), sizeof(bool));

  for(std::vector<Vertex>::iterator vi = vertices.begin();
      vi != vertices.end();
      vi++)
    {
#ifdef EXPENSIVE_CONSISTENCY_CHECKS
      for(unsigned int i = 0; i < vertices.size(); i++)
	assert(duplicate_array[i] == false);
#endif
      Vertex &v = *vi;
      v.remove_duplicate_edges(duplicate_array);
    }

  free(duplicate_array);
}





/*-------------------------------------------------------------------------
 *
 * Partition independent invariants
 *
 *-------------------------------------------------------------------------*/


unsigned int Graph::label_invariant(Graph *g, unsigned int v)
{
  DEBUG_ASSERT(v < g->vertices.size());
  return g->vertices[v].label;
}


unsigned int Graph::degree_invariant(Graph *g, unsigned int v)
{
  DEBUG_ASSERT(v < g->vertices.size());
  DEBUG_ASSERT(g->vertices[v].edges.size() ==
	       g->vertices[v].nof_edges);
  return g->vertices[v].nof_edges;
}







/*-------------------------------------------------------------------------
 *
 * Refine the partition p according to a partition independent invariant
 *
 *-------------------------------------------------------------------------*/

bool Graph::refine_according_to_invariant(unsigned int (*inv)(Graph * const g, unsigned int v))
{
  bool refined = false;

  for(Cell *cell = p.first_cell; cell; )
    {
      assert(cell->max_ival == 0);
      assert(cell->max_ival_count == 0);
      
      Cell * const next_cell = cell->next;

      if(cell->length == 1)
	{
	  cell = next_cell;
	  continue;
	}
      
      const unsigned int *ep = p.elements + cell->first;
      for(unsigned int i = cell->length; i > 0; i--, ep++)
	{
	  unsigned int ival = inv(this, *ep);
	  p.invariant_values[*ep] = ival;
	  if(ival > cell->max_ival) {
	    cell->max_ival = ival;
	    cell->max_ival_count = 1;
	  }
	  else if(ival == cell->max_ival) {
	    cell->max_ival_count++;
	  }
	}
      Cell * const last_new_cell = p.zplit_cell(cell, true);
      refined = (last_new_cell != cell);
      cell = next_cell;
    }

  return refined;
}





/*-------------------------------------------------------------------------
 *
 * Split the neighbourhood of a cell according to the equitable invariant
 *
 *-------------------------------------------------------------------------*/


void Graph::split_neighbourhood_of_cell(Cell * const cell)
{
  DEBUG_ASSERT(neighbour_heap.is_empty());
  DEBUG_ASSERT(cell->length > 1);

  eqref_hash.update(cell->first);
  eqref_hash.update(cell->length);

  unsigned int *ep = p.elements + cell->first;
  for(unsigned int i = cell->length; i > 0; i--)
    {
      const Vertex &v = vertices[*ep];
      ep++;
      
      std::vector<unsigned int>::const_iterator ei = v.edges.begin();
      for(unsigned int j = v.nof_edges; j > 0; j--)
	{
	  const unsigned int dest_vertex = *ei++;
	  Cell * const neighbour_cell = p.element_to_cell_map[dest_vertex];
	  if(neighbour_cell->length == 1)
	    continue;
	  const unsigned int ival = p.invariant_values[dest_vertex] + 1;
	  p.invariant_values[dest_vertex] = ival;
	  if(ival > neighbour_cell->max_ival) {
	    neighbour_cell->max_ival = ival;
	    neighbour_cell->max_ival_count = 1;
	  }
	  else if(ival == neighbour_cell->max_ival) {
	    neighbour_cell->max_ival_count++;
	  }
	  if(!neighbour_cell->in_neighbour_heap) {
	    neighbour_cell->in_neighbour_heap = true;
	    neighbour_heap.insert(neighbour_cell->first);
	  }
	}
    }

  while(!neighbour_heap.is_empty())
    {
      const unsigned int start = neighbour_heap.remove();
      Cell * const neighbour_cell = p.element_to_cell_map[p.elements[start]];
      DEBUG_ASSERT(neighbour_cell->first == start);
      DEBUG_ASSERT(neighbour_cell->in_neighbour_heap);
      neighbour_cell->in_neighbour_heap = false;
      
      DEBUG_ASSERT(neighbour_cell->length > 1);
      DEBUG_ASSERT(neighbour_cell->max_ival >= 1);
      DEBUG_ASSERT(neighbour_cell->max_ival_count >= 1);
      
      eqref_hash.update(neighbour_cell->first);
      eqref_hash.update(neighbour_cell->length);
      eqref_hash.update(neighbour_cell->max_ival);
      eqref_hash.update(neighbour_cell->max_ival_count);

      Cell * const last_new_cell = p.zplit_cell(neighbour_cell, true);
      /* Update hash */
      const Cell *c = neighbour_cell;
      while(1)
	{
	  eqref_hash.update(c->first);
	  eqref_hash.update(c->length);
	  if(c == last_new_cell)
	    break;
	  c = c->next;
	}
    }
}


bool Graph::split_neighbourhood_of_unit_cell(Cell * const unit_cell)
{
  DEBUG_ASSERT(neighbour_heap.is_empty());

  DEBUG_ASSERT(unit_cell->length == 1);
  DEBUG_ASSERT(p.element_to_cell_map[p.elements[unit_cell->first]] == unit_cell);
  DEBUG_ASSERT(p.in_pos[p.elements[unit_cell->first]] ==
               p.elements + unit_cell->first);

  eqref_hash.update(0x87654321);
  eqref_hash.update(unit_cell->first);
  eqref_hash.update(1);

  const Vertex &v = vertices[p.elements[unit_cell->first]];

  std::vector<unsigned int>::const_iterator ei = v.edges.begin();
  for(unsigned int j = v.nof_edges; j > 0; j--)
    {
      const unsigned int dest_vertex = *ei++;
      Cell * const neighbour_cell = p.element_to_cell_map[dest_vertex];
      DEBUG_ASSERT(*p.in_pos[dest_vertex] == dest_vertex);
      
      if(neighbour_cell->length == 1) {
	DEBUG_ASSERT(!neighbour_cell->in_neighbour_heap);
	if(in_search) {
	  neighbour_cell->in_neighbour_heap = true;
	  neighbour_heap.insert(neighbour_cell->first);
	}
	continue;
      }
      if(!neighbour_cell->in_neighbour_heap) {
	neighbour_cell->in_neighbour_heap = true;
	neighbour_heap.insert(neighbour_cell->first);
      }
      neighbour_cell->max_ival_count++;
      DEBUG_ASSERT(neighbour_cell->max_ival_count <= neighbour_cell->length);
      
      unsigned int * const swap_position =
	p.elements + neighbour_cell->first + neighbour_cell->length -
	neighbour_cell->max_ival_count;
      DEBUG_ASSERT(p.in_pos[dest_vertex] <= swap_position);
      *p.in_pos[dest_vertex] = *swap_position;
      p.in_pos[*swap_position] = p.in_pos[dest_vertex];
      *swap_position = dest_vertex;
      p.in_pos[dest_vertex] = swap_position;
    }

  while(!neighbour_heap.is_empty())
    {
      const unsigned int start = neighbour_heap.remove();
      Cell *neighbour_cell = p.element_to_cell_map[p.elements[start]];
      DEBUG_ASSERT(neighbour_cell->in_neighbour_heap);
      neighbour_cell->in_neighbour_heap = false;

#ifdef DEBUG
      assert(neighbour_cell->first == start);
      if(neighbour_cell->length == 1) {
	assert(neighbour_cell->max_ival_count == 0);
      } else {
	assert(neighbour_cell->max_ival_count > 0);
	assert(neighbour_cell->max_ival_count <= neighbour_cell->length);
      }
#endif

      eqref_hash.update(neighbour_cell->first);
      eqref_hash.update(neighbour_cell->length);
      eqref_hash.update(neighbour_cell->max_ival_count);

      if(neighbour_cell->length > 1 &&
	 neighbour_cell->max_ival_count != neighbour_cell->length) {

	p.consistency_check();

	Cell * const new_cell = p.aux_split_in_two(neighbour_cell, neighbour_cell->length - neighbour_cell->max_ival_count);
	unsigned int *ep = p.elements + new_cell->first;
	unsigned int * const lp = p.elements+new_cell->first+new_cell->length;
	while(ep < lp) {
	  DEBUG_ASSERT(p.in_pos[*ep] == ep);
	  p.element_to_cell_map[*ep] = new_cell;
	  ep++;
	}
	neighbour_cell->max_ival_count = 0;

	p.consistency_check();

	/* update hash */
	eqref_hash.update(neighbour_cell->first);
	eqref_hash.update(neighbour_cell->length);
	eqref_hash.update(0);
	eqref_hash.update(new_cell->first);
	eqref_hash.update(new_cell->length);
	eqref_hash.update(1);

	/* Add cells in splitting_queue */
	DEBUG_ASSERT(!new_cell->in_splitting_queue);
	if(neighbour_cell->in_splitting_queue) {
	  /* Both cells must be included in splitting_queue in order
	     to have refinement to equitable partition */
	  p.add_in_splitting_queue(new_cell);
	} else {
	  Cell *min_cell, *max_cell;
	  if(neighbour_cell->length <= new_cell->length) {
	    min_cell = neighbour_cell;
	    max_cell = new_cell;
	  } else {
	    min_cell = new_cell;
	    max_cell = neighbour_cell;
	  }
	  /* Put the smaller cell in splitting_queue */
	  p.add_in_splitting_queue(min_cell);
	  if(max_cell->length == 1) {
	    /* Put the "larger" cell also in splitting_queue */
	    p.add_in_splitting_queue(max_cell);
	  }
	}
	/* Update pointer for certificate generation */
	neighbour_cell = new_cell;
      }
      else
	neighbour_cell->max_ival_count = 0;

      /*
       * Build certificate if required
       */
      if(in_search)
	{
	  for(unsigned int i = neighbour_cell->first,
		j = neighbour_cell->length,
		c_index = certificate_current_path.size();
	      j > 0;
	      j--, i++, c_index += 2)
	    {
	      if(refine_compare_certificate)
		{
		  if(refine_equal_to_first)
		    {
		      if(c_index >= refine_first_path_subcertificate_end)
			refine_equal_to_first = false;
		      else if(certificate_first_path[c_index] !=
			      unit_cell->first)
			refine_equal_to_first = false;
		      else if(certificate_first_path[c_index+1] != i)
			refine_equal_to_first = false;
		    }
		  if(refine_cmp_to_best == 0)
		    {
		      if(c_index >= refine_best_path_subcertificate_end)
			{
			  refine_cmp_to_best = 1;
			}
		      else if(unit_cell->first>certificate_best_path[c_index])
			{
			  refine_cmp_to_best = 1;
			}
		      else if(unit_cell->first<certificate_best_path[c_index])
			{
			  refine_cmp_to_best = -1;
			}
		      else if(i > certificate_best_path[c_index+1])
			{
			  refine_cmp_to_best = 1;
			}
		      else if(i < certificate_best_path[c_index+1])
			{
			  refine_cmp_to_best = -1;
			}
		    }
		  if((refine_equal_to_first == false) &&
		     (refine_cmp_to_best < 0))
		    goto worse_exit;
		}
	      certificate_current_path.push_back(unit_cell->first);
	      certificate_current_path.push_back(i);
	    }
	} /* if(in_search) */
    } /* while(!neighbour_heap.is_empty()) */
  
  return false;

 worse_exit:
  while(!neighbour_heap.is_empty())
    {
      const unsigned int start = neighbour_heap.remove();
      Cell * const neighbour_cell = p.element_to_cell_map[p.elements[start]];
      DEBUG_ASSERT(neighbour_cell->in_neighbour_heap);
      neighbour_cell->in_neighbour_heap = false;
      neighbour_cell->max_ival_count = 0;
    }
  return true;
}





/*-------------------------------------------------------------------------
 *
 * Check whether the current partition p is equitable
 * Slow: use only for debugging purposes
 * Side effect: resets max_ival and max_ival_count fields in cells
 *
 *-------------------------------------------------------------------------*/

bool Graph::is_equitable()
{
  bool result = true;

  /*
   * Max ival and max_ival_count are used for counting purposes,
   * they should be reset...
   */
  for(Cell *cell = p.first_cell; cell; cell = cell->next)
    {
      assert(cell->prev_next_ptr && *(cell->prev_next_ptr) == cell);
      assert(cell->max_ival == 0);
      assert(cell->max_ival_count == 0);
    }


  for(Cell *cell = p.first_cell; cell; cell = cell->next)
    {
      if(cell->length == 1)
	continue;

      unsigned int *ep = p.elements + cell->first;
      Vertex &first_vertex = vertices[*ep++];

      /* Count edges of the first vertex for cells in max_ival */
      std::vector<unsigned int>::const_iterator ei = first_vertex.edges.begin();
      for(unsigned int j = first_vertex.nof_edges; j > 0; j--)
	{
	  p.element_to_cell_map[*ei++]->max_ival++;
	}

      /* Count and compare edges of the other vertices */
      for(unsigned int i = cell->length; i > 1; i--)
	{
	  Vertex &vertex = vertices[*ep++];
	  std::vector<unsigned int>::const_iterator ei = vertex.edges.begin();
	  for(unsigned int j = vertex.nof_edges; j > 0; j--)
	    {
	      p.element_to_cell_map[*ei++]->max_ival_count++;
	    }
	  for(Cell *cell2 = p.first_cell; cell2; cell2 = cell2->next)
	    {
	      if(cell2->max_ival != cell2->max_ival_count)
		{
		  result = false;
		  goto done;
		}
	      cell2->max_ival_count = 0;
	    }
	}
      /* Reset max_ival */
      for(Cell *cell2 = p.first_cell; cell2; cell2 = cell2->next)
	{
	  cell2->max_ival = 0;
	  assert(cell2->max_ival_count == 0);
	}
    }

 done:

  for(Cell *cell = p.first_cell; cell; cell = cell->next)
    {
      cell->max_ival = 0;
      cell->max_ival_count = 0;
    }

  return result;
}





/*-------------------------------------------------------------------------
 *
 * Build the initial equitable partition
 *
 *-------------------------------------------------------------------------*/

void Graph::make_initial_equitable_partition()
{
  refine_according_to_invariant(&label_invariant);
  p.clear_splitting_queue();
  //p.print_signature(stderr); fprintf(stderr, "\n");

  refine_according_to_invariant(&degree_invariant);
  p.clear_splitting_queue();
  //p.print_signature(stderr); fprintf(stderr, "\n");

  /* To do: add loop invariant */

  refine_to_equitable();
  p.refinement_stack.clean();
  //p.print_signature(stderr); fprintf(stderr, "\n");
}





/*-------------------------------------------------------------------------
 *
 * Find the next cell to be splitted
 *
 *-------------------------------------------------------------------------*/

Cell *Graph::find_next_cell_to_be_splitted(Cell *cell)
{
  assert(!p.is_discrete());
  switch(sh) {
  case sh_f:
    return sh_first(cell);
  case sh_fs:
    return sh_first_smallest(cell);
  case sh_fl:
    return sh_first_largest(cell);
  case sh_fm:
    return sh_first_max_neighbours(cell);
  case sh_fsm:
    return sh_first_smallest_max_neighbours(cell);
  case sh_flm:
    return sh_first_largest_max_neighbours(cell);
  default:
    assert(false && "Unknown splitting heuristics");
    return 0;
  }
}

/* First nonsingleton cell */
Cell *Graph::sh_first(Cell *cell)
{
  return p.first_nonsingleton_cell;
}

/* First smallest nonsingleton cell. */
Cell *Graph::sh_first_smallest(Cell *cell)
{
  Cell *best_cell = 0;
  unsigned int best_size = UINT_MAX;
  for(cell = p.first_nonsingleton_cell; cell; cell = cell->next_nonsingleton)
    {
      assert(cell->length > 1);
      if(cell->length < best_size)
	{
	  best_size = cell->length;
	  best_cell = cell;
	}
    }
  assert(best_cell);
  return best_cell;
}

/* First largest nonsingleton cell. */
Cell *Graph::sh_first_largest(Cell *cell)
{
  Cell *best_cell = 0;
  unsigned int best_size = 0;
  for(cell = p.first_nonsingleton_cell; cell; cell = cell->next_nonsingleton)
    {
      assert(cell->length > 1);
      if(cell->length > best_size)
	{
	  best_size = cell->length;
	  best_cell = cell;
	}
    }
  assert(best_cell);
  return best_cell;
}

/* First nonsingleton cell with max number of neighbouring
 * nonsingleton cells.
 * Assumes that the partition p is equitable.
 * Messes up in_neighbour_heap and max_ival fields of cells
 * (assumes they are false/0).
 */
Cell *Graph::sh_first_max_neighbours(Cell *cell)
{
  Cell *best_cell = 0;
  int best_value = -1;
  for(cell = p.first_nonsingleton_cell; cell; cell = cell->next_nonsingleton)
    {
      assert(cell->length > 1);
	
      const Vertex &v = vertices[p.elements[cell->first]];
      std::vector<unsigned int>::const_iterator ei = v.edges.begin();
      std::list<Cell*> neighbour_cells_visited;
      for(unsigned int j = v.nof_edges; j > 0; j--)
	{
	  const unsigned int dest_vertex = *ei++;
	  Cell * const neighbour_cell = p.element_to_cell_map[dest_vertex];
	  if(neighbour_cell->length == 1)
	    continue;
	  neighbour_cell->max_ival++;
	  if(neighbour_cell->in_neighbour_heap)
	    continue;
	  neighbour_cell->in_neighbour_heap = true;
	  neighbour_cells_visited.push_back(neighbour_cell);
	}
      int value = 0;
      while(!neighbour_cells_visited.empty())
	{
	  Cell * const neighbour_cell = neighbour_cells_visited.front();
	  neighbour_cells_visited.pop_front();
	  assert(neighbour_cell->in_neighbour_heap);
	  neighbour_cell->in_neighbour_heap = false;
	  if(neighbour_cell->max_ival != neighbour_cell->length)
	    value++;
	  neighbour_cell->max_ival = 0;
	}
      if(value > best_value)
	{
	  best_value = value;
	  best_cell = cell;
	}
    }
  assert(best_cell);
  return best_cell;
}
/* First smallest nonsingleton cell with max number of neighbouring
 * nonsingleton cells.
 * Assumes that the partition p is equitable.
 * Messes up in_neighbour_heap and max_ival fields of cells
 * (assumes they are false).
 */
Cell *Graph::sh_first_smallest_max_neighbours(Cell *cell)
{
  Cell *best_cell = 0;
  int best_value = -1;
  int best_size = INT_MAX;
  for(cell = p.first_nonsingleton_cell; cell; cell = cell->next_nonsingleton)
    {
      assert(cell->length > 1);
	
      const Vertex &v = vertices[p.elements[cell->first]];
      std::vector<unsigned int>::const_iterator ei = v.edges.begin();
      std::list<Cell*> neighbour_cells_visited;
      for(unsigned int j = v.nof_edges; j > 0; j--)
	{
	  const unsigned int dest_vertex = *ei++;
	  Cell * const neighbour_cell = p.element_to_cell_map[dest_vertex];
	  if(neighbour_cell->length == 1)
	    continue;
	  neighbour_cell->max_ival++;
	  if(neighbour_cell->in_neighbour_heap)
	    continue;
	  neighbour_cell->in_neighbour_heap = true;
	  neighbour_cells_visited.push_back(neighbour_cell);
	}
      int value = 0;
      while(!neighbour_cells_visited.empty())
	{
	  Cell * const neighbour_cell = neighbour_cells_visited.front();
	  neighbour_cells_visited.pop_front();
	  assert(neighbour_cell->in_neighbour_heap);
	  neighbour_cell->in_neighbour_heap = false;
	  if(neighbour_cell->max_ival != neighbour_cell->length)
	    value++;
	  neighbour_cell->max_ival = 0;
	}
      if((value > best_value) ||
	 (value == best_value && (int)cell->length < best_size))
	{
	  best_value = value;
	  best_size = cell->length;
	  best_cell = cell;
	}
    }
  assert(best_cell);
  return best_cell;
}
/* First largest nonsingleton cell with max number of neighbouring
 * nonsingleton cells.
 * Assumes that the partition p is equitable.
 * Messes up in_neighbour_heap and max_ival fields of cells
 * (assumes they are false/0).
 */
Cell *Graph::sh_first_largest_max_neighbours(Cell *cell)
{
  Cell *best_cell = 0;
  int best_value = -1;
  int best_size = -1;
  for(cell = p.first_nonsingleton_cell; cell; cell = cell->next_nonsingleton)
    {
      assert(cell->length > 1);
	
      const Vertex &v = vertices[p.elements[cell->first]];
      std::vector<unsigned int>::const_iterator ei = v.edges.begin();
      std::list<Cell*> neighbour_cells_visited;
      for(unsigned int j = v.nof_edges; j > 0; j--)
	{
	  const unsigned int dest_vertex = *ei++;
	  Cell * const neighbour_cell = p.element_to_cell_map[dest_vertex];
	  if(neighbour_cell->length == 1)
	    continue;
	  neighbour_cell->max_ival++;
	  if(neighbour_cell->in_neighbour_heap)
	    continue;
	  neighbour_cell->in_neighbour_heap = true;
	  neighbour_cells_visited.push_back(neighbour_cell);
	}
      int value = 0;
      while(!neighbour_cells_visited.empty())
	{
	  Cell * const neighbour_cell = neighbour_cells_visited.front();
	  neighbour_cells_visited.pop_front();
	  assert(neighbour_cell->in_neighbour_heap);
	  neighbour_cell->in_neighbour_heap = false;
	  if(neighbour_cell->max_ival != neighbour_cell->length)
	    value++;
	  neighbour_cell->max_ival = 0;
	}
      if((value > best_value) ||
	 (value == best_value && (int)cell->length > best_size))
	{
	  best_value = value;
	  best_size = cell->length;
	  best_cell = cell;
	}
    }
  assert(best_cell);
  return best_cell;
}





/*-------------------------------------------------------------------------
 *
 * Initialize the certificate size and memory
 *
 *-------------------------------------------------------------------------*/

void Graph::initialize_certificate()
{
  certificate_size = 0;
  for(Cell *cell = p.first_cell; cell; cell = cell->next)
    {
      if(cell->length > 1) {
	certificate_size +=
	  vertices[p.elements[cell->first]].nof_edges * 2 * cell->length;
      }
    }
  //if(certificate)
  //  free(certificate);
  //certificate = (unsigned int*)malloc(certificate_size * sizeof(unsigned int));
  certificate_index = 0;

  certificate_current_path.clear();
  certificate_first_path.clear();
  certificate_best_path.clear();
}





/*-------------------------------------------------------------------------
 *
 * Check whether perm is an automorphism
 *
 *-------------------------------------------------------------------------*/

bool Graph::is_automorphism(unsigned int * const perm)
{
  std::set<unsigned int, std::less<unsigned int> > edges1;
  std::set<unsigned int, std::less<unsigned int> > edges2;

  bool result = true;

  for(unsigned int i = 0; i < vertices.size(); i++)
    {
      Vertex &v1 = vertices[i];
      edges1.clear();
      for(std::vector<unsigned int>::iterator ei = v1.edges.begin();
	  ei != v1.edges.end();
	  ei++)
	edges1.insert(perm[*ei]);
      
      Vertex &v2 = vertices[perm[i]];
      edges2.clear();
      for(std::vector<unsigned int>::iterator ei = v2.edges.begin();
	  ei != v2.edges.end();
	  ei++)
	edges2.insert(*ei);

      if(!(edges1 == edges2))
	{
	  result = false;
	  goto done;
	}
    }

 done:

  return result;
}

}
