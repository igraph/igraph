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

#include <cassert>
#include <vector>
#include <list>
#include "bliss_graph.hh"
#include "bliss_partition.hh"

using namespace std;

namespace igraph {

static const bool should_not_happen = false;

Partition::Partition()
{
  elements = 0;
  in_pos = 0;
  invariant_values = 0;
  cells = 0;
  free_cells = 0;
  element_to_cell_map = 0;
  graph = 0;
}

Partition::~Partition()
{
  if(elements) {
    free(elements); elements = 0; }
  if(cells) {
    free(cells); cells = 0; }
  if(element_to_cell_map) {
    free(element_to_cell_map); element_to_cell_map = 0; }
  if(in_pos) {
    free(in_pos); in_pos = 0; }
  if(invariant_values) {
    free(invariant_values); invariant_values = 0; }
}

void Partition::init(const unsigned int N)
{
  assert(N > 0);

  if(elements)
    free(elements);
  elements = (unsigned int*)malloc(N * sizeof(unsigned int));
  for(unsigned int i = 0; i < N; i++)
    elements[i] = i;

  if(in_pos)
    free(in_pos);
  in_pos = (unsigned int**)malloc(N * sizeof(unsigned int*));
  for(unsigned int i = 0; i < N; i++)
    in_pos[i] = elements + i;

  if(invariant_values)
    free(invariant_values);
  invariant_values = (unsigned int*)malloc(N * sizeof(unsigned int));
  for(unsigned int i = 0; i < N; i++)
    invariant_values[i] = 0;

  if(cells)
    free(cells);
  cells = (Cell*)malloc(N * sizeof(Cell));

  cells[0].first = 0;
  cells[0].length = N;
  cells[0].max_ival = 0;
  cells[0].max_ival_count = 0;
  cells[0].in_splitting_queue = false;
  cells[0].in_neighbour_heap = false;
  cells[0].next = 0;
  cells[0].prev_next_ptr = &first_cell;
  cells[0].next_nonsingleton = 0;
  cells[0].prev_nonsingleton = 0;
  cells[0].split_level = 0;
  first_cell = &cells[0];
  if(N == 1)
    first_nonsingleton_cell = 0;
  else
    first_nonsingleton_cell = &cells[0];

  for(unsigned int i = 1; i < N; i++)
    {
      cells[i].first = 0;
      cells[i].length = 0;
      cells[i].max_ival = 0;
      cells[i].max_ival_count = 0;
      cells[i].in_splitting_queue = false;
      cells[i].in_neighbour_heap = false;
      cells[i].next = (i < N-1)?&cells[i+1]:0;
      cells[i].prev_next_ptr = (i == 1)?&free_cells:&(cells[i-1].next);
      cells[i].next_nonsingleton = 0;
      cells[i].prev_nonsingleton = 0;
    }
  if(N > 1)
    free_cells = &cells[1];
  else
    free_cells = 0;

  if(element_to_cell_map)
    free(element_to_cell_map);
  element_to_cell_map = (Cell **)malloc(N * sizeof(Cell *));
  for(unsigned int i = 0; i < N; i++)
    element_to_cell_map[i] = first_cell;

  splitting_queue.init(N);
  refinement_stack.init(N);
  level = 0;
}

/*
 * For debugging purposes only.
 * Checks that the ordered list of nonsingleton cells is consistent.
 */
void Partition::consistency_check()
{
#ifdef DEBUG
  for(const Cell *cell = first_cell; cell; cell = cell->next)
    {
      assert(cell->prev_next_ptr && *(cell->prev_next_ptr) == cell);
    }
  const bool do_print = false;
  if(do_print)
    {
      fprintf(stderr, "\nRef stack: ");
      for(unsigned int j = 0; j < refinement_stack.size(); j++)
	{
	  const RefInfo i = refinement_stack.element_at(j);
	  fprintf(stderr, "f%u,%d,%d ",
	      i.split_cell_first,
		  i.prev_nonsingleton_first,
		  i.next_nonsingleton_first);
	}
      fprintf(stderr, "\n");
      for(const Cell *cell = first_nonsingleton_cell; cell;
	  cell = cell->next_nonsingleton)
	{
	  fprintf(stderr, "%u:%u->", cell->first, cell->length);
	  if(cell->next_nonsingleton)
	    assert(cell->first < cell->next_nonsingleton->first);
	}
      fprintf(stderr, "\n");
    }
  const Cell *next_nonsingleton = first_nonsingleton_cell;
  const Cell *prev_nonsingleton = 0;
  if(next_nonsingleton)
    assert(next_nonsingleton->prev_nonsingleton == 0);
  for(const Cell *cell = first_cell; cell; cell = cell->next)
    {
      assert(!cell->next || cell->next->prev_next_ptr == &(cell->next));
      if(cell->length > 1)
	{
	  if(do_print)
	    fprintf(stderr, "%u:%u=>", cell->first, cell->length);
	  assert(cell == next_nonsingleton); 
	  assert(cell->prev_nonsingleton == prev_nonsingleton);
	  next_nonsingleton = cell->next_nonsingleton;
	  prev_nonsingleton = cell;
	  if(next_nonsingleton)
	    assert(next_nonsingleton->first > cell->first);
	}
      else
	{
	  assert(cell != next_nonsingleton);
	  assert(cell->next_nonsingleton == 0);
	  assert(cell->prev_nonsingleton == 0);
	}
    }
  assert(next_nonsingleton == 0);
  if(do_print)
    fprintf(stderr, "\n");
#endif
}

Cell *Partition::aux_split_in_two(Cell * const cell,
				  const unsigned int first_half_size)
{
  RefInfo i;

  DEBUG_ASSERT(first_half_size > 0);
  DEBUG_ASSERT(first_half_size < cell->length);

  /* (Pseudo)allocate new cell */
  Cell * const new_cell = free_cells;
  DEBUG_ASSERT(new_cell);
  free_cells = new_cell->next;
  /* Update new cell parameters */
  new_cell->first = cell->first + first_half_size;
  new_cell->length = cell->length - first_half_size;
  new_cell->next = cell->next;
  if(new_cell->next)
    new_cell->next->prev_next_ptr = &(new_cell->next);
  new_cell->prev_next_ptr = &(cell->next);
  new_cell->split_level = cell->split_level;
  /* Update old, splitted cell parameters */
  cell->length = first_half_size;
  cell->next = new_cell;
  cell->split_level = level;
  
  /* Add cell in refinement_stack for backtracking */
  i.split_cell_first = cell->first;
  if(cell->prev_nonsingleton)
    i.prev_nonsingleton_first = cell->prev_nonsingleton->first;
  else
    i.prev_nonsingleton_first = -1;
  if(cell->next_nonsingleton)
    i.next_nonsingleton_first = cell->next_nonsingleton->first;
  else
    i.next_nonsingleton_first = -1;
  refinement_stack.push(i);

  /* Modify nonsingleton cell list */
  if(new_cell->length > 1)
    {
      new_cell->prev_nonsingleton = cell;
      new_cell->next_nonsingleton = cell->next_nonsingleton;
      if(new_cell->next_nonsingleton)
	new_cell->next_nonsingleton->prev_nonsingleton = new_cell;
      cell->next_nonsingleton = new_cell;
    }
  else
    {
      new_cell->next_nonsingleton = 0;
      new_cell->prev_nonsingleton = 0;
    }

  if(cell->length == 1)
    {
      if(cell->prev_nonsingleton)
	cell->prev_nonsingleton->next_nonsingleton = cell->next_nonsingleton;
      else
	first_nonsingleton_cell = cell->next_nonsingleton;
      if(cell->next_nonsingleton)
	cell->next_nonsingleton->prev_nonsingleton = cell->prev_nonsingleton;
      cell->next_nonsingleton = 0;
      cell->prev_nonsingleton = 0;
    }

  return new_cell;
} 

void Partition::add_in_splitting_queue(Cell * const cell)
{
  static const unsigned int smallish_cell_threshold = 1;
  assert(!cell->in_splitting_queue);
  cell->in_splitting_queue = true;
  if(cell->length <= smallish_cell_threshold)
    splitting_queue.push_front(cell);
  else
    splitting_queue.push_back(cell);    
}


void Partition::clear_splitting_queue()
{
  while(!splitting_queue.is_empty())
    {
      Cell * const cell = splitting_queue.pop_front();
      assert(cell->in_splitting_queue);
      cell->in_splitting_queue = false;
    }
}



/*
 * Assumes that the invariant values are NOT the same
 * and that the cell contains more than one element
 */
Cell *Partition::sort_and_split_cell1(Cell *cell)
{
#if defined(EXPENSIVE_CONSISTENCY_CHECKS)
  assert(cell->length > 1);
  assert(cell->first + cell->length <= graph->get_nof_vertices());
  bool found0 = false, found1 = false;
  for(unsigned int i = 0; i < cell->length; i++)
    {
      if(invariant_values[elements[cell->first + i]] == 0)
	found0 = true;
      else if(invariant_values[elements[cell->first + i]] == 1)
	found1 = true;
      else
	assert(should_not_happen);
    }
  assert(found0);
  assert(found1);
#endif

  consistency_check();

  /* Allocate new cell */
  Cell *new_cell = free_cells;
  DEBUG_ASSERT(new_cell);
  free_cells = new_cell->next;
  if(free_cells) free_cells->prev_next_ptr = &(free_cells);

  /* Sort vertices in the cell according to the invariant values */
  unsigned int *ep0 = elements + cell->first;
  unsigned int *ep1 = ep0 + cell->length;
  while(ep1 > ep0)
    {
      const unsigned int element = *ep0;
      const unsigned int ival = invariant_values[element];
      invariant_values[element] = 0;
      DEBUG_ASSERT(ival <= 1);
      DEBUG_ASSERT(element_to_cell_map[element] == cell);
      DEBUG_ASSERT(in_pos[element] == ep0);
      if(ival == 0)
	{
	  ep0++;
	}
      else
	{
	  ep1--;
	  *ep0 = *ep1;
	  *ep1 = element;
	  element_to_cell_map[element] = new_cell;
	  in_pos[element] = ep1;
	  in_pos[*ep0] = ep0;
	}
    }

  DEBUG_ASSERT(ep1 != elements + cell->first);
  DEBUG_ASSERT(ep0 != elements + cell->first + cell->length);

  /* Update new cell parameters */
  new_cell->first = ep1 - elements;
  new_cell->length = cell->length - (new_cell->first - cell->first);
  new_cell->next = cell->next;
  if(new_cell->next)
    new_cell->next->prev_next_ptr = &(new_cell->next);
  new_cell->prev_next_ptr = &(cell->next);
  new_cell->split_level = cell->split_level;
  /* Update old, splitted cell parameters */
  cell->length = new_cell->first - cell->first;
  cell->next = new_cell;
  cell->split_level = level;

  /* Add cell in refinement stack for backtracking */
  {
    RefInfo i;
    i.split_cell_first = cell->first;
    if(cell->prev_nonsingleton)
      i.prev_nonsingleton_first = cell->prev_nonsingleton->first;
    else
      i.prev_nonsingleton_first = -1;
    if(cell->next_nonsingleton)
      i.next_nonsingleton_first = cell->next_nonsingleton->first;
    else
      i.next_nonsingleton_first = -1;
    /* Modify nonsingleton cell list */
    if(new_cell->length > 1)
      {
	new_cell->prev_nonsingleton = cell;
	new_cell->next_nonsingleton = cell->next_nonsingleton;
	if(new_cell->next_nonsingleton)
	  new_cell->next_nonsingleton->prev_nonsingleton = new_cell;
	cell->next_nonsingleton = new_cell;
      }
    else
      {
	new_cell->next_nonsingleton = 0;
	new_cell->prev_nonsingleton = 0;
      }
    if(cell->length == 1)
      {
	if(cell->prev_nonsingleton)
	  cell->prev_nonsingleton->next_nonsingleton = cell->next_nonsingleton;
	else
	  first_nonsingleton_cell = cell->next_nonsingleton;
	if(cell->next_nonsingleton)
	  cell->next_nonsingleton->prev_nonsingleton = cell->prev_nonsingleton;
	cell->next_nonsingleton = 0;
	cell->prev_nonsingleton = 0;
      }
    refinement_stack.push(i);
  }


  /* Add cells in splitting queue */
  DEBUG_ASSERT(!new_cell->in_splitting_queue);
  if(cell->in_splitting_queue) {
    /* Both cells must be included in splitting_queue in order to have
       refinement to equitable partition */
    add_in_splitting_queue(new_cell);
  } else {
    Cell *min_cell, *max_cell;
    if(cell->length <= new_cell->length) {
      min_cell = cell;
      max_cell = new_cell;
    } else {
      min_cell = new_cell;
      max_cell = cell;
    }
    /* Put the smaller cell in splitting_queue */
    add_in_splitting_queue(min_cell);
    if(max_cell->length == 1) {
      /* Put the "larger" cell also in splitting_queue */
      add_in_splitting_queue(max_cell);
    }
  }

  consistency_check();

  return new_cell;
}



/*
 * Tables and a subroutine for distribution count sorting
 */
static IGRAPH_THREAD_LOCAL unsigned int count[256] = {
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
};
static IGRAPH_THREAD_LOCAL unsigned int start[256];

/*
 * Build start array so that start[0] = 0 and start[i+1] = start[i] + count[i]
 */
static void cumulate_count(const unsigned int max) 
{
  DEBUG_ASSERT(max <= 255);
  unsigned int *count_p = count;
  unsigned int *start_p = start;
  unsigned int sum = 0;
  for(unsigned int i = max+1; i > 0; i--)
    {
      *start_p = sum;
      start_p++;
      sum += *count_p;
      count_p++;
    }
}


/*
 * Distribution count sorting of cells with invariant values less than 256 
 */
Cell *Partition::sort_and_split_cell255(Cell * const cell,
					const unsigned int max_ival)
{
  //DEBUG_ASSERT(cell->first + cell->length <= graph->vertices.size());

  DEBUG_ASSERT(max_ival <= 255);

  if(cell->length == 1)
    {
      /* Reset invariant value */
      invariant_values[elements[cell->first]] = 0;
      return cell;
    }
  
#ifdef CONSISTENCY_CHECKS
  for(unsigned int i = 0; i < 256; i++)
    assert(count[i] == 0);
#endif

  /*
   * Compute the distribution of invariant values to the count array
   */
  {
    const unsigned int *ep = elements + cell->first;
    DEBUG_ASSERT(element_to_cell_map[*ep] == cell);
    const unsigned int ival = invariant_values[*ep];
    DEBUG_ASSERT(ival <= 255);
    count[ival]++;
    ep++;
#ifdef CONSISTENCY_CHECKS
    bool equal_invariant_values = true;
#endif
    for(unsigned int i = cell->length - 1; i > 0; i--)
      {
	DEBUG_ASSERT(element_to_cell_map[*ep] == cell);
	const unsigned int ival2 = invariant_values[*ep];
	DEBUG_ASSERT(ival2 <= 255);
	DEBUG_ASSERT(ival2 <= max_ival);
	count[ival2]++;
#ifdef CONSISTENCY_CHECKS
	if(ival2 != ival) {
	  equal_invariant_values = false;
	}
#endif
	ep++;
      }
#ifdef CONSISTENCY_CHECKS
    DEBUG_ASSERT(!equal_invariant_values);
    if(equal_invariant_values) {
      DEBUG_ASSERT(count[ival] == cell->length);
      count[ival] = 0;
      clear_ivs(cell);
      return cell;
    }
#endif
  }

  /* Build start array */
  cumulate_count(max_ival);

  //DEBUG_ASSERT(start[255] + count[255] == cell->length);
  DEBUG_ASSERT(start[max_ival] + count[max_ival] == cell->length);

  /* Do the sorting */
  for(unsigned int i = 0; i <= max_ival; i++)
    {
      unsigned int *ep = elements + cell->first + start[i];
      for(unsigned int j = count[i]; j > 0; j--)
	{
	  while(true)
	    {
	      const unsigned int element = *ep;
	      const unsigned int ival = invariant_values[element];
	      if(ival == i)
		break;
	      DEBUG_ASSERT(ival > i);
	      DEBUG_ASSERT(count[ival] > 0);
	      *ep = elements[cell->first + start[ival]];
	      elements[cell->first + start[ival]] = element;
	      start[ival]++;
	      count[ival]--;
	    }
	  ep++;
	}
      count[i] = 0;
    }

#if defined(CONSISTENCY_CHECKS)
  for(unsigned int i = 0; i < 256; i++)
    assert(count[i] == 0);
#endif

#if defined(VERBOSEDEBUG)
  {
    const unsigned int *ep = elements + cell->first;
    fprintf(stderr, "\n");
    for(unsigned int i = cell->length; i > 0; i--, ep++)
      fprintf(stderr, "%u ", invariant_values[*ep]);
    fprintf(stderr, "\n");
  }
#endif

  /* split cell */
  Cell * const cell2 = split_cell(cell);
  DEBUG_ASSERT(cell2 != cell);
  return cell2;
}


/*
 * Sort the elements in a cell according to their invariant values
 * The invariant values are not cleared
 * Warning: the in_pos array is left in incorrect state
 */
bool Partition::shellsort_cell(Cell *cell)
{
  unsigned int h;
  unsigned int *ep;

  //DEBUG_ASSERT(cell->first + cell->length <= graph->vertices.size());

  if(cell->length == 1)
    return false;

  /* Check whether all the elements have the same invariant value */
  bool equal_invariant_values = true;
  {
    ep = elements + cell->first;
    const unsigned int ival = invariant_values[*ep];
    DEBUG_ASSERT(element_to_cell_map[*ep] == cell);
    ep++;
    for(unsigned int i = cell->length - 1; i > 0; i--)
      {
	DEBUG_ASSERT(element_to_cell_map[*ep] == cell);
	if(invariant_values[*ep] != ival) {
	  equal_invariant_values = false;
	  break;
	}
	ep++;
      }
  }
  if(equal_invariant_values)
    return false;

  ep = elements + cell->first;

  for(h = 1; h <= cell->length/9; h = 3*h + 1)
    ;
  for( ; h > 0; h = h/3) {
    for(unsigned int i = h; i < cell->length; i++) {
      const unsigned int element = ep[i];
      const unsigned int ival = invariant_values[element];
      unsigned int j = i;
      while(j >= h && invariant_values[ep[j-h]] > ival) {
        ep[j] = ep[j-h];
        j -= h;
      }
      ep[j] = element;
    }
  }
  return true;
}


void Partition::clear_ivs(Cell * const cell)
{
  unsigned int *ep = elements + cell->first;
  for(unsigned int i = cell->length; i > 0; i--, ep++)
    invariant_values[*ep] = 0;
}


/*
 * Assumes that the elements in the cell are sorted according to their
 * invariant values
 */
Cell *Partition::split_cell(Cell * const original_cell)
{
  Cell *cell = original_cell;
  const bool original_cell_was_in_splitting_queue =
    original_cell->in_splitting_queue;
  Cell *largest_new_cell = 0;

  consistency_check();

  while(true) 
    {
      unsigned int *ep = elements + cell->first;
      const unsigned int * const lp = ep + cell->length;
      const unsigned int ival = invariant_values[*ep];
      invariant_values[*ep] = 0;
      element_to_cell_map[*ep] = cell;
      in_pos[*ep] = ep;
      ep++;
      while(ep < lp)
	{
	  const unsigned int e = *ep;
	  if(invariant_values[e] != ival)
	    break;
	  invariant_values[e] = 0;
	  in_pos[e] = ep;
	  ep++;
	  element_to_cell_map[e] = cell;
	}
      if(ep == lp)
	break;
      
      Cell * const new_cell = aux_split_in_two(cell, (ep - elements) - cell->first);
      
      if(graph->in_search)
	{
	  graph->eqref_hash.update(new_cell->first);
	  graph->eqref_hash.update(new_cell->length);
	  graph->eqref_hash.update(ival);
	}
      
      /* Add cells in splitting_queue */
      assert(!new_cell->in_splitting_queue);
      if(original_cell_was_in_splitting_queue)
	{
	  /* In this case, all new cells are inserted in splitting_queue */
	  assert(cell->in_splitting_queue);
	  add_in_splitting_queue(new_cell);
	}
      else
	{
	  /* Otherwise, we can omit one new cell from splitting_queue */
	  assert(!cell->in_splitting_queue);
	  if(largest_new_cell == 0) {
	    largest_new_cell = cell;
	  } else {
	    assert(!largest_new_cell->in_splitting_queue);
	    if(cell->length > largest_new_cell->length) {
	      add_in_splitting_queue(largest_new_cell);
	      largest_new_cell = cell;
	    } else {
	      add_in_splitting_queue(cell);
	    }
	  }
	}
      /* Process the rest of the cell */
      cell = new_cell;
    }

  consistency_check();
  
  if(original_cell == cell) {
    /* All the elements in cell had the same invariant value */
    return cell;
  }

  /* Add cells in splitting_queue */
  if(!original_cell_was_in_splitting_queue)
    {
      /* Also consider the last new cell */
      assert(largest_new_cell);
      if(cell->length > largest_new_cell->length)
	{
	  add_in_splitting_queue(largest_new_cell);
	  largest_new_cell = cell;
	}
      else
	{
	  add_in_splitting_queue(cell);
	}
      if(largest_new_cell->length == 1)
	{
	  /* Needed in certificate computation */
	  add_in_splitting_queue(largest_new_cell);
	}
    }

  return cell;
}


Cell *Partition::zplit_cell(Cell * const cell, const bool max_ival_info_ok)
{
  assert(cell);

  Cell *last_new_cell = cell;

  if(!max_ival_info_ok)
    {
      /* Compute max_ival info */
      assert(cell->max_ival == 0);
      assert(cell->max_ival_count == 0);
      unsigned int *ep = elements + cell->first;
      for(unsigned int i = cell->length; i > 0; i--, ep++)
	{
	  const unsigned int ival = invariant_values[*ep];
	  if(ival > cell->max_ival)
	    {
	      cell->max_ival = ival;
	      cell->max_ival_count = 1;
	    }
	  else if(ival == cell->max_ival)
	    {
	      cell->max_ival_count++;
	    }
	}
    }

#ifdef CONSISTENCY_CHECKS
  /* Verify max_ival info */
  {
    unsigned int max_ival = 0;
    unsigned int max_ival_count = 0;
    unsigned int *ep = elements + cell->first;
    for(unsigned int i = cell->length; i > 0; i--, ep++)
      {
	const unsigned int ival = invariant_values[*ep];
	if(ival > max_ival)
	  {
	    max_ival = ival;
	    max_ival_count = 1;
	  }
	else if(ival == max_ival)
	  {
	    max_ival_count++;
	  }
      }
    assert(max_ival == cell->max_ival);
    assert(max_ival_count == cell->max_ival_count);
  }
#endif

  /* max_ival info has been computed */

  if(cell->max_ival_count == cell->length)
    {
      /* All invariant values are the same */
      if(cell->max_ival > 0)
	clear_ivs(cell);
      goto done;
    }

  /* All invariant values are not the same */
  if(cell->max_ival == 1)
    {
      /* Specialized splitting for cells with binary invariant values */
      last_new_cell = sort_and_split_cell1(cell);
      goto done;
    }
  if(cell->max_ival < 256)
    {
      /* Specialized splitting for cells with invariant values < 256 */
      last_new_cell = sort_and_split_cell255(cell, cell->max_ival);
      goto done;
    }
  {
    /* Generic sorting and splitting */
    const bool sorted = shellsort_cell(cell);
    assert(sorted);
    last_new_cell = split_cell(cell);
    goto done;
  }

 done:
  cell->max_ival = 0;
  cell->max_ival_count = 0;
  return last_new_cell;
}






void Partition::unrefine(unsigned int dest_split_level,
			 unsigned int dest_refinement_stack_size)
{
  assert(refinement_stack.size() >= dest_refinement_stack_size);
  while(refinement_stack.size() > dest_refinement_stack_size)
    {
      RefInfo i = refinement_stack.pop();
      const unsigned int first = i.split_cell_first;
      //const unsigned int first = refinement_stack.pop();
      Cell *cell = element_to_cell_map[elements[first]];

      if(cell->first != first) {
	assert(cell->split_level <= dest_split_level);
	goto done;
      }
      if(cell->split_level <= dest_split_level) {
	goto done;
      }

      {
	const unsigned int new_first = cell->first;
	do
	  {
	    Cell * const next_cell = cell->next;
	    assert(next_cell);
	    /* (Pseudo)free cell */
	    cell->first = 0;
	    cell->length = 0;
	    cell->next->prev_next_ptr = cell->prev_next_ptr;
	    *(cell->prev_next_ptr) = cell->next;
	    cell->next = free_cells;
	    if(cell->next)
	      cell->next->prev_next_ptr = &(cell->next);
	    cell->prev_next_ptr = &free_cells;
	    free_cells = cell;
	    
	    cell = next_cell;
	  }
	while(cell->split_level > dest_split_level);
	/* Update element_to_cell_map values of elements added in cell */
	unsigned int *ep = elements + new_first;
	unsigned int * const lp = elements + cell->first;
	while(ep < lp) {
	  element_to_cell_map[*ep] = cell;
	  ep++;
	}
	/* Update cell parameters */
	cell->length = (cell->first + cell->length) - new_first;
	cell->first = new_first;
      }

    done:
      if(i.prev_nonsingleton_first >= 0)
	{
	  Cell * const prev_cell = element_to_cell_map[elements[i.prev_nonsingleton_first]];
	  DEBUG_ASSERT(prev_cell->length > 1);
	  cell->prev_nonsingleton = prev_cell;
	  prev_cell->next_nonsingleton = cell;
	}
      else
	{
	  //assert(cell->prev_nonsingleton == 0);
	  cell->prev_nonsingleton = 0;
	  first_nonsingleton_cell = cell;
	}

      if(i.next_nonsingleton_first >= 0)
	{
	  Cell * const next_cell = element_to_cell_map[elements[i.next_nonsingleton_first]];
	  DEBUG_ASSERT(next_cell->length > 1);
	  cell->next_nonsingleton = next_cell;
	  next_cell->prev_nonsingleton = cell;
	}
      else
	{
	  //assert(cell->next_nonsingleton == 0);
	  cell->next_nonsingleton = 0;
	}
    }

  consistency_check();
}

}
