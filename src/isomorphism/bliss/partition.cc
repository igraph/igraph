#include <new>
#include <cassert>

#include "graph.hh"
#include "partition.hh"

/* Allow using 'and' instead of '&&' with MSVC */
#if _MSC_VER
#include <ciso646>
#endif

/*
  Copyright (c) 2003-2021 Tommi Junttila
  Released under the GNU Lesser General Public License version 3.

  This file is part of bliss.

  bliss is free software: you can redistribute it and/or modify
  it under the terms of the GNU Lesser General Public License as published by
  the Free Software Foundation, version 3 of the License.

  bliss is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public License
  along with bliss.  If not, see <http://www.gnu.org/licenses/>.
*/

namespace bliss {

Partition::Partition()
{
  N = 0;
  elements = 0;
  in_pos = 0;
  invariant_values = 0;
  cells = 0;
  free_cells = 0;
  element_to_cell_map = 0;
  graph = 0;
  discrete_cell_count = 0;
  /* Initialize a distribution count sorting array. */
  for(unsigned int i = 0; i < 256; i++)
    dcs_count[i] = 0;

  cr_enabled = false;
  cr_cells = 0;
  cr_levels = 0;
}



Partition::~Partition()
{
  delete[] elements; elements = nullptr;
  delete[] cells; cells = nullptr;
  delete[] element_to_cell_map; element_to_cell_map = nullptr;
  delete[] in_pos; in_pos = nullptr;
  delete[] invariant_values; invariant_values = nullptr;
  N = 0;
}



void Partition::init(const unsigned int M)
{
  assert(M > 0);
  N = M;

  delete[] elements;
  elements = new unsigned int[N];
  for(unsigned int i = 0; i < N; i++)
    elements[i] = i;

  delete[] in_pos;
  in_pos = new unsigned int*[N];
  for(unsigned int i = 0; i < N; i++)
    in_pos[i] = elements + i;

  delete[] invariant_values;
  invariant_values = new unsigned int[N];
  for(unsigned int i = 0; i < N; i++)
    invariant_values[i] = 0;

  delete[] cells;
  cells = new Cell[N];

  cells[0].first = 0;
  cells[0].length = N;
  cells[0].max_ival = 0;
  cells[0].max_ival_count = 0;
  cells[0].in_splitting_queue = false;
  cells[0].in_neighbour_heap = false;
  cells[0].prev = 0;
  cells[0].next = 0;
  cells[0].next_nonsingleton = 0;
  cells[0].prev_nonsingleton = 0;
  cells[0].split_level = 0;
  first_cell = &cells[0];
  if(N == 1)
    {
      first_nonsingleton_cell = 0;
      discrete_cell_count = 1;
    }
  else
    {
      first_nonsingleton_cell = &cells[0];
      discrete_cell_count = 0;
    }

  for(unsigned int i = 1; i < N; i++)
    {
      cells[i].first = 0;
      cells[i].length = 0;
      cells[i].max_ival = 0;
      cells[i].max_ival_count = 0;
      cells[i].in_splitting_queue = false;
      cells[i].in_neighbour_heap = false;
      cells[i].prev = 0;
      cells[i].next = (i < N-1)?&cells[i+1]:0;
      cells[i].next_nonsingleton = 0;
      cells[i].prev_nonsingleton = 0;
    }
  if(N > 1)
    free_cells = &cells[1];
  else
    free_cells = 0;

  delete[] element_to_cell_map;
  element_to_cell_map = new Cell*[N];
  for(unsigned int i = 0; i < N; i++)
    element_to_cell_map[i] = first_cell;

  splitting_queue.init(N);
  refinement_stack.init(N);

  /* Reset the main backtracking stack */
  bt_stack.clear();
}






Partition::BacktrackPoint
Partition::set_backtrack_point()
{
  BacktrackInfo info;
  info.refinement_stack_size = refinement_stack.size();
  if(cr_enabled)
    info.cr_backtrack_point = cr_get_backtrack_point();
  BacktrackPoint p = bt_stack.size();
  bt_stack.push_back(info);
  return p;
}



void
Partition::goto_backtrack_point(BacktrackPoint p)
{
  assert(p < bt_stack.size());
  BacktrackInfo info = bt_stack[p];
  bt_stack.resize(p);

  if(cr_enabled)
    cr_goto_backtrack_point(info.cr_backtrack_point);

  const unsigned int dest_refinement_stack_size = info.refinement_stack_size;

  assert(refinement_stack.size() >= dest_refinement_stack_size);
  while(refinement_stack.size() > dest_refinement_stack_size)
    {
      RefInfo i = refinement_stack.pop();
      const unsigned int first = i.split_cell_first;
      Cell* cell = get_cell(elements[first]);

      if(cell->first != first)
        {
          assert(cell->first < first);
          assert(cell->split_level <= dest_refinement_stack_size);
          goto done;
        }
      assert(cell->split_level > dest_refinement_stack_size);

      while(cell->split_level > dest_refinement_stack_size)
        {
          assert(cell->prev);
          cell = cell->prev;
        }
      while(cell->next and
            cell->next->split_level > dest_refinement_stack_size)
        {
          /* Merge next cell */
          Cell* const next_cell = cell->next;
          if(cell->length == 1)
            discrete_cell_count--;
          if(next_cell->length == 1)
            discrete_cell_count--;
          /* Update element_to_cell_map values of elements added in cell */
          unsigned int* ep = elements + next_cell->first;
          unsigned int* const lp = ep + next_cell->length;
          for( ; ep < lp; ep++)
            element_to_cell_map[*ep] = cell;
          /* Update cell parameters */
          cell->length += next_cell->length;
          if(next_cell->next)
            next_cell->next->prev = cell;
          cell->next = next_cell->next;
          /* (Pseudo)free next_cell */
          next_cell->first = 0;
          next_cell->length = 0;
          next_cell->prev = 0;
          next_cell->next = free_cells;
          free_cells = next_cell;
        }

    done:
      if(i.prev_nonsingleton_first >= 0)
        {
          Cell* const prev_cell = get_cell(elements[i.prev_nonsingleton_first]);
          assert(prev_cell->length > 1);
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
          Cell* const next_cell = get_cell(elements[i.next_nonsingleton_first]);
          assert(next_cell->length > 1);
          cell->next_nonsingleton = next_cell;
          next_cell->prev_nonsingleton = cell;
        }
      else
        {
          //assert(cell->next_nonsingleton == 0);
          cell->next_nonsingleton = 0;
        }
    }

}



Partition::Cell*
Partition::individualize(Partition::Cell * const cell,
                         const unsigned int element)
{
  assert(!cell->is_unit());

  unsigned int * const pos = in_pos[element];
  assert((unsigned int)(pos - elements) >= cell->first);
  assert((unsigned int)(pos - elements) < cell->first + cell->length);
  assert(*pos == element);

  const unsigned int last = cell->first + cell->length - 1;
  *pos = elements[last];
  in_pos[*pos] = pos;
  elements[last] = element;
  in_pos[element] = elements + last;

  Partition::Cell * const new_cell = aux_split_in_two(cell, cell->length-1);
  assert(elements[new_cell->first] == element);
  element_to_cell_map[element] = new_cell;

  return new_cell;
}



Partition::Cell*
Partition::aux_split_in_two(Partition::Cell* const cell,
                            const unsigned int first_half_size)
{
  RefInfo i;

  assert(0 < first_half_size && first_half_size < cell->length);

  /* (Pseudo)allocate new cell */
  Cell * const new_cell = free_cells;
  assert(new_cell != 0);
  free_cells = new_cell->next;
  /* Update new cell parameters */
  new_cell->first = cell->first + first_half_size;
  new_cell->length = cell->length - first_half_size;
  new_cell->next = cell->next;
  if(new_cell->next)
    new_cell->next->prev = new_cell;
  new_cell->prev = cell;
  new_cell->split_level = refinement_stack.size()+1;
  /* Update old, splitted cell parameters */
  cell->length = first_half_size;
  cell->next = new_cell;
  /* CR */
  if(cr_enabled)
    cr_create_at_level_trailed(new_cell->first, cr_get_level(cell->first));

  /* Add cell in refinement_stack for backtracking */
  i.split_cell_first = new_cell->first;
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
      discrete_cell_count++;
    }

  if(cell->is_unit())
    {
      if(cell->prev_nonsingleton)
        cell->prev_nonsingleton->next_nonsingleton = cell->next_nonsingleton;
      else
        first_nonsingleton_cell = cell->next_nonsingleton;
      if(cell->next_nonsingleton)
        cell->next_nonsingleton->prev_nonsingleton = cell->prev_nonsingleton;
      cell->next_nonsingleton = 0;
      cell->prev_nonsingleton = 0;
      discrete_cell_count++;
    }

  return new_cell;
}



void
Partition::splitting_queue_add(Cell* const cell)
{
  static const unsigned int smallish_cell_threshold = 1;
  assert(!cell->in_splitting_queue);
  cell->in_splitting_queue = true;
  if(cell->length <= smallish_cell_threshold)
    splitting_queue.push_front(cell);
  else
    splitting_queue.push_back(cell);
}



void
Partition::splitting_queue_clear()
{
  while(!splitting_queue_is_empty())
    splitting_queue_pop();
}





/*
 * Assumes that the invariant values are NOT the same
 * and that the cell contains more than one element
 */
Partition::Cell*
Partition::sort_and_split_cell1(Partition::Cell* const cell)
{
#if defined(BLISS_EXPENSIVE_CONSISTENCY_CHECKS)
  assert(cell->length > 1);
  assert(cell->first + cell->length <= N);
  unsigned int nof_0_found = 0;
  unsigned int nof_1_found = 0;
  for(unsigned int i = cell->first; i < cell->first + cell->length; i++)
    {
      const unsigned int ival = invariant_values[elements[i]];
      assert(ival == 0 or ival == 1);
      if(ival == 0) nof_0_found++;
      else nof_1_found++;
    }
  assert(nof_0_found > 0);
  assert(nof_1_found > 0);
  assert(nof_1_found == cell->max_ival_count);
  assert(nof_0_found + nof_1_found == cell->length);
  assert(cell->max_ival == 1);
#endif


  /* (Pseudo)allocate new cell */
  Cell* const new_cell = free_cells;
  assert(new_cell != 0);
  free_cells = new_cell->next;

#define NEW_SORT1
#ifdef NEW_SORT1
      unsigned int *ep0 = elements + cell->first;
      unsigned int *ep1 = ep0 + cell->length - cell->max_ival_count;
      if(cell->max_ival_count > cell->length / 2)
        {
          /* There are more ones than zeros, only move zeros */
          unsigned int * const end = ep0 + cell->length;
          while(ep1 < end)
            {
              while(invariant_values[*ep1] == 0)
                {
                  const unsigned int tmp = *ep1;
                  *ep1 = *ep0;
                  *ep0 = tmp;
                  in_pos[tmp] = ep0;
                  in_pos[*ep1] = ep1;
                  ep0++;
                }
              element_to_cell_map[*ep1] = new_cell;
              invariant_values[*ep1] = 0;
              ep1++;
            }
        }
      else
        {
          /* There are more zeros than ones, only move ones */
          unsigned int * const end = ep1;
          while(ep0 < end)
            {
              while(invariant_values[*ep0] != 0)
                {
                  const unsigned int tmp = *ep0;
                  *ep0 = *ep1;
                  *ep1 = tmp;
                  in_pos[tmp] = ep1;
                  in_pos[*ep0] = ep0;
                  ep1++;
                }
              ep0++;
            }
          ep1 = end;
          while(ep1 < elements + cell->first + cell->length)
            {
              element_to_cell_map[*ep1] = new_cell;
              invariant_values[*ep1] = 0;
              ep1++;
            }
        }
  /* Update new cell parameters */
  new_cell->first = cell->first + cell->length - cell->max_ival_count;
  new_cell->length = cell->length - (new_cell->first - cell->first);
  new_cell->next = cell->next;
  if(new_cell->next)
    new_cell->next->prev = new_cell;
  new_cell->prev = cell;
  new_cell->split_level = refinement_stack.size()+1;
  /* Update old, splitted cell parameters */
  cell->length = new_cell->first - cell->first;
  cell->next = new_cell;
  /* CR */
  if(cr_enabled)
    cr_create_at_level_trailed(new_cell->first, cr_get_level(cell->first));

#else
  /* Sort vertices in the cell according to the invariant values */
  unsigned int *ep0 = elements + cell->first;
  unsigned int *ep1 = ep0 + cell->length;
  while(ep1 > ep0)
    {
      const unsigned int element = *ep0;
      const unsigned int ival = invariant_values[element];
      invariant_values[element] = 0;
      assert(ival <= 1);
      assert(element_to_cell_map[element] == cell);
      assert(in_pos[element] == ep0);
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

  assert(ep1 != elements + cell->first);
  assert(ep0 != elements + cell->first + cell->length);

  /* Update new cell parameters */
  new_cell->first = ep1 - elements;
  new_cell->length = cell->length - (new_cell->first - cell->first);
  new_cell->next = cell->next;
  if(new_cell->next)
    new_cell->next->prev = new_cell;
  new_cell->prev = cell;
  new_cell->split_level = cell->split_level;
  /* Update old, splitted cell parameters */
  cell->length = new_cell->first - cell->first;
  cell->next = new_cell;
  cell->split_level = refinement_stack.size()+1;
  /* CR */
  if(cr_enabled)
    cr_create_at_level_trailed(new_cell->first, cr_get_level(cell->first));

#endif /* ifdef NEW_SORT1*/

  /* Add cell in refinement stack for backtracking */
  {
    RefInfo i;
    i.split_cell_first = new_cell->first;
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
        discrete_cell_count++;
      }
    if(cell->is_unit())
      {
        if(cell->prev_nonsingleton)
          cell->prev_nonsingleton->next_nonsingleton = cell->next_nonsingleton;
        else
          first_nonsingleton_cell = cell->next_nonsingleton;
        if(cell->next_nonsingleton)
          cell->next_nonsingleton->prev_nonsingleton = cell->prev_nonsingleton;
        cell->next_nonsingleton = 0;
        cell->prev_nonsingleton = 0;
        discrete_cell_count++;
      }
    refinement_stack.push(i);
  }


  /* Add cells in splitting queue */
  assert(!new_cell->in_splitting_queue);
  if(cell->in_splitting_queue) {
    /* Both cells must be included in splitting_queue in order to have
       refinement to equitable partition */
    splitting_queue_add(new_cell);
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
    splitting_queue_add(min_cell);
    if(max_cell->is_unit()) {
      /* Put the "larger" cell also in splitting_queue */
      splitting_queue_add(max_cell);
    }
  }


  return new_cell;
}





/**
 * An auxiliary function for distribution count sorting.
 * Build start array so that
 * dcs_start[0] = 0 and dcs_start[i+1] = dcs_start[i] + dcs_count[i].
 */
void
Partition::dcs_cumulate_count(const unsigned int max)
{
  assert(max <= 255);
  unsigned int* count_p = dcs_count;
  unsigned int* start_p = dcs_start;
  unsigned int sum = 0;
  for(unsigned int i = max+1; i > 0; i--)
    {
      *start_p = sum;
      start_p++;
      sum += *count_p;
      count_p++;
    }
}


/**
 * Distribution count sorting of cells with invariant values less than 256.
 */
Partition::Cell*
Partition::sort_and_split_cell255(Partition::Cell* const cell,
                                  const unsigned int max_ival)
{
  assert(max_ival <= 255);

  if(cell->is_unit())
    {
      /* Reset invariant value */
      invariant_values[elements[cell->first]] = 0;
      return cell;
    }

#ifdef BLISS_CONSISTENCY_CHECKS
  for(unsigned int i = 0; i < 256; i++)
    assert(dcs_count[i] == 0);
#endif

  /*
   * Compute the distribution of invariant values to the count array
   */
  {
    const unsigned int *ep = elements + cell->first;
    assert(element_to_cell_map[*ep] == cell);
    const unsigned int ival = invariant_values[*ep];
    assert(ival <= 255);
    dcs_count[ival]++;
    ep++;
#if defined(BLISS_CONSISTENCY_CHECKS)
    bool equal_invariant_values = true;
#endif
    for(unsigned int i = cell->length - 1; i != 0; i--)
      {
        assert(element_to_cell_map[*ep] == cell);
        const unsigned int ival2 = invariant_values[*ep];
        assert(ival2 <= 255);
        assert(ival2 <= max_ival);
        dcs_count[ival2]++;
#if defined(BLISS_CONSISTENCY_CHECKS)
        if(ival2 != ival) {
          equal_invariant_values = false;
        }
#endif
        ep++;
      }
#if defined(BLISS_CONSISTENCY_CHECKS)
    assert(!equal_invariant_values);
    if(equal_invariant_values) {
      assert(dcs_count[ival] == cell->length);
      dcs_count[ival] = 0;
      clear_ivs(cell);
      return cell;
    }
#endif
  }

  /* Build start array */
  dcs_cumulate_count(max_ival);

  //assert(dcs_start[255] + dcs_count[255] == cell->length);
  assert(dcs_start[max_ival] + dcs_count[max_ival] == cell->length);

  /* Do the sorting */
  for(unsigned int i = 0; i <= max_ival; i++)
    {
      unsigned int *ep = elements + cell->first + dcs_start[i];
      for(unsigned int j = dcs_count[i]; j > 0; j--)
        {
          while(true)
            {
              const unsigned int element = *ep;
              const unsigned int ival = invariant_values[element];
              if(ival == i)
                break;
              assert(ival > i);
              assert(dcs_count[ival] > 0);
              *ep = elements[cell->first + dcs_start[ival]];
              elements[cell->first + dcs_start[ival]] = element;
              dcs_start[ival]++;
              dcs_count[ival]--;
            }
          ep++;
        }
      dcs_count[i] = 0;
    }

#if defined(BLISS_CONSISTENCY_CHECKS)
  for(unsigned int i = 0; i < 256; i++)
    assert(dcs_count[i] == 0);
#endif

  /* split cell */
  Cell* const new_cell = split_cell(cell);
  assert(new_cell != cell);
  return new_cell;
}



/*
 * Sort the elements in a cell according to their invariant values.
 * The invariant values are not cleared.
 * Warning: the in_pos array is left in incorrect state.
 */
bool
Partition::shellsort_cell(Partition::Cell* const cell)
{
  unsigned int h;
  unsigned int* ep;

  //assert(cell->first + cell->length <= N);

  if(cell->is_unit())
    return false;

  /* Check whether all the elements have the same invariant value */
  bool equal_invariant_values = true;
  {
    ep = elements + cell->first;
    const unsigned int ival = invariant_values[*ep];
    assert(element_to_cell_map[*ep] == cell);
    ep++;
    for(unsigned int i = cell->length - 1; i > 0; i--)
      {
        assert(element_to_cell_map[*ep] == cell);
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
      while(j >= h and invariant_values[ep[j-h]] > ival) {
        ep[j] = ep[j-h];
        j -= h;
      }
      ep[j] = element;
    }
  }
  return true;
}



void
Partition::clear_ivs(Cell* const cell)
{
  unsigned int* ep = elements + cell->first;
  for(unsigned int i = cell->length; i > 0; i--, ep++)
    invariant_values[*ep] = 0;
}


/*
 * Assumes that the elements in the cell are sorted according to their
 * invariant values.
 */
Partition::Cell*
Partition::split_cell(Partition::Cell* const original_cell)
{
  Cell* cell = original_cell;
  const bool original_cell_was_in_splitting_queue =
    original_cell->in_splitting_queue;
  Cell* largest_new_cell = 0;

  while(true)
    {
      unsigned int* ep = elements + cell->first;
      const unsigned int* const lp = ep + cell->length;
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

      Cell* const new_cell = aux_split_in_two(cell,
                                              (ep - elements) - cell->first);

      if(graph and graph->compute_eqref_hash)
        {
          graph->eqref_hash.update(new_cell->first);
          graph->eqref_hash.update(new_cell->length);
          graph->eqref_hash.update(ival);
        }

      /* Add cells in splitting_queue */
      assert(!new_cell->is_in_splitting_queue());
      if(original_cell_was_in_splitting_queue)
        {
          /* In this case, all new cells are inserted in splitting_queue */
          assert(cell->is_in_splitting_queue());
          splitting_queue_add(new_cell);
        }
      else
        {
          /* Otherwise, we can omit one new cell from splitting_queue */
          assert(!cell->is_in_splitting_queue());
          if(largest_new_cell == 0) {
            largest_new_cell = cell;
          } else {
            assert(!largest_new_cell->is_in_splitting_queue());
            if(cell->length > largest_new_cell->length) {
              splitting_queue_add(largest_new_cell);
              largest_new_cell = cell;
            } else {
              splitting_queue_add(cell);
            }
          }
        }
      /* Process the rest of the cell */
      cell = new_cell;
    }


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
          splitting_queue_add(largest_new_cell);
          largest_new_cell = cell;
        }
      else
        {
          splitting_queue_add(cell);
        }
      if(largest_new_cell->is_unit())
        {
          /* Needed in certificate computation */
          splitting_queue_add(largest_new_cell);
        }
    }

  return cell;
}


Partition::Cell*
Partition::zplit_cell(Partition::Cell* const cell,
                      const bool max_ival_info_ok)
{
  assert(cell != 0);

  Cell* last_new_cell = cell;

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

#ifdef BLISS_CONSISTENCY_CHECKS
  /* Verify max_ival info */
  {
    unsigned int nof_zeros = 0;
    unsigned int max_ival = 0;
    unsigned int max_ival_count = 0;
    unsigned int *ep = elements + cell->first;
    for(unsigned int i = cell->length; i > 0; i--, ep++)
      {
        const unsigned int ival = invariant_values[*ep];
        if(ival == 0)
          nof_zeros++;
        if(ival > max_ival)
          {
            max_ival = ival;
            max_ival_count = 1;
          }
        else if(ival == max_ival)
          max_ival_count++;
      }
    assert(max_ival == cell->max_ival);
    assert(max_ival_count == cell->max_ival_count);
  }
#endif

  /* max_ival info has been computed */

  if(cell->max_ival_count == cell->length)
    {
      /* All invariant values are the same, clear 'em */
      if(cell->max_ival > 0)
        clear_ivs(cell);
    }
  else
    {
      /* All invariant values are not the same */
      if(cell->max_ival == 1)
        {
          /* Specialized splitting for cells with binary invariant values */
          last_new_cell = sort_and_split_cell1(cell);
        }
      else if(cell->max_ival < 256)
        {
          /* Specialized splitting for cells with invariant values < 256 */
          last_new_cell = sort_and_split_cell255(cell, cell->max_ival);
        }
      else
        {
          /* Generic sorting and splitting */
          const bool sorted = shellsort_cell(cell);
          assert(sorted);
          last_new_cell = split_cell(cell);
        }
    }
  cell->max_ival = 0;
  cell->max_ival_count = 0;
  return last_new_cell;
}



/*
 *
 * Component recursion specific code
 *
 */
void
Partition::cr_init()
{
  assert(bt_stack.empty());

  cr_enabled = true;

  delete[] cr_cells;
  cr_cells = new CRCell[N];

  delete[] cr_levels;
  cr_levels = new CRCell*[N];

  for(unsigned int i = 0; i < N; i++) {
    cr_levels[i] = 0;
    cr_cells[i].level = UINT_MAX;
    cr_cells[i].next = 0;
    cr_cells[i].prev_next_ptr = 0;
  }

  for(const Cell *cell = first_cell; cell; cell = cell->next)
    cr_create_at_level_trailed(cell->first, 0);

  cr_max_level = 0;
}


void
Partition::cr_free()
{
  delete[] cr_cells; cr_cells = nullptr;
  delete[] cr_levels; cr_levels = nullptr;

  cr_created_trail.clear();
  cr_splitted_level_trail.clear();
  cr_bt_info.clear();
  cr_max_level = 0;

  cr_enabled = false;
}


unsigned int
Partition::cr_split_level(const unsigned int level,
                          const std::vector<unsigned int>& splitted_cells)
{
  assert(cr_enabled);
  assert(level <= cr_max_level);
  cr_levels[++cr_max_level] = 0;
  cr_splitted_level_trail.push_back(level);

  for(unsigned int i = 0; i < splitted_cells.size(); i++)
    {
      const unsigned int cell_index = splitted_cells[i];
      assert(cell_index < N);
      CRCell& cr_cell = cr_cells[cell_index];
      assert(cr_cell.level == level);
      cr_cell.detach();
      cr_create_at_level(cell_index, cr_max_level);
    }

  return cr_max_level;
}


unsigned int
Partition::cr_get_backtrack_point()
{
  assert(cr_enabled);
  CR_BTInfo info;
  info.created_trail_index = cr_created_trail.size();
  info.splitted_level_trail_index = cr_splitted_level_trail.size();
  cr_bt_info.push_back(info);
  return cr_bt_info.size()-1;
}


void
Partition::cr_goto_backtrack_point(const unsigned int btpoint)
{
  assert(cr_enabled);
  assert(btpoint < cr_bt_info.size());
  while(cr_created_trail.size() > cr_bt_info[btpoint].created_trail_index)
    {
      const unsigned int cell_index = cr_created_trail.back();
      cr_created_trail.pop_back();
      CRCell& cr_cell = cr_cells[cell_index];
      assert(cr_cell.level != UINT_MAX);
      assert(cr_cell.prev_next_ptr);
      cr_cell.detach();
    }

  while(cr_splitted_level_trail.size() >
        cr_bt_info[btpoint].splitted_level_trail_index)
    {
      const unsigned int dest_level = cr_splitted_level_trail.back();
      cr_splitted_level_trail.pop_back();
      assert(cr_max_level > 0);
      assert(dest_level < cr_max_level);
      while(cr_levels[cr_max_level]) {
        CRCell *cr_cell = cr_levels[cr_max_level];
        cr_cell->detach();
        cr_create_at_level(cr_cell - cr_cells, dest_level);
      }
      cr_max_level--;
    }
  cr_bt_info.resize(btpoint);
}


void
Partition::cr_create_at_level(const unsigned int cell_index,
                              const unsigned int level)
{
  assert(cr_enabled);
  assert(cell_index < N);
  assert(level < N);
  CRCell& cr_cell = cr_cells[cell_index];
  assert(cr_cell.level == UINT_MAX);
  assert(cr_cell.next == 0);
  assert(cr_cell.prev_next_ptr == 0);
  if(cr_levels[level])
    cr_levels[level]->prev_next_ptr = &(cr_cell.next);
  cr_cell.next = cr_levels[level];
  cr_levels[level] = &cr_cell;
  cr_cell.prev_next_ptr = &cr_levels[level];
  cr_cell.level = level;
}


void
Partition::cr_create_at_level_trailed(const unsigned int cell_index,
                                      const unsigned int level)
{
  assert(cr_enabled);
  cr_create_at_level(cell_index, level);
  cr_created_trail.push_back(cell_index);
}


} // namespace bliss
