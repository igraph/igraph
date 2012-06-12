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

#ifndef BLISS_PARTITION_HH
#define BLISS_PARTITION_HH

namespace igraph {

class Cell;
class Partition;

}

#include <cstdlib>
#include "bliss_kstack.hh"
#include "bliss_kqueue.hh"
#include "bliss_heap.hh"
#include "bliss_orbit.hh"
#include "bliss_graph.hh"

namespace igraph {

class Cell
{
public:
  /* Index of the first element of the cell in the Partition::elements array */
  unsigned int first;
  unsigned int length;
  unsigned int max_ival;
  unsigned int max_ival_count;
  bool in_neighbour_heap;
  bool in_splitting_queue;
  Cell *next;
  Cell **prev_next_ptr;
  Cell *next_nonsingleton;
  Cell *prev_nonsingleton;
  unsigned int split_level;
};


class Partition
{
public:
  AbstractGraph *graph;

  /* Used during equitable partition refinement */
  KQueue<Cell *> splitting_queue;
  void add_in_splitting_queue(Cell * const cell);
  void clear_splitting_queue();

  class RefInfo {
  public:
    unsigned int split_cell_first;
    int prev_nonsingleton_first;
    int next_nonsingleton_first;
  };

  /* Used for unrefinement */
  KStack<RefInfo> refinement_stack;

  Cell *aux_split_in_two(Cell * const cell,
			 const unsigned int first_half_size);

  /* The current search level */
  unsigned int level;

  void unrefine(unsigned int dest_level, unsigned int dest_split_stack_size);

  void consistency_check();
public:
  Cell *cells;
  Cell *free_cells;
  Cell *first_cell;
  Cell *first_nonsingleton_cell;
  unsigned int *elements;
  unsigned int *invariant_values;
  /* element_to_cell_map[e] gives the cell of element e */
  Cell **element_to_cell_map;
  /* in_pos[e] points to the elements array s.t. *in_pos[e] = e  */
  unsigned int **in_pos;

  Partition();
  ~Partition();
  void init(const unsigned int);

  bool is_discrete() const {return(free_cells == 0); }

  /*
   * Splits "cell" into [cell_1,...,cell_n] so that &cell_1 == &cell
   *  according to the invariant_values of elements in the cell
   * Returns cell_n which is different from "cell" iff
   *   the cell was actually splitted.
   * max_ival_info_ok indicates that the max_ival and max_ival_count fields
   *  in cell have proper values
   * Clears the invariant values of elements in cell and cell's max_ival and
   *   max_ival_count fields
   */
  Cell *zplit_cell(Cell * const cell, const bool max_ival_info_ok);

private:
  /* Auxiliary routines for sorting and splitting cells */
  void clear_ivs(Cell * const cell);
  Cell *sort_and_split_cell1(Cell *cell);
  Cell *sort_and_split_cell255(Cell * const cell, const unsigned int max_ival);
  bool shellsort_cell(Cell *cell);
  Cell *split_cell(Cell * const cell);
};

}

#endif
