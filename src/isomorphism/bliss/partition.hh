#ifndef BLISS_PARTITION_HH
#define BLISS_PARTITION_HH

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
  class Partition;
}

#include <vector>
#include <climits>
#include "kstack.hh"
#include "kqueue.hh"
#include "graph.hh"


namespace bliss {

/**
 * \brief A class for refinable, backtrackable ordered partitions.
 *
 * This is rather a data structure with some helper functions than
 * a proper self-contained class.
 * That is, for efficiency reasons the fields of this class are directly
 * manipulated from bliss::AbstractGraph and its subclasses.
 * Conversely, some methods of this class modify the fields of
 * bliss::AbstractGraph, too.
 */
class Partition
{
public:
  /**
   * \brief Data structure for holding information about a cell in a Partition.
   */
  class Cell
  {
    friend class Partition;
  public:
    unsigned int length;
    /* Index of the first element of the cell in
       the Partition::elements array */
    unsigned int first;
    unsigned int max_ival;
    unsigned int max_ival_count;
  private:
    bool in_splitting_queue;
  public:
    bool in_neighbour_heap;
    /* Pointer to the next cell, null if this is the last one. */
    Cell* next;
    Cell* prev;
    Cell* next_nonsingleton;
    Cell* prev_nonsingleton;
    unsigned int split_level;
    /** Is this a unit cell? */
    bool is_unit() const {return(length == 1); }
    /** Is this cell in splitting queue? */
    bool is_in_splitting_queue() const {return(in_splitting_queue); }
  };


private:

  /** \internal
   * Data structure for remembering information about splits in order to
   * perform efficient backtracking over the splits.
   */
  class RefInfo {
  public:
    unsigned int split_cell_first;
    int prev_nonsingleton_first;
    int next_nonsingleton_first;
  };
  /** \internal
   * A stack for remembering the splits, used for backtracking.
   */
  KStack<RefInfo> refinement_stack;

  class BacktrackInfo {
  public:
    unsigned int refinement_stack_size;
    unsigned int cr_backtrack_point;
  };

  /** \internal
   * The main stack for enabling backtracking.
   */
  std::vector<BacktrackInfo> bt_stack;

public:
  AbstractGraph* graph;

  /* Used during equitable partition refinement */
  KQueue<Cell*> splitting_queue;
  void  splitting_queue_add(Cell* const cell);
  Cell* splitting_queue_pop();
  bool  splitting_queue_is_empty() const;
  void  splitting_queue_clear();


  /** Type for backtracking points. */
  typedef unsigned int BacktrackPoint;

  /**
   * Get a new backtrack point for the current partition
   */
  BacktrackPoint set_backtrack_point();

  /**
   * Backtrack to the point \a p and remove it.
   */
  void goto_backtrack_point(BacktrackPoint p);

  /**
   * Split the non-unit Cell \a cell = {\a element,e1,e2,...,en} containing
   * the element \a element in two:
   * \a cell = {e1,...,en} and \a newcell = {\a element}.
   * @param cell     a non-unit Cell
   * @param element  an element in \a cell
   * @return         the new unit Cell \a newcell
   */
  Cell* individualize(Cell* const cell,
                      const unsigned int element);

  Cell* aux_split_in_two(Cell* const cell,
                         const unsigned int first_half_size);


private:
  unsigned int N;
  Cell* cells;
  Cell* free_cells;
  unsigned int discrete_cell_count;
public:
  Cell* first_cell;
  Cell* first_nonsingleton_cell;
  unsigned int *elements;
  /* invariant_values[e] gives the invariant value of the element e */
  unsigned int *invariant_values;
  /* element_to_cell_map[e] gives the cell of the element e */
  Cell **element_to_cell_map;
  /** Get the cell of the element \a e */
  Cell* get_cell(const unsigned int e) const {
    assert(e < N);
    return element_to_cell_map[e];
  }
  /* in_pos[e] points to the elements array s.t. *in_pos[e] = e  */
  unsigned int **in_pos;

  Partition();
  ~Partition();

  /**
   * Initialize the partition to the unit partition (all elements in one cell)
   * over the \a N > 0 elements {0,...,\a N-1}.
   */
  void init(const unsigned int N);

  /**
   * Returns true iff the partition is discrete, meaning that all
   * the elements are in their own cells.
   */
  bool is_discrete() const {return(free_cells == 0); }

  unsigned int nof_discrete_cells() const {return(discrete_cell_count); }

  /*
   * Splits the Cell \a cell into [cell_1,...,cell_n]
   * according to the invariant_values of the elements in \a cell.
   * After splitting, cell_1 == \a cell.
   * Returns the pointer to the Cell cell_n;
   * cell_n != cell iff the Cell \a cell was actually splitted.
   * The flag \a max_ival_info_ok indicates whether the max_ival and
   * max_ival_count fields of the Cell \a cell have consistent values
   * when the method is called.
   * Clears the invariant values of elements in the Cell \a cell as well as
   * the max_ival and max_ival_count fields of the Cell \a cell.
   */
  Cell *zplit_cell(Cell * const cell, const bool max_ival_info_ok);

  /*
   * Routines for component recursion
   */
  void cr_init();
  void cr_free();
  unsigned int cr_get_level(const unsigned int cell_index) const;
  unsigned int cr_split_level(const unsigned int level,
                              const std::vector<unsigned int>& cells);

  /** Clear the invariant_values of the elements in the Cell \a cell. */
  void clear_ivs(Cell* const cell);

private:
  /*
   * Component recursion data structures
   */

  /* Is component recursion support in use? */
  bool cr_enabled;

  class CRCell {
  public:
    unsigned int level;
    CRCell* next;
    CRCell** prev_next_ptr;
    void detach() {
      if(next)
        next->prev_next_ptr = prev_next_ptr;
      *(prev_next_ptr) = next;
      level = UINT_MAX;
      next = 0;
      prev_next_ptr = 0;
    }
  };
  CRCell* cr_cells;
  CRCell** cr_levels;
  class CR_BTInfo {
  public:
    unsigned int created_trail_index;
    unsigned int splitted_level_trail_index;
  };
  std::vector<unsigned int> cr_created_trail;
  std::vector<unsigned int> cr_splitted_level_trail;
  std::vector<CR_BTInfo> cr_bt_info;
  unsigned int cr_max_level;
  void cr_create_at_level(const unsigned int cell_index, unsigned int level);
  void cr_create_at_level_trailed(const unsigned int cell_index, unsigned int level);
  unsigned int cr_get_backtrack_point();
  void cr_goto_backtrack_point(const unsigned int btpoint);


  /*
   *
   * Auxiliary routines for sorting and splitting cells
   *
   */
  Cell* sort_and_split_cell1(Cell* cell);
  Cell* sort_and_split_cell255(Cell* const cell, const unsigned int max_ival);
  bool shellsort_cell(Cell* cell);
  Cell* split_cell(Cell* const cell);

  /*
   * Some auxiliary stuff needed for distribution count sorting.
   * To make the code thread-safe (modulo the requirement that each graph is
   * only accessed in one thread at a time), the arrays are owned by
   * the partition instance, not statically defined.
   */
  unsigned int dcs_count[256];
  unsigned int dcs_start[256];
  void dcs_cumulate_count(const unsigned int max);
};


inline Partition::Cell*
Partition::splitting_queue_pop()
{
  assert(!splitting_queue.is_empty());
  Cell* const cell = splitting_queue.pop_front();
  assert(cell->in_splitting_queue);
  cell->in_splitting_queue = false;
  return cell;
}

inline bool
Partition::splitting_queue_is_empty() const
{
  return splitting_queue.is_empty();
}


inline unsigned int
Partition::cr_get_level(const unsigned int cell_index) const
{
  assert(cr_enabled);
  assert(cell_index < N);
  assert(cr_cells[cell_index].level != UINT_MAX);
  return(cr_cells[cell_index].level);
}

} // namespace bliss

#endif // BLISS_PARTITION_HH
