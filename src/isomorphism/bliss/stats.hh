#ifndef BLISS_STATS_HH
#define BLISS_STATS_HH

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

#include "graph.hh"
#include "bignum.hh"

namespace bliss {

/**
 * \brief Statistics returned by the bliss search algorithm.
 */
class Stats
{
  friend class AbstractGraph;
  /** \internal The size of the automorphism group. */
  BigNum group_size;
  /** \internal An approximation (due to possible overflows) of
   * the size of the automorphism group. */
  long double group_size_approx;
  /** \internal The number of nodes in the search tree. */
  long unsigned int nof_nodes;
  /** \internal The number of leaf nodes in the search tree. */
  long unsigned int nof_leaf_nodes;
  /** \internal The number of bad nodes in the search tree. */
  long unsigned int nof_bad_nodes;
  /** \internal The number of canonical representative updates. */
  long unsigned int nof_canupdates;
  /** \internal The number of generator permutations. */
  long unsigned int nof_generators;
  /** \internal The maximal depth of the search tree. */
  unsigned long int max_level;
  /** \internal Reset the statistics. */
  void reset()
  {
    group_size.assign(1);
    group_size_approx = 1.0;
    nof_nodes = 0;
    nof_leaf_nodes = 0;
    nof_bad_nodes = 0;
    nof_canupdates = 0;
    nof_generators = 0;
    max_level = 0;
  }
public:
  Stats() { reset(); }

  /** The size of the automorphism group. */
  const BigNum& get_group_size() const {return group_size;}
  /** An approximation (due to possible overflows/rounding errors) of
   * the size of the automorphism group. */
  long double get_group_size_approx() const {return group_size_approx;}
  /** The number of nodes in the search tree. */
  long unsigned int get_nof_nodes() const {return nof_nodes;}
  /** The number of leaf nodes in the search tree. */
  long unsigned int get_nof_leaf_nodes() const {return nof_leaf_nodes;}
  /** The number of bad nodes in the search tree. */
  long unsigned int get_nof_bad_nodes() const {return nof_bad_nodes;}
  /** The number of canonical representative updates. */
  long unsigned int get_nof_canupdates() const {return nof_canupdates;}
  /** The number of generator permutations. */
  long unsigned int get_nof_generators() const {return nof_generators;}
  /** The maximal depth of the search tree. */
  unsigned long int get_max_level() const {return max_level;}
};

} // namespace bliss

#endif // BLISS_STATS_HH
