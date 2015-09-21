#ifndef BLISS_ORBIT_HH
#define BLISS_ORBIT_HH

/*
  Copyright (c) 2003-2015 Tommi Junttila
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

/** \internal
 * \brief A class for representing orbit information.
 *
 * Given a set {0,...,N-1} of N elements, represent equivalence
 * classes (that is, unordered partitions) of the elements.
 * Supports only equivalence class merging, not splitting.
 * Merging two classes requires time O(k), where k is the number of
 * the elements in the smaller of the merged classes.
 * Getting the smallest representative in a class (and thus testing
 * whether two elements belong to the same class) is a constant time operation.
 */
class Orbit
{
  class OrbitEntry
  {
  public:
    unsigned int element;
    OrbitEntry *next;
    unsigned int size;
  };

  OrbitEntry *orbits;
  OrbitEntry **in_orbit;
  unsigned int nof_elements;
  unsigned int _nof_orbits;
  void merge_orbits(OrbitEntry *o1, OrbitEntry *o2);

public:
  /**
   * Create a new orbit information object.
   * The init() function must be called next to actually initialize
   * the object.
   */
  Orbit();
  ~Orbit();

  /**
   * Initialize the orbit information to consider sets of \a N elements.
   * It is required that \a N > 0.
   * The orbit information is reset so that each element forms
   * an orbit of its own.
   * Time complexity is O(N).
   * \sa reset()
   */
  void init(const unsigned int N);

  /**
   * Reset the orbits so that each element forms an orbit of its own.
   * Time complexity is O(N).
   */
  void reset();

  /**
   * Merge the orbits of the elements \a e1 and \a e2.
   * Time complexity is O(k), where k is the number of elements in
   * the smaller of the merged orbits.
   */
  void merge_orbits(unsigned int e1, unsigned int e2);

  /**
   * Is the element \a e the smallest element in its orbit?
   * Time complexity is O(1).
   */
  bool is_minimal_representative(unsigned int e) const;

  /**
   * Get the smallest element in the orbit of the element \a e.
   * Time complexity is O(1).
   */
  unsigned int get_minimal_representative(unsigned int e) const;

  /**
   * Get the number of elements in the orbit of the element \a e.
   * Time complexity is O(1).
   */
  unsigned int orbit_size(unsigned int e) const;

  /**
   * Get the number of orbits.
   * Time complexity is O(1).
   */
  unsigned int nof_orbits() const {return _nof_orbits; }
};

} // namespace bliss

#endif
