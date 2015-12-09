#ifndef BLISS_UINTSEQHASH_HH
#define BLISS_UINTSEQHASH_HH

#include <cstdio>

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
 * \brief A hash for sequences of unsigned ints.
 */
class UintSeqHash
{
protected:
  unsigned int h;
public:
  UintSeqHash() {h = 0; }
  UintSeqHash(const UintSeqHash &other) {h = other.h; }
  UintSeqHash& operator=(const UintSeqHash &other) {h = other.h; return *this; }
  
  /** Reset the hash value. */
  void reset() {h = 0; }

  /** Add the unsigned int \a n to the sequence. */
  void update(unsigned int n);

  /** Get the hash value of the sequence seen so far. */
  unsigned int get_value() const {return h; }

  /** Compare the hash values of this and \a other.
   * Return -1/0/1 if the value of this is smaller/equal/greater than
   * that of \a other. */
  int cmp(const UintSeqHash &other) const {
    return (h < other.h)?-1:((h == other.h)?0:1);
  }
  /** An abbreviation for cmp(other) < 0 */
  bool is_lt(const UintSeqHash &other) const {return(cmp(other) < 0); }
  /** An abbreviation for cmp(other) <= 0 */
  bool is_le(const UintSeqHash &other) const {return(cmp(other) <= 0); }
  /** An abbreviation for cmp(other) == 0 */
  bool is_equal(const UintSeqHash &other) const {return(cmp(other) == 0); }
};


} // namespace bliss

#endif
