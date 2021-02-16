#ifndef BLISS_UTILS_HH
#define BLISS_UTILS_HH

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

/**
 * \file
 * \brief Some small utilities.
 */

#include <vector>

namespace bliss {

/**
 * Check whether \a perm is a valid permutation on {0,...,N-1}.
 * Slow, mainly for debugging and validation purposes.
 */
bool is_permutation(const unsigned int N, const unsigned int* perm);

/**
 * Check whether \a perm is a valid permutation on {0,...,N-1}.
 * Slow, mainly for debugging and validation purposes.
 */
bool is_permutation(const std::vector<unsigned int>& perm);

} // namespace bliss

#endif // BLISS_UTILS_HH
