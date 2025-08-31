/*
  igraph library.
  Copyright (C) 2025 The igraph development team <igraph@igraph.org>

  This program is free software; you can redistribute it and/or modify it under
  the terms of the GNU General Public License as published by the Free Software
  Foundation; either version 2 of the License, or (at your option) any later
  version.

  This program is distributed in the hope that it will be useful, but WITHOUT
  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
  FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

  You should have received a copy of the GNU General Public License along with
  this program. If not, see <https://www.gnu.org/licenses/>.
*/

#include "random/random_internal.h"

#include <ctime>
#include <random>

/* This function attempts to produce an unpredictable number suitable for
 * seeding the RNG upon startup. It makes an effort to avoid producing the
 * same number when called within different processes, even if called
 * approximately at the same time. However, it cannot guarantee this on
 * all systems under all circumstancs.
 *
 * This function cannot fail.
 */
igraph_uint_t igraph_i_get_random_seed(void) {
    try {
        // Try to use C++'s std::random_device. This may fail, either because
        // the systems's random device ran out of entropy or because the system
        // does not offer this functionality.
        return std::random_device()();
    } catch (...) {
        // If random_device fails, use a time-based seed. A combination of
        // time() (calendar time, low resolution) and clock() (processor time
        // used, higher resolution) reduces the risk of seed collisions.
        return igraph_uint_t(std::clock()) + igraph_uint_t(std::time(NULL));
    }
}
