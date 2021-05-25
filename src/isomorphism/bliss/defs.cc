#include "igraph_error.h"

#include "defs.hh"

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

void
fatal_error(const char* reason)
{
  IGRAPH_FATAL(reason);
}

}
