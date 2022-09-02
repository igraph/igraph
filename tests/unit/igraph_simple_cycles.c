/*
   IGraph library.
   Copyright (C) 2021  The igraph development team <igraph@igraph.org>

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

#include <igraph.h>
#include "test_utilities.h"

int main()
{
  igraph_t g_cycle, g_no_cycle;

  igraph_star(&g_no_cycle, 7, IGRAPH_STAR_UNDIRECTED, 1);
  // TODO: call cycles finder, expect 0 cycle to be found

  igraph_ring(&g_cycle, 10, false, true, true);
  // TODO: call cycles finder, expect 1 cycle to be found
}
