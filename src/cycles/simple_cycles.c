/*
   IGraph library.
   Copyright (C) 2021-2022  The igraph development team <igraph@igraph.org>

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

#include "igraph_cycles.h"

#include "igraph_adjlist.h"
#include "igraph_components.h"
#include "igraph_dqueue.h"
#include "igraph_error.h"
#include "igraph_interface.h"

#include "core/interruption.h"

#include <algorithm>
#include <iostream>
#include <list>
#include <vector>

// Johnson's cycle detection algorithm


igraph_error_t igraph_simple_cycles_unblock(igraph_integer_t V, igraph_vector_bool_t *blocked, igraph_adjlist_t *B)
{
  VECTOR(blocked)
  [V] = false;

  while (!igraph_vector_empty(B->adjs[V]))
  {
    igraph_integer_t W = igraph_vector_pop_front(B->adjs[V]);

    if (VECTOR(blocked)[W])
    {
      igraph_simple_cycles_unblock(W, blocked, B);
    }
  }
}

igraph_error_t igraph_simple_cycles_circuit(igraph_adjlist_t *AK, igraph_adjlist_t *B, igraph_integer_t V, igraph_integer_t S, igraph_stack_t *stack, igraph_vector_bool_t *blocked, bool *found)
{
  bool localFound = false;
  igraph_stack_push(stack, F);
  VECTOR(blocked)
  [V] = true;

  for (igraph_integer_t i = 0; i < igraph_vector_size(AK->adjs[V]); ++i)
  {
    igraph_integer_t W = VECTOR(AK->adjs[V])[i];
    if (W == S)
    {
      // found a loop -> TODO: return stack!
      localFound = true;
    }
    else if (!blocked[W])
    {
      bool found = false;
      igraph_simple_cycles_circuit(graph, W, S, stack, blocked, &found);
      if (found)
      {
        localFound = true;
      }
    }
    IGRAPH_VIT_NEXT(vit);
  }

  if (localFound)
  {
    igraph_simple_cycles_unblock(V)
  }
  else
  {
    for (igraph_integer_t i = 0; i < igraph_vector_size(AK->adjs[V]); ++i)
    {
      igraph_integer_t W = VECTOR(AK->adjs[V])[i];
      igraph_integer_t pos;
      if (!igraph_vector_search(B->adjs[W], 0, V, &pos))
      {
        igraph_vector_push_back(B->adjs[W], V);
      }
    }
  }

  igraph_stack_pop_back(&stack);
  // return result
  found = localFound;
}

igraph_error_t igraph_simple_cycles(
    const igraph_t *graph,
    igraph_vector_int_list_t *result,
    igraph_integer_t bfs_cutoff)
{
  igraph_integer_t N = igraph_vcount(graph);
  igraph_stack_t stack;
  igraph_stack_init(&stack, N); // maximum size per cycle.
  igraph_vector_bool_t Blocked;
  igraph_vector_bool_init(&Blocked, N);

  igraph_adjlist_t AK;
  igraph_adjlist_init(graph, &AK, IGRAPH_ALL, IGRAPH_LOOPS_ONCE, IGRAPH_MULTIPLE); // TODO: understand what we actually want to include

  igraph_adjlist_t B;
  igraph_adjlist_empty(&B, N);

  // TODO: depending on the graph, it is rather unreasonable to search cycles from each and every node
  for (igraph_integer_t s = 0; s < N; ++s)
  {
    for (igraph_integer_t i = s; i < N; ++i)
    {
      VECTOR(Blocked)
      [i] = false;
      igraph_vector_int_clear(&B.adjs[i]);
    }

    bool found = false;
    // TODO: here, handle found cases again
    igraph_simple_cycles_circuit(graph, s, s, &stack, &blocked, &found);

    for (igraph_integer_t i = s + 1; s < N; ++i)
    {
      // we want to remove the element with value s, not at position s
      igraph_integer_t pos;
      if (igraph_vector_search(B.adjs[i], 0, s, &pos))
        igraph_vector_int_remove(&B.adjs[i], pos);
    }
  }
}

//
igraph_adjlist_destroy(&AK);
igraph_adjlist_destroy(&B);
igraph_stack_destroy(&stack);
}
