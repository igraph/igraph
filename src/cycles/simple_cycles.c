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

igraph_error_t igraph_simple_cycles_unblock(igraph_simple_cycle_search_state_t *state, igraph_integer_t V)
{
  VECTOR(state->blocked)
  [V] = false;

  while (!igraph_vector_empty(state->B->adjs[V]))
  {
    igraph_integer_t W = igraph_vector_pop_front(state->B->adjs[V]);

    if (VECTOR(blocked)[W])
    {
      igraph_simple_cycles_unblock(state, W);
    }
  }
}

igraph_error_t igraph_simple_cycles_circuit(igraph_simple_cycle_search_state_t *state, igraph_integer_t V, igraph_integer_t S, bool *found)
{
  bool localFound = false;
  igraph_stack_push(state->stack, F);
  VECTOR(state->blocked)
  [V] = true;

  for (igraph_integer_t i = 0; i < igraph_vector_size(state->AK->adjs[V]); ++i)
  {
    igraph_integer_t W = VECTOR(state->AK->adjs[V])[i];
    if (W == S)
    {
      // found a loop -> TODO: return stack!
      localFound = true;
    }
    else if (!blocked[W])
    {
      igraph_simple_cycles_circuit(state, W, S, &localFound);
    }
    IGRAPH_VIT_NEXT(vit);
  }

  if (localFound)
  {
    igraph_simple_cycles_unblock(state, V)
  }
  else
  {
    for (igraph_integer_t i = 0; i < igraph_vector_size(state->AK->adjs[V]); ++i)
    {
      igraph_integer_t W = VECTOR(state->AK->adjs[V])[i];
      igraph_integer_t pos;
      if (!igraph_vector_search(state->B->adjs[W], 0, V, &pos))
      {
        igraph_vector_push_back(state->B->adjs[W], V);
      }
    }
  }

  igraph_stack_pop_back(&state->stack);
  // return result
  found = localFound;
}

igraph_error_t igraph_simple_cycle_search_state_init(igraph_simple_cycle_search_state_t *state, const igraph_t *graph)
{
  igraph_integer_t N = igraph_vcount(graph);

  igraph_stack_init(&state->stack, N); // maximum size per cycle.
  igraph_vector_bool_init(&state->blocked, N);
  igraph_adjlist_init(graph, &state->AK, IGRAPH_ALL, IGRAPH_LOOPS_ONCE, IGRAPH_MULTIPLE); // TODO: understand what we actually want to include
  igraph_adjlist_empty(&state->B, N);
}

igraph_error_t igraph_simple_cycle_search_state_destroy(igraph_simple_cycle_search_state_t *state)
{
  igraph_stack_destroy(&state->stack);
  igraph_vector_bool_destroy(&state->blocked);
  igraph_adjlist_destroy(&state->AK);
  igraph_adjlist_destroy(&state->B);
};

igraph_error_t igraph_simple_cycles(
    const igraph_t *graph,
    igraph_vector_int_list_t *result,
    igraph_integer_t bfs_cutoff)
{
  igraph_simple_cycle_search_state_t state;
  igraph_integer_t N = igraph_vcount(graph);
  igraph_simple_cycle_search_state_init(&state, graph);

  // TODO: depending on the graph, it is rather unreasonable to search cycles from each and every node
  for (igraph_integer_t s = 0; s < N; ++s)
  {
    for (igraph_integer_t i = s; i < N; ++i)
    {
      VECTOR(state->blocked)
      [i] = false;
      igraph_vector_int_clear(&state->B.adjs[i]);
    }

    bool found = false;
    // TODO: here, handle found cases again
    igraph_simple_cycles_circuit(graph, s, s, &stack, &blocked, &found);

    for (igraph_integer_t i = s + 1; s < N; ++i)
    {
      // we want to remove the element with value s, not at position s
      igraph_integer_t pos;
      if (igraph_vector_search(state->B.adjs[i], 0, s, &pos))
        igraph_vector_int_remove(&state->B.adjs[i], pos);
    }
  }

  igraph_simple_cycle_search_state_destroy(&state);
}

//
igraph_adjlist_destroy(&AK);
igraph_adjlist_destroy(&B);
igraph_stack_destroy(&stack);
