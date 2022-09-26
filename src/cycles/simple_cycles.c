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

#define IGRAPH_CYCLE_FOUND 0
#define IGRAPH_ERROR_NO_CYCLE_FOUND 3

// Johnson's cycle detection algorithm

igraph_error_t igraph_simple_cycles_unblock(igraph_simple_cycle_search_state_t *state, igraph_integer_t V)
{
  VECTOR(state->blocked)
  [V] = false;

  while (!igraph_vector_int_empty(&state->B.adjs[V]))
  {
    igraph_integer_t W = igraph_vector_int_pop_front(state->B.adjs[V]);

    if (VECTOR(state->blocked)[W])
    {
      igraph_simple_cycles_unblock(state, W);
    }
  }
}

igraph_error_t igraph_simple_cycles_circuit(igraph_simple_cycle_search_state_t *state, igraph_integer_t V, igraph_integer_t S, bool *found)
{
  bool localFound = false;
  igraph_stack_int_push(&state->stack, V);
  VECTOR(state->blocked)
  [V] = true;

  for (igraph_integer_t i = 0; i < igraph_vector_int_size(&state->AK.adjs[V]); ++i)
  {
    igraph_integer_t W = VECTOR(state->AK.adjs[V])[i];
    if (W == S)
    {
      // found a loop -> TODO: return stack!
      localFound = true;
    }
    else if (!(VECTOR(state->blocked)[W]))
    {
      igraph_simple_cycles_circuit(state, W, S, &localFound);
    }
  }

  if (localFound)
  {
    igraph_simple_cycles_unblock(state, V);
  }
  else
  {
    for (igraph_integer_t i = 0; i < igraph_vector_int_size(&state->AK.adjs[V]); ++i)
    {
      igraph_integer_t W = VECTOR(state->AK.adjs[V])[i];
      igraph_integer_t pos;
      if (!igraph_vector_int_search(&state->B.adjs[W], 0, V, &pos))
      {
        igraph_vector_int_push_back(&state->B.adjs[W], V);
      }
    }
  }

  igraph_stack_int_pop(&state->stack); // _back
  // return result
  *found = localFound;
}

igraph_error_t igraph_simple_cycle_search_state_init(igraph_simple_cycle_search_state_t *state, const igraph_t *graph)
{
  igraph_integer_t N = igraph_vcount(graph);

  state->N = N;
  igraph_stack_int_init(&state->stack, N); // maximum size per cycle.
  igraph_vector_bool_init(&state->blocked, N);
  igraph_adjlist_init(graph, &state->AK, IGRAPH_ALL, IGRAPH_LOOPS_ONCE, IGRAPH_MULTIPLE); // TODO: understand what we actually want to include
  igraph_adjlist_init_empty(&state->B, N);
}

igraph_error_t igraph_simple_cycle_search_state_destroy(igraph_simple_cycle_search_state_t *state)
{
  igraph_stack_int_destroy(&state->stack);
  igraph_vector_bool_destroy(&state->blocked);
  igraph_adjlist_destroy(&state->AK);
  igraph_adjlist_destroy(&state->B);
}

igraph_error_t igraph_simple_cycles_search_one(
    igraph_simple_cycle_search_state_t *state, igraph_integer_t s, igraph_vector_int_t *result)
{
  for (igraph_integer_t i = s; i < state->N; ++i)
  {
    VECTOR(state->blocked)
    [i] = false;
    igraph_vector_int_clear(&state->B.adjs[i]);
  }

  bool found = false;
  igraph_simple_cycles_circuit(state, s, s, &found);

  for (igraph_integer_t i = s + 1; s < state->N; ++i)
  {
    // we want to remove the element with value s, not at position s
    igraph_integer_t pos;
    if (igraph_vector_int_search(&state->B.adjs[i], 0, s, &pos))
      igraph_vector_int_remove(&state->B.adjs[i], pos);
  }

  if (found)
  {
    // return stack // TODO: currently, only the nodes are returned
    igraph_vector_int_resize(result, igraph_stack_int_size(&state->stack));
    for (igraph_integer_t i = igraph_stack_int_size(&state->stack) - 1; i >= 0; --i)
    {
      VECTOR(*result)
      [i] = igraph_stack_int_pop(&state->stack); // igraph_stack_int_get(&state->stack, i);
    }
    return IGRAPH_CYCLE_FOUND;
  }
  else
  {
    return IGRAPH_ERROR_NO_CYCLE_FOUND;
  }
}

igraph_error_t igraph_simple_cycles_search_all(
    const igraph_t *graph,
    igraph_vector_int_list_t *result)
{
  igraph_simple_cycle_search_state_t state;
  igraph_simple_cycle_search_state_init(&state, graph);

  // int nfound = 0;
  igraph_vector_int_list_init(result, 0); // state.N);

  // TODO: depending on the graph, it is rather unreasonable to search cycles from each and every node
  for (igraph_integer_t s = 0; s < state.N; ++s)
  {
    igraph_vector_int_t res;
    igraph_vector_int_init(&res, 0);
    if (igraph_simple_cycles_search_one(&state, s, &res) == IGRAPH_CYCLE_FOUND)
    {
      // nfound += 1;
      igraph_vector_int_list_push_back_copy(result, &res);
      igraph_vector_int_empty(&res);
    }
  }

  // igraph_vector_int_list_resize(result, nfound);
  igraph_simple_cycle_search_state_destroy(&state);
}
