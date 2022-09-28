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

#include <stdlib.h>
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

  while (!igraph_vector_int_empty(igraph_adjlist_get(&state->B, V)))
  {
    // NOTE: pop front is highly inefficient. Check frequency of call, possibly replace adjlist
    // probably with a vector of lists
    igraph_integer_t W = igraph_vector_int_pop_front(igraph_adjlist_get(&state->B, V));

    if (VECTOR(state->blocked)[W])
    {
      igraph_simple_cycles_unblock(state, W);
    }
  }

  return IGRAPH_SUCCESS;
}

igraph_error_t igraph_simple_cycles_circuit(igraph_simple_cycle_search_state_t *state, igraph_integer_t V, igraph_integer_t S, igraph_vector_int_list_t *results, bool *found)
{
  bool localFound = false;
  igraph_stack_int_push(&state->stack, V);
  printf("Pushing to stack, stack size is %lld\n", igraph_stack_int_size(&state->stack));
  VECTOR(state->blocked)
  [V] = true;

  for (igraph_integer_t i = 0; i < igraph_vector_int_size(igraph_adjlist_get(&state->AK, V)); ++i)
  {
    igraph_integer_t W = VECTOR(*igraph_adjlist_get(&state->AK, V))[i];
    if (W == S)
    {
      localFound = true;
      // copy output: from stack to vector
      igraph_integer_t i = 0;
      igraph_vector_int_t res;
      igraph_vector_int_init(&res, igraph_stack_int_size(&state->stack));
      while (!igraph_stack_int_empty(&state->stack))
      {
        printf("Setting result value %lld\n", i);
        VECTOR(res)
        [i] = igraph_stack_int_pop(&state->stack); // igraph_stack_int_get(&state->stack, i);
        printf("Set result value %lld\n", VECTOR(res)
                                              [i]);
        i += 1;
      }
      igraph_vector_int_list_push_back_copy(results, &res);
      igraph_vector_int_destroy(&res);
    }
    else if (W > S && !(VECTOR(state->blocked)[W]))
    {
      igraph_simple_cycles_circuit(state, W, S, results, &localFound);
    }
  }
  *found = localFound;

  if (localFound)
  {
    igraph_simple_cycles_unblock(state, V);
  }
  else
  {
    for (igraph_integer_t i = 0; i < igraph_vector_int_size(igraph_adjlist_get(&state->AK, V)); ++i)
    {
      igraph_integer_t W = VECTOR(*igraph_adjlist_get(&state->AK, V))[i];
      igraph_integer_t pos;
      if (!igraph_vector_int_search(igraph_adjlist_get(&state->B, W), 0, V, &pos))
      {
        igraph_vector_int_push_back(igraph_adjlist_get(&state->B, W), V);
      }
    }
  }

  if (igraph_stack_int_size(&state->stack) > 0)
  {
    igraph_stack_int_pop(&state->stack); // _back
  }
  // return result

  return IGRAPH_SUCCESS;
}

igraph_error_t igraph_simple_cycle_search_state_init(igraph_simple_cycle_search_state_t *state, const igraph_t *graph)
{
  igraph_integer_t N = igraph_vcount(graph);

  state->N = N;
  igraph_stack_int_init(&state->stack, N); // maximum size per cycle.
  igraph_vector_bool_init(&state->blocked, N);
  igraph_adjlist_init(graph, &state->AK, IGRAPH_ALL, IGRAPH_LOOPS_ONCE, IGRAPH_MULTIPLE); // TODO: understand what we actually want to include
  igraph_adjlist_init_empty(&state->B, N);

  return IGRAPH_SUCCESS;
}

igraph_error_t igraph_simple_cycle_search_state_destroy(igraph_simple_cycle_search_state_t *state)
{
  igraph_stack_int_destroy(&state->stack);
  igraph_vector_bool_destroy(&state->blocked);
  igraph_adjlist_destroy(&state->AK);
  igraph_adjlist_destroy(&state->B);

  return IGRAPH_SUCCESS;
}

igraph_error_t igraph_simple_cycles_search_one(
    igraph_simple_cycle_search_state_t *state, igraph_integer_t s, igraph_vector_int_list_t *results)
{
  for (igraph_integer_t i = s; i < state->N; ++i)
  {
    VECTOR(state->blocked)
    [i] = false;
    igraph_vector_int_clear(igraph_adjlist_get(&state->B, i));
  }

  bool found = false;
  igraph_simple_cycles_circuit(state, s, s, results, &found);

  // for (igraph_integer_t i = s; i < state->N; ++i)
  // {
  //   // we want to remove the vertex with value s, not at position s
  //   igraph_integer_t pos;
  //   if (igraph_vector_int_search(igraph_adjlist_get(&state->AK, i), 0, s, &pos))
  //   {
  //     igraph_vector_int_remove(igraph_adjlist_get(&state->AK, i), pos);
  //   }
  // }

  if (found)
  {
    // TODO: currently, only the nodes are returned
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
    igraph_vector_int_list_t res;
    igraph_vector_int_list_init(&res, 0);
    if (igraph_simple_cycles_search_one(&state, s, &res) == IGRAPH_CYCLE_FOUND)
    {
      // nfound += 1;
      printf("Found %lld. adding to all results.\n", igraph_vector_int_list_size(&res));
      for (igraph_integer_t i = 0; i < igraph_vector_int_list_size(&res); i++)
      {
        igraph_vector_int_list_push_back_copy(result, igraph_vector_int_list_get_ptr(&res, i));
      }
    }
    igraph_vector_int_list_clear(&res);
  }

  // igraph_vector_int_list_resize(result, nfound);
  igraph_simple_cycle_search_state_destroy(&state);

  return IGRAPH_SUCCESS;
}
