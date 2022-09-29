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
#include "igraph_structural.h"
#include "../../tests/unit/test_utilities.h"

#include "core/interruption.h"

#define IGRAPH_CYCLE_FOUND 0
#define IGRAPH_ERROR_NO_CYCLE_FOUND 3

// Johnson's cycle detection algorithm

igraph_error_t igraph_simple_cycles_unblock(igraph_simple_cycle_search_state_t *state, igraph_integer_t u)
{
  VECTOR(state->blocked)
  [u] = false;

  while (!igraph_vector_int_empty(igraph_adjlist_get(&state->B, u)))
  {
    igraph_integer_t w = igraph_vector_int_pop_back(igraph_adjlist_get(&state->B, u));

    if (VECTOR(state->blocked)[w])
    {
      igraph_simple_cycles_unblock(state, w);
    }
  }

  return IGRAPH_SUCCESS;
}

igraph_error_t igraph_simple_cycles_circuit(igraph_simple_cycle_search_state_t *state, igraph_integer_t V, igraph_integer_t S, igraph_vector_int_list_t *results, bool *found, igraph_integer_t depth)
{
  bool local_found = false;
  // stack v
  if (V != S)
  {
    igraph_stack_int_push(&state->stack, V);
  }
  printf("Pushing %lld to stack, stack size is %lld, result size is %lld, depth %lld\n", V, igraph_stack_int_size(&state->stack), igraph_vector_int_list_size(results), depth);
  VECTOR(state->blocked)
  [V] = true;

  // L1
  for (igraph_integer_t i = 0; i < igraph_vector_int_size(igraph_adjlist_get(&state->AK, V)); ++i)
  {
    igraph_integer_t W = VECTOR(*igraph_adjlist_get(&state->AK, V))[i];
    // NOTE: possibly dangerous fix for undirected graphs,
    // disabling finding any two-vertex-loops
    if (W == S && (state->directed || depth > 1)) // or depth > 1
    {
      local_found = true;
      depth = 0;
      // output circuit composed of stack followed by s
      igraph_stack_int_push(&state->stack, S);
      printf("Found cycle with size %lld\n", igraph_stack_int_size(&state->stack));

      // copy output: from stack to vector
      igraph_integer_t res_idx = 0;
      igraph_vector_int_t res;
      igraph_vector_int_init(&res, igraph_stack_int_size(&state->stack));
      while (!igraph_stack_int_empty(&state->stack))
      {
        printf("Setting result value %lld\n", res_idx);
        VECTOR(res)
        [res_idx] = igraph_stack_int_pop(&state->stack); // igraph_stack_int_get(&state->stack, i);
        printf("Set result value %lld\n", VECTOR(res)
                                              [res_idx]);
        res_idx += 1;
      }
      igraph_vector_int_sort(&res);
      // undirected graphs lead to every cycle being found twice.
      // this is our naïve filter for now
      if (!state->directed)
      {
        igraph_bool_t duplicate_found = false;
        for (igraph_integer_t results_idx = 0; results_idx < igraph_vector_int_list_size(results); ++results_idx)
        {
          if (igraph_vector_int_size(igraph_vector_int_list_get_ptr(results, results_idx)) != igraph_vector_int_size(&res))
          {
            continue;
          }
          igraph_bool_t discrepancy_found = false;
          for (igraph_integer_t res_idx = 0; res_idx < igraph_vector_int_size(&res); ++res_idx)
          {
            if (igraph_vector_int_get(&res, res_idx) != igraph_vector_int_get(igraph_vector_int_list_get_ptr(results, results_idx), res_idx))
            {
              discrepancy_found = true;
              break;
            }
          }
          if (!discrepancy_found)
          {
            // found this loop already.
            duplicate_found = true;
            break;
          }
        }
        if (duplicate_found)
        {
          igraph_vector_int_destroy(&res);
          continue;
        }
      }
      // end filter
      igraph_vector_int_list_push_back_copy(results, &res);
      igraph_vector_int_destroy(&res);
    }
    else if (!(VECTOR(state->blocked)[W]))
    {
      igraph_simple_cycles_circuit(state, W, S, results, &local_found, depth + 1);
    }
  }
  *found = local_found;

  // L2
  if (local_found)
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

  if (V != S && igraph_stack_int_size(&state->stack) > 0)
  {
    // unstack v
    printf("Unstacking %lld\n", V);
    igraph_stack_int_pop(&state->stack); // _back
  }
  // return result

  return IGRAPH_SUCCESS;
}

/**
 * @brief Initialize the cycle search state
 *
 * @param state
 * @param graph
 * @return igraph_error_t
 *
 * Time complexity: O(|V|*|E|*log(|V|*|E|))
 */
igraph_error_t igraph_simple_cycle_search_state_init(igraph_simple_cycle_search_state_t *state, const igraph_t *graph)
{
  printf("\n\nInitialization\n");
  igraph_integer_t N = igraph_vcount(graph);

  state->N = N;
  igraph_stack_int_init(&state->stack, N); // maximum size per cycle.
  igraph_vector_bool_init(&state->blocked, N);
  igraph_adjlist_init(graph, &state->AK, IGRAPH_OUT, IGRAPH_LOOPS_ONCE, IGRAPH_MULTIPLE); // TODO: understand what we actually want to include
  // igraph_adjlist_reverse_sort(&state->AK);                                                 // here is where the time complexity comes in
  // // remove all neighbours where the indices do not follow the expected order
  // for (igraph_integer_t i = N - 1; i >= 0; --i)
  // {
  //   igraph_vector_int_t *curr_adjs = igraph_adjlist_get(&state->AK, i);
  //   igraph_integer_t nCurr_adjs = igraph_vector_int_size(curr_adjs);
  //   for (igraph_integer_t j = nCurr_adjs - 1; j >= 0; --j)
  //   {
  //     if (igraph_vector_int_get(curr_adjs, j) < i)
  //     {
  //       igraph_vector_int_push_back(igraph_adjlist_get(&state->AK, igraph_vector_int_pop_back(curr_adjs)), i);
  //     }
  //     else
  //     {
  //       // as they are sorted, we can break out of the loop
  //       break;
  //     }
  //   }
  // }

  igraph_adjlist_sort(&state->AK);
  state->directed = igraph_is_directed(graph);
  // if (!igraph_is_directed(graph))
  // {
  //   // filter out all back-and-forth edges
  //   for (igraph_integer_t i = 0; i < N; ++i)
  //   {
  //     igraph_vector_int_t *neighs = igraph_adjlist_get(&state->AK, i);
  //     igraph_integer_t nCurr_neighs = igraph_vector_int_size(neighs);
  //     igraph_integer_t pos;
  //     for (igraph_integer_t j = nCurr_neighs - 1; j >= 0; --j)
  //     {
  //       igraph_integer_t neigh_val = igraph_vector_int_get(neighs, j);
  //       if (neigh_val < i)
  //       {
  //         // as they are sorted, we can break out of the loop
  //         break;
  //       }
  //       if (igraph_vector_int_search(igraph_adjlist_get(&state->AK, neigh_val), 0, i, &pos))
  //       {
  //         igraph_vector_int_remove(igraph_adjlist_get(&state->AK, neigh_val), pos);
  //       }
  //     }
  //   }
  // }

  for (igraph_integer_t j = 0; j < igraph_adjlist_size(&state->AK); ++j)
  {
    print_vector_int(igraph_adjlist_get(&state->AK, j));
  }
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
  // L3:
  for (igraph_integer_t i = 0; i < state->N; ++i)
  {
    VECTOR(state->blocked)
    [i] = false;
    igraph_vector_int_clear(igraph_adjlist_get(&state->B, i));
  }

  bool found = false;
  igraph_simple_cycles_circuit(state, s, s, results, &found, 0);

  for (igraph_integer_t i = 0; i < state->N; ++i)
  {
    // we want to remove the vertex with value s, not at position s
    igraph_integer_t pos;
    if (igraph_vector_int_search(igraph_adjlist_get(&state->AK, i), 0, s, &pos))
    {
      igraph_vector_int_remove(igraph_adjlist_get(&state->AK, i), pos);
    }
  }
  igraph_vector_int_clear(igraph_adjlist_get(&state->AK, s));

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

/**
 * @see https://en.wikipedia.org/wiki/Johnson%27s_algorithm
 * @see https://stackoverflow.com/a/35922906/3909202
 * @see https://epubs.siam.org/doi/epdf/10.1137/0204007
 */
igraph_error_t igraph_simple_cycles_search_all(
    const igraph_t *graph,
    igraph_vector_int_list_t *result)
{
  igraph_simple_cycle_search_state_t state;
  igraph_simple_cycle_search_state_init(&state, graph);

  // int nfound = 0;
  // igraph_vector_int_list_init(result, 0); // state.N);

  // TODO: depending on the graph, it is rather unreasonable to search cycles from each and every node
  igraph_integer_t s = 0;
  while (s < state.N)
  {
    // if (!igraph_adjlist_empty(&state.AK))
    // {
    if (!igraph_vector_int_empty(igraph_adjlist_get(&state.AK, s)))
    {
      igraph_simple_cycles_search_one(&state, s, result);
    }
    s += 1;
    // }
    // else
    // {
    //   s = state.N;
    // }
  }

  // igraph_vector_int_list_resize(result, nfound);
  igraph_simple_cycle_search_state_destroy(&state);

  return IGRAPH_SUCCESS;
}
