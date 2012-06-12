/*
Copyright (C) 2003-2006 Tommi Junttila

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License version 2
 as published by the Free Software Foundation.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA.
*/

/* FSF address fixed in the above notice on 1 Oct 2009 by Tamas Nepusz */

#ifndef BLISS_GRAPH_HH
#define BLISS_GRAPH_HH

namespace igraph {
  class AbstractGraph;
}

#include <cstdio>
#include <vector>
#include "bliss_kstack.hh"
#include "bliss_kqueue.hh"
#include "bliss_heap.hh"
#include "bliss_orbit.hh"
#include "bliss_partition.hh"
#include "bliss_bignum.hh"
#include "bliss_eqrefhash.hh"

#include "igraph_datatype.h"

namespace igraph {

typedef struct t_Stats
{
  BigNum group_size;
  long unsigned int nof_nodes;
  long unsigned int nof_leaf_nodes;
  long unsigned int nof_bad_nodes;
  long unsigned int nof_canupdates;
  long unsigned int nof_generators;
  unsigned long int max_level;
} Stats;

class AbstractGraph
{
  friend class Partition;
protected:
  AbstractGraph();
  virtual ~AbstractGraph();
  
  Partition p;

  bool in_search;

  bool refine_compare_certificate;
  bool refine_equal_to_first;
  unsigned int refine_first_path_subcertificate_end;
  int refine_cmp_to_best;
  unsigned int refine_best_path_subcertificate_end;

  /* Max mem used by long prune in megabytes */
  static const unsigned int long_prune_options_max_mem = 20;
  static const unsigned int long_prune_options_max_stored_auts = 50;
  unsigned int long_prune_max_stored_autss;
  std::vector<std::vector<bool> *> long_prune_fixed;
  std::vector<std::vector<bool> *> long_prune_mcrs;
  std::vector<bool> long_prune_temp;
  unsigned int long_prune_begin;
  unsigned int long_prune_end;
  void long_prune_init();
  void long_prune_add_automorphism(const unsigned int *aut);
  std::vector<bool> &long_prune_get_fixed(const unsigned int index);
  std::vector<bool> &long_prune_get_mcrs(const unsigned int index);
  void long_prune_swap(const unsigned int i, const unsigned int j);

  /*
   * Data structures and routines for refining the partition p into equitable
   */
  Heap neighbour_heap;
  virtual bool split_neighbourhood_of_unit_cell(Cell *) = 0;
  virtual void split_neighbourhood_of_cell(Cell * const) = 0;
  void refine_to_equitable();
  void refine_to_equitable(Cell *cell);
  void refine_to_equitable(Cell *cell1, Cell *cell2);
  void do_refine_to_equitable();
  unsigned int eqref_max_certificate_index;
  bool eqref_worse_than_certificate;
  //void eqref_update_hash(unsigned int i);
  EqrefHash eqref_hash;

  /* For debugging purposes only */
  virtual bool is_equitable() {assert(0); return false; }

  void print_permutation(FILE *, const unsigned int *perm);

  unsigned int *first_path_labeling;
  unsigned int *first_path_labeling_inv;
  Orbit         first_path_orbits;
  unsigned int *first_path_automorphism;

  unsigned int *best_path_labeling;
  unsigned int *best_path_labeling_inv;
  Orbit         best_path_orbits;
  unsigned int *best_path_automorphism;

  void update_labeling(unsigned int * const lab);
  void update_labeling_and_its_inverse(unsigned int * const lab,
				       unsigned int * const lab_inv);
  void update_orbit_information(Orbit &o, const unsigned int *perm);

  void reset_permutation(unsigned int *perm);

  /* Mainly for debugging purposes */
  virtual bool is_automorphism(unsigned int * const perm);

  std::vector<unsigned int> certificate_current_path;
  std::vector<unsigned int> certificate_first_path;
  std::vector<unsigned int> certificate_best_path;

  //unsigned int *certificate;
  unsigned int certificate_size;
  unsigned int certificate_index;
  virtual void initialize_certificate() = 0;

  virtual void remove_duplicate_edges() = 0;
  virtual void make_initial_equitable_partition() = 0;
  virtual Cell *find_next_cell_to_be_splitted(Cell *cell) = 0;

  void search(const bool canonical, Stats &stats);

#ifdef PRINT_SEARCH_TREE_DOT
  FILE *dotty_output;
#endif

public:
  virtual unsigned int get_nof_vertices() = 0;

  void find_automorphisms(Stats &stats);
  const unsigned int *canonical_form(Stats &stats);
};




/*
 * Undirected, vertex labeled graph
 * Multiple edges between vertices are not allowed (i.e., will be ignored)
 */
class Graph : public AbstractGraph
{
  class Vertex {
  public:
    Vertex();
    ~Vertex();
    void add_edge(const unsigned int other_vertex);
    void remove_duplicate_edges(bool * const);

    unsigned int label;
    unsigned int nof_edges;
    std::vector<unsigned int> edges;
  };
  std::vector<Vertex> vertices;
  void remove_duplicate_edges();

  /*
   * Partition independent invariants for this graph class
   */
  static unsigned int label_invariant(Graph *g, unsigned int v);
  static unsigned int degree_invariant(Graph *g, unsigned int v);

  bool refine_according_to_invariant(unsigned int (*inv)(Graph * g,
							 unsigned int v));

  /*
   * Routines needed when refining the partition p into equitable
   */
  bool split_neighbourhood_of_unit_cell(Cell *);
  void split_neighbourhood_of_cell(Cell * const);

  /* For debugging purposes only */
  bool is_equitable();

  Cell *find_next_cell_to_be_splitted(Cell *cell);
  /* Splitting heuristics */
  unsigned int sh;
  Cell *sh_first(Cell *cell);
  Cell *sh_first_smallest(Cell *cell);
  Cell *sh_first_largest(Cell *cell);
  Cell *sh_first_max_neighbours(Cell *cell);
  Cell *sh_first_smallest_max_neighbours(Cell *cell);
  Cell *sh_first_largest_max_neighbours(Cell *cell);

  void make_initial_equitable_partition();
  void initialize_certificate();

  bool is_automorphism(unsigned int * const perm);

public:
  Graph(const unsigned int nof_vertices = 0);
  ~Graph();

  static const unsigned int sh_f   = 0;
  static const unsigned int sh_fs  = 1;
  static const unsigned int sh_fl  = 2;
  static const unsigned int sh_fm  = 3;
  static const unsigned int sh_fsm = 4;
  static const unsigned int sh_flm = 5;

  static Graph *read_dimacs(FILE *);
  static Graph *from_igraph(const igraph_t *graph);
  void print_dimacs(FILE *);

  void to_dot(FILE *fp);
  void to_dot(const char *file_name);

  unsigned int get_nof_vertices() {return vertices.size(); }
  Graph *permute(const unsigned int *perm);
  
  unsigned int add_vertex(const unsigned int label = 1);
  void add_edge(const unsigned int vertex1, const unsigned int vertex2);
  void change_label(const unsigned int vertex, const unsigned int new_label);

  void set_splitting_heuristics(unsigned int shs) {sh = shs; }
};

}

#endif
