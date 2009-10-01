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

#ifndef BLISS_EQREFHASH_HH
#define BLISS_EQREFHASH_HH

#include <vector>

#define EqrefHash BuzzHash
//#define EqrefHash PerfectHash

namespace igraph {

class BuzzHash
{
protected:
  unsigned int h;
public:
  void reset() {h = 0; }
  void update(unsigned int);
  int cmp(const BuzzHash &other);
  bool is_lt(const BuzzHash &other) {return(cmp(other) < 0); }
  bool is_le(const BuzzHash &other) {return(cmp(other) <= 0); }
  bool is_equal(const BuzzHash &other) {return(cmp(other) == 0); }
};

class PerfectHash
{
protected:
  std::vector<unsigned int> h;
public:
  void reset() {h.clear(); }
  void update(unsigned int i) {h.push_back(i); }
  int cmp(const PerfectHash &other);
  bool is_lt(const PerfectHash &other) {return(cmp(other) < 0); }
  bool is_le(const PerfectHash &other) {return(cmp(other) <= 0); }
  bool is_equal(const PerfectHash &other) {return(cmp(other) == 0); }
};

}

#endif
