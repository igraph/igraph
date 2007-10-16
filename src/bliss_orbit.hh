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
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#ifndef ORBIT_HH
#define ORBIT_HH

namespace igraph {

class Orbit
{
protected:
  class OrbitEntry {
  public:
    unsigned int element;
    OrbitEntry *next;
    unsigned int size;
  };
  OrbitEntry *orbits;
  OrbitEntry **in_orbit;
  unsigned int nof_elements;
  unsigned int _nof_orbits;
  void merge_orbits(OrbitEntry *o1, OrbitEntry *o2);
public:
  Orbit();
  ~Orbit();
  void init(const unsigned int n);
  void reset();
  void merge_orbits(unsigned int e1, unsigned int e2);
  bool is_minimal_representative(unsigned int e);
  unsigned int get_minimal_representative(unsigned int e);
  unsigned int orbit_size(unsigned int e);
  unsigned int nof_orbits() const {return _nof_orbits; }
};

}

#endif
