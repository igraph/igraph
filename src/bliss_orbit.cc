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

#include <cstdlib>
#include <cassert>
#include "bliss_defs.hh"
#include "bliss_orbit.hh"

using namespace std;

namespace igraph {

Orbit::Orbit()
{
  orbits = 0;
  in_orbit = 0;
  nof_elements = 0;
}


Orbit::~Orbit()
{
  if(orbits)
    {
      free(orbits);
      orbits = 0;
    }
  if(in_orbit)
    {
      free(in_orbit);
      in_orbit = 0;
    }
  nof_elements = 0;
}


void Orbit::init(const unsigned int n)
{
  assert(n > 0);
  if(orbits) free(orbits);
  orbits = (OrbitEntry*)malloc(n * sizeof(OrbitEntry));
  if(in_orbit) free(in_orbit);
  in_orbit = (OrbitEntry**)malloc(n * sizeof(OrbitEntry*));
  nof_elements = n;

  reset();
}


void Orbit::reset()
{
  assert(orbits);
  assert(in_orbit);

  for(unsigned int i = 0; i < nof_elements; i++)
    {
      orbits[i].element = i;
      orbits[i].next = 0;
      orbits[i].size = 1;
      in_orbit[i] = &orbits[i];
    }
  _nof_orbits = nof_elements;
}


void Orbit::merge_orbits(OrbitEntry *orbit1, OrbitEntry *orbit2)
{
  DEBUG_ASSERT((orbit1 == orbit2) == (orbit1->element == orbit2->element));
  DEBUG_ASSERT(orbit1->element < nof_elements);
  DEBUG_ASSERT(orbit2->element < nof_elements);

  if(orbit1 != orbit2)
    {
      _nof_orbits--;
      /* Only update the elements in the smaller orbit */
      if(orbit1->size > orbit2->size)
	{
	  OrbitEntry * const temp = orbit2;
	  orbit2 = orbit1;
	  orbit1 = temp;
	}
      /* Link the elements of orbit1 to the almost beginning of orbit2 */
      OrbitEntry *e = orbit1;
      while(e->next)
	{
	  in_orbit[e->element] = orbit2;
	  e = e->next;
	}
      in_orbit[e->element] = orbit2;
      e->next = orbit2->next;
      orbit2->next = orbit1;
      /* Keep the minimal orbit representative in the beginning */
      if(orbit1->element < orbit2->element)
	{
	  const unsigned int temp = orbit1->element;
	  orbit1->element = orbit2->element;
	  orbit2->element = temp;
	}
      orbit2->size += orbit1->size;
    }
}


void Orbit::merge_orbits(unsigned int e1, unsigned int e2)
{
  DEBUG_ASSERT(e1 < nof_elements);
  DEBUG_ASSERT(e2 < nof_elements);

  merge_orbits(in_orbit[e1], in_orbit[e2]);
}


bool Orbit::is_minimal_representative(unsigned int element)
{
  return(get_minimal_representative(element) == element);
}


unsigned int Orbit::get_minimal_representative(unsigned int element)
{
  DEBUG_ASSERT(element < nof_elements);

  OrbitEntry * const orbit = in_orbit[element];

  DEBUG_ASSERT(orbit->element <= element);
  return(orbit->element);
}


unsigned int Orbit::orbit_size(unsigned int element)
{
  DEBUG_ASSERT(element < nof_elements);
  
  return(in_orbit[element]->size);
}

}
