#include <cassert>
#include "orbit.hh"

/*
  Copyright (c) 2003-2021 Tommi Junttila
  Released under the GNU Lesser General Public License version 3.

  This file is part of bliss.

  bliss is free software: you can redistribute it and/or modify
  it under the terms of the GNU Lesser General Public License as published by
  the Free Software Foundation, version 3 of the License.

  bliss is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public License
  along with bliss.  If not, see <http://www.gnu.org/licenses/>.
*/

namespace bliss {

Orbit::Orbit()
{
  orbits = 0;
  in_orbit = 0;
  nof_elements = 0;
}


Orbit::~Orbit()
{
  delete[] orbits;
  orbits = 0;
  /*
  if(orbits)
    {
      free(orbits);
      orbits = 0;
    }
  */
  delete[] in_orbit;
  in_orbit = 0;
  /*
  if(in_orbit)
    {
      free(in_orbit);
      in_orbit = 0;
    }
  */
  nof_elements = 0;
  _nof_orbits = 0;
}


void Orbit::init(const unsigned int n)
{
  assert(n > 0);
  if(orbits) delete[] orbits;
  orbits = new OrbitEntry[n];
  delete[] in_orbit;
  in_orbit = new OrbitEntry*[n];
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

  merge_orbits(in_orbit[e1], in_orbit[e2]);
}


bool Orbit::is_minimal_representative(unsigned int element) const
{
  return(get_minimal_representative(element) == element);
}


unsigned int Orbit::get_minimal_representative(unsigned int element) const
{

  OrbitEntry * const orbit = in_orbit[element];

  return(orbit->element);
}


unsigned int Orbit::orbit_size(unsigned int element) const
{

  return(in_orbit[element]->size);
}


} // namespace bliss
