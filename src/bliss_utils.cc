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

#include "bliss_utils.hh"
#include "bliss_bignum.hh"

using namespace std;

#include <cstring>

#if defined(BLISS_USE_GMP)

namespace igraph {

int BigNum::tostring(char **str) { 
  *str=igraph_Calloc(mpz_sizeinbase(v, 10)+2, char);
  if (! *str) { 
    IGRAPH_ERROR("Cannot convert big number to string", IGRAPH_ENOMEM);
  }
  mpz_get_str(*str, 10, v);
  return 0;
}

}

#else

namespace igraph {

int BigNum::tostring(char **str) {
  int size=static_cast<int>( (log(abs(v))/log(10.0))+4 );
  *str=igraph_Calloc(size, char );
  if (! *str) {
    IGRAPH_ERROR("Cannot convert big number to string", IGRAPH_ENOMEM);
  }
  std::stringstream ss;
  ss << v;
  strncpy(*str, ss.str().c_str(), size);
  return 0;
}

}

#endif
