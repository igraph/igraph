#include <cassert>
#include <vector>
#include "utils.hh"

/*
  Copyright (c) 2003-2015 Tommi Junttila
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

#if 0
void
print_permutation(FILE* const fp,
		  const unsigned int N,
		  const unsigned int* perm,
		  const unsigned int offset)
{
  assert(N > 0);
  assert(perm);
  for(unsigned int i = 0; i < N; i++) {
    unsigned int j = perm[i];
    if(j == i)
      continue;
    bool is_first = true;
    while(j != i) {
      if(j < i) {
        is_first = false;
        break;
      }
      j = perm[j];
    }
    if(!is_first)
      continue;
    fprintf(fp, "(%u,", i+offset);
    j = perm[i];
    while(j != i) {
      fprintf(fp, "%u", j+offset);
      j = perm[j];
      if(j != i)
        fprintf(fp, ",");
    }
    fprintf(fp, ")");
  }
}

void
print_permutation(FILE* const fp,
		  const std::vector<unsigned int>& perm,
		  const unsigned int offset)
{
  const unsigned int N = perm.size();
  for(unsigned int i = 0; i < N; i++) {
    unsigned int j = perm[i];
    if(j == i)
      continue;
    bool is_first = true;
    while(j != i) {
      if(j < i) {
        is_first = false;
        break;
      }
      j = perm[j];
    }
    if(!is_first)
      continue;
    fprintf(fp, "(%u,", i+offset);
    j = perm[i];
    while(j != i) {
      fprintf(fp, "%u", j+offset);
      j = perm[j];
      if(j != i)
        fprintf(fp, ",");
    }
    fprintf(fp, ")");
  }
}
#endif

bool
is_permutation(const unsigned int N, const unsigned int* perm)
{
  if(N == 0)
    return true;
  std::vector<bool> m(N, false);
  for(unsigned int i = 0; i < N; i++) {
    if(perm[i] >= N) return false;
    if(m[perm[i]]) return false;
    m[perm[i]] = true;
  }
  return true;
}

bool
is_permutation(const std::vector<unsigned int>& perm)
{
  const unsigned int N = perm.size();
  if(N == 0)
    return true;
  std::vector<bool> m(N, false);
  for(unsigned int i = 0; i < N; i++) {
    if(perm[i] >= N) return false;
    if(m[perm[i]]) return false;
    m[perm[i]] = true;
  }
  return true;
}


} // namespace bliss
