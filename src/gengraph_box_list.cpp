#include "gengraph_box_list.h"
#include <assert.h>

namespace gengraph {

void box_list::insert(int v) {
  register int d = deg[v];
  if(d<1) return;
  if(d>dmax) dmax=d;
  int yo = list[d-1];
  list[d-1] = v;
  prev[v] = -1;
  next[v] = yo;
  if(yo>=0) prev[yo]=v;
}

void box_list::pop(int v) {
  register int p = prev[v];
  register int n = next[v];
  if(p<0) {
    register int d = deg[v];
    assert(list[d-1]==v);
    list[d-1] = n;
    if(d==dmax && n<0) do dmax--; while(dmax>0 && list[dmax-1]<0);
  }
  else next[p] = n;
  if(n>=0) prev[n] = p;
}

box_list::box_list(int n0, int *deg0) : n(n0), deg(deg0) {
  next = new int[n];
  prev = new int[n];
  dmax = -1;
  int i;
  for(i=0; i<n; i++) if(deg[i]>dmax) dmax=deg[i];
  list = new int[dmax];
  for(i=0; i<dmax; i++) list[i]=-1;
  for(i=0; i<n; i++) insert(i);
}

box_list::~box_list() {
  delete[] prev;
  delete[] next;
  delete[] list;
}

void box_list::pop_vertex(int v, int **neigh) {
  int k=deg[v];
  if(k<1) return;
  pop(v);
  int *w = neigh[v];
  while(k--) {
    int v2 = *(w++);
    register int *w2 = neigh[v2];
    while(*w2 != v) w2++;
    register int *w3 = neigh[v2]+(deg[v2]-1);
    assert(w2<=w3);
    register int tmp = *w3;
    *w3 = *w2;
    *w2 = tmp;
    pop(v2);
    deg[v2]--;
    insert(v2);
  }
}

} // namespace gengraph
