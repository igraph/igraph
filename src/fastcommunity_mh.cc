////////////////////////////////////////////////////////////////////////
// --- COPYRIGHT NOTICE ---------------------------------------------
// FastCommunityMH - infers community structure of networks
// Copyright (C) 2004 Aaron Clauset
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
// 
// See http://www.gnu.org/licenses/gpl.txt for more details.
// 
////////////////////////////////////////////////////////////////////////
// Author       : Aaron Clauset  (aaron@cs.unm.edu)				//
// Location     : U. Michigan, U. New Mexico						//
// Time         : January-August 2004							//
// Collaborators: Dr. Cris Moore (moore@cs.unm.edu)				//
//              : Dr. Mark Newman (mejn@umich.edu)				//
////////////////////////////////////////////////////////////////////////
// --- DEEPER DESCRIPTION ---------------------------------------------
//  see http://www.arxiv.org/abs/cond-mat/0408187 for more information
// 
//  - read network structure from data file (see below for constraints)
//  - builds dQ, H and a data structures
//  - runs new fast community structure inference algorithm
//  - records Q(t) function to file
//  - (optional) records community structure (at t==cutstep)
//  - (optional) records the list of members in each community (at t==cutstep)
//
////////////////////////////////////////////////////////////////////////
// --- PROGRAM USAGE NOTES ---------------------------------------------
// This program is rather complicated and requires a specific kind of input,
// so some notes on how to use it are in order. Mainly, the program requires
// a specific structure input file (.pairs) that has the following characteristics:
//  
//  1. .pairs is a list of tab-delimited pairs of numeric indices, e.g.,
//		"54\t91\n"
//  2. the network described is a SINGLE COMPONENT
//  3. there are NO SELF-LOOPS or MULTI-EDGES in the file; you can use
//     the 'netstats' utility to extract the giantcomponent (-gcomp.pairs)
//     and then use that file as input to this program
//  4. the MINIMUM NODE ID = 0 in the input file; the maximum can be
//     anything (the program will infer it from the input file)
// 
// Description of commandline arguments
// -f <filename>    give the target .pairs file to be processed
// -l <text>		the text label for this run; used to build output filenames
// -t <int>		timer period for reporting progress of file input to screen
// -s			calculate and record the support of the dQ matrix
// -v --v ---v		differing levels of screen output verbosity
// -c <int>		record the aglomerated network at step <int>
// 
////////////////////////////////////////////////////////////////////////
// Change Log:
// 2006-02-06: 1) modified readInputFile to be more descriptive of its actions
//             2) removed -o functionality; program writes output to directory
//             of input file. (also removed -h option of command line)
// 2006-07-31: 3) changed -c behavior so input should be exactly same number
//             as given in the .info file (previously it differed by 1)
// 2006-10-13: 4) Janne Aukia (jaukia@cc.hut.fi) suggested changes to the 
//             mergeCommunities() function here (see comments in that function),
//             and an indexing adjustment in printHeapTop10() in maxheap.h.
//
////////////////////////////////////////////////////////////////////////

#include <string>
#include "stdlib.h"
#include "time.h"
#include "math.h"

#include "fastcommunity_maxheap.h"
#include "fastcommunity_vektor.h"

#include "igraph.h"

using namespace std;

// ------------------------------------------------------------------------------------
// Edge object - defined by a pair of vertex indices and *edge pointer to next in linked-list
class edge {
public:
  int     so;					// originating node
  int     si;					// terminating node
  edge    *next;					// pointer for linked list of edges
	
  edge();						// default constructor
  ~edge();						// default destructor
};
edge::edge()  { so = 0; si = 0; next = NULL; }
edge::~edge() {}

// ------------------------------------------------------------------------------------
// Nodenub object - defined by a *node pointer and *node pointer 
struct nodenub {
  tuple	*heap_ptr;			// pointer to node(max,i,j) in max-heap of row maxes
  vektor    *v;					// pointer stored vector for (i,j)
};

// ------------------------------------------------------------------------------------
// tuple object - defined by an real value and (row,col) indices
#if !defined(TUPLE_INCLUDED)
#define TUPLE_INCLUDED
struct tuple {
  double    m;					// stored value
  int		i;					// row index
  int		j;					// column index
  int		k;					// heap index
};
#endif

// ordered pair structures (handy in the program)
struct apair { int x; int y; };
#if !defined(DPAIR_INCLUDED)
#define DPAIR_INCLUDED
class dpair {
public:
  int x; double y; dpair *next;
  dpair(); ~dpair();
};
dpair::dpair()  { x = 0; y = 0.0; next = NULL; }
dpair::~dpair() {}
#endif

// ------------------------------------------------------------------------------------
// List object - simple linked list of integers
class list {
public:
  int		index;				// node index
  list		*next;				// pointer to next element in linked list
  list();   ~list();
};
list::list()  { index= 0; next = NULL; }
list::~list() {}

// ------------------------------------------------------------------------------------
// Community stub object - stub for a community list
class stub {
public:
  bool		valid;				// is this community valid?
  int		size;				// size of community
  list		*members;				// pointer to list of community members
  list		*last;				// pointer to end of list
  stub();   ~stub();
};
stub::stub()  { valid = false; size = 0; members = NULL; last = NULL; }
stub::~stub() {
  list *current;
  if (members != NULL) {
    current = members;
    while (current != NULL) { members = current->next; delete current; current = members; }
  }
}

// ------------------------------------------------------------------------------------
// PROGRAM PARAMETERS -----------------------------------------------------------------

typedef struct netparameters {
  int			n;				// number of nodes in network
  int			m;				// number of edges in network
  int			maxid;			// maximum node id
  int			minid;			// minimum node id
} igraph_i_netparameters_t;

typedef struct groupstats {
  int			numgroups;		// number of groups
  double		meansize;			// mean size of groups
  int			maxsize;			// size of largest group
  int			minsize;			// size of smallest group
  double		*sizehist;		// distribution of sizes
} igraph_i_groupstats_t;

// ------------------------------------------------------------------------------------
// FUNCTION DECLARATIONS --------------------------------------------------------------

void igraph_i_buildDeltaQMatrix(igraph_i_netparameters_t &gparm, igraph_i_groupstats_t &gstats,
		       double *&a, edge *&e, double *&Q, nodenub *&dq, maxheap *&h,
		       edge *&elist);
void igraph_i_mergeCommunities(igraph_i_netparameters_t &gparm, igraph_i_groupstats_t &gstats,
		      nodenub *&dq, maxheap *&h, double *&a,
		      int i, int j);
void igraph_i_readInputFile(igraph_i_netparameters_t &gparm, igraph_i_groupstats_t &gstats,
		   edge *&elist, edge *&e, 
		   const igraph_t *graph);

// ------------------------------------------------------------------------------------
// ----------------------------------- GLOBAL VARIABLES -------------------------------

enum {NONE};

int igraph_community_fastgreedy(const igraph_t *graph,
				igraph_matrix_t *merges,
				igraph_vector_t *modularity) {

  long int no_of_nodes=(long int)igraph_vcount(graph);
  igraph_i_netparameters_t gparm;
  igraph_i_groupstats_t gstats;

  maxheap   *h=0;				// heap of values from max_i{dQ_ij}
  nodenub *dq=0;
  edge	*e=0;				// initial adjacency matrix (sparse)
  edge *elist=0;

  igraph_i_readInputFile(gparm, gstats, elist, e, graph);
  
  if (merges) {
    IGRAPH_CHECK(igraph_matrix_resize(merges, no_of_nodes-1, 2));
  }
  if (modularity) {
    IGRAPH_CHECK(igraph_vector_resize(modularity, no_of_nodes));
  }
  
  // ----------------------------------------------------------------------
  // Allocate data structures for main loop
  double *a     = new double [gparm.maxid];
  double *Q     = new double [gparm.n+1];
  apair *joins = new apair  [gparm.n+1];
  for (int i=0; i<gparm.maxid; i++) { a[i] = 0.0; }
  for (int i=0; i<gparm.n+1;   i++) { Q[i] = 0.0; joins[i].x = 0; joins[i].y = 0; }
  int t = 1;
  
  // builds dQ[] and h
  igraph_i_buildDeltaQMatrix(gparm, gstats, a, e, Q, dq, h, elist);

  if (modularity) {
    VECTOR(*modularity)[0]=Q[0];
  }

  // ----------------------------------------------------------------------
  // Start FastCommunity algorithm
  tuple  dQmax;
  int isupport, jsupport;
  while (h->heapSize() >= 2) {
		
    // ---------------------------------
    // Find largest dQ
    dQmax = h->popMaximum();					// select maximum dQ_ij // convention: insert i into j
    if (dQmax.m < -4000000000.0) { break; }		// no more joins possible
		
    // ---------------------------------
    // Merge the chosen communities
    isupport = dq[dQmax.i].v->returnNodecount();
    jsupport = dq[dQmax.j].v->returnNodecount();
    if (isupport < jsupport) {
      igraph_i_mergeCommunities(gparm, gstats, dq, h, a, dQmax.i, dQmax.j);	// merge community i into community j
      joins[t].x = dQmax.i;				// record merge of i(x) into j(y)
      joins[t].y = dQmax.j;				// 
    } else {								// 
      dq[dQmax.i].heap_ptr = dq[dQmax.j].heap_ptr; // take community j's heap pointer
      dq[dQmax.i].heap_ptr->i = dQmax.i;			//   mark it as i's
      dq[dQmax.i].heap_ptr->j = dQmax.j;			//   mark it as i's
      igraph_i_mergeCommunities(gparm, gstats, dq, h, a, dQmax.j, dQmax.i);	// merge community j into community i
      joins[t].x = dQmax.j;				// record merge of j(x) into i(y)
      joins[t].y = dQmax.i;				// 
    }									// 
    Q[t] = dQmax.m + Q[t-1];					// record Q(t)
		
    // ---------------------------------
    // Record join to file
    if (merges) {
      MATRIX(*merges, t-1, 0)=joins[t].x-1;
      MATRIX(*merges, t-1, 1)=joins[t].y-1;
    }
    if (modularity) {
      VECTOR(*modularity)[t]=Q[t];
    }

    t++;									// increment time
  } // ------------- end community merging loop
	
  return 0;  
}

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //
// FUNCTION DEFINITIONS --------------------------------------------------------------- //

void igraph_i_buildDeltaQMatrix(igraph_i_netparameters_t &gparm, igraph_i_groupstats_t &gstats,
		       double *&a, edge *&e, double *&Q, nodenub *&dq, maxheap *&h,
		       edge *&elist) {

  // Given that we've now populated a sparse (unordered) adjacency matrix e (e), 
  // we now need to construct the intial dQ matrix according to the definition of dQ
  // which may be derived from the definition of modularity Q:
  //    Q(t) = \sum_{i} (e_{ii} - a_{i}^2) = Tr(e) - ||e^2||
  // thus dQ is
  //    dQ_{i,j} = 2* ( e_{i,j} - a_{i}a_{j} )
  //    where a_{i} = \sum_{j} e_{i,j} (i.e., the sum over the ith row)
  // To create dQ, we must insert each value of dQ_{i,j} into a binary search tree,
  // for the jth column. That is, dQ is simply an array of such binary search trees,
  // each of which represents the dQ_{x,j} adjacency vector. Having created dQ as
  // such, we may happily delete the matrix e in order to free up more memory.
  // The next step is to create a max-heap data structure, which contains the entries
  // of the following form (value, s, t), where the heap-key is 'value'. Accessing the
  // root of the heap gives us the next dQ value, and the indices (s,t) of the vectors
  // in dQ which need to be updated as a result of the merge.
	
	
  // First we compute e_{i,j}, and the compute+store the a_{i} values. These will be used
  // shortly when we compute each dQ_{i,j}.
  edge   *current;
  double  eij = (double)(0.5/gparm.m);				// intially each e_{i,j} = 1/m
  for (int i=1; i<gparm.maxid; i++) {				// for each row
    a[i] = 0.0;								// 
    if (e[i].so != 0) {							//    ensure it exists
      current = &e[i];						//    grab first edge
      a[i] = eij;							//    initialize a[i]
      while (current->next != NULL) {			//    loop through remaining edges
	a[i] += eij;						//       add another eij
	current = current->next;				//
      }
      Q[0] += -1.0*a[i]*a[i];					// calculate initial value of Q
    }
  }

  // now we create an empty (ordered) sparse matrix dq[]
  dq = new nodenub [gparm.maxid];						// initialize dq matrix
  for (int i=0; i<gparm.maxid; i++) {					// 
    dq[i].heap_ptr = NULL;							// no pointer in the heap at first
    if (e[i].so != 0) { dq[i].v = new vektor(2+(int)floor(gparm.m*a[i])); }
    else {			dq[i].v = NULL; }
  }
  h = new maxheap(gparm.n);						// allocate max-heap of size = number of nodes
	
  // Now we do all the work, which happens as we compute and insert each dQ_{i,j} into 
  // the corresponding (ordered) sparse vector dq[i]. While computing each dQ for a
  // row i, we track the maximum dQmax of the row and its (row,col) indices (i,j). Upon
  // finishing all dQ's for a row, we insert the tuple into the max-heap hQmax. That
  // insertion returns the itemaddress, which we then store in the nodenub heap_ptr for 
  // that row's vector.
  double    dQ;
  tuple	dQmax;										// for heaping the row maxes

  for (int i=1; i<gparm.maxid; i++) {
    if (e[i].so != 0) {
      current = &e[i];								// grab first edge
      dQ      = 2.0*(eij-(a[current->so]*a[current->si]));   // compute its dQ
      dQmax.m = dQ;									// assume it is maximum so far
      dQmax.i = current->so;							// store its (row,col)
      dQmax.j = current->si;							// 
      dq[i].v->insertItem(current->si, dQ);				// insert its dQ
      while (current->next != NULL) {					// 
	current = current->next;						// step to next edge
	dQ = 2.0*(eij-(a[current->so]*a[current->si]));	// compute new dQ
	if (dQ > dQmax.m) {							// if dQ larger than current max
	  dQmax.m = dQ;							//    replace it as maximum so far
	  dQmax.j = current->si;					//    and store its (col)
	}
	dq[i].v->insertItem(current->si, dQ);			// insert it into vector[i]
      }
      dq[i].heap_ptr = h->insertItem(dQmax);				// store the pointer to its loc in heap
    }
  }

  delete [] elist;								// free-up adjacency matrix memory in two shots
  delete [] e;									// 
  return;
}

// ------------------------------------------------------------------------------------ //

void igraph_i_mergeCommunities(igraph_i_netparameters_t &gparm, igraph_i_groupstats_t &gstats,
		      nodenub *&dq, maxheap *&h, double *&a,
		      int i, int j) {
	
  // To do the join operation for a pair of communities (i,j), we must update the dQ
  // values which correspond to any neighbor of either i or j to reflect the change.
  // In doing this, there are three update rules (cases) to follow:
  //  1. jix-triangle, in which the community x is a neighbor of both i and j
  //  2. jix-chain, in which community x is a neighbor of i but not j
  //  3. ijx-chain, in which community x is a neighbor of j but not i
  //
  // For the first two cases, we may make these updates by simply getting a list of
  // the elements (x,dQ) of [i] and stepping through them. If x==j, then we can ignore
  // that value since it corresponds to an edge which is being absorbed by the joined
  // community (i,j). If [j] also contains an element (x,dQ), then we have a triangle;
  // if it does not, then we have a jix-chain.
  //
  // The last case requires that we step through the elements (x,dQ) of [j] and update each
  // if [i] does not have an element (x,dQ), since that implies a ijx-chain.
  // 
  // Let d([i]) be the degree of the vector [i], and let k = d([i]) + d([j]). The running
  // time of this operation is O(k log k)
  //
  // Essentially, we do most of the following operations for each element of
  // dq[i]_x where x \not= j
  //  1.  add dq[i]_x to dq[j]_x (2 cases)
  //  2.  remove dq[x]_i
  //  3.  update maxheap[x]
  //  4.  add dq[i]_x to dq[x]_j (2 cases)
  //  5.  remove dq[j]_i
  //  6.  update maxheap[j]
  //  7.  update a[j] and a[i]
  //  8.  delete dq[i]
	
  dpair *list, *current, *temp;
  tuple newMax;
  int t = 1;
	
  // -- Working with the community being inserted (dq[i])
  // The first thing we must do is get a list of the elements (x,dQ) in dq[i]. With this 
  // list, we can then insert each into dq[j].

  //	dq[i].v->printTree();
  list    = dq[i].v->returnTreeAsList();			// get a list of items in dq[i].v
  current = list;							// store ptr to head of list
	
  // ---------------------------------------------------------------------------------
  // SEARCHING FOR JIX-TRIANGLES AND JIX-CHAINS --------------------------------------
  // Now that we have a list of the elements of [i], we can step through them to check if
  // they correspond to an jix-triangle, a jix-chain, or the edge (i,j), and do the appropriate
  // operation depending.
	
  while (current!=NULL) {						// insert list elements appropriately
		
    // If the element (x,dQ) is actually (j,dQ), then we can ignore it, since it will 
    // correspond to an edge internal to the joined community (i,j) after the join.
    if (current->x != j) {

      // Now we must decide if we have a jix-triangle or a jix-chain by asking if
      // [j] contains some element (x,dQ). If the following conditional is TRUE,
      // then we have a jix-triangle, ELSE it is a jix-chain.
			
      if (dq[j].v->findItem(current->x)) {
	// CASE OF JIX-TRIANGLE
				
	// We first add (x,dQ) from [i] to [x] as (j,dQ), since [x] essentially now has
	// two connections to the joined community [j].
	dq[current->x].v->insertItem(j,current->y);			// (step 1)
	
	// Then we need to delete the element (i,dQ) in [x], since [i] is now a
	// part of [j] and [x] must reflect this connectivity.
	dq[current->x].v->deleteItem(i);					// (step 2)
	
	// After deleting an item, the tree may now have a new maximum element in [x],
	// so we need to check it against the old maximum element. If it's new, then
	// we need to update that value in the heap and reheapify.
	newMax = dq[current->x].v->returnMaxStored();		// (step 3)
	//				if (newMax.m > dq[current->x].heap_ptr->m || dq[current->x].heap_ptr->j==i) {
	h->updateItem(dq[current->x].heap_ptr, newMax);
	//				}
	// Change suggested by Janne Aukia (jaukia@cc.hut.fi) on 12 Oct 2006
				
	// Finally, we must insert (x,dQ) into [j] to note that [j] essentially now
	// has two connections with its neighbor [x].
	dq[j].v->insertItem(current->x,current->y);		// (step 4)
	
      } else {
	// CASE OF JIX-CHAIN
				
	// The first thing we need to do is calculate the adjustment factor (+) for updating elements.
	double axaj = -2.0*a[current->x]*a[j];
				
	// Then we insert a new element (j,dQ+) of [x] to represent that [x] has
	// acquired a connection to [j], which was [x]'d old connection to [i]
	dq[current->x].v->insertItem(j,current->y + axaj);	// (step 1)
				
	// Now the deletion of the connection from [x] to [i], since [i] is now
	// a part of [j]
	dq[current->x].v->deleteItem(i);					// (step 2)
				
	// Deleting that element may have changed the maximum element for [x], so we
	// need to check if the maximum of [x] is new (checking it against the value
	// in the heap) and then update the maximum in the heap if necessary.
	newMax = dq[current->x].v->returnMaxStored();		// (step 3)
	//				if (newMax.m > dq[current->x].heap_ptr->m || dq[current->x].heap_ptr->j==i) {
	h->updateItem(dq[current->x].heap_ptr, newMax);
	//				}
	// Change suggested by Janne Aukia (jaukia@cc.hut.fi) on 12 Oct 2006
					
	// Finally, we insert a new element (x,dQ+) of [j] to represent [j]'s new
	// connection to [x]
	dq[j].v->insertItem(current->x,current->y + axaj);	// (step 4)

      }    // if (dq[j].v->findItem(current->x))
			
    }    // if (current->x != j)
		
    temp    = current;
    current = current->next;						// move to next element
    delete temp;
    temp = NULL;
    t++;
  }    // while (current!=NULL)

  // We've now finished going through all of [i]'s connections, so we need to delete the element
  // of [j] which represented the connection to [i]

  dq[j].v->deleteItem(i);						// (step 5)
	
  // We can be fairly certain that the maximum element of [j] was also the maximum
  // element of [i], so we need to check to update the maximum value of [j] that
  // is in the heap.
  newMax = dq[j].v->returnMaxStored();			// (step 6)
  h->updateItem(dq[j].heap_ptr, newMax);
	
  // ---------------------------------------------------------------------------------
  // SEARCHING FOR IJX-CHAINS --------------------------------------------------------
  // So far, we've treated all of [i]'s previous connections, and updated the elements
  // of dQ[] which corresponded to neighbors of [i] (which may also have been neighbors
  // of [j]. Now we need to update the neighbors of [j] (as necessary)
	
  // Again, the first thing we do is get a list of the elements of [j], so that we may
  // step through them and determine if that element constitutes an ijx-chain which
  // would require some action on our part.
  list = dq[j].v->returnTreeAsList();			// get a list of items in dq[j].v
  current = list;							// store ptr to head of list
  t       = 1;

  while (current != NULL) {					// insert list elements appropriately

    // If the element (x,dQ) of [j] is not also (i,dQ) (which it shouldn't be since we've
    // already deleted it previously in this function), and [i] does not also have an
    // element (x,dQ), then we have an ijx-chain.
    if ((current->x != i) && (!dq[i].v->findItem(current->x))) {
      // CASE OF IJX-CHAIN
			
      // First we must calculate the adjustment factor (+).
      double axai = -2.0*a[current->x]*a[i];
			
      // Now we must add an element (j,+) to [x], since [x] has essentially now acquired
      // a new connection to [i] (via [j] absorbing [i]).
      dq[current->x].v->insertItem(j, axai);			// (step 1)

      // This new item may have changed the maximum dQ of [x], so we must update it.
      newMax = dq[current->x].v->returnMaxStored();	// (step 3)
      h->updateItem(dq[current->x].heap_ptr, newMax);

      // And we must add an element (x,+) to [j], since [j] as acquired a new connection
      // to [x] (via absorbing [i]).
      dq[j].v->insertItem(current->x, axai);			// (step 4)
      newMax = dq[j].v->returnMaxStored();			// (step 6)
      h->updateItem(dq[j].heap_ptr, newMax);

    }    //  (current->x != i && !dq[i].v->findItem(current->x))
		
    temp    = current;
    current = current->next;						// move to next element
    delete temp;
    temp = NULL;
    t++;
  }    // while (current!=NULL)
	
  // Now that we've updated the connections values for all of [i]'s and [j]'s neighbors, 
  // we need to update the a[] vector to reflect the change in fractions of edges after
  // the join operation.
  a[j] += a[i];								// (step 7)
  a[i] = 0.0;
	
  // ---------------------------------------------------------------------------------
  // Finally, now we need to clean up by deleting the vector [i] since we'll never
  // need it again, and it'll conserve memory. For safety, we also set the pointers
  // to be NULL to prevent inadvertent access to the deleted data later on.

  delete dq[i].v;							// (step 8)
  dq[i].v        = NULL;						// (step 8)
  dq[i].heap_ptr = NULL;						//
	
  return;
 
}

// ------------------------------------------------------------------------------------

void igraph_i_readInputFile(igraph_i_netparameters_t &gparm, igraph_i_groupstats_t &gstats,
		   edge *&elist, edge *&e, 
		   const igraph_t *graph) {
	
  // temporary variables for this function
  long int no_of_edges=(long int)igraph_ecount(graph);
  igraph_integer_t from ,to;
  long int i;
  long int numnodes = (long int)igraph_vcount(graph)-1;
  long int numlinks = no_of_edges;
  long int s,f,t;
  edge **last;
  edge *newedge;
  edge *current;								// pointer for checking edge existence
  bool existsFlag;							// flag for edge existence
	
  gparm.maxid = numnodes+2;					// store maximum index
  elist = new edge [2*numlinks];				// create requisite number of edges
  long int ecounter = 0;							// index of next edge of elist to be used

  // Now that we know numnodes, we can allocate the space for the sparse matrix, and
  // then reparse the file, adding edges as necessary.
  e        = new  edge [gparm.maxid];			// (unordered) sparse adjacency matrix
  last     = new edge* [gparm.maxid];			// list of pointers to the last edge in each row
  numnodes = 0;								// numnodes now counts number of actual used node ids
  numlinks = 0;								// numlinks now counts number of bi-directional edges created
	
  for (i=0; i<no_of_edges; i++) {
    igraph_edge(graph, i, &from, &to);
    s=(long int)from+1; f=(long int)to+1;
    if (f < s) { t = s; s = f; f = t; }		// guarantee s < f
    numlinks++;							// increment link count (preemptive)
    if (e[s].so == 0) {						// if first edge with s, add s and (s,f)
      e[s].so = s;						// 
      e[s].si = f;						// 
      last[s] = &e[s];					//    point last[s] at self
      numnodes++;						//    increment node count
    } else {								//    try to add (s,f) to s-edgelist
      current = &e[s];					// 
      existsFlag = false;					// 
      while (current != NULL) {			// check if (s,f) already in edgelist
	if (current->si==f) {			// 
	  existsFlag = true;			//    link already exists
	  numlinks--;				//    adjust link-count downward
	  break;					// 
	}							// 
	current = current->next;			//    look at next edge
      }								// 
      if (!existsFlag) {					// if not already exists, append it
	newedge = &elist[ecounter++];		//    grab next-free-edge
	newedge -> so = s;				// 
	newedge -> si = f;				// 
	last[s] -> next = newedge;		//    append newedge to [s]'s list
	last[s]         = newedge;		//    point last[s] to newedge
      }								// 
    }									// 
	  
    if (e[f].so == 0) {						// if first edge with f, add f and (f,s)
      e[f].so = f;						// 
      e[f].si = s;						// 
      last[f] = &e[f];					//    point last[s] at self
      numnodes++;						//    increment node count
    } else {								// try to add (f,s) to f-edgelist
      if (!existsFlag) {					//    if (s,f) wasn't in s-edgelist, then
	newedge = &elist[ecounter++];		//       (f,s) not in f-edgelist
	newedge -> so = f;				// 
	newedge -> si = s;				// 
	last[f] -> next = newedge;		//    append newedge to [f]'s list
	last[f]		 = newedge;		//    point last[f] to newedge
      }								// 
    }									
    existsFlag = false;						// reset existsFlag
	  
  }
	
  gparm.m = numlinks;							// store actual number of edges created
  gparm.n = numnodes;							// store actual number of nodes used
  return;
}
