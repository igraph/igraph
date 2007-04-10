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

#if !defined(TUPLE_INCLUDED)
#define TUPLE_INCLUDED
struct tuple {
	double    m;					// stored value
	int		i;					// row index
	int		j;					// column index
	int		k;					// heap index
};
#endif

#if !defined(MAXHEAP_INCLUDED)
#define MAXHEAP_INCLUDED

/*   Because using this heap requires us to be able to modify an arbitrary element's
	data in constant O(1) time, I use to tricky tactic of having elements in an array-
	based heap only contain addresses to the data, rather than the data itself. In this
	manner, an external program may freely modify an element of the heap, given that it
	possesses a pointer to said element (in an array-based heap, the addresses and the 
	value in that address are not bound and thus may change during the heapify() operation).
*/

struct hnode { tuple     *d; };
const int heapmin = 3;
//const double tiny = -4294967296.0;

class maxheap {
private:
	hnode    *A;					// maxheap array
	int       heaplimit;			// first unused element of heap array
	int		arraysize;			// size of array
	bool		isempty;				// T if heap is empty; F otherwise

	int		downsift(int i);		// sift A[i] down in heap
	int		upsift  (int i);		// sift A[i] up in heap
	int		left    (int i);		// returns index of left child
	int		right   (int i);		// returns index of right child
	int		parent  (int i);		// returns index of parent
	void		grow();				// increase size of array A
	void		shrink();				// decrease size of array A
  
public:
	maxheap();					// default constructor
	maxheap(int size);				// default constructor
	~maxheap();					// default destructor
	
	int		heapSize();							// returns heaplimit value
	bool		heapIsEmpty();							// returns isempty value
	tuple	*insertItem(const tuple newData);			// heap-inserts newData, returns the address of it
	tuple	popMaximum();							// removes and returns heap max, reheapifies
	tuple	returnMaximum();						// returns the heap max; no change to heap
	void		updateItem(tuple *address, tuple newData);   // updates the value of the tuple at address
	void		updateItem(tuple *address, double newStored);// update only the stored value of tuple at address
	void		deleteItem(tuple *address);				// remove an item from the heap
	int		returnArraysize();						// 
	int		returnHeaplimit();						// 
	
};

// ------------------------------------------------------------------------------------
// Max Heap methods

// Constructor + Destructor -----------------------------------------------------------
maxheap::maxheap() {
	tuple *newtuple;
	heaplimit = 1;							// first free location is A[1]
	arraysize = heapmin+1;					// 
	isempty   = true;						// heap is initially empty
	A = new hnode [arraysize];				// allocate array for heap
	for (int i=0; i<arraysize; i++) {			// initialize heap values
		newtuple = new tuple;				// 
		A[i].d = newtuple;					// 
		A[i].d->m = -4294967296.0;			// a very negative value; unlikely to be a valid dQ
		A[i].d->i = 0;						// 
		A[i].d->j = 0;						// 
		A[i].d->k = i;						// 
	}
}

maxheap::maxheap(int size) {
	tuple *newtuple;
	heaplimit = 1;							// first free location is A[1]
	arraysize = size+1;						// 
	isempty   = true;						// heap is initially empty
	A = new hnode [arraysize];				// allocate array for heap
	for (int i=0; i<arraysize; i++) {			// initialize heap values
		newtuple = new tuple;				// 
		A[i].d = newtuple;					// 
		A[i].d->m = -4294967296.0;			// a very negative value; unlikely to be a valid dQ
		A[i].d->i = 0;						// 
		A[i].d->j = 0;						// 
		A[i].d->k = i;						// 
	}
}

maxheap::~maxheap() {
	for (int i=0; i<arraysize; i++) { delete A[i].d; }	// first delete the containers pointed to
	delete [] A;									// then delete the list of pointers to containers
}

// Pop Maximum ------------------------------------------------------------------------
tuple maxheap::popMaximum() {				// O(log k) time
	tuple temp = returnMaximum();
	deleteItem(A[1].d);
	return temp;
}

// Return Maximum --------------------------------------------------------------------
tuple maxheap::returnMaximum() {			// O(1) time
	tuple temp;
	temp.m = A[1].d->m;					// 
	temp.i = A[1].d->i;					// 
	temp.j = A[1].d->j;					// 
	temp.k = A[1].d->k;					// grab A's data
	return temp;
}

int maxheap::returnArraysize() { return arraysize; }
int maxheap::returnHeaplimit() { return heaplimit; }


// Heapification functions ------------------------------------------------------------
int maxheap::downsift(int index) {			// O(log k) time
	bool stopFlag  = false;
	int L		= left(index);
	int R		= right(index);
	int swap;
	tuple* temp;
	
	while (!stopFlag) {
		// check that both children are within the array boundaries
		if ((L < heaplimit) && (R < heaplimit)) {
			if (A[L].d->m > A[R].d->m) { swap = L; } else { swap = R; } // first choose larger of the children
		} else { if (L < heaplimit) { swap = L; } else { break; } }		// only one child to consider
		
		// now decide if need to exchange A[index] with A[swap]
		if (A[index].d->m < A[swap].d->m) {
			temp          = A[index].d;   // exchange pointers A[index] and A[swap]
			A[index].d    = A[swap].d;    // 
			A[index].d->k = index;		// note A[index].d's change of array location
			A[swap].d     = temp;		// 
			A[swap].d->k  = swap;		// note A[swap].d's change in array location
			
			index = swap;				// update indices for next pass
			L     = left(index);		// 
			R     = right(index);		// 			
		} else { stopFlag = true; }
	}
	return index;						// return the new index location of downsifted element
}

int maxheap::upsift(int index) {			// O(log k) time
	bool stopFlag = false;
	int P         = parent(index);
	tuple* temp;
	while (!stopFlag) {
		// decide if A[index] needs to move up in tree
		if ((P > 0) && (A[index].d->m > A[P].d->m)) {
			temp          = A[index].d;   // exchange A[index] and A[P]
			A[index].d    = A[P].d;		// 
			A[index].d->k = index;		// note A[index].d's change of array location
			A[P].d        = temp;		// 
			A[P].d->k     = P;			// note A[P].d's change of array location
			
			index = P;				// update indices for next pass
			P     = parent(index);		// 
		} else { stopFlag = true; }
	}
	return index;
}


int  maxheap::left  (int index) { return 2*index;		}
int  maxheap::right (int index) { return 2*index+1;    }
int  maxheap::parent(int index) { return (int)index/2; }

void maxheap::grow() {								// O(k) time
	tuple	*newtuple;
	hnode	*B;									// scratch space for expansion of A
	B = new hnode [arraysize];						// 
	for (int i=0; i<arraysize; i++) { B[i].d = A[i].d; }   // copy A into B
	delete [] A;									// delete old array of addresses
	A = new hnode [2*arraysize];						// grow A by factor of 2
	for (int i=0; i<arraysize; i++) { A[i].d = B[i].d; }   // copy B into first half of A
	for (int i=arraysize; i<(2*arraysize); i++) {		// initialize new heap values
		newtuple  = new tuple;						// 
		A[i].d    = newtuple;						// 
		A[i].d->m = -4294967296.0;					// 
		A[i].d->i = 0;								// 
		A[i].d->j = 0;								// 
		A[i].d->k = i;								// 
	}										 
	delete [] B;									// delete scratch space B
	arraysize = 2*arraysize;							// update size of array
	return;
}

void maxheap::shrink() {								// O(k) time
	tuple	*newtuple;
	hnode	*B;									// scratch space for contraction of A
	B = new hnode [heaplimit];						// 
	for (int i=0; i<heaplimit; i++) { B[i].d = A[i].d; }   // copy A into B
	delete [] A;									// delete old array of addresses
	A = new hnode [arraysize/2];						// shrink A by factor of 2
	for (int i=0; i<heaplimit; i++) { A[i].d = B[i].d; }   // copy B into A
	for (int i=heaplimit; i<(arraysize/2); i++) {		// initialize new heap values
		newtuple  = new tuple;						// 
		A[i].d    = newtuple;						// 
		A[i].d->m = -4294967296.0;					// 
		A[i].d->i = 0;								// 
		A[i].d->j = 0;								// 
		A[i].d->k = i;								// 
	}										 
	delete [] B;									// delete scratch space B
	arraysize = arraysize/2;							// update size of array
	return;
}

// Insert + Update Functions --------------------------------------------------------------

tuple* maxheap::insertItem(const tuple newData) {			// O(log k) time
	int index;
	tuple *pointer;
	if (heaplimit >= (arraysize-1)) { grow(); }  // if heap is full, grow by factor of 2

	index		= heaplimit;			// 
	A[index].d->m  = newData.m;			// copy newData onto the bottom of the heap
	A[index].d->i  = newData.i;			// 
	A[index].d->j  = newData.j;			// 
	A[index].d->k  = index;				//
	pointer		= A[index].d;			// store pointer to container
	heaplimit++;						// note the larger heap
	upsift(index);						// upsift new item up to proper spot
	if  (heaplimit > 1) { isempty = false; }
	
	return pointer;
}

void maxheap::updateItem(tuple *address, tuple newData) {   // O(log k) time
	double oldm = address->m;
	address->m = newData.m;						//   udpate the dQ stored
	address->j = newData.j;						//   udpate the community to which to join
	if (oldm > newData.m) { downsift(address->k); }   //   downsift if old value was larger
	else				  { upsift(address->k);   }   //   upsift otherwise
	return;
}

void maxheap::updateItem(tuple *address, double newStored) {	// O(log k) time
	double oldm    = address->m;
	address->m	= newStored;						// udpate the dQ stored
	if (oldm > newStored) { downsift(address->k); }		// downsift if old value was larger
	else				  { upsift(address->k);   }		// upsift otherwise
	return;
}

void maxheap::deleteItem(tuple *address) {
	tuple *swap;
	int small = heaplimit-1;
	int index = address->k;

	if (heaplimit==2) {						// check if deleting last item in heap
		A[1].d->m		= -4294967296.0;		// zero out root data
		A[1].d->i		= 0;					// 
		A[1].d->j		= 0;					// 
		A[1].d->k		= 1;					// 
		isempty		= true;				// note the heap's emptiness
		heaplimit--;						// decrement size of heap to empty
	} else {
		A[index].d->m  = -4294967296.0;		// zero out the deleted item's data
		A[index].d->i	= 0;					// 
		A[index].d->j  = 0;					// 
		swap			= A[index].d;			// swap this element with item at end of heap
		A[index].d	= A[small].d;			// 
		A[small].d	= swap;				// 
		A[index].d->k  = index;				// note the change in locations
		A[small].d->k  = small;				// 
		heaplimit--;						// note change in heap size
		downsift(index);					// downsift moved element to new location; O(log k)
		if ((heaplimit/2 > heapmin) && (heaplimit == (arraysize/2))) { shrink(); } // shrink by factor of 2 if necessary
	}
	return;
}

bool maxheap::heapIsEmpty() { return isempty; }
int  maxheap::heapSize()    { return heaplimit-1; }

// ------------------------------------------------------------------------------------

#endif
