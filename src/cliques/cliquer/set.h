
/*
 * This file contains the set handling routines.
 *
 * Copyright (C) 2002 Sampo Niskanen, Patric Östergård.
 * Licensed under the GNU GPL, read the file LICENSE for details.
 */

#ifndef CLIQUER_SET_H
#define CLIQUER_SET_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include "misc.h"

/*
 * Sets are arrays of setelement's (typically unsigned long int's) with
 * representative bits for each value they can contain.  The values
 * are numbered 0,...,n-1.
 */


/*** Variable types and constants. ***/


/*
 * If setelement hasn't been declared:
 *   - use "unsigned long int" as setelement
 *   - try to deduce size from ULONG_MAX
 */

#ifndef ELEMENTSIZE
typedef unsigned long int setelement;
# if (ULONG_MAX == 65535)
#  define ELEMENTSIZE 16
# elif (ULONG_MAX == 4294967295)
#  define ELEMENTSIZE 32
# else
#  define ELEMENTSIZE 64
# endif
#endif  /* !ELEMENTSIZE */

typedef setelement * set_t;


/*** Counting amount of 1 bits in a setelement ***/

/* Array for amount of 1 bits in a byte. */
static int set_bit_count[256] = {
	0,1,1,2,1,2,2,3,1,2,2,3,2,3,3,4,
	1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,
	1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,
	2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
	1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,
	2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
	2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
	3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
	1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,
	2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
	2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
	3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
	2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
	3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
	3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
	4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8 };

/* The following macros assume that all higher bits are 0.
 * They may in some cases be useful also on with other ELEMENTSIZE's,
 * so we define them all.  */
#define SET_ELEMENT_BIT_COUNT_8(a)  (set_bit_count[(a)])
#define SET_ELEMENT_BIT_COUNT_16(a) (set_bit_count[(a)>>8] + \
				     set_bit_count[(a)&0xFF])
#define SET_ELEMENT_BIT_COUNT_32(a) (set_bit_count[(a)>>24] + \
				     set_bit_count[((a)>>16)&0xFF] + \
				     set_bit_count[((a)>>8)&0xFF] + \
				     set_bit_count[(a)&0xFF])
#define SET_ELEMENT_BIT_COUNT_64(a) (set_bit_count[(a)>>56] + \
				     set_bit_count[((a)>>48)&0xFF] + \
				     set_bit_count[((a)>>40)&0xFF] + \
				     set_bit_count[((a)>>32)&0xFF] + \
				     set_bit_count[((a)>>24)&0xFF] + \
				     set_bit_count[((a)>>16)&0xFF] + \
				     set_bit_count[((a)>>8)&0xFF] + \
				     set_bit_count[(a)&0xFF])
#if (ELEMENTSIZE==64)
# define SET_ELEMENT_BIT_COUNT(a) SET_ELEMENT_BIT_COUNT_64(a)
# define FULL_ELEMENT ((setelement)0xFFFFFFFFFFFFFFFF)
#elif (ELEMENTSIZE==32)
# define SET_ELEMENT_BIT_COUNT(a) SET_ELEMENT_BIT_COUNT_32(a)
# define FULL_ELEMENT ((setelement)0xFFFFFFFF)
#elif (ELEMENTSIZE==16)
# define SET_ELEMENT_BIT_COUNT(a) SET_ELEMENT_BIT_COUNT_16(a)
# define FULL_ELEMENT ((setelement)0xFFFF)
#else
# error "SET_ELEMENT_BIT_COUNT(a) not defined for current ELEMENTSIZE"
#endif



/*** Macros and functions ***/

/*
 * Gives a value with bit x (counting from lsb up) set.
 *
 * Making this as a table might speed up things on some machines
 * (though on most modern machines it's faster to shift instead of
 * using memory).  Making it a macro makes it easy to change.
 */
#define SET_BIT_MASK(x) ((setelement)1<<(x))



/* Set element handling macros */

#define SET_ELEMENT_INTERSECT(a,b)  ((a)&(b))
#define SET_ELEMENT_UNION(a,b)      ((a)|(b))
#define SET_ELEMENT_DIFFERENCE(a,b) ((a)&(~(b)))
#define SET_ELEMENT_CONTAINS(e,v)   ((e)&SET_BIT_MASK(v))


/* Set handling macros */

#define SET_ADD_ELEMENT(s,a) \
                       ((s)[(a)/ELEMENTSIZE] |= SET_BIT_MASK((a)%ELEMENTSIZE))
#define SET_DEL_ELEMENT(s,a) \
                       ((s)[(a)/ELEMENTSIZE] &= ~SET_BIT_MASK((a)%ELEMENTSIZE))
#define SET_CONTAINS_FAST(s,a) (SET_ELEMENT_CONTAINS((s)[(a)/ELEMENTSIZE], \
						      (a)%ELEMENTSIZE))
#define SET_CONTAINS(s,a) (((a)<SET_MAX_SIZE(s))?SET_CONTAINS_FAST(s,a):FALSE)

/* Sets can hold values between 0,...,SET_MAX_SIZE(s)-1 */
#define SET_MAX_SIZE(s) ((s)[-1])
/* Sets consist of an array of SET_ARRAY_LENGTH(s) setelements */
#define SET_ARRAY_LENGTH(s) (((s)[-1]+ELEMENTSIZE-1)/ELEMENTSIZE)


/*
 * set_new()
 *
 * Create a new set that can hold values in the range 0,...,size-1.
 */
UNUSED_FUNCTION
static set_t set_new(int size) {
	int n;
	set_t s;

	ASSERT(size>0);

	n=(size/ELEMENTSIZE+1)+1;
	s=calloc(n,sizeof(setelement));
	s[0]=size;

	return &(s[1]);
}

/*
 * set_free()
 *
 * Free the memory associated with set s.
 */
UNUSED_FUNCTION INLINE
static void set_free(set_t s) {
	ASSERT(s!=NULL);
	free(&(s[-1]));
}

/*
 * set_resize()
 *
 * Resizes set s to given size.  If the size is less than SET_MAX_SIZE(s),
 * the last elements are dropped.
 *
 * Returns a pointer to the new set.
 */
UNUSED_FUNCTION INLINE
static set_t set_resize(set_t s, unsigned int size) {
	unsigned int n;

	ASSERT(size>0);

	n=(size/ELEMENTSIZE+1);
	s=((setelement *)realloc(s-1,(n+1)*sizeof(setelement)))+1;

	if (n>SET_ARRAY_LENGTH(s))
		memset(s+SET_ARRAY_LENGTH(s),0,
		       (n-SET_ARRAY_LENGTH(s))*sizeof(setelement));
	if (size < SET_MAX_SIZE(s))
		s[(size-1)/ELEMENTSIZE] &= (FULL_ELEMENT >>
					    (ELEMENTSIZE-size%ELEMENTSIZE));
	s[-1]=size;

	return s;
}

/*
 * set_size()
 *
 * Returns the number of elements in set s.
 */
UNUSED_FUNCTION INLINE
static int set_size(set_t s) {
	int count=0;
	setelement *c;

	for (c=s; c < s+SET_ARRAY_LENGTH(s); c++)
		count+=SET_ELEMENT_BIT_COUNT(*c);
	return count;
}

/*
 * set_duplicate()
 *
 * Returns a newly allocated duplicate of set s.
 */
UNUSED_FUNCTION INLINE
static set_t set_duplicate(set_t s) {
	set_t new;

	new=set_new(SET_MAX_SIZE(s));
	memcpy(new,s,SET_ARRAY_LENGTH(s)*sizeof(setelement));
	return new;
}

/*
 * set_copy()
 *
 * Copies set src to dest.  If dest is NULL, is equal to set_duplicate.
 * If dest smaller than src, it is freed and a new set of the same size as
 * src is returned.
 */
UNUSED_FUNCTION INLINE
static set_t set_copy(set_t dest,set_t src) {
	if (dest==NULL)
		return set_duplicate(src);
	if (SET_MAX_SIZE(dest)<SET_MAX_SIZE(src)) {
		set_free(dest);
		return set_duplicate(src);
	}
	memcpy(dest,src,SET_ARRAY_LENGTH(src)*sizeof(setelement));
	memset(dest+SET_ARRAY_LENGTH(src),0,((SET_ARRAY_LENGTH(dest) -
					      SET_ARRAY_LENGTH(src)) *
					     sizeof(setelement)));
	return dest;
}

/*
 * set_empty()
 *
 * Removes all elements from the set s.
 */
UNUSED_FUNCTION INLINE
static void set_empty(set_t s) {
	memset(s,0,SET_ARRAY_LENGTH(s)*sizeof(setelement));
	return;
}

/*
 * set_intersection()
 *
 * Store the intersection of sets a and b into res.  If res is NULL,
 * a new set is created and the result is written to it.  If res is
 * smaller than the larger one of a and b, it is freed and a new set
 * is created and the result is returned.
 *
 * Returns either res or a new set that has been allocated in its stead.
 *
 * Note:  res may not be a or b.
 */
UNUSED_FUNCTION INLINE
static set_t set_intersection(set_t res,set_t a,set_t b) {
	int i,max;

	if (res==NULL) {
		res = set_new(MAX(SET_MAX_SIZE(a),SET_MAX_SIZE(b)));
	} else if (SET_MAX_SIZE(res) < MAX(SET_MAX_SIZE(a),SET_MAX_SIZE(b))) {
		set_free(res);
		res = set_new(MAX(SET_MAX_SIZE(a),SET_MAX_SIZE(b)));
	} else {
		set_empty(res);
	}

	max=MIN(SET_ARRAY_LENGTH(a),SET_ARRAY_LENGTH(b));
	for (i=0; i<max; i++) {
		res[i]=SET_ELEMENT_INTERSECT(a[i],b[i]);
	}

	return res;
}

/*
 * set_union()
 *
 * Store the union of sets a and b into res.  If res is NULL, a new set
 * is created and the result is written to it.  If res is smaller than
 * the larger one of a and b, it is freed and a new set is created and
 * the result is returned.
 *
 * Returns either res or a new set that has been allocated in its stead.
 *
 * Note:  res may not be a or b.
 */
UNUSED_FUNCTION INLINE
static set_t set_union(set_t res,set_t a,set_t b) {
	int i,max;

	if (res==NULL) {
		res = set_new(MAX(SET_MAX_SIZE(a),SET_MAX_SIZE(b)));
	} else if (SET_MAX_SIZE(res) < MAX(SET_MAX_SIZE(a),SET_MAX_SIZE(b))) {
		set_free(res);
		res = set_new(MAX(SET_MAX_SIZE(a),SET_MAX_SIZE(b)));
	} else {
		set_empty(res);
	}

	max=MAX(SET_ARRAY_LENGTH(a),SET_ARRAY_LENGTH(b));
	for (i=0; i<max; i++) {
		res[i]=SET_ELEMENT_UNION(a[i],b[i]);
	}

	return res;
}


/*
 * set_return_next()
 *
 * Returns the smallest value in set s which is greater than n, or -1 if
 * such a value does not exist.
 *
 * Can be used to iterate through all values of s:
 *
 * int i=-1;
 * while ((i=set_return_next(s,i))>=0) {
 *         // i is in set s
 * }
 */
UNUSED_FUNCTION INLINE
static int set_return_next(set_t s, unsigned int n) {
	n++;
	if (n >= SET_MAX_SIZE(s))
		return -1;

	while (n%ELEMENTSIZE) {
		if (SET_CONTAINS(s,n))
			return n;
		n++;
		if (n >= SET_MAX_SIZE(s))
			return -1;
	}

	while (s[n/ELEMENTSIZE]==0) {
		n+=ELEMENTSIZE;
		if (n >= SET_MAX_SIZE(s))
			return -1;
	}
	while (!SET_CONTAINS(s,n)) {
		n++;
		if (n >= SET_MAX_SIZE(s))
			return -1;
	}
	return n;
}


/*
 * set_print()
 *
 * Prints the size and contents of set s to stdout.
 * Mainly useful for debugging purposes and trivial output.
 */
/*
UNUSED_FUNCTION
static void set_print(set_t s) {
	int i;
	printf("size=%d(max %d)",set_size(s),(int)SET_MAX_SIZE(s));
	for (i=0; i<SET_MAX_SIZE(s); i++)
		if (SET_CONTAINS(s,i))
			printf(" %d",i);
	printf("\n");
	return;
}
*/

#endif /* !CLIQUER_SET_H */
