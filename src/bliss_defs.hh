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

#ifndef BLISS_DEFS_HH
#define BLISS_DEFS_HH

#include "config.h"
#include <cassert>

/* Define this if you have gmp and want to have exact group sizes.
 * Remember to include -lgmp in LIB in Makefile. */
#if HAVE_GMP == 1
#  define BLISS_USE_GMP
#endif

#if defined(DEBUG)
#define CONSISTENCY_CHECKS
#define EXPENSIVE_CONSISTENCY_CHECKS
#endif
//#define PRINT_SEARCH_TREE_DOT

/* Force a check that the found automorphisms are valid */
#if defined(CONSISTENCY_CHECKS)
#define VERIFY_AUTOMORPHISMS
/* Force a check that the generated partitions are equitable */
#define VERIFY_EQUITABLEDNESS
#endif


#if defined(CONSISTENCY_CHECKS)
#define DEBUG_ASSERT(a) assert(a)
#else
#define DEBUG_ASSERT(a) ;
#endif

#endif
