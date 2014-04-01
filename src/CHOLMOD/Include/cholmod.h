/* ========================================================================== */
/* === Include/cholmod.h ==================================================== */
/* ========================================================================== */

/* -----------------------------------------------------------------------------
 * CHOLMOD/Include/cholmod.h.
 * Copyright (C) 2005-2013, Univ. of Florida.  Author: Timothy A. Davis
 * CHOLMOD/Include/cholmod.h is licensed under Version 2.1 of the GNU
 * Lesser General Public License.  See lesser.txt for a text of the license.
 * CHOLMOD is also available under other licenses; contact authors for details.
 * http://www.suitesparse.com
 *
 * Portions of CHOLMOD (the Core and Partition Modules) are copyrighted by the
 * University of Florida.  The Modify Module is co-authored by William W.
 * Hager, Univ. of Florida.
 *
 * Acknowledgements:  this work was supported in part by the National Science
 * Foundation (NFS CCR-0203270 and DMS-9803599), and a grant from Sandia
 * National Laboratories (Dept. of Energy) which supported the development of
 * CHOLMOD's Partition Module.
 * -------------------------------------------------------------------------- */

/* CHOLMOD include file, for inclusion user programs.
 *
 * The include files listed below include a short description of each user-
 * callable routine.  Each routine in CHOLMOD has a consistent interface.
 * More details about the CHOLMOD data types is in the cholmod_core.h file.
 *
 * Naming convention:
 * ------------------
 *
 *	All routine names, data types, and CHOLMOD library files use the
 *	cholmod_ prefix.  All macros and other #define's use the CHOLMOD
 *	prefix.
 * 
 * Return value:
 * -------------
 *
 *	Most CHOLMOD routines return an int (TRUE (1) if successful, or FALSE
 *	(0) otherwise.  A SuiteSparse_long or double return value is >= 0 if
 *	successful, or -1 otherwise.  A size_t return value is > 0 if
 *	successful, or 0 otherwise.
 *
 *	If a routine returns a pointer, it is a pointer to a newly allocated
 *	object or NULL if a failure occured, with one exception.  cholmod_free
 *	always returns NULL.
 *
 * "Common" parameter:
 * ------------------
 *
 *	The last parameter in all CHOLMOD routines is a pointer to the CHOLMOD
 *	"Common" object.  This contains control parameters, statistics, and
 *	workspace used between calls to CHOLMOD.  It is always an input/output
 *	parameter.
 *
 * Input, Output, and Input/Output parameters:
 * -------------------------------------------
 *
 *	Input parameters are listed first.  They are not modified by CHOLMOD.
 *
 *	Input/output are listed next.  They must be defined on input, and
 *	are modified on output.
 *
 *	Output parameters are listed next.  If they are pointers, they must
 *	point to allocated space on input, but their contents are not defined
 *	on input.
 *
 *	Workspace parameters appear next.  They are used in only two routines
 *	in the Supernodal module.
 *
 *	The cholmod_common *Common parameter always appears as the last
 *	parameter.  It is always an input/output parameter.
 */

#ifndef CHOLMOD_H
#define CHOLMOD_H

/* make it easy for C++ programs to include CHOLMOD */
#ifdef __cplusplus
extern "C" {
#endif

/* assume large file support.  If problems occur, compile with -DNLARGEFILE */
#include "cholmod_io64.h"

#include "SuiteSparse_config.h"

#include "cholmod_config.h"

/* CHOLMOD always includes the Core module. */
#include "cholmod_core.h"

#ifndef NCHECK
#include "cholmod_check.h"
#endif

#ifndef NCHOLESKY
#include "cholmod_cholesky.h"
#endif

#ifndef NMATRIXOPS
#include "cholmod_matrixops.h"
#endif

#ifndef NMODIFY
#include "cholmod_modify.h"
#endif

#ifndef NCAMD
#include "cholmod_camd.h"
#endif

#ifndef NPARTITION
#include "cholmod_partition.h"
#endif

#ifndef NSUPERNODAL
#include "cholmod_supernodal.h"
#endif

#ifdef __cplusplus
}
#endif

#endif
