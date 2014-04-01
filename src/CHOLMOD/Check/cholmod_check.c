/* ========================================================================== */
/* === Check/cholmod_check ================================================== */
/* ========================================================================== */

/* -----------------------------------------------------------------------------
 * CHOLMOD/Check Module.  Copyright (C) 2005-2013, Timothy A. Davis
 * The CHOLMOD/Check Module is licensed under Version 2.1 of the GNU
 * Lesser General Public License.  See lesser.txt for a text of the license.
 * CHOLMOD is also available under other licenses; contact authors for details.
 * -------------------------------------------------------------------------- */

/* Routines to check and print the contents of the 5 CHOLMOD objects:
 *
 * No CHOLMOD routine calls the check or print routines.  If a user wants to
 * check CHOLMOD's input parameters, a separate call to the appropriate check
 * routine should be used before calling other CHOLMOD routines.
 *
 * cholmod_check_common		check statistics and workspace in Common
 * cholmod_check_sparse		check sparse matrix in compressed column form
 * cholmod_check_dense		check dense matrix
 * cholmod_check_factor		check factorization
 * cholmod_check_triplet	check sparse matrix in triplet form
 *
 * cholmod_print_common		print statistics in Common
 * cholmod_print_sparse		print sparse matrix in compressed column form
 * cholmod_print_dense		print dense matrix
 * cholmod_print_factor		print factorization
 * cholmod_print_triplet	print sparse matrix in triplet form
 *
 * In addition, this file contains routines to check and print three types of
 * integer vectors:
 * 
 * cholmod_check_perm		check a permutation of 0:n-1 (no duplicates)
 * cholmod_check_subset		check a subset of 0:n-1 (duplicates OK)
 * cholmod_check_parent		check an elimination tree
 *
 * cholmod_print_perm		print a permutation
 * cholmod_print_subset		print a subset
 * cholmod_print_parent		print an elimination tree
 *
 * Each Common->print level prints the items at or below the given level:
 *
 *	0: print nothing; just check the data structures and return TRUE/FALSE
 *	1: error messages
 *	2: warning messages
 *	3: one-line summary of each object printed
 *	4: short summary of each object (first and last few entries)
 *	5: entire contents of the object
 *
 * No CHOLMOD routine calls these routines, so no printing occurs unless
 * the user specifically calls a cholmod_print_* routine.  Thus, the default
 * print level is 3.
 *
 * Common->precise controls the # of digits printed for numerical entries
 * (5 if FALSE, 15 if TRUE).
 *
 * If Common->print_function is NULL, then no printing occurs.  The
 * cholmod_check_* and cholmod_print_* routines still check their inputs and
 * return TRUE/FALSE if the object is valid or not.
 *
 * This file also includes debugging routines that are enabled only when
 * NDEBUG is defined in cholmod_internal.h (cholmod_dump_*).
 */

#ifndef NCHECK

#include "cholmod_internal.h"
#include "cholmod_check.h"

/* ========================================================================== */
/* === printing definitions ================================================= */
/* ========================================================================== */

#ifdef LONG
#define I8 "%8ld"
#define I_8 "%-8ld"
#else
#define I8 "%8d"
#define I_8 "%-8d"
#endif

#define PR(i,format,arg) \
{ \
    if (print >= i && Common->print_function != NULL) \
    { \
	(Common->print_function) (format, arg) ; \
    } \
}

#define P1(format,arg) PR(1,format,arg)
#define P2(format,arg) PR(2,format,arg)
#define P3(format,arg) PR(3,format,arg)
#define P4(format,arg) PR(4,format,arg)

#define ERR(msg) \
{ \
    P1 ("\nCHOLMOD ERROR: %s: ", type) ; \
    if (name != NULL) \
    { \
	P1 ("%s", name) ; \
    } \
    P1 (": %s\n", msg) ; \
    ERROR (CHOLMOD_INVALID, "invalid") ; \
    return (FALSE) ; \
}

/* print a numerical value */
#define PRINTVALUE(value) \
{ \
    if (Common->precise) \
    { \
	P4 (" %23.15e", value) ; \
    } \
    else \
    { \
	P4 (" %.5g", value) ; \
    } \
}

/* start printing */
#define ETC_START(count,limit) \
{ \
    count = (init_print == 4) ? (limit) : (-1) ; \
}

/* re-enable printing if condition is met */
#define ETC_ENABLE(condition,count,limit) \
{ \
    if ((condition) && init_print == 4) \
    { \
	count = limit ; \
	print = 4 ; \
    } \
}

/* turn off printing if limit is reached */
#define ETC_DISABLE(count) \
{ \
    if ((count >= 0) && (count-- == 0) && print == 4) \
    { \
	P4 ("%s", "    ...\n")  ; \
	print = 3 ; \
    } \
}

/* re-enable printing, or turn if off after limit is reached */
#define ETC(condition,count,limit) \
{ \
    ETC_ENABLE (condition, count, limit) ; \
    ETC_DISABLE (count) ; \
}

#define BOOLSTR(x) ((x) ? "true " : "false")

/* ========================================================================== */
/* === print_value ========================================================== */
/* ========================================================================== */

static void print_value
(
    Int print,
    Int xtype,
    double *Xx,
    double *Xz,
    Int p,
    cholmod_common *Common)
{
    if (xtype == CHOLMOD_REAL)
    {
	PRINTVALUE (Xx [p]) ;
    }
    else if (xtype == CHOLMOD_COMPLEX)
    {
	P4 ("%s", "(") ;
	PRINTVALUE (Xx [2*p  ]) ;
	P4 ("%s", " , ") ;
	PRINTVALUE (Xx [2*p+1]) ;
	P4 ("%s", ")") ;
    }
    else if (xtype == CHOLMOD_ZOMPLEX)
    {
	P4 ("%s", "(") ;
	PRINTVALUE (Xx [p]) ;
	P4 ("%s", " , ") ;
	PRINTVALUE (Xz [p]) ;
	P4 ("%s", ")") ;
    }
}

/* ========================================================================== */
/* === cholmod_check_common ================================================= */
/* ========================================================================== */

/* Print and verify the contents of Common */

static int check_common
(
    Int print,
    const char *name,
    cholmod_common *Common
)
{
    double fl, lnz ;
    double *Xwork ;
    Int *Flag, *Head ;
    SuiteSparse_long mark ;
    Int i, nrow, nmethods, ordering, xworksize, amd_backup, init_print ;
    const char *type = "common" ;

    /* ---------------------------------------------------------------------- */
    /* print control parameters and statistics */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (FALSE) ;
    init_print = print ;

    P2 ("%s", "\n") ;

    P1 ("CHOLMOD version %d", CHOLMOD_MAIN_VERSION) ;
    P1 (".%d", CHOLMOD_SUB_VERSION) ;
    P1 (".%d", CHOLMOD_SUBSUB_VERSION) ;
    P1 (", %s: ", CHOLMOD_DATE) ;

    if (name != NULL)
    {
	P1 ("%s: ", name) ;
    }
    switch (Common->status)
    {

	case CHOLMOD_OK:
	    P1 ("%s", "status: OK\n") ;
	    break ;

	case CHOLMOD_OUT_OF_MEMORY:
	    P1 ("%s", "status: ERROR, out of memory\n") ;
	    break ;

	case CHOLMOD_INVALID:
	    P1 ("%s", "status: ERROR, invalid parameter\n") ;
	    break ;

	case CHOLMOD_TOO_LARGE:
	    P1 ("%s", "status: ERROR, problem too large\n") ;
	    break ;

	case CHOLMOD_NOT_INSTALLED:
	    P1 ("%s", "status: ERROR, method not installed\n") ;
	    break ;

#if GPU_BLAS
	case CHOLMOD_GPU_PROBLEM:
	    P1 ("%s", "status: ERROR, GPU had a fatal error\n") ;
	    break ;
#endif

	case CHOLMOD_NOT_POSDEF:
	    P1 ("%s", "status: warning, matrix not positive definite\n") ;
	    break ;

	case CHOLMOD_DSMALL:
	    P1 ("%s", "status: warning, diagonal entry has tiny abs. value\n") ;
	    break ;

	default:
	    ERR ("unknown status") ;
    }

    P2 ("  Architecture: %s\n", CHOLMOD_ARCHITECTURE) ;
    P3 ("    sizeof(int):      %d\n", (int) sizeof (int)) ;
    P3 ("    sizeof(SuiteSparse_long):  %d\n", (int) sizeof (SuiteSparse_long));
    P3 ("    sizeof(void *):   %d\n", (int) sizeof (void *)) ;
    P3 ("    sizeof(double):   %d\n", (int) sizeof (double)) ;
    P3 ("    sizeof(Int):      %d (CHOLMOD's basic integer)\n", (int) sizeof (Int)) ;
    P3 ("    sizeof(BLAS_INT): %d (integer used in the BLAS)\n",
	    (int) sizeof (BLAS_INT)) ;

    if (Common->fl != EMPTY)
    {
	P2 ("%s", "  Results from most recent analysis:\n") ;
	P2 ("    Cholesky flop count: %.5g\n", Common->fl) ;
	P2 ("    Nonzeros in L:       %.5g\n", Common->lnz) ;
    }
    if (Common->modfl != EMPTY)
    {
	P2 ("    Update/downdate flop count: %.5g\n", Common->modfl) ;
    }

    P2 ("  memory blocks in use:    %8.0f\n", (double) (Common->malloc_count)) ;
    P2 ("  memory in use (MB):      %8.1f\n", 
	(double) (Common->memory_inuse) / 1048576.) ;
    P2 ("  peak memory usage (MB):  %8.1f\n", 
	(double) (Common->memory_usage) / 1048576.) ;

    /* ---------------------------------------------------------------------- */
    /* primary control parameters and related ordering statistics */
    /* ---------------------------------------------------------------------- */

    P3 ("  maxrank:    update/downdate rank:   "ID"\n",
	    (Int) CHOLMOD(maxrank) (0, Common)) ;
    P3 ("  supernodal control: %d", Common->supernodal) ;
    P3 (" %g ", Common->supernodal_switch) ;
    if (Common->supernodal <= CHOLMOD_SIMPLICIAL)
    {
	P3 ("%s", "(always do simplicial)\n") ;
    }
    else if (Common->supernodal == CHOLMOD_AUTO)
    {
	P3 ("(supernodal if flops/lnz >= %g)\n", Common->supernodal_switch) ;
    }
    else
    {
	P3 ("%s", "(always do supernodal)\n") ;
    }

    nmethods = MIN (Common->nmethods, CHOLMOD_MAXMETHODS) ;
    nmethods = MAX (0, nmethods) ;

    if (nmethods > 0)
    {
	P3 ("%s", "  nmethods:   number of ordering methods to try: ") ;
	P3 (""ID"\n", nmethods) ;
        amd_backup = (nmethods > 1) || (nmethods == 1 &&
            (Common->method [0].ordering == CHOLMOD_METIS ||
             Common->method [0].ordering == CHOLMOD_NESDIS)) ;
    }
    else
    {
	P3 ("%s", "  nmethods=0: default strategy:  Try user permutation if "
		"given.  Try AMD.\n") ;
#ifndef NPARTITION
	if (Common->default_nesdis)
	{
	    P3 ("%s", "    Try NESDIS if AMD reports flops/nnz(L) >= 500 and "
		"nnz(L)/nnz(A) >= 5.\n") ;
	}
	else
	{
	    P3 ("%s", "    Try METIS if AMD reports flops/nnz(L) >= 500 and "
		"nnz(L)/nnz(A) >= 5.\n") ;
	}
#endif
	P3 ("%s", "    Select best ordering tried.\n") ;
	Common->method [0].ordering = CHOLMOD_GIVEN ;
	Common->method [1].ordering = CHOLMOD_AMD ;
	Common->method [2].ordering = 
            (Common->default_nesdis ? CHOLMOD_NESDIS : CHOLMOD_METIS) ;
        amd_backup = FALSE ;
#ifndef NPARTITION
	nmethods = 3 ;
#else
	nmethods = 2 ;
#endif
    }

    for (i = 0 ; i < nmethods ; i++)
    {
	P3 ("    method "ID": ", i) ;
	ordering = Common->method [i].ordering ;
	fl = Common->method [i].fl ;
	lnz = Common->method [i].lnz ;
	switch (ordering)
	{

	    case CHOLMOD_NATURAL:
		P3 ("%s", "natural\n") ;
		break ;

	    case CHOLMOD_GIVEN:
		P3 ("%s", "user permutation (if given)\n") ;
		break ;

	    case CHOLMOD_AMD:
		P3 ("%s", "AMD (or COLAMD if factorizing AA')\n") ;
		amd_backup = FALSE ;
		break ;

	    case CHOLMOD_COLAMD:
		P3 ("%s", "AMD if factorizing A, COLAMD if factorizing AA')\n");
		amd_backup = FALSE ;
		break ;

	    case CHOLMOD_METIS:
		P3 ("%s", "METIS_NodeND nested dissection\n") ;
		break ;

	    case CHOLMOD_NESDIS:
		P3 ("%s", "CHOLMOD nested dissection\n") ;

		P3 ("        nd_small: # nodes in uncut subgraph: "ID"\n",
			(Int) (Common->method [i].nd_small)) ;
		P3 ("        nd_compress: compress the graph:     %s\n",
			BOOLSTR (Common->method [i].nd_compress)) ;
		P3 ("        nd_camd: use constrained min degree: %s\n",
			BOOLSTR (Common->method [i].nd_camd)) ;
		break ;

	    default:
		P3 (ID, ordering) ;
		ERR ("unknown ordering method") ;
		break ;

	}

	if (!(ordering == CHOLMOD_NATURAL || ordering == CHOLMOD_GIVEN))
	{
	    if (Common->method [i].prune_dense < 0)
	    {
		P3 ("        prune_dense: for pruning dense nodes:   %s\n",
			" none pruned") ;
	    }
	    else
	    {
		P3 ("        prune_dense: for pruning dense nodes:   "
		    "%.5g\n",
		    Common->method [i].prune_dense) ;
		P3 ("        a dense node has degree "
			">= max(16,(%.5g)*sqrt(n))\n",
		    Common->method [i].prune_dense) ;
	    }
	}

	if (ordering == CHOLMOD_COLAMD || ordering == CHOLMOD_NESDIS)
	{
	    if (Common->method [i].prune_dense2 < 0)
	    {
		P3 ("        prune_dense2: for pruning dense rows for AA':"
			"  %s\n", " none pruned") ;
	    }
	    else
	    {
		P3 ("        prune_dense2: for pruning dense rows for AA':"
		    " %.5g\n", Common->method [i].prune_dense2) ;
		P3 ("        a dense row has degree "
			">= max(16,(%.5g)*sqrt(ncol))\n",
		    Common->method [i].prune_dense2) ;
	    }
	}

	if (fl  != EMPTY) P3 ("        flop count: %.5g\n", fl) ;
	if (lnz != EMPTY) P3 ("        nnz(L):     %.5g\n", lnz) ;
    }

    /* backup AMD results, if any */
    if (amd_backup)
    {
	P3 ("%s", "    backup method: ") ;
	P3 ("%s", "AMD (or COLAMD if factorizing AA')\n") ;
	fl = Common->method [nmethods].fl ;
	lnz = Common->method [nmethods].lnz ;
	if (fl  != EMPTY) P3 ("        AMD flop count: %.5g\n", fl) ;
	if (lnz != EMPTY) P3 ("        AMD nnz(L):     %.5g\n", lnz) ;
    }

    /* ---------------------------------------------------------------------- */
    /* arcane control parameters */
    /* ---------------------------------------------------------------------- */

    if (Common->final_asis)
    {
	P4 ("%s", "  final_asis: TRUE, leave as is\n") ;
    }
    else
    {
	P4 ("%s", "  final_asis: FALSE, convert when done\n") ;
	if (Common->final_super)
	{
	    P4 ("%s", "  final_super: TRUE, leave in supernodal form\n") ;
	}
	else
	{
	    P4 ("%s", "  final_super: FALSE, convert to simplicial form\n") ;
	}
	if (Common->final_ll)
	{
	    P4 ("%s", "  final_ll: TRUE, convert to LL' form\n") ;
	}
	else
	{
	    P4 ("%s", "  final_ll: FALSE, convert to LDL' form\n") ;
	}
	if (Common->final_pack)
	{
	    P4 ("%s", "  final_pack: TRUE, pack when done\n") ;
	}
	else
	{
	    P4 ("%s", "  final_pack: FALSE, do not pack when done\n") ;
	}
	if (Common->final_monotonic)
	{
	    P4 ("%s", "  final_monotonic: TRUE, ensure L is monotonic\n") ;
	}
	else
	{
	    P4 ("%s",
		"  final_monotonic: FALSE, do not ensure L is monotonic\n") ;
	}
	P4 ("  final_resymbol: remove zeros from amalgamation: %s\n",
		BOOLSTR (Common->final_resymbol)) ;
    }

    P4 ("  dbound:  LDL' diagonal threshold: % .5g\n    Entries with abs. value"
	    " less than dbound are replaced with +/- dbound.\n",
	    Common->dbound) ;

    P4 ("  grow0: memory reallocation: % .5g\n", Common->grow0) ;
    P4 ("  grow1: memory reallocation: % .5g\n", Common->grow1) ;
    P4 ("  grow2: memory reallocation: %g\n", (double) (Common->grow2)) ;

    P4 ("%s", "  nrelax, zrelax:  supernodal amalgamation rule:\n") ;
    P4 ("%s", "    s = # columns in two adjacent supernodes\n") ;
    P4 ("%s", "    z = % of zeros in new supernode if they are merged.\n") ;
    P4 ("%s", "    Two supernodes are merged if") ;
    P4 (" (s <= %g) or (no new zero entries) or\n",
	    (double) (Common->nrelax [0])) ;
    P4 ("    (s <= %g and ",  (double) (Common->nrelax [1])) ;
    P4 ("z < %.5g%%) or",      Common->zrelax [0] * 100) ;
    P4 (" (s <= %g and ",     (double) (Common->nrelax [2])) ;
    P4 ("z < %.5g%%) or",      Common->zrelax [1] * 100) ;
    P4 (" (z < %.5g%%)\n",     Common->zrelax [2] * 100) ;

    /* ---------------------------------------------------------------------- */
    /* check workspace */
    /* ---------------------------------------------------------------------- */

    mark = Common->mark ;
    nrow = Common->nrow ;
    Flag = Common->Flag ;
    Head = Common->Head ;
    if (nrow > 0)
    {
	if (mark < 0 || Flag == NULL || Head == NULL)
	{
	    ERR ("workspace corrupted (Flag and/or Head missing)") ;
	}
	for (i = 0 ; i < nrow ; i++)
	{
	    if (Flag [i] >= mark)
	    {
		PRINT0 (("Flag ["ID"]="ID", mark = %ld\n", i, Flag [i], mark)) ;
		ERR ("workspace corrupted (Flag)") ;
	    }
	}
	for (i = 0 ; i <= nrow ; i++)
	{
	    if (Head [i] != EMPTY)
	    {
		PRINT0 (("Head ["ID"] = "ID",\n", i, Head [i])) ;
		ERR ("workspace corrupted (Head)") ;
	    }
	}
    }
    xworksize = Common->xworksize ;
    Xwork = Common->Xwork ;
    if (xworksize > 0)
    {
	if (Xwork == NULL)
	{
	    ERR ("workspace corrupted (Xwork missing)") ;
	}
	for (i = 0 ; i < xworksize ; i++)
	{
	    if (Xwork [i] != 0.)
	    {
		PRINT0 (("Xwork ["ID"] = %g\n", i, Xwork [i])) ;
		ERR ("workspace corrupted (Xwork)") ;
	    }
	}
    }

    /* workspace and parameters are valid */
    P3 ("%s", "  OK\n") ;
    P4 ("%s", "\n") ;
    return (TRUE) ;
}


int CHOLMOD(check_common)
(
    cholmod_common *Common
)
{
    return (check_common (0, NULL, Common)) ;
}


int CHOLMOD(print_common)
(
    /* ---- input ---- */
    const char *name,		/* printed name of Common object */
    /* --------------- */
    cholmod_common *Common
)
{
    Int print = (Common == NULL) ? 3 : (Common->print) ;
    return (check_common (print, name, Common)) ;
}


/* ========================================================================== */
/* === cholmod_gpu_stats ==================================================== */
/* ========================================================================== */

/* Print CPU / GPU statistics.  If the timer is not installed, the times are
   reported as zero, but this function still works.  Likewise, the function
   still works if the GPU BLAS is not installed. */

int CHOLMOD(gpu_stats)
(
    cholmod_common *Common      /* input */
)
{
    double cpu_time, gpu_time ;
    int print ;

    RETURN_IF_NULL_COMMON (FALSE) ;
    print = Common->print ;

    P2 ("%s", "\nCHOLMOD GPU/CPU statistics:\n") ;
    P2 ("SYRK  CPU calls %12.0f", (double) Common->CHOLMOD_CPU_SYRK_CALLS) ;
    P2 (" time %12.4e\n", Common->CHOLMOD_CPU_SYRK_TIME) ;
    P2 ("      GPU calls %12.0f", (double) Common->CHOLMOD_GPU_SYRK_CALLS) ;
    P2 (" time %12.4e\n", Common->CHOLMOD_GPU_SYRK_TIME) ;
    P2 ("GEMM  CPU calls %12.0f", (double) Common->CHOLMOD_CPU_GEMM_CALLS) ;
    P2 (" time %12.4e\n", Common->CHOLMOD_CPU_GEMM_TIME) ;
    P2 ("      GPU calls %12.0f", (double) Common->CHOLMOD_GPU_GEMM_CALLS) ;
    P2 (" time %12.4e\n", Common->CHOLMOD_GPU_GEMM_TIME) ;
    P2 ("POTRF CPU calls %12.0f", (double) Common->CHOLMOD_CPU_POTRF_CALLS) ;
    P2 (" time %12.4e\n", Common->CHOLMOD_CPU_POTRF_TIME) ;
    P2 ("      GPU calls %12.0f", (double) Common->CHOLMOD_GPU_POTRF_CALLS) ;
    P2 (" time %12.4e\n", Common->CHOLMOD_GPU_POTRF_TIME) ;
    P2 ("TRSM  CPU calls %12.0f", (double) Common->CHOLMOD_CPU_TRSM_CALLS) ;
    P2 (" time %12.4e\n", Common->CHOLMOD_CPU_TRSM_TIME) ;
    P2 ("      GPU calls %12.0f", (double) Common->CHOLMOD_GPU_TRSM_CALLS) ;
    P2 (" time %12.4e\n", Common->CHOLMOD_GPU_TRSM_TIME) ;

    cpu_time = Common->CHOLMOD_CPU_SYRK_TIME + Common->CHOLMOD_CPU_TRSM_TIME +
               Common->CHOLMOD_CPU_GEMM_TIME + Common->CHOLMOD_CPU_POTRF_TIME ;

    gpu_time = Common->CHOLMOD_GPU_SYRK_TIME + Common->CHOLMOD_GPU_TRSM_TIME +
               Common->CHOLMOD_GPU_GEMM_TIME + Common->CHOLMOD_GPU_POTRF_TIME ;

    P2 ("time in the BLAS: CPU %12.4e", cpu_time) ;
    P2 (" GPU %12.4e", gpu_time) ;
    P2 (" total: %12.4e\n", cpu_time + gpu_time) ;

    P2 ("assembly time %12.4e", Common->CHOLMOD_ASSEMBLE_TIME) ;
    P2 ("  %12.4e\n", Common->CHOLMOD_ASSEMBLE_TIME2) ;
    return (TRUE) ;
}


/* ========================================================================== */
/* === cholmod_check_sparse ================================================= */
/* ========================================================================== */

/* Ensure that a sparse matrix in column-oriented form is valid, and optionally
 * print it.  Returns the number of entries on the diagonal or -1 if error.
 *
 * workspace: Iwork (nrow)
 */

static SuiteSparse_long check_sparse
(
    Int *Wi,
    Int print,
    const char *name,
    cholmod_sparse *A,
    SuiteSparse_long *nnzdiag,
    cholmod_common *Common
)
{
    double *Ax, *Az ;
    Int *Ap, *Ai, *Anz ;
    Int nrow, ncol, nzmax, sorted, packed, j, p, pend, i, nz, ilast,
	space, init_print, dnz, count, xtype ;
    const char *type = "sparse" ;

    /* ---------------------------------------------------------------------- */
    /* print header information */
    /* ---------------------------------------------------------------------- */

    P4 ("%s", "\n") ;
    P3 ("%s", "CHOLMOD sparse:  ") ;
    if (name != NULL)
    {
	P3 ("%s: ", name) ;
    }

    if (A == NULL)
    {
	ERR ("null") ;
    }

    nrow = A->nrow ;
    ncol = A->ncol ;
    nzmax = A->nzmax ;
    sorted = A->sorted ;
    packed = A->packed ;
    xtype = A->xtype ;
    Ap = A->p ;
    Ai = A->i ;
    Ax = A->x ;
    Az = A->z ;
    Anz = A->nz ;
    nz = CHOLMOD(nnz) (A, Common) ;

    P3 (" "ID"", nrow) ;
    P3 ("-by-"ID", ", ncol) ;
    P3 ("nz "ID",", nz) ;
    if (A->stype > 0)
    {
	P3 ("%s", " upper.") ;
    }
    else if (A->stype < 0)
    {
	P3 ("%s", " lower.") ;
    }
    else
    {
	P3 ("%s", " up/lo.") ;
    }

    P4 ("\n  nzmax "ID", ", nzmax) ;
    if (nz > nzmax)
    {
	ERR ("nzmax too small") ;
    }
    if (!sorted)
    {
	P4 ("%s", "un") ;
    }
    P4 ("%s", "sorted, ") ;
    if (!packed)
    {
	P4 ("%s", "un") ;
    }
    P4 ("%s", "packed, ") ;

    switch (A->itype)
    {
	case CHOLMOD_INT:     P4 ("%s", "\n  scalar types: int, ") ; break ;
	case CHOLMOD_INTLONG: ERR ("mixed int/long type unsupported") ;
	case CHOLMOD_LONG:    P4 ("%s", "\n  scalar types: SuiteSparse_long, ");
        break ;
	default:	      ERR ("unknown itype") ;
    }

    switch (A->xtype)
    {
	case CHOLMOD_PATTERN: P4 ("%s", "pattern") ;	break ;
	case CHOLMOD_REAL:    P4 ("%s", "real") ;	break ;
	case CHOLMOD_COMPLEX: P4 ("%s", "complex") ;	break ;
	case CHOLMOD_ZOMPLEX: P4 ("%s", "zomplex") ;	break ;
	default:	      ERR ("unknown xtype") ;
    }

    switch (A->dtype)
    {
	case CHOLMOD_DOUBLE:  P4 ("%s", ", double\n") ;	       break ;
	case CHOLMOD_SINGLE:  ERR ("float unsupported") ;
	default:	      ERR ("unknown dtype") ;
    }

    if (A->itype != ITYPE || A->dtype != DTYPE)
    {
	ERR ("integer and real type must match routine") ;
    }

    if (A->stype && nrow != ncol)
    {
	ERR ("symmetric but not square") ;
    }

    /* check for existence of Ap, Ai, Anz, Ax, and Az arrays */
    if (Ap == NULL)
    {
	ERR ("p array not present") ;
    }
    if (Ai == NULL)
    {
	ERR ("i array not present") ;
    }
    if (!packed && Anz == NULL)
    {
	ERR ("nz array not present") ;
    }
    if (xtype != CHOLMOD_PATTERN && Ax == NULL)
    {
	ERR ("x array not present") ;
    }
    if (xtype == CHOLMOD_ZOMPLEX && Az == NULL)
    {
	ERR ("z array not present") ;
    }

    /* packed matrices must start at Ap [0] = 0 */
    if (packed && Ap [0] != 0)
    {
	ERR ("p [0] must be zero") ;
    }
    if (packed && (Ap [ncol] < Ap [0] || Ap [ncol] > nzmax))
    {
	ERR ("p [ncol] invalid") ;
    }

    /* ---------------------------------------------------------------------- */
    /* allocate workspace if needed */
    /* ---------------------------------------------------------------------- */

    if (!sorted)
    {
	if (Wi == NULL)
	{
	    CHOLMOD(allocate_work) (0, nrow, 0, Common) ;
	    Wi = Common->Iwork ;	/* size nrow, (i/i/l) */
	}
	if (Common->status < CHOLMOD_OK)
	{
	    return (FALSE) ;	    /* out of memory */
	}
	for (i = 0 ; i < nrow ; i++)
	{
	    Wi [i] = EMPTY ;
	}
    }

    /* ---------------------------------------------------------------------- */
    /* check and print each column */
    /* ---------------------------------------------------------------------- */

    init_print = print ;
    dnz = 0 ;
    ETC_START (count, 8) ;

    for (j = 0 ; j < ncol ; j++)
    {
	ETC (j == ncol-1, count, 4) ;
	p = Ap [j] ;
	if (packed)
	{
	    pend = Ap [j+1] ;
	    nz = pend - p ;
	}
	else
	{
	    /* Note that Anz [j] < 0 is treated as zero */
	    nz = MAX (0, Anz [j]) ;
	    pend = p + nz ;
	}
	/* Note that space can be negative if the matrix is non-monotonic */
	space = Ap [j+1] - p ;
	P4 ("  col "ID":", j) ;
	P4 (" nz "ID"", nz) ;
	P4 (" start "ID"", p) ;
	P4 (" end "ID"", pend) ;
	if (!packed)
	{
	    P4 (" space "ID"", space) ;
	}
	P4 ("%s", ":\n") ;
	if (p < 0 || pend > nzmax)
	{
	    ERR ("pointer invalid") ;
	}
	if (nz < 0 || nz > nrow)
	{
	    ERR ("nz invalid") ;
	}
	ilast = EMPTY ;

	for ( ; p < pend ; p++)
	{
	    ETC (j == ncol-1 && p >= pend-4, count, -1) ;
	    i = Ai [p] ;
	    P4 ("  "I8":", i) ;

	    print_value (print, xtype, Ax, Az, p, Common) ;

	    if (i == j)
	    {
		dnz++ ;
	    }
	    if (i < 0 || i >= nrow)
	    {
		ERR ("row index out of range") ;
	    }
	    if (sorted && i <= ilast)
	    {
		ERR ("row indices out of order") ;
	    }
	    if (!sorted && Wi [i] == j)
	    {
		ERR ("duplicate row index") ;
	    }
	    P4 ("%s", "\n") ;
	    ilast = i ;
	    if (!sorted)
	    {
		Wi [i] = j ;
	    }
	}
    }

    /* matrix is valid */
    P4 ("  nnz on diagonal: "ID"\n", dnz) ;
    P3 ("%s", "  OK\n") ;
    P4 ("%s", "\n") ;
    *nnzdiag = dnz ;
    return (TRUE) ;
}


int CHOLMOD(check_sparse)
(
    /* ---- input ---- */
    cholmod_sparse *A,	/* sparse matrix to check */
    /* --------------- */
    cholmod_common *Common
)
{
    SuiteSparse_long nnzdiag ;
    RETURN_IF_NULL_COMMON (FALSE) ;
    Common->status = CHOLMOD_OK ;
    return (check_sparse (NULL, 0, NULL, A, &nnzdiag, Common)) ;
}


int CHOLMOD(print_sparse)
(
    /* ---- input ---- */
    cholmod_sparse *A,	/* sparse matrix to print */
    const char *name,	/* printed name of sparse matrix */
    /* --------------- */
    cholmod_common *Common
)
{
    SuiteSparse_long nnzdiag ;
    RETURN_IF_NULL_COMMON (FALSE) ;
    Common->status = CHOLMOD_OK ;
    return (check_sparse (NULL, Common->print, name, A, &nnzdiag, Common)) ;
}


/* ========================================================================== */
/* === cholmod_check_dense ================================================== */
/* ========================================================================== */

/* Ensure a dense matrix is valid, and optionally print it. */

static int check_dense
(
    Int print,
    const char *name,
    cholmod_dense *X,
    cholmod_common *Common
)
{
    double *Xx, *Xz ;
    Int i, j, d, nrow, ncol, nzmax, nz, init_print, count, xtype ;
    const char *type = "dense" ;

    /* ---------------------------------------------------------------------- */
    /* print header information */
    /* ---------------------------------------------------------------------- */

    P4 ("%s", "\n") ;
    P3 ("%s", "CHOLMOD dense:   ") ;
    if (name != NULL)
    {
	P3 ("%s: ", name) ;
    }

    if (X == NULL)
    {
	ERR ("null") ;
    }

    nrow = X->nrow ;
    ncol = X->ncol ;
    nzmax = X->nzmax ;
    d = X->d ;
    Xx = X->x ;
    Xz = X->z ;
    xtype = X->xtype ;

    P3 (" "ID"", nrow) ;
    P3 ("-by-"ID", ", ncol) ;
    P4 ("\n  leading dimension "ID", ", d) ;
    P4 ("nzmax "ID", ", nzmax) ;
    if (d * ncol > nzmax)
    {
	ERR ("nzmax too small") ;
    }
    if (d < nrow)
    {
	ERR ("leading dimension must be >= # of rows") ;
    }
    if (Xx == NULL)
    {
	ERR ("null") ;
    }

    switch (X->xtype)
    {
	case CHOLMOD_PATTERN: ERR ("pattern unsupported") ;  break ;
	case CHOLMOD_REAL:    P4 ("%s", "real") ;	break ;
	case CHOLMOD_COMPLEX: P4 ("%s", "complex") ;	break ;
	case CHOLMOD_ZOMPLEX: P4 ("%s", "zomplex") ;	break ;
	default:	      ERR ("unknown xtype") ;
    }

    switch (X->dtype)
    {
	case CHOLMOD_DOUBLE:  P4 ("%s", ", double\n") ;	       break ;
	case CHOLMOD_SINGLE:  ERR ("single unsupported") ;
	default:	      ERR ("unknown dtype") ;
    }

    /* ---------------------------------------------------------------------- */
    /* check and print each entry */
    /* ---------------------------------------------------------------------- */

    if (print >= 4)
    {
	init_print = print ;
	ETC_START (count, 9) ;
	nz = nrow * ncol ;
	for (j = 0 ; j < ncol ; j++)
	{
	    ETC (j == ncol-1, count, 5) ;
	    P4 ("  col "ID":\n", j) ;
	    for (i = 0 ; i < nrow ; i++)
	    {
		ETC (j == ncol-1 && i >= nrow-4, count, -1) ;
		P4 ("  "I8":", i) ;

		print_value (print, xtype, Xx, Xz, i+j*d, Common) ;

		P4 ("%s", "\n") ;
	    }
	}
    }

    /* dense  is valid */
    P3 ("%s", "  OK\n") ;
    P4 ("%s", "\n") ;
    return (TRUE) ;
}


int CHOLMOD(check_dense)
(
    /* ---- input ---- */
    cholmod_dense *X,	/* dense matrix to check */
    /* --------------- */
    cholmod_common *Common
)
{
    RETURN_IF_NULL_COMMON (FALSE) ;
    Common->status = CHOLMOD_OK ;
    return (check_dense (0, NULL, X, Common)) ;
}


int CHOLMOD(print_dense)
(
    /* ---- input ---- */
    cholmod_dense *X,	/* dense matrix to print */
    const char *name,	/* printed name of dense matrix */
    /* --------------- */
    cholmod_common *Common
)
{
    RETURN_IF_NULL_COMMON (FALSE) ;
    Common->status = CHOLMOD_OK ;
    return (check_dense (Common->print, name, X, Common)) ;
}


/* ========================================================================== */
/* === cholmod_check_subset ================================================= */
/* ========================================================================== */

/* Ensure S (0:len-1) is a subset of 0:n-1.  Duplicates are allowed.  S may be
 * NULL.  A negative len denotes the set 0:n-1.
 *
 * To check the rset and cset for A(rset,cset), where nc and nr are the length
 * of cset and rset respectively:
 *
 *	cholmod_check_subset (cset, nc, A->ncol, Common) ;
 *	cholmod_check_subset (rset, nr, A->nrow, Common) ;
 *
 * workspace: none
 */

static int check_subset
(
    Int *S,
    SuiteSparse_long len,
    size_t n,
    Int print,
    const char *name,
    cholmod_common *Common
)
{
    Int i, k, init_print, count ;
    const char *type = "subset" ;

    init_print = print ;

    if (S == NULL)
    {
	/* zero len denotes S = [ ], negative len denotes S = 0:n-1 */
	len = (len < 0) ? (-1) : 0 ;
    }

    P4 ("%s", "\n") ;
    P3 ("%s", "CHOLMOD subset:  ") ;
    if (name != NULL)
    {
	P3 ("%s: ", name) ;
    }

    P3 (" len: %ld ", len) ;
    if (len < 0)
    {
	P3 ("%s", "(denotes 0:n-1) ") ;
    }
    P3 ("n: "ID"", (Int) n) ;
    P4 ("%s", "\n") ;

    if (len <= 0 || S == NULL)
    {
	P3 ("%s", "  OK\n") ;
	P4 ("%s", "\n") ;
	return (TRUE) ;
    }

    if (print >= 4)
    {
	ETC_START (count, 8) ;
	for (k = 0 ; k < ((Int) len) ; k++)
	{
	    ETC (k == ((Int) len) - 4, count, -1) ;
	    i = S [k] ;
	    P4 ("  "I8":", k) ;
	    P4 (" "ID"\n", i) ;
	    if (i < 0 || i >= ((Int) n))
	    {
		ERR ("entry out range") ;
	    }
	}
    }
    else
    {
	for (k = 0 ; k < ((Int) len) ; k++)
	{
	    i = S [k] ;
	    if (i < 0 || i >= ((Int) n))
	    {
		ERR ("entry out range") ;
	    }
	}
    }
    P3 ("%s", "  OK\n") ;
    P4 ("%s", "\n") ;
    return (TRUE) ;
}


int CHOLMOD(check_subset)
(
    /* ---- input ---- */
    Int *Set,		/* Set [0:len-1] is a subset of 0:n-1.  Duplicates OK */
    SuiteSparse_long len, /* size of Set (an integer array), or < 0 if 0:n-1 */
    size_t n,		/* 0:n-1 is valid range */
    /* --------------- */
    cholmod_common *Common
)
{
    RETURN_IF_NULL_COMMON (FALSE) ;
    Common->status = CHOLMOD_OK ;
    return (check_subset (Set, len, n, 0, NULL, Common)) ;
}


int CHOLMOD(print_subset)
(
    /* ---- input ---- */
    Int *Set,		/* Set [0:len-1] is a subset of 0:n-1.  Duplicates OK */
    SuiteSparse_long len, /* size of Set (an integer array), or < 0 if 0:n-1 */
    size_t n,		/* 0:n-1 is valid range */
    const char *name,	/* printed name of Set */
    /* --------------- */
    cholmod_common *Common
)
{
    RETURN_IF_NULL_COMMON (FALSE) ;
    Common->status = CHOLMOD_OK ;
    return (check_subset (Set, len, n, Common->print, name, Common)) ;
}


/* ========================================================================== */
/* === cholmod_check_perm =================================================== */
/* ========================================================================== */

/* Ensure that Perm [0..len-1] is a permutation of a subset of 0:n-1.  Perm
 * may be NULL, which is interpreted as the identity permutation.  There can
 * be no duplicate entries (len must be <= n).
 *
 * If n <= Common->nrow, then this routine takes O(len) time and does not
 * allocate any memory, by using Common->Flag.  Otherwise, it takes O(n) time
 * and ensures that Common->Iwork is at least n*sizeof(Int) in size.
 *
 * To check the fset:	    cholmod_check_perm (fset, fsize, ncol, Common) ;
 * To check a permutation:  cholmod_check_perm (Perm, n, n, Common) ;
 *
 * workspace:  Flag (n) if n <= Common->nrow, Iwork (n) otherwise.
 */

static int check_perm
(
    Int *Wi,
    Int print,
    const char *name,
    Int *Perm,
    size_t len,
    size_t n,
    cholmod_common *Common
)
{
    Int *Flag ;
    Int i, k, mark, init_print, count ;
    const char *type = "perm" ;

    /* ---------------------------------------------------------------------- */
    /* checks that take O(1) time */
    /* ---------------------------------------------------------------------- */

    if (Perm == NULL || n == 0)
    {
	/* Perm is valid implicit identity, or empty */
	return (TRUE) ;
    }

    /* ---------------------------------------------------------------------- */
    /* checks that take O(n) time or require memory allocation */
    /* ---------------------------------------------------------------------- */

    init_print = print ;
    ETC_START (count, 8) ;

    if (Wi == NULL && n <= Common->nrow)
    {
	/* use the Common->Flag array if it's big enough */
	mark = CHOLMOD(clear_flag) (Common) ;
	Flag = Common->Flag ;
	ASSERT (CHOLMOD(dump_work) (TRUE, FALSE, 0, Common)) ;
	if (print >= 4)
	{
	    for (k = 0 ; k < ((Int) len) ; k++)
	    {
		ETC (k >= ((Int) len) - 4, count, -1) ;
		i = Perm [k] ;
		P4 ("  "I8":", k) ;
		P4 (""ID"\n", i) ;
		if (i < 0 || i >= ((Int) n) || Flag [i] == mark)
		{
		    CHOLMOD(clear_flag) (Common) ;
		    ERR ("invalid permutation") ;
		}
		Flag [i] = mark ;
	    }
	}
	else
	{
	    for (k = 0 ; k < ((Int) len) ; k++)
	    {
		i = Perm [k] ;
		if (i < 0 || i >= ((Int) n) || Flag [i] == mark)
		{
		    CHOLMOD(clear_flag) (Common) ;
		    ERR ("invalid permutation") ;
		}
		Flag [i] = mark ;
	    }
	}
	CHOLMOD(clear_flag) (Common) ;
	ASSERT (CHOLMOD(dump_work) (TRUE, FALSE, 0, Common)) ;
    }
    else
    {
	if (Wi == NULL)
	{
	    /* use Common->Iwork instead, but initialize it first */
	    CHOLMOD(allocate_work) (0, n, 0, Common) ;
	    Wi = Common->Iwork ;		    /* size n, (i/i/i) is OK */
	}
	if (Common->status < CHOLMOD_OK)
	{
	    return (FALSE) ;	    /* out of memory */
	}
	for (i = 0 ; i < ((Int) n) ; i++)
	{
	    Wi [i] = FALSE ;
	}
	if (print >= 4)
	{
	    for (k = 0 ; k < ((Int) len) ; k++)
	    {
		ETC (k >= ((Int) len) - 4, count, -1) ;
		i = Perm [k] ;
		P4 ("  "I8":", k) ;
		P4 (""ID"\n", i) ;
		if (i < 0 || i >= ((Int) n) || Wi [i])
		{
		    ERR ("invalid permutation") ;
		}
		Wi [i] = TRUE ;
	    }
	}
	else
	{
	    for (k = 0 ; k < ((Int) len) ; k++)
	    {
		i = Perm [k] ;
		if (i < 0 || i >= ((Int) n) || Wi [i])
		{
		    ERR ("invalid permutation") ;
		}
		Wi [i] = TRUE ;
	    }
	}
    }

    /* perm is valid */
    return (TRUE) ;
}


int CHOLMOD(check_perm)
(
    /* ---- input ---- */
    Int *Perm,		/* Perm [0:len-1] is a permutation of subset of 0:n-1 */
    size_t len,		/* size of Perm (an integer array) */
    size_t n,		/* 0:n-1 is valid range */
    /* --------------- */
    cholmod_common *Common
)
{
    RETURN_IF_NULL_COMMON (FALSE) ;
    Common->status = CHOLMOD_OK ;
    return (check_perm (NULL, 0, NULL, Perm, len, n, Common)) ;
}


int CHOLMOD(print_perm)
(
    /* ---- input ---- */
    Int *Perm,		/* Perm [0:len-1] is a permutation of subset of 0:n-1 */
    size_t len,		/* size of Perm (an integer array) */
    size_t n,		/* 0:n-1 is valid range */
    const char *name,	/* printed name of Perm */
    /* --------------- */
    cholmod_common *Common
)
{
    Int ok, print ;
    RETURN_IF_NULL_COMMON (FALSE) ;
    Common->status = CHOLMOD_OK ;
    print = Common->print ;
    P4 ("%s", "\n") ;
    P3 ("%s", "CHOLMOD perm:    ") ;
    if (name != NULL)
    {
	P3 ("%s: ", name) ;
    }
    P3 (" len: "ID"", (Int) len) ;
    P3 (" n: "ID"", (Int) n) ;
    P4 ("%s", "\n") ;
    ok = check_perm (NULL, print, name, Perm, len, n, Common) ;
    if (ok)
    {
	P3 ("%s", "  OK\n") ;
	P4 ("%s", "\n") ;
    }
    return (ok) ;
}


/* ========================================================================== */
/* === cholmod_check_parent ================================================= */
/* ========================================================================== */

/* Ensure that Parent is a valid elimination tree of nodes 0 to n-1.
 * If j is a root of the tree then Parent [j] is EMPTY (-1).
 *
 * NOTE: this check will fail if applied to the component tree (CParent) in
 * cholmod_nested_dissection, unless it has been postordered and renumbered.
 *
 * workspace: none
 */

static int check_parent
(
    Int *Parent,
    size_t n,
    Int print,
    const char *name,
    cholmod_common *Common
)
{
    Int j, p, init_print, count ;
    const char *type = "parent" ;

    init_print = print ;

    P4 ("%s", "\n") ;
    P3 ("%s", "CHOLMOD parent:  ") ;
    if (name != NULL)
    {
	P3 ("%s: ", name) ;
    }

    P3 (" n: "ID"", (Int) n) ;
    P4 ("%s", "\n") ;

    if (Parent == NULL)
    {
	ERR ("null") ;
    }

    /* ---------------------------------------------------------------------- */
    /* checks that take O(n) time */
    /* ---------------------------------------------------------------------- */

    ETC_START (count, 8) ;
    for (j = 0 ; j < ((Int) n) ; j++)
    {
	ETC (j == ((Int) n) - 4, count, -1) ;
	p = Parent [j] ;
	P4 ("  "I8":", j) ;
	P4 (" "ID"\n", p) ;
	if (!(p == EMPTY || p > j))
	{
	    ERR ("invalid") ;
	}
    }
    P3 ("%s", "  OK\n") ;
    P4 ("%s", "\n") ;
    return (TRUE) ;
}


int CHOLMOD(check_parent)
(
    /* ---- input ---- */
    Int *Parent,	/* Parent [0:n-1] is an elimination tree */
    size_t n,		/* size of Parent */
    /* --------------- */
    cholmod_common *Common
)
{
    RETURN_IF_NULL_COMMON (FALSE) ;
    Common->status = CHOLMOD_OK ;
    return (check_parent (Parent, n, 0, NULL, Common)) ;
}


int CHOLMOD(print_parent)
(
    /* ---- input ---- */
    Int *Parent,	/* Parent [0:n-1] is an elimination tree */
    size_t n,		/* size of Parent */
    const char *name,	/* printed name of Parent */
    /* --------------- */
    cholmod_common *Common
)
{
    RETURN_IF_NULL_COMMON (FALSE) ;
    Common->status = CHOLMOD_OK ;
    return (check_parent (Parent, n, Common->print, name, Common)) ;
}


/* ========================================================================== */
/* === cholmod_check_factor ================================================= */
/* ========================================================================== */

static int check_factor
(
    Int *Wi,
    Int print,
    const char *name,
    cholmod_factor *L,
    cholmod_common *Common
)
{
    double *Lx, *Lz ;
    Int *Lp, *Li, *Lnz, *Lnext, *Lprev, *Perm, *ColCount, *Lpi, *Lpx, *Super,
	*Ls ;
    Int n, nzmax, j, p, pend, i, nz, ordering, space, is_monotonic, minor,
	count, precise, init_print, ilast, lnz, head, tail, jprev, plast,
	jnext, examine_super, nsuper, s, k1, k2, psi, psend, psx, nsrow, nscol,
	ps2, psxend, ssize, xsize, maxcsize, maxesize, nsrow2, jj, ii, xtype ;
    Int for_cholesky ;
    const char *type = "factor" ;

    /* ---------------------------------------------------------------------- */
    /* print header information */
    /* ---------------------------------------------------------------------- */

    P4 ("%s", "\n") ;
    P3 ("%s", "CHOLMOD factor:  ") ;
    if (name != NULL)
    {
	P3 ("%s: ", name) ;
    }

    if (L == NULL)
    {
	ERR ("null") ;
    }

    n = L->n ;
    minor = L->minor ;
    ordering = L->ordering ;
    xtype = L->xtype ;

    Perm = L->Perm ;
    ColCount = L->ColCount ;
    lnz = 0 ;

    precise = Common->precise ;

    P3 (" "ID"", n) ;
    P3 ("-by-"ID"", n) ;

    if (minor < n)
    {
	P3 (" not positive definite (column "ID")", minor) ;
    }

    switch (L->itype)
    {
	case CHOLMOD_INT:     P4 ("%s", "\n  scalar types: int, ") ; break ;
	case CHOLMOD_INTLONG: ERR ("mixed int/long type unsupported") ;
	case CHOLMOD_LONG:    P4 ("%s", "\n  scalar types: SuiteSparse_long, ");
        break ;
	default:	      ERR ("unknown itype") ;
    }

    switch (L->xtype)
    {
	case CHOLMOD_PATTERN: P4 ("%s", "pattern") ;	break ;
	case CHOLMOD_REAL:    P4 ("%s", "real") ;	break ;
	case CHOLMOD_COMPLEX: P4 ("%s", "complex") ;	break ;
	case CHOLMOD_ZOMPLEX: P4 ("%s", "zomplex") ;	break ;
	default:	      ERR ("unknown xtype") ;
    }

    switch (L->dtype)
    {
	case CHOLMOD_DOUBLE:  P4 ("%s", ", double\n") ;	       break ;
	case CHOLMOD_SINGLE:  ERR ("single unsupported") ;
	default:	      ERR ("unknown dtype") ;
    }

    if (L->itype != ITYPE || L->dtype != DTYPE)
    {
	ERR ("integer and real type must match routine") ;
    }

    if (L->is_super)
    {
	P3 ("%s", "  supernodal") ;
    }
    else
    {
	P3 ("%s", "  simplicial") ;
    }

    if (L->is_ll)
    {
	P3 ("%s", ", LL'.") ;
    }
    else
    {
	P3 ("%s", ", LDL'.") ;
    }

    P4 ("%s", "\n  ordering method used: ") ;
    switch (L->ordering)
    {
	case CHOLMOD_POSTORDERED:P4 ("%s", "natural (postordered)") ;	 break ;
	case CHOLMOD_NATURAL:	P4 ("%s", "natural") ;			 break ;
	case CHOLMOD_GIVEN:	P4 ("%s", "user-provided") ;		 break ;
	case CHOLMOD_AMD:	P4 ("%s", "AMD") ;			 break ;
	case CHOLMOD_COLAMD:	P4 ("%s", "AMD for A, COLAMD for A*A'") ;break ;
#ifndef NPARTITION
	case CHOLMOD_METIS:	P4 ("%s", "METIS NodeND") ;		 break ;
	case CHOLMOD_NESDIS:	P4 ("%s", "CHOLMOD nested dissection") ; break ;
#endif
	default:		ERR ("unknown ordering") ;
    }

    P4 ("%s", "\n") ;

    init_print = print ;

    if (L->is_super && L->xtype == CHOLMOD_ZOMPLEX)
    {
	ERR ("Supernodal zomplex L not supported") ;
    }

    /* ---------------------------------------------------------------------- */
    /* check L->Perm */
    /* ---------------------------------------------------------------------- */

    if (!check_perm (Wi, print, name, Perm, n, n, Common))
    {
	return (FALSE) ;
    }

    /* ---------------------------------------------------------------------- */
    /* check L->ColCount */
    /* ---------------------------------------------------------------------- */

    if (ColCount == NULL)
    {
	ERR ("ColCount vector invalid") ;
    }

    ETC_START (count, 8) ;
    for (j = 0 ; j < n ; j++)
    {
	ETC (j >= n-4, count, -1) ;
	P4 ("  col: "ID" ", j) ;
	nz = ColCount [j] ;
	P4 ("colcount: "ID"\n", nz) ;
	if (nz < 0 || nz > n-j)
	{
	    ERR ("ColCount out of range") ;
	}
    }

    /* ---------------------------------------------------------------------- */
    /* check factor */
    /* ---------------------------------------------------------------------- */

    if (L->xtype == CHOLMOD_PATTERN && !(L->is_super))
    {

	/* ------------------------------------------------------------------ */
	/* check simplicial symbolic factor */
	/* ------------------------------------------------------------------ */

	/* nothing else to do */ ;

    }
    else if (L->xtype != CHOLMOD_PATTERN && !(L->is_super))
    {

	/* ------------------------------------------------------------------ */
	/* check simplicial numerical factor */
	/* ------------------------------------------------------------------ */

	P4 ("monotonic: %d\n", L->is_monotonic) ;
	nzmax = L->nzmax ;
	P3 (" nzmax "ID".", nzmax) ;
	P4 ("%s", "\n") ;
	Lp = L->p ;
	Li = L->i ;
	Lx = L->x ;
	Lz = L->z ;
	Lnz = L->nz ;
	Lnext = L->next ;
	Lprev = L->prev ;

	/* check for existence of Lp, Li, Lnz, Lnext, Lprev, and Lx arrays */
	if (Lp == NULL)
	{
	    ERR ("p array not present") ;
	}
	if (Li == NULL)
	{
	    ERR ("i array not present") ;
	}
	if (Lnz == NULL)
	{
	    ERR ("nz array not present") ;
	}
	if (Lx == NULL)
	{
	    ERR ("x array not present") ;
	}
	if (xtype == CHOLMOD_ZOMPLEX && Lz == NULL)
	{
	    ERR ("z array not present") ;
	}
	if (Lnext == NULL)
	{
	    ERR ("next array not present") ;
	}
	if (Lprev == NULL)
	{
	    ERR ("prev array not present") ;
	}

	ETC_START (count, 8) ;

	/* check each column of L */
	plast = 0 ;
	is_monotonic = TRUE ;
	for (j = 0 ; j < n ; j++)
	{
	    ETC (j >= n-3, count, -1) ;
	    p = Lp [j] ;
	    nz = Lnz [j] ;
	    pend = p + nz ;
	    lnz += nz ;

	    P4 ("  col "ID":", j) ;
	    P4 (" nz "ID"", nz) ;
	    P4 (" start "ID"", p) ;
	    P4 (" end "ID"", pend) ;

	    if (Lnext [j] < 0 || Lnext [j] > n)
	    {
		ERR ("invalid link list")  ;
	    }
	    space = Lp [Lnext [j]] - p ;

	    P4 (" space "ID"", space) ;
	    P4 (" free "ID":\n", space - nz) ;

	    if (p < 0 || pend > nzmax || space < 1)
	    {
		ERR ("pointer invalid") ;
	    }
	    if (nz < 1 || nz > (n-j) || nz > space)
	    {
		ERR ("nz invalid") ;
	    }
	    ilast = j-1 ;

	    if (p < plast)
	    {
		is_monotonic = FALSE ;
	    }
	    plast = p ;

	    i = Li [p] ;
	    P4 ("  "I8":", i) ;
	    if (i != j)
	    {
		ERR ("diagonal missing") ;
	    }

	    print_value (print, xtype, Lx, Lz, p, Common) ;

	    P4 ("%s", "\n") ;
	    ilast = j ;
	    for (p++ ; p < pend ; p++)
	    {
		ETC_DISABLE (count) ;
		i = Li [p] ;
		P4 ("  "I8":", i) ;
		if (i < j || i >= n)
		{
		    ERR ("row index out of range") ;
		}
		if (i <= ilast)
		{
		    ERR ("row indices out of order") ;
		}

		print_value (print, xtype, Lx, Lz, p, Common) ;

		P4 ("%s", "\n") ;
		ilast = i ;
	    }
	}

	if (L->is_monotonic && !is_monotonic)
	{
	    ERR ("columns not monotonic") ;
	}

	/* check the link list */
	head = n+1 ;
	tail = n ;
	j = head ;
	jprev = EMPTY ;
	count = 0 ;
	for ( ; ; )
	{
	    if (j < 0 || j > n+1 || count > n+2)
	    {
		ERR ("invalid link list") ;
	    }
	    jnext = Lnext [j] ;
	    if (j >= 0 && j < n)
	    {
		if (jprev != Lprev [j])
		{
		    ERR ("invalid link list") ;
		}
	    }
	    count++ ;
	    if (j == tail)
	    {
		break ;
	    }
	    jprev = j ;
	    j = jnext ;
	}
	if (Lnext [tail] != EMPTY || count != n+2)
	{
	    ERR ("invalid link list") ;
	}

    }
    else
    {

	/* ------------------------------------------------------------------ */
	/* check supernodal numeric or symbolic factor */
	/* ------------------------------------------------------------------ */

	nsuper = L->nsuper ;
	ssize = L->ssize ;
	xsize = L->xsize ;
	maxcsize = L->maxcsize ;
	maxesize = L->maxesize ;
	Ls = L->s ;
	Lpi = L->pi ;
	Lpx = L->px ;
	Super = L->super ;
	Lx = L->x ;
	ETC_START (count, 8) ;

	P4 ("  ssize "ID" ", ssize) ;
	P4 ("xsize "ID" ", xsize) ;
	P4 ("maxcsize "ID" ", maxcsize) ;
	P4 ("maxesize "ID"\n", maxesize) ;

	if (Ls == NULL)
	{
	    ERR ("invalid: L->s missing") ;
	}
	if (Lpi == NULL)
	{
	    ERR ("invalid: L->pi missing") ;
	}
	if (Lpx == NULL)
	{
	    ERR ("invalid: L->px missing") ;
	}
	if (Super == NULL)
	{
	    ERR ("invalid: L->super missing") ;
	}

	if (L->xtype != CHOLMOD_PATTERN)
	{
	    /* numerical supernodal factor */
	    if (Lx == NULL)
	    {
		ERR ("invalid: L->x missing") ;
	    }
	    if (Ls [0] == EMPTY)
	    {
		ERR ("invalid: L->s not defined") ;
	    }
	    examine_super = TRUE ;
	}
	else
	{
	    /* symbolic supernodal factor, but only if it has been computed */
	    examine_super = (Ls [0] != EMPTY) ;
	}

	if (examine_super)
	{
	    if (Lpi [0] != 0 || MAX (1, Lpi [nsuper]) != ssize)
	    {
		PRINT0 (("Lpi [0] "ID", Lpi [nsuper = "ID"] = "ID"\n",
			Lpi [0], nsuper, Lpi [nsuper])) ;
		ERR ("invalid: L->pi invalid") ;
	    }

            for_cholesky = (Lpx [0] != 123456) ;
	    if (for_cholesky && (Lpx [0] != 0 || MAX (1, Lpx[nsuper]) != xsize))
	    {
		ERR ("invalid: L->px invalid") ;
	    }

	    /* check and print each supernode */
	    for (s = 0 ; s < nsuper ; s++)
	    {
		k1 = Super [s] ;
		k2 = Super [s+1] ;
		psi = Lpi [s] ;
		psend = Lpi [s+1] ;
		nsrow = psend - psi ;
		nscol = k2 - k1 ;
		nsrow2 = nsrow - nscol ;
		ps2 = psi + nscol ;

                if (for_cholesky)
                {
                    psx = Lpx [s] ;
                    psxend = Lpx [s+1] ;
                }

		ETC (s == nsuper-1, count, 4) ;

		P4 ("  supernode "ID", ", s) ;
		P4 ("col "ID" ", k1) ;
		P4 ("to "ID". ", k2-1) ;
		P4 ("nz in first col: "ID".\n", nsrow) ;

                if (for_cholesky)
                {
                    P4 ("  values start "ID", ", psx) ;
                    P4 ("end "ID"\n", psxend) ;
                }

		if (k1 > k2 || k1 < 0 || k2 > n || nsrow < nscol || nsrow2 < 0
                    || (for_cholesky && psxend - psx != nsrow * nscol))
		{
		    ERR ("invalid supernode") ;
		}

		lnz += nscol * nsrow - (nscol*nscol - nscol)/2 ;

		if (L->xtype != CHOLMOD_PATTERN)
		{
		    /* print each column of the supernode */
		    for (jj = 0 ; jj < nscol ; jj++)
		    {
			ETC_ENABLE (s == nsuper-1 && jj >= nscol-3, count, -1) ;
			j = k1 + jj ;
			P4 ("  col "ID"\n", j) ;
			ilast = j ;
			i = Ls [psi + jj] ;
			P4 ("  "I8":", i) ;
			if (i != j)
			{
			    ERR ("row index invalid") ;
			}

			/* PRINTVALUE (Lx [psx + jj + jj*nsrow]) ; */
			print_value (print, xtype, Lx, NULL,
				psx + jj + jj*nsrow, Common) ;

			P4 ("%s", "\n") ;
			for (ii = jj + 1 ; ii < nsrow ; ii++)
			{
			    ETC_DISABLE (count) ;
			    i = Ls [psi + ii] ;
			    P4 ("  "I8":", i) ;
			    if (i <= ilast || i > n)
			    {
				ERR ("row index out of range") ;
			    }

			    /* PRINTVALUE (Lx [psx + ii + jj*nsrow]) ; */
			    print_value (print, xtype, Lx, NULL,
				    psx + ii + jj*nsrow, Common) ;

			    P4 ("%s", "\n") ;
			    ilast = i ;
			}
		    }
		}
		else
		{
		    /* just print the leading column of the supernode */
		    P4 ("  col "ID"\n", k1) ;
		    for (jj = 0 ; jj < nscol ; jj++)
		    {
			ETC (s == nsuper-1 && jj >= nscol-3, count, -1) ;
			j = k1 + jj ;
			i = Ls [psi + jj] ;
			P4 ("  "I8"", i) ;
			if (i != j)
			{
			    ERR ("row index invalid") ;
			}
			P4 ("%s", "\n") ;
		    }
		    ilast = j ;
		    for (ii = nscol ; ii < nsrow ; ii++)
		    {
			ETC_DISABLE (count) ;
			i = Ls [psi + ii] ;
			P4 ("  "I8"", i) ;
			if (i <= ilast || i > n)
			{
			    ERR ("row index out of range") ;
			}
			P4 ("%s", "\n") ;
			ilast = i ;
		    }
		}
	    }
	}
    }

    /* factor is valid */
    P3 ("  nz "ID"", lnz) ;
    P3 ("%s", "  OK\n") ;
    P4 ("%s", "\n") ;
    return (TRUE) ;
}


int CHOLMOD(check_factor)
(
    /* ---- input ---- */
    cholmod_factor *L,	/* factor to check */
    /* --------------- */
    cholmod_common *Common
)
{
    RETURN_IF_NULL_COMMON (FALSE) ;
    Common->status = CHOLMOD_OK ;
    return (check_factor (NULL, 0, NULL, L, Common)) ;
}


int CHOLMOD(print_factor)
(
    /* ---- input ---- */
    cholmod_factor *L,	/* factor to print */
    const char *name,	/* printed name of factor */
    /* --------------- */
    cholmod_common *Common
)
{
    RETURN_IF_NULL_COMMON (FALSE) ;
    Common->status = CHOLMOD_OK ;
    return (check_factor (NULL, Common->print, name, L, Common)) ;
}


/* ========================================================================== */
/* === cholmod_check_triplet ================================================ */
/* ========================================================================== */

/* Ensure a triplet matrix is valid, and optionally print it. */

static int check_triplet
(
    Int print,
    const char *name,
    cholmod_triplet *T,
    cholmod_common *Common
)
{
    double *Tx, *Tz ;
    Int *Ti, *Tj ;
    Int i, j, p, nrow, ncol, nzmax, nz, xtype, init_print, count ;
    const char *type = "triplet" ;

    /* ---------------------------------------------------------------------- */
    /* print header information */
    /* ---------------------------------------------------------------------- */

    P4 ("%s", "\n") ;
    P3 ("%s", "CHOLMOD triplet: ") ;
    if (name != NULL)
    {
	P3 ("%s: ", name) ;
    }

    if (T == NULL)
    {
	ERR ("null") ;
    }

    nrow = T->nrow ;
    ncol = T->ncol ;
    nzmax = T->nzmax ;
    nz = T->nnz ;
    Ti = T->i ;
    Tj = T->j ;
    Tx = T->x ;
    Tz = T->z ;
    xtype = T->xtype ;


    P3 (" "ID"", nrow) ;
    P3 ("-by-"ID", ", ncol) ;
    P3 ("nz "ID",", nz) ;
    if (T->stype > 0)
    {
	P3 ("%s", " upper.") ;
    }
    else if (T->stype < 0)
    {
	P3 ("%s", " lower.") ;
    }
    else
    {
	P3 ("%s", " up/lo.") ;
    }

    P4 ("\n  nzmax "ID", ", nzmax) ;
    if (nz > nzmax)
    {
	ERR ("nzmax too small") ;
    }

    switch (T->itype)
    {
	case CHOLMOD_INT:     P4 ("%s", "\n  scalar types: int, ") ; break ;
	case CHOLMOD_INTLONG: ERR ("mixed int/long type unsupported") ;
	case CHOLMOD_LONG:    P4 ("%s", "\n  scalar types: SuiteSparse_long, ");
        break ;
	default:	      ERR ("unknown itype") ;
    }

    switch (T->xtype)
    {
	case CHOLMOD_PATTERN: P4 ("%s", "pattern") ;	break ;
	case CHOLMOD_REAL:    P4 ("%s", "real") ;	break ;
	case CHOLMOD_COMPLEX: P4 ("%s", "complex") ;	break ;
	case CHOLMOD_ZOMPLEX: P4 ("%s", "zomplex") ;	break ;
	default:	      ERR ("unknown xtype") ;
    }

    switch (T->dtype)
    {
	case CHOLMOD_DOUBLE:  P4 ("%s", ", double\n") ;	       break ;
	case CHOLMOD_SINGLE:  ERR ("single unsupported") ;
	default:	      ERR ("unknown dtype") ;
    }

    if (T->itype != ITYPE || T->dtype != DTYPE)
    {
	ERR ("integer and real type must match routine") ;
    }

    if (T->stype && nrow != ncol)
    {
	ERR ("symmetric but not square") ;
    }

    /* check for existence of Ti, Tj, Tx arrays */
    if (Tj == NULL)
    {
	ERR ("j array not present") ;
    }
    if (Ti == NULL)
    {
	ERR ("i array not present") ;
    }

    if (xtype != CHOLMOD_PATTERN && Tx == NULL)
    {
	ERR ("x array not present") ;
    }
    if (xtype == CHOLMOD_ZOMPLEX && Tz == NULL)
    {
	ERR ("z array not present") ;
    }

    /* ---------------------------------------------------------------------- */
    /* check and print each entry */
    /* ---------------------------------------------------------------------- */

    init_print = print ;
    ETC_START (count, 8) ;

    for (p = 0 ; p < nz ; p++)
    {
	ETC (p >= nz-4, count, -1) ;
	i = Ti [p] ;
	P4 ("  "I8":", p) ;
	P4 (" "I_8"", i) ;
	if (i < 0 || i >= nrow)
	{
	    ERR ("row index out of range") ;
	}
	j = Tj [p] ;
	P4 (" "I_8"", j) ;
	if (j < 0 || j >= ncol)
	{
	    ERR ("column index out of range") ;
	}

	print_value (print, xtype, Tx, Tz, p, Common) ;

	P4 ("%s", "\n") ;
    }

    /* triplet matrix is valid */
    P3 ("%s", "  OK\n") ;
    P4 ("%s", "\n") ;
    return (TRUE) ;
}



int CHOLMOD(check_triplet)
(
    /* ---- input ---- */
    cholmod_triplet *T,	/* triplet matrix to check */
    /* --------------- */
    cholmod_common *Common
)
{
    RETURN_IF_NULL_COMMON (FALSE) ;
    Common->status = CHOLMOD_OK ;
    return (check_triplet (0, NULL, T, Common)) ;
}


int CHOLMOD(print_triplet)
(
    /* ---- input ---- */
    cholmod_triplet *T,	/* triplet matrix to print */
    const char *name,	/* printed name of triplet matrix */
    /* --------------- */
    cholmod_common *Common
)
{
    RETURN_IF_NULL_COMMON (FALSE) ;
    Common->status = CHOLMOD_OK ;
    return (check_triplet (Common->print, name, T, Common)) ;
}



/* ========================================================================== */
/* === CHOLMOD debugging routines =========================================== */
/* ========================================================================== */

#ifndef NDEBUG

/* The global variables present only when debugging enabled. */
int CHOLMOD(dump) = 0 ;
int CHOLMOD(dump_malloc) = -1 ;

/* workspace: no debug routines use workspace in Common */

/* ========================================================================== */
/* === cholmod_dump_init ==================================================== */
/* ========================================================================== */

void CHOLMOD(dump_init) (const char *s, cholmod_common *Common)
{
    int i = 0 ;
    FILE *f ;
    f = fopen ("debug", "r") ;
    CHOLMOD(dump) = 0 ;
    if (f != NULL)
    {
	i = fscanf (f, "%d", &CHOLMOD(dump)) ;
	fclose (f) ;
    }
    PRINT1 (("%s: cholmod_dump_init, D = %d\n", s, CHOLMOD(dump))) ;
}


/* ========================================================================== */
/* === cholmod_dump_sparse ================================================== */
/* ========================================================================== */

/* returns nnz (diag (A)) or EMPTY if error */

SuiteSparse_long CHOLMOD(dump_sparse)
(
    cholmod_sparse *A,
    const char *name,
    cholmod_common *Common
)
{
    Int *Wi ;
    SuiteSparse_long nnzdiag ;
    Int ok ;

    if (CHOLMOD(dump) < -1)
    {
	/* no checks if debug level is -2 or less */
	return (0) ;
    }

    RETURN_IF_NULL_COMMON (FALSE) ;
    RETURN_IF_NULL (A, FALSE) ;
    Wi = malloc (MAX (1, A->nrow) * sizeof (Int)) ;
    ok = check_sparse (Wi, CHOLMOD(dump), name, A, &nnzdiag, Common) ;
    if (Wi != NULL) free (Wi) ;
    return (ok ? nnzdiag : EMPTY) ;
}


/* ========================================================================== */
/* === cholmod_dump_factor ================================================== */
/* ========================================================================== */

int CHOLMOD(dump_factor)
(
    cholmod_factor *L,
    const char *name,
    cholmod_common *Common
)
{
    Int *Wi ;
    int ok ;

    if (CHOLMOD(dump) < -1)
    {
	/* no checks if debug level is -2 or less */
	return (TRUE) ;
    }
    RETURN_IF_NULL_COMMON (FALSE) ;
    RETURN_IF_NULL (L, FALSE) ;
    Wi = malloc (MAX (1, L->n) * sizeof (Int)) ;
    ok = check_factor (Wi, CHOLMOD(dump), name, L, Common) ;
    if (Wi != NULL) free (Wi) ;
    return (ok) ;
}


/* ========================================================================== */
/* === cholmod_dump_perm ==================================================== */
/* ========================================================================== */

int CHOLMOD(dump_perm)
(
    Int *Perm,
    size_t len,
    size_t n,
    const char *name,
    cholmod_common *Common
)
{
    Int *Wi ;
    int ok ;

    if (CHOLMOD(dump) < -1)
    {
	/* no checks if debug level is -2 or less */
	return (TRUE) ;
    }
    RETURN_IF_NULL_COMMON (FALSE) ;
    Wi = malloc (MAX (1, n) * sizeof (Int)) ;
    ok = check_perm (Wi, CHOLMOD(dump), name, Perm, len, n,Common) ;
    if (Wi != NULL) free (Wi) ;
    return (ok) ;
}


/* ========================================================================== */
/* === cholmod_dump_dense =================================================== */
/* ========================================================================== */

int CHOLMOD(dump_dense)
(
    cholmod_dense *X,
    const char *name,
    cholmod_common *Common
)
{
    if (CHOLMOD(dump) < -1)
    {
	/* no checks if debug level is -2 or less */
	return (TRUE) ;
    }
    RETURN_IF_NULL_COMMON (FALSE) ;
    return (check_dense (CHOLMOD(dump), name, X, Common)) ;
}


/* ========================================================================== */
/* === cholmod_dump_triplet ================================================= */
/* ========================================================================== */

int CHOLMOD(dump_triplet)
(
    cholmod_triplet *T,
    const char *name,
    cholmod_common *Common
)
{
    if (CHOLMOD(dump) < -1)
    {
	/* no checks if debug level is -2 or less */
	return (TRUE) ;
    }
    RETURN_IF_NULL_COMMON (FALSE) ;
    return (check_triplet (CHOLMOD(dump), name, T, Common)) ;
}


/* ========================================================================== */
/* === cholmod_dump_subset ================================================== */
/* ========================================================================== */

int CHOLMOD(dump_subset)
(
    Int *S,
    size_t len,
    size_t n,
    const char *name,
    cholmod_common *Common
)
{
    if (CHOLMOD(dump) < -1)
    {
	/* no checks if debug level is -2 or less */
	return (TRUE) ;
    }
    RETURN_IF_NULL_COMMON (FALSE) ;
    return (check_subset (S, len, n, CHOLMOD(dump), name, Common)) ;
}


/* ========================================================================== */
/* === cholmod_dump_parent ================================================== */
/* ========================================================================== */

int CHOLMOD(dump_parent)
(
    Int *Parent,
    size_t n,
    const char *name,
    cholmod_common *Common
)
{
    if (CHOLMOD(dump) < -1)
    {
	/* no checks if debug level is -2 or less */
	return (TRUE) ;
    }
    RETURN_IF_NULL_COMMON (FALSE) ;
    return (check_parent (Parent, n, CHOLMOD(dump), name, Common)) ;
}


/* ========================================================================== */
/* === cholmod_dump_real ==================================================== */
/* ========================================================================== */

void CHOLMOD(dump_real)
(
    const char *name,
    Real *X, SuiteSparse_long nrow, SuiteSparse_long ncol, int lower,
    int xentry, cholmod_common *Common
)
{
    /* dump an nrow-by-ncol real dense matrix */
    SuiteSparse_long i, j ;
    double x, z ;
    if (CHOLMOD(dump) < -1)
    {
	/* no checks if debug level is -2 or less */
	return ;
    }
    PRINT1 (("%s: dump_real, nrow: %ld ncol: %ld lower: %d\n",
		name, nrow, ncol, lower)) ;
    for (j = 0 ; j < ncol ; j++)
    {
	PRINT2 (("    col %ld\n", j)) ;
	for (i = 0 ; i < nrow ; i++)
	{
	    /* X is stored in column-major form */
	    if (lower && i < j)
	    {
		PRINT2 (("        %5ld: -", i)) ;
	    }
	    else
	    {
		x = *X ;
		PRINT2 (("        %5ld: %e", i, x)) ;
		if (xentry == 2)
		{
		    z = *(X+1) ;
		    PRINT2 ((", %e", z)) ;
		}
	    }
	    PRINT2 (("\n")) ;
	    X += xentry ;
	}
    }
}


/* ========================================================================== */
/* === cholmod_dump_super =================================================== */
/* ========================================================================== */

void CHOLMOD(dump_super)
(
    SuiteSparse_long s,
    Int *Super, Int *Lpi, Int *Ls, Int *Lpx, double *Lx,
    int xentry,
    cholmod_common *Common
)
{
    Int k1, k2, do_values, psi, psx, nsrow, nscol, psend, ilast, p, i ;
    if (CHOLMOD(dump) < -1)
    {
	/* no checks if debug level is -2 or less */
	return ;
    }
    k1 = Super [s] ;
    k2 = Super [s+1] ;
    nscol = k2 - k1 ;
    do_values = (Lpx != NULL) && (Lx != NULL) ;
    psi = Lpi [s] ;
    psend = Lpi [s+1] ;
    nsrow = psend - psi ;
    PRINT1 (("\nSuper %ld, columns "ID" to "ID", "ID" rows "ID" cols\n",
		s, k1, k2-1, nsrow, nscol)) ;
    ilast = -1 ;
    for (p = psi ; p < psend ; p++)
    {
	i = Ls [p] ;
	PRINT2 (("  "ID" : p-psi "ID"\n", i, p-psi)) ;
	ASSERT (IMPLIES (p-psi < nscol, i == k1 + (p-psi))) ;
	if (p-psi == nscol-1) PRINT2 (("------\n")) ;
	ASSERT (i > ilast) ;
	ilast = i ;
    }
    if (do_values)
    {
	psx = Lpx [s] ;
	CHOLMOD(dump_real) ("Supernode", Lx + xentry*psx, nsrow, nscol, TRUE, 
		xentry, Common) ;
    }
}


/* ========================================================================== */
/* === cholmod_dump_mem ===================================================== */
/* ========================================================================== */

int CHOLMOD(dump_mem)
(
    const char *where,
    SuiteSparse_long should,
    cholmod_common *Common
)
{
    SuiteSparse_long diff = should - Common->memory_inuse ;
    if (diff != 0)
    {
	PRINT0 (("mem: %-15s peak %10g inuse %10g should %10g\n",
	    where, (double) Common->memory_usage, (double) Common->memory_inuse,
	    (double) should)) ;
	PRINT0 (("mem: %s diff %ld !\n", where, diff)) ;
    }
    return (diff == 0) ;
}


/* ========================================================================== */
/* === cholmod_dump_partition =============================================== */
/* ========================================================================== */

/* make sure we have a proper separator (for debugging only)
 *
 * workspace: none
 */

int CHOLMOD(dump_partition)
(
    SuiteSparse_long n,
    Int *Cp,
    Int *Ci,
    Int *Cnw,
    Int *Part,
    SuiteSparse_long sepsize,
    cholmod_common *Common
)
{
    Int chek [3], which, ok, i, j, p ;
    PRINT1 (("bisect sepsize %ld\n", sepsize)) ;
    ok = TRUE ;
    chek [0] = 0 ;
    chek [1] = 0 ;
    chek [2] = 0 ;
    for (j = 0 ; j < n ; j++)
    {
	PRINT2 (("--------j "ID" in part "ID" nw "ID"\n", j, Part [j], Cnw[j]));
	which = Part [j] ;
	for (p = Cp [j] ; p < Cp [j+1] ; p++)
	{
	    i = Ci [p] ;
	    PRINT3 (("i "ID", part "ID"\n", i, Part [i])) ;
	    if (which == 0)
	    {
		if (Part [i] == 1)
		{
		    PRINT0 (("Error! "ID" "ID"\n", i, j)) ;
		    ok = FALSE ;
		}
	    }
	    else if (which == 1)
	    {
		if (Part [i] == 0)
		{
		    PRINT0 (("Error! "ID" "ID"\n", i, j)) ;
		    ok = FALSE ;
		}
	    }
	}
	if (which < 0 || which > 2)
	{
	    PRINT0 (("Part out of range\n")) ;
	    ok = FALSE ;
	}
	chek [which] += Cnw [j] ;
    }
    PRINT1 (("sepsize %ld check "ID" "ID" "ID"\n",
		sepsize, chek[0], chek[1],chek[2]));
    if (sepsize != chek[2])
    {
	PRINT0 (("mismatch!\n")) ;
	ok = FALSE ;
    }
    return (ok) ;
}


/* ========================================================================== */
/* === cholmod_dump_work ==================================================== */
/* ========================================================================== */

int CHOLMOD(dump_work) (int flag, int head, SuiteSparse_long wsize,
    cholmod_common *Common)
{
    double *W ;
    Int *Flag, *Head ;
    Int k, nrow, mark ;

    if (CHOLMOD(dump) < -1)
    {
	/* no checks if debug level is -2 or less */
	return (TRUE) ;
    }

    RETURN_IF_NULL_COMMON (FALSE) ;
    nrow = Common->nrow ;
    Flag = Common->Flag ;
    Head = Common->Head ;
    W = Common->Xwork ;
    mark = Common->mark ;

    if (wsize < 0)
    {
	/* check all of Xwork */
	wsize = Common->xworksize ;
    }
    else
    {
	/* check on the first wsize doubles in Xwork */
	wsize = MIN (wsize, (Int) (Common->xworksize)) ;
    }

    if (flag)
    {
	for (k = 0 ; k < nrow ; k++)
	{
	    if (Flag [k] >= mark)
	    {
		PRINT0 (("Flag invalid, Flag ["ID"] = "ID", mark = "ID"\n",
			    k, Flag [k], mark)) ;
		ASSERT (0) ;
		return (FALSE) ;
	    }
	}
    }

    if (head)
    {
	for (k = 0 ; k < nrow ; k++)
	{
	    if (Head [k] != EMPTY)
	    {
		PRINT0 (("Head invalid, Head ["ID"] = "ID"\n", k, Head [k])) ;
		ASSERT (0) ;
		return (FALSE) ;
	    }
	}
    }

    for (k = 0 ; k < wsize ; k++)
    {
	if (W [k] != 0.)
	{
	    PRINT0 (("W invalid, W ["ID"] = %g\n", k, W [k])) ;
	    ASSERT (0) ;
	    return (FALSE) ;
	}
    }

    return (TRUE) ;
}
#endif
#endif
