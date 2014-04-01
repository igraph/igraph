/* ========================================================================== */
/* === Check/cholmod_write ================================================== */
/* ========================================================================== */

/* Write a matrix to a file in Matrix Market form.
 *
 * A can be sparse or full.
 *
 * If present and non-empty, A and Z must have the same dimension.  Z contains
 * the explicit zero entries in the matrix (which MATLAB drops).  The entries
 * of Z appear as explicit zeros in the output file.  Z is optional.  If it is
 * an empty matrix it is ignored.  Z must be sparse or empty, if present.
 * It is ignored if A is full.
 *
 * filename is the name of the output file.  comments is file whose
 * contents are include after the Matrix Market header and before the first
 * data line.  Ignored if an empty string or not present.
 *
 * Except for the workspace used by cholmod_symmetry (ncol integers) for
 * the sparse case, these routines use no workspace at all.
 */

#ifndef NCHECK

#include "cholmod_internal.h"
#include "cholmod_check.h"
#include "cholmod_matrixops.h"
#include <string.h>
#include <ctype.h>

#define MMLEN 1024
#define MAXLINE MMLEN+6

/* ========================================================================== */
/* === include_comments ===================================================== */
/* ========================================================================== */

/* Read in the comments file, if it exists, and copy it to the Matrix Market
 * file.  A "%" is prepended to each line.  Returns TRUE if successful, FALSE
 * otherwise.
 */

static int include_comments (FILE *f, const char *comments)
{
    FILE *cf = NULL ;
    char buffer [MAXLINE] ;
    int ok = TRUE ;
    if (comments != NULL && comments [0] != '\0')
    {
	cf = fopen (comments, "r") ;
	if (cf == NULL)
	{
	    return (FALSE) ;
	}
	while (ok && fgets (buffer, MAXLINE, cf) != NULL)
	{
	    /* ensure the line is not too long */
	    buffer [MMLEN-1] = '\0' ;
	    buffer [MMLEN-2] = '\n' ;
	    ok = ok && (fprintf (f, "%%%s", buffer) > 0) ;
	}
	fclose (cf) ;
    }
    return (ok) ;
}


/* ========================================================================== */
/* === get_value ============================================================ */
/* ========================================================================== */

/* Get the pth value in the matrix. */

static void get_value
(
    double *Ax,	    /* real values, or real/imag. for CHOLMOD_COMPLEX type */
    double *Az,	    /* imaginary values for CHOLMOD_ZOMPLEX type */
    Int p,	    /* get the pth entry */
    Int xtype,	    /* A->xtype: pattern, real, complex, or zomplex */
    double *x,	    /* the real part */
    double *z	    /* the imaginary part */
)
{
    switch (xtype)
    {
	case CHOLMOD_PATTERN:
	    *x = 1 ;
	    *z = 0 ;
	    break ;

	case CHOLMOD_REAL:
	    *x = Ax [p] ;
	    *z = 0 ;
	    break ;

	case CHOLMOD_COMPLEX:
	    *x = Ax [2*p] ;
	    *z = Ax [2*p+1] ;
	    break ;

	case CHOLMOD_ZOMPLEX:
	    *x = Ax [p] ;
	    *z = Az [p] ;
	    break ;
    }
}


/* ========================================================================== */
/* === print_value ========================================================== */
/* ========================================================================== */

/* Print a numeric value to the file, using the shortest format that ensures
 * the value is written precisely.  Returns TRUE if successful, FALSE otherwise.
 */ 

static int print_value
(
    FILE *f,	    /* file to print to */
    double x,	    /* value to print */
    Int is_integer  /* TRUE if printing as an integer */
)
{
    double y ;
    char s [MAXLINE], *p ;
    Int i, dest = 0, src = 0 ;
    int width, ok ;

    if (is_integer)
    {
	i = (Int) x ;
	ok = (fprintf (f, ID, i) > 0) ;
	return (ok) ;
    }

    /* ---------------------------------------------------------------------- */
    /* handle Inf and NaN */
    /* ---------------------------------------------------------------------- */

    /* change -inf to -HUGE_DOUBLE, and change +inf and nan to +HUGE_DOUBLE */
    if (CHOLMOD_IS_NAN (x) || x >= HUGE_DOUBLE)
    {
	x = HUGE_DOUBLE ;
    }
    else if (x <= -HUGE_DOUBLE)
    {
	x = -HUGE_DOUBLE ;
    }

    /* ---------------------------------------------------------------------- */
    /* find the smallest acceptable precision */
    /* ---------------------------------------------------------------------- */

    for (width = 6 ; width < 20 ; width++)
    {
	sprintf (s, "%.*g", width, x) ;
	sscanf (s, "%lg", &y) ;
	if (x == y) break ;
    }

    /* ---------------------------------------------------------------------- */
    /* shorten the string */
    /* ---------------------------------------------------------------------- */

    /* change "e+0" to "e", change "e+" to "e", and change "e-0" to "e-" */
    for (i = 0 ; i < MAXLINE && s [i] != '\0' ; i++)
    {
	if (s [i] == 'e')
	{
	    if (s [i+1] == '+')
	    {
		dest = i+1 ;
		if (s [i+2] == '0')
		{
		    /* delete characters s[i+1] and s[i+2] */
		    src = i+3 ;
		}
		else
		{
		    /* delete characters s[i+1] */
		    src = i+2 ;
		}
	    }
	    else if (s [i+1] == '-')
	    {
		dest = i+2 ;
		if (s [i+2] == '0')
		{
		    /* delete character s[i+2] */
		    src = i+3 ;
		}
		else
		{
		    /* no change */
		    break ;
		}
	    }
	    while (s [src] != '\0')
	    {
		s [dest++] = s [src++] ;
	    }
	    s [dest] = '\0' ;
	    break ;
	}
    }

    /* delete the leading "0" if present and not necessary */
    p = s ;
    s [MAXLINE-1] = '\0' ;
    i = strlen (s) ;
    if (i > 2 && s [0] == '0' && s [1] == '.')
    {
	/* change "0.x" to ".x" */
	p = s + 1 ;
    }
    else if (i > 3 && s [0] == '-' && s [1] == '0' && s [2] == '.')
    {
	/* change "-0.x" to "-.x" */
	s [1] = '-' ;
	p = s + 1 ;
    }

#if 0
    /* double-check */
    i = sscanf (p, "%lg", &z) ;
    if (i != 1 || y != z)
    {
	/* oops! something went wrong in the "e+0" edit, above. */
	/* this "cannot" happen */
	sprintf (s, "%.*g", width, x) ;
	p = s ;
    }
#endif

    /* ---------------------------------------------------------------------- */
    /* print the value to the file */
    /* ---------------------------------------------------------------------- */

    ok = (fprintf (f, "%s", p) > 0) ;
    return (ok) ;
}


/* ========================================================================== */
/* === print_triplet ======================================================== */
/* ========================================================================== */

/* Print a triplet, converting it to one-based.  Returns TRUE if successful,
 * FALSE otherwise.
 */

static int print_triplet
(
    FILE *f,		/* file to print to */
    Int is_binary,	/* TRUE if file is "pattern" */
    Int is_complex,	/* TRUE if file is "complex" */
    Int is_integer,	/* TRUE if file is "integer" */
    Int i,		/* row index (zero-based) */
    Int j,		/* column index (zero-based) */
    double x,		/* real part */
    double z		/* imaginary part */
)
{
    int ok ; 
    ok = (fprintf (f, ID " " ID, 1+i, 1+j) > 0) ;
    if (!is_binary)
    {
	fprintf (f, " ") ;
	ok = ok && print_value (f, x, is_integer) ;
	if (is_complex)
	{
	    fprintf (f, " ") ;
	    ok = ok && print_value (f, z, is_integer) ;
	}
    }
    ok = ok && (fprintf (f, "\n") > 0) ;
    return (ok) ;
}


/* ========================================================================== */
/* === ntriplets ============================================================ */
/* ========================================================================== */

/* Compute the number of triplets that will be printed to the file
 * from the matrix A. */

static Int ntriplets
(
    cholmod_sparse *A,	    /* matrix that will be printed */
    Int is_sym		    /* TRUE if the file is symmetric (lower part only)*/
)
{
    Int *Ap, *Ai, *Anz, packed, i, j, p, pend, ncol, stype, nz = 0 ;
    if (A == NULL)
    {
	/* the Z matrix is NULL */
	return (0) ;
    }
    stype = A->stype ;
    Ap = A->p ;
    Ai = A->i ;
    Anz = A->nz ;
    packed = A->packed ;
    ncol = A->ncol ;
    for (j = 0 ; j < ncol ; j++)
    {
	p = Ap [j] ;
	pend = (packed) ? Ap [j+1] : p + Anz [j] ;
	for ( ; p < pend ; p++)
	{
	    i = Ai [p] ;
	    if ((stype < 0 && i >= j) || (stype == 0 && (i >= j || !is_sym)))
	    {
		/* CHOLMOD matrix is symmetric-lower (and so is the file);
		 * or CHOLMOD matrix is unsymmetric and either A(i,j) is in
		 * the lower part or the file is unsymmetric. */
		nz++ ;
	    }
	    else if (stype > 0 && i <= j)
	    {
		/* CHOLMOD matrix is symmetric-upper, but the file is
		 * symmetric-lower.  Need to transpose the entry. */
		nz++ ;
	    }
	}
    }
    return (nz) ;
}


/* ========================================================================== */
/* === cholmod_write_sparse ================================================= */
/* ========================================================================== */

/* Write a sparse matrix to a file in Matrix Market format.   Optionally include
 * comments, and print explicit zero entries given by the pattern of the Z
 * matrix.  If not NULL, the Z matrix must have the same dimensions and stype
 * as A.
 *
 * Returns the symmetry in which the matrix was printed (1 to 7, see the
 * CHOLMOD_MM_* codes in CHOLMOD/Include/cholmod_core.h), or -1 on failure.
 *
 * If A and Z are sorted on input, and either unsymmetric (stype = 0) or
 * symmetric-lower (stype < 0), and if A and Z do not overlap, then the triplets
 * are sorted, first by column and then by row index within each column, with
 * no duplicate entries.  If all the above holds except stype > 0, then the
 * triplets are sorted by row first and then column.
 */

int CHOLMOD(write_sparse)
(
    /* ---- input ---- */
    FILE *f,		    /* file to write to, must already be open */
    cholmod_sparse *A,	    /* matrix to print */
    cholmod_sparse *Z,	    /* optional matrix with pattern of explicit zeros */
    const char *comments,    /* optional filename of comments to include */
    /* --------------- */
    cholmod_common *Common
)
{
    double x = 0, z = 0 ;
    double *Ax, *Az ;
    Int *Ap, *Ai, *Anz, *Zp, *Zi, *Znz ;
    Int nrow, ncol, is_complex, symmetry, i, j, q, iz, p, nz, is_binary, stype,
	is_integer, asym, is_sym, xtype, apacked, zpacked, pend, qend, zsym ;
    int ok ;

    /* ---------------------------------------------------------------------- */
    /* check inputs */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (EMPTY) ;
    RETURN_IF_NULL (f, EMPTY) ;
    RETURN_IF_NULL (A, EMPTY) ;
    RETURN_IF_XTYPE_INVALID (A, CHOLMOD_PATTERN, CHOLMOD_ZOMPLEX, EMPTY) ;
    if (Z != NULL && (Z->nrow == 0 || Z->ncol == 0))
    {
	/* Z is non-NULL but empty, so treat it as a NULL matrix */
	Z = NULL ;
    }
    if (Z != NULL)
    {
	RETURN_IF_XTYPE_INVALID (Z, CHOLMOD_PATTERN, CHOLMOD_ZOMPLEX, EMPTY) ;
	if (Z->nrow != A->nrow || Z->ncol != A->ncol || Z->stype != A->stype)
	{
	    ERROR (CHOLMOD_INVALID, "dimension or type of A and Z mismatch") ;
	    return (EMPTY) ;
	}
    }
    Common->status = CHOLMOD_OK ;

    /* ---------------------------------------------------------------------- */
    /* get the A matrix */
    /* ---------------------------------------------------------------------- */

    Ap = A->p ;
    Ai = A->i ;
    Ax = A->x ;
    Az = A->z ;
    Anz = A->nz ;
    nrow = A->nrow ;
    ncol = A->ncol ;
    xtype = A->xtype ;
    apacked = A->packed ;

    if (xtype == CHOLMOD_PATTERN)
    {
	/* a CHOLMOD pattern matrix is printed as "pattern" in the file */
	is_binary = TRUE ;
	is_integer = FALSE ;
	is_complex = FALSE ;
    }
    else if (xtype == CHOLMOD_REAL)
    {
	/* determine if a real matrix is in fact binary or integer */
	is_binary = TRUE ;
	is_integer = TRUE ;
	is_complex = FALSE ;
	for (j = 0 ; (is_binary || is_integer) && j < ncol ; j++)
	{
	    p = Ap [j] ;
	    pend = (apacked) ? Ap [j+1] : p + Anz [j] ;
	    for ( ; (is_binary || is_integer) && p < pend ; p++)
	    {
		x = Ax [p] ;
		if (x != 1)
		{
		    is_binary = FALSE ;
		}
		/* convert to Int and then back to double */
		i = (Int) x ;
		z = (double) i ;
		if (z != x)
		{
		    is_integer = FALSE ;
		}
	    }
	}
    }
    else
    {
	/* a CHOLMOD complex matrix is printed as "complex" in the file */
	is_binary = FALSE ;
	is_integer = FALSE ;
	is_complex = TRUE ;
    }

    /* ---------------------------------------------------------------------- */
    /* get the Z matrix (only consider the pattern) */
    /* ---------------------------------------------------------------------- */

    Zp = NULL ;
    Zi = NULL ;
    Znz = NULL ;
    zpacked = TRUE ;
    if (Z != NULL)
    {
	Zp = Z->p ;
	Zi = Z->i ;
	Znz = Z->nz ;
	zpacked = Z->packed ;
    }

    /* ---------------------------------------------------------------------- */
    /* determine the symmetry of A and Z */
    /* ---------------------------------------------------------------------- */

    stype = A->stype ;
    if (A->nrow != A->ncol)
    {
	asym = CHOLMOD_MM_RECTANGULAR ;
    }
    else if (stype != 0)
    {
	/* CHOLMOD's A and Z matrices have a symmetric (and matching) stype.
	 * Note that the diagonal is not checked. */
	asym = is_complex ? CHOLMOD_MM_HERMITIAN : CHOLMOD_MM_SYMMETRIC ;
    }
    else if (!A->sorted)
    {
	/* A is in unsymmetric storage, but unsorted */
	asym = CHOLMOD_MM_UNSYMMETRIC ;
    }
    else
    {
	/* CHOLMOD's stype is zero (stored in unsymmetric form) */
	asym = EMPTY ;
	zsym = EMPTY ;

#ifndef NMATRIXOPS
	/* determine if the matrices are in fact symmetric or Hermitian */
	asym = CHOLMOD(symmetry) (A, 1, NULL, NULL, NULL, NULL, Common) ;
	zsym = (Z == NULL) ? 999 :
	       CHOLMOD(symmetry) (Z, 1, NULL, NULL, NULL, NULL, Common) ;
#endif

	if (asym == EMPTY || zsym <= CHOLMOD_MM_UNSYMMETRIC)
	{
	    /* not computed, out of memory, or Z is unsymmetric */
	    asym = CHOLMOD_MM_UNSYMMETRIC ;
	}
    }

    /* ---------------------------------------------------------------------- */
    /* write the Matrix Market header */
    /* ---------------------------------------------------------------------- */

    ok = fprintf (f, "%%%%MatrixMarket matrix coordinate") > 0 ;

    if (is_complex)
    {
	ok = ok && (fprintf (f, " complex") > 0) ;
    }
    else if (is_binary)
    {
	ok = ok && (fprintf (f, " pattern") > 0) ;
    }
    else if (is_integer)
    {
	ok = ok && (fprintf (f, " integer") > 0) ;
    }
    else
    {
	ok = ok && (fprintf (f, " real") > 0) ;
    }

    is_sym = FALSE ;

    switch (asym)
    {
	case CHOLMOD_MM_RECTANGULAR:
	case CHOLMOD_MM_UNSYMMETRIC:
	    /* A is rectangular or unsymmetric */
	    ok = ok && (fprintf (f, " general\n") > 0) ;
	    is_sym = FALSE ;
	    symmetry = CHOLMOD_MM_UNSYMMETRIC ;
	    break ;

	case CHOLMOD_MM_SYMMETRIC:
	case CHOLMOD_MM_SYMMETRIC_POSDIAG:
	    /* A is symmetric */
	    ok = ok && (fprintf (f, " symmetric\n") > 0) ;
	    is_sym = TRUE ;
	    symmetry = CHOLMOD_MM_SYMMETRIC ;
	    break ;

	case CHOLMOD_MM_HERMITIAN:
	case CHOLMOD_MM_HERMITIAN_POSDIAG:
	    /* A is Hermitian */
	    ok = ok && (fprintf (f, " Hermitian\n") > 0) ;
	    is_sym = TRUE ;
	    symmetry = CHOLMOD_MM_HERMITIAN ;
	    break ;

	case CHOLMOD_MM_SKEW_SYMMETRIC:
	    /* A is skew symmetric */
	    ok = ok && (fprintf (f, " skew-symmetric\n") > 0) ;
	    is_sym = TRUE ;
	    symmetry = CHOLMOD_MM_SKEW_SYMMETRIC ;
	    break ;
    }

    /* ---------------------------------------------------------------------- */
    /* include the comments if present */
    /* ---------------------------------------------------------------------- */

    ok = ok && include_comments (f, comments) ;

    /* ---------------------------------------------------------------------- */
    /* write a sparse matrix (A and Z) */
    /* ---------------------------------------------------------------------- */

    nz = ntriplets (A, is_sym) + ntriplets (Z, is_sym) ;

    /* write the first data line, with nrow, ncol, and # of triplets */
    ok = ok && (fprintf (f, ID " " ID " " ID "\n", nrow, ncol, nz) > 0) ;

    for (j = 0 ; ok && j < ncol ; j++)
    {
	/* merge column of A and Z */
	p = Ap [j] ;
	pend = (apacked) ? Ap [j+1] : p + Anz [j] ;
	q = (Z == NULL) ? 0 : Zp [j] ;
	qend = (Z == NULL) ? 0 : ((zpacked) ? Zp [j+1] : q + Znz [j]) ;
	while (ok)
	{
	    /* get the next row index from A and Z */
	    i  = (p < pend) ? Ai [p] : (nrow+1) ;
	    iz = (q < qend) ? Zi [q] : (nrow+2) ;
	    if (i <= iz)
	    {
		/* get A(i,j), or quit if both A and Z are exhausted */
		if (i == nrow+1) break ;
		get_value (Ax, Az, p, xtype, &x, &z) ;
		p++ ;
	    }
	    else
	    {
		/* get Z(i,j) */
		i = iz ;
		x = 0 ;
		z = 0 ;
		q++ ;
	    }
	    if ((stype < 0 && i >= j) || (stype == 0 && (i >= j || !is_sym)))
	    {
		/* CHOLMOD matrix is symmetric-lower (and so is the file);
		 * or CHOLMOD matrix is unsymmetric and either A(i,j) is in
		 * the lower part or the file is unsymmetric. */
		ok = ok && print_triplet (f, is_binary, is_complex, is_integer,
		    i,j, x,z) ;
	    }
	    else if (stype > 0 && i <= j)
	    {
		/* CHOLMOD matrix is symmetric-upper, but the file is
		 * symmetric-lower.  Need to transpose the entry.   If the
		 * matrix is real, the complex part is ignored.  If the matrix
		 * is complex, it Hermitian.
		 */
		ASSERT (IMPLIES (is_complex, asym == CHOLMOD_MM_HERMITIAN)) ;
		if (z != 0)
		{
		    z = -z ;
		}
		ok = ok && print_triplet (f, is_binary, is_complex, is_integer,
		    j,i, x,z) ;
	    }
	}
    }

    if (!ok)
    {
	ERROR (CHOLMOD_INVALID, "error reading/writing file") ;
	return (EMPTY) ;
    }

    return (asym) ;
}


/* ========================================================================== */
/* === cholmod_write_dense ================================================== */
/* ========================================================================== */

/* Write a dense matrix to a file in Matrix Market format.   Optionally include
 * comments.  Returns > 0 if successful, -1 otherwise (1 if rectangular, 2 if
 * square).  Future versions may return 1 to 7 on success (a CHOLMOD_MM_* code,
 * just as cholmod_write_sparse does).
 *
 * A dense matrix is written in "general" format; symmetric formats in the
 * Matrix Market standard are not exploited.
 */

int CHOLMOD(write_dense)
(
    /* ---- input ---- */
    FILE *f,		    /* file to write to, must already be open */
    cholmod_dense *X,	    /* matrix to print */
    const char *comments,    /* optional filename of comments to include */
    /* --------------- */
    cholmod_common *Common
)
{
    double x = 0, z = 0 ;
    double *Xx, *Xz ;
    Int nrow, ncol, is_complex, i, j, xtype, p ;
    int ok ;

    /* ---------------------------------------------------------------------- */
    /* check inputs */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (EMPTY) ;
    RETURN_IF_NULL (f, EMPTY) ;
    RETURN_IF_NULL (X, EMPTY) ;
    RETURN_IF_XTYPE_INVALID (X, CHOLMOD_REAL, CHOLMOD_ZOMPLEX, EMPTY) ;
    Common->status = CHOLMOD_OK ;

    /* ---------------------------------------------------------------------- */
    /* get the X matrix */
    /* ---------------------------------------------------------------------- */

    Xx = X->x ;
    Xz = X->z ;
    nrow = X->nrow ;
    ncol = X->ncol ;
    xtype = X->xtype ;
    is_complex = (xtype == CHOLMOD_COMPLEX) || (xtype == CHOLMOD_ZOMPLEX) ;

    /* ---------------------------------------------------------------------- */
    /* write the Matrix Market header */
    /* ---------------------------------------------------------------------- */

    ok = (fprintf (f, "%%%%MatrixMarket matrix array") > 0) ;
    if (is_complex)
    {
	ok = ok && (fprintf (f, " complex general\n") > 0) ;
    }
    else
    {
	ok = ok && (fprintf (f, " real general\n") > 0) ;
    }

    /* ---------------------------------------------------------------------- */
    /* include the comments if present */
    /* ---------------------------------------------------------------------- */

    ok = ok && include_comments (f, comments) ;

    /* ---------------------------------------------------------------------- */
    /* write a dense matrix */
    /* ---------------------------------------------------------------------- */

    /* write the first data line, with nrow and ncol */
    ok = ok && (fprintf (f, ID " " ID "\n", nrow, ncol) > 0) ;

    Xx = X->x ;
    Xz = X->z ;
    for (j = 0 ; ok && j < ncol ; j++)
    {
	for (i = 0 ; ok && i < nrow ; i++)
	{
	    p = i + j*nrow ;
	    get_value (Xx, Xz, p, xtype, &x, &z) ;
	    ok = ok && print_value (f, x, FALSE) ;
	    if (is_complex)
	    {
		ok = ok && (fprintf (f, " ") > 0) ;
		ok = ok && print_value (f, z, FALSE) ;
	    }
	    ok = ok && (fprintf (f, "\n") > 0) ;
	}
    }

    if (!ok)
    {
	ERROR (CHOLMOD_INVALID, "error reading/writing file") ;
	return (EMPTY) ;
    }

    return ((nrow == ncol) ? CHOLMOD_MM_UNSYMMETRIC : CHOLMOD_MM_RECTANGULAR) ;
}
#endif
