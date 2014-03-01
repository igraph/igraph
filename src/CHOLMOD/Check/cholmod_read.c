/* ========================================================================== */
/* === Check/cholmod_read =================================================== */
/* ========================================================================== */

/* -----------------------------------------------------------------------------
 * CHOLMOD/Check Module.  Copyright (C) 2005-2006, Timothy A. Davis.
 * The CHOLMOD/Check Module is licensed under Version 2.1 of the GNU
 * Lesser General Public License.  See lesser.txt for a text of the license.
 * CHOLMOD is also available under other licenses; contact authors for details.
 * -------------------------------------------------------------------------- */

/* Read a sparse matrix in triplet or dense form.  A triplet matrix can be
 * returned as compressed-column sparse matrix.  The file format is compatible
 * with all variations of the Matrix Market "coordinate" and "array" format
 * (http://www.nist.gov/MatrixMarket).  The format supported by these routines
 * also allow other formats, where the Matrix Market header is optional.
 *
 * Although the Matrix Market header is optional, I recommend that users stick
 * with the strict Matrix Market format.  The optional format appears here to
 * support the reading of symmetric matrices stored with just their upper
 * triangular parts present, for testing and development of the A->stype > 0
 * format in CHOLMOD.  That format is not included in the Matrix Market format.
 *
 * If the first line of the file starts with %%MatrixMarket, then it is
 * interpretted as a file in Matrix Market format.  This line must have
 * the following format:
 *
 *	%%MatrixMarket matrix <fmt> <type> <storage>
 *
 *	<fmt> is one of: coordinate or array.  The former is a sparse matrix in
 *	triplet form.  The latter is a dense matrix in column-major form.
 *
 *	<type> is one of: real, complex, pattern, or integer.
 *	The functions here convert the  "integer" and "pattern" types to real.  
 *
 *	<storage> is one of: general, hermitian, symmetric, or skew-symmetric
 *
 *	The strings are case-insensitive.  Only the first character is
 *	significant (or the first two for skew-symmetric).
 *
 *	<type> is ignored for all matrices; the actual type (real, complex,
 *	or pattern) is inferred from the number of tokens in each line of the
 *	file.  For a "coordinate" matrix: 2: pattern, 3: real, 4: complex; for
 *	a dense "array" matrix: 1: real, 2: complex.  This is compatible with
 *	the Matrix Market format, since pattern matrices must have two tokens
 *	per line, real matrices must have 3, and complex matrices must have 4.
 *	A storage of "general" implies an stype of zero (see below). 
 *	"symmetric" and "hermitian" imply an stype of -1. Skew-symmetric and
 *	complex symmetric matrices are always returned with both upper and lower
 *	triangular parts present, with an stype of zero, since CHOLMOD does not
 *	have a method for representing skew-symmetric and complex symmetric
 *	matrices.  Real symmetric and complex Hermitian matrices may optionally
 *	be returned with both parts present.
 *
 * Any other lines starting with "%" are treated as comments, and are ignored.
 * Blank lines are ignored.  The Matrix Market header is optional in this
 * routine (it is not optional in the Matrix Market format).
 *
 * Note that complex matrices are always returned in CHOLMOD_COMPLEX format,
 * not CHOLMOD_ZOMPLEX.
 *
 * -----------------------------------------------------------------------------
 * Triplet matrices:
 * -----------------------------------------------------------------------------
 *
 * The first data line of a triplet matrix contains 3 or 4 integers:
 *
 *	nrow ncol nnz stype
 *
 * where stype is optional (stype does not appear in the Matrix Market format).
 * The matrix is nrow-by-ncol.  The following nnz lines (excluding comments
 * and blank lines) each contain a single entry.  Duplicates are permitted,
 * and are summed in the output matrix.
 *
 * The stype is first derived from the Matrix Market header.  If the stype
 * appears as the fourth integer in the first data line, it is determined from
 * that line.
 * 
 * If stype is present, it denotes the storage format for the matrix.
 * stype = 0 denotes an unsymmetric matrix (same as Matrix Market "general").
 * stype = -1 denotes a real symmetric or complex Hermitian matrix whose lower
 *	triangular entries are stored.  Entries may be present in the upper
 *	triangular part, but these are ignored (same as Matrix Market
 *	"real symmetric" and "complex Hermitian").
 * stype = 1 denotes a real symmetric or complex Hermitian matrix whose upper
 *	triangular entries are stored.  Entries may be present in the lower
 *	triangular part, but these are ignored.  This option is not present
 *	in the Matrix Market format.
 *
 * If stype is not present (no Matrix Market header and not in the first data
 * line) it is inferred from the rest of the data.  If the matrix is
 * rectangular, or has entries in both the upper and lower triangular parts,
 * then it is assumed to be unsymmetric (stype=0).  If only entries in the
 * lower triangular part are present, the matrix is assumed to have stype = -1.
 * If only entries in the upper triangular part are present, the matrix is
 * assumed to have stype = 1.
 *
 * After the first data line (with nrow, ncol, nnz, and optionally stype),
 * each nonzero consists of one line with 2, 3, or 4 entries.  All lines must
 * have the same number of entries.  The first two entries are the row and
 * column indices of the nonzero.  If 3 entries are present, the 3rd entry is
 * the numerical value, and the matrix is real.  If 4 entries are present,
 * the 3rd and 4th entries in the line are the real and imaginary parts of
 * a complex value.
 *
 * The matrix can be either 0-based or 1-based.  It is first assumed to be
 * one-based (all matrices in the Matrix Market are one-based), with row indices
 * in the range 1 to ncol and column indices in the range 1 to nrow.  If a row
 * or column index of zero is found, the matrix is assumed to be zero-based
 * (with row indices in the range 0 to ncol-1 and column indices in the range 0
 * to nrow-1).
 *
 * If Common->prefer_binary is set to its default value of FALSE, then
 * for symmetric pattern-only matrices, the kth diagonal (if present) is set to
 * one plus the degree of the row/column k, and the off-diagonal entries are set
 * to -1.  A symmetric pattern-only matrix with a zero-free diagonal is thus
 * converted into a symmetric positive definite matrix.  All entries are set to
 * one for an unsymmetric pattern-only matrix.  This differs from the
 * Matrix Market format (A = mmread ('file') returns a binary pattern for A for
 * symmetric pattern-only matrices).  If Common->prefer_binary is TRUE, then
 * this function returns a binary matrix (just like mmread('file')).
 *
 * -----------------------------------------------------------------------------
 * Dense matrices:
 * -----------------------------------------------------------------------------
 *
 * A dense matrix is specified by the Matrix Market "array" format.  The
 * Matrix Market header is optional; if not present, the matrix is assumed to
 * be in the Matrix Market "general" format.  The first data line contains just
 * two integers:
 *
 *	nrow ncol
 *
 * The <type> can be real, integer, or complex (not pattern).  These functions
 * convert an integer type to real.  The entries in the matrix are stored in
 * column-major format, with one line per entry.  Two entries are present in
 * each line for complex matrices, one for real and integer matrices.  In
 * rectangular and unsymmetric matrices, all entries are present.  For real
 * symmetric or complex Hermitian matrices, only entries in the lower triangular
 * part appear.  For skew-symmetric matrices, only entries in the strictly
 * lower triangular part appear.
 *
 * Since CHOLMOD does not have a data structure for presenting dense symmetric/
 * Hermitian matrices, these functions always return a dense matrix in its
 * general form, with both upper and lower parts present.
 */

#ifndef NCHECK

#include "cholmod_internal.h"
#include "cholmod_check.h"
#include <string.h>
#include <ctype.h>

/* The MatrixMarket format specificies a maximum line length of 1024 */
#define MAXLINE 1030

/* ========================================================================== */
/* === get_line ============================================================= */
/* ========================================================================== */

/* Read one line of the file, return TRUE if successful, FALSE if EOF. */

static int get_line (FILE *f, char *buf)
{
    buf [0] = '\0' ;
    buf [1] = '\0' ;
    buf [MAXLINE] = '\0' ;
    return (fgets (buf, MAXLINE, f) != NULL) ;
}

/* ========================================================================== */
/* === fix_inf ============================================================== */
/* ========================================================================== */

/* Replace huge values with +/- Inf's, since scanf and printf don't deal
 * with Inf's properly.
 */

static double fix_inf (double x)
{
    if ((x >= HUGE_DOUBLE) || (x <= -HUGE_DOUBLE))
    {
	/* treat this as +/- Inf (assume 2*x leads to overflow) */
	x = 2*x ;
    }
    return (x) ;
}

/* ========================================================================== */
/* === is_blank_line ======================================================== */
/* ========================================================================== */

/* TRUE if s is a blank line or comment, FALSE otherwise */

static int is_blank_line
(
    char *s
)
{
    int c, k ;
    if (s [0] == '%')
    {
	/* a comment line */
	return (TRUE) ;
    }
    for (k = 0 ; k <= MAXLINE ; k++)
    {
	c = s [k] ;
	if (c == '\0')
	{
	    /* end of line */
	    break ;
	}
	if (!isspace (c))
	{
	    /* non-space character */
	    return (FALSE) ;
	}
    }
    return (TRUE) ;
}


/* ========================================================================== */
/* === read_header ========================================================== */
/* ========================================================================== */

/* Read the header.  This consists of zero or more comment lines (blank, or
 * starting with a "%" in the first column), followed by a single data line
 * containing up to four numerical values.
 *
 * The first line may optionally be a Matrix Market header line, of the form
 *
 *	%%MatrixMarket matrix <fmt> <type> <storage>
 *
 * The first data line of a sparse matrix in triplet form consists of 3 or 4
 * numerical values:
 *
 *	nrow ncol nnz stype
 *
 * where stype is optional (it does not appear in the Matrix Market file
 * format).  The first line of a dense matrix in column-major form consists of
 * two numerical values:
 *
 *	nrow ncol
 *
 * The stype of the matrix is determine either from the Matrix Market header,
 * or (optionally) from the first data line.  stypes of 0 to -3 directly
 * correlate with the Matrix Market format; stype = 1 is an extension to that
 * format.
 *
 *	999: unknown (will be inferred from the data)
 *	1: real symmetric or complex Hermitian with upper part stored
 *		(not in the Matrix Market format)
 *	0: unsymmetric (same as Matrix Market "general")
 *	-1: real symmetric or complex Hermitian, with lower part stored
 *		(Matrix Market "real symmetric" or "complex hermitian")
 *	-2: real or complex skew symmetric (lower part stored, can only be
 *		specified by Matrix Market header)
 *	-3: complex symmetric (lower part stored)
 *		specified by Matrix Market header)
 *
 * The Matrix Market header is optional.  If stype appears in the first data
 * line, it is determine by that data line.  Otherwise, if the Matrix Market
 * header appears, stype is determined from that header.  If stype does not
 * appear, it is set to "unknown" (999).
 */

#define STYPE_UNKNOWN 999
#define STYPE_SYMMETRIC_UPPER 1
#define STYPE_UNSYMMETRIC 0
#define STYPE_SYMMETRIC_LOWER -1
#define STYPE_SKEW_SYMMETRIC -2
#define STYPE_COMPLEX_SYMMETRIC_LOWER -3

static int read_header	/* returns TRUE if successful, FALSE on error */
(
    /* ---- input ---- */
    FILE *f,		/* file to read from */
    /* ---- output --- */
    char *buf,		/* a character array of size MAXLINE+1 */
    int *mtype,		/* CHOLMOD_TRIPLET or CHOLMOD_DENSE */
    size_t *nrow,	/* number of rows in the matrix */
    size_t *ncol,	/* number of columns in the matrix */
    size_t *nnz,	/* number of entries in a triplet matrix (0 for dense)*/
    int *stype		/* stype (see above) */
)
{
    char *p ;
    int first = TRUE, got_mm_header = FALSE, c, c2, is_complex, nitems ;
    double l1, l2, l3, l4 ;

    *mtype = CHOLMOD_TRIPLET ;
    *nrow = 0 ;
    *ncol = 0 ;
    *nnz = 0 ;
    *stype = STYPE_UNKNOWN ;

    for ( ; ; )
    {

	/* ------------------------------------------------------------------ */
	/* get the next line */
	/* ------------------------------------------------------------------ */

	if (!get_line (f, buf))
	{
	    /* premature end of file */
	    return (FALSE) ;
	}

	if (first && (strncmp (buf, "%%MatrixMarket", 14) == 0))
	{

	    /* -------------------------------------------------------------- */
	    /* read a Matrix Market header */
	    /* -------------------------------------------------------------- */

	    got_mm_header = TRUE ;
	    p = buf ;

	    /* -------------------------------------------------------------- */
	    /* get "matrix" token */
	    /* -------------------------------------------------------------- */

	    while (*p && !isspace (*p)) p++ ;
	    while (*p &&  isspace (*p)) p++ ;
	    c = tolower (*p) ;
	    if (c != 'm')
	    {
		/* bad format */
		return (FALSE) ;
	    }

	    /* -------------------------------------------------------------- */
	    /* get the fmt token ("coord" or "array") */
	    /* -------------------------------------------------------------- */

	    while (*p && !isspace (*p)) p++ ;
	    while (*p &&  isspace (*p)) p++ ;
	    c = tolower (*p) ;
	    if (c == 'c')
	    {
		*mtype = CHOLMOD_TRIPLET  ;
	    }
	    else if (c == 'a')
	    {
		*mtype = CHOLMOD_DENSE  ;
	    }
	    else
	    {
		/* bad format, neither "coordinate" nor "array" */
		return (FALSE) ;
	    }

	    /* -------------------------------------------------------------- */
	    /* get type token (real, pattern, complex, integer) */
	    /* -------------------------------------------------------------- */

	    while (*p && !isspace (*p)) p++ ;
	    while (*p &&  isspace (*p)) p++ ;
	    c = tolower (*p) ;
	    if (!(c == 'r' || c == 'p' || c == 'c' || c == 'i'))
	    {
		/* bad format */
		return (FALSE) ;
	    }
	    is_complex = (c == 'c') ;

	    /* -------------------------------------------------------------- */
	    /* get storage token (general, hermitian, symmetric, skew) */
	    /* -------------------------------------------------------------- */

	    while (*p && !isspace (*p)) p++ ;
	    while (*p &&  isspace (*p)) p++ ;
	    c = tolower (*p) ;
	    c2 = tolower (*(p+1)) ;
	    if (c == 'g')
	    {
		/* "general" storage (unsymmetric matrix), both parts present */
		*stype = STYPE_UNSYMMETRIC ;
	    }
	    else if (c == 's' && c2 == 'y')
	    {
		/* "symmetric" */
		if (is_complex)
		{
		    /* complex symmetric, lower triangular part present */
		    *stype = STYPE_COMPLEX_SYMMETRIC_LOWER ;
		}
		else
		{
		    /* real symmetric, lower triangular part present */
		    *stype = STYPE_SYMMETRIC_LOWER ;
		}
	    }
	    else if (c == 'h')
	    {
		/* "hermitian" matrix, lower triangular part present */
		*stype = STYPE_SYMMETRIC_LOWER ;
	    }
	    else if (c == 's' && c2 == 'k')
	    {
		/* "skew-symmetric" (real or complex), lower part present */
		*stype = STYPE_SKEW_SYMMETRIC ;
	    }
	    else
	    {
		/* bad format */
		return (FALSE) ;
	    }

	}
	else if (is_blank_line (buf))
	{

	    /* -------------------------------------------------------------- */
	    /* blank line or comment line */
	    /* -------------------------------------------------------------- */

	    continue ;

	}
	else
	{

	    /* -------------------------------------------------------------- */
	    /* read the first data line and return */
	    /* -------------------------------------------------------------- */

	    /* format: nrow ncol nnz stype */
	    l1 = EMPTY ;
	    l2 = EMPTY ;
	    l3 = 0 ;
	    l4 = 0 ;
	    nitems = sscanf (buf, "%lg %lg %lg %lg\n", &l1, &l2, &l3, &l4) ;
	    if (nitems < 2 || nitems > 4 || l1 > Int_max || l2 > Int_max)
	    {
		/* invalid matrix */
		return (FALSE) ;
	    }
	    *nrow = l1 ;
	    *ncol = l2 ;
	    if (nitems == 2)
	    {
		/* a dense matrix */
		if (!got_mm_header)
		{
		    *mtype = CHOLMOD_DENSE ;
		    *stype = STYPE_UNSYMMETRIC ;
		}
	    }
	    if (nitems == 3 || nitems == 4)
	    {
		/* a sparse triplet matrix */
		*nnz = l3 ;
		if (!got_mm_header)
		{
		    *mtype = CHOLMOD_TRIPLET ;
		}
	    }
	    if (nitems == 4)
	    {
		/* an stype specified here can only be 1, 0, or -1 */
		if (l4 < 0)
		{
		    *stype = STYPE_SYMMETRIC_LOWER ;
		}
		else if (l4 > 0)
		{
		    *stype = STYPE_SYMMETRIC_UPPER ;
		}
		else
		{
		    *stype = STYPE_UNSYMMETRIC ;
		}
	    }
	    if (*nrow != *ncol)
	    {
		/* a rectangular matrix must be unsymmetric */
		*stype = STYPE_UNSYMMETRIC ;
	    }
	    return (TRUE) ;
	}

	first = FALSE ;
    }
}


/* ========================================================================== */
/* === read_triplet ========================================================= */
/* ========================================================================== */

/* Header has already been read in, including first line (nrow ncol nnz stype).
 * Read the triplets. */

static cholmod_triplet *read_triplet
(
    /* ---- input ---- */
    FILE *f,		    /* file to read from, must already be open */
    size_t nrow,	    /* number of rows */
    size_t ncol,	    /* number of columns */
    size_t nnz,		    /* number of triplets in file to read */
    int stype,		    /* stype from header, or "unknown" */
    int prefer_unsym,	    /* if TRUE, always return T->stype of zero */
    /* ---- workspace */
    char *buf,		    /* of size MAXLINE+1 */
    /* --------------- */
    cholmod_common *Common
)
{
    double x, z ;
    double *Tx ;
    Int *Ti, *Tj, *Rdeg, *Cdeg ;
    cholmod_triplet *T ;
    double l1, l2 ;
    Int nitems, xtype, unknown, k, nshould, is_lower, is_upper, one_based, i, j,
	imax, jmax, skew_symmetric, p, complex_symmetric ;
    size_t s, nnz2, extra ;
    int ok = TRUE ;

    /* ---------------------------------------------------------------------- */
    /* quick return for empty matrix */
    /* ---------------------------------------------------------------------- */

    if (nrow == 0 || ncol == 0 || nnz == 0)
    {
	/* return an empty matrix */
	return (CHOLMOD(allocate_triplet) (nrow, ncol, 0, 0, CHOLMOD_REAL,
		    Common)) ;
    }

    /* ---------------------------------------------------------------------- */
    /* special stype cases: unknown, skew symmetric, and complex symmetric  */
    /* ---------------------------------------------------------------------- */

    unknown = (stype == STYPE_UNKNOWN) ;
    skew_symmetric = (stype == STYPE_SKEW_SYMMETRIC) ;
    complex_symmetric = (stype == STYPE_COMPLEX_SYMMETRIC_LOWER) ;

    extra = 0 ;
    if (stype < STYPE_SYMMETRIC_LOWER
	|| (prefer_unsym && stype != STYPE_UNSYMMETRIC))
    {
	/* 999: unknown might be converted to unsymmetric */
	/*  1:  symmetric upper converted to unsym. if prefer_unsym is TRUE */
	/* -1:  symmetric lower converted to unsym. if prefer_unsym is TRUE */
	/* -2:  real or complex skew symmetric converted to unsymmetric */
	/* -3:  complex symmetric converted to unsymmetric */
	stype = STYPE_UNSYMMETRIC ;
	extra = nnz ;
    }
    nnz2 = CHOLMOD(add_size_t) (nnz, extra, &ok) ;

    /* ---------------------------------------------------------------------- */
    /* allocate workspace */
    /* ---------------------------------------------------------------------- */

    /* s = nrow + ncol */
    s = CHOLMOD(add_size_t) (nrow, ncol, &ok) ;
    if (!ok || nrow > Int_max || ncol > Int_max || nnz > Int_max)
    {
	ERROR (CHOLMOD_TOO_LARGE, "problem too large") ;
	return (NULL) ;
    }

    CHOLMOD(allocate_work) (0, s, 0, Common) ;
    Rdeg = Common->Iwork ;	/* size nrow */
    Cdeg = Rdeg + nrow ;	/* size ncol */

    /* ---------------------------------------------------------------------- */
    /* read the triplets */
    /* ---------------------------------------------------------------------- */

    is_lower = TRUE ;
    is_upper = TRUE ;
    one_based = TRUE ;
    imax = 0 ;
    jmax = 0 ;

    Tx = NULL ;
    Ti = NULL ;
    Tj = NULL ;
    xtype = 999 ;
    nshould = 0 ;

    for (k = 0 ; k < (Int) nnz ; k++)
    {

	/* ------------------------------------------------------------------ */
	/* get the next triplet, skipping blank lines and comment lines */
	/* ------------------------------------------------------------------ */

	l1 = EMPTY ;
	l2 = EMPTY ;
	x = 0 ;
	z = 0 ;

	for ( ; ; )
	{
	    if (!get_line (f, buf))
	    {
		/* premature end of file - not enough triplets read in */
		ERROR (CHOLMOD_INVALID, "premature EOF") ;
		return (NULL) ;
	    }
	    if (is_blank_line (buf))
	    {
		/* blank line or comment */
		continue ;
	    }
	    nitems = sscanf (buf, "%lg %lg %lg %lg\n", &l1, &l2, &x, &z) ;
	    x = fix_inf (x) ;
	    z = fix_inf (z) ;
	    break ;
	}

	nitems = (nitems == EOF) ? 0 : nitems ;
	i = l1 ;
	j = l2 ;

	/* ------------------------------------------------------------------ */
	/* for first triplet: determine type and allocate triplet matrix */
	/* ------------------------------------------------------------------ */

	if (k == 0)
	{
	    if (nitems < 2 || nitems > 4)
	    {
		/* invalid matrix */
		ERROR (CHOLMOD_INVALID, "invalid format") ;
		return (NULL) ;
	    }
	    else if (nitems == 2)
	    {
		/* this will be converted into a real matrix later */
		xtype = CHOLMOD_PATTERN ;
	    }
	    else if (nitems == 3)
	    {
		xtype = CHOLMOD_REAL ;
	    }
	    else if (nitems == 4)
	    {
		xtype = CHOLMOD_COMPLEX ;
	    }

	    /* the rest of the lines should have the same number of entries */
	    nshould = nitems ;

	    /* allocate triplet matrix */
	    T = CHOLMOD(allocate_triplet) (nrow, ncol, nnz2, stype,
		    (xtype == CHOLMOD_PATTERN ? CHOLMOD_REAL : xtype), Common) ;
	    if (Common->status < CHOLMOD_OK)
	    {
		/* out of memory */
		return (NULL) ;
	    }
	    Ti = T->i ;
	    Tj = T->j ;
	    Tx = T->x ;
	    T->nnz = nnz ;
	}

	/* ------------------------------------------------------------------ */
	/* save the entry in the triplet matrix */
	/* ------------------------------------------------------------------ */

	if (nitems != nshould || i < 0 || j < 0)
	{
	    /* wrong format, premature end-of-file, or negative indices */
	    CHOLMOD(free_triplet) (&T, Common) ;
	    ERROR (CHOLMOD_INVALID, "invalid matrix file") ;
	    return (NULL) ;
	}

	Ti [k] = i ;
	Tj [k] = j ;

	if (i < j)
	{
	    /* this entry is in the upper triangular part */
	    is_lower = FALSE ;
	}
	if (i > j)
	{
	    /* this entry is in the lower triangular part */
	    is_upper = FALSE ;
	}

	if (xtype == CHOLMOD_REAL)
	{
	    Tx [k] = x ;
	}
	else if (xtype == CHOLMOD_COMPLEX)
	{
	    Tx [2*k  ] = x ;	/* real part */
	    Tx [2*k+1] = z ;	/* imaginary part */
	}

	if (i == 0 || j == 0)
	{
	    one_based = FALSE ;
	}

	imax = MAX (i, imax) ;
	jmax = MAX (j, jmax) ;
    }

    /* ---------------------------------------------------------------------- */
    /* convert to zero-based */
    /* ---------------------------------------------------------------------- */

    if (one_based)
    {
	/* input matrix is one-based; convert matrix to zero-based */
	for (k = 0 ; k < (Int) nnz ; k++)
	{
	    Ti [k]-- ;
	    Tj [k]-- ;
	}
    }

    if (one_based ?
	(imax >  (Int) nrow || jmax >  (Int) ncol) :
	(imax >= (Int) nrow || jmax >= (Int) ncol))
    {
	/* indices out of range */
	CHOLMOD(free_triplet) (&T, Common) ;
	ERROR (CHOLMOD_INVALID, "indices out of range") ;
	return (NULL) ;
    }

    /* ---------------------------------------------------------------------- */
    /* determine the stype, if not yet known */
    /* ---------------------------------------------------------------------- */

    if (unknown)
    {
	if (is_lower && is_upper)
	{
	    /* diagonal matrix, symmetric with upper part present */
	    stype = STYPE_SYMMETRIC_UPPER ;
	}
	else if (is_lower && !is_upper)
	{
	    /* symmetric, lower triangular part present */
	    stype = STYPE_SYMMETRIC_LOWER ;
	}
	else if (!is_lower && is_upper)
	{
	    /* symmetric, upper triangular part present */
	    stype = STYPE_SYMMETRIC_UPPER ;
	}
	else
	{
	    /* unsymmetric */
	    stype = STYPE_UNSYMMETRIC ;
	    extra = 0 ;
	}
    }

    /* ---------------------------------------------------------------------- */
    /* add the remainder of symmetric, skew-symmetric or Hermitian matrices */
    /* ---------------------------------------------------------------------- */

    /* note that this step is not done for real symmetric or complex Hermitian
     * matrices, unless prefer_unsym is TRUE */
    if (extra > 0)
    {
	p = nnz ;
	for (k = 0 ; k < (Int) nnz ; k++)
	{
	    i = Ti [k] ;
	    j = Tj [k] ;
	    if (i != j)
	    {
		Ti [p] = j ;
		Tj [p] = i ;
		if (xtype == CHOLMOD_REAL)
		{
		    if (skew_symmetric)
		    {
			Tx [p] = -Tx [k] ;
		    }
		    else
		    {
			Tx [p] =  Tx [k] ;
		    }
		}
		else if (xtype == CHOLMOD_COMPLEX)
		{
		    if (skew_symmetric)
		    {
			Tx [2*p  ] = -Tx [2*k ] ;
			Tx [2*p+1] = -Tx [2*k+1] ;
		    }
		    else if (complex_symmetric)
		    {
			Tx [2*p  ] =  Tx [2*k ] ;
			Tx [2*p+1] =  Tx [2*k+1] ;
		    }
		    else /* Hermitian */
		    {
			Tx [2*p  ] =  Tx [2*k ] ;
			Tx [2*p+1] = -Tx [2*k+1] ;
		    }
		}
		p++ ;
	    }
	}
	T->nnz = p ;
	nnz = p ;
    }

    T->stype = stype ;

    /* ---------------------------------------------------------------------- */
    /* create values for a pattern-only matrix */
    /* ---------------------------------------------------------------------- */

    if (xtype == CHOLMOD_PATTERN)
    {
	if (stype == STYPE_UNSYMMETRIC || Common->prefer_binary)
	{
	    /* unsymmetric case, or binary case */
	    for (k = 0 ; k < (Int) nnz ; k++)
	    {
		Tx [k] = 1 ;
	    }
	}
	else
	{
	    /* compute the row and columm degrees (excluding the diagonal) */
	    for (i = 0 ; i < (Int) nrow ; i++)
	    {
		Rdeg [i] = 0 ;
	    }
	    for (j = 0 ; j < (Int) ncol ; j++)
	    {
		Cdeg [j] = 0 ;
	    }
	    for (k = 0 ; k < (Int) nnz ; k++)
	    {
		i = Ti [k] ;
		j = Tj [k] ;
		if ((stype < 0 && i > j) || (stype > 0 && i < j))
		{
		    /* both a(i,j) and a(j,i) appear in the matrix */
		    Rdeg [i]++ ;
		    Cdeg [j]++ ;
		    Rdeg [j]++ ;
		    Cdeg [i]++ ;
		}
	    }
	    /* assign the numerical values */
	    for (k = 0 ; k < (Int) nnz ; k++)
	    {
		i = Ti [k] ;
		j = Tj [k] ;
		Tx [k] = (i == j) ? (1 + MAX (Rdeg [i], Cdeg [j])) : (-1) ;
	    }
	}
    }

    /* ---------------------------------------------------------------------- */
    /* return the new triplet matrix */
    /* ---------------------------------------------------------------------- */

    return (T) ;
}


/* ========================================================================== */
/* === read_dense =========================================================== */
/* ========================================================================== */

/* Header has already been read in, including first line (nrow ncol).
 * Read a dense matrix. */

static cholmod_dense *read_dense
(
    /* ---- input ---- */
    FILE *f,		    /* file to read from, must already be open */
    size_t nrow,	    /* number of rows */
    size_t ncol,	    /* number of columns */
    int stype,		    /* stype from header */
    /* ---- workspace */
    char *buf,		    /* of size MAXLINE+1 */
    /* --------------- */
    cholmod_common *Common
)
{
    double x, z ;
    double *Xx = NULL ;
    cholmod_dense *X ;
    Int nitems, xtype = -1, nshould = 0, i, j, k, kup, first ;

    /* ---------------------------------------------------------------------- */
    /* quick return for empty matrix */
    /* ---------------------------------------------------------------------- */

    if (nrow == 0 || ncol == 0)
    {
	/* return an empty dense matrix */
	return (CHOLMOD(zeros) (nrow, ncol, CHOLMOD_REAL, Common)) ;
    }

    /* ---------------------------------------------------------------------- */
    /* read the entries */
    /* ---------------------------------------------------------------------- */

    first = TRUE ;

    for (j = 0 ; j < (Int) ncol ; j++)
    {

	/* ------------------------------------------------------------------ */
	/* get the row index of the first entry in the file for column j */
	/* ------------------------------------------------------------------ */

	if (stype == STYPE_UNSYMMETRIC)
	{
	    i = 0 ;
	}
	else if (stype == STYPE_SKEW_SYMMETRIC)
	{
	    i = j+1 ;
	}
	else /* real symmetric or complex Hermitian lower */
	{
	    i = j ;
	}

	/* ------------------------------------------------------------------ */
	/* get column j */
	/* ------------------------------------------------------------------ */

	for ( ; i < (Int) nrow ; i++)
	{

	    /* -------------------------------------------------------------- */
	    /* get the next entry, skipping blank lines and comment lines */
	    /* -------------------------------------------------------------- */

	    x = 0 ;
	    z = 0 ;
	    for ( ; ; )
	    {

		if (!get_line (f, buf))
		{
		    /* premature end of file - not enough entries read in */
		    ERROR (CHOLMOD_INVALID, "premature EOF") ;
		    return (NULL) ;
		}

		if (is_blank_line (buf))
		{
		    /* blank line or comment */
		    continue ;
		}
		nitems = sscanf (buf, "%lg %lg\n", &x, &z) ;
		x = fix_inf (x) ;
		z = fix_inf (z) ;
		break ;
	    }

	    nitems = (nitems == EOF) ? 0 : nitems ;

	    /* -------------------------------------------------------------- */
	    /* for first entry: determine type and allocate dense matrix */
	    /* -------------------------------------------------------------- */

	    if (first)
	    {
		first = FALSE ;

		if (nitems < 1 || nitems > 2)
		{
		    /* invalid matrix */
		    ERROR (CHOLMOD_INVALID, "invalid format") ;
		    return (NULL) ;
		}
		else if (nitems == 1)
		{
		    /* a real matrix */
		    xtype = CHOLMOD_REAL ;
		}
		else if (nitems == 2)
		{
		    /* a complex matrix */
		    xtype = CHOLMOD_COMPLEX ;
		}

		/* the rest of the lines should have same number of entries */
		nshould = nitems ;

		/* allocate the result */
		X = CHOLMOD(zeros) (nrow, ncol, xtype, Common) ;
		if (Common->status < CHOLMOD_OK)
		{
		    /* out of memory */
		    return (NULL) ;
		}
		Xx = X->x ;
	    }

	    /* -------------------------------------------------------------- */
	    /* save the entry in the dense matrix */
	    /* -------------------------------------------------------------- */

	    if (nitems != nshould)
	    {
		/* wrong format or premature end-of-file */
		CHOLMOD(free_dense) (&X, Common) ;
		ERROR (CHOLMOD_INVALID, "invalid matrix file") ;
		return (NULL) ;
	    }

	    k = i + j*nrow ;
	    kup = j + i*nrow ;

	    if (xtype == CHOLMOD_REAL)
	    {
		/* real matrix */
		Xx [k] = x ;
		if (k != kup)
		{
		    if (stype == STYPE_SYMMETRIC_LOWER)
		    {
			/* real symmetric matrix */
			Xx [kup] = x ;
		    }
		    else if (stype == STYPE_SKEW_SYMMETRIC)
		    {
			/* real skew symmetric matrix */
			Xx [kup] = -x ;
		    }
		}
	    }
	    else if (xtype == CHOLMOD_COMPLEX)
	    {
		Xx [2*k  ] = x ;	    /* real part */
		Xx [2*k+1] = z ;	    /* imaginary part */
		if (k != kup)
		{
		    if (stype == STYPE_SYMMETRIC_LOWER)
		    {
			/* complex Hermitian */
			Xx [2*kup  ] = x ;	    /* real part */
			Xx [2*kup+1] = -z ;	    /* imaginary part */
		    }
		    else if (stype == STYPE_SKEW_SYMMETRIC)
		    {
			/* complex skew symmetric */
			Xx [2*kup  ] = -x ;	    /* real part */
			Xx [2*kup+1] = -z ;	    /* imaginary part */
		    }
		    if (stype == STYPE_COMPLEX_SYMMETRIC_LOWER)
		    {
			/* complex symmetric */
			Xx [2*kup  ] = x ;	    /* real part */
			Xx [2*kup+1] = z ;	    /* imaginary part */
		    }
		}
	    }
	}
    }

    /* ---------------------------------------------------------------------- */
    /* return the new dense matrix */
    /* ---------------------------------------------------------------------- */

    return (X) ;
}


/* ========================================================================== */
/* === cholmod_read_triplet ================================================= */
/* ========================================================================== */

/* Read in a triplet matrix from a file. */

cholmod_triplet *CHOLMOD(read_triplet)
(
    /* ---- input ---- */
    FILE *f,		/* file to read from, must already be open */
    /* --------------- */
    cholmod_common *Common
)
{
    char buf [MAXLINE+1] ;
    size_t nrow, ncol, nnz ;
    int stype, mtype ;

    /* ---------------------------------------------------------------------- */
    /* check inputs */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (NULL) ;
    RETURN_IF_NULL (f, NULL) ;
    Common->status = CHOLMOD_OK ;

    /* ---------------------------------------------------------------------- */
    /* read the header and first data line */
    /* ---------------------------------------------------------------------- */

    if (!read_header (f, buf, &mtype, &nrow, &ncol, &nnz, &stype) ||
	mtype != CHOLMOD_TRIPLET)
    {
	/* invalid matrix - this function can only read in a triplet matrix */
	ERROR (CHOLMOD_INVALID, "invalid format") ;
	return (NULL) ;
    }

    /* ---------------------------------------------------------------------- */
    /* read the triplet matrix */
    /* ---------------------------------------------------------------------- */

    return (read_triplet (f, nrow, ncol, nnz, stype, FALSE, buf, Common)) ;
}


/* ========================================================================== */
/* === cholmod_read_sparse ================================================== */
/* ========================================================================== */

/* Read a sparse matrix from a file.  See cholmod_read_triplet for a discussion
 * of the file format.
 *
 * If Common->prefer_upper is TRUE (the default case), a symmetric matrix is
 * returned stored in upper-triangular form (A->stype == 1).
 */

cholmod_sparse *CHOLMOD(read_sparse)
(
    /* ---- input ---- */
    FILE *f,		/* file to read from, must already be open */
    /* --------------- */
    cholmod_common *Common
)
{
    cholmod_sparse *A, *A2 ;
    cholmod_triplet *T ;

    /* ---------------------------------------------------------------------- */
    /* check inputs */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (NULL) ;
    RETURN_IF_NULL (f, NULL) ;
    Common->status = CHOLMOD_OK ;

    /* ---------------------------------------------------------------------- */
    /* convert to a sparse matrix in compressed-column form */
    /* ---------------------------------------------------------------------- */

    T = CHOLMOD(read_triplet) (f, Common) ;
    A = CHOLMOD(triplet_to_sparse) (T, 0, Common) ;
    CHOLMOD(free_triplet) (&T, Common) ;

    if (Common->prefer_upper && A != NULL && A->stype == -1)
    {
	/* A=A' */
	A2 = CHOLMOD(transpose) (A, 2, Common) ;
	CHOLMOD(free_sparse) (&A, Common) ;
	A = A2 ;
    }
    return (A) ;
}


/* ========================================================================== */
/* === cholmod_read_dense =================================================== */
/* ========================================================================== */

/* Read a dense matrix from a file. */

cholmod_dense *CHOLMOD(read_dense)
(
    /* ---- input ---- */
    FILE *f,		/* file to read from, must already be open */
    /* --------------- */
    cholmod_common *Common
)
{
    char buf [MAXLINE+1] ;
    size_t nrow, ncol, nnz ;
    int stype, mtype ;

    /* ---------------------------------------------------------------------- */
    /* check inputs */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (NULL) ;
    RETURN_IF_NULL (f, NULL) ;
    Common->status = CHOLMOD_OK ;

    /* ---------------------------------------------------------------------- */
    /* read the header and first data line */
    /* ---------------------------------------------------------------------- */

    if (!read_header (f, buf, &mtype, &nrow, &ncol, &nnz, &stype) ||
	mtype != CHOLMOD_DENSE)
    {
	/* invalid matrix - this function can only read in a dense matrix */
	ERROR (CHOLMOD_INVALID, "invalid format") ;
	return (NULL) ;
    }

    /* ---------------------------------------------------------------------- */
    /* read the dense matrix */
    /* ---------------------------------------------------------------------- */

    return (read_dense (f, nrow, ncol, stype, buf, Common)) ;
}


/* ========================================================================== */
/* === cholmod_read_matrix ================================================== */
/* ========================================================================== */

/* Read a triplet matrix, sparse matrix or a dense matrix from a file.  Returns
 * a void pointer to either a cholmod_triplet, cholmod_sparse, or cholmod_dense
 * object.  The type of object is passed back to the caller as the mtype
 * argument. */

void *CHOLMOD(read_matrix)
(
    /* ---- input ---- */
    FILE *f,		/* file to read from, must already be open */
    int prefer,		/* If 0, a sparse matrix is always return as a
			 *	cholmod_triplet form.  It can have any stype
			 *	(symmetric-lower, unsymmetric, or
			 *	symmetric-upper).
			 * If 1, a sparse matrix is returned as an unsymmetric
			 *	cholmod_sparse form (A->stype == 0), with both
			 *	upper and lower triangular parts present.
			 *	This is what the MATLAB mread mexFunction does,
			 *	since MATLAB does not have an stype.
			 * If 2, a sparse matrix is returned with an stype of 0
			 *	or 1 (unsymmetric, or symmetric with upper part
			 *	stored).
			 * This argument has no effect for dense matrices.
			 */
    /* ---- output---- */
    int *mtype,		/* CHOLMOD_TRIPLET, CHOLMOD_SPARSE or CHOLMOD_DENSE */
    /* --------------- */
    cholmod_common *Common
)
{
    void *G = NULL ;
    cholmod_sparse *A, *A2 ;
    cholmod_triplet *T ;
    char buf [MAXLINE+1] ;
    size_t nrow, ncol, nnz ;
    int stype ;

    /* ---------------------------------------------------------------------- */
    /* check inputs */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (NULL) ;
    RETURN_IF_NULL (f, NULL) ;
    RETURN_IF_NULL (mtype, NULL) ;
    Common->status = CHOLMOD_OK ;

    /* ---------------------------------------------------------------------- */
    /* read the header to determine the mtype */
    /* ---------------------------------------------------------------------- */

    if (!read_header (f, buf, mtype, &nrow, &ncol, &nnz, &stype))
    {
	/* invalid matrix */
	ERROR (CHOLMOD_INVALID, "invalid format") ;
	return (NULL) ;
    }

    /* ---------------------------------------------------------------------- */
    /* read a matrix */
    /* ---------------------------------------------------------------------- */

    if (*mtype == CHOLMOD_TRIPLET)
    {
	/* read in the triplet matrix, converting to unsymmetric format if
	 * prefer == 1 */
	T = read_triplet (f, nrow, ncol, nnz, stype, prefer == 1, buf, Common) ;
	if (prefer == 0)
	{
	    /* return matrix in its original triplet form */
	    G = T ;
	}
	else
	{
	    /* return matrix in a compressed-column form */
	    A = CHOLMOD(triplet_to_sparse) (T, 0, Common) ;
	    CHOLMOD(free_triplet) (&T, Common) ;
	    if (A != NULL && prefer == 2 && A->stype == -1)
	    {
		/* convert A from symmetric-lower to symmetric-upper */
		A2 = CHOLMOD(transpose) (A, 2, Common) ;
		CHOLMOD(free_sparse) (&A, Common) ;
		A = A2 ;
	    }
	    *mtype = CHOLMOD_SPARSE ;
	    G = A ;
	}
    }
    else if (*mtype == CHOLMOD_DENSE)
    {
	/* return a dense matrix */
	G = read_dense (f, nrow, ncol, stype, buf, Common) ;
    }
    return (G) ;
}
#endif
