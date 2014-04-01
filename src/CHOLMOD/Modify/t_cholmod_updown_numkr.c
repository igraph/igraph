/* ========================================================================== */
/* === Modify/t_cholmod_updown_numkr ======================================== */
/* ========================================================================== */

/* -----------------------------------------------------------------------------
 * CHOLMOD/Modify Module.  Copyright (C) 2005-2006,
 * Timothy A. Davis and William W. Hager.
 * The CHOLMOD/Modify Module is licensed under Version 2.0 of the GNU
 * General Public License.  See gpl.txt for a text of the license.
 * CHOLMOD is also available under other licenses; contact authors for details.
 * http://www.suitesparse.com
 * -------------------------------------------------------------------------- */

/* Supernodal numerical update/downdate of rank K = RANK, along a single path.
 * This routine operates on a simplicial factor, but operates on adjacent
 * columns of L that would fit within a single supernode.  "Adjacent" means
 * along a single path in the elimination tree; they may or may not be
 * adjacent in the matrix L.
 *
 * external defines:  NUMERIC, WDIM, RANK.
 *
 * WDIM is 1, 2, 4, or 8.  RANK can be 1 to WDIM.
 *
 * A simple method is included (#define SIMPLE).  The code works, but is slow.
 * It is meant only to illustrate what this routine is doing.
 *
 * A rank-K update proceeds along a single path, using single-column, dual-
 * column, or quad-column updates of L.  If a column j and the next column
 * in the path (its parent) do not have the same nonzero pattern, a single-
 * column update is used.  If they do, but the 3rd and 4th column from j do
 * not have the same pattern, a dual-column update is used, in which the two
 * columns are treated as if they were a single supernode of two columns.  If
 * there are 4 columns in the path that all have the same nonzero pattern, then
 * a quad-column update is used.  All three kinds of updates can be used along
 * a single path, in a single call to this function.
 *
 * Single-column update:
 *
 *	When updating a single column of L, each iteration of the for loop,
 *	below, processes four rows of W (all columns involved) and one column
 *	of L.  Suppose we have a rank-5 update, and columns 2 through 6 of W
 *	are involved.  In this case, W in this routine is a pointer to column
 *	2 of the matrix W in the caller.  W (in the caller, shown as 'W') is
 *	held in row-major order, and is 8-by-n (a dense matrix storage format),
 *	but shown below in column form to match the column of L.  Suppose there
 *	are 13 nonzero entries in column 27 of L, with row indices 27 (the
 *	diagonal, D), 28, 30, 31, 42, 43, 44, 50, 51, 67, 81, 83, and 84.  This
 *	pattern is held in Li [Lp [27] ... Lp [27 + Lnz [27] - 1], where
 *	Lnz [27] = 13.   The modification of the current column j of L is done
 *	in the following order.  A dot (.) means the entry of W is not accessed.
 *
 *	W0 points to row 27 of W, and G is a 1-by-8 temporary vector.
 *
 *		     G[0]        G[4]
 *	G	      x  x  x  x  x  .  .  .
 *
 *		      W0
 *		      |
 *		      v
 *	27      .  .  x  x  x  x  x  .   W0 points to W (27,2)
 *
 *
 *	row    'W'    W                      column j = 27
 *	|       |     |                      of L
 *	v       v     v                      |
 *	first iteration of for loop:         v
 *
 *	28      .  .  1  5  9 13 17  .       x
 *	30      .  .  2  6 10 14 18  .       x
 *	31      .  .  3  7 11 15 19  .       x
 *	42      .  .  4  8 12 16 20  .       x
 *
 *	second iteration of for loop:
 *
 *	43      .  .  1  5  9 13 17  .       x
 *	44      .  .  2  6 10 14 18  .       x
 *	50      .  .  3  7 11 15 19  .       x
 *	51      .  .  4  8 12 16 20  .       x
 *
 *	third iteration of for loop:
 *
 *	67      .  .  1  5  9 13 17  .       x
 *	81      .  .  2  6 10 14 18  .       x
 *	83      .  .  3  7 11 15 19  .       x
 *	84      .  .  4  8 12 16 20  .       x
 *
 *	If the number of offdiagonal nonzeros in column j of L is not divisible
 *	by 4, then the switch-statement does the work for the first nz % 4 rows.
 *
 * Dual-column update:
 *
 *	In this case, two columns of L that are adjacent in the path are being
 *	updated, by 1 to 8 columns of W.  Suppose columns j=27 and j=28 are
 *	adjacent columns in the path (they need not be j and j+1).  Two rows
 *	of G and W are used as coefficients during the update: (G0, G1) and
 *	(W0, W1).
 *
 *	G0	      x  x  x  x  x  .  .  .
 *	G1	      x  x  x  x  x  .  .  .
 *
 *	27      .  .  x  x  x  x  x  .   W0 points to W (27,2)
 *	28      .  .  x  x  x  x  x  .   W1 points to W (28,2)
 *
 *
 *	row    'W'    W0,W1                  column j = 27
 *	|       |     |                      of L
 *	v       v     v                      |
 *					     | |-- column j = 28 of L
 *					     v v
 *	update L (j1,j):
 *
 *	28      .  .  1  2  3  4  5  .       x -    ("-" is not stored in L)
 *
 *	cleanup iteration since length is odd:
 *
 *	30      .  .  1  2  3  4  5  .       x x
 *
 *	then each iteration does two rows of both columns of L:
 *
 *	31      .  .  1  3  5  7  9  .       x x
 *	42      .  .  2  4  6  8 10  .       x x
 *
 *	43      .  .  1  3  5  7  9  .       x x
 *	44      .  .  2  4  6  8 10  .       x x
 *
 *	50      .  .  1  3  5  7  9  .       x x
 *	51      .  .  2  4  6  8 10  .       x x
 *
 *	67      .  .  1  3  5  7  9  .       x x
 *	81      .  .  2  4  6  8 10  .       x x
 *
 *	83      .  .  1  3  5  7  9  .       x x
 *	84      .  .  2  4  6  8 10  .       x x
 *
 *	If the number of offdiagonal nonzeros in column j of L is not even,
 *	then the cleanup iteration does the work for the first row.
 *
 * Quad-column update:
 *
 *	In this case, four columns of L that are adjacent in the path are being
 *	updated, by 1 to 8 columns of W.  Suppose columns j=27, 28, 30, and 31
 *	are adjacent columns in the path (they need not be j, j+1, ...).  Four
 *	rows of G and W are used as coefficients during the update: (G0 through
 *	G3) and (W0 through W3).  j=27, j1=28, j2=30, and j3=31.
 *
 *	G0	      x  x  x  x  x  .  .  .
 *	G1	      x  x  x  x  x  .  .  .
 *	G3	      x  x  x  x  x  .  .  .
 *	G4	      x  x  x  x  x  .  .  .
 *
 *	27      .  .  x  x  x  x  x  .   W0 points to W (27,2)
 *	28      .  .  x  x  x  x  x  .   W1 points to W (28,2)
 *	30      .  .  x  x  x  x  x  .   W2 points to W (30,2)
 *	31      .  .  x  x  x  x  x  .   W3 points to W (31,2)
 *
 *
 *	row    'W'    W0,W1,..               column j = 27
 *	|       |     |                      of L
 *	v       v     v                      |
 *					     | |-- column j = 28 of L
 *					     | | |-- column j = 30 of L
 *					     | | | |-- column j =  31 of L
 *					     v v v v
 *	update L (j1,j):
 *	28      .  .  1  2  3  4  5  .       x - - -
 *
 *	update L (j2,j):
 *	30      .  .  1  2  3  4  5  .       # x - -	 (# denotes modified)
 *
 *	update L (j2,j1)
 *	30      .  .  1  2  3  4  5  .       x # - -
 *
 *	update L (j3,j)
 *	31      .  .  1  2  3  4  5  .       # x x -
 *
 *	update L (j3,j1)
 *	31      .  .  1  2  3  4  5  .       x # x -
 *
 *	update L (j3,j2)
 *	31      .  .  1  2  3  4  5  .       x x # -
 *
 *	cleanup iteration since length is odd:
 *	42      .  .  1  2  3  4  5  .       x x x x
 *
 *
 * ----- CHOLMOD v1.1.1 did the following --------------------------------------
 *	then each iteration does two rows of all four colummns of L:
 *
 *	43      .  .  1  3  5  7  9  .       x x x x
 *	44      .  .  2  4  6  8 10  .       x x x x
 *
 *	50      .  .  1  3  5  7  9  .       x x x x
 *	51      .  .  2  4  6  8 10  .       x x x x
 *
 *	67      .  .  1  3  5  7  9  .       x x x x
 *	81      .  .  2  4  6  8 10  .       x x x x
 *
 *	83      .  .  1  3  5  7  9  .       x x x x
 *	84      .  .  2  4  6  8 10  .       x x x x
 *
 * ----- CHOLMOD v1.2.0 does the following -------------------------------------
 *	then each iteration does one rows of all four colummns of L:
 *
 *	43      .  .  1  2  3  4  5  .       x x x x
 *	44      .  .  1  2  3  4  5  .       x x x x
 *	50      .  .  1  3  5  4  5  .       x x x x
 *	51      .  .  1  2  3  4  5  .       x x x x
 *	67      .  .  1  3  5  4  5  .       x x x x
 *	81      .  .  1  2  3  4  5  .       x x x x
 *	83      .  .  1  3  5  4  5  .       x x x x
 *	84      .  .  1  2  3  4  5  .       x x x x
 *
 * This file is included in t_cholmod_updown.c, only.
 * It is not compiled separately.  It contains no user-callable routines.
 *
 * workspace: Xwork (WDIM*nrow)
 */

/* ========================================================================== */
/* === loop unrolling macros ================================================ */
/* ========================================================================== */

#undef RANK1
#undef RANK2
#undef RANK3
#undef RANK4
#undef RANK5
#undef RANK6
#undef RANK7
#undef RANK8

#define RANK1(statement) statement

#if RANK < 2
#define RANK2(statement)
#else
#define RANK2(statement) statement
#endif

#if RANK < 3
#define RANK3(statement)
#else
#define RANK3(statement) statement
#endif

#if RANK < 4
#define RANK4(statement)
#else
#define RANK4(statement) statement
#endif

#if RANK < 5
#define RANK5(statement)
#else
#define RANK5(statement) statement
#endif

#if RANK < 6
#define RANK6(statement)
#else
#define RANK6(statement) statement
#endif

#if RANK < 7
#define RANK7(statement)
#else
#define RANK7(statement) statement
#endif

#if RANK < 8
#define RANK8(statement)
#else
#define RANK8(statement) statement
#endif

#define FOR_ALL_K \
    RANK1 (DO (0)) \
    RANK2 (DO (1)) \
    RANK3 (DO (2)) \
    RANK4 (DO (3)) \
    RANK5 (DO (4)) \
    RANK6 (DO (5)) \
    RANK7 (DO (6)) \
    RANK8 (DO (7))

/* ========================================================================== */
/* === alpha/gamma ========================================================== */
/* ========================================================================== */

#undef ALPHA_GAMMA

#define ALPHA_GAMMA(Dj,Alpha,Gamma,W) \
{ \
    double dj = Dj ; \
    if (update) \
    { \
	for (k = 0 ; k < RANK ; k++) \
	{ \
	    double w = W [k] ; \
	    double alpha = Alpha [k] ; \
	    double a = alpha + (w * w) / dj ; \
	    dj *= a ; \
	    Alpha [k] = a ; \
	    Gamma [k] = (- w / dj) ; \
	    dj /= alpha ; \
	} \
    } \
    else \
    { \
	for (k = 0 ; k < RANK ; k++) \
	{ \
	    double w = W [k] ; \
	    double alpha = Alpha [k] ; \
	    double a = alpha - (w * w) / dj ; \
	    dj *= a ; \
	    Alpha [k] = a ; \
	    Gamma [k] = w / dj ; \
	    dj /= alpha ; \
	} \
    } \
    Dj = ((use_dbound) ? (CHOLMOD(dbound) (dj, Common)) : (dj)) ; \
}

/* ========================================================================== */
/* === numeric update/downdate along one path =============================== */
/* ========================================================================== */

static void NUMERIC (WDIM, RANK)
(
    int update,		/* TRUE for update, FALSE for downdate */
    Int j,		/* first column in the path */
    Int e,		/* last column in the path */
    double Alpha [ ],	/* alpha, for each column of W */
    double W [ ],	/* W is an n-by-WDIM array, stored in row-major order */
    cholmod_factor *L,	/* with unit diagonal (diagonal not stored) */
    cholmod_common *Common
)
{

#ifdef SIMPLE
#define w(row,col) W [WDIM*(row) + (col)]

    /* ---------------------------------------------------------------------- */
    /* concise but slow version for illustration only */
    /* ---------------------------------------------------------------------- */

    double Gamma [WDIM] ;
    double *Lx ;
    Int *Li, *Lp, *Lnz ;
    Int p, k ;
    Int use_dbound = IS_GT_ZERO (Common->dbound) ;

    Li = L->i ;
    Lx = L->x ;
    Lp = L->p ;
    Lnz = L->nz ;

    /* walk up the etree from node j to its ancestor e */
    for ( ; j <= e ; j = (Lnz [j] > 1) ? (Li [Lp [j] + 1]) : Int_max)
    {
	/* update the diagonal entry D (j,j) with each column of W */
	ALPHA_GAMMA (Lx [Lp [j]], Alpha, Gamma, (&(w (j,0)))) ;
	/* update column j of L */
	for (p = Lp [j] + 1 ; p < Lp [j] + Lnz [j] ; p++)
	{
	    /* update row Li [p] of column j of L with each column of W */
	    Int i = Li [p] ;
	    for (k = 0 ; k < RANK ; k++)
	    {
		w (i,k) -= w (j,k) * Lx [p] ;
		Lx [p] -= Gamma [k] * w (i,k) ;
	    }
	}
	/* clear workspace W */
	for (k = 0 ; k < RANK ; k++)
	{
	    w (j,k) = 0 ;
	}
    }

#else

    /* ---------------------------------------------------------------------- */
    /* dynamic supernodal version: supernodes detected dynamically */
    /* ---------------------------------------------------------------------- */

    double G0 [RANK], G1 [RANK], G2 [RANK], G3 [RANK] ;
    double Z0 [RANK], Z1 [RANK], Z2 [RANK], Z3 [RANK] ;
    double *W0, *W1, *W2, *W3, *Lx ;
    Int *Li, *Lp, *Lnz ;
    Int j1, j2, j3, p0, p1, p2, p3, parent, lnz, pend, k ;
    Int use_dbound = IS_GT_ZERO (Common->dbound) ;

    Li = L->i ;
    Lx = L->x ;
    Lp = L->p ;
    Lnz = L->nz ;

    /* walk up the etree from node j to its ancestor e */
    for ( ; j <= e ; j = parent)
    {

	p0 = Lp [j] ;		/* col j is Li,Lx [p0 ... p0+lnz-1] */
	lnz = Lnz [j] ;

	W0 = W + WDIM * j ;	/* pointer to row j of W */
	pend = p0 + lnz ;

	/* for k = 0 to RANK-1 do: */
	#define DO(k) Z0 [k] = W0 [k] ;
	FOR_ALL_K
	#undef DO

	/* for k = 0 to RANK-1 do: */
	#define DO(k) W0 [k] = 0 ;
	FOR_ALL_K
	#undef DO

	/* update D (j,j) */
	ALPHA_GAMMA (Lx [p0], Alpha, G0, Z0) ;
	p0++ ;

	/* determine how many columns of L to update at the same time */
	parent = (lnz > 1) ? (Li [p0]) : Int_max ;
	if (parent <= e && lnz == Lnz [parent] + 1)
	{

	    /* -------------------------------------------------------------- */
	    /* node j and its parent j1 can be updated at the same time */
	    /* -------------------------------------------------------------- */

	    j1 = parent ;
	    j2 = (lnz > 2) ? (Li [p0+1]) : Int_max ;
	    j3 = (lnz > 3) ? (Li [p0+2]) : Int_max ;
	    W1 = W + WDIM * j1 ;	/* pointer to row j1 of W */
	    p1 = Lp [j1] ;

	    /* for k = 0 to RANK-1 do: */
	    #define DO(k) Z1 [k] = W1 [k] ;
	    FOR_ALL_K
	    #undef DO

	    /* for k = 0 to RANK-1 do: */
	    #define DO(k) W1 [k] = 0 ;
	    FOR_ALL_K
	    #undef DO

	    /* update L (j1,j) */
	    {
		double lx = Lx [p0] ;

		/* for k = 0 to RANK-1 do: */
		#define DO(k) \
		    Z1 [k] -= Z0 [k] * lx ; \
		    lx -= G0 [k] * Z1 [k] ;
		FOR_ALL_K
		#undef DO

		Lx [p0++] = lx ;
	    }

	    /* update D (j1,j1) */
	    ALPHA_GAMMA (Lx [p1], Alpha, G1, Z1) ;
	    p1++ ;

	    /* -------------------------------------------------------------- */
	    /* update 2 or 4 columns of L */
	    /* -------------------------------------------------------------- */

	    if ((j2 <= e) &&			/* j2 in the current path */
		(j3 <= e) &&			/* j3 in the current path */
		(lnz == Lnz [j2] + 2) &&	/* column j2 matches */
		(lnz == Lnz [j3] + 3))		/* column j3 matches */
	    {

		/* ---------------------------------------------------------- */
		/* update 4 columns of L */
		/* ---------------------------------------------------------- */

		/* p0 and p1 currently point to row j2 in cols j and j1 of L */

		parent = (lnz > 4) ? (Li [p0+2]) : Int_max ;
		W2 = W + WDIM * j2 ;	    /* pointer to row j2 of W */
		W3 = W + WDIM * j3 ;	    /* pointer to row j3 of W */
		p2 = Lp [j2] ;
		p3 = Lp [j3] ;

		/* for k = 0 to RANK-1 do: */
		#define DO(k) Z2 [k] = W2 [k] ;
		FOR_ALL_K
		#undef DO

		/* for k = 0 to RANK-1 do: */
		#define DO(k) Z3 [k] = W3 [k] ;
		FOR_ALL_K
		#undef DO

		/* for k = 0 to RANK-1 do: */
		#define DO(k) W2 [k] = 0 ;
		FOR_ALL_K
		#undef DO

		/* for k = 0 to RANK-1 do: */
		#define DO(k) W3 [k] = 0 ;
		FOR_ALL_K
		#undef DO

		/* update L (j2,j) and update L (j2,j1) */
		{
		    double lx [2] ;
		    lx [0] = Lx [p0] ;
		    lx [1] = Lx [p1] ;

		    /* for k = 0 to RANK-1 do: */
		    #define DO(k) \
		    Z2 [k] -= Z0 [k] * lx [0] ; lx [0] -= G0 [k] * Z2 [k] ; \
		    Z2 [k] -= Z1 [k] * lx [1] ; lx [1] -= G1 [k] * Z2 [k] ;
		    FOR_ALL_K
		    #undef DO

		    Lx [p0++] = lx [0] ;
		    Lx [p1++] = lx [1] ;
		}

		/* update D (j2,j2) */
		ALPHA_GAMMA (Lx [p2], Alpha, G2, Z2) ;
		p2++ ;

		/* update L (j3,j), L (j3,j1), and L (j3,j2) */
		{
		    double lx [3] ;
		    lx [0] = Lx [p0] ;
		    lx [1] = Lx [p1] ;
		    lx [2] = Lx [p2] ;

		    /* for k = 0 to RANK-1 do: */
		    #define DO(k) \
		    Z3 [k] -= Z0 [k] * lx [0] ; lx [0] -= G0 [k] * Z3 [k] ; \
		    Z3 [k] -= Z1 [k] * lx [1] ; lx [1] -= G1 [k] * Z3 [k] ; \
		    Z3 [k] -= Z2 [k] * lx [2] ; lx [2] -= G2 [k] * Z3 [k] ;
		    FOR_ALL_K
		    #undef DO

		    Lx [p0++] = lx [0] ;
		    Lx [p1++] = lx [1] ;
		    Lx [p2++] = lx [2] ;
		}

		/* update D (j3,j3) */
		ALPHA_GAMMA (Lx [p3], Alpha, G3, Z3) ;
		p3++ ;

		/* each iteration updates L (i, [j j1 j2 j3]) */
		for ( ; p0 < pend ; p0++, p1++, p2++, p3++)
		{
		    double lx [4], *w0 ;
		    lx [0] = Lx [p0] ;
		    lx [1] = Lx [p1] ;
		    lx [2] = Lx [p2] ;
		    lx [3] = Lx [p3] ;
		    w0 = W + WDIM * Li [p0] ;

		    /* for k = 0 to RANK-1 do: */
		    #define DO(k) \
		    w0 [k] -= Z0 [k] * lx [0] ; lx [0] -= G0 [k] * w0 [k] ; \
		    w0 [k] -= Z1 [k] * lx [1] ; lx [1] -= G1 [k] * w0 [k] ; \
		    w0 [k] -= Z2 [k] * lx [2] ; lx [2] -= G2 [k] * w0 [k] ; \
		    w0 [k] -= Z3 [k] * lx [3] ; lx [3] -= G3 [k] * w0 [k] ;
		    FOR_ALL_K
		    #undef DO

		    Lx [p0] = lx [0] ;
		    Lx [p1] = lx [1] ;
		    Lx [p2] = lx [2] ;
		    Lx [p3] = lx [3] ;
		}
	    }
	    else
	    {

		/* ---------------------------------------------------------- */
		/* update 2 columns of L */
		/* ---------------------------------------------------------- */

		parent = j2 ;

		/* cleanup iteration if length is odd */
		if ((lnz - 2) % 2)
		{
		    double lx [2] , *w0 ;
		    lx [0] = Lx [p0] ;
		    lx [1] = Lx [p1] ;
		    w0 = W + WDIM * Li [p0] ;

		    /* for k = 0 to RANK-1 do: */
		    #define DO(k) \
		    w0 [k] -= Z0 [k] * lx [0] ; lx [0] -= G0 [k] * w0 [k] ; \
		    w0 [k] -= Z1 [k] * lx [1] ; lx [1] -= G1 [k] * w0 [k] ;
		    FOR_ALL_K
		    #undef DO

		    Lx [p0++] = lx [0] ;
		    Lx [p1++] = lx [1] ;
		}

		for ( ; p0 < pend ; p0 += 2, p1 += 2)
		{
		    double lx [2][2], w [2], *w0, *w1 ;
		    lx [0][0] = Lx [p0  ] ;
		    lx [1][0] = Lx [p0+1] ;
		    lx [0][1] = Lx [p1  ] ;
		    lx [1][1] = Lx [p1+1] ;
		    w0 = W + WDIM * Li [p0  ] ;
		    w1 = W + WDIM * Li [p0+1] ;

		    /* for k = 0 to RANK-1 do: */
		    #define DO(k) \
		    w [0] = w0 [k] - Z0 [k] * lx [0][0] ; \
		    w [1] = w1 [k] - Z0 [k] * lx [1][0] ; \
		    lx [0][0] -= G0 [k] * w [0] ; \
		    lx [1][0] -= G0 [k] * w [1] ; \
		    w0 [k] = w [0] -= Z1 [k] * lx [0][1] ; \
		    w1 [k] = w [1] -= Z1 [k] * lx [1][1] ; \
		    lx [0][1] -= G1 [k] * w [0] ; \
		    lx [1][1] -= G1 [k] * w [1] ;
		    FOR_ALL_K
		    #undef DO

		    Lx [p0  ] = lx [0][0] ;
		    Lx [p0+1] = lx [1][0] ;
		    Lx [p1  ] = lx [0][1] ;
		    Lx [p1+1] = lx [1][1] ;
		}
	    }
	}
	else
	{

	    /* -------------------------------------------------------------- */
	    /* update one column of L */
	    /* -------------------------------------------------------------- */

	    /* cleanup iteration if length is not a multiple of 4 */
	    switch ((lnz - 1) % 4)
	    {
		case 1:
		{
		    double lx , *w0 ;
		    lx = Lx [p0] ;
		    w0 = W + WDIM * Li [p0] ;

		    /* for k = 0 to RANK-1 do: */
		    #define DO(k) \
		    w0 [k] -= Z0 [k] * lx ; lx -= G0 [k] * w0 [k] ;
		    FOR_ALL_K
		    #undef DO

		    Lx [p0++] = lx ;
		}
		break ;

		case 2:
		{
		    double lx [2], *w0, *w1 ;
		    lx [0] = Lx [p0  ] ;
		    lx [1] = Lx [p0+1] ;
		    w0 = W + WDIM * Li [p0  ] ;
		    w1 = W + WDIM * Li [p0+1] ;

		    /* for k = 0 to RANK-1 do: */
		    #define DO(k) \
		    w0 [k] -= Z0 [k] * lx [0] ; \
		    w1 [k] -= Z0 [k] * lx [1] ; \
		    lx [0] -= G0 [k] * w0 [k] ; \
		    lx [1] -= G0 [k] * w1 [k] ;
		    FOR_ALL_K
		    #undef DO

		    Lx [p0++] = lx [0] ;
		    Lx [p0++] = lx [1] ;
		}
		break ;

		case 3:
		{
		    double lx [3], *w0, *w1, *w2 ;
		    lx [0] = Lx [p0  ] ;
		    lx [1] = Lx [p0+1] ;
		    lx [2] = Lx [p0+2] ;
		    w0 = W + WDIM * Li [p0  ] ;
		    w1 = W + WDIM * Li [p0+1] ;
		    w2 = W + WDIM * Li [p0+2] ;

		    /* for k = 0 to RANK-1 do: */
		    #define DO(k) \
		    w0 [k] -= Z0 [k] * lx [0] ; \
		    w1 [k] -= Z0 [k] * lx [1] ; \
		    w2 [k] -= Z0 [k] * lx [2] ; \
		    lx [0] -= G0 [k] * w0 [k] ; \
		    lx [1] -= G0 [k] * w1 [k] ; \
		    lx [2] -= G0 [k] * w2 [k] ;
		    FOR_ALL_K
		    #undef DO

		    Lx [p0++] = lx [0] ;
		    Lx [p0++] = lx [1] ;
		    Lx [p0++] = lx [2] ;
		}
	    }

	    for ( ; p0 < pend ; p0 += 4)
	    {
		double lx [4], *w0, *w1, *w2, *w3 ;
		lx [0] = Lx [p0  ] ;
		lx [1] = Lx [p0+1] ;
		lx [2] = Lx [p0+2] ;
		lx [3] = Lx [p0+3] ;
		w0 = W + WDIM * Li [p0  ] ;
		w1 = W + WDIM * Li [p0+1] ;
		w2 = W + WDIM * Li [p0+2] ;
		w3 = W + WDIM * Li [p0+3] ;

		/* for k = 0 to RANK-1 do: */
		#define DO(k) \
		w0 [k] -= Z0 [k] * lx [0] ; \
		w1 [k] -= Z0 [k] * lx [1] ; \
		w2 [k] -= Z0 [k] * lx [2] ; \
		w3 [k] -= Z0 [k] * lx [3] ; \
		lx [0] -= G0 [k] * w0 [k] ; \
		lx [1] -= G0 [k] * w1 [k] ; \
		lx [2] -= G0 [k] * w2 [k] ; \
		lx [3] -= G0 [k] * w3 [k] ;
		FOR_ALL_K
		#undef DO

		Lx [p0  ] = lx [0] ;
		Lx [p0+1] = lx [1] ;
		Lx [p0+2] = lx [2] ;
		Lx [p0+3] = lx [3] ;
	    }
	}
    }
#endif
}
/* prepare this file for another inclusion in t_cholmod_updown.c: */
#undef RANK
