/* ========================================================================== */
/* === Core/t_cholmod_triplet =============================================== */
/* ========================================================================== */

/* -----------------------------------------------------------------------------
 * CHOLMOD/Core Module.  Copyright (C) 2005-2006,
 * Univ. of Florida.  Author: Timothy A. Davis
 * The CHOLMOD/Core Module is licensed under Version 2.1 of the GNU
 * Lesser General Public License.  See lesser.txt for a text of the license.
 * CHOLMOD is also available under other licenses; contact authors for details.
 * -------------------------------------------------------------------------- */

/* Template routine for cholmod_triplet.  All xtypes supported */

#include "cholmod_template.h"

/* ========================================================================== */
/* === t_cholmod_triplet_to_sparse ========================================== */
/* ========================================================================== */

static size_t TEMPLATE (cholmod_triplet_to_sparse)
(
    /* ---- input ---- */
    cholmod_triplet *T,	/* matrix to copy */
    /* ---- in/out --- */
    cholmod_sparse *R,	/* output matrix */
    /* --------------- */
    cholmod_common *Common
)
{
    double *Rx, *Rz, *Tx, *Tz ;
    Int *Wj, *Rp, *Ri, *Rnz, *Ti, *Tj  ;
    Int i, j, p, p1, p2, pdest, pj, k, stype, nrow, ncol, nz ;
    size_t anz ;

    /* ---------------------------------------------------------------------- */
    /* get inputs */
    /* ---------------------------------------------------------------------- */

    /* Wj contains a copy of Rp on input [ */
    Wj = Common->Iwork ;	/* size MAX (nrow,ncol). (i/l/l) */

    Rp = R->p ;
    Ri = R->i ;
    Rnz = R->nz ;
    Rx = R->x ;
    Rz = R->z ;

    Ti = T->i ;
    Tj = T->j ;
    Tx = T->x ;
    Tz = T->z ;
    nz = T->nnz ;
    nrow = T->nrow ;
    ncol = T->ncol ;
    stype = SIGN (T->stype) ;

    /* ---------------------------------------------------------------------- */
    /* construct the row form */
    /* ---------------------------------------------------------------------- */

    /* if Ti is jumbled, this part dominates the run time */

    if (stype > 0)
    {
	for (k = 0 ; k < nz ; k++)
	{
	    i = Ti [k] ;
	    j = Tj [k] ;
	    if (i < j)
	    {
		/* place triplet (j,i,x) in column i of R */
		p = Wj [i]++ ;
		Ri [p] = j ;
	    }
	    else
	    {
		/* place triplet (i,j,x) in column j of R */
		p = Wj [j]++ ;
		Ri [p] = i ;
	    }
	    ASSIGN (Rx, Rz, p, Tx, Tz, k) ;
	}
    }
    else if (stype < 0)
    {
	for (k = 0 ; k < nz ; k++)
	{
	    i = Ti [k] ;
	    j = Tj [k] ;
	    if (i > j)
	    {
		/* place triplet (j,i,x) in column i of R */
		p = Wj [i]++ ;
		Ri [p] = j ;
	    }
	    else
	    {
		/* place triplet (i,j,x) in column j of R */
		p = Wj [j]++ ;
		Ri [p] = i ;
	    }
	    ASSIGN (Rx, Rz, p, Tx, Tz, k) ;
	}
    }
    else
    {
	for (k = 0 ; k < nz ; k++)
	{
	    /* place triplet (i,j,x) in column i of R */
	    p = Wj [Ti [k]]++ ;
	    Ri [p] = Tj [k] ;
	    ASSIGN (Rx, Rz, p, Tx, Tz, k) ;
	}
    }

    /* done using Wj (i/l/l) as temporary row pointers ] */

    /* ---------------------------------------------------------------------- */
    /* sum up duplicates */
    /* ---------------------------------------------------------------------- */

    /* use Wj (i/l/l) of size ncol to keep track of duplicates in each row [ */
    for (j = 0 ; j < ncol ; j++)
    {
	Wj [j] = EMPTY ;
    }

    anz = 0 ;
    for (i = 0 ; i < nrow ; i++)
    {
	p1 = Rp [i] ;
	p2 = Rp [i+1] ;
	pdest = p1 ;
	/* at this point Wj [j] < p1 holds true for all columns j, because
	 * Ri/Rx is stored in row oriented manner */
	for (p = p1 ; p < p2 ; p++)
	{
	    j = Ri [p] ;
	    pj = Wj [j] ;
	    if (pj >= p1)
	    {
		/* this column index j is already in row i at position pj;
		 * sum up the duplicate entry */
		/* Rx [pj] += Rx [p] ; */
		ASSEMBLE (Rx, Rz, pj, Rx, Rz, p) ;
	    }
	    else
	    {
		/* keep the entry and keep track in Wj [j] for case above */
		Wj [j] = pdest ;
		if (pdest != p)
		{
		    Ri [pdest] = j ;
		    ASSIGN (Rx, Rz, pdest, Rx, Rz, p) ;
		}
		pdest++ ;
	    }
	}
	Rnz [i] = pdest - p1 ;
	anz += (pdest - p1) ;
    }
    /* done using Wj to keep track of duplicate entries in each row ] */

    /* ---------------------------------------------------------------------- */
    /* return number of entries after summing up duplicates */
    /* ---------------------------------------------------------------------- */

    return (anz) ;
}

#undef PATTERN
#undef REAL
#undef COMPLEX
#undef ZOMPLEX
