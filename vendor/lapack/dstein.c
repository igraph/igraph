/*  -- translated by f2c (version 20240504).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#include "f2c.h"

/* Table of constant values */

static integer c__2 = 2;
static integer c__1 = 1;
static integer c_n1 = -1;

/* > \brief \b DSTEIN   

    =========== DOCUMENTATION ===========   

   Online html documentation available at   
              http://www.netlib.org/lapack/explore-html/   

   > \htmlonly   
   > Download DSTEIN + dependencies   
   > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dstein.
f">   
   > [TGZ]</a>   
   > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dstein.
f">   
   > [ZIP]</a>   
   > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dstein.
f">   
   > [TXT]</a>   
   > \endhtmlonly   

    Definition:   
    ===========   

         SUBROUTINE DSTEIN( N, D, E, M, W, IBLOCK, ISPLIT, Z, LDZ, WORK,   
                            IWORK, IFAIL, INFO )   

         INTEGER            INFO, LDZ, M, N   
         INTEGER            IBLOCK( * ), IFAIL( * ), ISPLIT( * ),   
        $                   IWORK( * )   
         DOUBLE PRECISION   D( * ), E( * ), W( * ), WORK( * ), Z( LDZ, * )   


   > \par Purpose:   
    =============   
   >   
   > \verbatim   
   >   
   > DSTEIN computes the eigenvectors of a real symmetric tridiagonal   
   > matrix T corresponding to specified eigenvalues, using inverse   
   > iteration.   
   >   
   > The maximum number of iterations allowed for each eigenvector is   
   > specified by an internal parameter MAXITS (currently set to 5).   
   > \endverbatim   

    Arguments:   
    ==========   

   > \param[in] N   
   > \verbatim   
   >          N is INTEGER   
   >          The order of the matrix.  N >= 0.   
   > \endverbatim   
   >   
   > \param[in] D   
   > \verbatim   
   >          D is DOUBLE PRECISION array, dimension (N)   
   >          The n diagonal elements of the tridiagonal matrix T.   
   > \endverbatim   
   >   
   > \param[in] E   
   > \verbatim   
   >          E is DOUBLE PRECISION array, dimension (N-1)   
   >          The (n-1) subdiagonal elements of the tridiagonal matrix   
   >          T, in elements 1 to N-1.   
   > \endverbatim   
   >   
   > \param[in] M   
   > \verbatim   
   >          M is INTEGER   
   >          The number of eigenvectors to be found.  0 <= M <= N.   
   > \endverbatim   
   >   
   > \param[in] W   
   > \verbatim   
   >          W is DOUBLE PRECISION array, dimension (N)   
   >          The first M elements of W contain the eigenvalues for   
   >          which eigenvectors are to be computed.  The eigenvalues   
   >          should be grouped by split-off block and ordered from   
   >          smallest to largest within the block.  ( The output array   
   >          W from DSTEBZ with ORDER = 'B' is expected here. )   
   > \endverbatim   
   >   
   > \param[in] IBLOCK   
   > \verbatim   
   >          IBLOCK is INTEGER array, dimension (N)   
   >          The submatrix indices associated with the corresponding   
   >          eigenvalues in W; IBLOCK(i)=1 if eigenvalue W(i) belongs to   
   >          the first submatrix from the top, =2 if W(i) belongs to   
   >          the second submatrix, etc.  ( The output array IBLOCK   
   >          from DSTEBZ is expected here. )   
   > \endverbatim   
   >   
   > \param[in] ISPLIT   
   > \verbatim   
   >          ISPLIT is INTEGER array, dimension (N)   
   >          The splitting points, at which T breaks up into submatrices.   
   >          The first submatrix consists of rows/columns 1 to   
   >          ISPLIT( 1 ), the second of rows/columns ISPLIT( 1 )+1   
   >          through ISPLIT( 2 ), etc.   
   >          ( The output array ISPLIT from DSTEBZ is expected here. )   
   > \endverbatim   
   >   
   > \param[out] Z   
   > \verbatim   
   >          Z is DOUBLE PRECISION array, dimension (LDZ, M)   
   >          The computed eigenvectors.  The eigenvector associated   
   >          with the eigenvalue W(i) is stored in the i-th column of   
   >          Z.  Any vector which fails to converge is set to its current   
   >          iterate after MAXITS iterations.   
   > \endverbatim   
   >   
   > \param[in] LDZ   
   > \verbatim   
   >          LDZ is INTEGER   
   >          The leading dimension of the array Z.  LDZ >= max(1,N).   
   > \endverbatim   
   >   
   > \param[out] WORK   
   > \verbatim   
   >          WORK is DOUBLE PRECISION array, dimension (5*N)   
   > \endverbatim   
   >   
   > \param[out] IWORK   
   > \verbatim   
   >          IWORK is INTEGER array, dimension (N)   
   > \endverbatim   
   >   
   > \param[out] IFAIL   
   > \verbatim   
   >          IFAIL is INTEGER array, dimension (M)   
   >          On normal exit, all elements of IFAIL are zero.   
   >          If one or more eigenvectors fail to converge after   
   >          MAXITS iterations, then their indices are stored in   
   >          array IFAIL.   
   > \endverbatim   
   >   
   > \param[out] INFO   
   > \verbatim   
   >          INFO is INTEGER   
   >          = 0: successful exit.   
   >          < 0: if INFO = -i, the i-th argument had an illegal value   
   >          > 0: if INFO = i, then i eigenvectors failed to converge   
   >               in MAXITS iterations.  Their indices are stored in   
   >               array IFAIL.   
   > \endverbatim   

   > \par Internal Parameters:   
    =========================   
   >   
   > \verbatim   
   >  MAXITS  INTEGER, default = 5   
   >          The maximum number of iterations performed.   
   >   
   >  EXTRA   INTEGER, default = 2   
   >          The number of iterations performed after norm growth   
   >          criterion is satisfied, should be at least 1.   
   > \endverbatim   

    Authors:   
    ========   

   > \author Univ. of Tennessee   
   > \author Univ. of California Berkeley   
   > \author Univ. of Colorado Denver   
   > \author NAG Ltd.   

   > \date November 2011   

   > \ingroup doubleOTHERcomputational   

    =====================================================================   
   Subroutine */ int igraphdstein_(integer *n, doublereal *d__, doublereal *e, 
	integer *m, doublereal *w, integer *iblock, integer *isplit, 
	doublereal *z__, integer *ldz, doublereal *work, integer *iwork, 
	integer *ifail, integer *info)
{
    /* System generated locals */
    integer z_dim1, z_offset, i__1, i__2, i__3;
    doublereal d__1, d__2, d__3, d__4, d__5;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    integer i__, j, b1, j1, bn;
    doublereal xj, scl, eps, sep, nrm, tol;
    integer its;
    doublereal xjm, ztr, eps1;
    integer jblk, nblk;
    extern doublereal igraphddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    integer jmax;
    extern doublereal igraphdnrm2_(integer *, doublereal *, integer *);
    extern /* Subroutine */ int igraphdscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    integer iseed[4], gpind, iinfo;
    extern doublereal igraphdasum_(integer *, doublereal *, integer *);
    extern /* Subroutine */ int igraphdcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), igraphdaxpy_(integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *);
    doublereal ortol;
    integer indrv1, indrv2, indrv3, indrv4, indrv5;
    extern doublereal igraphdlamch_(char *);
    extern /* Subroutine */ int igraphdlagtf_(integer *, doublereal *, doublereal *,
	     doublereal *, doublereal *, doublereal *, doublereal *, integer *
	    , integer *);
    extern integer igraphidamax_(integer *, doublereal *, integer *);
    extern /* Subroutine */ int igraphxerbla_(char *, integer *, ftnlen), igraphdlagts_(
	    integer *, integer *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *);
    integer nrmchk;
    extern /* Subroutine */ int igraphdlarnv_(integer *, integer *, integer *, 
	    doublereal *);
    integer blksiz;
    doublereal onenrm, dtpcrt, pertol;


/*  -- LAPACK computational routine (version 3.4.0) --   
    -- LAPACK is a software package provided by Univ. of Tennessee,    --   
    -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--   
       November 2011   


    =====================================================================   


       Test the input parameters.   

       Parameter adjustments */
    --d__;
    --e;
    --w;
    --iblock;
    --isplit;
    z_dim1 = *ldz;
    z_offset = 1 + z_dim1;
    z__ -= z_offset;
    --work;
    --iwork;
    --ifail;

    /* Function Body */
    *info = 0;
    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ifail[i__] = 0;
/* L10: */
    }

    if (*n < 0) {
	*info = -1;
    } else if (*m < 0 || *m > *n) {
	*info = -4;
    } else if (*ldz < max(1,*n)) {
	*info = -9;
    } else {
	i__1 = *m;
	for (j = 2; j <= i__1; ++j) {
	    if (iblock[j] < iblock[j - 1]) {
		*info = -6;
		goto L30;
	    }
	    if (iblock[j] == iblock[j - 1] && w[j] < w[j - 1]) {
		*info = -5;
		goto L30;
	    }
/* L20: */
	}
L30:
	;
    }

    if (*info != 0) {
	i__1 = -(*info);
	igraphxerbla_("DSTEIN", &i__1, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible */

    if (*n == 0 || *m == 0) {
	return 0;
    } else if (*n == 1) {
	z__[z_dim1 + 1] = 1.;
	return 0;
    }

/*     Get machine constants. */

    eps = igraphdlamch_("Precision");

/*     Initialize seed for random number generator DLARNV. */

    for (i__ = 1; i__ <= 4; ++i__) {
	iseed[i__ - 1] = 1;
/* L40: */
    }

/*     Initialize pointers. */

    indrv1 = 0;
    indrv2 = indrv1 + *n;
    indrv3 = indrv2 + *n;
    indrv4 = indrv3 + *n;
    indrv5 = indrv4 + *n;

/*     Compute eigenvectors of matrix blocks. */

    j1 = 1;
    i__1 = iblock[*m];
    for (nblk = 1; nblk <= i__1; ++nblk) {

/*        Find starting and ending indices of block nblk. */

	if (nblk == 1) {
	    b1 = 1;
	} else {
	    b1 = isplit[nblk - 1] + 1;
	}
	bn = isplit[nblk];
	blksiz = bn - b1 + 1;
	if (blksiz == 1) {
	    goto L60;
	}
	gpind = b1;

/*        Compute reorthogonalization criterion and stopping criterion. */

	onenrm = (d__1 = d__[b1], abs(d__1)) + (d__2 = e[b1], abs(d__2));
/* Computing MAX */
	d__3 = onenrm, d__4 = (d__1 = d__[bn], abs(d__1)) + (d__2 = e[bn - 1],
		 abs(d__2));
	onenrm = max(d__3,d__4);
	i__2 = bn - 1;
	for (i__ = b1 + 1; i__ <= i__2; ++i__) {
/* Computing MAX */
	    d__4 = onenrm, d__5 = (d__1 = d__[i__], abs(d__1)) + (d__2 = e[
		    i__ - 1], abs(d__2)) + (d__3 = e[i__], abs(d__3));
	    onenrm = max(d__4,d__5);
/* L50: */
	}
	ortol = onenrm * .001;

	dtpcrt = sqrt(.1 / blksiz);

/*        Loop through eigenvalues of block nblk. */

L60:
	jblk = 0;
	i__2 = *m;
	for (j = j1; j <= i__2; ++j) {
	    if (iblock[j] != nblk) {
		j1 = j;
		goto L160;
	    }
	    ++jblk;
	    xj = w[j];

/*           Skip all the work if the block size is one. */

	    if (blksiz == 1) {
		work[indrv1 + 1] = 1.;
		goto L120;
	    }

/*           If eigenvalues j and j-1 are too close, add a relatively   
             small perturbation. */

	    if (jblk > 1) {
		eps1 = (d__1 = eps * xj, abs(d__1));
		pertol = eps1 * 10.;
		sep = xj - xjm;
		if (sep < pertol) {
		    xj = xjm + pertol;
		}
	    }

	    its = 0;
	    nrmchk = 0;

/*           Get random starting vector. */

	    igraphdlarnv_(&c__2, iseed, &blksiz, &work[indrv1 + 1]);

/*           Copy the matrix T so it won't be destroyed in factorization. */

	    igraphdcopy_(&blksiz, &d__[b1], &c__1, &work[indrv4 + 1], &c__1);
	    i__3 = blksiz - 1;
	    igraphdcopy_(&i__3, &e[b1], &c__1, &work[indrv2 + 2], &c__1);
	    i__3 = blksiz - 1;
	    igraphdcopy_(&i__3, &e[b1], &c__1, &work[indrv3 + 1], &c__1);

/*           Compute LU factors with partial pivoting  ( PT = LU ) */

	    tol = 0.;
	    igraphdlagtf_(&blksiz, &work[indrv4 + 1], &xj, &work[indrv2 + 2], &work[
		    indrv3 + 1], &tol, &work[indrv5 + 1], &iwork[1], &iinfo);

/*           Update iteration count. */

L70:
	    ++its;
	    if (its > 5) {
		goto L100;
	    }

/*           Normalize and scale the righthand side vector Pb.   

   Computing MAX */
	    d__2 = eps, d__3 = (d__1 = work[indrv4 + blksiz], abs(d__1));
	    scl = blksiz * onenrm * max(d__2,d__3) / igraphdasum_(&blksiz, &work[
		    indrv1 + 1], &c__1);
	    igraphdscal_(&blksiz, &scl, &work[indrv1 + 1], &c__1);

/*           Solve the system LU = Pb. */

	    igraphdlagts_(&c_n1, &blksiz, &work[indrv4 + 1], &work[indrv2 + 2], &
		    work[indrv3 + 1], &work[indrv5 + 1], &iwork[1], &work[
		    indrv1 + 1], &tol, &iinfo);

/*           Reorthogonalize by modified Gram-Schmidt if eigenvalues are   
             close enough. */

	    if (jblk == 1) {
		goto L90;
	    }
	    if ((d__1 = xj - xjm, abs(d__1)) > ortol) {
		gpind = j;
	    }
	    if (gpind != j) {
		i__3 = j - 1;
		for (i__ = gpind; i__ <= i__3; ++i__) {
		    ztr = -igraphddot_(&blksiz, &work[indrv1 + 1], &c__1, &z__[b1 + 
			    i__ * z_dim1], &c__1);
		    igraphdaxpy_(&blksiz, &ztr, &z__[b1 + i__ * z_dim1], &c__1, &
			    work[indrv1 + 1], &c__1);
/* L80: */
		}
	    }

/*           Check the infinity norm of the iterate. */

L90:
	    jmax = igraphidamax_(&blksiz, &work[indrv1 + 1], &c__1);
	    nrm = (d__1 = work[indrv1 + jmax], abs(d__1));

/*           Continue for additional iterations after norm reaches   
             stopping criterion. */

	    if (nrm < dtpcrt) {
		goto L70;
	    }
	    ++nrmchk;
	    if (nrmchk < 3) {
		goto L70;
	    }

	    goto L110;

/*           If stopping criterion was not satisfied, update info and   
             store eigenvector number in array ifail. */

L100:
	    ++(*info);
	    ifail[*info] = j;

/*           Accept iterate as jth eigenvector. */

L110:
	    scl = 1. / igraphdnrm2_(&blksiz, &work[indrv1 + 1], &c__1);
	    jmax = igraphidamax_(&blksiz, &work[indrv1 + 1], &c__1);
	    if (work[indrv1 + jmax] < 0.) {
		scl = -scl;
	    }
	    igraphdscal_(&blksiz, &scl, &work[indrv1 + 1], &c__1);
L120:
	    i__3 = *n;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		z__[i__ + j * z_dim1] = 0.;
/* L130: */
	    }
	    i__3 = blksiz;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		z__[b1 + i__ - 1 + j * z_dim1] = work[indrv1 + i__];
/* L140: */
	    }

/*           Save the shift to check eigenvalue spacing at next   
             iteration. */

	    xjm = xj;

/* L150: */
	}
L160:
	;
    }

    return 0;

/*     End of DSTEIN */

} /* igraphdstein_ */

