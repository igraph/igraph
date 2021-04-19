/*  -- translated by f2c (version 20191129).
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

static logical c_true = TRUE_;
static integer c__1 = 1;
static doublereal c_b18 = 1.;
static doublereal c_b20 = 0.;

/* -----------------------------------------------------------------------   
   \BeginDoc   

   \Name: dneigh   

   \Description:   
    Compute the eigenvalues of the current upper Hessenberg matrix   
    and the corresponding Ritz estimates given the current residual norm.   

   \Usage:   
    call dneigh   
       ( RNORM, N, H, LDH, RITZR, RITZI, BOUNDS, Q, LDQ, WORKL, IERR )   

   \Arguments   
    RNORM   Double precision scalar.  (INPUT)   
            Residual norm corresponding to the current upper Hessenberg   
            matrix H.   

    N       Integer.  (INPUT)   
            Size of the matrix H.   

    H       Double precision N by N array.  (INPUT)   
            H contains the current upper Hessenberg matrix.   

    LDH     Integer.  (INPUT)   
            Leading dimension of H exactly as declared in the calling   
            program.   

    RITZR,  Double precision arrays of length N.  (OUTPUT)   
    RITZI   On output, RITZR(1:N) (resp. RITZI(1:N)) contains the real   
            (respectively imaginary) parts of the eigenvalues of H.   

    BOUNDS  Double precision array of length N.  (OUTPUT)   
            On output, BOUNDS contains the Ritz estimates associated with   
            the eigenvalues RITZR and RITZI.  This is equal to RNORM   
            times the last components of the eigenvectors corresponding   
            to the eigenvalues in RITZR and RITZI.   

    Q       Double precision N by N array.  (WORKSPACE)   
            Workspace needed to store the eigenvectors of H.   

    LDQ     Integer.  (INPUT)   
            Leading dimension of Q exactly as declared in the calling   
            program.   

    WORKL   Double precision work array of length N**2 + 3*N.  (WORKSPACE)   
            Private (replicated) array on each PE or array allocated on   
            the front end.  This is needed to keep the full Schur form   
            of H and also in the calculation of the eigenvectors of H.   

    IERR    Integer.  (OUTPUT)   
            Error exit flag from dlaqrb or dtrevc.   

   \EndDoc   

   -----------------------------------------------------------------------   

   \BeginLib   

   \Local variables:   
       xxxxxx  real   

   \Routines called:   
       dlaqrb  ARPACK routine to compute the real Schur form of an   
               upper Hessenberg matrix and last row of the Schur vectors.   
       second  ARPACK utility routine for timing.   
       dmout   ARPACK utility routine that prints matrices   
       dvout   ARPACK utility routine that prints vectors.   
       dlacpy  LAPACK matrix copy routine.   
       dlapy2  LAPACK routine to compute sqrt(x**2+y**2) carefully.   
       dtrevc  LAPACK routine to compute the eigenvectors of a matrix   
               in upper quasi-triangular form   
       dgemv   Level 2 BLAS routine for matrix vector multiplication.   
       dcopy   Level 1 BLAS that copies one vector to another .   
       dnrm2   Level 1 BLAS that computes the norm of a vector.   
       dscal   Level 1 BLAS that scales a vector.   


   \Author   
       Danny Sorensen               Phuong Vu   
       Richard Lehoucq              CRPC / Rice University   
       Dept. of Computational &     Houston, Texas   
       Applied Mathematics   
       Rice University   
       Houston, Texas   

   \Revision history:   
       xx/xx/92: Version ' 2.1'   

   \SCCS Information: @(#)   
   FILE: neigh.F   SID: 2.3   DATE OF SID: 4/20/96   RELEASE: 2   

   \Remarks   
       None   

   \EndLib   

   -----------------------------------------------------------------------   

   Subroutine */ int igraphdneigh_(doublereal *rnorm, integer *n, doublereal *h__, 
	integer *ldh, doublereal *ritzr, doublereal *ritzi, doublereal *
	bounds, doublereal *q, integer *ldq, doublereal *workl, integer *ierr)
{
    /* System generated locals */
    integer h_dim1, h_offset, q_dim1, q_offset, i__1;
    doublereal d__1, d__2;

    /* Local variables */
    integer i__;
    IGRAPH_F77_SAVE real t0, t1;
    doublereal vl[1], temp;
    extern doublereal igraphdnrm2_(integer *, doublereal *, integer *);
    extern /* Subroutine */ int igraphdscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    integer iconj;
    extern /* Subroutine */ int igraphdgemv_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *), igraphdmout_(integer *, 
	    integer *, integer *, doublereal *, integer *, integer *, char *, 
	    ftnlen), igraphdvout_(integer *, integer *, doublereal *, integer *, 
	    char *, ftnlen);
    extern doublereal igraphdlapy2_(doublereal *, doublereal *);
    extern /* Subroutine */ int igraphdlaqrb_(logical *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *);
    integer mneigh = 0;
    extern /* Subroutine */ int igraphsecond_(real *), igraphdlacpy_(char *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, integer *);
    integer logfil, ndigit;
    logical select[1];
    real tneigh = 0.;
    extern /* Subroutine */ int igraphdtrevc_(char *, char *, logical *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, integer *, integer *, doublereal *, integer *);
    integer msglvl;


/*     %----------------------------------------------------%   
       | Include files for debugging and timing information |   
       %----------------------------------------------------%   


       %------------------%   
       | Scalar Arguments |   
       %------------------%   


       %-----------------%   
       | Array Arguments |   
       %-----------------%   


       %------------%   
       | Parameters |   
       %------------%   


       %------------------------%   
       | Local Scalars & Arrays |   
       %------------------------%   


       %----------------------%   
       | External Subroutines |   
       %----------------------%   


       %--------------------%   
       | External Functions |   
       %--------------------%   


       %---------------------%   
       | Intrinsic Functions |   
       %---------------------%   


       %-----------------------%   
       | Executable Statements |   
       %-----------------------%   


       %-------------------------------%   
       | Initialize timing statistics  |   
       | & message level for debugging |   
       %-------------------------------%   

       Parameter adjustments */
    --workl;
    --bounds;
    --ritzi;
    --ritzr;
    h_dim1 = *ldh;
    h_offset = 1 + h_dim1;
    h__ -= h_offset;
    q_dim1 = *ldq;
    q_offset = 1 + q_dim1;
    q -= q_offset;

    /* Function Body */
    igraphsecond_(&t0);
    msglvl = mneigh;

    if (msglvl > 2) {
	igraphdmout_(&logfil, n, n, &h__[h_offset], ldh, &ndigit, "_neigh: Enterin"
		"g upper Hessenberg matrix H ", (ftnlen)43);
    }

/*     %-----------------------------------------------------------%   
       | 1. Compute the eigenvalues, the last components of the    |   
       |    corresponding Schur vectors and the full Schur form T  |   
       |    of the current upper Hessenberg matrix H.              |   
       | dlaqrb returns the full Schur form of H in WORKL(1:N**2)  |   
       | and the last components of the Schur vectors in BOUNDS.   |   
       %-----------------------------------------------------------% */

    igraphdlacpy_("All", n, n, &h__[h_offset], ldh, &workl[1], n);
    igraphdlaqrb_(&c_true, n, &c__1, n, &workl[1], n, &ritzr[1], &ritzi[1], &bounds[
	    1], ierr);
    if (*ierr != 0) {
	goto L9000;
    }

    if (msglvl > 1) {
	igraphdvout_(&logfil, n, &bounds[1], &ndigit, "_neigh: last row of the Sch"
		"ur matrix for H", (ftnlen)42);
    }

/*     %-----------------------------------------------------------%   
       | 2. Compute the eigenvectors of the full Schur form T and  |   
       |    apply the last components of the Schur vectors to get  |   
       |    the last components of the corresponding eigenvectors. |   
       | Remember that if the i-th and (i+1)-st eigenvalues are    |   
       | complex conjugate pairs, then the real & imaginary part   |   
       | of the eigenvector components are split across adjacent   |   
       | columns of Q.                                             |   
       %-----------------------------------------------------------% */

    igraphdtrevc_("R", "A", select, n, &workl[1], n, vl, n, &q[q_offset], ldq, n, n,
	     &workl[*n * *n + 1], ierr);

    if (*ierr != 0) {
	goto L9000;
    }

/*     %------------------------------------------------%   
       | Scale the returning eigenvectors so that their |   
       | euclidean norms are all one. LAPACK subroutine |   
       | dtrevc returns each eigenvector normalized so  |   
       | that the element of largest magnitude has      |   
       | magnitude 1; here the magnitude of a complex   |   
       | number (x,y) is taken to be |x| + |y|.         |   
       %------------------------------------------------% */

    iconj = 0;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if ((d__1 = ritzi[i__], abs(d__1)) <= 0.) {

/*           %----------------------%   
             | Real eigenvalue case |   
             %----------------------% */

	    temp = igraphdnrm2_(n, &q[i__ * q_dim1 + 1], &c__1);
	    d__1 = 1. / temp;
	    igraphdscal_(n, &d__1, &q[i__ * q_dim1 + 1], &c__1);
	} else {

/*           %-------------------------------------------%   
             | Complex conjugate pair case. Note that    |   
             | since the real and imaginary part of      |   
             | the eigenvector are stored in consecutive |   
             | columns, we further normalize by the      |   
             | square root of two.                       |   
             %-------------------------------------------% */

	    if (iconj == 0) {
		d__1 = igraphdnrm2_(n, &q[i__ * q_dim1 + 1], &c__1);
		d__2 = igraphdnrm2_(n, &q[(i__ + 1) * q_dim1 + 1], &c__1);
		temp = igraphdlapy2_(&d__1, &d__2);
		d__1 = 1. / temp;
		igraphdscal_(n, &d__1, &q[i__ * q_dim1 + 1], &c__1);
		d__1 = 1. / temp;
		igraphdscal_(n, &d__1, &q[(i__ + 1) * q_dim1 + 1], &c__1);
		iconj = 1;
	    } else {
		iconj = 0;
	    }
	}
/* L10: */
    }

    igraphdgemv_("T", n, n, &c_b18, &q[q_offset], ldq, &bounds[1], &c__1, &c_b20, &
	    workl[1], &c__1);

    if (msglvl > 1) {
	igraphdvout_(&logfil, n, &workl[1], &ndigit, "_neigh: Last row of the eige"
		"nvector matrix for H", (ftnlen)48);
    }

/*     %----------------------------%   
       | Compute the Ritz estimates |   
       %----------------------------% */

    iconj = 0;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if ((d__1 = ritzi[i__], abs(d__1)) <= 0.) {

/*           %----------------------%   
             | Real eigenvalue case |   
             %----------------------% */

	    bounds[i__] = *rnorm * (d__1 = workl[i__], abs(d__1));
	} else {

/*           %-------------------------------------------%   
             | Complex conjugate pair case. Note that    |   
             | since the real and imaginary part of      |   
             | the eigenvector are stored in consecutive |   
             | columns, we need to take the magnitude    |   
             | of the last components of the two vectors |   
             %-------------------------------------------% */

	    if (iconj == 0) {
		bounds[i__] = *rnorm * igraphdlapy2_(&workl[i__], &workl[i__ + 1]);
		bounds[i__ + 1] = bounds[i__];
		iconj = 1;
	    } else {
		iconj = 0;
	    }
	}
/* L20: */
    }

    if (msglvl > 2) {
	igraphdvout_(&logfil, n, &ritzr[1], &ndigit, "_neigh: Real part of the eig"
		"envalues of H", (ftnlen)41);
	igraphdvout_(&logfil, n, &ritzi[1], &ndigit, "_neigh: Imaginary part of th"
		"e eigenvalues of H", (ftnlen)46);
	igraphdvout_(&logfil, n, &bounds[1], &ndigit, "_neigh: Ritz estimates for "
		"the eigenvalues of H", (ftnlen)47);
    }

    igraphsecond_(&t1);
    tneigh += t1 - t0;

L9000:
    return 0;

/*     %---------------%   
       | End of dneigh |   
       %---------------% */

} /* igraphdneigh_ */

