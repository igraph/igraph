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

static doublereal c_b3 = .66666666666666663;
static integer c__1 = 1;
static integer c__0 = 0;
static integer c__4 = 4;
static logical c_true = TRUE_;
static integer c__2 = 2;

/* \BeginDoc   

   \Name: dnaup2   

   \Description:   
    Intermediate level interface called by dnaupd.   

   \Usage:   
    call dnaup2   
       ( IDO, BMAT, N, WHICH, NEV, NP, TOL, RESID, MODE, IUPD,   
         ISHIFT, MXITER, V, LDV, H, LDH, RITZR, RITZI, BOUNDS,   
         Q, LDQ, WORKL, IPNTR, WORKD, INFO )   

   \Arguments   

    IDO, BMAT, N, WHICH, NEV, TOL, RESID: same as defined in dnaupd.   
    MODE, ISHIFT, MXITER: see the definition of IPARAM in dnaupd.   

    NP      Integer.  (INPUT/OUTPUT)   
            Contains the number of implicit shifts to apply during   
            each Arnoldi iteration.   
            If ISHIFT=1, NP is adjusted dynamically at each iteration   
            to accelerate convergence and prevent stagnation.   
            This is also roughly equal to the number of matrix-vector   
            products (involving the operator OP) per Arnoldi iteration.   
            The logic for adjusting is contained within the current   
            subroutine.   
            If ISHIFT=0, NP is the number of shifts the user needs   
            to provide via reverse comunication. 0 < NP < NCV-NEV.   
            NP may be less than NCV-NEV for two reasons. The first, is   
            to keep complex conjugate pairs of "wanted" Ritz values   
            together. The second, is that a leading block of the current   
            upper Hessenberg matrix has split off and contains "unwanted"   
            Ritz values.   
            Upon termination of the IRA iteration, NP contains the number   
            of "converged" wanted Ritz values.   

    IUPD    Integer.  (INPUT)   
            IUPD .EQ. 0: use explicit restart instead implicit update.   
            IUPD .NE. 0: use implicit update.   

    V       Double precision N by (NEV+NP) array.  (INPUT/OUTPUT)   
            The Arnoldi basis vectors are returned in the first NEV   
            columns of V.   

    LDV     Integer.  (INPUT)   
            Leading dimension of V exactly as declared in the calling   
            program.   

    H       Double precision (NEV+NP) by (NEV+NP) array.  (OUTPUT)   
            H is used to store the generated upper Hessenberg matrix   

    LDH     Integer.  (INPUT)   
            Leading dimension of H exactly as declared in the calling   
            program.   

    RITZR,  Double precision arrays of length NEV+NP.  (OUTPUT)   
    RITZI   RITZR(1:NEV) (resp. RITZI(1:NEV)) contains the real (resp.   
            imaginary) part of the computed Ritz values of OP.   

    BOUNDS  Double precision array of length NEV+NP.  (OUTPUT)   
            BOUNDS(1:NEV) contain the error bounds corresponding to   
            the computed Ritz values.   

    Q       Double precision (NEV+NP) by (NEV+NP) array.  (WORKSPACE)   
            Private (replicated) work array used to accumulate the   
            rotation in the shift application step.   

    LDQ     Integer.  (INPUT)   
            Leading dimension of Q exactly as declared in the calling   
            program.   

    WORKL   Double precision work array of length at least   
            (NEV+NP)**2 + 3*(NEV+NP).  (INPUT/WORKSPACE)   
            Private (replicated) array on each PE or array allocated on   
            the front end.  It is used in shifts calculation, shifts   
            application and convergence checking.   

            On exit, the last 3*(NEV+NP) locations of WORKL contain   
            the Ritz values (real,imaginary) and associated Ritz   
            estimates of the current Hessenberg matrix.  They are   
            listed in the same order as returned from dneigh.   

            If ISHIFT .EQ. O and IDO .EQ. 3, the first 2*NP locations   
            of WORKL are used in reverse communication to hold the user   
            supplied shifts.   

    IPNTR   Integer array of length 3.  (OUTPUT)   
            Pointer to mark the starting locations in the WORKD for   
            vectors used by the Arnoldi iteration.   
            -------------------------------------------------------------   
            IPNTR(1): pointer to the current operand vector X.   
            IPNTR(2): pointer to the current result vector Y.   
            IPNTR(3): pointer to the vector B * X when used in the   
                      shift-and-invert mode.  X is the current operand.   
            -------------------------------------------------------------   

    WORKD   Double precision work array of length 3*N.  (WORKSPACE)   
            Distributed array to be used in the basic Arnoldi iteration   
            for reverse communication.  The user should not use WORKD   
            as temporary workspace during the iteration !!!!!!!!!!   
            See Data Distribution Note in DNAUPD.   

    INFO    Integer.  (INPUT/OUTPUT)   
            If INFO .EQ. 0, a randomly initial residual vector is used.   
            If INFO .NE. 0, RESID contains the initial residual vector,   
                            possibly from a previous run.   
            Error flag on output.   
            =     0: Normal return.   
            =     1: Maximum number of iterations taken.   
                     All possible eigenvalues of OP has been found.   
                     NP returns the number of converged Ritz values.   
            =     2: No shifts could be applied.   
            =    -8: Error return from LAPACK eigenvalue calculation;   
                     This should never happen.   
            =    -9: Starting vector is zero.   
            = -9999: Could not build an Arnoldi factorization.   
                     Size that was built in returned in NP.   

   \EndDoc   

   -----------------------------------------------------------------------   

   \BeginLib   

   \Local variables:   
       xxxxxx  real   

   \References:   
    1. D.C. Sorensen, "Implicit Application of Polynomial Filters in   
       a k-Step Arnoldi Method", SIAM J. Matr. Anal. Apps., 13 (1992),   
       pp 357-385.   
    2. R.B. Lehoucq, "Analysis and Implementation of an Implicitly   
       Restarted Arnoldi Iteration", Rice University Technical Report   
       TR95-13, Department of Computational and Applied Mathematics.   

   \Routines called:   
       dgetv0  ARPACK initial vector generation routine.   
       dnaitr  ARPACK Arnoldi factorization routine.   
       dnapps  ARPACK application of implicit shifts routine.   
       dnconv  ARPACK convergence of Ritz values routine.   
       dneigh  ARPACK compute Ritz values and error bounds routine.   
       dngets  ARPACK reorder Ritz values and error bounds routine.   
       dsortc  ARPACK sorting routine.   
       ivout   ARPACK utility routine that prints integers.   
       second  ARPACK utility routine for timing.   
       dmout   ARPACK utility routine that prints matrices   
       dvout   ARPACK utility routine that prints vectors.   
       dlamch  LAPACK routine that determines machine constants.   
       dlapy2  LAPACK routine to compute sqrt(x**2+y**2) carefully.   
       dcopy   Level 1 BLAS that copies one vector to another .   
       ddot    Level 1 BLAS that computes the scalar product of two vectors.   
       dnrm2   Level 1 BLAS that computes the norm of a vector.   
       dswap   Level 1 BLAS that swaps two vectors.   

   \Author   
       Danny Sorensen               Phuong Vu   
       Richard Lehoucq              CRPC / Rice University   
       Dept. of Computational &     Houston, Texas   
       Applied Mathematics   
       Rice University   
       Houston, Texas   

   \SCCS Information: @(#)   
   FILE: naup2.F   SID: 2.4   DATE OF SID: 7/30/96   RELEASE: 2   

   \Remarks   
       1. None   

   \EndLib   

   -----------------------------------------------------------------------   

   Subroutine */ int igraphdnaup2_(integer *ido, char *bmat, integer *n, char *
	which, integer *nev, integer *np, doublereal *tol, doublereal *resid, 
	integer *mode, integer *iupd, integer *ishift, integer *mxiter, 
	doublereal *v, integer *ldv, doublereal *h__, integer *ldh, 
	doublereal *ritzr, doublereal *ritzi, doublereal *bounds, doublereal *
	q, integer *ldq, doublereal *workl, integer *ipntr, doublereal *workd,
	 integer *info)
{
    /* System generated locals */
    integer h_dim1, h_offset, q_dim1, q_offset, v_dim1, v_offset, i__1, i__2;
    doublereal d__1, d__2;

    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *);
    integer s_cmp(char *, char *, ftnlen, ftnlen);
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    double sqrt(doublereal);

    /* Local variables */
    IGRAPH_F77_SAVE integer j;
    IGRAPH_F77_SAVE real t0, t1, t2, t3;
    IGRAPH_F77_SAVE integer kp[4], np0, nbx, nev0;
    extern doublereal igraphddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    IGRAPH_F77_SAVE doublereal eps23;
    IGRAPH_F77_SAVE integer ierr, iter;
    IGRAPH_F77_SAVE doublereal temp;
    extern doublereal igraphdnrm2_(integer *, doublereal *, integer *);
    IGRAPH_F77_SAVE logical getv0, cnorm;
    extern /* Subroutine */ int igraphdcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    IGRAPH_F77_SAVE integer nconv;
    extern /* Subroutine */ int igraphdmout_(integer *, integer *, integer *, 
	    doublereal *, integer *, integer *, char *, ftnlen);
    IGRAPH_F77_SAVE logical initv;
    IGRAPH_F77_SAVE doublereal rnorm;
    IGRAPH_F77_SAVE real tmvbx;
    extern /* Subroutine */ int igraphdvout_(integer *, integer *, doublereal *, 
	    integer *, char *, ftnlen), igraphivout_(integer *, integer *, integer *
	    , integer *, char *, ftnlen), igraphdgetv0_(integer *, char *, integer *
	    , logical *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *);
    extern doublereal igraphdlapy2_(doublereal *, doublereal *);
    IGRAPH_F77_SAVE integer mnaup2;
    IGRAPH_F77_SAVE real tnaup2;
    extern doublereal igraphdlamch_(char *);
    extern /* Subroutine */ int igraphdneigh_(doublereal *, integer *, doublereal *,
	     integer *, doublereal *, doublereal *, doublereal *, doublereal *
	    , integer *, doublereal *, integer *);
    IGRAPH_F77_SAVE integer nevbef;
    extern /* Subroutine */ int igraphsecond_(real *);
    IGRAPH_F77_SAVE integer logfil, ndigit;
    extern /* Subroutine */ int igraphdnaitr_(integer *, char *, integer *, integer 
	    *, integer *, integer *, doublereal *, doublereal *, doublereal *,
	     integer *, doublereal *, integer *, integer *, doublereal *, 
	    integer *);
    IGRAPH_F77_SAVE logical update;
    extern /* Subroutine */ int igraphdngets_(integer *, char *, integer *, integer 
	    *, doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *), igraphdnapps_(integer *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, doublereal *,
	     integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    doublereal *), igraphdnconv_(integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *), igraphdsortc_(char *, logical *,
	     integer *, doublereal *, doublereal *, doublereal *);
    IGRAPH_F77_SAVE logical ushift;
    IGRAPH_F77_SAVE char wprime[2];
    IGRAPH_F77_SAVE integer msglvl, nptemp, numcnv, kplusp;


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


       %---------------%   
       | Local Scalars |   
       %---------------%   


       %-----------------------%   
       | Local array arguments |   
       %-----------------------%   


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

       Parameter adjustments */
    --workd;
    --resid;
    --workl;
    --bounds;
    --ritzi;
    --ritzr;
    v_dim1 = *ldv;
    v_offset = 1 + v_dim1;
    v -= v_offset;
    h_dim1 = *ldh;
    h_offset = 1 + h_dim1;
    h__ -= h_offset;
    q_dim1 = *ldq;
    q_offset = 1 + q_dim1;
    q -= q_offset;
    --ipntr;

    /* Function Body */
    if (*ido == 0) {

	igraphsecond_(&t0);

	msglvl = mnaup2;

/*        %-------------------------------------%   
          | Get the machine dependent constant. |   
          %-------------------------------------% */

	eps23 = igraphdlamch_("Epsilon-Machine");
	eps23 = pow_dd(&eps23, &c_b3);

	nev0 = *nev;
	np0 = *np;

/*        %-------------------------------------%   
          | kplusp is the bound on the largest  |   
          |        Lanczos factorization built. |   
          | nconv is the current number of      |   
          |        "converged" eigenvlues.      |   
          | iter is the counter on the current  |   
          |      iteration step.                |   
          %-------------------------------------% */

	kplusp = *nev + *np;
	nconv = 0;
	iter = 0;

/*        %---------------------------------------%   
          | Set flags for computing the first NEV |   
          | steps of the Arnoldi factorization.   |   
          %---------------------------------------% */

	getv0 = TRUE_;
	update = FALSE_;
	ushift = FALSE_;
	cnorm = FALSE_;

	if (*info != 0) {

/*           %--------------------------------------------%   
             | User provides the initial residual vector. |   
             %--------------------------------------------% */

	    initv = TRUE_;
	    *info = 0;
	} else {
	    initv = FALSE_;
	}
    }

/*     %---------------------------------------------%   
       | Get a possibly random starting vector and   |   
       | force it into the range of the operator OP. |   
       %---------------------------------------------%   

   L10: */

    if (getv0) {
	igraphdgetv0_(ido, bmat, &c__1, &initv, n, &c__1, &v[v_offset], ldv, &resid[
		1], &rnorm, &ipntr[1], &workd[1], info);

	if (*ido != 99) {
	    goto L9000;
	}

	if (rnorm == 0.) {

/*           %-----------------------------------------%   
             | The initial vector is zero. Error exit. |   
             %-----------------------------------------% */

	    *info = -9;
	    goto L1100;
	}
	getv0 = FALSE_;
	*ido = 0;
    }

/*     %-----------------------------------%   
       | Back from reverse communication : |   
       | continue with update step         |   
       %-----------------------------------% */

    if (update) {
	goto L20;
    }

/*     %-------------------------------------------%   
       | Back from computing user specified shifts |   
       %-------------------------------------------% */

    if (ushift) {
	goto L50;
    }

/*     %-------------------------------------%   
       | Back from computing residual norm   |   
       | at the end of the current iteration |   
       %-------------------------------------% */

    if (cnorm) {
	goto L100;
    }

/*     %----------------------------------------------------------%   
       | Compute the first NEV steps of the Arnoldi factorization |   
       %----------------------------------------------------------% */

    igraphdnaitr_(ido, bmat, n, &c__0, nev, mode, &resid[1], &rnorm, &v[v_offset], 
	    ldv, &h__[h_offset], ldh, &ipntr[1], &workd[1], info);

/*     %---------------------------------------------------%   
       | ido .ne. 99 implies use of reverse communication  |   
       | to compute operations involving OP and possibly B |   
       %---------------------------------------------------% */

    if (*ido != 99) {
	goto L9000;
    }

    if (*info > 0) {
	*np = *info;
	*mxiter = iter;
	*info = -9999;
	goto L1200;
    }

/*     %--------------------------------------------------------------%   
       |                                                              |   
       |           M A I N  ARNOLDI  I T E R A T I O N  L O O P       |   
       |           Each iteration implicitly restarts the Arnoldi     |   
       |           factorization in place.                            |   
       |                                                              |   
       %--------------------------------------------------------------% */

L1000:

    ++iter;

    if (msglvl > 0) {
	igraphivout_(&logfil, &c__1, &iter, &ndigit, "_naup2: **** Start of major "
		"iteration number ****", (ftnlen)49);
    }

/*        %-----------------------------------------------------------%   
          | Compute NP additional steps of the Arnoldi factorization. |   
          | Adjust NP since NEV might have been updated by last call  |   
          | to the shift application routine dnapps.                  |   
          %-----------------------------------------------------------% */

    *np = kplusp - *nev;

    if (msglvl > 1) {
	igraphivout_(&logfil, &c__1, nev, &ndigit, "_naup2: The length of the curr"
		"ent Arnoldi factorization", (ftnlen)55);
	igraphivout_(&logfil, &c__1, np, &ndigit, "_naup2: Extend the Arnoldi fact"
		"orization by", (ftnlen)43);
    }

/*        %-----------------------------------------------------------%   
          | Compute NP additional steps of the Arnoldi factorization. |   
          %-----------------------------------------------------------% */

    *ido = 0;
L20:
    update = TRUE_;

    igraphdnaitr_(ido, bmat, n, nev, np, mode, &resid[1], &rnorm, &v[v_offset], ldv,
	     &h__[h_offset], ldh, &ipntr[1], &workd[1], info);

/*        %---------------------------------------------------%   
          | ido .ne. 99 implies use of reverse communication  |   
          | to compute operations involving OP and possibly B |   
          %---------------------------------------------------% */

    if (*ido != 99) {
	goto L9000;
    }

    if (*info > 0) {
	*np = *info;
	*mxiter = iter;
	*info = -9999;
	goto L1200;
    }
    update = FALSE_;

    if (msglvl > 1) {
	igraphdvout_(&logfil, &c__1, &rnorm, &ndigit, "_naup2: Corresponding B-nor"
		"m of the residual", (ftnlen)44);
    }

/*        %--------------------------------------------------------%   
          | Compute the eigenvalues and corresponding error bounds |   
          | of the current upper Hessenberg matrix.                |   
          %--------------------------------------------------------% */

    igraphdneigh_(&rnorm, &kplusp, &h__[h_offset], ldh, &ritzr[1], &ritzi[1], &
	    bounds[1], &q[q_offset], ldq, &workl[1], &ierr);

    if (ierr != 0) {
	*info = -8;
	goto L1200;
    }

/*        %----------------------------------------------------%   
          | Make a copy of eigenvalues and corresponding error |   
          | bounds obtained from dneigh.                       |   
          %----------------------------------------------------%   

   Computing 2nd power */
    i__1 = kplusp;
    igraphdcopy_(&kplusp, &ritzr[1], &c__1, &workl[i__1 * i__1 + 1], &c__1);
/* Computing 2nd power */
    i__1 = kplusp;
    igraphdcopy_(&kplusp, &ritzi[1], &c__1, &workl[i__1 * i__1 + kplusp + 1], &c__1)
	    ;
/* Computing 2nd power */
    i__1 = kplusp;
    igraphdcopy_(&kplusp, &bounds[1], &c__1, &workl[i__1 * i__1 + (kplusp << 1) + 1]
	    , &c__1);

/*        %---------------------------------------------------%   
          | Select the wanted Ritz values and their bounds    |   
          | to be used in the convergence test.               |   
          | The wanted part of the spectrum and corresponding |   
          | error bounds are in the last NEV loc. of RITZR,   |   
          | RITZI and BOUNDS respectively. The variables NEV  |   
          | and NP may be updated if the NEV-th wanted Ritz   |   
          | value has a non zero imaginary part. In this case |   
          | NEV is increased by one and NP decreased by one.  |   
          | NOTE: The last two arguments of dngets are no     |   
          | longer used as of version 2.1.                    |   
          %---------------------------------------------------% */

    *nev = nev0;
    *np = np0;
    numcnv = *nev;
    igraphdngets_(ishift, which, nev, np, &ritzr[1], &ritzi[1], &bounds[1], &workl[
	    1], &workl[*np + 1]);
    if (*nev == nev0 + 1) {
	numcnv = nev0 + 1;
    }

/*        %-------------------%   
          | Convergence test. |   
          %-------------------% */

    igraphdcopy_(nev, &bounds[*np + 1], &c__1, &workl[(*np << 1) + 1], &c__1);
    igraphdnconv_(nev, &ritzr[*np + 1], &ritzi[*np + 1], &workl[(*np << 1) + 1], 
	    tol, &nconv);

    if (msglvl > 2) {
	kp[0] = *nev;
	kp[1] = *np;
	kp[2] = numcnv;
	kp[3] = nconv;
	igraphivout_(&logfil, &c__4, kp, &ndigit, "_naup2: NEV, NP, NUMCNV, NCONV "
		"are", (ftnlen)34);
	igraphdvout_(&logfil, &kplusp, &ritzr[1], &ndigit, "_naup2: Real part of t"
		"he eigenvalues of H", (ftnlen)41);
	igraphdvout_(&logfil, &kplusp, &ritzi[1], &ndigit, "_naup2: Imaginary part"
		" of the eigenvalues of H", (ftnlen)46);
	igraphdvout_(&logfil, &kplusp, &bounds[1], &ndigit, "_naup2: Ritz estimate"
		"s of the current NCV Ritz values", (ftnlen)53);
    }

/*        %---------------------------------------------------------%   
          | Count the number of unwanted Ritz values that have zero |   
          | Ritz estimates. If any Ritz estimates are equal to zero |   
          | then a leading block of H of order equal to at least    |   
          | the number of Ritz values with zero Ritz estimates has  |   
          | split off. None of these Ritz values may be removed by  |   
          | shifting. Decrease NP the number of shifts to apply. If |   
          | no shifts may be applied, then prepare to exit          |   
          %---------------------------------------------------------% */

    nptemp = *np;
    i__1 = nptemp;
    for (j = 1; j <= i__1; ++j) {
	if (bounds[j] == 0.) {
	    --(*np);
	    ++(*nev);
	}
/* L30: */
    }

    if (nconv >= numcnv || iter > *mxiter || *np == 0) {

	if (msglvl > 4) {
/* Computing 2nd power */
	    i__1 = kplusp;
	    igraphdvout_(&logfil, &kplusp, &workl[i__1 * i__1 + 1], &ndigit, "_nau"
		    "p2: Real part of the eig computed by _neigh:", (ftnlen)48)
		    ;
/* Computing 2nd power */
	    i__1 = kplusp;
	    igraphdvout_(&logfil, &kplusp, &workl[i__1 * i__1 + kplusp + 1], &
		    ndigit, "_naup2: Imag part of the eig computed by _neigh:"
		    , (ftnlen)48);
/* Computing 2nd power */
	    i__1 = kplusp;
	    igraphdvout_(&logfil, &kplusp, &workl[i__1 * i__1 + (kplusp << 1) + 1], 
		    &ndigit, "_naup2: Ritz eistmates computed by _neigh:", (
		    ftnlen)42);
	}

/*           %------------------------------------------------%   
             | Prepare to exit. Put the converged Ritz values |   
             | and corresponding bounds in RITZ(1:NCONV) and  |   
             | BOUNDS(1:NCONV) respectively. Then sort. Be    |   
             | careful when NCONV > NP                        |   
             %------------------------------------------------%   

             %------------------------------------------%   
             |  Use h( 3,1 ) as storage to communicate  |   
             |  rnorm to _neupd if needed               |   
             %------------------------------------------% */
	h__[h_dim1 + 3] = rnorm;

/*           %----------------------------------------------%   
             | To be consistent with dngets, we first do a  |   
             | pre-processing sort in order to keep complex |   
             | conjugate pairs together.  This is similar   |   
             | to the pre-processing sort used in dngets    |   
             | except that the sort is done in the opposite |   
             | order.                                       |   
             %----------------------------------------------% */

	if (s_cmp(which, "LM", (ftnlen)2, (ftnlen)2) == 0) {
	    s_copy(wprime, "SR", (ftnlen)2, (ftnlen)2);
	}
	if (s_cmp(which, "SM", (ftnlen)2, (ftnlen)2) == 0) {
	    s_copy(wprime, "LR", (ftnlen)2, (ftnlen)2);
	}
	if (s_cmp(which, "LR", (ftnlen)2, (ftnlen)2) == 0) {
	    s_copy(wprime, "SM", (ftnlen)2, (ftnlen)2);
	}
	if (s_cmp(which, "SR", (ftnlen)2, (ftnlen)2) == 0) {
	    s_copy(wprime, "LM", (ftnlen)2, (ftnlen)2);
	}
	if (s_cmp(which, "LI", (ftnlen)2, (ftnlen)2) == 0) {
	    s_copy(wprime, "SM", (ftnlen)2, (ftnlen)2);
	}
	if (s_cmp(which, "SI", (ftnlen)2, (ftnlen)2) == 0) {
	    s_copy(wprime, "LM", (ftnlen)2, (ftnlen)2);
	}

	igraphdsortc_(wprime, &c_true, &kplusp, &ritzr[1], &ritzi[1], &bounds[1]);

/*           %----------------------------------------------%   
             | Now sort Ritz values so that converged Ritz  |   
             | values appear within the first NEV locations |   
             | of ritzr, ritzi and bounds, and the most     |   
             | desired one appears at the front.            |   
             %----------------------------------------------% */

	if (s_cmp(which, "LM", (ftnlen)2, (ftnlen)2) == 0) {
	    s_copy(wprime, "SM", (ftnlen)2, (ftnlen)2);
	}
	if (s_cmp(which, "SM", (ftnlen)2, (ftnlen)2) == 0) {
	    s_copy(wprime, "LM", (ftnlen)2, (ftnlen)2);
	}
	if (s_cmp(which, "LR", (ftnlen)2, (ftnlen)2) == 0) {
	    s_copy(wprime, "SR", (ftnlen)2, (ftnlen)2);
	}
	if (s_cmp(which, "SR", (ftnlen)2, (ftnlen)2) == 0) {
	    s_copy(wprime, "LR", (ftnlen)2, (ftnlen)2);
	}
	if (s_cmp(which, "LI", (ftnlen)2, (ftnlen)2) == 0) {
	    s_copy(wprime, "SI", (ftnlen)2, (ftnlen)2);
	}
	if (s_cmp(which, "SI", (ftnlen)2, (ftnlen)2) == 0) {
	    s_copy(wprime, "LI", (ftnlen)2, (ftnlen)2);
	}

	igraphdsortc_(wprime, &c_true, &kplusp, &ritzr[1], &ritzi[1], &bounds[1]);

/*           %--------------------------------------------------%   
             | Scale the Ritz estimate of each Ritz value       |   
             | by 1 / max(eps23,magnitude of the Ritz value).   |   
             %--------------------------------------------------% */

	i__1 = nev0;
	for (j = 1; j <= i__1; ++j) {
/* Computing MAX */
	    d__1 = eps23, d__2 = igraphdlapy2_(&ritzr[j], &ritzi[j]);
	    temp = max(d__1,d__2);
	    bounds[j] /= temp;
/* L35: */
	}

/*           %----------------------------------------------------%   
             | Sort the Ritz values according to the scaled Ritz  |   
             | esitmates.  This will push all the converged ones  |   
             | towards the front of ritzr, ritzi, bounds          |   
             | (in the case when NCONV < NEV.)                    |   
             %----------------------------------------------------% */

	s_copy(wprime, "LR", (ftnlen)2, (ftnlen)2);
	igraphdsortc_(wprime, &c_true, &nev0, &bounds[1], &ritzr[1], &ritzi[1]);

/*           %----------------------------------------------%   
             | Scale the Ritz estimate back to its original |   
             | value.                                       |   
             %----------------------------------------------% */

	i__1 = nev0;
	for (j = 1; j <= i__1; ++j) {
/* Computing MAX */
	    d__1 = eps23, d__2 = igraphdlapy2_(&ritzr[j], &ritzi[j]);
	    temp = max(d__1,d__2);
	    bounds[j] *= temp;
/* L40: */
	}

/*           %------------------------------------------------%   
             | Sort the converged Ritz values again so that   |   
             | the "threshold" value appears at the front of  |   
             | ritzr, ritzi and bound.                        |   
             %------------------------------------------------% */

	igraphdsortc_(which, &c_true, &nconv, &ritzr[1], &ritzi[1], &bounds[1]);

	if (msglvl > 1) {
	    igraphdvout_(&logfil, &kplusp, &ritzr[1], &ndigit, "_naup2: Sorted rea"
		    "l part of the eigenvalues", (ftnlen)43);
	    igraphdvout_(&logfil, &kplusp, &ritzi[1], &ndigit, "_naup2: Sorted ima"
		    "ginary part of the eigenvalues", (ftnlen)48);
	    igraphdvout_(&logfil, &kplusp, &bounds[1], &ndigit, "_naup2: Sorted ri"
		    "tz estimates.", (ftnlen)30);
	}

/*           %------------------------------------%   
             | Max iterations have been exceeded. |   
             %------------------------------------% */

	if (iter > *mxiter && nconv < numcnv) {
	    *info = 1;
	}

/*           %---------------------%   
             | No shifts to apply. |   
             %---------------------% */

	if (*np == 0 && nconv < numcnv) {
	    *info = 2;
	}

	*np = nconv;
	goto L1100;

    } else if (nconv < numcnv && *ishift == 1) {

/*           %-------------------------------------------------%   
             | Do not have all the requested eigenvalues yet.  |   
             | To prevent possible stagnation, adjust the size |   
             | of NEV.                                         |   
             %-------------------------------------------------% */

	nevbef = *nev;
/* Computing MIN */
	i__1 = nconv, i__2 = *np / 2;
	*nev += min(i__1,i__2);
	if (*nev == 1 && kplusp >= 6) {
	    *nev = kplusp / 2;
	} else if (*nev == 1 && kplusp > 3) {
	    *nev = 2;
	}
	*np = kplusp - *nev;

/*           %---------------------------------------%   
             | If the size of NEV was just increased |   
             | resort the eigenvalues.               |   
             %---------------------------------------% */

	if (nevbef < *nev) {
	    igraphdngets_(ishift, which, nev, np, &ritzr[1], &ritzi[1], &bounds[1], 
		    &workl[1], &workl[*np + 1]);
	}

    }

    if (msglvl > 0) {
	igraphivout_(&logfil, &c__1, &nconv, &ndigit, "_naup2: no. of \"converge"
		"d\" Ritz values at this iter.", (ftnlen)52);
	if (msglvl > 1) {
	    kp[0] = *nev;
	    kp[1] = *np;
	    igraphivout_(&logfil, &c__2, kp, &ndigit, "_naup2: NEV and NP are", (
		    ftnlen)22);
	    igraphdvout_(&logfil, nev, &ritzr[*np + 1], &ndigit, "_naup2: \"wante"
		    "d\" Ritz values -- real part", (ftnlen)41);
	    igraphdvout_(&logfil, nev, &ritzi[*np + 1], &ndigit, "_naup2: \"wante"
		    "d\" Ritz values -- imag part", (ftnlen)41);
	    igraphdvout_(&logfil, nev, &bounds[*np + 1], &ndigit, "_naup2: Ritz es"
		    "timates of the \"wanted\" values ", (ftnlen)46);
	}
    }

    if (*ishift == 0) {

/*           %-------------------------------------------------------%   
             | User specified shifts: reverse comminucation to       |   
             | compute the shifts. They are returned in the first    |   
             | 2*NP locations of WORKL.                              |   
             %-------------------------------------------------------% */

	ushift = TRUE_;
	*ido = 3;
	goto L9000;
    }

L50:

/*        %------------------------------------%   
          | Back from reverse communication;   |   
          | User specified shifts are returned |   
          | in WORKL(1:2*NP)                   |   
          %------------------------------------% */

    ushift = FALSE_;

    if (*ishift == 0) {

/*            %----------------------------------%   
              | Move the NP shifts from WORKL to |   
              | RITZR, RITZI to free up WORKL    |   
              | for non-exact shift case.        |   
              %----------------------------------% */

	igraphdcopy_(np, &workl[1], &c__1, &ritzr[1], &c__1);
	igraphdcopy_(np, &workl[*np + 1], &c__1, &ritzi[1], &c__1);
    }

    if (msglvl > 2) {
	igraphivout_(&logfil, &c__1, np, &ndigit, "_naup2: The number of shifts to"
		" apply ", (ftnlen)38);
	igraphdvout_(&logfil, np, &ritzr[1], &ndigit, "_naup2: Real part of the sh"
		"ifts", (ftnlen)31);
	igraphdvout_(&logfil, np, &ritzi[1], &ndigit, "_naup2: Imaginary part of t"
		"he shifts", (ftnlen)36);
	if (*ishift == 1) {
	    igraphdvout_(&logfil, np, &bounds[1], &ndigit, "_naup2: Ritz estimates"
		    " of the shifts", (ftnlen)36);
	}
    }

/*        %---------------------------------------------------------%   
          | Apply the NP implicit shifts by QR bulge chasing.       |   
          | Each shift is applied to the whole upper Hessenberg     |   
          | matrix H.                                               |   
          | The first 2*N locations of WORKD are used as workspace. |   
          %---------------------------------------------------------% */

    igraphdnapps_(n, nev, np, &ritzr[1], &ritzi[1], &v[v_offset], ldv, &h__[
	    h_offset], ldh, &resid[1], &q[q_offset], ldq, &workl[1], &workd[1]
	    );

/*        %---------------------------------------------%   
          | Compute the B-norm of the updated residual. |   
          | Keep B*RESID in WORKD(1:N) to be used in    |   
          | the first step of the next call to dnaitr.  |   
          %---------------------------------------------% */

    cnorm = TRUE_;
    igraphsecond_(&t2);
    if (*(unsigned char *)bmat == 'G') {
	++nbx;
	igraphdcopy_(n, &resid[1], &c__1, &workd[*n + 1], &c__1);
	ipntr[1] = *n + 1;
	ipntr[2] = 1;
	*ido = 2;

/*           %----------------------------------%   
             | Exit in order to compute B*RESID |   
             %----------------------------------% */

	goto L9000;
    } else if (*(unsigned char *)bmat == 'I') {
	igraphdcopy_(n, &resid[1], &c__1, &workd[1], &c__1);
    }

L100:

/*        %----------------------------------%   
          | Back from reverse communication; |   
          | WORKD(1:N) := B*RESID            |   
          %----------------------------------% */

    if (*(unsigned char *)bmat == 'G') {
	igraphsecond_(&t3);
	tmvbx += t3 - t2;
    }

    if (*(unsigned char *)bmat == 'G') {
	rnorm = igraphddot_(n, &resid[1], &c__1, &workd[1], &c__1);
	rnorm = sqrt((abs(rnorm)));
    } else if (*(unsigned char *)bmat == 'I') {
	rnorm = igraphdnrm2_(n, &resid[1], &c__1);
    }
    cnorm = FALSE_;

    if (msglvl > 2) {
	igraphdvout_(&logfil, &c__1, &rnorm, &ndigit, "_naup2: B-norm of residual "
		"for compressed factorization", (ftnlen)55);
	igraphdmout_(&logfil, nev, nev, &h__[h_offset], ldh, &ndigit, "_naup2: Com"
		"pressed upper Hessenberg matrix H", (ftnlen)44);
    }

    goto L1000;

/*     %---------------------------------------------------------------%   
       |                                                               |   
       |  E N D     O F     M A I N     I T E R A T I O N     L O O P  |   
       |                                                               |   
       %---------------------------------------------------------------% */

L1100:

    *mxiter = iter;
    *nev = numcnv;

L1200:
    *ido = 99;

/*     %------------%   
       | Error Exit |   
       %------------% */

    igraphsecond_(&t1);
    tnaup2 = t1 - t0;

L9000:

/*     %---------------%   
       | End of dnaup2 |   
       %---------------% */

    return 0;
} /* igraphdnaup2_ */

