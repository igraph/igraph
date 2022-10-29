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
static integer c__3 = 3;
static logical c_true = TRUE_;
static integer c__2 = 2;

/* -----------------------------------------------------------------------   
   \BeginDoc   

   \Name: dsaup2   

   \Description:   
    Intermediate level interface called by dsaupd.   

   \Usage:   
    call dsaup2   
       ( IDO, BMAT, N, WHICH, NEV, NP, TOL, RESID, MODE, IUPD,   
         ISHIFT, MXITER, V, LDV, H, LDH, RITZ, BOUNDS, Q, LDQ, WORKL,   
         IPNTR, WORKD, INFO )   

   \Arguments   

    IDO, BMAT, N, WHICH, NEV, TOL, RESID: same as defined in dsaupd.   
    MODE, ISHIFT, MXITER: see the definition of IPARAM in dsaupd.   

    NP      Integer.  (INPUT/OUTPUT)   
            Contains the number of implicit shifts to apply during   
            each Arnoldi/Lanczos iteration.   
            If ISHIFT=1, NP is adjusted dynamically at each iteration   
            to accelerate convergence and prevent stagnation.   
            This is also roughly equal to the number of matrix-vector   
            products (involving the operator OP) per Arnoldi iteration.   
            The logic for adjusting is contained within the current   
            subroutine.   
            If ISHIFT=0, NP is the number of shifts the user needs   
            to provide via reverse comunication. 0 < NP < NCV-NEV.   
            NP may be less than NCV-NEV since a leading block of the current   
            upper Tridiagonal matrix has split off and contains "unwanted"   
            Ritz values.   
            Upon termination of the IRA iteration, NP contains the number   
            of "converged" wanted Ritz values.   

    IUPD    Integer.  (INPUT)   
            IUPD .EQ. 0: use explicit restart instead implicit update.   
            IUPD .NE. 0: use implicit update.   

    V       Double precision N by (NEV+NP) array.  (INPUT/OUTPUT)   
            The Lanczos basis vectors.   

    LDV     Integer.  (INPUT)   
            Leading dimension of V exactly as declared in the calling   
            program.   

    H       Double precision (NEV+NP) by 2 array.  (OUTPUT)   
            H is used to store the generated symmetric tridiagonal matrix   
            The subdiagonal is stored in the first column of H starting   
            at H(2,1).  The main diagonal is stored in the second column   
            of H starting at H(1,2). If dsaup2 converges store the   
            B-norm of the final residual vector in H(1,1).   

    LDH     Integer.  (INPUT)   
            Leading dimension of H exactly as declared in the calling   
            program.   

    RITZ    Double precision array of length NEV+NP.  (OUTPUT)   
            RITZ(1:NEV) contains the computed Ritz values of OP.   

    BOUNDS  Double precision array of length NEV+NP.  (OUTPUT)   
            BOUNDS(1:NEV) contain the error bounds corresponding to RITZ.   

    Q       Double precision (NEV+NP) by (NEV+NP) array.  (WORKSPACE)   
            Private (replicated) work array used to accumulate the   
            rotation in the shift application step.   

    LDQ     Integer.  (INPUT)   
            Leading dimension of Q exactly as declared in the calling   
            program.   

    WORKL   Double precision array of length at least 3*(NEV+NP).  (INPUT/WORKSPACE)   
            Private (replicated) array on each PE or array allocated on   
            the front end.  It is used in the computation of the   
            tridiagonal eigenvalue problem, the calculation and   
            application of the shifts and convergence checking.   
            If ISHIFT .EQ. O and IDO .EQ. 3, the first NP locations   
            of WORKL are used in reverse communication to hold the user   
            supplied shifts.   

    IPNTR   Integer array of length 3.  (OUTPUT)   
            Pointer to mark the starting locations in the WORKD for   
            vectors used by the Lanczos iteration.   
            -------------------------------------------------------------   
            IPNTR(1): pointer to the current operand vector X.   
            IPNTR(2): pointer to the current result vector Y.   
            IPNTR(3): pointer to the vector B * X when used in one of   
                      the spectral transformation modes.  X is the current   
                      operand.   
            -------------------------------------------------------------   

    WORKD   Double precision work array of length 3*N.  (REVERSE COMMUNICATION)   
            Distributed array to be used in the basic Lanczos iteration   
            for reverse communication.  The user should not use WORKD   
            as temporary workspace during the iteration !!!!!!!!!!   
            See Data Distribution Note in dsaupd.   

    INFO    Integer.  (INPUT/OUTPUT)   
            If INFO .EQ. 0, a randomly initial residual vector is used.   
            If INFO .NE. 0, RESID contains the initial residual vector,   
                            possibly from a previous run.   
            Error flag on output.   
            =     0: Normal return.   
            =     1: All possible eigenvalues of OP has been found.   
                     NP returns the size of the invariant subspace   
                     spanning the operator OP.   
            =     2: No shifts could be applied.   
            =    -8: Error return from trid. eigenvalue calculation;   
                     This should never happen.   
            =    -9: Starting vector is zero.   
            = -9999: Could not build an Lanczos factorization.   
                     Size that was built in returned in NP.   

   \EndDoc   

   -----------------------------------------------------------------------   

   \BeginLib   

   \References:   
    1. D.C. Sorensen, "Implicit Application of Polynomial Filters in   
       a k-Step Arnoldi Method", SIAM J. Matr. Anal. Apps., 13 (1992),   
       pp 357-385.   
    2. R.B. Lehoucq, "Analysis and Implementation of an Implicitly   
       Restarted Arnoldi Iteration", Rice University Technical Report   
       TR95-13, Department of Computational and Applied Mathematics.   
    3. B.N. Parlett, "The Symmetric Eigenvalue Problem". Prentice-Hall,   
       1980.   
    4. B.N. Parlett, B. Nour-Omid, "Towards a Black Box Lanczos Program",   
       Computer Physics Communications, 53 (1989), pp 169-179.   
    5. B. Nour-Omid, B.N. Parlett, T. Ericson, P.S. Jensen, "How to   
       Implement the Spectral Transformation", Math. Comp., 48 (1987),   
       pp 663-673.   
    6. R.G. Grimes, J.G. Lewis and H.D. Simon, "A Shifted Block Lanczos   
       Algorithm for Solving Sparse Symmetric Generalized Eigenproblems",   
       SIAM J. Matr. Anal. Apps.,  January (1993).   
    7. L. Reichel, W.B. Gragg, "Algorithm 686: FORTRAN Subroutines   
       for Updating the QR decomposition", ACM TOMS, December 1990,   
       Volume 16 Number 4, pp 369-377.   

   \Routines called:   
       dgetv0  ARPACK initial vector generation routine.   
       dsaitr  ARPACK Lanczos factorization routine.   
       dsapps  ARPACK application of implicit shifts routine.   
       dsconv  ARPACK convergence of Ritz values routine.   
       dseigt  ARPACK compute Ritz values and error bounds routine.   
       dsgets  ARPACK reorder Ritz values and error bounds routine.   
       dsortr  ARPACK sorting routine.   
       ivout   ARPACK utility routine that prints integers.   
       second  ARPACK utility routine for timing.   
       dvout   ARPACK utility routine that prints vectors.   
       dlamch  LAPACK routine that determines machine constants.   
       dcopy   Level 1 BLAS that copies one vector to another.   
       ddot    Level 1 BLAS that computes the scalar product of two vectors.   
       dnrm2   Level 1 BLAS that computes the norm of a vector.   
       dscal   Level 1 BLAS that scales a vector.   
       dswap   Level 1 BLAS that swaps two vectors.   

   \Author   
       Danny Sorensen               Phuong Vu   
       Richard Lehoucq              CRPC / Rice University   
       Dept. of Computational &     Houston, Texas   
       Applied Mathematics   
       Rice University   
       Houston, Texas   

   \Revision history:   
       12/15/93: Version ' 2.4'   
       xx/xx/95: Version ' 2.4'.  (R.B. Lehoucq)   

   \SCCS Information: @(#)   
   FILE: saup2.F   SID: 2.6   DATE OF SID: 8/16/96   RELEASE: 2   

   \EndLib   

   -----------------------------------------------------------------------   

   Subroutine */ int igraphdsaup2_(integer *ido, char *bmat, integer *n, char *
	which, integer *nev, integer *np, doublereal *tol, doublereal *resid, 
	integer *mode, integer *iupd, integer *ishift, integer *mxiter, 
	doublereal *v, integer *ldv, doublereal *h__, integer *ldh, 
	doublereal *ritz, doublereal *bounds, doublereal *q, integer *ldq, 
	doublereal *workl, integer *ipntr, doublereal *workd, integer *info)
{
    /* System generated locals */
    integer h_dim1, h_offset, q_dim1, q_offset, v_dim1, v_offset, i__1, i__2, 
	    i__3;
    doublereal d__1, d__2, d__3;

    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *);
    integer s_cmp(char *, char *, ftnlen, ftnlen);
    /* Subroutine */ void s_copy(char *, char *, ftnlen, ftnlen);
    double sqrt(doublereal);

    /* Local variables */
    integer j;
    IGRAPH_F77_SAVE real t0, t1, t2, t3;
    integer kp[3];
    IGRAPH_F77_SAVE integer np0;
    integer nbx = 0;
    IGRAPH_F77_SAVE integer nev0;
    extern doublereal igraphddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    IGRAPH_F77_SAVE doublereal eps23;
    integer ierr;
    IGRAPH_F77_SAVE integer iter;
    doublereal temp;
    integer nevd2;
    extern doublereal igraphdnrm2_(integer *, doublereal *, integer *);
    IGRAPH_F77_SAVE logical getv0;
    integer nevm2;
    IGRAPH_F77_SAVE logical cnorm;
    extern /* Subroutine */ int igraphdcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), igraphdswap_(integer *, doublereal *, integer 
	    *, doublereal *, integer *);
    IGRAPH_F77_SAVE integer nconv;
    IGRAPH_F77_SAVE logical initv;
    IGRAPH_F77_SAVE doublereal rnorm;
    real tmvbx = 0.0;
    extern /* Subroutine */ int igraphdvout_(integer *, integer *, doublereal *, 
	    integer *, char *, ftnlen), igraphivout_(integer *, integer *, integer *
	    , integer *, char *, ftnlen), igraphdgetv0_(integer *, char *, integer *
	    , logical *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *);
    integer msaup2 = 0;
    real tsaup2;
    extern doublereal igraphdlamch_(char *);
    integer nevbef;
    extern /* Subroutine */ int igraphsecond_(real *);
    integer logfil = 0, ndigit;
    extern /* Subroutine */ int igraphdseigt_(doublereal *, integer *, doublereal *,
	     integer *, doublereal *, doublereal *, doublereal *, integer *);
    IGRAPH_F77_SAVE logical update;
    extern /* Subroutine */ int igraphdsaitr_(integer *, char *, integer *, integer 
	    *, integer *, integer *, doublereal *, doublereal *, doublereal *,
	     integer *, doublereal *, integer *, integer *, doublereal *, 
	    integer *), igraphdsgets_(integer *, char *, integer *, integer 
	    *, doublereal *, doublereal *, doublereal *), igraphdsapps_(
	    integer *, integer *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *), igraphdsconv_(integer *, doublereal *, 
	    doublereal *, doublereal *, integer *);
    IGRAPH_F77_SAVE logical ushift;
    char wprime[2];
    IGRAPH_F77_SAVE integer msglvl;
    integer nptemp;
    extern /* Subroutine */ int igraphdsortr_(char *, logical *, integer *, 
	    doublereal *, doublereal *);
    IGRAPH_F77_SAVE integer kplusp;


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
    --ritz;
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

/*        %-------------------------------%   
          | Initialize timing statistics  |   
          | & message level for debugging |   
          %-------------------------------% */

	igraphsecond_(&t0);
	msglvl = msaup2;

/*        %---------------------------------%   
          | Set machine dependent constant. |   
          %---------------------------------% */

	eps23 = igraphdlamch_("Epsilon-Machine");
	eps23 = pow_dd(&eps23, &c_b3);

/*        %-------------------------------------%   
          | nev0 and np0 are integer variables  |   
          | hold the initial values of NEV & NP |   
          %-------------------------------------% */

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

	kplusp = nev0 + np0;
	nconv = 0;
	iter = 0;

/*        %--------------------------------------------%   
          | Set flags for computing the first NEV steps |   
          | of the Lanczos factorization.              |   
          %--------------------------------------------% */

	getv0 = TRUE_;
	update = FALSE_;
	ushift = FALSE_;
	cnorm = FALSE_;

	if (*info != 0) {

/*        %--------------------------------------------%   
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
	    goto L1200;
	}
	getv0 = FALSE_;
	*ido = 0;
    }

/*     %------------------------------------------------------------%   
       | Back from reverse communication: continue with update step |   
       %------------------------------------------------------------% */

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
       | Compute the first NEV steps of the Lanczos factorization |   
       %----------------------------------------------------------% */

    igraphdsaitr_(ido, bmat, n, &c__0, &nev0, mode, &resid[1], &rnorm, &v[v_offset],
	     ldv, &h__[h_offset], ldh, &ipntr[1], &workd[1], info);

/*     %---------------------------------------------------%   
       | ido .ne. 99 implies use of reverse communication  |   
       | to compute operations involving OP and possibly B |   
       %---------------------------------------------------% */

    if (*ido != 99) {
	goto L9000;
    }

    if (*info > 0) {

/*        %-----------------------------------------------------%   
          | dsaitr was unable to build an Lanczos factorization |   
          | of length NEV0. INFO is returned with the size of   |   
          | the factorization built. Exit main loop.            |   
          %-----------------------------------------------------% */

	*np = *info;
	*mxiter = iter;
	*info = -9999;
	goto L1200;
    }

/*     %--------------------------------------------------------------%   
       |                                                              |   
       |           M A I N  LANCZOS  I T E R A T I O N  L O O P       |   
       |           Each iteration implicitly restarts the Lanczos     |   
       |           factorization in place.                            |   
       |                                                              |   
       %--------------------------------------------------------------% */

L1000:

    ++iter;

    if (msglvl > 0) {
	igraphivout_(&logfil, &c__1, &iter, &ndigit, "_saup2: **** Start of major "
		"iteration number ****", (ftnlen)49);
    }
    if (msglvl > 1) {
	igraphivout_(&logfil, &c__1, nev, &ndigit, "_saup2: The length of the curr"
		"ent Lanczos factorization", (ftnlen)55);
	igraphivout_(&logfil, &c__1, np, &ndigit, "_saup2: Extend the Lanczos fact"
		"orization by", (ftnlen)43);
    }

/*        %------------------------------------------------------------%   
          | Compute NP additional steps of the Lanczos factorization. |   
          %------------------------------------------------------------% */

    *ido = 0;
L20:
    update = TRUE_;

    igraphdsaitr_(ido, bmat, n, nev, np, mode, &resid[1], &rnorm, &v[v_offset], ldv,
	     &h__[h_offset], ldh, &ipntr[1], &workd[1], info);

/*        %---------------------------------------------------%   
          | ido .ne. 99 implies use of reverse communication  |   
          | to compute operations involving OP and possibly B |   
          %---------------------------------------------------% */

    if (*ido != 99) {
	goto L9000;
    }

    if (*info > 0) {

/*           %-----------------------------------------------------%   
             | dsaitr was unable to build an Lanczos factorization |   
             | of length NEV0+NP0. INFO is returned with the size  |   
             | of the factorization built. Exit main loop.         |   
             %-----------------------------------------------------% */

	*np = *info;
	*mxiter = iter;
	*info = -9999;
	goto L1200;
    }
    update = FALSE_;

    if (msglvl > 1) {
	igraphdvout_(&logfil, &c__1, &rnorm, &ndigit, "_saup2: Current B-norm of r"
		"esidual for factorization", (ftnlen)52);
    }

/*        %--------------------------------------------------------%   
          | Compute the eigenvalues and corresponding error bounds |   
          | of the current symmetric tridiagonal matrix.           |   
          %--------------------------------------------------------% */

    igraphdseigt_(&rnorm, &kplusp, &h__[h_offset], ldh, &ritz[1], &bounds[1], &
	    workl[1], &ierr);

    if (ierr != 0) {
	*info = -8;
	goto L1200;
    }

/*        %----------------------------------------------------%   
          | Make a copy of eigenvalues and corresponding error |   
          | bounds obtained from _seigt.                       |   
          %----------------------------------------------------% */

    igraphdcopy_(&kplusp, &ritz[1], &c__1, &workl[kplusp + 1], &c__1);
    igraphdcopy_(&kplusp, &bounds[1], &c__1, &workl[(kplusp << 1) + 1], &c__1);

/*        %---------------------------------------------------%   
          | Select the wanted Ritz values and their bounds    |   
          | to be used in the convergence test.               |   
          | The selection is based on the requested number of |   
          | eigenvalues instead of the current NEV and NP to  |   
          | prevent possible misconvergence.                  |   
          | * Wanted Ritz values := RITZ(NP+1:NEV+NP)         |   
          | * Shifts := RITZ(1:NP) := WORKL(1:NP)             |   
          %---------------------------------------------------% */

    *nev = nev0;
    *np = np0;
    igraphdsgets_(ishift, which, nev, np, &ritz[1], &bounds[1], &workl[1]);

/*        %-------------------%   
          | Convergence test. |   
          %-------------------% */

    igraphdcopy_(nev, &bounds[*np + 1], &c__1, &workl[*np + 1], &c__1);
    igraphdsconv_(nev, &ritz[*np + 1], &workl[*np + 1], tol, &nconv);

    if (msglvl > 2) {
	kp[0] = *nev;
	kp[1] = *np;
	kp[2] = nconv;
	igraphivout_(&logfil, &c__3, kp, &ndigit, "_saup2: NEV, NP, NCONV are", (
		ftnlen)26);
	igraphdvout_(&logfil, &kplusp, &ritz[1], &ndigit, "_saup2: The eigenvalues"
		" of H", (ftnlen)28);
	igraphdvout_(&logfil, &kplusp, &bounds[1], &ndigit, "_saup2: Ritz estimate"
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

    if (nconv >= nev0 || iter > *mxiter || *np == 0) {

/*           %------------------------------------------------%   
             | Prepare to exit. Put the converged Ritz values |   
             | and corresponding bounds in RITZ(1:NCONV) and  |   
             | BOUNDS(1:NCONV) respectively. Then sort. Be    |   
             | careful when NCONV > NP since we don't want to |   
             | swap overlapping locations.                    |   
             %------------------------------------------------% */

	if (s_cmp(which, "BE", (ftnlen)2, (ftnlen)2) == 0) {

/*              %-----------------------------------------------------%   
                | Both ends of the spectrum are requested.            |   
                | Sort the eigenvalues into algebraically decreasing  |   
                | order first then swap low end of the spectrum next  |   
                | to high end in appropriate locations.               |   
                | NOTE: when np < floor(nev/2) be careful not to swap |   
                | overlapping locations.                              |   
                %-----------------------------------------------------% */

	    s_copy(wprime, "SA", (ftnlen)2, (ftnlen)2);
	    igraphdsortr_(wprime, &c_true, &kplusp, &ritz[1], &bounds[1])
		    ;
	    nevd2 = *nev / 2;
	    nevm2 = *nev - nevd2;
	    if (*nev > 1) {
		i__1 = min(nevd2,*np);
/* Computing MAX */
		i__2 = kplusp - nevd2 + 1, i__3 = kplusp - *np + 1;
		igraphdswap_(&i__1, &ritz[nevm2 + 1], &c__1, &ritz[max(i__2,i__3)], 
			&c__1);
		i__1 = min(nevd2,*np);
/* Computing MAX */
		i__2 = kplusp - nevd2 + 1, i__3 = kplusp - *np;
		igraphdswap_(&i__1, &bounds[nevm2 + 1], &c__1, &bounds[max(i__2,
			i__3) + 1], &c__1);
	    }

	} else {

/*              %--------------------------------------------------%   
                | LM, SM, LA, SA case.                             |   
                | Sort the eigenvalues of H into the an order that |   
                | is opposite to WHICH, and apply the resulting    |   
                | order to BOUNDS.  The eigenvalues are sorted so  |   
                | that the wanted part are always within the first |   
                | NEV locations.                                   |   
                %--------------------------------------------------% */

	    if (s_cmp(which, "LM", (ftnlen)2, (ftnlen)2) == 0) {
		s_copy(wprime, "SM", (ftnlen)2, (ftnlen)2);
	    }
	    if (s_cmp(which, "SM", (ftnlen)2, (ftnlen)2) == 0) {
		s_copy(wprime, "LM", (ftnlen)2, (ftnlen)2);
	    }
	    if (s_cmp(which, "LA", (ftnlen)2, (ftnlen)2) == 0) {
		s_copy(wprime, "SA", (ftnlen)2, (ftnlen)2);
	    }
	    if (s_cmp(which, "SA", (ftnlen)2, (ftnlen)2) == 0) {
		s_copy(wprime, "LA", (ftnlen)2, (ftnlen)2);
	    }

	    igraphdsortr_(wprime, &c_true, &kplusp, &ritz[1], &bounds[1])
		    ;

	}

/*           %--------------------------------------------------%   
             | Scale the Ritz estimate of each Ritz value       |   
             | by 1 / max(eps23,magnitude of the Ritz value).   |   
             %--------------------------------------------------% */

	i__1 = nev0;
	for (j = 1; j <= i__1; ++j) {
/* Computing MAX */
	    d__2 = eps23, d__3 = (d__1 = ritz[j], abs(d__1));
	    temp = max(d__2,d__3);
	    bounds[j] /= temp;
/* L35: */
	}

/*           %----------------------------------------------------%   
             | Sort the Ritz values according to the scaled Ritz  |   
             | esitmates.  This will push all the converged ones  |   
             | towards the front of ritzr, ritzi, bounds          |   
             | (in the case when NCONV < NEV.)                    |   
             %----------------------------------------------------% */

	s_copy(wprime, "LA", (ftnlen)2, (ftnlen)2);
	igraphdsortr_(wprime, &c_true, &nev0, &bounds[1], &ritz[1]);

/*           %----------------------------------------------%   
             | Scale the Ritz estimate back to its original |   
             | value.                                       |   
             %----------------------------------------------% */

	i__1 = nev0;
	for (j = 1; j <= i__1; ++j) {
/* Computing MAX */
	    d__2 = eps23, d__3 = (d__1 = ritz[j], abs(d__1));
	    temp = max(d__2,d__3);
	    bounds[j] *= temp;
/* L40: */
	}

/*           %--------------------------------------------------%   
             | Sort the "converged" Ritz values again so that   |   
             | the "threshold" values and their associated Ritz |   
             | estimates appear at the appropriate position in  |   
             | ritz and bound.                                  |   
             %--------------------------------------------------% */

	if (s_cmp(which, "BE", (ftnlen)2, (ftnlen)2) == 0) {

/*              %------------------------------------------------%   
                | Sort the "converged" Ritz values in increasing |   
                | order.  The "threshold" values are in the      |   
                | middle.                                        |   
                %------------------------------------------------% */

	    s_copy(wprime, "LA", (ftnlen)2, (ftnlen)2);
	    igraphdsortr_(wprime, &c_true, &nconv, &ritz[1], &bounds[1]);

	} else {

/*              %----------------------------------------------%   
                | In LM, SM, LA, SA case, sort the "converged" |   
                | Ritz values according to WHICH so that the   |   
                | "threshold" value appears at the front of    |   
                | ritz.                                        |   
                %----------------------------------------------% */
	    igraphdsortr_(which, &c_true, &nconv, &ritz[1], &bounds[1]);

	}

/*           %------------------------------------------%   
             |  Use h( 1,1 ) as storage to communicate  |   
             |  rnorm to _seupd if needed               |   
             %------------------------------------------% */

	h__[h_dim1 + 1] = rnorm;

	if (msglvl > 1) {
	    igraphdvout_(&logfil, &kplusp, &ritz[1], &ndigit, "_saup2: Sorted Ritz"
		    " values.", (ftnlen)27);
	    igraphdvout_(&logfil, &kplusp, &bounds[1], &ndigit, "_saup2: Sorted ri"
		    "tz estimates.", (ftnlen)30);
	}

/*           %------------------------------------%   
             | Max iterations have been exceeded. |   
             %------------------------------------% */

	if (iter > *mxiter && nconv < *nev) {
	    *info = 1;
	}

/*           %---------------------%   
             | No shifts to apply. |   
             %---------------------% */

	if (*np == 0 && nconv < nev0) {
	    *info = 2;
	}

	*np = nconv;
	goto L1100;

    } else if (nconv < *nev && *ishift == 1) {

/*           %---------------------------------------------------%   
             | Do not have all the requested eigenvalues yet.    |   
             | To prevent possible stagnation, adjust the number |   
             | of Ritz values and the shifts.                    |   
             %---------------------------------------------------% */

	nevbef = *nev;
/* Computing MIN */
	i__1 = nconv, i__2 = *np / 2;
	*nev += min(i__1,i__2);
	if (*nev == 1 && kplusp >= 6) {
	    *nev = kplusp / 2;
	} else if (*nev == 1 && kplusp > 2) {
	    *nev = 2;
	}
	*np = kplusp - *nev;

/*           %---------------------------------------%   
             | If the size of NEV was just increased |   
             | resort the eigenvalues.               |   
             %---------------------------------------% */

	if (nevbef < *nev) {
	    igraphdsgets_(ishift, which, nev, np, &ritz[1], &bounds[1], &workl[1]);
	}

    }

    if (msglvl > 0) {
	igraphivout_(&logfil, &c__1, &nconv, &ndigit, "_saup2: no. of \"converge"
		"d\" Ritz values at this iter.", (ftnlen)52);
	if (msglvl > 1) {
	    kp[0] = *nev;
	    kp[1] = *np;
	    igraphivout_(&logfil, &c__2, kp, &ndigit, "_saup2: NEV and NP are", (
		    ftnlen)22);
	    igraphdvout_(&logfil, nev, &ritz[*np + 1], &ndigit, "_saup2: \"wante"
		    "d\" Ritz values.", (ftnlen)29);
	    igraphdvout_(&logfil, nev, &bounds[*np + 1], &ndigit, "_saup2: Ritz es"
		    "timates of the \"wanted\" values ", (ftnlen)46);
	}
    }

    if (*ishift == 0) {

/*           %-----------------------------------------------------%   
             | User specified shifts: reverse communication to     |   
             | compute the shifts. They are returned in the first  |   
             | NP locations of WORKL.                              |   
             %-----------------------------------------------------% */

	ushift = TRUE_;
	*ido = 3;
	goto L9000;
    }

L50:

/*        %------------------------------------%   
          | Back from reverse communication;   |   
          | User specified shifts are returned |   
          | in WORKL(1:*NP)                   |   
          %------------------------------------% */

    ushift = FALSE_;


/*        %---------------------------------------------------------%   
          | Move the NP shifts to the first NP locations of RITZ to |   
          | free up WORKL.  This is for the non-exact shift case;   |   
          | in the exact shift case, dsgets already handles this.   |   
          %---------------------------------------------------------% */

    if (*ishift == 0) {
	igraphdcopy_(np, &workl[1], &c__1, &ritz[1], &c__1);
    }

    if (msglvl > 2) {
	igraphivout_(&logfil, &c__1, np, &ndigit, "_saup2: The number of shifts to"
		" apply ", (ftnlen)38);
	igraphdvout_(&logfil, np, &workl[1], &ndigit, "_saup2: shifts selected", (
		ftnlen)23);
	if (*ishift == 1) {
	    igraphdvout_(&logfil, np, &bounds[1], &ndigit, "_saup2: corresponding "
		    "Ritz estimates", (ftnlen)36);
	}
    }

/*        %---------------------------------------------------------%   
          | Apply the NP0 implicit shifts by QR bulge chasing.      |   
          | Each shift is applied to the entire tridiagonal matrix. |   
          | The first 2*N locations of WORKD are used as workspace. |   
          | After dsapps is done, we have a Lanczos                 |   
          | factorization of length NEV.                            |   
          %---------------------------------------------------------% */

    igraphdsapps_(n, nev, np, &ritz[1], &v[v_offset], ldv, &h__[h_offset], ldh, &
	    resid[1], &q[q_offset], ldq, &workd[1]);

/*        %---------------------------------------------%   
          | Compute the B-norm of the updated residual. |   
          | Keep B*RESID in WORKD(1:N) to be used in    |   
          | the first step of the next call to dsaitr.  |   
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
/* L130: */

    if (msglvl > 2) {
	igraphdvout_(&logfil, &c__1, &rnorm, &ndigit, "_saup2: B-norm of residual "
		"for NEV factorization", (ftnlen)48);
	igraphdvout_(&logfil, nev, &h__[(h_dim1 << 1) + 1], &ndigit, "_saup2: main"
		" diagonal of compressed H matrix", (ftnlen)44);
	i__1 = *nev - 1;
	igraphdvout_(&logfil, &i__1, &h__[h_dim1 + 2], &ndigit, "_saup2: subdiagon"
		"al of compressed H matrix", (ftnlen)42);
    }

    goto L1000;

/*     %---------------------------------------------------------------%   
       |                                                               |   
       |  E N D     O F     M A I N     I T E R A T I O N     L O O P  |   
       |                                                               |   
       %---------------------------------------------------------------% */

L1100:

    *mxiter = iter;
    *nev = nconv;

L1200:
    *ido = 99;

/*     %------------%   
       | Error exit |   
       %------------% */

    igraphsecond_(&t1);
    tsaup2 = t1 - t0;

L9000:
    return 0;

/*     %---------------%   
       | End of dsaup2 |   
       %---------------% */

} /* igraphdsaup2_ */

