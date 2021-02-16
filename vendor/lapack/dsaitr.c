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

static integer c__1 = 1;
static logical c_false = FALSE_;
static doublereal c_b24 = 1.;
static doublereal c_b49 = 0.;
static doublereal c_b57 = -1.;
static integer c__2 = 2;

/* -----------------------------------------------------------------------   
   \BeginDoc   

   \Name: dsaitr   

   \Description:   
    Reverse communication interface for applying NP additional steps to   
    a K step symmetric Arnoldi factorization.   

    Input:  OP*V_{k}  -  V_{k}*H = r_{k}*e_{k}^T   

            with (V_{k}^T)*B*V_{k} = I, (V_{k}^T)*B*r_{k} = 0.   

    Output: OP*V_{k+p}  -  V_{k+p}*H = r_{k+p}*e_{k+p}^T   

            with (V_{k+p}^T)*B*V_{k+p} = I, (V_{k+p}^T)*B*r_{k+p} = 0.   

    where OP and B are as in dsaupd.  The B-norm of r_{k+p} is also   
    computed and returned.   

   \Usage:   
    call dsaitr   
       ( IDO, BMAT, N, K, NP, MODE, RESID, RNORM, V, LDV, H, LDH,   
         IPNTR, WORKD, INFO )   

   \Arguments   
    IDO     Integer.  (INPUT/OUTPUT)   
            Reverse communication flag.   
            -------------------------------------------------------------   
            IDO =  0: first call to the reverse communication interface   
            IDO = -1: compute  Y = OP * X  where   
                      IPNTR(1) is the pointer into WORK for X,   
                      IPNTR(2) is the pointer into WORK for Y.   
                      This is for the restart phase to force the new   
                      starting vector into the range of OP.   
            IDO =  1: compute  Y = OP * X  where   
                      IPNTR(1) is the pointer into WORK for X,   
                      IPNTR(2) is the pointer into WORK for Y,   
                      IPNTR(3) is the pointer into WORK for B * X.   
            IDO =  2: compute  Y = B * X  where   
                      IPNTR(1) is the pointer into WORK for X,   
                      IPNTR(2) is the pointer into WORK for Y.   
            IDO = 99: done   
            -------------------------------------------------------------   
            When the routine is used in the "shift-and-invert" mode, the   
            vector B * Q is already available and does not need to be   
            recomputed in forming OP * Q.   

    BMAT    Character*1.  (INPUT)   
            BMAT specifies the type of matrix B that defines the   
            semi-inner product for the operator OP.  See dsaupd.   
            B = 'I' -> standard eigenvalue problem A*x = lambda*x   
            B = 'G' -> generalized eigenvalue problem A*x = lambda*M*x   

    N       Integer.  (INPUT)   
            Dimension of the eigenproblem.   

    K       Integer.  (INPUT)   
            Current order of H and the number of columns of V.   

    NP      Integer.  (INPUT)   
            Number of additional Arnoldi steps to take.   

    MODE    Integer.  (INPUT)   
            Signifies which form for "OP". If MODE=2 then   
            a reduction in the number of B matrix vector multiplies   
            is possible since the B-norm of OP*x is equivalent to   
            the inv(B)-norm of A*x.   

    RESID   Double precision array of length N.  (INPUT/OUTPUT)   
            On INPUT:  RESID contains the residual vector r_{k}.   
            On OUTPUT: RESID contains the residual vector r_{k+p}.   

    RNORM   Double precision scalar.  (INPUT/OUTPUT)   
            On INPUT the B-norm of r_{k}.   
            On OUTPUT the B-norm of the updated residual r_{k+p}.   

    V       Double precision N by K+NP array.  (INPUT/OUTPUT)   
            On INPUT:  V contains the Arnoldi vectors in the first K   
            columns.   
            On OUTPUT: V contains the new NP Arnoldi vectors in the next   
            NP columns.  The first K columns are unchanged.   

    LDV     Integer.  (INPUT)   
            Leading dimension of V exactly as declared in the calling   
            program.   

    H       Double precision (K+NP) by 2 array.  (INPUT/OUTPUT)   
            H is used to store the generated symmetric tridiagonal matrix   
            with the subdiagonal in the first column starting at H(2,1)   
            and the main diagonal in the second column.   

    LDH     Integer.  (INPUT)   
            Leading dimension of H exactly as declared in the calling   
            program.   

    IPNTR   Integer array of length 3.  (OUTPUT)   
            Pointer to mark the starting locations in the WORK for   
            vectors used by the Arnoldi iteration.   
            -------------------------------------------------------------   
            IPNTR(1): pointer to the current operand vector X.   
            IPNTR(2): pointer to the current result vector Y.   
            IPNTR(3): pointer to the vector B * X when used in the   
                      shift-and-invert mode.  X is the current operand.   
            -------------------------------------------------------------   

    WORKD   Double precision work array of length 3*N.  (REVERSE COMMUNICATION)   
            Distributed array to be used in the basic Arnoldi iteration   
            for reverse communication.  The calling program should not   
            use WORKD as temporary workspace during the iteration !!!!!!   
            On INPUT, WORKD(1:N) = B*RESID where RESID is associated   
            with the K step Arnoldi factorization. Used to save some   
            computation at the first step.   
            On OUTPUT, WORKD(1:N) = B*RESID where RESID is associated   
            with the K+NP step Arnoldi factorization.   

    INFO    Integer.  (OUTPUT)   
            = 0: Normal exit.   
            > 0: Size of an invariant subspace of OP is found that is   
                 less than K + NP.   

   \EndDoc   

   -----------------------------------------------------------------------   

   \BeginLib   

   \Local variables:   
       xxxxxx  real   

   \Routines called:   
       dgetv0  ARPACK routine to generate the initial vector.   
       ivout   ARPACK utility routine that prints integers.   
       dmout   ARPACK utility routine that prints matrices.   
       dvout   ARPACK utility routine that prints vectors.   
       dlamch  LAPACK routine that determines machine constants.   
       dlascl  LAPACK routine for careful scaling of a matrix.   
       dgemv   Level 2 BLAS routine for matrix vector multiplication.   
       daxpy   Level 1 BLAS that computes a vector triad.   
       dscal   Level 1 BLAS that scales a vector.   
       dcopy   Level 1 BLAS that copies one vector to another .   
       ddot    Level 1 BLAS that computes the scalar product of two vectors.   
       dnrm2   Level 1 BLAS that computes the norm of a vector.   

   \Author   
       Danny Sorensen               Phuong Vu   
       Richard Lehoucq              CRPC / Rice University   
       Dept. of Computational &     Houston, Texas   
       Applied Mathematics   
       Rice University   
       Houston, Texas   

   \Revision history:   
       xx/xx/93: Version ' 2.4'   

   \SCCS Information: @(#)   
   FILE: saitr.F   SID: 2.6   DATE OF SID: 8/28/96   RELEASE: 2   

   \Remarks   
    The algorithm implemented is:   

    restart = .false.   
    Given V_{k} = [v_{1}, ..., v_{k}], r_{k};   
    r_{k} contains the initial residual vector even for k = 0;   
    Also assume that rnorm = || B*r_{k} || and B*r_{k} are already   
    computed by the calling program.   

    betaj = rnorm ; p_{k+1} = B*r_{k} ;   
    For  j = k+1, ..., k+np  Do   
       1) if ( betaj < tol ) stop or restart depending on j.   
          if ( restart ) generate a new starting vector.   
       2) v_{j} = r(j-1)/betaj;  V_{j} = [V_{j-1}, v_{j}];   
          p_{j} = p_{j}/betaj   
       3) r_{j} = OP*v_{j} where OP is defined as in dsaupd   
          For shift-invert mode p_{j} = B*v_{j} is already available.   
          wnorm = || OP*v_{j} ||   
       4) Compute the j-th step residual vector.   
          w_{j} =  V_{j}^T * B * OP * v_{j}   
          r_{j} =  OP*v_{j} - V_{j} * w_{j}   
          alphaj <- j-th component of w_{j}   
          rnorm = || r_{j} ||   
          betaj+1 = rnorm   
          If (rnorm > 0.717*wnorm) accept step and go back to 1)   
       5) Re-orthogonalization step:   
          s = V_{j}'*B*r_{j}   
          r_{j} = r_{j} - V_{j}*s;  rnorm1 = || r_{j} ||   
          alphaj = alphaj + s_{j};   
       6) Iterative refinement step:   
          If (rnorm1 > 0.717*rnorm) then   
             rnorm = rnorm1   
             accept step and go back to 1)   
          Else   
             rnorm = rnorm1   
             If this is the first time in step 6), go to 5)   
             Else r_{j} lies in the span of V_{j} numerically.   
                Set r_{j} = 0 and rnorm = 0; go to 1)   
          EndIf   
    End Do   

   \EndLib   

   -----------------------------------------------------------------------   

   Subroutine */ int igraphdsaitr_(integer *ido, char *bmat, integer *n, integer *k,
	 integer *np, integer *mode, doublereal *resid, doublereal *rnorm, 
	doublereal *v, integer *ldv, doublereal *h__, integer *ldh, integer *
	ipntr, doublereal *workd, integer *info)
{
    /* Initialized data */

    IGRAPH_F77_SAVE logical first = TRUE_;

    /* System generated locals */
    integer h_dim1, h_offset, v_dim1, v_offset, i__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    integer i__;
    IGRAPH_F77_SAVE integer j;
    IGRAPH_F77_SAVE real t0, t1, t2, t3, t4, t5;
    integer jj;
    IGRAPH_F77_SAVE integer ipj, irj;
    integer nbx = 0;
    IGRAPH_F77_SAVE integer ivj;
    extern doublereal igraphddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    IGRAPH_F77_SAVE integer ierr, iter;
    integer nopx = 0;
    IGRAPH_F77_SAVE integer itry;
    extern doublereal igraphdnrm2_(integer *, doublereal *, integer *);
    doublereal temp1;
    IGRAPH_F77_SAVE logical orth1, orth2, step3, step4;
    extern /* Subroutine */ int igraphdscal_(integer *, doublereal *, doublereal *, 
	    integer *), igraphdgemv_(char *, integer *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *);
    integer infol;
    extern /* Subroutine */ int igraphdcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    doublereal xtemp[2];
    real tmvbx = 0;
    extern /* Subroutine */ int igraphdvout_(integer *, integer *, doublereal *, 
	    integer *, char *, ftnlen);
    IGRAPH_F77_SAVE doublereal wnorm;
    extern /* Subroutine */ int igraphivout_(integer *, integer *, integer *, 
	    integer *, char *, ftnlen), igraphdgetv0_(integer *, char *, integer *, 
	    logical *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *);
    IGRAPH_F77_SAVE doublereal rnorm1;
    extern doublereal igraphdlamch_(char *);
    extern /* Subroutine */ int igraphdlascl_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, integer *), igraphsecond_(real *);
    integer logfil;
    IGRAPH_F77_SAVE doublereal safmin;
    integer ndigit = 0, nitref = 0;
    real titref = 0;
    integer msaitr = 0;
    IGRAPH_F77_SAVE integer msglvl;
    real tsaitr = 0;
    integer nrorth = 0;
    IGRAPH_F77_SAVE logical rstart;
    integer nrstrt = 0;
    real tmvopx = 0;


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
       | Local Array Arguments |   
       %-----------------------%   


       %----------------------%   
       | External Subroutines |   
       %----------------------%   


       %--------------------%   
       | External Functions |   
       %--------------------%   


       %-----------------%   
       | Data statements |   
       %-----------------%   

       Parameter adjustments */
    --workd;
    --resid;
    v_dim1 = *ldv;
    v_offset = 1 + v_dim1;
    v -= v_offset;
    h_dim1 = *ldh;
    h_offset = 1 + h_dim1;
    h__ -= h_offset;
    --ipntr;

    /* Function Body   

       %-----------------------%   
       | Executable Statements |   
       %-----------------------% */

    if (first) {
	first = FALSE_;

/*        %--------------------------------%   
          | safmin = safe minimum is such  |   
          | that 1/sfmin does not overflow |   
          %--------------------------------% */

	safmin = igraphdlamch_("safmin");
    }

    if (*ido == 0) {

/*        %-------------------------------%   
          | Initialize timing statistics  |   
          | & message level for debugging |   
          %-------------------------------% */

	igraphsecond_(&t0);
	msglvl = msaitr;

/*        %------------------------------%   
          | Initial call to this routine |   
          %------------------------------% */

	*info = 0;
	step3 = FALSE_;
	step4 = FALSE_;
	rstart = FALSE_;
	orth1 = FALSE_;
	orth2 = FALSE_;

/*        %--------------------------------%   
          | Pointer to the current step of |   
          | the factorization to build     |   
          %--------------------------------% */

	j = *k + 1;

/*        %------------------------------------------%   
          | Pointers used for reverse communication  |   
          | when using WORKD.                        |   
          %------------------------------------------% */

	ipj = 1;
	irj = ipj + *n;
	ivj = irj + *n;
    }

/*     %-------------------------------------------------%   
       | When in reverse communication mode one of:      |   
       | STEP3, STEP4, ORTH1, ORTH2, RSTART              |   
       | will be .true.                                  |   
       | STEP3: return from computing OP*v_{j}.          |   
       | STEP4: return from computing B-norm of OP*v_{j} |   
       | ORTH1: return from computing B-norm of r_{j+1}  |   
       | ORTH2: return from computing B-norm of          |   
       |        correction to the residual vector.       |   
       | RSTART: return from OP computations needed by   |   
       |         dgetv0.                                 |   
       %-------------------------------------------------% */

    if (step3) {
	goto L50;
    }
    if (step4) {
	goto L60;
    }
    if (orth1) {
	goto L70;
    }
    if (orth2) {
	goto L90;
    }
    if (rstart) {
	goto L30;
    }

/*     %------------------------------%   
       | Else this is the first step. |   
       %------------------------------%   

       %--------------------------------------------------------------%   
       |                                                              |   
       |        A R N O L D I     I T E R A T I O N     L O O P       |   
       |                                                              |   
       | Note:  B*r_{j-1} is already in WORKD(1:N)=WORKD(IPJ:IPJ+N-1) |   
       %--------------------------------------------------------------% */

L1000:

    if (msglvl > 2) {
	igraphivout_(&logfil, &c__1, &j, &ndigit, "_saitr: generating Arnoldi vect"
		"or no.", (ftnlen)37);
	igraphdvout_(&logfil, &c__1, rnorm, &ndigit, "_saitr: B-norm of the curren"
		"t residual =", (ftnlen)40);
    }

/*        %---------------------------------------------------------%   
          | Check for exact zero. Equivalent to determing whether a |   
          | j-step Arnoldi factorization is present.                |   
          %---------------------------------------------------------% */

    if (*rnorm > 0.) {
	goto L40;
    }

/*           %---------------------------------------------------%   
             | Invariant subspace found, generate a new starting |   
             | vector which is orthogonal to the current Arnoldi |   
             | basis and continue the iteration.                 |   
             %---------------------------------------------------% */

    if (msglvl > 0) {
	igraphivout_(&logfil, &c__1, &j, &ndigit, "_saitr: ****** restart at step "
		"******", (ftnlen)37);
    }

/*           %---------------------------------------------%   
             | ITRY is the loop variable that controls the |   
             | maximum amount of times that a restart is   |   
             | attempted. NRSTRT is used by stat.h         |   
             %---------------------------------------------% */

    ++nrstrt;
    itry = 1;
L20:
    rstart = TRUE_;
    *ido = 0;
L30:

/*           %--------------------------------------%   
             | If in reverse communication mode and |   
             | RSTART = .true. flow returns here.   |   
             %--------------------------------------% */

    igraphdgetv0_(ido, bmat, &itry, &c_false, n, &j, &v[v_offset], ldv, &resid[1], 
	    rnorm, &ipntr[1], &workd[1], &ierr);
    if (*ido != 99) {
	goto L9000;
    }
    if (ierr < 0) {
	++itry;
	if (itry <= 3) {
	    goto L20;
	}

/*              %------------------------------------------------%   
                | Give up after several restart attempts.        |   
                | Set INFO to the size of the invariant subspace |   
                | which spans OP and exit.                       |   
                %------------------------------------------------% */

	*info = j - 1;
	igraphsecond_(&t1);
	tsaitr += t1 - t0;
	*ido = 99;
	goto L9000;
    }

L40:

/*        %---------------------------------------------------------%   
          | STEP 2:  v_{j} = r_{j-1}/rnorm and p_{j} = p_{j}/rnorm  |   
          | Note that p_{j} = B*r_{j-1}. In order to avoid overflow |   
          | when reciprocating a small RNORM, test against lower    |   
          | machine bound.                                          |   
          %---------------------------------------------------------% */

    igraphdcopy_(n, &resid[1], &c__1, &v[j * v_dim1 + 1], &c__1);
    if (*rnorm >= safmin) {
	temp1 = 1. / *rnorm;
	igraphdscal_(n, &temp1, &v[j * v_dim1 + 1], &c__1);
	igraphdscal_(n, &temp1, &workd[ipj], &c__1);
    } else {

/*            %-----------------------------------------%   
              | To scale both v_{j} and p_{j} carefully |   
              | use LAPACK routine SLASCL               |   
              %-----------------------------------------% */

	igraphdlascl_("General", &i__, &i__, rnorm, &c_b24, n, &c__1, &v[j * v_dim1 
		+ 1], n, &infol);
	igraphdlascl_("General", &i__, &i__, rnorm, &c_b24, n, &c__1, &workd[ipj], 
		n, &infol);
    }

/*        %------------------------------------------------------%   
          | STEP 3:  r_{j} = OP*v_{j}; Note that p_{j} = B*v_{j} |   
          | Note that this is not quite yet r_{j}. See STEP 4    |   
          %------------------------------------------------------% */

    step3 = TRUE_;
    ++nopx;
    igraphsecond_(&t2);
    igraphdcopy_(n, &v[j * v_dim1 + 1], &c__1, &workd[ivj], &c__1);
    ipntr[1] = ivj;
    ipntr[2] = irj;
    ipntr[3] = ipj;
    *ido = 1;

/*        %-----------------------------------%   
          | Exit in order to compute OP*v_{j} |   
          %-----------------------------------% */

    goto L9000;
L50:

/*        %-----------------------------------%   
          | Back from reverse communication;  |   
          | WORKD(IRJ:IRJ+N-1) := OP*v_{j}.   |   
          %-----------------------------------% */

    igraphsecond_(&t3);
    tmvopx += t3 - t2;

    step3 = FALSE_;

/*        %------------------------------------------%   
          | Put another copy of OP*v_{j} into RESID. |   
          %------------------------------------------% */

    igraphdcopy_(n, &workd[irj], &c__1, &resid[1], &c__1);

/*        %-------------------------------------------%   
          | STEP 4:  Finish extending the symmetric   |   
          |          Arnoldi to length j. If MODE = 2 |   
          |          then B*OP = B*inv(B)*A = A and   |   
          |          we don't need to compute B*OP.   |   
          | NOTE: If MODE = 2 WORKD(IVJ:IVJ+N-1) is   |   
          | assumed to have A*v_{j}.                  |   
          %-------------------------------------------% */

    if (*mode == 2) {
	goto L65;
    }
    igraphsecond_(&t2);
    if (*(unsigned char *)bmat == 'G') {
	++nbx;
	step4 = TRUE_;
	ipntr[1] = irj;
	ipntr[2] = ipj;
	*ido = 2;

/*           %-------------------------------------%   
             | Exit in order to compute B*OP*v_{j} |   
             %-------------------------------------% */

	goto L9000;
    } else if (*(unsigned char *)bmat == 'I') {
	igraphdcopy_(n, &resid[1], &c__1, &workd[ipj], &c__1);
    }
L60:

/*        %-----------------------------------%   
          | Back from reverse communication;  |   
          | WORKD(IPJ:IPJ+N-1) := B*OP*v_{j}. |   
          %-----------------------------------% */

    if (*(unsigned char *)bmat == 'G') {
	igraphsecond_(&t3);
	tmvbx += t3 - t2;
    }

    step4 = FALSE_;

/*        %-------------------------------------%   
          | The following is needed for STEP 5. |   
          | Compute the B-norm of OP*v_{j}.     |   
          %-------------------------------------% */

L65:
    if (*mode == 2) {

/*           %----------------------------------%   
             | Note that the B-norm of OP*v_{j} |   
             | is the inv(B)-norm of A*v_{j}.   |   
             %----------------------------------% */

	wnorm = igraphddot_(n, &resid[1], &c__1, &workd[ivj], &c__1);
	wnorm = sqrt((abs(wnorm)));
    } else if (*(unsigned char *)bmat == 'G') {
	wnorm = igraphddot_(n, &resid[1], &c__1, &workd[ipj], &c__1);
	wnorm = sqrt((abs(wnorm)));
    } else if (*(unsigned char *)bmat == 'I') {
	wnorm = igraphdnrm2_(n, &resid[1], &c__1);
    }

/*        %-----------------------------------------%   
          | Compute the j-th residual corresponding |   
          | to the j step factorization.            |   
          | Use Classical Gram Schmidt and compute: |   
          | w_{j} <-  V_{j}^T * B * OP * v_{j}      |   
          | r_{j} <-  OP*v_{j} - V_{j} * w_{j}      |   
          %-----------------------------------------%   


          %------------------------------------------%   
          | Compute the j Fourier coefficients w_{j} |   
          | WORKD(IPJ:IPJ+N-1) contains B*OP*v_{j}.  |   
          %------------------------------------------% */

    if (*mode != 2) {
	igraphdgemv_("T", n, &j, &c_b24, &v[v_offset], ldv, &workd[ipj], &c__1, &
		c_b49, &workd[irj], &c__1);
    } else if (*mode == 2) {
	igraphdgemv_("T", n, &j, &c_b24, &v[v_offset], ldv, &workd[ivj], &c__1, &
		c_b49, &workd[irj], &c__1);
    }

/*        %--------------------------------------%   
          | Orthgonalize r_{j} against V_{j}.    |   
          | RESID contains OP*v_{j}. See STEP 3. |   
          %--------------------------------------% */

    igraphdgemv_("N", n, &j, &c_b57, &v[v_offset], ldv, &workd[irj], &c__1, &c_b24, 
	    &resid[1], &c__1);

/*        %--------------------------------------%   
          | Extend H to have j rows and columns. |   
          %--------------------------------------% */

    h__[j + (h_dim1 << 1)] = workd[irj + j - 1];
    if (j == 1 || rstart) {
	h__[j + h_dim1] = 0.;
    } else {
	h__[j + h_dim1] = *rnorm;
    }
    igraphsecond_(&t4);

    orth1 = TRUE_;
    iter = 0;

    igraphsecond_(&t2);
    if (*(unsigned char *)bmat == 'G') {
	++nbx;
	igraphdcopy_(n, &resid[1], &c__1, &workd[irj], &c__1);
	ipntr[1] = irj;
	ipntr[2] = ipj;
	*ido = 2;

/*           %----------------------------------%   
             | Exit in order to compute B*r_{j} |   
             %----------------------------------% */

	goto L9000;
    } else if (*(unsigned char *)bmat == 'I') {
	igraphdcopy_(n, &resid[1], &c__1, &workd[ipj], &c__1);
    }
L70:

/*        %---------------------------------------------------%   
          | Back from reverse communication if ORTH1 = .true. |   
          | WORKD(IPJ:IPJ+N-1) := B*r_{j}.                    |   
          %---------------------------------------------------% */

    if (*(unsigned char *)bmat == 'G') {
	igraphsecond_(&t3);
	tmvbx += t3 - t2;
    }

    orth1 = FALSE_;

/*        %------------------------------%   
          | Compute the B-norm of r_{j}. |   
          %------------------------------% */

    if (*(unsigned char *)bmat == 'G') {
	*rnorm = igraphddot_(n, &resid[1], &c__1, &workd[ipj], &c__1);
	*rnorm = sqrt((abs(*rnorm)));
    } else if (*(unsigned char *)bmat == 'I') {
	*rnorm = igraphdnrm2_(n, &resid[1], &c__1);
    }

/*        %-----------------------------------------------------------%   
          | STEP 5: Re-orthogonalization / Iterative refinement phase |   
          | Maximum NITER_ITREF tries.                                |   
          |                                                           |   
          |          s      = V_{j}^T * B * r_{j}                     |   
          |          r_{j}  = r_{j} - V_{j}*s                         |   
          |          alphaj = alphaj + s_{j}                          |   
          |                                                           |   
          | The stopping criteria used for iterative refinement is    |   
          | discussed in Parlett's book SEP, page 107 and in Gragg &  |   
          | Reichel ACM TOMS paper; Algorithm 686, Dec. 1990.         |   
          | Determine if we need to correct the residual. The goal is |   
          | to enforce ||v(:,1:j)^T * r_{j}|| .le. eps * || r_{j} ||  |   
          %-----------------------------------------------------------% */

    if (*rnorm > wnorm * .717f) {
	goto L100;
    }
    ++nrorth;

/*        %---------------------------------------------------%   
          | Enter the Iterative refinement phase. If further  |   
          | refinement is necessary, loop back here. The loop |   
          | variable is ITER. Perform a step of Classical     |   
          | Gram-Schmidt using all the Arnoldi vectors V_{j}  |   
          %---------------------------------------------------% */

L80:

    if (msglvl > 2) {
	xtemp[0] = wnorm;
	xtemp[1] = *rnorm;
	igraphdvout_(&logfil, &c__2, xtemp, &ndigit, "_saitr: re-orthonalization ;"
		" wnorm and rnorm are", (ftnlen)48);
    }

/*        %----------------------------------------------------%   
          | Compute V_{j}^T * B * r_{j}.                       |   
          | WORKD(IRJ:IRJ+J-1) = v(:,1:J)'*WORKD(IPJ:IPJ+N-1). |   
          %----------------------------------------------------% */

    igraphdgemv_("T", n, &j, &c_b24, &v[v_offset], ldv, &workd[ipj], &c__1, &c_b49, 
	    &workd[irj], &c__1);

/*        %----------------------------------------------%   
          | Compute the correction to the residual:      |   
          | r_{j} = r_{j} - V_{j} * WORKD(IRJ:IRJ+J-1).  |   
          | The correction to H is v(:,1:J)*H(1:J,1:J) + |   
          | v(:,1:J)*WORKD(IRJ:IRJ+J-1)*e'_j, but only   |   
          | H(j,j) is updated.                           |   
          %----------------------------------------------% */

    igraphdgemv_("N", n, &j, &c_b57, &v[v_offset], ldv, &workd[irj], &c__1, &c_b24, 
	    &resid[1], &c__1);

    if (j == 1 || rstart) {
	h__[j + h_dim1] = 0.;
    }
    h__[j + (h_dim1 << 1)] += workd[irj + j - 1];

    orth2 = TRUE_;
    igraphsecond_(&t2);
    if (*(unsigned char *)bmat == 'G') {
	++nbx;
	igraphdcopy_(n, &resid[1], &c__1, &workd[irj], &c__1);
	ipntr[1] = irj;
	ipntr[2] = ipj;
	*ido = 2;

/*           %-----------------------------------%   
             | Exit in order to compute B*r_{j}. |   
             | r_{j} is the corrected residual.  |   
             %-----------------------------------% */

	goto L9000;
    } else if (*(unsigned char *)bmat == 'I') {
	igraphdcopy_(n, &resid[1], &c__1, &workd[ipj], &c__1);
    }
L90:

/*        %---------------------------------------------------%   
          | Back from reverse communication if ORTH2 = .true. |   
          %---------------------------------------------------% */

    if (*(unsigned char *)bmat == 'G') {
	igraphsecond_(&t3);
	tmvbx += t3 - t2;
    }

/*        %-----------------------------------------------------%   
          | Compute the B-norm of the corrected residual r_{j}. |   
          %-----------------------------------------------------% */

    if (*(unsigned char *)bmat == 'G') {
	rnorm1 = igraphddot_(n, &resid[1], &c__1, &workd[ipj], &c__1);
	rnorm1 = sqrt((abs(rnorm1)));
    } else if (*(unsigned char *)bmat == 'I') {
	rnorm1 = igraphdnrm2_(n, &resid[1], &c__1);
    }

    if (msglvl > 0 && iter > 0) {
	igraphivout_(&logfil, &c__1, &j, &ndigit, "_saitr: Iterative refinement fo"
		"r Arnoldi residual", (ftnlen)49);
	if (msglvl > 2) {
	    xtemp[0] = *rnorm;
	    xtemp[1] = rnorm1;
	    igraphdvout_(&logfil, &c__2, xtemp, &ndigit, "_saitr: iterative refine"
		    "ment ; rnorm and rnorm1 are", (ftnlen)51);
	}
    }

/*        %-----------------------------------------%   
          | Determine if we need to perform another |   
          | step of re-orthogonalization.           |   
          %-----------------------------------------% */

    if (rnorm1 > *rnorm * .717f) {

/*           %--------------------------------%   
             | No need for further refinement |   
             %--------------------------------% */

	*rnorm = rnorm1;

    } else {

/*           %-------------------------------------------%   
             | Another step of iterative refinement step |   
             | is required. NITREF is used by stat.h     |   
             %-------------------------------------------% */

	++nitref;
	*rnorm = rnorm1;
	++iter;
	if (iter <= 1) {
	    goto L80;
	}

/*           %-------------------------------------------------%   
             | Otherwise RESID is numerically in the span of V |   
             %-------------------------------------------------% */

	i__1 = *n;
	for (jj = 1; jj <= i__1; ++jj) {
	    resid[jj] = 0.;
/* L95: */
	}
	*rnorm = 0.;
    }

/*        %----------------------------------------------%   
          | Branch here directly if iterative refinement |   
          | wasn't necessary or after at most NITER_REF  |   
          | steps of iterative refinement.               |   
          %----------------------------------------------% */

L100:

    rstart = FALSE_;
    orth2 = FALSE_;

    igraphsecond_(&t5);
    titref += t5 - t4;

/*        %----------------------------------------------------------%   
          | Make sure the last off-diagonal element is non negative  |   
          | If not perform a similarity transformation on H(1:j,1:j) |   
          | and scale v(:,j) by -1.                                  |   
          %----------------------------------------------------------% */

    if (h__[j + h_dim1] < 0.) {
	h__[j + h_dim1] = -h__[j + h_dim1];
	if (j < *k + *np) {
	    igraphdscal_(n, &c_b57, &v[(j + 1) * v_dim1 + 1], &c__1);
	} else {
	    igraphdscal_(n, &c_b57, &resid[1], &c__1);
	}
    }

/*        %------------------------------------%   
          | STEP 6: Update  j = j+1;  Continue |   
          %------------------------------------% */

    ++j;
    if (j > *k + *np) {
	igraphsecond_(&t1);
	tsaitr += t1 - t0;
	*ido = 99;

	if (msglvl > 1) {
	    i__1 = *k + *np;
	    igraphdvout_(&logfil, &i__1, &h__[(h_dim1 << 1) + 1], &ndigit, "_saitr"
		    ": main diagonal of matrix H of step K+NP.", (ftnlen)47);
	    if (*k + *np > 1) {
		i__1 = *k + *np - 1;
		igraphdvout_(&logfil, &i__1, &h__[h_dim1 + 2], &ndigit, "_saitr: s"
			"ub diagonal of matrix H of step K+NP.", (ftnlen)46);
	    }
	}

	goto L9000;
    }

/*        %--------------------------------------------------------%   
          | Loop back to extend the factorization by another step. |   
          %--------------------------------------------------------% */

    goto L1000;

/*     %---------------------------------------------------------------%   
       |                                                               |   
       |  E N D     O F     M A I N     I T E R A T I O N     L O O P  |   
       |                                                               |   
       %---------------------------------------------------------------% */

L9000:
    return 0;

/*     %---------------%   
       | End of dsaitr |   
       %---------------% */

} /* igraphdsaitr_ */

