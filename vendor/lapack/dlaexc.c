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

static integer c__1 = 1;
static integer c__4 = 4;
static logical c_false = FALSE_;
static integer c_n1 = -1;
static integer c__2 = 2;
static integer c__3 = 3;

/* > \brief \b DLAEXC swaps adjacent diagonal blocks of a real upper quasi-triangular matrix in Schur canonica
l form, by an orthogonal similarity transformation.   

    =========== DOCUMENTATION ===========   

   Online html documentation available at   
              http://www.netlib.org/lapack/explore-html/   

   > \htmlonly   
   > Download DLAEXC + dependencies   
   > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlaexc.
f">   
   > [TGZ]</a>   
   > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlaexc.
f">   
   > [ZIP]</a>   
   > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlaexc.
f">   
   > [TXT]</a>   
   > \endhtmlonly   

    Definition:   
    ===========   

         SUBROUTINE DLAEXC( WANTQ, N, T, LDT, Q, LDQ, J1, N1, N2, WORK,   
                            INFO )   

         LOGICAL            WANTQ   
         INTEGER            INFO, J1, LDQ, LDT, N, N1, N2   
         DOUBLE PRECISION   Q( LDQ, * ), T( LDT, * ), WORK( * )   


   > \par Purpose:   
    =============   
   >   
   > \verbatim   
   >   
   > DLAEXC swaps adjacent diagonal blocks T11 and T22 of order 1 or 2 in   
   > an upper quasi-triangular matrix T by an orthogonal similarity   
   > transformation.   
   >   
   > T must be in Schur canonical form, that is, block upper triangular   
   > with 1-by-1 and 2-by-2 diagonal blocks; each 2-by-2 diagonal block   
   > has its diagonal elemnts equal and its off-diagonal elements of   
   > opposite sign.   
   > \endverbatim   

    Arguments:   
    ==========   

   > \param[in] WANTQ   
   > \verbatim   
   >          WANTQ is LOGICAL   
   >          = .TRUE. : accumulate the transformation in the matrix Q;   
   >          = .FALSE.: do not accumulate the transformation.   
   > \endverbatim   
   >   
   > \param[in] N   
   > \verbatim   
   >          N is INTEGER   
   >          The order of the matrix T. N >= 0.   
   > \endverbatim   
   >   
   > \param[in,out] T   
   > \verbatim   
   >          T is DOUBLE PRECISION array, dimension (LDT,N)   
   >          On entry, the upper quasi-triangular matrix T, in Schur   
   >          canonical form.   
   >          On exit, the updated matrix T, again in Schur canonical form.   
   > \endverbatim   
   >   
   > \param[in] LDT   
   > \verbatim   
   >          LDT is INTEGER   
   >          The leading dimension of the array T. LDT >= max(1,N).   
   > \endverbatim   
   >   
   > \param[in,out] Q   
   > \verbatim   
   >          Q is DOUBLE PRECISION array, dimension (LDQ,N)   
   >          On entry, if WANTQ is .TRUE., the orthogonal matrix Q.   
   >          On exit, if WANTQ is .TRUE., the updated matrix Q.   
   >          If WANTQ is .FALSE., Q is not referenced.   
   > \endverbatim   
   >   
   > \param[in] LDQ   
   > \verbatim   
   >          LDQ is INTEGER   
   >          The leading dimension of the array Q.   
   >          LDQ >= 1; and if WANTQ is .TRUE., LDQ >= N.   
   > \endverbatim   
   >   
   > \param[in] J1   
   > \verbatim   
   >          J1 is INTEGER   
   >          The index of the first row of the first block T11.   
   > \endverbatim   
   >   
   > \param[in] N1   
   > \verbatim   
   >          N1 is INTEGER   
   >          The order of the first block T11. N1 = 0, 1 or 2.   
   > \endverbatim   
   >   
   > \param[in] N2   
   > \verbatim   
   >          N2 is INTEGER   
   >          The order of the second block T22. N2 = 0, 1 or 2.   
   > \endverbatim   
   >   
   > \param[out] WORK   
   > \verbatim   
   >          WORK is DOUBLE PRECISION array, dimension (N)   
   > \endverbatim   
   >   
   > \param[out] INFO   
   > \verbatim   
   >          INFO is INTEGER   
   >          = 0: successful exit   
   >          = 1: the transformed matrix T would be too far from Schur   
   >               form; the blocks are not swapped and T and Q are   
   >               unchanged.   
   > \endverbatim   

    Authors:   
    ========   

   > \author Univ. of Tennessee   
   > \author Univ. of California Berkeley   
   > \author Univ. of Colorado Denver   
   > \author NAG Ltd.   

   > \date September 2012   

   > \ingroup doubleOTHERauxiliary   

    =====================================================================   
   Subroutine */ int igraphdlaexc_(logical *wantq, integer *n, doublereal *t, 
	integer *ldt, doublereal *q, integer *ldq, integer *j1, integer *n1, 
	integer *n2, doublereal *work, integer *info)
{
    /* System generated locals */
    integer q_dim1, q_offset, t_dim1, t_offset, i__1;
    doublereal d__1, d__2, d__3;

    /* Local variables */
    doublereal d__[16]	/* was [4][4] */;
    integer k;
    doublereal u[3], x[4]	/* was [2][2] */;
    integer j2, j3, j4;
    doublereal u1[3], u2[3];
    integer nd;
    doublereal cs, t11, t22, t33, sn, wi1, wi2, wr1, wr2, eps, tau, tau1, 
	    tau2;
    integer ierr;
    doublereal temp;
    extern /* Subroutine */ int igraphdrot_(integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *);
    doublereal scale, dnorm, xnorm;
    extern /* Subroutine */ int igraphdlanv2_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *), igraphdlasy2_(
	    logical *, logical *, integer *, integer *, integer *, doublereal 
	    *, integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *);
    extern doublereal igraphdlamch_(char *), igraphdlange_(char *, integer *, 
	    integer *, doublereal *, integer *, doublereal *);
    extern /* Subroutine */ int igraphdlarfg_(integer *, doublereal *, doublereal *,
	     integer *, doublereal *), igraphdlacpy_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *), 
	    igraphdlartg_(doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *), igraphdlarfx_(char *, integer *, integer *, doublereal *,
	     doublereal *, doublereal *, integer *, doublereal *);
    doublereal thresh, smlnum;


/*  -- LAPACK auxiliary routine (version 3.4.2) --   
    -- LAPACK is a software package provided by Univ. of Tennessee,    --   
    -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--   
       September 2012   


    =====================================================================   


       Parameter adjustments */
    t_dim1 = *ldt;
    t_offset = 1 + t_dim1;
    t -= t_offset;
    q_dim1 = *ldq;
    q_offset = 1 + q_dim1;
    q -= q_offset;
    --work;

    /* Function Body */
    *info = 0;

/*     Quick return if possible */

    if (*n == 0 || *n1 == 0 || *n2 == 0) {
	return 0;
    }
    if (*j1 + *n1 > *n) {
	return 0;
    }

    j2 = *j1 + 1;
    j3 = *j1 + 2;
    j4 = *j1 + 3;

    if (*n1 == 1 && *n2 == 1) {

/*        Swap two 1-by-1 blocks. */

	t11 = t[*j1 + *j1 * t_dim1];
	t22 = t[j2 + j2 * t_dim1];

/*        Determine the transformation to perform the interchange. */

	d__1 = t22 - t11;
	igraphdlartg_(&t[*j1 + j2 * t_dim1], &d__1, &cs, &sn, &temp);

/*        Apply transformation to the matrix T. */

	if (j3 <= *n) {
	    i__1 = *n - *j1 - 1;
	    igraphdrot_(&i__1, &t[*j1 + j3 * t_dim1], ldt, &t[j2 + j3 * t_dim1], 
		    ldt, &cs, &sn);
	}
	i__1 = *j1 - 1;
	igraphdrot_(&i__1, &t[*j1 * t_dim1 + 1], &c__1, &t[j2 * t_dim1 + 1], &c__1, 
		&cs, &sn);

	t[*j1 + *j1 * t_dim1] = t22;
	t[j2 + j2 * t_dim1] = t11;

	if (*wantq) {

/*           Accumulate transformation in the matrix Q. */

	    igraphdrot_(n, &q[*j1 * q_dim1 + 1], &c__1, &q[j2 * q_dim1 + 1], &c__1, 
		    &cs, &sn);
	}

    } else {

/*        Swapping involves at least one 2-by-2 block.   

          Copy the diagonal block of order N1+N2 to the local array D   
          and compute its norm. */

	nd = *n1 + *n2;
	igraphdlacpy_("Full", &nd, &nd, &t[*j1 + *j1 * t_dim1], ldt, d__, &c__4);
	dnorm = igraphdlange_("Max", &nd, &nd, d__, &c__4, &work[1]);

/*        Compute machine-dependent threshold for test for accepting   
          swap. */

	eps = igraphdlamch_("P");
	smlnum = igraphdlamch_("S") / eps;
/* Computing MAX */
	d__1 = eps * 10. * dnorm;
	thresh = max(d__1,smlnum);

/*        Solve T11*X - X*T22 = scale*T12 for X. */

	igraphdlasy2_(&c_false, &c_false, &c_n1, n1, n2, d__, &c__4, &d__[*n1 + 1 + 
		(*n1 + 1 << 2) - 5], &c__4, &d__[(*n1 + 1 << 2) - 4], &c__4, &
		scale, x, &c__2, &xnorm, &ierr);

/*        Swap the adjacent diagonal blocks. */

	k = *n1 + *n1 + *n2 - 3;
	switch (k) {
	    case 1:  goto L10;
	    case 2:  goto L20;
	    case 3:  goto L30;
	}

L10:

/*        N1 = 1, N2 = 2: generate elementary reflector H so that:   

          ( scale, X11, X12 ) H = ( 0, 0, * ) */

	u[0] = scale;
	u[1] = x[0];
	u[2] = x[2];
	igraphdlarfg_(&c__3, &u[2], u, &c__1, &tau);
	u[2] = 1.;
	t11 = t[*j1 + *j1 * t_dim1];

/*        Perform swap provisionally on diagonal block in D. */

	igraphdlarfx_("L", &c__3, &c__3, u, &tau, d__, &c__4, &work[1]);
	igraphdlarfx_("R", &c__3, &c__3, u, &tau, d__, &c__4, &work[1]);

/*        Test whether to reject swap.   

   Computing MAX */
	d__2 = abs(d__[2]), d__3 = abs(d__[6]), d__2 = max(d__2,d__3), d__3 = 
		(d__1 = d__[10] - t11, abs(d__1));
	if (max(d__2,d__3) > thresh) {
	    goto L50;
	}

/*        Accept swap: apply transformation to the entire matrix T. */

	i__1 = *n - *j1 + 1;
	igraphdlarfx_("L", &c__3, &i__1, u, &tau, &t[*j1 + *j1 * t_dim1], ldt, &
		work[1]);
	igraphdlarfx_("R", &j2, &c__3, u, &tau, &t[*j1 * t_dim1 + 1], ldt, &work[1]);

	t[j3 + *j1 * t_dim1] = 0.;
	t[j3 + j2 * t_dim1] = 0.;
	t[j3 + j3 * t_dim1] = t11;

	if (*wantq) {

/*           Accumulate transformation in the matrix Q. */

	    igraphdlarfx_("R", n, &c__3, u, &tau, &q[*j1 * q_dim1 + 1], ldq, &work[
		    1]);
	}
	goto L40;

L20:

/*        N1 = 2, N2 = 1: generate elementary reflector H so that:   

          H (  -X11 ) = ( * )   
            (  -X21 ) = ( 0 )   
            ( scale ) = ( 0 ) */

	u[0] = -x[0];
	u[1] = -x[1];
	u[2] = scale;
	igraphdlarfg_(&c__3, u, &u[1], &c__1, &tau);
	u[0] = 1.;
	t33 = t[j3 + j3 * t_dim1];

/*        Perform swap provisionally on diagonal block in D. */

	igraphdlarfx_("L", &c__3, &c__3, u, &tau, d__, &c__4, &work[1]);
	igraphdlarfx_("R", &c__3, &c__3, u, &tau, d__, &c__4, &work[1]);

/*        Test whether to reject swap.   

   Computing MAX */
	d__2 = abs(d__[1]), d__3 = abs(d__[2]), d__2 = max(d__2,d__3), d__3 = 
		(d__1 = d__[0] - t33, abs(d__1));
	if (max(d__2,d__3) > thresh) {
	    goto L50;
	}

/*        Accept swap: apply transformation to the entire matrix T. */

	igraphdlarfx_("R", &j3, &c__3, u, &tau, &t[*j1 * t_dim1 + 1], ldt, &work[1]);
	i__1 = *n - *j1;
	igraphdlarfx_("L", &c__3, &i__1, u, &tau, &t[*j1 + j2 * t_dim1], ldt, &work[
		1]);

	t[*j1 + *j1 * t_dim1] = t33;
	t[j2 + *j1 * t_dim1] = 0.;
	t[j3 + *j1 * t_dim1] = 0.;

	if (*wantq) {

/*           Accumulate transformation in the matrix Q. */

	    igraphdlarfx_("R", n, &c__3, u, &tau, &q[*j1 * q_dim1 + 1], ldq, &work[
		    1]);
	}
	goto L40;

L30:

/*        N1 = 2, N2 = 2: generate elementary reflectors H(1) and H(2) so   
          that:   

          H(2) H(1) (  -X11  -X12 ) = (  *  * )   
                    (  -X21  -X22 )   (  0  * )   
                    ( scale    0  )   (  0  0 )   
                    (    0  scale )   (  0  0 ) */

	u1[0] = -x[0];
	u1[1] = -x[1];
	u1[2] = scale;
	igraphdlarfg_(&c__3, u1, &u1[1], &c__1, &tau1);
	u1[0] = 1.;

	temp = -tau1 * (x[2] + u1[1] * x[3]);
	u2[0] = -temp * u1[1] - x[3];
	u2[1] = -temp * u1[2];
	u2[2] = scale;
	igraphdlarfg_(&c__3, u2, &u2[1], &c__1, &tau2);
	u2[0] = 1.;

/*        Perform swap provisionally on diagonal block in D. */

	igraphdlarfx_("L", &c__3, &c__4, u1, &tau1, d__, &c__4, &work[1])
		;
	igraphdlarfx_("R", &c__4, &c__3, u1, &tau1, d__, &c__4, &work[1])
		;
	igraphdlarfx_("L", &c__3, &c__4, u2, &tau2, &d__[1], &c__4, &work[1]);
	igraphdlarfx_("R", &c__4, &c__3, u2, &tau2, &d__[4], &c__4, &work[1]);

/*        Test whether to reject swap.   

   Computing MAX */
	d__1 = abs(d__[2]), d__2 = abs(d__[6]), d__1 = max(d__1,d__2), d__2 = 
		abs(d__[3]), d__1 = max(d__1,d__2), d__2 = abs(d__[7]);
	if (max(d__1,d__2) > thresh) {
	    goto L50;
	}

/*        Accept swap: apply transformation to the entire matrix T. */

	i__1 = *n - *j1 + 1;
	igraphdlarfx_("L", &c__3, &i__1, u1, &tau1, &t[*j1 + *j1 * t_dim1], ldt, &
		work[1]);
	igraphdlarfx_("R", &j4, &c__3, u1, &tau1, &t[*j1 * t_dim1 + 1], ldt, &work[
		1]);
	i__1 = *n - *j1 + 1;
	igraphdlarfx_("L", &c__3, &i__1, u2, &tau2, &t[j2 + *j1 * t_dim1], ldt, &
		work[1]);
	igraphdlarfx_("R", &j4, &c__3, u2, &tau2, &t[j2 * t_dim1 + 1], ldt, &work[1]
		);

	t[j3 + *j1 * t_dim1] = 0.;
	t[j3 + j2 * t_dim1] = 0.;
	t[j4 + *j1 * t_dim1] = 0.;
	t[j4 + j2 * t_dim1] = 0.;

	if (*wantq) {

/*           Accumulate transformation in the matrix Q. */

	    igraphdlarfx_("R", n, &c__3, u1, &tau1, &q[*j1 * q_dim1 + 1], ldq, &
		    work[1]);
	    igraphdlarfx_("R", n, &c__3, u2, &tau2, &q[j2 * q_dim1 + 1], ldq, &work[
		    1]);
	}

L40:

	if (*n2 == 2) {

/*           Standardize new 2-by-2 block T11 */

	    igraphdlanv2_(&t[*j1 + *j1 * t_dim1], &t[*j1 + j2 * t_dim1], &t[j2 + *
		    j1 * t_dim1], &t[j2 + j2 * t_dim1], &wr1, &wi1, &wr2, &
		    wi2, &cs, &sn);
	    i__1 = *n - *j1 - 1;
	    igraphdrot_(&i__1, &t[*j1 + (*j1 + 2) * t_dim1], ldt, &t[j2 + (*j1 + 2) 
		    * t_dim1], ldt, &cs, &sn);
	    i__1 = *j1 - 1;
	    igraphdrot_(&i__1, &t[*j1 * t_dim1 + 1], &c__1, &t[j2 * t_dim1 + 1], &
		    c__1, &cs, &sn);
	    if (*wantq) {
		igraphdrot_(n, &q[*j1 * q_dim1 + 1], &c__1, &q[j2 * q_dim1 + 1], &
			c__1, &cs, &sn);
	    }
	}

	if (*n1 == 2) {

/*           Standardize new 2-by-2 block T22 */

	    j3 = *j1 + *n2;
	    j4 = j3 + 1;
	    igraphdlanv2_(&t[j3 + j3 * t_dim1], &t[j3 + j4 * t_dim1], &t[j4 + j3 * 
		    t_dim1], &t[j4 + j4 * t_dim1], &wr1, &wi1, &wr2, &wi2, &
		    cs, &sn);
	    if (j3 + 2 <= *n) {
		i__1 = *n - j3 - 1;
		igraphdrot_(&i__1, &t[j3 + (j3 + 2) * t_dim1], ldt, &t[j4 + (j3 + 2)
			 * t_dim1], ldt, &cs, &sn);
	    }
	    i__1 = j3 - 1;
	    igraphdrot_(&i__1, &t[j3 * t_dim1 + 1], &c__1, &t[j4 * t_dim1 + 1], &
		    c__1, &cs, &sn);
	    if (*wantq) {
		igraphdrot_(n, &q[j3 * q_dim1 + 1], &c__1, &q[j4 * q_dim1 + 1], &
			c__1, &cs, &sn);
	    }
	}

    }
    return 0;

/*     Exit with INFO = 1 if swap was rejected. */

L50:
    *info = 1;
    return 0;

/*     End of DLAEXC */

} /* igraphdlaexc_ */

