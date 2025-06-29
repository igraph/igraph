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
static integer c_n1 = -1;
static doublereal c_b12 = 0.;
static doublereal c_b13 = 1.;
static logical c_true = TRUE_;

/* > \brief \b DLAQR2 performs the orthogonal similarity transformation of a Hessenberg matrix to detect and d
eflate fully converged eigenvalues from a trailing principal submatrix (aggressive early deflation). 
  

    =========== DOCUMENTATION ===========   

   Online html documentation available at   
              http://www.netlib.org/lapack/explore-html/   

   > \htmlonly   
   > Download DLAQR2 + dependencies   
   > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlaqr2.
f">   
   > [TGZ]</a>   
   > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlaqr2.
f">   
   > [ZIP]</a>   
   > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlaqr2.
f">   
   > [TXT]</a>   
   > \endhtmlonly   

    Definition:   
    ===========   

         SUBROUTINE DLAQR2( WANTT, WANTZ, N, KTOP, KBOT, NW, H, LDH, ILOZ,   
                            IHIZ, Z, LDZ, NS, ND, SR, SI, V, LDV, NH, T,   
                            LDT, NV, WV, LDWV, WORK, LWORK )   

         INTEGER            IHIZ, ILOZ, KBOT, KTOP, LDH, LDT, LDV, LDWV,   
        $                   LDZ, LWORK, N, ND, NH, NS, NV, NW   
         LOGICAL            WANTT, WANTZ   
         DOUBLE PRECISION   H( LDH, * ), SI( * ), SR( * ), T( LDT, * ),   
        $                   V( LDV, * ), WORK( * ), WV( LDWV, * ),   
        $                   Z( LDZ, * )   


   > \par Purpose:   
    =============   
   >   
   > \verbatim   
   >   
   >    DLAQR2 is identical to DLAQR3 except that it avoids   
   >    recursion by calling DLAHQR instead of DLAQR4.   
   >   
   >    Aggressive early deflation:   
   >   
   >    This subroutine accepts as input an upper Hessenberg matrix   
   >    H and performs an orthogonal similarity transformation   
   >    designed to detect and deflate fully converged eigenvalues from   
   >    a trailing principal submatrix.  On output H has been over-   
   >    written by a new Hessenberg matrix that is a perturbation of   
   >    an orthogonal similarity transformation of H.  It is to be   
   >    hoped that the final version of H has many zero subdiagonal   
   >    entries.   
   > \endverbatim   

    Arguments:   
    ==========   

   > \param[in] WANTT   
   > \verbatim   
   >          WANTT is LOGICAL   
   >          If .TRUE., then the Hessenberg matrix H is fully updated   
   >          so that the quasi-triangular Schur factor may be   
   >          computed (in cooperation with the calling subroutine).   
   >          If .FALSE., then only enough of H is updated to preserve   
   >          the eigenvalues.   
   > \endverbatim   
   >   
   > \param[in] WANTZ   
   > \verbatim   
   >          WANTZ is LOGICAL   
   >          If .TRUE., then the orthogonal matrix Z is updated so   
   >          so that the orthogonal Schur factor may be computed   
   >          (in cooperation with the calling subroutine).   
   >          If .FALSE., then Z is not referenced.   
   > \endverbatim   
   >   
   > \param[in] N   
   > \verbatim   
   >          N is INTEGER   
   >          The order of the matrix H and (if WANTZ is .TRUE.) the   
   >          order of the orthogonal matrix Z.   
   > \endverbatim   
   >   
   > \param[in] KTOP   
   > \verbatim   
   >          KTOP is INTEGER   
   >          It is assumed that either KTOP = 1 or H(KTOP,KTOP-1)=0.   
   >          KBOT and KTOP together determine an isolated block   
   >          along the diagonal of the Hessenberg matrix.   
   > \endverbatim   
   >   
   > \param[in] KBOT   
   > \verbatim   
   >          KBOT is INTEGER   
   >          It is assumed without a check that either   
   >          KBOT = N or H(KBOT+1,KBOT)=0.  KBOT and KTOP together   
   >          determine an isolated block along the diagonal of the   
   >          Hessenberg matrix.   
   > \endverbatim   
   >   
   > \param[in] NW   
   > \verbatim   
   >          NW is INTEGER   
   >          Deflation window size.  1 .LE. NW .LE. (KBOT-KTOP+1).   
   > \endverbatim   
   >   
   > \param[in,out] H   
   > \verbatim   
   >          H is DOUBLE PRECISION array, dimension (LDH,N)   
   >          On input the initial N-by-N section of H stores the   
   >          Hessenberg matrix undergoing aggressive early deflation.   
   >          On output H has been transformed by an orthogonal   
   >          similarity transformation, perturbed, and the returned   
   >          to Hessenberg form that (it is to be hoped) has some   
   >          zero subdiagonal entries.   
   > \endverbatim   
   >   
   > \param[in] LDH   
   > \verbatim   
   >          LDH is integer   
   >          Leading dimension of H just as declared in the calling   
   >          subroutine.  N .LE. LDH   
   > \endverbatim   
   >   
   > \param[in] ILOZ   
   > \verbatim   
   >          ILOZ is INTEGER   
   > \endverbatim   
   >   
   > \param[in] IHIZ   
   > \verbatim   
   >          IHIZ is INTEGER   
   >          Specify the rows of Z to which transformations must be   
   >          applied if WANTZ is .TRUE.. 1 .LE. ILOZ .LE. IHIZ .LE. N.   
   > \endverbatim   
   >   
   > \param[in,out] Z   
   > \verbatim   
   >          Z is DOUBLE PRECISION array, dimension (LDZ,N)   
   >          IF WANTZ is .TRUE., then on output, the orthogonal   
   >          similarity transformation mentioned above has been   
   >          accumulated into Z(ILOZ:IHIZ,ILO:IHI) from the right.   
   >          If WANTZ is .FALSE., then Z is unreferenced.   
   > \endverbatim   
   >   
   > \param[in] LDZ   
   > \verbatim   
   >          LDZ is integer   
   >          The leading dimension of Z just as declared in the   
   >          calling subroutine.  1 .LE. LDZ.   
   > \endverbatim   
   >   
   > \param[out] NS   
   > \verbatim   
   >          NS is integer   
   >          The number of unconverged (ie approximate) eigenvalues   
   >          returned in SR and SI that may be used as shifts by the   
   >          calling subroutine.   
   > \endverbatim   
   >   
   > \param[out] ND   
   > \verbatim   
   >          ND is integer   
   >          The number of converged eigenvalues uncovered by this   
   >          subroutine.   
   > \endverbatim   
   >   
   > \param[out] SR   
   > \verbatim   
   >          SR is DOUBLE PRECISION array, dimension (KBOT)   
   > \endverbatim   
   >   
   > \param[out] SI   
   > \verbatim   
   >          SI is DOUBLE PRECISION array, dimension (KBOT)   
   >          On output, the real and imaginary parts of approximate   
   >          eigenvalues that may be used for shifts are stored in   
   >          SR(KBOT-ND-NS+1) through SR(KBOT-ND) and   
   >          SI(KBOT-ND-NS+1) through SI(KBOT-ND), respectively.   
   >          The real and imaginary parts of converged eigenvalues   
   >          are stored in SR(KBOT-ND+1) through SR(KBOT) and   
   >          SI(KBOT-ND+1) through SI(KBOT), respectively.   
   > \endverbatim   
   >   
   > \param[out] V   
   > \verbatim   
   >          V is DOUBLE PRECISION array, dimension (LDV,NW)   
   >          An NW-by-NW work array.   
   > \endverbatim   
   >   
   > \param[in] LDV   
   > \verbatim   
   >          LDV is integer scalar   
   >          The leading dimension of V just as declared in the   
   >          calling subroutine.  NW .LE. LDV   
   > \endverbatim   
   >   
   > \param[in] NH   
   > \verbatim   
   >          NH is integer scalar   
   >          The number of columns of T.  NH.GE.NW.   
   > \endverbatim   
   >   
   > \param[out] T   
   > \verbatim   
   >          T is DOUBLE PRECISION array, dimension (LDT,NW)   
   > \endverbatim   
   >   
   > \param[in] LDT   
   > \verbatim   
   >          LDT is integer   
   >          The leading dimension of T just as declared in the   
   >          calling subroutine.  NW .LE. LDT   
   > \endverbatim   
   >   
   > \param[in] NV   
   > \verbatim   
   >          NV is integer   
   >          The number of rows of work array WV available for   
   >          workspace.  NV.GE.NW.   
   > \endverbatim   
   >   
   > \param[out] WV   
   > \verbatim   
   >          WV is DOUBLE PRECISION array, dimension (LDWV,NW)   
   > \endverbatim   
   >   
   > \param[in] LDWV   
   > \verbatim   
   >          LDWV is integer   
   >          The leading dimension of W just as declared in the   
   >          calling subroutine.  NW .LE. LDV   
   > \endverbatim   
   >   
   > \param[out] WORK   
   > \verbatim   
   >          WORK is DOUBLE PRECISION array, dimension (LWORK)   
   >          On exit, WORK(1) is set to an estimate of the optimal value   
   >          of LWORK for the given values of N, NW, KTOP and KBOT.   
   > \endverbatim   
   >   
   > \param[in] LWORK   
   > \verbatim   
   >          LWORK is integer   
   >          The dimension of the work array WORK.  LWORK = 2*NW   
   >          suffices, but greater efficiency may result from larger   
   >          values of LWORK.   
   >   
   >          If LWORK = -1, then a workspace query is assumed; DLAQR2   
   >          only estimates the optimal workspace size for the given   
   >          values of N, NW, KTOP and KBOT.  The estimate is returned   
   >          in WORK(1).  No error message related to LWORK is issued   
   >          by XERBLA.  Neither H nor Z are accessed.   
   > \endverbatim   

    Authors:   
    ========   

   > \author Univ. of Tennessee   
   > \author Univ. of California Berkeley   
   > \author Univ. of Colorado Denver   
   > \author NAG Ltd.   

   > \date September 2012   

   > \ingroup doubleOTHERauxiliary   

   > \par Contributors:   
    ==================   
   >   
   >       Karen Braman and Ralph Byers, Department of Mathematics,   
   >       University of Kansas, USA   
   >   
    =====================================================================   
   Subroutine */ int igraphdlaqr2_(logical *wantt, logical *wantz, integer *n, 
	integer *ktop, integer *kbot, integer *nw, doublereal *h__, integer *
	ldh, integer *iloz, integer *ihiz, doublereal *z__, integer *ldz, 
	integer *ns, integer *nd, doublereal *sr, doublereal *si, doublereal *
	v, integer *ldv, integer *nh, doublereal *t, integer *ldt, integer *
	nv, doublereal *wv, integer *ldwv, doublereal *work, integer *lwork)
{
    /* System generated locals */
    integer h_dim1, h_offset, t_dim1, t_offset, v_dim1, v_offset, wv_dim1, 
	    wv_offset, z_dim1, z_offset, i__1, i__2, i__3, i__4;
    doublereal d__1, d__2, d__3, d__4, d__5, d__6;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    integer i__, j, k;
    doublereal s, aa, bb, cc, dd, cs, sn;
    integer jw;
    doublereal evi, evk, foo;
    integer kln;
    doublereal tau, ulp;
    integer lwk1, lwk2;
    doublereal beta;
    integer kend, kcol, info, ifst, ilst, ltop, krow;
    extern /* Subroutine */ int igraphdlarf_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *), igraphdgemm_(char *, char *, integer *, integer *
	    , integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *);
    logical bulge;
    extern /* Subroutine */ int igraphdcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    integer infqr, kwtop;
    extern /* Subroutine */ int igraphdlanv2_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *), igraphdlabad_(
	    doublereal *, doublereal *);
    extern doublereal igraphdlamch_(char *);
    extern /* Subroutine */ int igraphdgehrd_(integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    integer *), igraphdlarfg_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *), igraphdlahqr_(logical *, logical *, integer *,
	     integer *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, integer *, doublereal *, integer *, 
	    integer *), igraphdlacpy_(char *, integer *, integer *, doublereal *, 
	    integer *, doublereal *, integer *);
    doublereal safmin;
    extern /* Subroutine */ int igraphdlaset_(char *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *);
    doublereal safmax;
    extern /* Subroutine */ int igraphdtrexc_(char *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, integer *, integer *, 
	    doublereal *, integer *), igraphdormhr_(char *, char *, integer 
	    *, integer *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    integer *);
    logical sorted;
    doublereal smlnum;
    integer lwkopt;


/*  -- LAPACK auxiliary routine (version 3.4.2) --   
    -- LAPACK is a software package provided by Univ. of Tennessee,    --   
    -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--   
       September 2012   


    ================================================================   

       ==== Estimate optimal workspace. ====   

       Parameter adjustments */
    h_dim1 = *ldh;
    h_offset = 1 + h_dim1;
    h__ -= h_offset;
    z_dim1 = *ldz;
    z_offset = 1 + z_dim1;
    z__ -= z_offset;
    --sr;
    --si;
    v_dim1 = *ldv;
    v_offset = 1 + v_dim1;
    v -= v_offset;
    t_dim1 = *ldt;
    t_offset = 1 + t_dim1;
    t -= t_offset;
    wv_dim1 = *ldwv;
    wv_offset = 1 + wv_dim1;
    wv -= wv_offset;
    --work;

    /* Function Body   
   Computing MIN */
    i__1 = *nw, i__2 = *kbot - *ktop + 1;
    jw = min(i__1,i__2);
    if (jw <= 2) {
	lwkopt = 1;
    } else {

/*        ==== Workspace query call to DGEHRD ==== */

	i__1 = jw - 1;
	igraphdgehrd_(&jw, &c__1, &i__1, &t[t_offset], ldt, &work[1], &work[1], &
		c_n1, &info);
	lwk1 = (integer) work[1];

/*        ==== Workspace query call to DORMHR ==== */

	i__1 = jw - 1;
	igraphdormhr_("R", "N", &jw, &jw, &c__1, &i__1, &t[t_offset], ldt, &work[1],
		 &v[v_offset], ldv, &work[1], &c_n1, &info);
	lwk2 = (integer) work[1];

/*        ==== Optimal workspace ==== */

	lwkopt = jw + max(lwk1,lwk2);
    }

/*     ==== Quick return in case of workspace query. ==== */

    if (*lwork == -1) {
	work[1] = (doublereal) lwkopt;
	return 0;
    }

/*     ==== Nothing to do ...   
       ... for an empty active block ... ==== */
    *ns = 0;
    *nd = 0;
    work[1] = 1.;
    if (*ktop > *kbot) {
	return 0;
    }
/*     ... nor for an empty deflation window. ==== */
    if (*nw < 1) {
	return 0;
    }

/*     ==== Machine constants ==== */

    safmin = igraphdlamch_("SAFE MINIMUM");
    safmax = 1. / safmin;
    igraphdlabad_(&safmin, &safmax);
    ulp = igraphdlamch_("PRECISION");
    smlnum = safmin * ((doublereal) (*n) / ulp);

/*     ==== Setup deflation window ====   

   Computing MIN */
    i__1 = *nw, i__2 = *kbot - *ktop + 1;
    jw = min(i__1,i__2);
    kwtop = *kbot - jw + 1;
    if (kwtop == *ktop) {
	s = 0.;
    } else {
	s = h__[kwtop + (kwtop - 1) * h_dim1];
    }

    if (*kbot == kwtop) {

/*        ==== 1-by-1 deflation window: not much to do ==== */

	sr[kwtop] = h__[kwtop + kwtop * h_dim1];
	si[kwtop] = 0.;
	*ns = 1;
	*nd = 0;
/* Computing MAX */
	d__2 = smlnum, d__3 = ulp * (d__1 = h__[kwtop + kwtop * h_dim1], abs(
		d__1));
	if (abs(s) <= max(d__2,d__3)) {
	    *ns = 0;
	    *nd = 1;
	    if (kwtop > *ktop) {
		h__[kwtop + (kwtop - 1) * h_dim1] = 0.;
	    }
	}
	work[1] = 1.;
	return 0;
    }

/*     ==== Convert to spike-triangular form.  (In case of a   
       .    rare QR failure, this routine continues to do   
       .    aggressive early deflation using that part of   
       .    the deflation window that converged using INFQR   
       .    here and there to keep track.) ==== */

    igraphdlacpy_("U", &jw, &jw, &h__[kwtop + kwtop * h_dim1], ldh, &t[t_offset], 
	    ldt);
    i__1 = jw - 1;
    i__2 = *ldh + 1;
    i__3 = *ldt + 1;
    igraphdcopy_(&i__1, &h__[kwtop + 1 + kwtop * h_dim1], &i__2, &t[t_dim1 + 2], &
	    i__3);

    igraphdlaset_("A", &jw, &jw, &c_b12, &c_b13, &v[v_offset], ldv);
    igraphdlahqr_(&c_true, &c_true, &jw, &c__1, &jw, &t[t_offset], ldt, &sr[kwtop], 
	    &si[kwtop], &c__1, &jw, &v[v_offset], ldv, &infqr);

/*     ==== DTREXC needs a clean margin near the diagonal ==== */

    i__1 = jw - 3;
    for (j = 1; j <= i__1; ++j) {
	t[j + 2 + j * t_dim1] = 0.;
	t[j + 3 + j * t_dim1] = 0.;
/* L10: */
    }
    if (jw > 2) {
	t[jw + (jw - 2) * t_dim1] = 0.;
    }

/*     ==== Deflation detection loop ==== */

    *ns = jw;
    ilst = infqr + 1;
L20:
    if (ilst <= *ns) {
	if (*ns == 1) {
	    bulge = FALSE_;
	} else {
	    bulge = t[*ns + (*ns - 1) * t_dim1] != 0.;
	}

/*        ==== Small spike tip test for deflation ==== */

	if (! bulge) {

/*           ==== Real eigenvalue ==== */

	    foo = (d__1 = t[*ns + *ns * t_dim1], abs(d__1));
	    if (foo == 0.) {
		foo = abs(s);
	    }
/* Computing MAX */
	    d__2 = smlnum, d__3 = ulp * foo;
	    if ((d__1 = s * v[*ns * v_dim1 + 1], abs(d__1)) <= max(d__2,d__3))
		     {

/*              ==== Deflatable ==== */

		--(*ns);
	    } else {

/*              ==== Undeflatable.   Move it up out of the way.   
                .    (DTREXC can not fail in this case.) ==== */

		ifst = *ns;
		igraphdtrexc_("V", &jw, &t[t_offset], ldt, &v[v_offset], ldv, &ifst,
			 &ilst, &work[1], &info);
		++ilst;
	    }
	} else {

/*           ==== Complex conjugate pair ==== */

	    foo = (d__3 = t[*ns + *ns * t_dim1], abs(d__3)) + sqrt((d__1 = t[*
		    ns + (*ns - 1) * t_dim1], abs(d__1))) * sqrt((d__2 = t[*
		    ns - 1 + *ns * t_dim1], abs(d__2)));
	    if (foo == 0.) {
		foo = abs(s);
	    }
/* Computing MAX */
	    d__3 = (d__1 = s * v[*ns * v_dim1 + 1], abs(d__1)), d__4 = (d__2 =
		     s * v[(*ns - 1) * v_dim1 + 1], abs(d__2));
/* Computing MAX */
	    d__5 = smlnum, d__6 = ulp * foo;
	    if (max(d__3,d__4) <= max(d__5,d__6)) {

/*              ==== Deflatable ==== */

		*ns += -2;
	    } else {

/*              ==== Undeflatable. Move them up out of the way.   
                .    Fortunately, DTREXC does the right thing with   
                .    ILST in case of a rare exchange failure. ==== */

		ifst = *ns;
		igraphdtrexc_("V", &jw, &t[t_offset], ldt, &v[v_offset], ldv, &ifst,
			 &ilst, &work[1], &info);
		ilst += 2;
	    }
	}

/*        ==== End deflation detection loop ==== */

	goto L20;
    }

/*        ==== Return to Hessenberg form ==== */

    if (*ns == 0) {
	s = 0.;
    }

    if (*ns < jw) {

/*        ==== sorting diagonal blocks of T improves accuracy for   
          .    graded matrices.  Bubble sort deals well with   
          .    exchange failures. ==== */

	sorted = FALSE_;
	i__ = *ns + 1;
L30:
	if (sorted) {
	    goto L50;
	}
	sorted = TRUE_;

	kend = i__ - 1;
	i__ = infqr + 1;
	if (i__ == *ns) {
	    k = i__ + 1;
	} else if (t[i__ + 1 + i__ * t_dim1] == 0.) {
	    k = i__ + 1;
	} else {
	    k = i__ + 2;
	}
L40:
	if (k <= kend) {
	    if (k == i__ + 1) {
		evi = (d__1 = t[i__ + i__ * t_dim1], abs(d__1));
	    } else {
		evi = (d__3 = t[i__ + i__ * t_dim1], abs(d__3)) + sqrt((d__1 =
			 t[i__ + 1 + i__ * t_dim1], abs(d__1))) * sqrt((d__2 =
			 t[i__ + (i__ + 1) * t_dim1], abs(d__2)));
	    }

	    if (k == kend) {
		evk = (d__1 = t[k + k * t_dim1], abs(d__1));
	    } else if (t[k + 1 + k * t_dim1] == 0.) {
		evk = (d__1 = t[k + k * t_dim1], abs(d__1));
	    } else {
		evk = (d__3 = t[k + k * t_dim1], abs(d__3)) + sqrt((d__1 = t[
			k + 1 + k * t_dim1], abs(d__1))) * sqrt((d__2 = t[k + 
			(k + 1) * t_dim1], abs(d__2)));
	    }

	    if (evi >= evk) {
		i__ = k;
	    } else {
		sorted = FALSE_;
		ifst = i__;
		ilst = k;
		igraphdtrexc_("V", &jw, &t[t_offset], ldt, &v[v_offset], ldv, &ifst,
			 &ilst, &work[1], &info);
		if (info == 0) {
		    i__ = ilst;
		} else {
		    i__ = k;
		}
	    }
	    if (i__ == kend) {
		k = i__ + 1;
	    } else if (t[i__ + 1 + i__ * t_dim1] == 0.) {
		k = i__ + 1;
	    } else {
		k = i__ + 2;
	    }
	    goto L40;
	}
	goto L30;
L50:
	;
    }

/*     ==== Restore shift/eigenvalue array from T ==== */

    i__ = jw;
L60:
    if (i__ >= infqr + 1) {
	if (i__ == infqr + 1) {
	    sr[kwtop + i__ - 1] = t[i__ + i__ * t_dim1];
	    si[kwtop + i__ - 1] = 0.;
	    --i__;
	} else if (t[i__ + (i__ - 1) * t_dim1] == 0.) {
	    sr[kwtop + i__ - 1] = t[i__ + i__ * t_dim1];
	    si[kwtop + i__ - 1] = 0.;
	    --i__;
	} else {
	    aa = t[i__ - 1 + (i__ - 1) * t_dim1];
	    cc = t[i__ + (i__ - 1) * t_dim1];
	    bb = t[i__ - 1 + i__ * t_dim1];
	    dd = t[i__ + i__ * t_dim1];
	    igraphdlanv2_(&aa, &bb, &cc, &dd, &sr[kwtop + i__ - 2], &si[kwtop + i__ 
		    - 2], &sr[kwtop + i__ - 1], &si[kwtop + i__ - 1], &cs, &
		    sn);
	    i__ += -2;
	}
	goto L60;
    }

    if (*ns < jw || s == 0.) {
	if (*ns > 1 && s != 0.) {

/*           ==== Reflect spike back into lower triangle ==== */

	    igraphdcopy_(ns, &v[v_offset], ldv, &work[1], &c__1);
	    beta = work[1];
	    igraphdlarfg_(ns, &beta, &work[2], &c__1, &tau);
	    work[1] = 1.;

	    i__1 = jw - 2;
	    i__2 = jw - 2;
	    igraphdlaset_("L", &i__1, &i__2, &c_b12, &c_b12, &t[t_dim1 + 3], ldt);

	    igraphdlarf_("L", ns, &jw, &work[1], &c__1, &tau, &t[t_offset], ldt, &
		    work[jw + 1]);
	    igraphdlarf_("R", ns, ns, &work[1], &c__1, &tau, &t[t_offset], ldt, &
		    work[jw + 1]);
	    igraphdlarf_("R", &jw, ns, &work[1], &c__1, &tau, &v[v_offset], ldv, &
		    work[jw + 1]);

	    i__1 = *lwork - jw;
	    igraphdgehrd_(&jw, &c__1, ns, &t[t_offset], ldt, &work[1], &work[jw + 1]
		    , &i__1, &info);
	}

/*        ==== Copy updated reduced window into place ==== */

	if (kwtop > 1) {
	    h__[kwtop + (kwtop - 1) * h_dim1] = s * v[v_dim1 + 1];
	}
	igraphdlacpy_("U", &jw, &jw, &t[t_offset], ldt, &h__[kwtop + kwtop * h_dim1]
		, ldh);
	i__1 = jw - 1;
	i__2 = *ldt + 1;
	i__3 = *ldh + 1;
	igraphdcopy_(&i__1, &t[t_dim1 + 2], &i__2, &h__[kwtop + 1 + kwtop * h_dim1],
		 &i__3);

/*        ==== Accumulate orthogonal matrix in order update   
          .    H and Z, if requested.  ==== */

	if (*ns > 1 && s != 0.) {
	    i__1 = *lwork - jw;
	    igraphdormhr_("R", "N", &jw, ns, &c__1, ns, &t[t_offset], ldt, &work[1],
		     &v[v_offset], ldv, &work[jw + 1], &i__1, &info);
	}

/*        ==== Update vertical slab in H ==== */

	if (*wantt) {
	    ltop = 1;
	} else {
	    ltop = *ktop;
	}
	i__1 = kwtop - 1;
	i__2 = *nv;
	for (krow = ltop; i__2 < 0 ? krow >= i__1 : krow <= i__1; krow += 
		i__2) {
/* Computing MIN */
	    i__3 = *nv, i__4 = kwtop - krow;
	    kln = min(i__3,i__4);
	    igraphdgemm_("N", "N", &kln, &jw, &jw, &c_b13, &h__[krow + kwtop * 
		    h_dim1], ldh, &v[v_offset], ldv, &c_b12, &wv[wv_offset], 
		    ldwv);
	    igraphdlacpy_("A", &kln, &jw, &wv[wv_offset], ldwv, &h__[krow + kwtop * 
		    h_dim1], ldh);
/* L70: */
	}

/*        ==== Update horizontal slab in H ==== */

	if (*wantt) {
	    i__2 = *n;
	    i__1 = *nh;
	    for (kcol = *kbot + 1; i__1 < 0 ? kcol >= i__2 : kcol <= i__2; 
		    kcol += i__1) {
/* Computing MIN */
		i__3 = *nh, i__4 = *n - kcol + 1;
		kln = min(i__3,i__4);
		igraphdgemm_("C", "N", &jw, &kln, &jw, &c_b13, &v[v_offset], ldv, &
			h__[kwtop + kcol * h_dim1], ldh, &c_b12, &t[t_offset],
			 ldt);
		igraphdlacpy_("A", &jw, &kln, &t[t_offset], ldt, &h__[kwtop + kcol *
			 h_dim1], ldh);
/* L80: */
	    }
	}

/*        ==== Update vertical slab in Z ==== */

	if (*wantz) {
	    i__1 = *ihiz;
	    i__2 = *nv;
	    for (krow = *iloz; i__2 < 0 ? krow >= i__1 : krow <= i__1; krow +=
		     i__2) {
/* Computing MIN */
		i__3 = *nv, i__4 = *ihiz - krow + 1;
		kln = min(i__3,i__4);
		igraphdgemm_("N", "N", &kln, &jw, &jw, &c_b13, &z__[krow + kwtop * 
			z_dim1], ldz, &v[v_offset], ldv, &c_b12, &wv[
			wv_offset], ldwv);
		igraphdlacpy_("A", &kln, &jw, &wv[wv_offset], ldwv, &z__[krow + 
			kwtop * z_dim1], ldz);
/* L90: */
	    }
	}
    }

/*     ==== Return the number of deflations ... ==== */

    *nd = jw - *ns;

/*     ==== ... and the number of shifts. (Subtracting   
       .    INFQR from the spike length takes care   
       .    of the case of a rare QR failure while   
       .    calculating eigenvalues of the deflation   
       .    window.)  ==== */

    *ns -= infqr;

/*      ==== Return optimal workspace. ==== */

    work[1] = (doublereal) lwkopt;

/*     ==== End of DLAQR2 ==== */

    return 0;
} /* igraphdlaqr2_ */

