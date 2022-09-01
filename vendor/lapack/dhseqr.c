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

static doublereal c_b11 = 0.;
static doublereal c_b12 = 1.;
static integer c__12 = 12;
static integer c__2 = 2;
static integer c__49 = 49;

/* > \brief \b DHSEQR   

    =========== DOCUMENTATION ===========   

   Online html documentation available at   
              http://www.netlib.org/lapack/explore-html/   

   > \htmlonly   
   > Download DHSEQR + dependencies   
   > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dhseqr.
f">   
   > [TGZ]</a>   
   > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dhseqr.
f">   
   > [ZIP]</a>   
   > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dhseqr.
f">   
   > [TXT]</a>   
   > \endhtmlonly   

    Definition:   
    ===========   

         SUBROUTINE DHSEQR( JOB, COMPZ, N, ILO, IHI, H, LDH, WR, WI, Z,   
                            LDZ, WORK, LWORK, INFO )   

         INTEGER            IHI, ILO, INFO, LDH, LDZ, LWORK, N   
         CHARACTER          COMPZ, JOB   
         DOUBLE PRECISION   H( LDH, * ), WI( * ), WORK( * ), WR( * ),   
        $                   Z( LDZ, * )   


   > \par Purpose:   
    =============   
   >   
   > \verbatim   
   >   
   >    DHSEQR computes the eigenvalues of a Hessenberg matrix H   
   >    and, optionally, the matrices T and Z from the Schur decomposition   
   >    H = Z T Z**T, where T is an upper quasi-triangular matrix (the   
   >    Schur form), and Z is the orthogonal matrix of Schur vectors.   
   >   
   >    Optionally Z may be postmultiplied into an input orthogonal   
   >    matrix Q so that this routine can give the Schur factorization   
   >    of a matrix A which has been reduced to the Hessenberg form H   
   >    by the orthogonal matrix Q:  A = Q*H*Q**T = (QZ)*T*(QZ)**T.   
   > \endverbatim   

    Arguments:   
    ==========   

   > \param[in] JOB   
   > \verbatim   
   >          JOB is CHARACTER*1   
   >           = 'E':  compute eigenvalues only;   
   >           = 'S':  compute eigenvalues and the Schur form T.   
   > \endverbatim   
   >   
   > \param[in] COMPZ   
   > \verbatim   
   >          COMPZ is CHARACTER*1   
   >           = 'N':  no Schur vectors are computed;   
   >           = 'I':  Z is initialized to the unit matrix and the matrix Z   
   >                   of Schur vectors of H is returned;   
   >           = 'V':  Z must contain an orthogonal matrix Q on entry, and   
   >                   the product Q*Z is returned.   
   > \endverbatim   
   >   
   > \param[in] N   
   > \verbatim   
   >          N is INTEGER   
   >           The order of the matrix H.  N .GE. 0.   
   > \endverbatim   
   >   
   > \param[in] ILO   
   > \verbatim   
   >          ILO is INTEGER   
   > \endverbatim   
   >   
   > \param[in] IHI   
   > \verbatim   
   >          IHI is INTEGER   
   >   
   >           It is assumed that H is already upper triangular in rows   
   >           and columns 1:ILO-1 and IHI+1:N. ILO and IHI are normally   
   >           set by a previous call to DGEBAL, and then passed to ZGEHRD   
   >           when the matrix output by DGEBAL is reduced to Hessenberg   
   >           form. Otherwise ILO and IHI should be set to 1 and N   
   >           respectively.  If N.GT.0, then 1.LE.ILO.LE.IHI.LE.N.   
   >           If N = 0, then ILO = 1 and IHI = 0.   
   > \endverbatim   
   >   
   > \param[in,out] H   
   > \verbatim   
   >          H is DOUBLE PRECISION array, dimension (LDH,N)   
   >           On entry, the upper Hessenberg matrix H.   
   >           On exit, if INFO = 0 and JOB = 'S', then H contains the   
   >           upper quasi-triangular matrix T from the Schur decomposition   
   >           (the Schur form); 2-by-2 diagonal blocks (corresponding to   
   >           complex conjugate pairs of eigenvalues) are returned in   
   >           standard form, with H(i,i) = H(i+1,i+1) and   
   >           H(i+1,i)*H(i,i+1).LT.0. If INFO = 0 and JOB = 'E', the   
   >           contents of H are unspecified on exit.  (The output value of   
   >           H when INFO.GT.0 is given under the description of INFO   
   >           below.)   
   >   
   >           Unlike earlier versions of DHSEQR, this subroutine may   
   >           explicitly H(i,j) = 0 for i.GT.j and j = 1, 2, ... ILO-1   
   >           or j = IHI+1, IHI+2, ... N.   
   > \endverbatim   
   >   
   > \param[in] LDH   
   > \verbatim   
   >          LDH is INTEGER   
   >           The leading dimension of the array H. LDH .GE. max(1,N).   
   > \endverbatim   
   >   
   > \param[out] WR   
   > \verbatim   
   >          WR is DOUBLE PRECISION array, dimension (N)   
   > \endverbatim   
   >   
   > \param[out] WI   
   > \verbatim   
   >          WI is DOUBLE PRECISION array, dimension (N)   
   >   
   >           The real and imaginary parts, respectively, of the computed   
   >           eigenvalues. If two eigenvalues are computed as a complex   
   >           conjugate pair, they are stored in consecutive elements of   
   >           WR and WI, say the i-th and (i+1)th, with WI(i) .GT. 0 and   
   >           WI(i+1) .LT. 0. If JOB = 'S', the eigenvalues are stored in   
   >           the same order as on the diagonal of the Schur form returned   
   >           in H, with WR(i) = H(i,i) and, if H(i:i+1,i:i+1) is a 2-by-2   
   >           diagonal block, WI(i) = sqrt(-H(i+1,i)*H(i,i+1)) and   
   >           WI(i+1) = -WI(i).   
   > \endverbatim   
   >   
   > \param[in,out] Z   
   > \verbatim   
   >          Z is DOUBLE PRECISION array, dimension (LDZ,N)   
   >           If COMPZ = 'N', Z is not referenced.   
   >           If COMPZ = 'I', on entry Z need not be set and on exit,   
   >           if INFO = 0, Z contains the orthogonal matrix Z of the Schur   
   >           vectors of H.  If COMPZ = 'V', on entry Z must contain an   
   >           N-by-N matrix Q, which is assumed to be equal to the unit   
   >           matrix except for the submatrix Z(ILO:IHI,ILO:IHI). On exit,   
   >           if INFO = 0, Z contains Q*Z.   
   >           Normally Q is the orthogonal matrix generated by DORGHR   
   >           after the call to DGEHRD which formed the Hessenberg matrix   
   >           H. (The output value of Z when INFO.GT.0 is given under   
   >           the description of INFO below.)   
   > \endverbatim   
   >   
   > \param[in] LDZ   
   > \verbatim   
   >          LDZ is INTEGER   
   >           The leading dimension of the array Z.  if COMPZ = 'I' or   
   >           COMPZ = 'V', then LDZ.GE.MAX(1,N).  Otherwize, LDZ.GE.1.   
   > \endverbatim   
   >   
   > \param[out] WORK   
   > \verbatim   
   >          WORK is DOUBLE PRECISION array, dimension (LWORK)   
   >           On exit, if INFO = 0, WORK(1) returns an estimate of   
   >           the optimal value for LWORK.   
   > \endverbatim   
   >   
   > \param[in] LWORK   
   > \verbatim   
   >          LWORK is INTEGER   
   >           The dimension of the array WORK.  LWORK .GE. max(1,N)   
   >           is sufficient and delivers very good and sometimes   
   >           optimal performance.  However, LWORK as large as 11*N   
   >           may be required for optimal performance.  A workspace   
   >           query is recommended to determine the optimal workspace   
   >           size.   
   >   
   >           If LWORK = -1, then DHSEQR does a workspace query.   
   >           In this case, DHSEQR checks the input parameters and   
   >           estimates the optimal workspace size for the given   
   >           values of N, ILO and IHI.  The estimate is returned   
   >           in WORK(1).  No error message related to LWORK is   
   >           issued by XERBLA.  Neither H nor Z are accessed.   
   > \endverbatim   
   >   
   > \param[out] INFO   
   > \verbatim   
   >          INFO is INTEGER   
   >             =  0:  successful exit   
   >           .LT. 0:  if INFO = -i, the i-th argument had an illegal   
   >                    value   
   >           .GT. 0:  if INFO = i, DHSEQR failed to compute all of   
   >                the eigenvalues.  Elements 1:ilo-1 and i+1:n of WR   
   >                and WI contain those eigenvalues which have been   
   >                successfully computed.  (Failures are rare.)   
   >   
   >                If INFO .GT. 0 and JOB = 'E', then on exit, the   
   >                remaining unconverged eigenvalues are the eigen-   
   >                values of the upper Hessenberg matrix rows and   
   >                columns ILO through INFO of the final, output   
   >                value of H.   
   >   
   >                If INFO .GT. 0 and JOB   = 'S', then on exit   
   >   
   >           (*)  (initial value of H)*U  = U*(final value of H)   
   >   
   >                where U is an orthogonal matrix.  The final   
   >                value of H is upper Hessenberg and quasi-triangular   
   >                in rows and columns INFO+1 through IHI.   
   >   
   >                If INFO .GT. 0 and COMPZ = 'V', then on exit   
   >   
   >                  (final value of Z)  =  (initial value of Z)*U   
   >   
   >                where U is the orthogonal matrix in (*) (regard-   
   >                less of the value of JOB.)   
   >   
   >                If INFO .GT. 0 and COMPZ = 'I', then on exit   
   >                      (final value of Z)  = U   
   >                where U is the orthogonal matrix in (*) (regard-   
   >                less of the value of JOB.)   
   >   
   >                If INFO .GT. 0 and COMPZ = 'N', then Z is not   
   >                accessed.   
   > \endverbatim   

    Authors:   
    ========   

   > \author Univ. of Tennessee   
   > \author Univ. of California Berkeley   
   > \author Univ. of Colorado Denver   
   > \author NAG Ltd.   

   > \date November 2011   

   > \ingroup doubleOTHERcomputational   

   > \par Contributors:   
    ==================   
   >   
   >       Karen Braman and Ralph Byers, Department of Mathematics,   
   >       University of Kansas, USA   

   > \par Further Details:   
    =====================   
   >   
   > \verbatim   
   >   
   >             Default values supplied by   
   >             ILAENV(ISPEC,'DHSEQR',JOB(:1)//COMPZ(:1),N,ILO,IHI,LWORK).   
   >             It is suggested that these defaults be adjusted in order   
   >             to attain best performance in each particular   
   >             computational environment.   
   >   
   >            ISPEC=12: The DLAHQR vs DLAQR0 crossover point.   
   >                      Default: 75. (Must be at least 11.)   
   >   
   >            ISPEC=13: Recommended deflation window size.   
   >                      This depends on ILO, IHI and NS.  NS is the   
   >                      number of simultaneous shifts returned   
   >                      by ILAENV(ISPEC=15).  (See ISPEC=15 below.)   
   >                      The default for (IHI-ILO+1).LE.500 is NS.   
   >                      The default for (IHI-ILO+1).GT.500 is 3*NS/2.   
   >   
   >            ISPEC=14: Nibble crossover point. (See IPARMQ for   
   >                      details.)  Default: 14% of deflation window   
   >                      size.   
   >   
   >            ISPEC=15: Number of simultaneous shifts in a multishift   
   >                      QR iteration.   
   >   
   >                      If IHI-ILO+1 is ...   
   >   
   >                      greater than      ...but less    ... the   
   >                      or equal to ...      than        default is   
   >   
   >                           1               30          NS =   2(+)   
   >                          30               60          NS =   4(+)   
   >                          60              150          NS =  10(+)   
   >                         150              590          NS =  **   
   >                         590             3000          NS =  64   
   >                        3000             6000          NS = 128   
   >                        6000             infinity      NS = 256   
   >   
   >                  (+)  By default some or all matrices of this order   
   >                       are passed to the implicit double shift routine   
   >                       DLAHQR and this parameter is ignored.  See   
   >                       ISPEC=12 above and comments in IPARMQ for   
   >                       details.   
   >   
   >                 (**)  The asterisks (**) indicate an ad-hoc   
   >                       function of N increasing from 10 to 64.   
   >   
   >            ISPEC=16: Select structured matrix multiply.   
   >                      If the number of simultaneous shifts (specified   
   >                      by ISPEC=15) is less than 14, then the default   
   >                      for ISPEC=16 is 0.  Otherwise the default for   
   >                      ISPEC=16 is 2.   
   > \endverbatim   

   > \par References:   
    ================   
   >   
   >       K. Braman, R. Byers and R. Mathias, The Multi-Shift QR   
   >       Algorithm Part I: Maintaining Well Focused Shifts, and Level 3   
   >       Performance, SIAM Journal of Matrix Analysis, volume 23, pages   
   >       929--947, 2002.   
   > \n   
   >       K. Braman, R. Byers and R. Mathias, The Multi-Shift QR   
   >       Algorithm Part II: Aggressive Early Deflation, SIAM Journal   
   >       of Matrix Analysis, volume 23, pages 948--973, 2002.   

    =====================================================================   
   Subroutine */ int igraphdhseqr_(char *job, char *compz, integer *n, integer *ilo,
	 integer *ihi, doublereal *h__, integer *ldh, doublereal *wr, 
	doublereal *wi, doublereal *z__, integer *ldz, doublereal *work, 
	integer *lwork, integer *info)
{
    /* System generated locals */
    address a__1[2];
    integer h_dim1, h_offset, z_dim1, z_offset, i__1, i__2[2], i__3;
    doublereal d__1;
    char ch__1[2];

    /* Builtin functions   
       Subroutine */ void s_cat(char *, char **, integer *, integer *, ftnlen);

    /* Local variables */
    integer i__;
    doublereal hl[2401]	/* was [49][49] */;
    integer kbot, nmin;
    extern logical igraphlsame_(char *, char *);
    logical initz;
    doublereal workl[49];
    logical wantt, wantz;
    extern /* Subroutine */ int igraphdlaqr0_(logical *, logical *, integer *, 
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, integer *), igraphdlahqr_(logical *, logical *,
	     integer *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, integer *), igraphdlacpy_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *), 
	    igraphdlaset_(char *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *);
    extern integer igraphilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    extern /* Subroutine */ int igraphxerbla_(char *, integer *, ftnlen);
    logical lquery;


/*  -- LAPACK computational routine (version 3.4.0) --   
    -- LAPACK is a software package provided by Univ. of Tennessee,    --   
    -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--   
       November 2011   


    =====================================================================   


       ==== Matrices of order NTINY or smaller must be processed by   
       .    DLAHQR because of insufficient subdiagonal scratch space.   
       .    (This is a hard limit.) ====   

       ==== NL allocates some local workspace to help small matrices   
       .    through a rare DLAHQR failure.  NL .GT. NTINY = 11 is   
       .    required and NL .LE. NMIN = ILAENV(ISPEC=12,...) is recom-   
       .    mended.  (The default value of NMIN is 75.)  Using NL = 49   
       .    allows up to six simultaneous shifts and a 16-by-16   
       .    deflation window.  ====   

       ==== Decode and check the input parameters. ====   

       Parameter adjustments */
    h_dim1 = *ldh;
    h_offset = 1 + h_dim1;
    h__ -= h_offset;
    --wr;
    --wi;
    z_dim1 = *ldz;
    z_offset = 1 + z_dim1;
    z__ -= z_offset;
    --work;

    /* Function Body */
    wantt = igraphlsame_(job, "S");
    initz = igraphlsame_(compz, "I");
    wantz = initz || igraphlsame_(compz, "V");
    work[1] = (doublereal) max(1,*n);
    lquery = *lwork == -1;

    *info = 0;
    if (! igraphlsame_(job, "E") && ! wantt) {
	*info = -1;
    } else if (! igraphlsame_(compz, "N") && ! wantz) {
	*info = -2;
    } else if (*n < 0) {
	*info = -3;
    } else if (*ilo < 1 || *ilo > max(1,*n)) {
	*info = -4;
    } else if (*ihi < min(*ilo,*n) || *ihi > *n) {
	*info = -5;
    } else if (*ldh < max(1,*n)) {
	*info = -7;
    } else if (*ldz < 1 || wantz && *ldz < max(1,*n)) {
	*info = -11;
    } else if (*lwork < max(1,*n) && ! lquery) {
	*info = -13;
    }

    if (*info != 0) {

/*        ==== Quick return in case of invalid argument. ==== */

	i__1 = -(*info);
	igraphxerbla_("DHSEQR", &i__1, (ftnlen)6);
	return 0;

    } else if (*n == 0) {

/*        ==== Quick return in case N = 0; nothing to do. ==== */

	return 0;

    } else if (lquery) {

/*        ==== Quick return in case of a workspace query ==== */

	igraphdlaqr0_(&wantt, &wantz, n, ilo, ihi, &h__[h_offset], ldh, &wr[1], &wi[
		1], ilo, ihi, &z__[z_offset], ldz, &work[1], lwork, info);
/*        ==== Ensure reported workspace size is backward-compatible with   
          .    previous LAPACK versions. ====   
   Computing MAX */
	d__1 = (doublereal) max(1,*n);
	work[1] = max(d__1,work[1]);
	return 0;

    } else {

/*        ==== copy eigenvalues isolated by DGEBAL ==== */

	i__1 = *ilo - 1;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    wr[i__] = h__[i__ + i__ * h_dim1];
	    wi[i__] = 0.;
/* L10: */
	}
	i__1 = *n;
	for (i__ = *ihi + 1; i__ <= i__1; ++i__) {
	    wr[i__] = h__[i__ + i__ * h_dim1];
	    wi[i__] = 0.;
/* L20: */
	}

/*        ==== Initialize Z, if requested ==== */

	if (initz) {
	    igraphdlaset_("A", n, n, &c_b11, &c_b12, &z__[z_offset], ldz)
		    ;
	}

/*        ==== Quick return if possible ==== */

	if (*ilo == *ihi) {
	    wr[*ilo] = h__[*ilo + *ilo * h_dim1];
	    wi[*ilo] = 0.;
	    return 0;
	}

/*        ==== DLAHQR/DLAQR0 crossover point ====   

   Writing concatenation */
	i__2[0] = 1, a__1[0] = job;
	i__2[1] = 1, a__1[1] = compz;
	s_cat(ch__1, a__1, i__2, &c__2, (ftnlen)2);
	nmin = igraphilaenv_(&c__12, "DHSEQR", ch__1, n, ilo, ihi, lwork, (ftnlen)6,
		 (ftnlen)2);
	nmin = max(11,nmin);

/*        ==== DLAQR0 for big matrices; DLAHQR for small ones ==== */

	if (*n > nmin) {
	    igraphdlaqr0_(&wantt, &wantz, n, ilo, ihi, &h__[h_offset], ldh, &wr[1], 
		    &wi[1], ilo, ihi, &z__[z_offset], ldz, &work[1], lwork, 
		    info);
	} else {

/*           ==== Small matrix ==== */

	    igraphdlahqr_(&wantt, &wantz, n, ilo, ihi, &h__[h_offset], ldh, &wr[1], 
		    &wi[1], ilo, ihi, &z__[z_offset], ldz, info);

	    if (*info > 0) {

/*              ==== A rare DLAHQR failure!  DLAQR0 sometimes succeeds   
                .    when DLAHQR fails. ==== */

		kbot = *info;

		if (*n >= 49) {

/*                 ==== Larger matrices have enough subdiagonal scratch   
                   .    space to call DLAQR0 directly. ==== */

		    igraphdlaqr0_(&wantt, &wantz, n, ilo, &kbot, &h__[h_offset], 
			    ldh, &wr[1], &wi[1], ilo, ihi, &z__[z_offset], 
			    ldz, &work[1], lwork, info);

		} else {

/*                 ==== Tiny matrices don't have enough subdiagonal   
                   .    scratch space to benefit from DLAQR0.  Hence,   
                   .    tiny matrices must be copied into a larger   
                   .    array before calling DLAQR0. ==== */

		    igraphdlacpy_("A", n, n, &h__[h_offset], ldh, hl, &c__49);
		    hl[*n + 1 + *n * 49 - 50] = 0.;
		    i__1 = 49 - *n;
		    igraphdlaset_("A", &c__49, &i__1, &c_b11, &c_b11, &hl[(*n + 1) *
			     49 - 49], &c__49);
		    igraphdlaqr0_(&wantt, &wantz, &c__49, ilo, &kbot, hl, &c__49, &
			    wr[1], &wi[1], ilo, ihi, &z__[z_offset], ldz, 
			    workl, &c__49, info);
		    if (wantt || *info != 0) {
			igraphdlacpy_("A", n, n, hl, &c__49, &h__[h_offset], ldh);
		    }
		}
	    }
	}

/*        ==== Clear out the trash, if necessary. ==== */

	if ((wantt || *info != 0) && *n > 2) {
	    i__1 = *n - 2;
	    i__3 = *n - 2;
	    igraphdlaset_("L", &i__1, &i__3, &c_b11, &c_b11, &h__[h_dim1 + 3], ldh);
	}

/*        ==== Ensure reported workspace size is backward-compatible with   
          .    previous LAPACK versions. ====   

   Computing MAX */
	d__1 = (doublereal) max(1,*n);
	work[1] = max(d__1,work[1]);
    }

/*     ==== End of DHSEQR ==== */

    return 0;
} /* igraphdhseqr_ */

