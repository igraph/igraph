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

/* > \brief \b IPARMQ   

    =========== DOCUMENTATION ===========   

   Online html documentation available at   
              http://www.netlib.org/lapack/explore-html/   

   > \htmlonly   
   > Download IPARMQ + dependencies   
   > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/iparmq.
f">   
   > [TGZ]</a>   
   > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/iparmq.
f">   
   > [ZIP]</a>   
   > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/iparmq.
f">   
   > [TXT]</a>   
   > \endhtmlonly   

    Definition:   
    ===========   

         INTEGER FUNCTION IPARMQ( ISPEC, NAME, OPTS, N, ILO, IHI, LWORK )   

         INTEGER            IHI, ILO, ISPEC, LWORK, N   
         CHARACTER          NAME*( * ), OPTS*( * )   


   > \par Purpose:   
    =============   
   >   
   > \verbatim   
   >   
   >      This program sets problem and machine dependent parameters   
   >      useful for xHSEQR and its subroutines. It is called whenever   
   >      ILAENV is called with 12 <= ISPEC <= 16   
   > \endverbatim   

    Arguments:   
    ==========   

   > \param[in] ISPEC   
   > \verbatim   
   >          ISPEC is integer scalar   
   >              ISPEC specifies which tunable parameter IPARMQ should   
   >              return.   
   >   
   >              ISPEC=12: (INMIN)  Matrices of order nmin or less   
   >                        are sent directly to xLAHQR, the implicit   
   >                        double shift QR algorithm.  NMIN must be   
   >                        at least 11.   
   >   
   >              ISPEC=13: (INWIN)  Size of the deflation window.   
   >                        This is best set greater than or equal to   
   >                        the number of simultaneous shifts NS.   
   >                        Larger matrices benefit from larger deflation   
   >                        windows.   
   >   
   >              ISPEC=14: (INIBL) Determines when to stop nibbling and   
   >                        invest in an (expensive) multi-shift QR sweep.   
   >                        If the aggressive early deflation subroutine   
   >                        finds LD converged eigenvalues from an order   
   >                        NW deflation window and LD.GT.(NW*NIBBLE)/100,   
   >                        then the next QR sweep is skipped and early   
   >                        deflation is applied immediately to the   
   >                        remaining active diagonal block.  Setting   
   >                        IPARMQ(ISPEC=14) = 0 causes TTQRE to skip a   
   >                        multi-shift QR sweep whenever early deflation   
   >                        finds a converged eigenvalue.  Setting   
   >                        IPARMQ(ISPEC=14) greater than or equal to 100   
   >                        prevents TTQRE from skipping a multi-shift   
   >                        QR sweep.   
   >   
   >              ISPEC=15: (NSHFTS) The number of simultaneous shifts in   
   >                        a multi-shift QR iteration.   
   >   
   >              ISPEC=16: (IACC22) IPARMQ is set to 0, 1 or 2 with the   
   >                        following meanings.   
   >                        0:  During the multi-shift QR sweep,   
   >                            xLAQR5 does not accumulate reflections and   
   >                            does not use matrix-matrix multiply to   
   >                            update the far-from-diagonal matrix   
   >                            entries.   
   >                        1:  During the multi-shift QR sweep,   
   >                            xLAQR5 and/or xLAQRaccumulates reflections and uses   
   >                            matrix-matrix multiply to update the   
   >                            far-from-diagonal matrix entries.   
   >                        2:  During the multi-shift QR sweep.   
   >                            xLAQR5 accumulates reflections and takes   
   >                            advantage of 2-by-2 block structure during   
   >                            matrix-matrix multiplies.   
   >                        (If xTRMM is slower than xGEMM, then   
   >                        IPARMQ(ISPEC=16)=1 may be more efficient than   
   >                        IPARMQ(ISPEC=16)=2 despite the greater level of   
   >                        arithmetic work implied by the latter choice.)   
   > \endverbatim   
   >   
   > \param[in] NAME   
   > \verbatim   
   >          NAME is character string   
   >               Name of the calling subroutine   
   > \endverbatim   
   >   
   > \param[in] OPTS   
   > \verbatim   
   >          OPTS is character string   
   >               This is a concatenation of the string arguments to   
   >               TTQRE.   
   > \endverbatim   
   >   
   > \param[in] N   
   > \verbatim   
   >          N is integer scalar   
   >               N is the order of the Hessenberg matrix H.   
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
   >               It is assumed that H is already upper triangular   
   >               in rows and columns 1:ILO-1 and IHI+1:N.   
   > \endverbatim   
   >   
   > \param[in] LWORK   
   > \verbatim   
   >          LWORK is integer scalar   
   >               The amount of workspace available.   
   > \endverbatim   

    Authors:   
    ========   

   > \author Univ. of Tennessee   
   > \author Univ. of California Berkeley   
   > \author Univ. of Colorado Denver   
   > \author NAG Ltd.   

   > \date November 2011   

   > \ingroup auxOTHERauxiliary   

   > \par Further Details:   
    =====================   
   >   
   > \verbatim   
   >   
   >       Little is known about how best to choose these parameters.   
   >       It is possible to use different values of the parameters   
   >       for each of CHSEQR, DHSEQR, SHSEQR and ZHSEQR.   
   >   
   >       It is probably best to choose different parameters for   
   >       different matrices and different parameters at different   
   >       times during the iteration, but this has not been   
   >       implemented --- yet.   
   >   
   >   
   >       The best choices of most of the parameters depend   
   >       in an ill-understood way on the relative execution   
   >       rate of xLAQR3 and xLAQR5 and on the nature of each   
   >       particular eigenvalue problem.  Experiment may be the   
   >       only practical way to determine which choices are most   
   >       effective.   
   >   
   >       Following is a list of default values supplied by IPARMQ.   
   >       These defaults may be adjusted in order to attain better   
   >       performance in any particular computational environment.   
   >   
   >       IPARMQ(ISPEC=12) The xLAHQR vs xLAQR0 crossover point.   
   >                        Default: 75. (Must be at least 11.)   
   >   
   >       IPARMQ(ISPEC=13) Recommended deflation window size.   
   >                        This depends on ILO, IHI and NS, the   
   >                        number of simultaneous shifts returned   
   >                        by IPARMQ(ISPEC=15).  The default for   
   >                        (IHI-ILO+1).LE.500 is NS.  The default   
   >                        for (IHI-ILO+1).GT.500 is 3*NS/2.   
   >   
   >       IPARMQ(ISPEC=14) Nibble crossover point.  Default: 14.   
   >   
   >       IPARMQ(ISPEC=15) Number of simultaneous shifts, NS.   
   >                        a multi-shift QR iteration.   
   >   
   >                        If IHI-ILO+1 is ...   
   >   
   >                        greater than      ...but less    ... the   
   >                        or equal to ...      than        default is   
   >   
   >                                0               30       NS =   2+   
   >                               30               60       NS =   4+   
   >                               60              150       NS =  10   
   >                              150              590       NS =  **   
   >                              590             3000       NS =  64   
   >                             3000             6000       NS = 128   
   >                             6000             infinity   NS = 256   
   >   
   >                    (+)  By default matrices of this order are   
   >                         passed to the implicit double shift routine   
   >                         xLAHQR.  See IPARMQ(ISPEC=12) above.   These   
   >                         values of NS are used only in case of a rare   
   >                         xLAHQR failure.   
   >   
   >                    (**) The asterisks (**) indicate an ad-hoc   
   >                         function increasing from 10 to 64.   
   >   
   >       IPARMQ(ISPEC=16) Select structured matrix multiply.   
   >                        (See ISPEC=16 above for details.)   
   >                        Default: 3.   
   > \endverbatim   
   >   
    ===================================================================== */
integer igraphiparmq_(integer *ispec, char *name__, char *opts, integer *n, integer 
	*ilo, integer *ihi, integer *lwork)
{
    /* System generated locals */
    integer ret_val, i__1, i__2;
    real r__1;

    /* Builtin functions */
    double log(doublereal);
    integer i_nint(real *);

    /* Local variables */
    integer nh, ns;


/*  -- LAPACK auxiliary routine (version 3.4.0) --   
    -- LAPACK is a software package provided by Univ. of Tennessee,    --   
    -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--   
       November 2011   


    ================================================================ */
    if (*ispec == 15 || *ispec == 13 || *ispec == 16) {

/*        ==== Set the number simultaneous shifts ==== */

	nh = *ihi - *ilo + 1;
	ns = 2;
	if (nh >= 30) {
	    ns = 4;
	}
	if (nh >= 60) {
	    ns = 10;
	}
	if (nh >= 150) {
/* Computing MAX */
	    r__1 = log((real) nh) / log(2.f);
	    i__1 = 10, i__2 = nh / i_nint(&r__1);
	    ns = max(i__1,i__2);
	}
	if (nh >= 590) {
	    ns = 64;
	}
	if (nh >= 3000) {
	    ns = 128;
	}
	if (nh >= 6000) {
	    ns = 256;
	}
/* Computing MAX */
	i__1 = 2, i__2 = ns - ns % 2;
	ns = max(i__1,i__2);
    }

    if (*ispec == 12) {


/*        ===== Matrices of order smaller than NMIN get sent   
          .     to xLAHQR, the classic double shift algorithm.   
          .     This must be at least 11. ==== */

	ret_val = 75;

    } else if (*ispec == 14) {

/*        ==== INIBL: skip a multi-shift qr iteration and   
          .    whenever aggressive early deflation finds   
          .    at least (NIBBLE*(window size)/100) deflations. ==== */

	ret_val = 14;

    } else if (*ispec == 15) {

/*        ==== NSHFTS: The number of simultaneous shifts ===== */

	ret_val = ns;

    } else if (*ispec == 13) {

/*        ==== NW: deflation window size.  ==== */

	if (nh <= 500) {
	    ret_val = ns;
	} else {
	    ret_val = ns * 3 / 2;
	}

    } else if (*ispec == 16) {

/*        ==== IACC22: Whether to accumulate reflections   
          .     before updating the far-from-diagonal elements   
          .     and whether to use 2-by-2 block structure while   
          .     doing it.  A small amount of work could be saved   
          .     by making this choice dependent also upon the   
          .     NH=IHI-ILO+1. */

	ret_val = 0;
	if (ns >= 14) {
	    ret_val = 1;
	}
	if (ns >= 14) {
	    ret_val = 2;
	}

    } else {
/*        ===== invalid value of ispec ===== */
	ret_val = -1;

    }

/*     ==== End of IPARMQ ==== */

    return ret_val;
} /* igraphiparmq_ */

