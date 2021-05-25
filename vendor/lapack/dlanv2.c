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

static doublereal c_b4 = 1.;

/* > \brief \b DLANV2 computes the Schur factorization of a real 2-by-2 nonsymmetric matrix in standard form. 
  

    =========== DOCUMENTATION ===========   

   Online html documentation available at   
              http://www.netlib.org/lapack/explore-html/   

   > \htmlonly   
   > Download DLANV2 + dependencies   
   > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlanv2.
f">   
   > [TGZ]</a>   
   > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlanv2.
f">   
   > [ZIP]</a>   
   > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlanv2.
f">   
   > [TXT]</a>   
   > \endhtmlonly   

    Definition:   
    ===========   

         SUBROUTINE DLANV2( A, B, C, D, RT1R, RT1I, RT2R, RT2I, CS, SN )   

         DOUBLE PRECISION   A, B, C, CS, D, RT1I, RT1R, RT2I, RT2R, SN   


   > \par Purpose:   
    =============   
   >   
   > \verbatim   
   >   
   > DLANV2 computes the Schur factorization of a real 2-by-2 nonsymmetric   
   > matrix in standard form:   
   >   
   >      [ A  B ] = [ CS -SN ] [ AA  BB ] [ CS  SN ]   
   >      [ C  D ]   [ SN  CS ] [ CC  DD ] [-SN  CS ]   
   >   
   > where either   
   > 1) CC = 0 so that AA and DD are real eigenvalues of the matrix, or   
   > 2) AA = DD and BB*CC < 0, so that AA + or - sqrt(BB*CC) are complex   
   > conjugate eigenvalues.   
   > \endverbatim   

    Arguments:   
    ==========   

   > \param[in,out] A   
   > \verbatim   
   >          A is DOUBLE PRECISION   
   > \endverbatim   
   >   
   > \param[in,out] B   
   > \verbatim   
   >          B is DOUBLE PRECISION   
   > \endverbatim   
   >   
   > \param[in,out] C   
   > \verbatim   
   >          C is DOUBLE PRECISION   
   > \endverbatim   
   >   
   > \param[in,out] D   
   > \verbatim   
   >          D is DOUBLE PRECISION   
   >          On entry, the elements of the input matrix.   
   >          On exit, they are overwritten by the elements of the   
   >          standardised Schur form.   
   > \endverbatim   
   >   
   > \param[out] RT1R   
   > \verbatim   
   >          RT1R is DOUBLE PRECISION   
   > \endverbatim   
   >   
   > \param[out] RT1I   
   > \verbatim   
   >          RT1I is DOUBLE PRECISION   
   > \endverbatim   
   >   
   > \param[out] RT2R   
   > \verbatim   
   >          RT2R is DOUBLE PRECISION   
   > \endverbatim   
   >   
   > \param[out] RT2I   
   > \verbatim   
   >          RT2I is DOUBLE PRECISION   
   >          The real and imaginary parts of the eigenvalues. If the   
   >          eigenvalues are a complex conjugate pair, RT1I > 0.   
   > \endverbatim   
   >   
   > \param[out] CS   
   > \verbatim   
   >          CS is DOUBLE PRECISION   
   > \endverbatim   
   >   
   > \param[out] SN   
   > \verbatim   
   >          SN is DOUBLE PRECISION   
   >          Parameters of the rotation matrix.   
   > \endverbatim   

    Authors:   
    ========   

   > \author Univ. of Tennessee   
   > \author Univ. of California Berkeley   
   > \author Univ. of Colorado Denver   
   > \author NAG Ltd.   

   > \date September 2012   

   > \ingroup doubleOTHERauxiliary   

   > \par Further Details:   
    =====================   
   >   
   > \verbatim   
   >   
   >  Modified by V. Sima, Research Institute for Informatics, Bucharest,   
   >  Romania, to reduce the risk of cancellation errors,   
   >  when computing real eigenvalues, and to ensure, if possible, that   
   >  abs(RT1R) >= abs(RT2R).   
   > \endverbatim   
   >   
    =====================================================================   
   Subroutine */ int igraphdlanv2_(doublereal *a, doublereal *b, doublereal *c__, 
	doublereal *d__, doublereal *rt1r, doublereal *rt1i, doublereal *rt2r,
	 doublereal *rt2i, doublereal *cs, doublereal *sn)
{
    /* System generated locals */
    doublereal d__1, d__2;

    /* Builtin functions */
    double d_sign(doublereal *, doublereal *), sqrt(doublereal);

    /* Local variables */
    doublereal p, z__, aa, bb, cc, dd, cs1, sn1, sab, sac, eps, tau, temp, 
	    scale, bcmax, bcmis, sigma;
    extern doublereal igraphdlapy2_(doublereal *, doublereal *), igraphdlamch_(char *);


/*  -- LAPACK auxiliary routine (version 3.4.2) --   
    -- LAPACK is a software package provided by Univ. of Tennessee,    --   
    -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--   
       September 2012   


    ===================================================================== */


    eps = igraphdlamch_("P");
    if (*c__ == 0.) {
	*cs = 1.;
	*sn = 0.;
	goto L10;

    } else if (*b == 0.) {

/*        Swap rows and columns */

	*cs = 0.;
	*sn = 1.;
	temp = *d__;
	*d__ = *a;
	*a = temp;
	*b = -(*c__);
	*c__ = 0.;
	goto L10;
    } else if (*a - *d__ == 0. && d_sign(&c_b4, b) != d_sign(&c_b4, c__)) {
	*cs = 1.;
	*sn = 0.;
	goto L10;
    } else {

	temp = *a - *d__;
	p = temp * .5;
/* Computing MAX */
	d__1 = abs(*b), d__2 = abs(*c__);
	bcmax = max(d__1,d__2);
/* Computing MIN */
	d__1 = abs(*b), d__2 = abs(*c__);
	bcmis = min(d__1,d__2) * d_sign(&c_b4, b) * d_sign(&c_b4, c__);
/* Computing MAX */
	d__1 = abs(p);
	scale = max(d__1,bcmax);
	z__ = p / scale * p + bcmax / scale * bcmis;

/*        If Z is of the order of the machine accuracy, postpone the   
          decision on the nature of eigenvalues */

	if (z__ >= eps * 4.) {

/*           Real eigenvalues. Compute A and D. */

	    d__1 = sqrt(scale) * sqrt(z__);
	    z__ = p + d_sign(&d__1, &p);
	    *a = *d__ + z__;
	    *d__ -= bcmax / z__ * bcmis;

/*           Compute B and the rotation matrix */

	    tau = igraphdlapy2_(c__, &z__);
	    *cs = z__ / tau;
	    *sn = *c__ / tau;
	    *b -= *c__;
	    *c__ = 0.;
	} else {

/*           Complex eigenvalues, or real (almost) equal eigenvalues.   
             Make diagonal elements equal. */

	    sigma = *b + *c__;
	    tau = igraphdlapy2_(&sigma, &temp);
	    *cs = sqrt((abs(sigma) / tau + 1.) * .5);
	    *sn = -(p / (tau * *cs)) * d_sign(&c_b4, &sigma);

/*           Compute [ AA  BB ] = [ A  B ] [ CS -SN ]   
                     [ CC  DD ]   [ C  D ] [ SN  CS ] */

	    aa = *a * *cs + *b * *sn;
	    bb = -(*a) * *sn + *b * *cs;
	    cc = *c__ * *cs + *d__ * *sn;
	    dd = -(*c__) * *sn + *d__ * *cs;

/*           Compute [ A  B ] = [ CS  SN ] [ AA  BB ]   
                     [ C  D ]   [-SN  CS ] [ CC  DD ] */

	    *a = aa * *cs + cc * *sn;
	    *b = bb * *cs + dd * *sn;
	    *c__ = -aa * *sn + cc * *cs;
	    *d__ = -bb * *sn + dd * *cs;

	    temp = (*a + *d__) * .5;
	    *a = temp;
	    *d__ = temp;

	    if (*c__ != 0.) {
		if (*b != 0.) {
		    if (d_sign(&c_b4, b) == d_sign(&c_b4, c__)) {

/*                    Real eigenvalues: reduce to upper triangular form */

			sab = sqrt((abs(*b)));
			sac = sqrt((abs(*c__)));
			d__1 = sab * sac;
			p = d_sign(&d__1, c__);
			tau = 1. / sqrt((d__1 = *b + *c__, abs(d__1)));
			*a = temp + p;
			*d__ = temp - p;
			*b -= *c__;
			*c__ = 0.;
			cs1 = sab * tau;
			sn1 = sac * tau;
			temp = *cs * cs1 - *sn * sn1;
			*sn = *cs * sn1 + *sn * cs1;
			*cs = temp;
		    }
		} else {
		    *b = -(*c__);
		    *c__ = 0.;
		    temp = *cs;
		    *cs = -(*sn);
		    *sn = temp;
		}
	    }
	}

    }

L10:

/*     Store eigenvalues in (RT1R,RT1I) and (RT2R,RT2I). */

    *rt1r = *a;
    *rt2r = *d__;
    if (*c__ == 0.) {
	*rt1i = 0.;
	*rt2i = 0.;
    } else {
	*rt1i = sqrt((abs(*b))) * sqrt((abs(*c__)));
	*rt2i = -(*rt1i);
    }
    return 0;

/*     End of DLANV2 */

} /* igraphdlanv2_ */

