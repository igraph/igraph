/* igraphdlanv2.f -- translated by f2c (version 19991025).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"
#include "config.h"
#include "arpack_internal.h"

/* Table of constant values */

static doublereal c_b3 = 1.;

/* Subroutine */ int igraphdlanv2_(a, b, c__, d__, rt1r, rt1i, rt2r, rt2i, cs, sn)
doublereal *a, *b, *c__, *d__, *rt1r, *rt1i, *rt2r, *rt2i, *cs, *sn;
{
    /* System generated locals */
    doublereal d__1;

    /* Builtin functions */
    double igraphd_sign(), sqrt();

    /* Local variables */
    static doublereal temp, p, sigma;
    extern doublereal igraphdlapy2_();
    static doublereal aa, bb, cc, dd, cs1, sn1, sab, sac, tau;


/*  -- LAPACK auxiliary routine (version 2.0) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
/*     Courant Institute, Argonne National Lab, and Rice University */
/*     September 30, 1994 */

/*     .. Scalar Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  DLANV2 computes the Schur factorization of a real 2-by-2 nonsymmetric */
/*  matrix in standard form: */

/*       [ A  B ] = [ CS -SN ] [ AA  BB ] [ CS  SN ] */
/*       [ C  D ]   [ SN  CS ] [ CC  DD ] [-SN  CS ] */

/*  where either */
/*  1) CC = 0 so that AA and DD are real eigenvalues of the matrix, or */
/*  2) AA = DD and BB*CC < 0, so that AA + or - sqrt(BB*CC) are complex */
/*  conjugate eigenvalues. */

/*  Arguments */
/*  ========= */

/*  A       (input/output) DOUBLE PRECISION */
/*  B       (input/output) DOUBLE PRECISION */
/*  C       (input/output) DOUBLE PRECISION */
/*  D       (input/output) DOUBLE PRECISION */
/*          On entry, the elements of the input matrix. */
/*          On exit, they are overwritten by the elements of the */
/*          standardised Schur form. */

/*  RT1R    (output) DOUBLE PRECISION */
/*  RT1I    (output) DOUBLE PRECISION */
/*  RT2R    (output) DOUBLE PRECISION */
/*  RT2I    (output) DOUBLE PRECISION */
/*          The real and imaginary parts of the eigenvalues. If the */
/*          eigenvalues are both real, abs(RT1R) >= abs(RT2R); if the */
/*          eigenvalues are a complex conjugate pair, RT1I > 0. */

/*  CS      (output) DOUBLE PRECISION */
/*  SN      (output) DOUBLE PRECISION */
/*          Parameters of the rotation matrix. */

/*  ===================================================================== */

/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Initialize CS and SN */

    *cs = 1.;
    *sn = 0.;

    if (*c__ == 0.) {
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
    } else if (*a - *d__ == 0. && igraphd_sign(&c_b3, b) != igraphd_sign(&c_b3, c__)) {
	goto L10;
    } else {

/*        Make diagonal elements equal */

	temp = *a - *d__;
	p = temp * .5;
	sigma = *b + *c__;
	tau = igraphdlapy2_(&sigma, &temp);
	cs1 = sqrt((abs(sigma) / tau + 1.) * .5);
	sn1 = -(p / (tau * cs1)) * igraphd_sign(&c_b3, &sigma);

/*        Compute [ AA  BB ] = [ A  B ] [ CS1 -SN1 ] */
/*                [ CC  DD ]   [ C  D ] [ SN1  CS1 ] */

	aa = *a * cs1 + *b * sn1;
	bb = -(*a) * sn1 + *b * cs1;
	cc = *c__ * cs1 + *d__ * sn1;
	dd = -(*c__) * sn1 + *d__ * cs1;

/*        Compute [ A  B ] = [ CS1  SN1 ] [ AA  BB ] */
/*                [ C  D ]   [-SN1  CS1 ] [ CC  DD ] */

	*a = aa * cs1 + cc * sn1;
	*b = bb * cs1 + dd * sn1;
	*c__ = -aa * sn1 + cc * cs1;
	*d__ = -bb * sn1 + dd * cs1;

/*        Accumulate transformation */

	temp = *cs * cs1 - *sn * sn1;
	*sn = *cs * sn1 + *sn * cs1;
	*cs = temp;

	temp = (*a + *d__) * .5;
	*a = temp;
	*d__ = temp;

	if (*c__ != 0.) {
	    if (*b != 0.) {
		if (igraphd_sign(&c_b3, b) == igraphd_sign(&c_b3, c__)) {

/*                 Real eigenvalues: reduce to upper triangular form */

		    sab = sqrt((abs(*b)));
		    sac = sqrt((abs(*c__)));
		    d__1 = sab * sac;
		    p = igraphd_sign(&d__1, c__);
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

