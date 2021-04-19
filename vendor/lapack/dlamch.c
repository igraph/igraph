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

static doublereal c_b2 = 0.;

/* > \brief \b DLAMCH   

    =========== DOCUMENTATION ===========   

   Online html documentation available at   
              http://www.netlib.org/lapack/explore-html/   

    Definition:   
    ===========   

        DOUBLE PRECISION FUNCTION DLAMCH( CMACH )   


   > \par Purpose:   
    =============   
   >   
   > \verbatim   
   >   
   > DLAMCH determines double precision machine parameters.   
   > \endverbatim   

    Arguments:   
    ==========   

   > \param[in] CMACH   
   > \verbatim   
   >          Specifies the value to be returned by DLAMCH:   
   >          = 'E' or 'e',   DLAMCH := eps   
   >          = 'S' or 's ,   DLAMCH := sfmin   
   >          = 'B' or 'b',   DLAMCH := base   
   >          = 'P' or 'p',   DLAMCH := eps*base   
   >          = 'N' or 'n',   DLAMCH := t   
   >          = 'R' or 'r',   DLAMCH := rnd   
   >          = 'M' or 'm',   DLAMCH := emin   
   >          = 'U' or 'u',   DLAMCH := rmin   
   >          = 'L' or 'l',   DLAMCH := emax   
   >          = 'O' or 'o',   DLAMCH := rmax   
   >          where   
   >          eps   = relative machine precision   
   >          sfmin = safe minimum, such that 1/sfmin does not overflow   
   >          base  = base of the machine   
   >          prec  = eps*base   
   >          t     = number of (base) digits in the mantissa   
   >          rnd   = 1.0 when rounding occurs in addition, 0.0 otherwise   
   >          emin  = minimum exponent before (gradual) underflow   
   >          rmin  = underflow threshold - base**(emin-1)   
   >          emax  = largest exponent before overflow   
   >          rmax  = overflow threshold  - (base**emax)*(1-eps)   
   > \endverbatim   

    Authors:   
    ========   

   > \author Univ. of Tennessee   
   > \author Univ. of California Berkeley   
   > \author Univ. of Colorado Denver   
   > \author NAG Ltd.   

   > \date November 2011   

   > \ingroup auxOTHERauxiliary   

    ===================================================================== */
doublereal igraphdlamch_(char *cmach)
{
    /* System generated locals */
    doublereal ret_val;

    /* Local variables */
    extern doublereal radixdbl_(doublereal *), digitsdbl_(doublereal *), 
	    epsilondbl_(doublereal *);
    doublereal rnd, eps, rmach;
    extern logical igraphlsame_(char *, char *);
    doublereal small, sfmin;
    extern integer minexponentdbl_(doublereal *), maxexponentdbl_(doublereal *
	    );
    extern doublereal hugedbl_(doublereal *), tinydbl_(doublereal *);


/*  -- LAPACK auxiliary routine (version 3.4.0) --   
    -- LAPACK is a software package provided by Univ. of Tennessee,    --   
    -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--   
       November 2011   



   =====================================================================   



       Assume rounding, not chopping. Always. */

    rnd = 1.;

    if (1. == rnd) {
	eps = epsilondbl_(&c_b2) * .5f;
    } else {
	eps = epsilondbl_(&c_b2);
    }

    if (igraphlsame_(cmach, "E")) {
	rmach = eps;
    } else if (igraphlsame_(cmach, "S")) {
	sfmin = tinydbl_(&c_b2);
	small = 1. / hugedbl_(&c_b2);
	if (small >= sfmin) {

/*           Use SMALL plus a bit, to avoid the possibility of rounding   
             causing overflow when computing  1/sfmin. */

	    sfmin = small * (eps + 1.);
	}
	rmach = sfmin;
    } else if (igraphlsame_(cmach, "B")) {
	rmach = radixdbl_(&c_b2);
    } else if (igraphlsame_(cmach, "P")) {
	rmach = eps * radixdbl_(&c_b2);
    } else if (igraphlsame_(cmach, "N")) {
	rmach = digitsdbl_(&c_b2);
    } else if (igraphlsame_(cmach, "R")) {
	rmach = rnd;
    } else if (igraphlsame_(cmach, "M")) {
	rmach = (doublereal) minexponentdbl_(&c_b2);
    } else if (igraphlsame_(cmach, "U")) {
	rmach = tinydbl_(&c_b2);
    } else if (igraphlsame_(cmach, "L")) {
	rmach = (doublereal) maxexponentdbl_(&c_b2);
    } else if (igraphlsame_(cmach, "O")) {
	rmach = hugedbl_(&c_b2);
    } else {
	rmach = 0.;
    }

    ret_val = rmach;
    return ret_val;

/*     End of DLAMCH */

} /* igraphdlamch_   

   ***********************************************************************   
   > \brief \b DLAMC3   
   > \details   
   > \b Purpose:   
   > \verbatim   
   > DLAMC3  is intended to force  A  and  B  to be stored prior to doing   
   > the addition of  A  and  B ,  for use in situations where optimizers   
   > might hold one of these in a register.   
   > \endverbatim   
   > \author LAPACK is a software package provided by Univ. of Tennessee, Univ. of California Berkeley, Univ. 
of Colorado Denver and NAG Ltd..   
   > \date November 2011   
   > \ingroup auxOTHERauxiliary   
   >   
   > \param[in] A   
   > \verbatim   
   >          A is a DOUBLE PRECISION   
   > \endverbatim   
   >   
   > \param[in] B   
   > \verbatim   
   >          B is a DOUBLE PRECISION   
   >          The values A and B.   
   > \endverbatim   
   > */
doublereal igraphdlamc3_(doublereal *a, doublereal *b)
{
    /* System generated locals */
    doublereal ret_val;


/*  -- LAPACK auxiliary routine (version 3.4.0) --   
       Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..   
       November 2010   

   ===================================================================== */


    ret_val = *a + *b;

    return ret_val;

/*     End of DLAMC3 */

} /* igraphdlamc3_ */

