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

/* > \brief \b DLAISNAN tests input for NaN by comparing two arguments for inequality.   

    =========== DOCUMENTATION ===========   

   Online html documentation available at   
              http://www.netlib.org/lapack/explore-html/   

   > \htmlonly   
   > Download DLAISNAN + dependencies   
   > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlaisna
n.f">   
   > [TGZ]</a>   
   > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlaisna
n.f">   
   > [ZIP]</a>   
   > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlaisna
n.f">   
   > [TXT]</a>   
   > \endhtmlonly   

    Definition:   
    ===========   

         LOGICAL FUNCTION DLAISNAN( DIN1, DIN2 )   

         DOUBLE PRECISION   DIN1, DIN2   


   > \par Purpose:   
    =============   
   >   
   > \verbatim   
   >   
   > This routine is not for general use.  It exists solely to avoid   
   > over-optimization in DISNAN.   
   >   
   > DLAISNAN checks for NaNs by comparing its two arguments for   
   > inequality.  NaN is the only floating-point value where NaN != NaN   
   > returns .TRUE.  To check for NaNs, pass the same variable as both   
   > arguments.   
   >   
   > A compiler must assume that the two arguments are   
   > not the same variable, and the test will not be optimized away.   
   > Interprocedural or whole-program optimization may delete this   
   > test.  The ISNAN functions will be replaced by the correct   
   > Fortran 03 intrinsic once the intrinsic is widely available.   
   > \endverbatim   

    Arguments:   
    ==========   

   > \param[in] DIN1   
   > \verbatim   
   >          DIN1 is DOUBLE PRECISION   
   > \endverbatim   
   >   
   > \param[in] DIN2   
   > \verbatim   
   >          DIN2 is DOUBLE PRECISION   
   >          Two numbers to compare for inequality.   
   > \endverbatim   

    Authors:   
    ========   

   > \author Univ. of Tennessee   
   > \author Univ. of California Berkeley   
   > \author Univ. of Colorado Denver   
   > \author NAG Ltd.   

   > \date September 2012   

   > \ingroup auxOTHERauxiliary   

    ===================================================================== */
logical igraphdlaisnan_(doublereal *din1, doublereal *din2)
{
    /* System generated locals */
    logical ret_val;


/*  -- LAPACK auxiliary routine (version 3.4.2) --   
    -- LAPACK is a software package provided by Univ. of Tennessee,    --   
    -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--   
       September 2012   


    ===================================================================== */

    ret_val = *din1 != *din2;
    return ret_val;
} /* igraphdlaisnan_ */

