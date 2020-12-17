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

/* Subroutine */ int igraphsecond_(real *t)
{
    real t1;
    extern doublereal etime_(real *);
    real tarray[2];



/*  -- LAPACK auxiliary routine (preliminary version) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       July 26, 1991   

    Purpose   
    =======   

    SECOND returns the user time for a process in seconds.   
    This version gets the time from the system function ETIME. */


    t1 = etime_(tarray);
    *t = tarray[0];
    return 0;

/*     End of SECOND */

} /* igraphsecond_ */

