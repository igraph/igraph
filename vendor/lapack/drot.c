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

/* > \brief \b DROT   

    =========== DOCUMENTATION ===========   

   Online html documentation available at   
              http://www.netlib.org/lapack/explore-html/   

    Definition:   
    ===========   

         SUBROUTINE DROT(N,DX,INCX,DY,INCY,C,S)   

         DOUBLE PRECISION C,S   
         INTEGER INCX,INCY,N   
         DOUBLE PRECISION DX(*),DY(*)   


   > \par Purpose:   
    =============   
   >   
   > \verbatim   
   >   
   >    DROT applies a plane rotation.   
   > \endverbatim   

    Arguments:   
    ==========   

   > \param[in] N   
   > \verbatim   
   >          N is INTEGER   
   >         number of elements in input vector(s)   
   > \endverbatim   
   >   
   > \param[in,out] DX   
   > \verbatim   
   >          DX is DOUBLE PRECISION array, dimension ( 1 + ( N - 1 )*abs( INCX ) )   
   > \endverbatim   
   >   
   > \param[in] INCX   
   > \verbatim   
   >          INCX is INTEGER   
   >         storage spacing between elements of DX   
   > \endverbatim   
   >   
   > \param[in,out] DY   
   > \verbatim   
   >          DY is DOUBLE PRECISION array, dimension ( 1 + ( N - 1 )*abs( INCY ) )   
   > \endverbatim   
   >   
   > \param[in] INCY   
   > \verbatim   
   >          INCY is INTEGER   
   >         storage spacing between elements of DY   
   > \endverbatim   
   >   
   > \param[in] C   
   > \verbatim   
   >          C is DOUBLE PRECISION   
   > \endverbatim   
   >   
   > \param[in] S   
   > \verbatim   
   >          S is DOUBLE PRECISION   
   > \endverbatim   

    Authors:   
    ========   

   > \author Univ. of Tennessee   
   > \author Univ. of California Berkeley   
   > \author Univ. of Colorado Denver   
   > \author NAG Ltd.   

   > \date November 2017   

   > \ingroup double_blas_level1   

   > \par Further Details:   
    =====================   
   >   
   > \verbatim   
   >   
   >     jack dongarra, linpack, 3/11/78.   
   >     modified 12/3/93, array(1) declarations changed to array(*)   
   > \endverbatim   
   >   
    =====================================================================   
   Subroutine */ int igraphdrot_(integer *n, doublereal *dx, integer *incx, 
	doublereal *dy, integer *incy, doublereal *c__, doublereal *s)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    integer i__, ix, iy;
    doublereal dtemp;


/*  -- Reference BLAS level1 routine (version 3.8.0) --   
    -- Reference BLAS is a software package provided by Univ. of Tennessee,    --   
    -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--   
       November 2017   


    =====================================================================   

       Parameter adjustments */
    --dy;
    --dx;

    /* Function Body */
    if (*n <= 0) {
	return 0;
    }
    if (*incx == 1 && *incy == 1) {

/*       code for both increments equal to 1 */

	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    dtemp = *c__ * dx[i__] + *s * dy[i__];
	    dy[i__] = *c__ * dy[i__] - *s * dx[i__];
	    dx[i__] = dtemp;
	}
    } else {

/*       code for unequal increments or equal increments not equal   
           to 1 */

	ix = 1;
	iy = 1;
	if (*incx < 0) {
	    ix = (-(*n) + 1) * *incx + 1;
	}
	if (*incy < 0) {
	    iy = (-(*n) + 1) * *incy + 1;
	}
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    dtemp = *c__ * dx[ix] + *s * dy[iy];
	    dy[iy] = *c__ * dy[iy] - *s * dx[ix];
	    dx[ix] = dtemp;
	    ix += *incx;
	    iy += *incy;
	}
    }
    return 0;
} /* igraphdrot_ */

